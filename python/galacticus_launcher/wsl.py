"""Windows support: drive Galacticus through WSL 2.

There is no native Windows build of Galacticus (the model, its build system, and
its run-time tools all assume a POSIX environment), but the Linux build runs
under the Windows Subsystem for Linux (WSL 2).  On a Windows host the launcher
therefore delegates every command to a copy of itself installed inside a WSL
distribution, so that ``pip install galacticus`` on Windows gives the same
one-command experience as on Linux and macOS:

* ``galacticus install-wsl`` walks the user through installing WSL 2 (an
  administrator prompt and a reboot are needed when WSL itself is missing),
  installing the default Ubuntu distribution, and provisioning Galacticus inside
  it.  It always asks before changing the system unless ``--yes`` is passed.
* Every other command (``run``, ``install``, ``validate``, ...) is forwarded to
  the launcher inside the distribution, with Windows paths on the command line
  translated to their ``/mnt/<drive>/...`` equivalents.  The current directory
  carries over automatically (``wsl.exe`` starts in the translated directory),
  so output files land where they would on Linux.

The copy inside the distribution lives in a dedicated virtual environment
(``~/.local/share/galacticus/venv``), and the binary, datasets, and tools it
downloads live on the Linux filesystem.  That placement matters: file access
across the Windows/Linux boundary is far too slow for the datasets and the HDF5
output.

Only ``wsl.exe`` and ``powershell.exe`` are used, both of which ship with
Windows.  ``wsl.exe`` writes its own messages as UTF-16, while the output of a
Linux command run through it is passed through as UTF-8, so :func:`decode`
handles both.
"""

import os
import platform
import re
import shutil
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

from . import __version__

WSL_EXE = "wsl.exe"
POWERSHELL_EXE = "powershell.exe"

# The distribution `install-wsl` installs when none is present. Ubuntu is what
# `wsl --install` itself defaults to, and is the distribution the Linux binary
# is tested against in CI.
DEFAULT_DISTRIBUTION = "Ubuntu"

# Location of the launcher's virtual environment inside the distribution,
# relative to the Linux user's home directory.
VENV = ".local/share/galacticus/venv"

# Distributions that are the backing store of a container runtime rather than a
# usable Linux system. They are often the *default* distribution, so they are
# skipped explicitly rather than trusting the default marker.
_UTILITY_DISTRIBUTIONS = frozenset({
    "docker-desktop", "docker-desktop-data",
    "rancher-desktop", "rancher-desktop-data",
    "podman-machine-default",
})

# Launcher environment variables forwarded into the distribution (via WSLENV),
# so that e.g. `GALACTICUS_RELEASE_TAG=... galacticus run` behaves the same on
# Windows. The path variables are deliberately NOT forwarded: Windows paths mean
# nothing to the launcher inside the distribution.
_FORWARDED = ("GALACTICUS_RELEASE_TAG", "GALACTICUS_DATASETS_REF")

_DOCS = ("https://galacticus.readthedocs.io/en/latest/manuals/user-guide/"
         "installation/pip.html#windows")

# Well-known WSL failure codes and what to do about them.
_HINTS = (
    ("0x80370102", "virtualization is disabled: enable it (Intel VT-x / AMD-V) "
                   "in the firmware (BIOS/UEFI) settings, or, in a virtual "
                   "machine, enable nested virtualization"),
    ("0x80370114", "the Virtual Machine Platform feature is not enabled; run "
                   "`galacticus install-wsl` again (or `wsl --install`) as an "
                   "administrator"),
    ("0x800701bc", "the WSL 2 kernel is missing or outdated; run `wsl --update` "
                   "and try again"),
)

STATE_UNSUPPORTED = "unsupported"          # not Windows, or a Windows without wsl.exe
STATE_NOT_INSTALLED = "not-installed"      # wsl.exe present, WSL feature not installed
STATE_NO_DISTRIBUTION = "no-distribution"  # WSL installed, no usable distribution
STATE_READY = "ready"                      # a usable distribution exists

Status = namedtuple("Status", ["state", "distribution", "version", "detail"])
Distribution = namedtuple("Distribution", ["name", "state", "version", "default"])
Result = namedtuple("Result", ["returncode", "output"])


def is_windows():
    return platform.system() == "Windows"


# --- running wsl.exe -------------------------------------------------------

def decode(data):
    """Decode output from ``wsl.exe``.

    Messages from ``wsl.exe`` itself (``--list``, ``--status``, errors) are
    UTF-16LE; output of a Linux command run through it is UTF-8. UTF-16 text of
    any length has NUL bytes, so that is the discriminator.
    """
    if not data:
        return ""
    if b"\x00" in data:
        text = data.decode("utf-16-le", errors="replace")
    else:
        text = data.decode("utf-8", errors="replace")
    return text.lstrip("\ufeff").replace("\x00", "").replace("\r\n", "\n")


def _environment():
    """The environment for a ``wsl.exe`` child: forward the launcher's release
    overrides into the distribution through ``WSLENV``."""
    env = dict(os.environ)
    names = [name for name in _FORWARDED if env.get(name)]
    if names:
        current = [entry for entry in env.get("WSLENV", "").split(":") if entry]
        env["WSLENV"] = ":".join(current + [n for n in names if n not in current])
    return env


def _run(args, *, capture=True):
    """Run ``wsl.exe args...``. Returns :class:`Result`; output is empty when
    not captured (the child then shares this process's console)."""
    completed = subprocess.run([WSL_EXE, *args], capture_output=capture,
                               env=_environment())
    output = ""
    if capture:
        output = decode(completed.stdout) + decode(completed.stderr)
    return Result(completed.returncode, output)


def run(distribution, command, *, user=None, capture=True):
    """Run `command` (an argv list, no shell) inside `distribution`."""
    args = ["-d", distribution]
    if user:
        args += ["-u", user]
    args += ["--exec", *command]
    return _run(args, capture=capture)


def launcher_command(argv):
    """The argv that runs the launcher inside the distribution.

    ``$HOME`` is only known to the Linux side, so a shell expands it; the
    launcher's own arguments are passed as positional parameters and forwarded
    with ``"$@"``, which keeps every argument intact (no re-quoting).
    """
    script = f'exec "$HOME/{VENV}/bin/galacticus" "$@"'
    return ["/bin/sh", "-c", script, "galacticus", *argv]


# --- status ----------------------------------------------------------------

def parse_distribution_list(text):
    """Parse ``wsl --list --verbose`` output.

    The table looks like::

          NAME              STATE           VERSION
        * Ubuntu            Running         2
          docker-desktop    Stopped         2

    where ``*`` marks the default distribution.
    """
    distributions = []
    for line in text.splitlines():
        default = line.lstrip().startswith("*")
        fields = line.lstrip(" *").split()
        if len(fields) < 3 or fields[0].upper() == "NAME" or not fields[2].isdigit():
            continue
        distributions.append(Distribution(fields[0], fields[1], int(fields[2]), default))
    return distributions


def list_distributions():
    result = _run(["--list", "--verbose"])
    if result.returncode != 0:
        return []
    return parse_distribution_list(result.output)


def status():
    """Classify the host: see the ``STATE_*`` constants.

    The distribution reported is the one the launcher will use: the default
    one, unless that is a container runtime's utility distribution, in which
    case the first ordinary one.
    """
    if not is_windows():
        return Status(STATE_UNSUPPORTED, None, None, "not a Windows host")
    if shutil.which(WSL_EXE) is None:
        return Status(STATE_UNSUPPORTED, None, None,
                      "wsl.exe was not found. WSL 2 needs Windows 10 (build 19041 "
                      "or later) or Windows 11.")
    usable = [d for d in list_distributions() if d.name not in _UTILITY_DISTRIBUTIONS]
    if usable:
        chosen = next((d for d in usable if d.default), usable[0])
        return Status(STATE_READY, chosen.name, chosen.version, "")
    # No usable distribution. `wsl --status` succeeds once the WSL feature is
    # installed (it reports the default version), and fails on the stub that
    # Windows ships before it is.
    result = _run(["--status"])
    state = STATE_NO_DISTRIBUTION if result.returncode == 0 else STATE_NOT_INSTALLED
    return Status(state, None, None, result.output.strip())


def hint(output):
    """A remedy for a well-known WSL error code appearing in `output`, or None."""
    for code, remedy in _HINTS:
        if code in output.lower():
            return remedy
    return None


# --- path translation ------------------------------------------------------

def _path_and_suffix(argument):
    """Split a command-line path into (existing path to translate, suffix).

    ``wslpath`` resolves through the filesystem, so a not-yet-existing output
    file is translated via its (existing) parent directory, with the file name
    re-appended afterwards.
    """
    if os.path.exists(argument):
        return os.path.abspath(argument), ""
    parent, name = os.path.split(os.path.abspath(argument))
    if name and os.path.isdir(parent):
        return parent, "/" + name
    return None, ""


def translate_paths(distribution, windows_paths):
    """Translate Windows paths to their paths inside `distribution`.

    A single ``wslpath`` invocation handles the whole list (one process start),
    honoring however the distribution mounts the Windows drives.
    """
    if not windows_paths:
        return []
    script = 'for p in "$@"; do wslpath -a "$p" || exit 1; done'
    result = run(distribution, ["/bin/sh", "-c", script, "sh", *windows_paths])
    lines = [line for line in result.output.splitlines() if line.strip()]
    if result.returncode != 0 or len(lines) != len(windows_paths):
        raise RuntimeError(
            f"could not translate paths for WSL distribution {distribution}: "
            f"{result.output.strip() or 'no output from wslpath'}")
    return lines


def translate_arguments(distribution, argv):
    """Return `argv` (sub-command first) with Windows paths rewritten.

    Positional arguments naming an existing file or directory are translated,
    as is the value of ``-o``/``--output`` (which need not exist yet). Options
    and anything that does not resolve on the Windows side pass through
    unchanged, so the launcher inside the distribution resolves bundled paths
    such as ``parameters/quickTest.xml`` itself.
    """
    argv = list(argv)
    targets = {}   # index -> (path, suffix)
    expect_output = False
    for index, argument in enumerate(argv[1:], start=1):
        if expect_output:
            expect_output = False
            path, suffix = _path_and_suffix(argument)
        elif argument in ("-o", "--output"):
            expect_output = True
            continue
        elif argument.startswith("-") or not os.path.exists(argument):
            continue
        else:
            path, suffix = os.path.abspath(argument), ""
        if path is not None:
            targets[index] = (path, suffix)
    translated = translate_paths(distribution, [p for p, _ in targets.values()])
    for (index, (_, suffix)), wsl_path in zip(targets.items(), translated):
        argv[index] = wsl_path + suffix
    return argv


# --- provisioning inside the distribution ----------------------------------

def is_release(version=None):
    version = __version__ if version is None else version
    return version != "0.0.0" and re.fullmatch(r"\d+\.\d+\.\d+", version) is not None


def source_tree():
    """The repository this launcher was installed from in development mode
    (``pip install -e .``), or None for a normal (wheel) install."""
    root = Path(__file__).resolve().parents[2]
    return root if (root / "pyproject.toml").is_file() else None


def pip_target(distribution, version=None, tree=None):
    """Arguments telling ``pip install`` inside the distribution what to
    install: the same released version as this launcher; for a development
    install, the same source checkout (editable), else the newest release."""
    version = __version__ if version is None else version
    if is_release(version):
        return [f"galacticus=={version}"]
    tree = source_tree() if tree is None else tree
    if tree is not None:
        return ["-e", translate_paths(distribution, [str(tree)])[0]]
    return ["galacticus"]


# `pip` refuses to install into the system Python of recent Ubuntu releases, so
# the launcher gets its own virtual environment. A venv left half-made by an
# earlier failure (no pip) is rebuilt.
_PROVISION_SCRIPT = f'''
set -e
venv="$HOME/{VENV}"
if ! "$venv/bin/python" -m pip --version >/dev/null 2>&1; then
    rm -rf "$venv"
    mkdir -p "$(dirname "$venv")"
    python3 -m venv "$venv"
fi
"$venv/bin/python" -m pip install --quiet --upgrade "$@"
'''

# Ubuntu's WSL image lacks the `venv` module (and, on some images, Python
# itself). Installed as root, which WSL allows without a password.
_PACKAGES_SCRIPT = "apt-get update -qq && apt-get install -y -qq python3 python3-venv"


def installed(distribution):
    """True if the launcher inside `distribution` exists and, for a released
    version, matches this launcher's version."""
    result = run(distribution, launcher_command(["--version"]))
    if result.returncode != 0:
        return False
    if not is_release():
        return True
    return result.output.strip().endswith(__version__)


def provision(distribution, *, log=print):
    """Install (or update) the launcher inside `distribution`. Returns an exit code."""
    target = pip_target(distribution)
    log(f"Installing Galacticus ({' '.join(target)}) in WSL distribution {distribution} ...")
    result = run(distribution, ["/bin/sh", "-c", _PROVISION_SCRIPT, "sh", *target],
                 capture=False)
    if result.returncode != 0:
        log("Installing the Python packages the launcher needs (python3, "
            "python3-venv) ...")
        result = run(distribution, ["/bin/sh", "-c", _PACKAGES_SCRIPT], user="root",
                     capture=False)
        if result.returncode == 0:
            result = run(distribution, ["/bin/sh", "-c", _PROVISION_SCRIPT, "sh", *target],
                         capture=False)
    if result.returncode != 0:
        print(f"galacticus: installing the launcher in WSL distribution "
              f"{distribution} failed (exit code {result.returncode}). See {_DOCS}",
              file=sys.stderr)
    return result.returncode


# --- delegation ------------------------------------------------------------

def _not_ready_message(current):
    if current.state == STATE_UNSUPPORTED:
        return f"galacticus: {current.detail} See {_DOCS}"
    if current.state == STATE_NOT_INSTALLED:
        what = "WSL 2 is not installed on this Windows machine"
    else:
        what = "WSL 2 is installed but has no Linux distribution"
    return (f"galacticus: {what}. Galacticus runs on Windows through WSL 2; run "
            f"`galacticus install-wsl` to set it up. See {_DOCS}")


def dispatch(argv, *, log=print):
    """Run the launcher command line `argv` inside WSL. Returns an exit code."""
    current = status()
    if current.state != STATE_READY:
        print(_not_ready_message(current), file=sys.stderr)
        return 2
    if current.version != 2:
        print(f"galacticus: WSL distribution {current.distribution} runs under WSL "
              f"{current.version}, but Galacticus needs WSL 2. Run `galacticus "
              f"install-wsl` to convert it. See {_DOCS}", file=sys.stderr)
        return 2
    if not installed(current.distribution):
        code = provision(current.distribution, log=log)
        if code != 0:
            return code
    if argv and argv[0] == "info":
        log(f"Windows host: commands run inside WSL distribution "
            f"{current.distribution} (WSL {current.version}).")
    try:
        argv = translate_arguments(current.distribution, argv)
    except RuntimeError as error:
        print(f"galacticus: {error}", file=sys.stderr)
        return 1
    sys.stdout.flush()
    return run(current.distribution, launcher_command(argv), capture=False).returncode


# --- `galacticus install-wsl` ----------------------------------------------

def _confirm(question, yes):
    if yes:
        return True
    if not sys.stdin.isatty():
        print("galacticus: not running in a terminal; pass --yes to proceed without "
              "confirmation.", file=sys.stderr)
        return False
    return input(f"{question} [y/N] ").strip().lower() in ("y", "yes")


def enable_wsl():
    """Install the WSL feature (no distribution). Needs administrator rights, so
    it is started elevated through PowerShell, which shows the standard Windows
    consent prompt; the elevated ``wsl.exe`` runs in its own console window.
    Returns its exit code (non-zero also if the prompt was declined)."""
    script = ("$process = Start-Process -FilePath 'wsl.exe' "
              "-ArgumentList '--install','--no-distribution' "
              "-Verb RunAs -Wait -PassThru; exit $process.ExitCode")
    completed = subprocess.run([POWERSHELL_EXE, "-NoProfile", "-NonInteractive",
                                "-Command", script])
    return completed.returncode


def install_distribution(name=DEFAULT_DISTRIBUTION):
    """Install `name` and run its first-launch setup interactively in this
    console (it asks for a Linux user name and password, then opens a shell).
    Returns the exit code."""
    return _run(["--install", "-d", name], capture=False).returncode


def install(*, yes=False, tools=None, log=print):
    """The ``galacticus install-wsl`` command. Returns an exit code.

    Idempotent: each run does the next step that is missing (install WSL; reboot;
    install a distribution; provision the launcher; download the Galacticus
    assets), so after a reboot the user simply runs it again.
    """
    if not is_windows():
        log("galacticus install-wsl is only needed on Windows; nothing to do.")
        return 0
    current = status()
    if current.state == STATE_UNSUPPORTED:
        print(f"galacticus: {current.detail} See {_DOCS}", file=sys.stderr)
        return 2

    if current.state == STATE_NOT_INSTALLED:
        log("Galacticus runs on Windows through WSL 2 (the Windows Subsystem for "
            "Linux), which is not installed on this machine.\n"
            "Installing it enables the Virtual Machine Platform Windows feature "
            "and the WSL 2 kernel. Windows will ask for administrator approval, "
            "and a reboot is needed afterwards.")
        if not _confirm("Install WSL 2 now?", yes):
            return 1
        code = enable_wsl()
        if code != 0:
            print(f"galacticus: installing WSL did not complete (exit code {code}). "
                  f"If the administrator prompt was declined, run `galacticus "
                  f"install-wsl` again; otherwise see {_DOCS}", file=sys.stderr)
            return code
        log("WSL 2 is installed. Reboot Windows, then run `galacticus install-wsl` "
            "again to add the Linux distribution and Galacticus.")
        return 0

    if current.state == STATE_NO_DISTRIBUTION:
        log(f"WSL 2 is installed but has no Linux distribution. The next step "
            f"downloads {DEFAULT_DISTRIBUTION} (about 1 GB) and starts it once to "
            f"create your Linux user: it will ask for a user name and password, "
            f"then show a Linux prompt. Type `exit` at that prompt to continue.")
        if not _confirm(f"Install {DEFAULT_DISTRIBUTION} now?", yes):
            return 1
        code = install_distribution()
        current = status()
        if code != 0 or current.state != STATE_READY:
            remedy = hint(current.detail)
            print(f"galacticus: installing {DEFAULT_DISTRIBUTION} did not complete "
                  f"(exit code {code})." + (f" Likely cause: {remedy}." if remedy else "")
                  + f" See {_DOCS}", file=sys.stderr)
            return code or 1

    if current.version != 2:
        log(f"WSL distribution {current.distribution} runs under WSL "
            f"{current.version}; Galacticus needs WSL 2. Converting it can take "
            f"several minutes.")
        if not _confirm(f"Convert {current.distribution} to WSL 2 now?", yes):
            return 1
        result = _run(["--set-version", current.distribution, "2"], capture=False)
        if result.returncode != 0:
            print(f"galacticus: `wsl --set-version {current.distribution} 2` failed "
                  f"(exit code {result.returncode}). See {_DOCS}", file=sys.stderr)
            return result.returncode

    code = provision(current.distribution, log=log)
    if code != 0:
        return code
    # Fetch the binary, datasets and tools now, so the first `galacticus run` is
    # immediate; the tri-state --tools choice is passed straight through.
    inner = ["install"]
    if tools is True:
        inner.append("--tools")
    elif tools is False:
        inner.append("--no-tools")
    code = run(current.distribution, launcher_command(inner), capture=False).returncode
    if code == 0:
        log(f"Galacticus is ready in WSL distribution {current.distribution}. "
            f"Run a model with `galacticus run <parameter file>` from any Windows "
            f"command prompt.")
    return code
