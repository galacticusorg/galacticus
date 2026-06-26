"""``galacticus`` command-line entry point.

Sub-commands:

* ``install`` / ``update`` -- provision (or refresh) the managed install.
* ``run <params.xml> [args...]`` -- validate then dispatch to the executable.
* ``validate <params.xml>`` -- run validation only.
* ``clean`` -- purge the regenerable cache (never the durable install).
* ``info`` -- report the resolved install, paths, and environment.

``galacticus <params.xml>`` is accepted as shorthand for ``galacticus run``.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import platformdirs

from . import __version__, download, paths, platforms, validate as _validate

_COMMANDS = {"install", "update", "run", "validate", "clean", "info"}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Shorthand: `galacticus params.xml` -> `galacticus run params.xml`.
    if argv and argv[0] not in _COMMANDS and not argv[0].startswith("-"):
        argv = ["run"] + argv

    parser = argparse.ArgumentParser(
        prog="galacticus",
        description="Install and run the Galacticus galaxy-formation model.",
    )
    parser.add_argument("--version", action="version", version=f"galacticus {__version__}")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("install", help="download the binary, datasets, and tools")
    sub.add_parser("update", help="re-download the install for the current version")

    run_parser = sub.add_parser("run", help="run a parameter file")
    run_parser.add_argument("parameter_file")
    run_parser.add_argument("--no-validate", action="store_true",
                            help="skip pre-dispatch parameter validation")
    run_parser.add_argument("extra", nargs=argparse.REMAINDER,
                            help="arguments passed through to Galacticus.exe")

    validate_parser = sub.add_parser("validate", help="validate a parameter file")
    validate_parser.add_argument("parameter_file")
    validate_parser.add_argument("--structural", action="store_true",
                                 help="also run structural checks")

    clean_parser = sub.add_parser("clean", help="purge the regenerable cache")
    clean_parser.add_argument("--all", action="store_true",
                              help="remove the entire dynamic data cache")
    clean_parser.add_argument("--older-than", type=float, metavar="N",
                              help="only remove cache files older than N days")
    clean_parser.add_argument("--prune-installs", action="store_true",
                              help="also remove superseded per-version install dirs")
    clean_parser.add_argument("--force", action="store_true",
                              help="clean even when pre-built tools share the "
                                   "dynamic path (deletes them too)")
    clean_parser.add_argument("--dry-run", action="store_true",
                              help="report what would be freed; delete nothing")

    sub.add_parser("info", help="show the resolved install and environment")

    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 0

    try:
        return _dispatch(args)
    except platforms.UnsupportedPlatform as error:
        print(f"galacticus: {error}", file=sys.stderr)
        return 2
    except (RuntimeError, ValueError) as error:
        print(f"galacticus: {error}", file=sys.stderr)
        return 1


def _dispatch(args):
    if args.command == "clean":
        return _cmd_clean(args)
    if args.command == "info":
        return _cmd_info()

    install = paths.resolve()
    if args.command in ("install", "update"):
        return _cmd_install(install, force=(args.command == "update"))
    if args.command == "validate":
        _ensure(install)
        return _cmd_validate(install, args.parameter_file, args.structural)
    if args.command == "run":
        _ensure(install)
        return _cmd_run(install, args)
    raise ValueError(f"unknown command {args.command!r}")  # pragma: no cover


def _ensure(install):
    """Make sure `install` is runnable, provisioning a managed one on demand."""
    if install.managed:
        download.provision(install)
    elif install.binary is None or not Path(install.binary).is_file():
        raise RuntimeError(
            f"no Galacticus executable found for the {install.source} install at "
            f"{install.exec_path}. Build it, or unset GALACTICUS_EXEC_PATH/"
            "GALACTICUS_HOME to use a managed download."
        )


def _cmd_install(install, *, force):
    if not install.managed:
        print(f"Using existing {install.source} install at {install.exec_path}; "
              "nothing to download.")
        return 0
    done = download.provision(install, force=force)
    if done:
        print("Provisioned: " + ", ".join(done))
    else:
        print("Install already up to date.")
    print(f"Install location: {install.exec_path.parent}")
    return 0


def _cmd_validate(install, parameter_file, structural):
    result = _validate.validate(parameter_file, install, structural=structural)
    _report_findings(result)
    if result.ok:
        print(f"OK ({result.method}): {parameter_file}")
        return 0
    print(f"INVALID ({result.method}): {parameter_file}", file=sys.stderr)
    return 1


def _cmd_run(install, args):
    if not args.no_validate:
        result = _validate.validate(args.parameter_file, install)
        _report_findings(result)
        if not result.ok:
            print(f"galacticus: refusing to run; {args.parameter_file} failed "
                  "validation (use --no-validate to override).", file=sys.stderr)
            return 1
    if install.data_path is None:
        print("galacticus: warning: no datasets path resolved; the run may fail. "
              "Set GALACTICUS_DATA_PATH.", file=sys.stderr)
    env = install.environ()
    extra = [a for a in (args.extra or []) if a != "--"]
    command = [str(install.binary), str(args.parameter_file), *extra]
    sys.stdout.flush()
    os.execve(str(install.binary), command, env)  # replaces this process


def _cmd_info():
    install = paths.resolve()
    print(f"galacticus launcher {__version__}")
    print(f"install source : {install.source}")
    if install.tag:
        print(f"release tag    : {install.tag}")
    print(f"executable     : {install.binary}"
          f"{'' if install.binary and Path(install.binary).is_file() else '  (not present)'}")
    # Show the paths the binary will actually use. For a managed install these
    # are imposed by the launcher; for env/home they reflect the resolved
    # effective locations (dynamic defaults to <data>/dynamic when unset).
    rows = [
        ("GALACTICUS_EXEC_PATH", install.exec_path),
        ("GALACTICUS_DATA_PATH", install.data_path),
        ("GALACTICUS_TOOLS_PATH", install.tools_path),
        ("GALACTICUS_DYNAMIC_DATA_PATH", install.dynamic_path),
    ]
    for name, value in rows:
        suffix = ""
        if name == "GALACTICUS_TOOLS_PATH" and install.tools_path is None:
            value = "(unset; tools share the dynamic path)"
        elif name == "GALACTICUS_DYNAMIC_DATA_PATH" and not install.managed \
                and os.environ.get("GALACTICUS_DYNAMIC_DATA_PATH") is None:
            suffix = "  (default <data>/dynamic)"
        print(f"{name:<28} = {value if value is not None else '(unset)'}{suffix}")
    if install.dynamic_path is not None:
        size = _tree_size(install.dynamic_path)
        note = "" if install.tools_separated else "  (includes pre-built tools)"
        print(f"dynamic data   : {_human(size)} ({install.dynamic_path}){note}")
    return 0


def _cmd_clean(args):
    install = paths.resolve()
    target = install.dynamic_path
    if target is None or not Path(target).exists():
        print("Nothing to clean (no dynamic data path resolved).")
        return 0
    # Guard: when tools are NOT on a separate path they live under the dynamic
    # path, and purging it would delete pre-built tools a binary-only install
    # cannot rebuild. Refuse unless the user explicitly forces it.
    if not install.tools_separated and not args.force:
        size = _tree_size(target)
        print(
            f"galacticus: {target} ({_human(size)}) holds regenerable data AND "
            "pre-built tools (no separate GALACTICUS_TOOLS_PATH is set).\n"
            "Refusing to clean to avoid deleting tools that a binary-only install "
            "cannot rebuild. Set GALACTICUS_TOOLS_PATH to a separate location to "
            "keep tools out of the cache, or pass --force to clean everything here.",
            file=sys.stderr,
        )
        return 1

    cutoff = None
    if args.older_than is not None:
        cutoff = time.time() - args.older_than * 86400.0
    freed, removed = _purge_tree(Path(target), cutoff=cutoff, dry_run=args.dry_run)
    verb = "Would free" if args.dry_run else "Freed"
    print(f"{verb} {_human(freed)} from {target} ({removed} item(s)).")

    if args.prune_installs and install.managed:
        data_base = Path(platformdirs.user_data_dir(paths.APP_NAME))
        keep = install.tag
        for child in (data_base.iterdir() if data_base.is_dir() else []):
            if child.is_dir() and child.name != keep:
                size = _tree_size(child)
                print(f"{verb} {_human(size)} from install {child.name}.")
                if not args.dry_run:
                    _rmtree(child)
        print(f"Kept current install: {keep}")
    return 0


# --- helpers ---------------------------------------------------------------

def _report_findings(result):
    for finding in result.findings:
        print(f"  [{finding.level}/{finding.kind}] {finding.path}: {finding.message}")


def _purge_tree(base, *, cutoff, dry_run):
    """Delete files under `base` (older than `cutoff` if given). Returns
    ``(bytes_freed, items_removed)``.  Never deletes `base` itself."""
    if not base.is_dir():
        return 0, 0
    freed = 0
    removed = 0
    if cutoff is None:
        for child in base.iterdir():
            size = _tree_size(child)
            freed += size
            removed += 1
            if not dry_run:
                _rmtree(child) if child.is_dir() else child.unlink()
        return freed, removed
    # Age-based: walk files, remove those older than cutoff, then empty dirs.
    for dirpath, _dirnames, filenames in os.walk(base):
        for name in filenames:
            file_path = Path(dirpath) / name
            try:
                stat = file_path.stat()
            except OSError:
                continue
            if stat.st_mtime < cutoff:
                freed += stat.st_size
                removed += 1
                if not dry_run:
                    try:
                        file_path.unlink()
                    except OSError:
                        pass
    if not dry_run:
        _remove_empty_dirs(base)
    return freed, removed


def _remove_empty_dirs(base):
    for dirpath, dirnames, filenames in os.walk(base, topdown=False):
        if Path(dirpath) == base:
            continue
        if not dirnames and not filenames:
            try:
                Path(dirpath).rmdir()
            except OSError:
                pass


def _rmtree(path):
    import shutil
    shutil.rmtree(path, ignore_errors=True)


def _tree_size(path):
    if not Path(path).exists():
        return 0
    if Path(path).is_file():
        return Path(path).stat().st_size
    total = 0
    for dirpath, _dirnames, filenames in os.walk(path):
        for name in filenames:
            try:
                total += (Path(dirpath) / name).stat().st_size
            except OSError:
                pass
    return total


def _human(num_bytes):
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
