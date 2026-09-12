"""Map the host platform to the Galacticus release assets it needs.

The CI ``Deploy`` job publishes one executable and one tools archive per
platform to a GitHub release.  The names are fixed strings (there is no
platform suffix scheme to parse), so we map ``(system, machine)`` to them
explicitly.  Anything we do not recognise raises :class:`UnsupportedPlatform`
with an actionable message rather than guessing.
"""

import platform
from collections import namedtuple

# Describes the release assets for one platform.
#   binary        -- name of the executable asset (e.g. "Galacticus.exe").
#   tools         -- name of the run-time tools archive asset.
#   tools_format  -- "tar.zst", "tar.bz2", or "zip"; how to unpack `tools`.
#   tools_legacy  -- name of the tools archive published by releases predating
#                    the switch to zstd, or None if there is no earlier name.
#   tools_legacy_format -- how to unpack `tools_legacy`.
#   key           -- short human label for the platform (used in messages).
#   retired       -- None for a platform still built, or a message explaining
#                    that it is not, used when a release turns out to publish no
#                    binary for it (see `download._require_published_binary`).
#
# Two tools archives are named per platform because a release only ever carries
# the format that was current when it was cut: releases published before the
# switch have `tools_legacy` and nothing else, and those assets can not be
# regenerated after the fact.  Provisioning therefore prefers `tools` and falls
# back to `tools_legacy` (see `download._tools_archive`), which keeps every
# already-published version tag installable.
PlatformAssets = namedtuple(
    "PlatformAssets",
    ["binary", "tools", "tools_format", "tools_legacy", "tools_legacy_format", "key",
     "retired"],
    defaults=(None,),
)

# The last released version whose assets include a macOS Intel build.
MACOS_INTEL_FINAL_VERSION = "0.9.12"

# Why no newer release carries a macOS Intel binary.  A retired platform keeps
# its entry above so that the releases which *do* carry its assets stay
# installable; only a release which publishes none is refused, and then with
# this message rather than a bare download failure.
_MACOS_INTEL_RETIRED = (
    "Galacticus no longer builds macOS Intel (x86-64) binaries: GitHub Actions is "
    "retiring its Intel macOS runners, and Homebrew no longer publishes bottles "
    f"for that platform. Releases up to and including v{MACOS_INTEL_FINAL_VERSION} "
    f"still carry one -- `pip install 'galacticus=={MACOS_INTEL_FINAL_VERSION}'` "
    "installs it, and an existing install keeps working. Otherwise build from "
    "source: https://galacticus.readthedocs.io/en/latest/manuals/user-guide/"
    "installation/source-macos.html"
)


class UnsupportedPlatform(RuntimeError):
    """Raised when no pre-built binary exists for the host platform."""


# Local source builds always produce an executable called "Galacticus.exe"
# regardless of platform; only the *released* asset names differ.
LOCAL_BINARY_NAME = "Galacticus.exe"


def detect(system=None, machine=None):
    """Return the :class:`PlatformAssets` for the host (or the given override).

    `system`/`machine` default to :func:`platform.system` /
    :func:`platform.machine` and are accepted as arguments so the mapping can
    be unit-tested without monkey-patching.
    """
    system = (system if system is not None else platform.system()).strip()
    machine = (machine if machine is not None else platform.machine()).strip().lower()

    if system == "Linux":
        if machine in ("x86_64", "amd64"):
            return PlatformAssets("Galacticus.exe",
                                  "tools.tar.zst", "tar.zst",
                                  "tools.tar.bz2", "tar.bz2",
                                  "Linux x86-64")
        raise UnsupportedPlatform(
            f"Galacticus provides no pre-built Linux binary for machine '{machine}'. "
            "Build from source: https://galacticus.readthedocs.io/en/latest/"
            "manuals/user-guide/installation/source-linux.html"
        )
    if system == "Darwin":
        if machine in ("x86_64", "amd64"):
            return PlatformAssets("Galacticus_MacOS.exe",
                                  "toolsMacOS.tar.zst", "tar.zst",
                                  "toolsMacOS.zip", "zip",
                                  "macOS x86-64",
                                  _MACOS_INTEL_RETIRED)
        if machine in ("arm64", "aarch64"):
            return PlatformAssets("Galacticus_MacOS-M1.exe",
                                  "toolsMacOSM1.tar.zst", "tar.zst",
                                  "toolsMacOSM1.zip", "zip",
                                  "macOS Apple Silicon")
        raise UnsupportedPlatform(
            f"Galacticus provides no pre-built macOS binary for machine '{machine}'."
        )
    raise UnsupportedPlatform(
        f"Galacticus provides no pre-built binary for system '{system}'. "
        "On Windows, run `galacticus install-wsl` to set up WSL 2 and install the "
        "Linux build inside it. Otherwise build from source: "
        "https://galacticus.readthedocs.io/en/latest/manuals/user-guide/"
        "installation/index.html"
    )
