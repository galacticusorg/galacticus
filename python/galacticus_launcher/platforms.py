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
#   binary       -- name of the executable asset (e.g. "Galacticus.exe").
#   tools        -- name of the run-time tools archive asset.
#   tools_format -- "tar.bz2" or "zip"; how to unpack `tools`.
#   key          -- short human label for the platform (used in messages).
PlatformAssets = namedtuple(
    "PlatformAssets", ["binary", "tools", "tools_format", "key"]
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
            return PlatformAssets("Galacticus.exe", "tools.tar.bz2", "tar.bz2", "Linux x86-64")
        raise UnsupportedPlatform(
            f"Galacticus provides no pre-built Linux binary for machine '{machine}'. "
            "Build from source: https://galacticus.readthedocs.io/en/latest/"
            "manuals/user-guide/installation/source-linux.html"
        )
    if system == "Darwin":
        if machine in ("x86_64", "amd64"):
            return PlatformAssets("Galacticus_MacOS.exe", "toolsMacOS.zip", "zip", "macOS x86-64")
        if machine in ("arm64", "aarch64"):
            return PlatformAssets("Galacticus_MacOS-M1.exe", "toolsMacOSM1.zip", "zip", "macOS Apple Silicon")
        raise UnsupportedPlatform(
            f"Galacticus provides no pre-built macOS binary for machine '{machine}'."
        )
    raise UnsupportedPlatform(
        f"Galacticus provides no pre-built binary for system '{system}'. "
        "On Windows, use WSL 2 and install the Linux build. Otherwise build "
        "from source: https://galacticus.readthedocs.io/en/latest/manuals/"
        "user-guide/installation/index.html"
    )
