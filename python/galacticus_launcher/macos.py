"""macOS binary compatibility check.

The pre-built macOS binaries are built on a recent macOS (currently macOS 15)
with no lowered deployment target, so the Mach-O minimum-OS load command pins the
oldest macOS they will run on. Running one on an older host fails with a cryptic
``dyld`` error, so before dispatching we read that minimum directly from the
binary (pure Python -- no Xcode tools) and refuse with a clear message when the
host macOS is older.

Only the Mach-O header and load commands are read (a few KiB near the start), not
the whole multi-hundred-MB binary.
"""

import platform
import struct

# Mach-O magics (thin) and fat/universal magics (always big-endian on disk).
_MH_MAGIC_64 = 0xFEEDFACF
_MH_MAGIC_32 = 0xFEEDFACE
_FAT_MAGIC = 0xCAFEBABE
_FAT_MAGIC_64 = 0xCAFEBABF

# Load commands carrying a minimum-OS, and the macOS platform id.
_LC_VERSION_MIN_MACOSX = 0x24
_LC_BUILD_VERSION = 0x32
_PLATFORM_MACOS = 1

_DOCS = ("https://galacticus.readthedocs.io/en/latest/manuals/user-guide/"
         "installation/source-macos.html")


def _decode_version(value):
    """Decode an X.Y.Z nibble-packed Mach-O version into ``(major, minor)``."""
    return (value >> 16, (value >> 8) & 0xFF)


def _u32(fh, offset, endian="<"):
    fh.seek(offset)
    data = fh.read(4)
    if len(data) < 4:
        raise struct.error("short read")
    return struct.unpack(endian + "I", data)[0]


def parse_minimum_macos(path):
    """Return the binary's minimum macOS as ``(major, minor)``, or None if it is
    not a Mach-O / has no version load command / cannot be read."""
    try:
        with open(path, "rb") as fh:
            magic_be = _u32(fh, 0, ">")
            if magic_be in (_FAT_MAGIC, _FAT_MAGIC_64):
                return _parse_fat(fh, magic_be)
            return _parse_thin(fh, 0)
    except (OSError, struct.error):
        return None


def _parse_fat(fh, magic):
    """Scan the slices of a universal binary; return the first slice's minimum."""
    nfat = _u32(fh, 4, ">")
    wide = magic == _FAT_MAGIC_64
    entry_size = 32 if wide else 20  # fat_arch_64 vs fat_arch
    for i in range(nfat):
        base = 8 + i * entry_size
        # offset field follows cputype(4) + cpusubtype(4).
        if wide:
            fh.seek(base + 8)
            offset = struct.unpack(">Q", fh.read(8))[0]
        else:
            offset = _u32(fh, base + 8, ">")
        result = _parse_thin(fh, offset)
        if result is not None:
            return result
    return None


def _parse_thin(fh, base):
    raw = _u32(fh, base, "<")
    if raw in (_MH_MAGIC_64, _MH_MAGIC_32):
        endian, magic = "<", raw
    else:
        magic = _u32(fh, base, ">")
        if magic not in (_MH_MAGIC_64, _MH_MAGIC_32):
            return None
        endian = ">"
    is64 = magic == _MH_MAGIC_64
    ncmds = _u32(fh, base + 16, endian)          # field 4 of the mach_header
    offset = base + (32 if is64 else 28)         # header size (64- vs 32-bit)
    for _ in range(ncmds):
        cmd = _u32(fh, offset, endian)
        cmdsize = _u32(fh, offset + 4, endian)
        if cmdsize == 0:
            break
        if cmd == _LC_BUILD_VERSION:
            if _u32(fh, offset + 8, endian) == _PLATFORM_MACOS:   # platform
                return _decode_version(_u32(fh, offset + 12, endian))  # minos
        elif cmd == _LC_VERSION_MIN_MACOSX:
            return _decode_version(_u32(fh, offset + 8, endian))  # version
        offset += cmdsize
    return None


def host_version():
    """Return the host macOS version as ``(major, minor)``, or None."""
    release = platform.mac_ver()[0]
    if not release:
        return None
    parts = release.split(".")
    try:
        return (int(parts[0]), int(parts[1]) if len(parts) > 1 else 0)
    except ValueError:
        return None


def incompatible_reason(binary_path):
    """Return a human-readable message if the host macOS is older than
    `binary_path`'s minimum, else None.

    A no-op (None) on non-macOS hosts and whenever the minimum or host version
    cannot be determined -- the guard only fires on a definite mismatch.
    """
    if platform.system() != "Darwin":
        return None
    minimum = parse_minimum_macos(binary_path)
    host = host_version()
    if minimum is None or host is None or host >= minimum:
        return None
    return (
        f"this pre-built Galacticus binary requires macOS {minimum[0]}.{minimum[1]} "
        f"or newer, but this host is macOS {host[0]}.{host[1]}. Upgrade macOS, or "
        f"build Galacticus from source: {_DOCS}"
    )
