# Processes `sourceDigest` directives: emits a C-interop character array
# declaration bound to the per-build-target MD5 hash symbol (defined by
# the Perl build system's `Find_Hash` chain), plus the matching
# `use :: ISO_C_Binding, only : C_Char` import.
#
# Also provides `find_hash`, `hash_data_files`, `modification_time`, and
# `update_modification_time` — the Makefile-side helpers used by
# scripts/build/sourceDigests.* to compute per-file and composite MD5
# digests (ported from the same Perl module).
#
# Andrew Benson (ported to Python 2026)
#
# Mirrors `perl/Galacticus/Build/SourceTree/Process/SourceDigest.pm`.

import base64
import fcntl
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from build.fortran_utils                            import (
    get_matching_lines, read_file,
)
from Galacticus.Build.SourceTree                    import walk_tree
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import add_declarations
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses


# ---------------------------------------------------------------------------
# Shared state (mirrors the `our %digests`, `%compositeDigests`,
# `%modificationTimes` package globals in the Perl module).
# ---------------------------------------------------------------------------

_digests            = {}
_composite_digests  = {}
_modification_times = {}


def process_source_digests(tree, options):
    """Mirrors Process_SourceDigests() from SourceDigest.pm:24-58."""
    for node in walk_tree(tree):
        if node.get('type') != 'sourceDigest':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        name = directive['name']
        add_declarations(node['parent'], [{
            'intrinsic':     'character',
            'type':          'C_Char',
            'openMP':        False,
            'attributes':    [
                'dimension(23)',
                f'bind(C, name="{name}MD5")',
            ],
            'variables':     [name],
            'variableNames': [name],
        }])
        add_uses(node['parent'], {
            'moduleUse': {
                'ISO_C_Binding': {
                    'intrinsic': True,
                    'only':      {'C_Char': True},
                },
            },
            'moduleOrder': ['ISO_C_Binding'],
        })


register_process('sourceDigests', process_source_digests)


# ---------------------------------------------------------------------------
# Hash-computation helpers (ported from Find_Hash / Hash_Data_Files /
# modificationTime / updateModificationTime in SourceDigest.pm).
# ---------------------------------------------------------------------------

def _b64digest(hasher):
    """Return the 22-character base64 form of `hasher`, matching Perl's
    Digest::MD5::b64digest (no trailing '=' padding).
    """
    return base64.b64encode(hasher.digest()).rstrip(b'=').decode('ascii')


def modification_time(file_name):
    """Return the mtime of `file_name`, caching repeated lookups.

    Mirrors `modificationTime()` at SourceDigest.pm:287-293.  Perl's
    `stat()` on a missing file returns an empty list, quietly caching
    `undef`; do the same here by caching `None`.
    """
    if file_name not in _modification_times:
        try:
            _modification_times[file_name] = os.stat(file_name).st_mtime
        except FileNotFoundError:
            _modification_times[file_name] = None
    return _modification_times[file_name]


def update_modification_time(file_name):
    """Overwrite the cached mtime of `file_name` with the current on-disk
    value, or `None` if the path does not exist.

    Mirrors `updateModificationTime()` at SourceDigest.pm:295-299.
    """
    try:
        _modification_times[file_name] = os.stat(file_name).st_mtime
    except FileNotFoundError:
        _modification_times[file_name] = None


def hash_data_files(hasher, files):
    """Feed the bytes of each data file (found under `$GALACTICUS_DATA_PATH`)
    into `hasher`, skipping any that don't exist.

    Mirrors `Hash_Data_Files()` at SourceDigest.pm:272-285.
    """
    data_root = os.environ.get('GALACTICUS_DATA_PATH', '')
    for leaf in files:
        data_file = os.path.join(data_root, leaf)
        if os.path.exists(data_file):
            with open(data_file, 'rb') as fh:
                hasher.update(fh.read())


def find_hash(file_names, *, use_locks=True, include_files_excluded=None,
              report=False):
    """Compute a composite MD5 hash over `file_names` and all their Fortran
    dependencies, using (and updating) persistent per-file `.md5` / `.md5c`
    caches under `$BUILDPATH`.

    Mirrors `Find_Hash()` at SourceDigest.pm:68-270.

    Parameters
    ----------
    file_names : list of str
        The roots whose composite hash is sought (paths relative to the
        top-level build directory, e.g. `"source/foo.F90"`).
    use_locks : bool, default True
        If True, take an exclusive flock on each `<file>.md5*.lock` sidecar
        to serialise concurrent builds.
    include_files_excluded : list of str or None
        Names of `include` files that `read_file` should NOT follow.
    report : bool, default False
        If True, print verbose progress to stdout (mirrors the Perl `report`
        flag).

    Returns
    -------
    str
        The 22-character base64 MD5 digest over the union of composite
        digests of the supplied `file_names`.
    """
    if include_files_excluded is None:
        include_files_excluded = []

    build_path = os.environ['BUILDPATH']

    if report:
        print("=> Begin computing MD5 hash")

    hasher = hashlib.md5()

    for file_name in file_names:
        if report:
            print(f" => Process file: {file_name}")

        if file_name in _composite_digests:
            hasher.update(_composite_digests[file_name].encode('ascii'))
            if report:
                print(
                    f"  => Use pre-existing composite hash: "
                    f"{_composite_digests[file_name]}"
                )
            continue

        # Compute / load composite digest for this file.
        composite_hasher     = hashlib.md5()
        hash_file_name       = re.sub(
            r'\.F90$', '.md5c',
            os.path.join(build_path, file_name),
        )
        dependency_file_name = re.sub(
            r'\.F90$', '.d',
            os.path.join(build_path, file_name),
        )

        if not os.path.exists(dependency_file_name):
            continue

        if report:
            print("  => Processing dependencies")

        composite_lock = open(hash_file_name + '.lock', 'w')
        try:
            if use_locks:
                fcntl.flock(composite_lock, fcntl.LOCK_EX)

            use_stored_composite = os.path.exists(hash_file_name)

            # First pass: decide whether the stored composite is still
            # current against every source file it covers.
            if use_stored_composite:
                with open(dependency_file_name) as dep_fh:
                    for object_file_name in dep_fh:
                        object_file_name = object_file_name.rstrip('\n')
                        source_prefix = re.sub(
                            r'^' + re.escape(build_path) + r'/(.*)\.o$',
                            r'source/\1',
                            object_file_name,
                        )
                        if report:
                            print(f"   => Dependency file: {source_prefix}")
                        for suffix in ('F90', 'c', 'h', 'Inc', 'cpp'):
                            source_file_name = f"{source_prefix}.{suffix}"
                            if not os.path.exists(source_file_name):
                                continue
                            if report:
                                print(f"    => Dependency file: {source_file_name}")
                            # Mirror Perl: `$md5FileName = $sourceFileName; $md5FileName =~ s/^source//;
                            #              $md5FileName = $BUILDPATH.$md5FileName.".md5";`
                            stripped = re.sub(
                                r'^source', '', source_file_name,
                            )
                            md5_file_name = build_path + stripped + '.md5'
                            if source_file_name not in _digests:
                                if not (
                                    os.path.exists(md5_file_name)
                                    and modification_time(hash_file_name)
                                        > modification_time(md5_file_name)
                                    and modification_time(md5_file_name)
                                        > modification_time(source_file_name)
                                ):
                                    use_stored_composite = False
                                if report:
                                    if use_stored_composite:
                                        print(f"     => Can use stored hash: {source_file_name}")
                                    else:
                                        print(f"     => Can not use stored hash: {source_file_name}")
                            else:
                                if not (
                                    modification_time(hash_file_name)
                                    > modification_time(md5_file_name)
                                ):
                                    use_stored_composite = False
                                if report:
                                    if use_stored_composite:
                                        print("     => Digest pre-exists")
                                    else:
                                        print("     => Digest pre-exists but is outdated")

            # If the composite is still current, read and use it.
            if use_stored_composite:
                with open(hash_file_name) as md5_fh:
                    composite_digest = md5_fh.readline()
                if composite_digest:
                    _composite_digests[file_name] = composite_digest
                    hasher.update(composite_digest.encode('ascii'))
                    if report:
                        print(
                            f"   => Reading stored composite hash:\t"
                            f"{file_name}\t{hash_file_name}\t"
                            f"{_composite_digests[file_name]}"
                        )
                else:
                    use_stored_composite = False

            # Otherwise, recompute the composite from scratch.
            if not use_stored_composite:
                if report:
                    print("   => Computing composite hash")
                with open(dependency_file_name) as dep_fh:
                    for object_file_name in dep_fh:
                        object_file_name = object_file_name.rstrip('\n')
                        source_prefix = re.sub(
                            r'^' + re.escape(build_path) + r'/(.*)\.o$',
                            r'source/\1',
                            object_file_name,
                        )
                        for suffix in ('F90', 'c', 'h', 'Inc', 'cpp'):
                            source_file_name = f"{source_prefix}.{suffix}"
                            if not os.path.exists(source_file_name):
                                continue
                            if source_file_name not in _digests:
                                stripped = re.sub(
                                    r'^source', '', source_file_name,
                                )
                                md5_file_name = build_path + stripped + '.md5'
                                md5_lock = open(md5_file_name + '.lock', 'w')
                                try:
                                    if use_locks:
                                        fcntl.flock(md5_lock, fcntl.LOCK_EX)
                                    use_stored = False
                                    if (os.path.exists(md5_file_name)
                                            and modification_time(md5_file_name)
                                                > modification_time(source_file_name)):
                                        use_stored = True
                                        if suffix in ('F90', 'Inc'):
                                            for match in get_matching_lines(
                                                source_file_name,
                                                re.compile(
                                                    r'[\"\'](data/[a-zA-Z0-9_\.\-/]+\.(xml|hdf5))[\"\']'
                                                ),
                                            ):
                                                data_file_name = match['submatches'][0]
                                                if not (
                                                    os.path.exists(md5_file_name)
                                                    and modification_time(md5_file_name)
                                                        > modification_time(data_file_name)
                                                ):
                                                    use_stored = False
                                    if use_stored:
                                        with open(md5_file_name) as md5_fh:
                                            digest = md5_fh.readline()
                                        if digest:
                                            _digests[source_file_name] = digest
                                            if report:
                                                print(
                                                    f"   => Reading stored hash: "
                                                    f"{source_file_name}\t"
                                                    f"{_digests[source_file_name]}"
                                                )
                                        else:
                                            use_stored = False
                                    if not use_stored:
                                        if report:
                                            print("   => Computing hash")
                                        file_hasher = hashlib.md5()
                                        if suffix in ('F90', 'Inc'):
                                            text = read_file(
                                                source_file_name,
                                                state='raw',
                                                follow_includes=True,
                                                include_locations=['../source', '../' + build_path],
                                                include_files_excluded=include_files_excluded,
                                                strip_regex=re.compile(r'^\s*!(?!(!\[|\$)).*$'),
                                                strip_leading=True,
                                                strip_trailing=True,
                                                strip_empty=True,
                                            )
                                            file_hasher.update(text.encode('utf-8', errors='replace'))
                                            extra_files = [
                                                m['submatches'][1]
                                                for m in get_matching_lines(
                                                    source_file_name,
                                                    re.compile(
                                                        r'(char\s*\()??\s*galacticusPath\s*'
                                                        r'\(\s*pathTypeDataStatic\s*\)\s*\)??\/\/\]\s*'
                                                        r'[\"\']([a-zA-Z0-9_\.\-/]+\.(xml|hdf5))[\"\']'
                                                    ),
                                                )
                                            ]
                                            hash_data_files(file_hasher, extra_files)
                                            if report:
                                                print(
                                                    f"    => Computed hash from: "
                                                    f"{source_file_name} {{{', '.join(extra_files)}}}"
                                                )
                                        else:
                                            with open(source_file_name, 'rb') as raw_fh:
                                                file_hasher.update(raw_fh.read())
                                            if report:
                                                print(
                                                    f"    => Computed hash from: "
                                                    f"{source_file_name} {{RAW}}"
                                                )
                                        _digests[source_file_name] = _b64digest(file_hasher)
                                        with open(md5_file_name, 'w') as md5_fh:
                                            md5_fh.write(_digests[source_file_name])
                                        update_modification_time(md5_file_name)
                                        if report:
                                            print(
                                                f"   => Stored hash: {source_file_name}\t"
                                                f"{_digests[source_file_name]}"
                                            )
                                finally:
                                    md5_lock.close()
                            if source_file_name in _digests:
                                composite_hasher.update(
                                    _digests[source_file_name].encode('ascii')
                                )
                            else:
                                raise RuntimeError(
                                    "find_hash: failed to build digest for "
                                    f"'{source_file_name}'"
                                )
                _composite_digests[file_name] = _b64digest(composite_hasher)
                hasher.update(_composite_digests[file_name].encode('ascii'))
                with open(hash_file_name, 'w') as md5_fh:
                    md5_fh.write(_composite_digests[file_name])
                update_modification_time(file_name)
                if report:
                    print(
                        f"   => Composite hash stored\t{file_name}\t"
                        f"{hash_file_name}\t{_composite_digests[file_name]}"
                    )
        finally:
            composite_lock.close()

    final_hash = _b64digest(hasher)
    if report:
        print(f"=> MD5 hash: {final_hash}")
    return final_hash
