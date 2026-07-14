#!/usr/bin/env python3
"""Compute MD5 source digests for every `functionClass` / `sourceDigest`
type that participates in a given executable, and emit a C source file
binding each type's composite digest to a `<type>MD5` character array.

Walks every source file whose `.o` appears in `<target>.d`, extracts
`<sourceDigest>` / `<functionClass>` directives, and parses every
functionClass instance file's source tree to collect the derived-type
dependency graph.  For each type, `Find_Hash` (ported Phase-0) computes
an MD5 over the file and every source it `include`s, skipping a small
set of non-deterministic inputs (`fftw3.f03`, build-environment / version
includes).  Composite digests are then resolved by iteratively combining
each type's source digest with those of its declared dependencies.

Outputs:
  * `$BUILDPATH/<target>.md5s.c`  -- one `char <type>MD5[] = "..."` per type
                                     (slashes in b64 digests replaced by `@`
                                     so the identifier is valid in C).
  * `$BUILDPATH/<target>.md5.blob` -- pickle cache of per-file digest data.

Mirrors scripts/build/sourceDigests.pl.
Andrew Benson (ported to Python 2026).
"""

import base64
import hashlib
import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.Directives                             import extract_directives
from Galacticus.Build.ParallelScan                           import scan as parallel_scan, resolve_jobs
from Galacticus.Build.SourceTree                             import parse_file, walk_tree
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import class_dependencies
from Galacticus.Build.SourceTree.Process.SourceDigest        import find_hash, ensure_file_digest
from Galacticus._logging                                     import configure_default as _configure_default
from List.ExtraUtils                                         import as_array
from XML.Utils                                               import xml_to_dict
from Galacticus.Build.ScanCache import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
)

# Show INFO-level diagnostic output from the library modules (mirrors the
# verbose `print()`-driven output of the Perl-era driver).
_configure_default()


_EXCLUDED_INCLUDES = [
    'fftw3.f03',
    'output.build.environment.inc',
    'output.version.revision.inc',
]

# Version stamp for the `.md5.blob` cache schema. Bump when the per-file entry
# structure changes; caches with a different (or absent) version are discarded.
_BLOB_VERSION = 2


def _load_xml(path):
    if not os.path.exists(path):
        return None
    return xml_to_dict(ET.parse(path).getroot())


def _object_files_from_dep(dependency_file_name, build_path):
    """Extract bare `<name>.o` entries from `<target>.d`."""
    if not os.path.exists(dependency_file_name):
        return []
    prefix_re = re.compile(r'^' + re.escape(build_path) + r'/(.+\.o)$')
    out = []
    with open(dependency_file_name, 'r') as fh:
        for line in fh:
            m = prefix_re.match(line.rstrip('\n'))
            if m:
                out.append(m.group(1))
    return out


def _functionclass_names_and_modules(state_storables):
    """Return a list of functionClass *directive* names (i.e. the
    `<functionClass>Class` entries from `stateStorables.xml`, with the
    trailing `Class` stripped) plus the `functionClassInstances` names.

    Mirrors sourceDigests.pl:42-47.
    """
    fcs = (state_storables or {}).get('functionClasses') or {}
    if isinstance(fcs, dict):
        fc_keys = sorted(fcs.keys())
    elif isinstance(fcs, list):
        fc_keys = sorted(
            d['name'] for d in fcs if isinstance(d, dict) and 'name' in d
        )
    else:
        fc_keys = []
    function_class_names = [re.sub(r'Class$', '', k) for k in fc_keys]

    instances = (state_storables or {}).get('functionClassInstances') or []
    if isinstance(instances, str):
        instances = [instances]
    return function_class_names, list(instances)


def _append_dep_files(entry_files, file_name, build_path, source_root):
    """Read the `.d` sidecar for the source file `file_name` (a path relative
    to `source_root`) and record every source file
    (`.F90`/`.Inc`/`.cpp`/`.c`/`.h`) its object entries reference.

    The source tree is hierarchical, so both the `.d` sidecar
    (`$BUILDPATH/<subdir>/<name>.d`) and the object entries it lists
    (`$BUILDPATH/<subdir>/<name>.o`) carry sub-directory components that must
    be preserved when mapping back to source files.
    """
    if not file_name.endswith('.F90'):
        return
    dep_file = os.path.join(build_path, file_name[:-len('.F90')] + '.d')
    for obj in _object_files_from_dep(dep_file, build_path):
        root = os.path.join(source_root, obj[:-len('.o')]) + '.'
        for suffix in ('F90', 'Inc', 'cpp', 'c', 'h'):
            candidate = root + suffix
            if os.path.exists(candidate):
                entry_files.append(candidate)
                break


def _b64digest_no_pad(hasher):
    return base64.b64encode(hasher.digest()).rstrip(b'=').decode('ascii')


# ---------------------------------------------------------------------------
# Parallel pre-warm of the per-file `.md5` sidecars.
# ---------------------------------------------------------------------------
# The per-source-file digests are shared across every type's composite and
# dominate a cold run (all of it is single-threaded MD5 hashing of file
# contents and their includes). They are pre-computed in parallel below, one
# task per unique source file, BEFORE the serial per-type loop calls
# `find_hash` (which then just reads the warm sidecars). Read-only inputs are
# published to `_WORKER` so forked workers inherit them via copy-on-write.
_WORKER = {}


def _prewarm_worker(task):
    """Compute (and write the `.md5` sidecar for) one source file. The return
    value is unused — the parent re-reads each sidecar through `find_hash`."""
    source_file_name, suffix = task
    ensure_file_digest(
        source_file_name, suffix, _WORKER['build_path'],
        use_locks=_WORKER['use_locks'],
        include_files_excluded=_EXCLUDED_INCLUDES,
    )
    return None


def _prewarm_file_digests(object_files, build_path, use_locks):
    """Pre-compute, in parallel, the per-file `.md5` digest of every source in
    the target's dependency set.

    Each source file is a single task with a single writer, so the parallel
    pass needs no intra-run lock (the `flock` inside `ensure_file_digest`,
    taken when `use_locks`, still serialises concurrent *external* builds).
    This is a pure warm-up: any file it skips — including one that is stale
    only via a data-file dependency the cheap filter below ignores — is still
    computed, correctly and serially, by `find_hash`.
    """
    tasks = []
    seen  = set()
    for obj in sorted(object_files):
        source_prefix = 'source/' + re.sub(r'\.o$', '', obj)
        for suffix in ('F90', 'c', 'h', 'Inc', 'cpp'):
            source_file_name = source_prefix + '.' + suffix
            if source_file_name in seen or not os.path.exists(source_file_name):
                continue
            seen.add(source_file_name)
            # Cheap staleness pre-filter so a warm run does not fork a worker
            # pool for no-op tasks: skip files whose `.md5` sidecar is already
            # newer than the source. Over-selection is harmless
            # (`ensure_file_digest` re-checks precisely and returns fast).
            md5_file_name = (build_path
                             + re.sub(r'^source', '', source_file_name)
                             + '.md5')
            try:
                if os.stat(md5_file_name).st_mtime \
                        > os.stat(source_file_name).st_mtime:
                    continue
            except FileNotFoundError:
                pass
            tasks.append((source_file_name, suffix))

    if not tasks:
        return
    _WORKER['build_path'] = build_path
    _WORKER['use_locks']  = use_locks
    parallel_scan(tasks, _prewarm_worker, 'sourceDigests.py')


def _scan_file_worker(task):
    """Process one source file: collect its tracked files and the per-type
    source-MD5 records for its sourceDigest / functionClass / functionClass-
    instance directives.  Returns `(file_identifier, entry, type_records)`
    where `type_records` is an ordered list of `(hash_name, record)`.

    Runs in a forked worker (see `_WORKER`). The `find_hash` calls only READ
    the per-file `.md5` sidecars pre-warmed by `_prewarm_file_digests`, and
    each worker writes only its own file's `.md5c` composite — so distinct
    workers never write the same sidecar.
    """
    file_name, file_to_process, file_identifier = task
    build_path              = _WORKER['build_path']
    source_root             = _WORKER['source_root']
    use_locks               = _WORKER['use_locks']
    sd_files                = _WORKER['sd_files']
    fc_files                = _WORKER['fc_files']
    function_class_file_set = _WORKER['function_class_file_set']
    function_class_names    = _WORKER['function_class_names']

    entry = {'files': [file_to_process], 'typesProvided': []}
    _append_dep_files(entry['files'], file_name, build_path, source_root)
    type_records = []

    def _source_md5():
        return find_hash([file_name],
                         include_files_excluded=_EXCLUDED_INCLUDES,
                         use_locks=use_locks)

    # sourceDigest directives -> per-type source MD5 (no dependencies).
    if file_to_process in sd_files:
        for sd in extract_directives(file_to_process, 'sourceDigest'):
            hash_name = sd['name']
            type_records.append((hash_name,
                                 {'sourceMD5': _source_md5(),
                                  'dependencies': []}))
            entry['typesProvided'].append(hash_name)

    # functionClass directives -> per-type source MD5 + explicit base.
    if file_to_process in fc_files:
        for fc in extract_directives(file_to_process, 'functionClass'):
            hash_name = fc['name'] + 'Class'
            type_records.append((hash_name,
                                 {'sourceMD5': _source_md5(),
                                  'dependencies':
                                      [fc['extends']] if 'extends' in fc
                                      else ['functionClass']}))
            entry['typesProvided'].append(hash_name)

    # functionClass instance files: walk the AST and record a per-derived-
    # class hash + its inheritance dependencies.
    if file_to_process in function_class_file_set:
        tree = parse_file(file_to_process)
        for node in walk_tree(tree):
            ntype = node.get('type')
            if ntype not in function_class_names:
                continue
            hash_name = (node.get('directive') or {}).get('name')
            if not hash_name:
                continue
            class_record, deps = class_dependencies(node, ntype)
            class_type = (class_record or {}).get('type')
            deps = [d for d in deps if d != class_type]
            type_records.append((hash_name,
                                 {'sourceMD5': _source_md5(),
                                  'dependencies': list(deps)}))
            entry['typesProvided'].append(hash_name)

    return file_identifier, entry, type_records


def main(argv):
    if len(argv) != 4:
        print("Usage: sourceDigests.py <sourceDirectory> <target> <useLocks>",
              file=sys.stderr)
        sys.exit(1)
    source_directory_name = argv[1]
    target_name           = argv[2]
    use_locks_arg         = argv[3]
    use_locks             = use_locks_arg == 'yes' or use_locks_arg == '1'
    build_path            = os.environ['BUILDPATH']

    state_storables = _load_xml(
        os.path.join(build_path, 'stateStorables.xml'),
    )
    locations_path = os.path.join(build_path, 'directiveLocations.xml')
    if not os.path.exists(locations_path):
        sys.exit("Error: directiveLocations.xml not found")
    locations = xml_to_dict(ET.parse(locations_path).getroot())

    function_class_names, instance_names = _functionclass_names_and_modules(
        state_storables,
    )

    function_class_files = []
    for name in function_class_names:
        function_class_files.extend(
            as_array((locations.get(name) or {}).get('file'))
        )
    function_class_file_set = set(function_class_files)

    allowed_names = ['functionClass'] + function_class_names + instance_names

    dep_name = re.sub(r'\.(exe|o)$', '.d', target_name)
    blob_name = re.sub(r'\.(exe|o)$', '.md5.blob', target_name)
    out_name = re.sub(r'\.(exe|o)$', '.md5s.c', target_name)
    dependency_file_name = os.path.join(build_path, dep_name)
    blob_path            = os.path.join(build_path, blob_name)
    output_path          = os.path.join(build_path, out_name)

    object_files = set(_object_files_from_dep(dependency_file_name,
                                              build_path))

    digests_per_file, cache_mtime = _load_cache(blob_path)
    # Discard caches written by earlier code versions: entries now record the
    # types they provide (`typesProvided`), which the stale-type pruning below
    # relies on. A version mismatch costs one full rescan of this target.
    if digests_per_file.get('blobVersion') != _BLOB_VERSION:
        digests_per_file, cache_mtime = {}, None
    digests_per_file['blobVersion'] = _BLOB_VERSION
    digests_per_file.setdefault('types', {})

    updated_types     = {}
    seen_identifiers  = set()

    source_root = os.path.join(source_directory_name, 'source')
    if not os.path.isdir(source_root):
        sys.exit("sourceDigests.py: can not open the source directory: "
                 + source_root)
    source_names = []
    for dirpath, dirnames, filenames in os.walk(source_root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for f in filenames:
            source_names.append(
                os.path.relpath(os.path.join(dirpath, f), source_root))
    source_names.sort()

    # Warm the per-file `.md5` sidecars in parallel before the serial per-type
    # loop below reads them through `find_hash` (see `_prewarm_file_digests`).
    _prewarm_file_digests(object_files, build_path, use_locks)

    # Select the files to (re)scan, preserving the deterministic source_names
    # order so the parallel merge below reproduces the serial run's
    # last-writer-wins semantics for any shared type name. `seen_identifiers`
    # tracks every participating file (for stale-entry pruning), independent
    # of whether it needs a rescan.
    files_to_scan = []
    for file_name in source_names:
        if os.path.basename(file_name).startswith('.#'):
            continue
        if not re.search(r'\.(F90|cpp)$', file_name):
            continue
        object_file_name = re.sub(r'\.(F90|cpp)$', '.o', file_name)
        if object_file_name not in object_files:
            continue

        file_to_process = os.path.join(source_root, file_name)
        file_identifier = _file_identifier(file_to_process)
        seen_identifiers.add(file_identifier)

        rescan = True
        if file_identifier in digests_per_file:
            tracked = digests_per_file[file_identifier].get('files') or []
            stale = any(
                os.path.exists(t) and os.stat(t).st_mtime > cache_mtime
                for t in tracked
            )
            rescan = bool(stale)
        if not rescan:
            continue
        files_to_scan.append((file_name, file_to_process, file_identifier))

    # Compute the per-type source digests concurrently (one task per file).
    # `find_hash` inside each worker only reads the pre-warmed `.md5` sidecars
    # and writes its own file's `.md5c`, so distinct workers never contend on
    # a sidecar; locks are forced on when the scan actually runs in parallel
    # to guard the rare dependency `.md5` the pre-warm pass did not cover
    # (and any concurrent external build).
    _WORKER['build_path']              = build_path
    _WORKER['source_root']             = source_root
    _WORKER['use_locks']               = use_locks or resolve_jobs(len(files_to_scan)) > 1
    _WORKER['sd_files']                = set(
        as_array((locations.get('sourceDigest') or {}).get('file')))
    _WORKER['fc_files']                = set(
        as_array((locations.get('functionClass') or {}).get('file')))
    _WORKER['function_class_file_set'] = function_class_file_set
    _WORKER['function_class_names']    = function_class_names

    for file_identifier, entry, type_records in parallel_scan(
            files_to_scan, _scan_file_worker, 'sourceDigests.py'):
        digests_per_file.pop(file_identifier, None)
        digests_per_file[file_identifier] = entry
        for hash_name, record in type_records:
            digests_per_file['types'][hash_name] = record
            updated_types[hash_name] = 1

    # Manually add the base `functionClass` type (no dependencies).
    fc_base = 'objects/function_class.F90'
    fc_base_path = os.path.join(source_root, fc_base)
    if (cache_mtime is None or
            (os.path.exists(fc_base_path)
             and os.stat(fc_base_path).st_mtime > cache_mtime)):
        digests_per_file['types']['functionClass'] = {
            'sourceMD5': find_hash(
                [fc_base],
                include_files_excluded=_EXCLUDED_INCLUDES,
                use_locks=use_locks,
            ),
            'dependencies': [],
        }
        digests_per_file['types']['functionClass'].pop('compositeMD5', None)
        updated_types['functionClass'] = 1

    # Prune stale cache entries. Per-file entries for source files that no
    # longer participate in this target (deleted, renamed, or dropped from the
    # `.d` list) would otherwise persist forever, and any types they provided
    # would keep contributing (stale) digests. Pruned types are queued for
    # composite invalidation so that types depending on them are re-resolved.
    _RESERVED_KEYS = {'types', 'blobVersion'}
    for stale_id in [k for k in digests_per_file
                     if k not in _RESERVED_KEYS and k not in seen_identifiers]:
        del digests_per_file[stale_id]
    valid_types = {'functionClass'}
    for key, entry in digests_per_file.items():
        if key in _RESERVED_KEYS:
            continue
        valid_types.update(entry.get('typesProvided') or [])
    for stale_type in [t for t in digests_per_file['types']
                       if t not in valid_types]:
        del digests_per_file['types'][stale_type]
        updated_types[stale_type] = 1

    # Invert the dependency graph so we can transitively clear composite
    # hashes for every type depending on an updated one.
    dependencies_inverted = {}
    for hash_name, info in digests_per_file['types'].items():
        for dep in info.get('dependencies') or []:
            dependencies_inverted.setdefault(dep, []).append(hash_name)

    # Breadth-first invalidation of composite hashes.
    queue = dict(updated_types)
    while queue:
        reset_name, _ = queue.popitem()
        for dependent in dependencies_inverted.get(reset_name, []):
            queue[dependent] = 1
        info = digests_per_file['types'].get(reset_name)
        if info is not None:
            info.pop('compositeMD5', None)

    # Iteratively resolve composite hashes.
    allowed_set = set(allowed_names)
    resolved = False
    while not resolved:
        resolved = True
        updated  = False
        for hash_name in sorted(digests_per_file['types']):
            info = digests_per_file['types'][hash_name]
            if 'compositeMD5' in info:
                continue
            deps = info.get('dependencies') or []
            # Blocked if any dependency exists in the types table but
            # has no compositeMD5 yet.
            blocked = any(
                dep in digests_per_file['types']
                and 'compositeMD5' not in digests_per_file['types'][dep]
                for dep in deps
            )
            if blocked:
                resolved = False
                continue

            updated = True
            hasher  = hashlib.md5()
            hasher.update(info['sourceMD5'].encode('ascii'))
            for dep in deps:
                if dep not in allowed_set:
                    continue
                dep_info = digests_per_file['types'].get(dep)
                if dep_info is None:
                    continue
                hasher.update(dep_info['compositeMD5'].encode('ascii'))
            info['compositeMD5'] = _b64digest_no_pad(hasher)
        if not resolved and not updated:
            sys.exit("sourceDigest.py: failed to resolve dependencies")

    with open(blob_path, 'wb') as fh:
        pickle.dump(digests_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)

    with open(output_path, 'w') as fh:
        for hash_name in sorted(digests_per_file['types']):
            info = digests_per_file['types'][hash_name]
            digest = info.get('compositeMD5', '')
            # Replace `/` with `@` so the b64 digest is a valid C string
            # literal in a header-less context (mirrors Perl line 213).
            digest = digest.replace('/', '@')
            fh.write(f'char {hash_name}MD5[]="{digest}";\n')


if __name__ == '__main__':
    main(sys.argv)
