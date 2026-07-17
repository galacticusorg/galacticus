"""Shared helpers for the per-file pickle ("blob") caches used by the
`scripts/build/` catalog scripts.

Several build scripts (moduleDependencies, useDependencies, sourceDigests,
parameterDependencies, deepCopyActions, stateStorables, codeDirectivesParse,
buildCode) cache per-source-file scan results in a pickle blob under
``$BUILDPATH`` so that subsequent runs rescan only files whose mtimes have
advanced. They share a common protocol:

* Each file's cache entry is keyed by :func:`file_identifier` of its path.
* The blob's own mtime is the freshness reference: any file tracked by an
  entry whose mtime is newer than the blob forces a rescan of that entry.
* A missing, corrupt, or non-dict blob degrades gracefully to a full rescan
  (never fails the build) -- see :func:`load_cache`.
* Blobs are rewritten only-if-changed (via ``FileChanges.update``), which
  deliberately preserves the blob mtime when content is unchanged so the
  freshness window matches a serial run.

Historically each script carried its own copy of these helpers (a legacy of
the Perl-to-Python port); they live here so fixes land once.
"""

import os
import pickle
import re

__all__ = ['file_identifier', 'load_cache', 'prune']


def file_identifier(path):
    """Return the canonical blob-cache key for a source-file path.

    Slashes become underscores and a leading ``.`` (optionally followed by
    ``_``) is stripped -- mirrors the Perl
    ``(my $id = $path) =~ s/\\//_/g; $id =~ s/^\\._??//;`` idiom. The strip
    is a no-op for the absolute paths make passes; it matters only for
    hand-run invocations using ``./``-relative paths.
    """
    return re.sub(r'^\._?', '', path.replace('/', '_'))


def load_cache(blob_path):
    """Load a per-file pickle cache, returning ``(cache, mtime)``.

    Returns ``({}, None)`` when the blob is missing, unreadable, corrupt
    (e.g. a legacy Perl ``Storable`` blob), or not a dict -- forcing a full
    rescan rather than failing the build. Otherwise returns the cached dict
    together with the blob's mtime, which callers use as the freshness
    reference for per-entry staleness checks.
    """
    if not os.path.exists(blob_path):
        return {}, None
    try:
        with open(blob_path, 'rb') as fh:
            cache = pickle.load(fh)
    except (pickle.UnpicklingError, EOFError, AttributeError, ValueError,
            ImportError, ModuleNotFoundError):
        return {}, None
    if not isinstance(cache, dict):
        return {}, None
    return cache, os.stat(blob_path).st_mtime


def prune(cache, valid_identifiers, reserved=frozenset()):
    """Drop (in place) every cache entry whose key is not in
    ``valid_identifiers`` (nor in ``reserved``, for blobs that carry
    non-per-file bookkeeping keys such as ``blobVersion``). Returns the list
    of dropped keys.

    Callers must do this after computing their current file set and before
    any reduction over the whole cache: a rescan only ever *replaces* entries
    for files that still exist, so without pruning the entry for a deleted
    (or renamed, or directive-stripped) source file persists forever and its
    contributions are re-emitted into the generated catalogs and Makefile
    fragments — yielding phantom make dependencies on files that no longer
    exist, or stale generated code. Deleting the blob does not help: it is
    regenerated from the same stale cache.
    """
    stale = [fid for fid in cache
             if fid not in valid_identifiers and fid not in reserved]
    for fid in stale:
        del cache[fid]
    return stale
