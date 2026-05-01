# Replace a file only when its content has actually changed.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/File/Changes.pm -- the `Update(oldFile, newFile, proveUpdate =>)`
# helper that compares `old_file` with the freshly-written `new_file`, and
# either moves `new_file` over `old_file` (if different) or deletes
# `new_file` (if identical).  Preserving the old file's mtime in the
# identical case lets Make skip rebuilds of downstream targets.

import filecmp
import os
import shutil

__all__ = ['update']


def update(old_file, new_file, *, prove_update=False):
    """Atomically update `old_file` from `new_file`, leaving mtime untouched
    when nothing changed.

    * If `old_file` does not exist, `new_file` is moved into its place.
    * If the two files are byte-for-byte identical, `new_file` is deleted
      and `old_file` is left alone (mtime preserved).
    * Otherwise, `new_file` replaces `old_file`.

    When `prove_update=True`, an `<old_file>.up` sentinel is touched to
    record that work was done -- matching Perl `proveUpdate => 'yes'`.
    """
    if not os.path.exists(old_file):
        shutil.move(new_file, old_file)
    elif filecmp.cmp(old_file, new_file, shallow=False):
        os.unlink(new_file)
    else:
        shutil.move(new_file, old_file)

    if prove_update:
        sentinel = old_file + '.up'
        with open(sentinel, 'a'):
            os.utime(sentinel, None)
