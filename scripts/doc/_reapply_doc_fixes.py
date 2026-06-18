#!/usr/bin/env python3
"""Re-apply the manual docstring content fixes after a full re-conversion.

A full re-convert regenerates every ``source/*.F90`` from the pristine LaTeX, so
the hand edits committed on the branch (closing parens, a units rephrase, the
marsh citation case, the betaSatellite typo, …) are lost.  This re-applies them
by replaying the exact line replacements from their commits, asserting each old
text is present so any drift from converter changes surfaces immediately.

Run AFTER convertDocstringsToRST.py.  Idempotent-ish: it errors loudly rather
than silently skipping if an expected fragment is missing.
"""
import re
import subprocess
import sys

# Commits whose every source hunk is a manual content fix (no converter-driven
# change), applied oldest-first so dependent edits (e.g. the nBody high-pass
# line gains a paren, then "lowest"->"highest") replay in order.
PURE_FIX_COMMITS = ['9d05c2e00', 'cd9b29cec', 'b64d6569f']

# Mixed commits (also contain re-conversions): the specific manual edits only.
EXTRA = [
    # construct.read: closing paren after the descendantNode sentence.
    ('source/merger_trees.construct.read.F90',
     'descendant node (as specified by the ``descendantNode`` property. Where',
     'descendant node (as specified by the ``descendantNode`` property). Where'),
    # Behroozi table: the beta_sat row referenced the betaCut parameter.
    ('source/halo_model.conditional_mass_function.Behroozi2010.F90',
     '      * - :math:`\\beta_\\mathrm{sat}`\n'
     '        - 0.859\n'
     '        - ``[conditionalStellarMassFunctionBehrooziBetaCut]``',
     '      * - :math:`\\beta_\\mathrm{sat}`\n'
     '        - 0.859\n'
     '        - ``[conditionalStellarMassFunctionBehrooziBetaSatellite]``'),
]


def pairs_from_commit(commit):
    """Yield (path, old_line, new_line) from a commit's source diff, pairing
    each removed line with the immediately following added line."""
    diff = subprocess.check_output(
        ['git', 'show', commit, '--', 'source/*.F90']).decode('utf-8', 'replace')
    path = None
    lines = diff.split('\n')
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith('+++ b/'):
            path = ln[len('+++ b/'):]
        elif (ln.startswith('-') and not ln.startswith('---')
              and i + 1 < len(lines)
              and lines[i + 1].startswith('+')
              and not lines[i + 1].startswith('+++')):
            yield path, ln[1:], lines[i + 1][1:]
            i += 1
        i += 1


def main():
    fixes = []
    for commit in PURE_FIX_COMMITS:
        fixes.extend(pairs_from_commit(commit))
    fixes.extend(EXTRA)

    # A fix line captured from an older commit may contain ``\refClass`` output
    # in its previous form (inline ``code``) where the current converter now
    # emits a ``:galacticus-class:`` role.  When the literal text is not found,
    # retry with single-identifier literals rewritten to that role.
    def to_role(s):
        return re.sub(r'``([A-Za-z]\w*)``', r':galacticus-class:`\1`', s)

    applied = 0
    for path, old, new in fixes:
        with open(path, encoding='utf-8') as fh:
            text = fh.read()
        for o, n in ((old, new), (to_role(old), to_role(new))):
            if text.count(o) == 1:
                with open(path, 'w', encoding='utf-8') as fh:
                    fh.write(text.replace(o, n, 1))
                applied += 1
                break
        else:
            sys.exit(f'ERROR: {path}: fix target not found (literal or role '
                     f'form):\n  {old!r}')
    print(f're-applied {applied} content fixes across '
          f'{len({p for p, _, _ in fixes})} files')


if __name__ == '__main__':
    main()
