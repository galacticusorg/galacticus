"""Regression test for `Galacticus.Build.Dependencies.dependency_sort`.

Bug: `dependency_sort` mirrored Perl `Dependency_Sort`'s edge-building
pattern verbatim — for `X.after = Y` it pushed `dependencies[Y] += [X]`.
But Perl's `Sort::Topo` had a docstring/algorithm mismatch: the
docstring said `dependencies[X] = [Y]` meant "X depends on Y" (Y first),
while the actual algorithm emitted X *first*.  The two opposites
cancelled out so Perl callers got the right answer.

Our Python `Sort.Topo` wraps `graphlib.TopologicalSorter`, which
matches the *standard* "X depends on Y → Y first" convention with no
such inversion.  Copying the Perl edges verbatim therefore produced
results inverted from natural English: `<eventHookStatic>` directives
tagged `after="other"` ended up BEFORE `other` in the emitted call
list.  In Galacticus that surfaced as e.g.
`Stellar_Luminosities_Initializor` being called before
`Stellar_Luminosities_Initialize` even though it had
`after="Stellar_Luminosities_Initialize"` — and the resulting code
crashed at runtime.

Fix: build the edges the way the standard reading dictates —
`X.after = Y` → `deps[X] += [Y]` (X depends on Y, Y emitted first).
"""

from Galacticus.Build.Dependencies import dependency_sort


def test_after_places_X_after_Y():
    """`X.after = Y` → Y comes BEFORE X in the output.

    Mirrors the user-reported real-world failure: `Initializor` had
    `after="Initialize"`, so `Initialize` should be emitted first.
    """
    out = dependency_sort({
        'Stellar_Luminosities_Initialize':  {},
        'Stellar_Luminosities_Initializor': {
            'after': 'Stellar_Luminosities_Initialize'},
    })
    assert out == [
        'Stellar_Luminosities_Initialize',
        'Stellar_Luminosities_Initializor',
    ]


def test_before_places_X_before_Y():
    """`X.before = Y` → X emitted first, Y emitted later."""
    assert dependency_sort({'A': {'before': 'B'}, 'B': {}}) == ['A', 'B']


def test_unknown_reference_silently_dropped():
    """A name in `after`/`before` that isn't itself a task is ignored —
    matches Perl's `if (grep {$_ eq $dependent} @objects)` skip."""
    assert dependency_sort({'A': {'after': 'doesNotExist'}}) == ['A']


def test_list_of_after_references():
    """`after` may be a list of names; all of them must precede X."""
    out = dependency_sort({
        'A': {},
        'B': {},
        'C': {'after': ['A', 'B']},
    })
    assert out.index('A') < out.index('C')
    assert out.index('B') < out.index('C')


def test_no_constraints_alphabetical_baseline():
    """When nothing declares `after`/`before`, the order is the
    alphabetical sort `topo_sort` falls back on."""
    assert dependency_sort({'b': {}, 'a': {}, 'c': {}}) == ['a', 'b', 'c']


def test_after_and_before_combined():
    """`A.after = B` and `A.before = C` together pin A between B and C."""
    out = dependency_sort({
        'A': {'after': 'B', 'before': 'C'},
        'B': {},
        'C': {},
    })
    assert out.index('B') < out.index('A') < out.index('C')
