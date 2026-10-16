"""Regression tests for ``scripts/doc/extractContributors.py``.

Guards both halves of the contributor-marker contract:

* extraction — the names parsed out of ``!+`` markers, which become the
  Acknowledgments list in the documentation; and
* validation (``--check``) — the enforcement that a marker records *names
  only*.  Markers which described a contribution rather than naming its
  contributors used to leak fragments of prose ("assistance from Claude",
  "per-object OpenMP lock", "and reviewed") into that published list.
"""
import os
import sys

# extractDocsRST self-inserts its own directory on sys.path for its sibling
# imports; do the same so we can import this module directly.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
import extractContributors  # noqa: E402


CANONICAL = ('!+    Contributions to this file made by: '
             'Andrew Benson, Claude.\n')


# --- Extraction ------------------------------------------------------------

def test_extracts_names_from_canonical_marker():
    names = extractContributors.extract_contributors(
        CANONICAL, extractContributors._FORTRAN_MARKER)
    assert names == ['Andrew Benson', 'Claude']


def test_extraction_ignores_non_markers():
    assert extractContributors.extract_contributors(
        '  ! An ordinary comment.\n',
        extractContributors._FORTRAN_MARKER) == []


def test_render_contributors_rst():
    assert extractContributors.render_contributors_rst(
        ['Andrew Benson', 'Claude']) == 'Andrew Benson, Claude.\n'


# --- Name validation -------------------------------------------------------

def test_accepts_real_names():
    for name in ('Andrew Benson', 'Claude', 'Copilot', 'Matías Liempi',
                 'Stéphane Mangeon', 'Luiz Felippe S. Rodrigues',
                 'Ludwig van Beethoven'):
        assert extractContributors.name_problem(name) is None, name


def test_rejects_prose():
    for name in ('were drafted with assistance from Claude',
                 'per-object OpenMP lock',
                 'and reviewed and verified by',
                 'for issue #1353',
                 'Implemented with assistance from Claude (Anthropic)',
                 # Capitalized prose still ends a sentence first.
                 'Andrew Benson. The Fix For Issue'):
        assert extractContributors.name_problem(name) is not None, name


# --- Marker validation -----------------------------------------------------

def test_canonical_marker_passes():
    assert extractContributors.check_markers_in_text(CANONICAL) == []


def test_indented_marker_passes():
    assert extractContributors.check_markers_in_text(
        '  !+ Contributions to this file made by: Andrew Benson, Claude.\n') == []


def test_marker_wrapped_over_two_lines_passes():
    text = ('!+    Contributions to this file made by: Arya Farahi, '
            'Andrew Benson,\n'
            '!+    Christoph Behrens, Xiaolong Du, Claude.\n')
    assert extractContributors.check_markers_in_text(text) == []


def test_marker_without_the_canonical_label_is_reported():
    problems = extractContributors.check_markers_in_text(
        '!+ Implemented by Andrew Benson with assistance from Claude.\n')
    assert any('does not begin' in problem for _, _, problem in problems)


def test_marker_describing_the_contribution_is_reported():
    text = ('!+    Contributions to this file made by: Andrew Benson. The fix '
            'for issue #1321 was diagnosed and drafted with\n'
            '!+    assistance from Claude, and reviewed and verified by '
            'Andrew Benson.\n')
    problems = extractContributors.check_markers_in_text(text)
    assert problems, 'prose marker was not reported'
    assert all(path == '<text>' for path, _, _ in problems)


def test_description_appended_after_the_names_is_reported():
    problems = extractContributors.check_markers_in_text(
        '!+    Contributions to this file made by: Andrew Benson, Claude. '
        'Drafted for issue #1317.\n')
    assert problems


def test_ordinary_comments_are_not_markers():
    assert extractContributors.check_markers_in_text(
        '  ! Drafted with assistance from Claude, and verified.\n') == []


def test_problems_carry_the_path_and_line_number():
    text = ('module Test\n'
            '!+ This file was generated with Claude\n')
    problems = extractContributors.check_markers_in_text(text, 'test.F90')
    assert problems
    for path, line_number, _ in problems:
        assert path == 'test.F90'
        assert line_number == 2


# --- The source tree itself ------------------------------------------------

def test_source_tree_markers_are_well_formed():
    """The Galacticus source must itself satisfy the check enforced in CI."""
    root = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir,
                        os.pardir)
    problems = extractContributors.check_markers(root)
    assert problems == [], '\n'.join(f'{path}:{line}: {problem}'
                                     for path, line, problem in problems)
