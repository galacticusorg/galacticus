r"""Regression tests for ``scripts/doc/residualLatex.py``.

The migration of the embedded docstrings from LaTeX to RST left residue behind
*inside* regions that ``convertDocstringsToRST.py`` had already stamped as
converted.  Its ``--check`` mode cannot see that residue: the conversion is
idempotent, so it skips any region carrying ``!!{RST`` / ``docformat="rst"``,
and therefore only ever proves that a region was *visited*.  These tests pin
the complementary check — most importantly its two silences, since a checker
that cries wolf on legitimate RST (or on the LaTeX that is still correct by
design) would simply be switched off.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
import residualLatex  # noqa: E402


def kinds(text):
    return sorted(f.kind for f in residualLatex.scan_text(text, 'fixture.F90'))


# ---------------------------------------------------------------------------
# Silence on legitimate RST
# ---------------------------------------------------------------------------

def test_clean_rst_docstring_is_silent():
    """Roles, literals, explicit-markup blocks and RST's own backslash escapes
    all contain characters that look like LaTeX; none may be flagged."""
    text = (
        '  !!{RST\n'
        '  A rate of :math:`\\dot{M}_\\star` in ``[units]``, per\n'
        '  :cite:t:`smith_2020`.\n'
        '\n'
        '  .. math::\n'
        '\n'
        '     \\frac{\\mathrm{d}M}{\\mathrm{d}t} = -\\alpha M\n'
        '\n'
        '  A trailing ``literal``\\ s, and a `link <https://x/y#a>`_.\n'
        '  !!}\n'
    )
    assert kinds(text) == []


def test_openmp_sentinel_is_not_inline_math():
    """``!$omp`` quoted in prose is a Fortran sentinel, not a ``$…$`` pair."""
    text = ('  !!{RST\n'
            '  The read is made via !$omp atomic read, paired with the\n'
            '  !$omp atomic update on the counter.\n'
            '  !!}\n')
    assert kinds(text) == []


def test_code_block_contents_are_not_scanned():
    """A code-block may legitimately *show* LaTeX (or anything else)."""
    text = ('  !!{RST\n'
            '  For example:\n'
            '\n'
            '  .. code-block:: none\n'
            '\n'
            '     {\\normalfont \\ttfamily verbatim} $x$ \\emph{sample}\n'
            '\n'
            '  ...as shown above.\n'
            '  !!}\n')
    assert kinds(text) == []


# ---------------------------------------------------------------------------
# Silence on LaTeX that is still correct by design
# ---------------------------------------------------------------------------

def test_enumeration_entry_descriptions_are_exempt():
    r"""``<entry description="…">`` is run through ``latex_to_rst`` at build
    time by ``extractDocsRST``, so its ``\gls{…}`` must survive untouched."""
    text = ('  !![\n'
            '  <enumeration docformat="rst">\n'
            '   <entry label="a" description="In units of the \\gls{nsc} '
            'scale radius"/>\n'
            '  </enumeration>\n'
            '  !!]\n')
    assert kinds(text) == []


def test_constant_attributes_are_exempt():
    """``<constant>`` descriptions/units/references are converted at build
    time too, so their LaTeX math and citations are correct as-is."""
    text = ('  !![\n'
            '  <constant variable="x" docformat="rst" '
            'description="A ratio \\citep{hinkley_ratio_1969}" '
            'units="$\\mathrm{rad}/^\\circ$"/>\n'
            '  !!]\n')
    assert kinds(text) == []


def test_unconverted_block_is_left_to_the_conversion_check():
    """A block with no ``docformat="rst"`` is still LaTeX awaiting conversion;
    ``--check`` reports it, so flagging it here would double-report."""
    text = ('  !![\n'
            '  <description>Uses {\\normalfont \\ttfamily node}.</description>\n'
            '  !!]\n')
    assert kinds(text) == []


def test_unconverted_docstring_is_left_to_the_conversion_check():
    text = '  !!{\n  Uses {\\normalfont \\ttfamily node} and $x$.\n  !!}\n'
    assert kinds(text) == []


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def test_font_group_is_detected():
    text = '  !!{RST\n  Uses {\\normalfont \\ttfamily node} here.\n  !!}\n'
    assert 'font-group' in kinds(text)


def test_inline_math_reported_once_per_span():
    r"""``$\mathcal{O}(1)$`` is one defect, not one per macro inside it: the
    body is masked after the span is reported so the baseline stays legible."""
    text = '  !!{RST\n  An $\\mathcal{O}(1)$ correction factor.\n  !!}\n'
    assert kinds(text) == ['inline-math']


def test_unbalanced_dollar_is_still_reported():
    text = '  !!{RST\n  A stray $ sign.\n  !!}\n'
    assert kinds(text) == ['inline-math']


def test_environment_reported_without_duplicate_command():
    r"""``\begin{enumerate}`` is one finding, not also a ``\begin`` command."""
    text = ('  !!{RST\n'
            '  \\begin{enumerate}\n'
            '  \\item First.\n'
            '  \\end{enumerate}\n'
            '  !!}\n')
    assert kinds(text) == ['command:item',
                           'environment:enumerate', 'environment:enumerate']


def test_method_description_attribute_is_scanned():
    r"""``<method description="…">`` is the one attribute-form description
    that is RST rather than LaTeX."""
    text = ('  !![\n'
            '  <methods docformat="rst">\n'
            '   <method method="m" description="A value\\argin here"/>\n'
            '  </methods>\n'
            '  !!]\n')
    assert kinds(text) == ['command:argin']


def test_latex_interword_space_that_glues_words_is_detected():
    r"""``e.g.\ foo`` is LaTeX's inter-word space.  RST reads ``\ `` as an
    *escaped* space and deletes it, so this renders as "e.g.foo" — residue that
    corrupts prose silently instead of showing up as visible markup."""
    text = '  !!{RST\n  A halo scale (e.g.\\ virial radius) matters.\n  !!}\n'
    assert kinds(text) == ['control-space']


def test_interword_space_next_to_inline_markup_is_legitimate():
    r"""Next to inline markup, ``\ `` is the zero-width separator that lets
    docutils recognise the markup at all; flagging it would be wrong."""
    text = ('  !!{RST\n'
            '  In units of cm\\ :math:`^2` g\\ :math:`^{-1}`, plus a\n'
            '  ``literal``\\ s suffix, and a footnote\\ [#]_ marker.\n'
            '  !!}\n')
    assert kinds(text) == []


def test_latex_row_break_inside_math_is_not_a_control_space():
    r"""``\\`` inside display math is a LaTeX row break, not ``\ ``."""
    text = ('  !!{RST\n'
            '  .. math::\n'
            '\n'
            '     a = b \\\\ c = d\n'
            '  !!}\n')
    assert kinds(text) == []


def test_latex_character_escapes_are_detected():
    text = '  !!{RST\n  See issue \\#695 and 50\\% of foo\\_bar.\n  !!}\n'
    assert kinds(text) == ['escape:#', 'escape:%', 'escape:_']


def test_line_numbers_are_reported_against_the_source():
    text = ('module m\n'
            '  !!{RST\n'
            '  Fine.\n'
            '  Uses {\\normalfont \\ttfamily node}.\n'
            '  !!}\n')
    found = [f for f in residualLatex.scan_text(text, 'f.F90')
             if f.kind == 'font-group']
    assert [f.line for f in found] == [4]


# ---------------------------------------------------------------------------
# Baseline semantics
# ---------------------------------------------------------------------------

def _finding(path, kind):
    return residualLatex.Finding(path, 1, 'docstring', kind, '')


def test_baseline_key_ignores_line_number():
    """Keying on (path, kind) keeps the baseline from churning every time an
    unrelated edit shifts code up or down a file."""
    a = residualLatex.Finding('f.F90', 10, 'docstring', 'command:AA', '')
    b = residualLatex.Finding('f.F90', 99, 'docstring', 'command:AA', '')
    assert a.key == b.key


def test_compare_flags_a_rise_above_baseline():
    findings = [_finding('f.F90', 'command:AA')] * 3
    regressions, _ = residualLatex.compare(findings, {('f.F90', 'command:AA'): 2})
    assert regressions == [(('f.F90', 'command:AA'), 2, 3)]


def test_compare_flags_an_unlisted_pair():
    regressions, _ = residualLatex.compare([_finding('f.F90', 'font-group')], {})
    assert regressions == [(('f.F90', 'font-group'), 0, 1)]


def test_compare_accepts_counts_at_or_below_baseline():
    findings = [_finding('f.F90', 'command:AA')]
    regressions, improvements = residualLatex.compare(
        findings, {('f.F90', 'command:AA'): 2})
    assert regressions == []
    assert improvements == [(('f.F90', 'command:AA'), 2, 1)]


def test_compare_reports_a_fully_fixed_pair_as_improvement():
    _regressions, improvements = residualLatex.compare(
        [], {('f.F90', 'font-group'): 4})
    assert improvements == [(('f.F90', 'font-group'), 4, 0)]


def test_baseline_round_trips(tmp_path):
    findings = [_finding('a.F90', 'command:AA'), _finding('a.F90', 'command:AA'),
                _finding('b.F90', 'font-group')]
    path = str(tmp_path / 'baseline.txt')
    residualLatex.write_baseline(findings, path)
    assert residualLatex.load_baseline(path) == {('a.F90', 'command:AA'): 2,
                                                 ('b.F90', 'font-group'): 1}
    regressions, improvements = residualLatex.compare(
        findings, residualLatex.load_baseline(path))
    assert regressions == [] and improvements == []


def test_missing_baseline_is_an_empty_baseline(tmp_path):
    assert residualLatex.load_baseline(str(tmp_path / 'absent.txt')) == {}


# ---------------------------------------------------------------------------
# The committed baseline
# ---------------------------------------------------------------------------

def test_committed_baseline_matches_the_tree():
    """The whole point of the guard: ``source`` must not drift above the
    recorded baseline.  This mirrors the CI step, so a regression is caught
    by ``pytest`` too."""
    root = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir)
    source = os.path.join(root, 'source')
    if not os.path.isdir(source):                       # pragma: no cover
        return
    findings = residualLatex.scan_paths([source])
    baseline = residualLatex.load_baseline()
    # Paths in the baseline are recorded relative to the repository root.
    rebased = [residualLatex.Finding(
        os.path.relpath(f.path, root), f.line, f.context, f.kind, f.snippet)
        for f in findings]
    regressions, _ = residualLatex.compare(rebased, baseline)
    assert regressions == [], f'residual LaTeX above baseline: {regressions}'
