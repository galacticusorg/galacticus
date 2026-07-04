"""Regression tests for libraryInterfacesAudit.py's methods pass.

Two bugs made the audit's method-level report silently misleading — the
exact report a developer consults to prioritise library-interface coverage
work:

* ``_FC_BLOCK_RX`` required a *bare* ``<functionClass>`` open tag, but
  every real directive carries attributes (``<functionClass
  docformat="rst">``, …).  The block therefore matched nothing and the
  whole methods pass reported zero methods discovered — the coverage
  report looked "all clear" while actually seeing none of the source.

* ``_is_out_of_scope_reason`` decided whether a concrete-class blocker was
  a deferred non-fc hierarchy by searching the *entire* reason string for
  ``class(*)``.  But the explanatory tail of every concrete-class blocker
  reads "…and class(*) (with override) are supported", so the search
  always matched and every ``class(<non-fc>)`` blocker was mis-bucketed as
  in-scope — inflating the actionable worklist ~3x with work the team had
  explicitly deferred.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import libraryInterfacesAudit as audit  # noqa: E402


# ---------------------------------------------------------------------------
# _FC_BLOCK_RX — must match attributed functionClass open tags.
# ---------------------------------------------------------------------------

_ATTRIBUTED_BLOCK = (
    '<functionClass docformat="rst">\n'
    '  <name>transferFunction</name>\n'
    '  <method name="value">\n'
    '   <type>double precision</type>\n'
    '  </method>\n'
    '</functionClass>'
)


def test_fc_block_matches_attributed_open_tag():
    # The real source form carries attributes; a matcher requiring a bare
    # `<functionClass>` sees nothing and the methods pass reports zero.
    m = audit._FC_BLOCK_RX.search(_ATTRIBUTED_BLOCK)
    assert m is not None
    assert '<name>transferFunction</name>' in m.group('body')


def test_fc_block_still_matches_bare_open_tag():
    bare = '<functionClass>\n  <name>foo</name>\n</functionClass>'
    m = audit._FC_BLOCK_RX.search(bare)
    assert m is not None
    assert '<name>foo</name>' in m.group('body')


def test_discover_methods_finds_attributed_directive(tmp_path):
    src = tmp_path / 'thing.F90'
    src.write_text(
        'module Thing\n'
        '  !![\n' + _ATTRIBUTED_BLOCK + '\n  !!]\n'
        'end module Thing\n'
    )
    methods = audit.discover_methods(tmp_path, {'transferFunction'})
    assert [m['name'] for m in methods['transferFunction']] == ['value']


# ---------------------------------------------------------------------------
# _is_out_of_scope_reason — concrete non-fc class blockers are deferred;
# the genuine in-scope categories are not.
# ---------------------------------------------------------------------------

# The concrete-class blocker string classify_arg emits: note the tail
# literally contains "class(*)", which must NOT flip the verdict.
_CONCRETE_CLASS = (
    "class(coordinate) — only registered functionClasses and class(*) "
    "(with override) are supported (coordinates)"
)


def test_concrete_non_fc_class_is_out_of_scope():
    assert audit._is_out_of_scope_reason(_CONCRETE_CLASS) is True


def test_concrete_class_unregistered_fc_is_out_of_scope():
    reason = ("class(ompLockClass) — 'ompLock' is not a registered "
              "functionClass in libraryClasses.xml (lock)")
    assert audit._is_out_of_scope_reason(reason) is True


def test_class_star_blocker_is_in_scope():
    # class(*) just needs a method-side override hook — in-scope.
    reason = ("class(*) without a libraryClasses.xml override "
              "(<constructor><argument name=… type=… module=…/>) (obj)")
    assert audit._is_out_of_scope_reason(reason) is False


def test_class_non_fc_return_is_out_of_scope():
    reason = ("class(non-fc) return (class(fooClass) — 'foo' is not a "
              "registered functionClass)")
    assert audit._is_out_of_scope_reason(reason) is True


def test_scalar_derived_type_return_is_out_of_scope():
    assert audit._is_out_of_scope_reason(
        'scalar derived-type return type (type(abundances))') is True


def test_varying_string_is_in_scope():
    assert audit._is_out_of_scope_reason('type(varying_string)') is False


def test_output_array_arg_is_in_scope():
    reason = ('1D allocatable array argument (output arrays not yet '
              'supported) (wavenumbers)')
    assert audit._is_out_of_scope_reason(reason) is False


def test_procedure_pointer_arg_is_in_scope():
    reason = ('procedure(interruptTask) — procedure-pointer args are not '
              'supported (functionInterrupt)')
    assert audit._is_out_of_scope_reason(reason) is False


def test_complex_arg_is_in_scope():
    assert audit._is_out_of_scope_reason(
        'complex(c_double_complex) (windowFunction1)') is False
