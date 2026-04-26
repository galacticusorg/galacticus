# Regression test for `_build_assignment_method` in
# `Galacticus.Build.SourceTree.Process.FunctionClass`.
#
# Two related bugs:
#
#   1. Type-bound bindings (`procedure :: methodName => methodImpl`) are
#      `intrinsic == 'procedure'` declarations, but they are NOT data
#      members.  An earlier draft emitted `self%methodName=from%methodName`
#      for them — invalid Fortran since `methodName` resolves to the
#      bound procedure, not a slot.
#
#   2. Procedure-pointer DATA members (`procedure(intf), pointer :: foo`)
#      are ALSO `intrinsic == 'procedure'`, but they MUST be copied — with
#      `=>` since they are pointers.  An over-aggressive fix for #1
#      skipped every `intrinsic == 'procedure'` declaration, leaving
#      procedure-pointer slots null in the assignment, and producing
#      runtime null-procedure-pointer segfaults at the call site.
#
# The correct discriminator is the `pointer` attribute: pointer ⇒ data
# member (emit `self%X=>from%X`), no pointer ⇒ type-bound binding (skip).

from Galacticus.Build.SourceTree.Process.FunctionClass import (
    _build_assignment_method,
)


def _decl(intrinsic, type_, attributes, variables):
    return {
        'intrinsic':     intrinsic,
        'type':          type_,
        'attributes':    list(attributes),
        'variables':     list(variables),
        'variableNames': list(variables),
    }


def _class_with_declarations(class_name, declarations):
    """Build the minimal AST shape `_build_assignment_method` walks.

    The function looks up each non-abstract class's `tree.firstChild`,
    walks the sibling chain for a `type` node whose `name` matches the
    class, then walks that node's child chain for `declaration` nodes.
    """
    type_node = {
        'type':       'type',
        'name':       class_name,
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    declaration_node = {
        'type':         'declaration',
        'declarations': declarations,
        'firstChild':   None,
        'sibling':      None,
        'parent':       type_node,
    }
    type_node['firstChild'] = declaration_node
    tree = {'firstChild': type_node}
    return {
        'name':    class_name,
        'extends': '__top__',           # walk terminator (matches directive name)
        'tree':    tree,
        'node':    type_node,
    }


def _run(declarations):
    cls = _class_with_declarations('myClass', declarations)
    methods   = {}
    directive = {'name': '__top__', 'data': []}
    _build_assignment_method(
        directive,
        non_abstract_classes=[cls],
        classes={'myClass': cls},
        methods=methods,
        state_storables={},
    )
    return methods['assignment(=)']['code']


def test_procedure_pointer_data_member_uses_arrow_assignment():
    """`procedure(intf), pointer, nopass :: callback` is a data member and
    must be copied with `=>`."""
    code = _run([_decl(
        'procedure', 'myInterface',
        ['pointer', 'nopass'],
        ['callback'],
    )])
    assert 'self%callback=>from%callback' in code, code


def test_type_bound_procedure_binding_is_skipped():
    """`procedure :: bind => bindImpl` is a type-bound binding (no
    `pointer` attribute) — it must NOT appear in the generated assignment."""
    code = _run([_decl(
        'procedure', None,
        [],   # no `pointer`
        ['bind => bindImpl'],
    )])
    assert 'self%bind' not in code, code


def test_generic_and_final_bindings_are_skipped():
    code = _run([
        _decl('generic', None, [], ['accept => acceptInt, acceptReal']),
        _decl('final',   None, [], ['destructor']),
    ])
    assert 'self%accept'     not in code, code
    assert 'self%destructor' not in code, code


def test_plain_data_member_uses_equals_assignment():
    """A scalar `integer :: counter` is a regular data member; copy with `=`."""
    code = _run([_decl('integer', None, [], ['counter'])])
    assert 'self%counter=from%counter' in code, code


def test_pointer_and_type_bound_in_same_class():
    """The two kinds coexist in real classes (e.g. the spherical-Gnedin
    initializer that triggered the segfault) — make sure both rules fire
    on the same `<declaration>`-block walk."""
    code = _run([_decl(
        'procedure', 'initInterface',
        ['pointer', 'nopass'],
        ['initializationFunction'],
    ), _decl(
        'procedure', None,
        [],
        ['describe => describeImpl'],
    )])
    assert 'self%initializationFunction=>from%initializationFunction' in code
    assert 'self%describe'     not in code
