"""Tests for the `<methods>` metadata `bound_function_table` emits.

Each type-bound binding carries a return type and an argument list alongside its
description.  Both used to be dropped when the `<method …/>` element was
written, which left the generated documentation with nothing but a name and a
one-line description; these tests pin the emission.
"""

from Galacticus.Build.Components import bound_function_table


def _methods_block(bindings):
    """Return just the `<methods>` block of a rendered binding table."""
    out = bound_function_table('nodeComponent', bindings)
    assert '<methods docformat="rst">' in out
    return out[out.index('<methods'):out.index('</methods>')]


def test_return_type_and_arguments_are_emitted():
    block = _methods_block([{
        'type':        'procedure',
        'name':        'massDistribution',
        'function':    'Node_Component_Mass_Distribution',
        'description': "Return the mass distribution for this component.",
        'returnType':  r"``class(massDistribution)``",
        'arguments':   (r"``type(enumerationComponentTypeType)`` componentType "
                        r"(optional) [in], ``integer`` weightIndex (optional) [in]"),
    }])
    assert '<method method="massDistribution" description="Return the mass ' \
           'distribution for this component.">' in block
    assert '<type>``class(massDistribution)``</type>' in block
    assert '<argument>``type(enumerationComponentTypeType) componentType`` ' \
           '(optional) [in]</argument>' in block
    assert '<argument>``integer weightIndex`` (optional) [in]</argument>' in block


def test_argument_list_splits_at_arguments_not_at_every_comma():
    """An array shape may contain commas of its own."""
    block = _methods_block([{
        'type':        'procedure',
        'name':        'serialize',
        'function':    'Node_Component_Serialize',
        'description': "Serialize the component.",
        'returnType':  r"``void``",
        'arguments':   r"``double precision(:,:)`` array [out], ``integer`` count [in]",
    }])
    assert '<argument>``double precision(:,:) array`` [out]</argument>' in block
    assert '<argument>``integer count`` [in]</argument>' in block


def test_void_return_emits_no_type():
    """`void` names no return type — a subroutine has nothing to return."""
    block = _methods_block([{
        'type':        'procedure',
        'name':        'destroy',
        'function':    'Node_Component_Generic_Destroy',
        'description': "Destroy the object.",
        'returnType':  r"``void``",
        'arguments':   "",
    }])
    assert '<type>' not in block
    # With neither a type nor arguments, the short form is still used.
    assert '<method method="destroy" description="Destroy the object."/>' in block


def test_pointer_return_is_spelled_out():
    block = _methods_block([{
        'type':        'procedure',
        'name':        'host',
        'function':    'Node_Component_Host_Node',
        'description': "Return a pointer to the host node.",
        'returnType':  r"``*type(treeNode)``",
        'arguments':   "",
    }])
    assert '<type>``type(treeNode)`` (pointer)</type>' in block


def test_metadata_on_a_descriptor_is_emitted_too():
    """Property bindings carry their documentation on a `descriptor`."""
    block = _methods_block([{
        'type':       'procedure',
        'name':       'massGet',
        'descriptor': {
            'name':        'Node_Component_Disk_Mass_Get',
            'methodName':  'mass',
            'description': "Get the mass.",
            'returnType':  r"``double precision``",
            'arguments':   "",
        },
    }])
    assert '<method method="mass" description="Get the mass.">' in block
    assert '<type>``double precision``</type>' in block


def test_bindings_without_a_description_are_not_documented():
    """Only described bindings appear: an undescribed one is an implementation
    detail, and a block with none at all is not emitted."""
    out = bound_function_table('nodeComponent', [{
        'type':     'procedure',
        'name':     'internal',
        'function': 'Node_Component_Internal',
    }])
    assert '<methods' not in out
    assert 'procedure :: internal => Node_Component_Internal' in out
