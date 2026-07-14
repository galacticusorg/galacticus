"""Characterization tests for the legitimate-recursion machinery (issue #695).

Three `functionClass` families opt in to a *bounded* construction cycle via
``recursive="yes"`` on their directive (``darkMatterHaloScale``,
``virialDensityContrast``, and the ``starFormationHistory`` family). When the
parameters-driven build re-enters the node currently under construction, the
generated factory returns a lightweight "shim": an instance of the concrete
recursive type carrying ``isRecursive=.true.`` and a weak ``recursiveSelf``
pointer back to the instance already under construction. Every public method
of the hand-written class forwards to ``recursiveSelf`` when ``isRecursive``.

Detection is unified onto the shared object-build stack (issue #695 Phase 1):
the factory records the in-progress object on its build-stack entry
(``Input_Parameters_Build_Stack_Object_Set``) after allocation but before
dispatch, and a re-entrant build detects the cycle by querying the stack
(``Input_Parameters_Build_Stack_Recursive_Object``) rather than via the
former per-family thread-private ``RecursiveBuildNode``/``RecursiveBuildObject``
module variables. These tests pin that mechanism so a later change to it is
deliberate. They run against the generator source directly, without a full
Galacticus build.

See also ``test_object_builder_recursion_guard.py``, which pins the
complementary #397 *abort* path for unintended (unbounded) cycles.
"""

import inspect


def test_factory_uses_build_stack_not_module_variables():
    """Detection must go through the shared build stack, not the retired
    per-family ``RecursiveBuildNode``/``RecursiveBuildObject`` module
    variables (which are no longer emitted)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    # The old per-family thread-private module variables must no longer be
    # emitted (their declarations and threadprivate directive are gone).
    assert "RecursiveBuildNode => null()" not in src, src
    assert "threadprivate({directive_name}RecursiveBuildNode" not in src, src
    assert "RecursiveBuildObject => self" not in src, src
    # Detection now goes through the shared build stack.
    assert "Input_Parameters_Build_Stack_Recursive_Object" in src, src
    assert "Input_Parameters_Build_Stack_Object_Set" in src, src
    # Emission remains gated on any class in the family being recursive.
    assert "c.get('recursive') == 'yes' for c in classes_ordered" in src, src


def test_factory_records_object_under_construction_before_dispatch():
    """Before dispatching to the inner constructor of a recursive class, the
    factory must record the in-progress object on its build-stack entry, so
    that a re-entrant build can retrieve it. The record happens after
    ``allocate`` (self exists) but before the constructor dispatch."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    assert "Input_Parameters_Build_Stack_Object_Set(self)" in src, src
    # The record is guarded on the concrete class being recursive.
    assert "c.get('recursive') == 'yes'" in src, src
    # It is placed after allocation and before the dispatch select-type.
    obj_set = "Input_Parameters_Build_Stack_Object_Set(self)"
    dispatch = "self={c['name']}(subParameters)"
    assert src.index(obj_set) < src.index(dispatch), src


def test_factory_short_circuits_reentrant_build_to_shim():
    """On re-entry for the node currently under construction, the factory
    must query the build stack for the in-progress object, and if found
    allocate the concrete recursive type as a shim: set
    ``isRecursive=.true.``, point ``recursiveSelf`` at that object, and
    return immediately (before the main build ladder)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    # The short-circuit fires on a positive build-stack query for this node ...
    assert "Input_Parameters_Build_Stack_Recursive_Object" in src, src
    assert "if (associated(recursiveObject)) then" in src, src
    # ... allocating a shim of the concrete type with the recursion flag set
    # and the weak back-pointer wired to the in-progress object.
    assert "self%isRecursive=.true." in src, src
    assert "self%recursiveSelf => recursiveObject" in src, src
    # The short-circuit returns before reaching the main 'select case' ladder.
    short_circuit = "if (associated(recursiveObject)) then"
    ladder = "select case (char(instanceName))"
    assert src.index(short_circuit) < src.index(ladder), src


def test_state_store_forwards_for_recursive_shim():
    """A ``recursive="yes"`` class's generated ``stateStore``/``stateRestore``
    must short-circuit when ``isRecursive`` and forward to ``recursiveSelf``
    (the shim holds no real state of its own)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _build_state_store_methods,
    )
    src = inspect.getsource(_build_state_store_methods)
    assert "non_abstract.get('recursive') == 'yes'" in src, src
    assert "if (self%isRecursive) then" in src, src
    assert "self%recursiveSelf%stateStore" in src, src
    assert "self%recursiveSelf%stateRestore" in src, src
