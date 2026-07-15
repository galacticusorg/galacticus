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
    # ... allocating the generated shim type with the weak back-pointer wired
    # to the in-progress object.
    assert "allocate({shim_type} :: self)" in src, src
    assert "self%recursiveSelf => recursiveObject" in src, src
    # The short-circuit returns before reaching the main 'select case' ladder.
    short_circuit = "if (associated(recursiveObject)) then"
    ladder = "select case (char(instanceName))"
    assert src.index(short_circuit) < src.index(ladder), src


def test_state_store_forwards_for_recursive_shim():
    """State storage for a recursive re-entry is handled by the generated shim
    type overriding stateStore/stateRestore to forward to recursiveSelf; the
    concrete class's generated stateStore no longer carries an isRecursive
    short-circuit."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _build_state_store_methods, _generate_recursive_shim,
    )
    # The concrete-class generated stateStore no longer references isRecursive.
    store_src = inspect.getsource(_build_state_store_methods)
    assert "self%isRecursive" not in store_src, store_src
    # The shim forwards stateStore/stateRestore instead.
    shim_src = inspect.getsource(_generate_recursive_shim)
    assert "self%recursiveSelf%{store}" in shim_src, shim_src
    assert "for store in ('stateStore', 'stateRestore')" in shim_src, shim_src


# --- Phase 2: the generated shim type <name>Recursive (issue #695) ------------

def test_shim_type_is_generated_for_recursive_families():
    """A family containing a recursive="yes" class emits a dedicated shim type
    ``type, extends(<name>Class) :: <name>Recursive`` holding only a weak
    ``recursiveSelf`` back-pointer (plus the deferred-copy flag)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_recursive_shim,
    )
    src = inspect.getsource(_generate_recursive_shim)
    assert "extends({directive_name}Class) :: {shim_type}" in src, src
    assert "shim_type = directive_name + 'Recursive'" in src, src
    assert "recursiveSelf => null()" in src, src
    assert "logical :: parentDeferred" in src, src


def test_shim_forwards_methods_and_no_ops_hooks():
    """The shim forwards every physics method (and stateStore/descriptor/
    objectType) to recursiveSelf, and makes autoHook and allowedParameters
    no-ops (fixing the double-hook wart and avoiding contributing parameters
    the real object already owns)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_recursive_shim,
    )
    src = inspect.getsource(_generate_recursive_shim)
    # Physics methods forward to the real object.
    assert "self%recursiveSelf%{m}" in src, src
    # stateStore / descriptor / objectType forward.
    assert "self%recursiveSelf%objectType" in src, src
    assert "self%recursiveSelf%descriptor" in src, src
    assert "self%recursiveSelf%{store}" in src, src
    # autoHook and allowedParameters are no-ops (bodies just return).
    assert "AutoHook(self)" in src, src
    assert "AllowedParameters(self" in src, src


def test_shim_deepcopy_uses_copiedself_fixup():
    """The shim's deepCopy wires recursiveSelf to the copy of the real object
    via copiedSelf (D4), deferring to deepCopyFinalize when the real object is
    not yet copied."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_recursive_shim,
    )
    src = inspect.getsource(_generate_recursive_shim)
    assert "self%recursiveSelf%copiedSelf" in src, src
    assert "parentDeferred =  .true." in src, src
    assert "DeepCopyFinalize(self)" in src, src
    # A mold-allocated copy starts at referenceCount=0; the shim's deepCopy must
    # reset it to 1 (like the generated deepCopy) or the owner's objectDestructor
    # decrements to -1 and aborts.
    assert "destination%referenceCountReset()" in src, src


def test_deepcopy_finalize_does_not_null_copiedself():
    """D4: the generated deepCopyFinalize must NOT null copiedSelf (it must
    survive the finalize pass so a shim can resolve its deferred parent)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _build_deep_copy_methods,
    )
    src = inspect.getsource(_build_deep_copy_methods)
    # reset still nulls copiedSelf ...
    assert "deep_copy['resetCode']    += \"self%copiedSelf => null()\\n\"" in src, src
    # ... but finalize does not.
    assert (
        "deep_copy['finalizeCode'] += \"self%copiedSelf => null()\\n\""
        not in src
    ), src


def test_factory_short_circuit_allocates_shim_type():
    """The factory short-circuit allocates the generated shim type and wires
    recursiveSelf for every recursive family (the old concrete-type-with-
    isRecursive-flag shim is gone)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    assert "allocate({shim_type} :: self)" in src, src
    assert "self%recursiveSelf => recursiveObject" in src, src
    # The old concrete-type-with-flag path must be gone.
    assert "self%isRecursive=.true." not in src, src
