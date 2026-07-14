"""Characterization tests for the legitimate-recursion machinery (issue #695).

Three `functionClass` families opt in to a *bounded* construction cycle via
``recursive="yes"`` on their directive (``darkMatterHaloScale``,
``virialDensityContrast``, and the ``starFormationHistory`` family). When the
parameters-driven build re-enters the node currently under construction, the
generated factory returns a lightweight "shim": an instance of the concrete
recursive type carrying ``isRecursive=.true.`` and a weak ``recursiveSelf``
pointer back to the instance already under construction. Every public method
of the hand-written class forwards to ``recursiveSelf`` when ``isRecursive``.

These tests *pin today's mechanism* (the concrete-type-with-flag shim, the
per-family ``<name>RecursiveBuildNode``/``<name>RecursiveBuildObject``
thread-private module variables, and the ``stateStore``/``stateRestore``
forwarding) so that the planned refactor to a generated shim type (issue #695
Phases 1-3) is a deliberate, reviewed change rather than a silent one. They
run against the generator source directly, without a full Galacticus build.

See also ``test_object_builder_recursion_guard.py``, which pins the
complementary #397 *abort* path for unintended (unbounded) cycles.
"""

import inspect


def test_factory_emits_recursive_build_module_variables():
    """A family containing a ``recursive="yes"`` class must emit the
    thread-private ``<name>RecursiveBuildNode``/``<name>RecursiveBuildObject``
    module variables that carry the node + object currently under
    construction across the re-entrant build."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    # Emission is gated on any class in the family being recursive.
    assert "c.get('recursive') == 'yes' for c in classes_ordered" in src, src
    assert "RecursiveBuildNode => null()" in src, src
    assert "RecursiveBuildObject => null()" in src, src
    # ... and the pair is thread-private (one per OpenMP thread).
    assert "threadprivate({directive_name}RecursiveBuildNode," in src, src


def test_factory_records_object_under_construction_before_dispatch():
    """Before dispatching to the inner constructor of a recursive class, the
    factory must point ``<name>RecursiveBuildObject`` at ``self`` (and
    ``<name>RecursiveBuildNode`` at the parameter node) so that a re-entrant
    build can discover the instance already under construction, then clear
    them afterwards."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    assert "RecursiveBuildNode   => parameterNode" in src, src
    assert "RecursiveBuildObject => self" in src, src
    # The set is guarded on the concrete class being recursive.
    assert "c.get('recursive') == 'yes'" in src, src


def test_factory_short_circuits_reentrant_build_to_shim():
    """On re-entry for the node currently under construction, the factory
    must allocate the concrete recursive type as a shim: set
    ``isRecursive=.true.``, point ``recursiveSelf`` at the object under
    construction, and return immediately (before the main build ladder)."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    # The short-circuit fires when the node matches the one being built ...
    assert "associated(parameterNode," in src and \
           "RecursiveBuildNode)) then" in src, src
    # ... allocating a shim of the concrete type with the recursion flag set
    # and the weak back-pointer wired up.
    assert "self%isRecursive=.true." in src, src
    assert "self%recursiveSelf => " in src, src
    assert "RecursiveBuildObject" in src, src
    # The short-circuit returns before reaching the main 'select case' ladder.
    short_circuit = "RecursiveBuildNode)) then"
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
