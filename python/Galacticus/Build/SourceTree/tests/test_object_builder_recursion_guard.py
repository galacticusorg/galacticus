"""Regression tests for the object-build recursion guards (issue #397).

Two protections are generated into the code to prevent an unbounded
recursion when an object that composites a member of its own class
(directly, or via another, mutually-compositing class) is not provided
explicitly in the parameter file:

  1. ObjectBuilder.py emits an informative same-class guard at each
     ``objectBuilder`` site where ``self`` is in scope and is of the
     class being built (catches direct self-composition early, with a
     message telling the user to provide the object explicitly).

  2. FunctionClass factories push the parameter node being built onto a
     thread-private build stack (``Input_Parameters_Build_Stack_Push``)
     before dispatching to a constructor and pop it afterwards
     (``Input_Parameters_Build_Stack_Pop``).  Because every build --
     whether reached via an ``objectBuilder`` directive or via a
     ``functionGlobal`` wrapper -- routes through the factory, this is the
     general guard that also catches cross-class mutual recursion.

These are "pin the implementation" tests: if either guard is removed or
its message/identifier reverted, the corresponding assertion fires
immediately, without needing a full Galacticus build.
"""

import inspect


def test_object_builder_emits_informative_same_class_message():
    """The same-class guard in ObjectBuilder.py must abort with a message
    that tells the user the object must be provided explicitly -- not the
    old cryptic 'recursive build ... detected' text."""
    from Galacticus.Build.SourceTree.Process.ObjectBuilder import (
        _handle_object_builder,
    )
    src = inspect.getsource(_handle_object_builder)
    # The recursive-build guard is still keyed on the source%parent test ...
    assert "associated(parametersCurrent," in src, src
    # ... but now reports an actionable, explicit-provision message.
    assert "composites a member of its own class" in src, src
    assert "provide a [" in src.lower() or "provide a [" in src, src


def test_factory_pushes_and_pops_build_stack():
    """Every generated functionClass factory must bracket its constructor
    dispatch with the build-stack push/pop so that recursive (including
    cross-class and functionGlobal-mediated) builds are detected."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    assert "Input_Parameters_Build_Stack_Push(subParameters%parameters" in src, src
    assert "Input_Parameters_Build_Stack_Pop()" in src, src
    # The push/pop helpers must be imported into the factory's use list.
    assert "Input_Parameters_Build_Stack_Push" in src and \
           "Input_Parameters_Build_Stack_Pop" in src, src


def test_build_stack_push_is_placed_after_recursive_short_circuit():
    """The push must live in the dispatch paths (default-add branch and the
    main 'select case' ladder), i.e. after the ``recursive="yes"``
    short-circuit return, so legitimately self-recursive classes (which
    return a recursiveSelf reference before dispatching) are not flagged."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _generate_constructor,
    )
    src = inspect.getsource(_generate_constructor)
    push = "Input_Parameters_Build_Stack_Push(subParameters%parameters"
    ladder = "select case (char(instanceName))"
    short_circuit = "RecursiveBuildNode)) then"
    # The recursive short-circuit is generated before the main ladder, and
    # the ladder push follows the 'select case' opener.
    assert src.index(short_circuit) < src.index(ladder), src
    assert src.index(push) < src.index(ladder) and \
           src.rindex(push) > src.index(short_circuit), src
