# Regression test for `_build_allowed_parameters_method`'s
# objectBuilder → pointer-member matching.
#
# Bug: when the discovery pass walks each `<objectBuilder class="…"
# name="X_" .../>` directive in a class's constructor and tries to
# match it against the pointer declarations on the class type, the
# loop iterated `dec.get('variables')` directly.  But
# `parse_declaration` runs in `keep_qualifiers=True` mode here, so a
# pointer-with-initializer declaration like
#
#     class(galacticFilterClass), pointer :: galacticFilter_ => null()
#
# arrives as `variables = ['galacticfilter_=>null()']` — the bare
# pointer name with the `=>null()` tail still attached.  The loop
# compared `v.lower() == striplc(d.get('name'))` (where
# `d.get('name')` is the bare name `'galacticFilter_'`), so the
# `=>null()` tail made the comparison fail, the `objects`
# accumulator stayed empty, and the auto-generated
# `*AllowedParameters` method silently dropped its
# `if (associated(self%X_)) call self%X_%allowedParameters
# (allowedParameters,'parameters',.true.)` recursive forwarding
# lines.
#
# Fix: run each `v` through `strip_variable_name` to remove the
# `=>...` / `=...` tail before the comparison (and store the bare
# form in the accumulator).

import inspect
import re


def test_objectBuilder_loop_strips_variable_initializer():
    """Pin the strip in the objectBuilder match loop in the shipped
    source: any reversion that goes back to the bare `v.lower()`
    comparison instantly stops matching pointer-with-initializer
    declarations."""
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _build_allowed_parameters_method,
    )
    src = inspect.getsource(_build_allowed_parameters_method)
    # Look for the exact comparison line that now uses the stripped
    # bare name.  If anyone reverts the strip, this assertion fires.
    assert 'bare = strip_variable_name(v)' in src, src
    # And the comparison that follows uses `bare`, not raw `v`.
    assert 'bare.lower()' in src, src


def test_strip_variable_name_drops_pointer_initializer():
    """Demonstrate the underlying need: the variable name as it appears
    in the parsed declaration carries `=>null()`, but the directive's
    `name` attribute is the bare form."""
    from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
        strip_variable_name,
    )
    assert strip_variable_name('galacticfilter_=>null()') == 'galacticfilter_'
    assert strip_variable_name('counter=0')               == 'counter'
    assert strip_variable_name('plain_var')               == 'plain_var'
