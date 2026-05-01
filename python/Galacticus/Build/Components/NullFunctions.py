# Generate the on-demand "null" type-bound functions used by the
# components-build pipeline.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/NullFunctions.pm.  A null
# function is a placeholder type-bound implementation that satisfies a
# Fortran type's contract for a property method (get / set / rate / scale
# / analytic / inactive) when no real implementation is needed.  The
# generated function names are deduplicated by a fingerprint key so that
# any two callers asking for the same shape share one function.

import os
import sys


from Galacticus.Build.Components.Utils import (
    is_intrinsic,
    intrinsic_types,
    intrinsic_nulls,
)


# Module-level fingerprint cache.  Mirrors Perl `my %nullFunctionFingerprints`
# at NullFunctions.pm:20: a single build process generates each unique null
# function exactly once.
_null_function_fingerprints = set()


def create_null_function(build, descriptor):
    """Return the name of a null function matching `descriptor`, creating
    the function on `build['functions']` if it has not been emitted yet.

    `descriptor` is a dict with keys:

    * `selfType`  — `"generic"` or a component-class name; controls the
                    Fortran `class(nodeComponent…)` of the `self` argument.
    * `attribute` — one of `"get"`, `"set"`, `"rate"`, `"scale"`,
                    `"analytic"`, `"inactive"`.
    * `intent`    — `"in"` / `"inout"` / `"out"` for the `self` argument.
    * `property`  — sub-dict carrying at least `type` and `rank`.

    Mirrors Perl `createNullFunction`.
    """
    prop = descriptor['property']
    fingerprint = ":".join(str(descriptor[k]) for k in ('selfType', 'attribute', 'intent')) \
                + ":" + ":".join(str(prop[k]) for k in ('type', 'rank'))

    function_name = (
        "null"
        + ''.join(_ucfirst(str(descriptor[k])) for k in ('selfType', 'attribute', 'intent'))
        + ''.join(_ucfirst(str(prop      [k])) for k in ('type',     'rank'                ))
    )

    if fingerprint in _null_function_fingerprints:
        return function_name
    _null_function_fingerprints.add(fingerprint)

    self_type = "nodeComponent" + (
        "" if descriptor['selfType'] == "generic" else descriptor['selfType']
    )

    # Build the descriptor for the property argument.
    prop_descriptor = {}
    if is_intrinsic(prop['type']):
        prop_descriptor['intrinsic'] = intrinsic_types[prop['type']]
    else:
        prop_descriptor['intrinsic'] = 'type'
        prop_descriptor['type'     ] = prop['type']

    rank = int(prop['rank'])
    rank_attribute = (
        ['dimension(' + ','.join([':'] * rank) + ')'] if rank > 0 else []
    )

    self_var = {
        'intrinsic':  'class',
        'type':       self_type,
        'attributes': [f"intent({descriptor['intent']})"],
        'variables':  ['self'],
    }

    variables  = []
    modules    = []
    return_type = None
    attribute   = descriptor['attribute']

    if attribute == 'rate':
        prop_descriptor['variables']  = ['setValue']
        prop_descriptor['attributes'] = ['intent(in   )', *rank_attribute]
        return_type = 'void'
        variables = [
            self_var,
            prop_descriptor,
            {
                'intrinsic':  'logical',
                'attributes': ['intent(inout)', 'optional'],
                'variables':  ['interrupt'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'interruptTask',
                'attributes': ['intent(inout)', 'optional', 'pointer'],
                'variables':  ['interruptProcedure'],
            },
        ]
    elif attribute in ('set', 'scale'):
        prop_descriptor['variables']  = ['setValue']
        prop_descriptor['attributes'] = ['intent(in   )', *rank_attribute]
        return_type = 'void'
        variables = [self_var, prop_descriptor]
        modules   = ['Error']
    elif attribute == 'get':
        prop_descriptor['variables']  = ['setValue']
        prop_descriptor['attributes'] = ['intent(in   )', *rank_attribute]
        head = (f"type({prop_descriptor['type']})"
                if 'type' in prop_descriptor
                else prop_descriptor['intrinsic'])
        rank_text = (
            ", dimension(" + ",".join([":"] * rank) + "), allocatable"
            if rank > 0 else ""
        )
        return_type = f"{head}{rank_text} => getValue"
        variables = [self_var]
    elif attribute in ('analytic', 'inactive'):
        return_type = 'void'
        variables = [self_var]
    else:
        raise RuntimeError(
            f"createNullFunction: attribute '{attribute}' not supported"
        )

    # Mark every variable as unused for GCC, matching Perl's `!$GLC
    # attributes unused :: …` line.
    unused_names = []
    for v in variables:
        unused_names.extend(v.get('variables') or [])
    unused_line = (
        "!$GLC attributes unused :: " + ", ".join(unused_names) + "\n"
    )

    function = {
        'type':        return_type,
        'name':        function_name,
        'description': (
            f"A null {attribute} rate function for a rank {prop['rank']} "
            f"\\mono{{{self_type.lower()}}} class.\n"
        ),
        'variables':   variables,
        'content':     unused_line,
    }
    if modules:
        function['modules'] = modules

    if attribute == 'get':
        if is_intrinsic(prop['type']):
            null_value = intrinsic_nulls[prop['type']]
        else:
            null_value = f"null{_ucfirst(prop['type'])}{prop['rank']}d"
        function['content'] += f"getValue={null_value}\n"
    if attribute == 'set':
        function['content'] += (
            "call Error_Report('attempt to set value in null component'"
            "//{introspection:location})\n"
        )

    build.setdefault('functions', []).append(function)
    return function_name


def reset_fingerprints():
    """Clear the fingerprint cache.

    Useful for tests that exercise `create_null_function` repeatedly: since
    the cache lives at module scope, the second test would see a stale
    `_null_function_fingerprints` and skip emitting the function.
    """
    _null_function_fingerprints.clear()


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
