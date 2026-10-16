"""Generate the on-demand "null" type-bound functions used by the
components-build pipeline.

Andrew Benson (ported to Python 2026)

A null
function is a placeholder type-bound implementation that satisfies a
Fortran type's contract for a property method (get / set / rate / scale
/ analytic / inactive) when no real implementation is needed.  The
generated function names are deduplicated by a fingerprint key so that
any two callers asking for the same shape share one function.
"""



from Galacticus.Build.Components.Utils import (
    is_intrinsic,
    intrinsic_types,
    intrinsic_nulls,
)


# Module-level fingerprint cache: a single build process generates each
# unique null function exactly once. `set` functions are generated per-property
# (see `create_null_function`), so the name each fingerprint was given is
# recorded, along with the number of functions sharing each base name, from
# which the suffix distinguishing them is formed - the property name itself
# would exceed the 63 character limit on a Fortran name.
_null_function_fingerprints = {}
_null_function_counts       = {}


def create_null_function(build, descriptor):
    """Return the name of a null function matching `descriptor`, creating
    the function on `build['functions']` if it has not been emitted yet.

    `descriptor` is a dict with keys:

    * `selfType`  — `"generic"` or a component-class name; controls the
                    Fortran `class(nodeComponent…)` of the `self` argument.
    * `attribute` — one of `"get"`, `"set"`, `"rate"`, `"scale"`,
                    `"analytic"`, `"inactive"`.
    * `intent`    — `"in"` / `"inout"` / `"out"` for the `self` argument.
    * `property`  — sub-dict carrying at least `type` and `rank`, and, for a
                    `set` attribute, `name`: setting a value in a null
                    component is an error, so those functions are generated
                    per-property in order to name the property responsible.
    """
    prop = descriptor['property']
    # A `set` function reports an error naming the component and property, so it can not be shared between properties.
    named = descriptor['attribute'] == 'set' and prop.get('name') is not None
    fingerprint = ":".join(str(descriptor[k]) for k in ('selfType', 'attribute', 'intent')) \
                + ":" + ":".join(str(prop[k]) for k in ('type', 'rank')) \
                + (":" + prop['name'] if named else "")

    function_name = (
        "null"
        + ''.join(_ucfirst(str(descriptor[k])) for k in ('selfType', 'attribute', 'intent'))
        + ''.join(_ucfirst(str(prop      [k])) for k in ('type',     'rank'                ))
    )

    if fingerprint in _null_function_fingerprints:
        return _null_function_fingerprints[fingerprint]
    count = _null_function_counts.get(function_name, 0)
    _null_function_counts[function_name] = count+1
    if count > 0:
        function_name += f"_{count}"
    _null_function_fingerprints[fingerprint] = function_name

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

    # Mark every variable as unused for GCC via a `!$GLC attributes
    # unused :: …` line.
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
            f"``{self_type.lower()}`` class.\n"
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
        # Name the component class and property, so that the error identifies what must be created (or not set).
        if named:
            target = (f"the `{prop['name']}` property of the null "
                      f"`{descriptor['selfType']}` component")
        else:
            target = "a value in a null component"
        function['content'] += (
            f"call Error_Report('attempt to set {target}'"
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
    _null_function_counts      .clear()


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
