# `treeNode.map*` methods that apply a function across every component
# of a node and reduce the results.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/TreeNodes/Map.pm.  Five hooks
# on the `functions` phase:
#
#   Tree_Node_Map_Optimizations — enumeration constants for optimized
#                                 map dispatch.
#   Tree_Node_Map_Void          — simple void-returning map.
#   Tree_Node_Map_Double0       — rank-0 double map with sum/product
#                                 reduction.
#   Tree_Node_Map_Double1       — rank-1 double map with sum/product
#                                 reduction.
#   Tree_Node_Map_TensorR2D3    — rank-2 dimension-3 tensor map with
#                                 summation reduction.
#
# `Double0`, `Double1`, and `TensorR2D3` share a common shape: one
# optimized branch per (boundFunction, reduction) pair drawn from any
# bound function whose descriptor carries a `mappable` key, followed by
# a generic fallback that switches on the reduction at run-time.

import os
import sys


from Galacticus.Build.Components.Utils import register


# ---------------------------------------------------------------------------
# Reduction tables
# ---------------------------------------------------------------------------

_DOUBLE_IDENTITY = {'summation': '0.0d0',             'product': '1.0d0'}
_DOUBLE_OPERATOR = {'summation': '+',                 'product': '*'    }
_TENSOR_IDENTITY = {'summation': 'tensorNullR2D3Sym'                    }
_TENSOR_OPERATOR = {'summation': '+'                                    }


def Tree_Node_Map_Optimizations(build):
    """Emit `optimizeFor<Method><Reduction>` integer parameters.

    Mirrors `Tree_Node_Map_Optimizations`.  Generates one parameter per
    (mappable bound function, reduction) pair, numbered from 0 upwards.
    """
    base_type = build.get('types', {}).get('nodeComponent') or {}
    if not _has_mappable_bound_functions(base_type):
        return

    value = -1
    for bf in base_type.get('boundFunctions') or []:
        if 'mappable' not in bf:
            continue
        for reduction in bf['mappable'].split(':'):
            value += 1
            name = (
                f"optimizeFor{_ucfirst(bf['name'])}{_ucfirst(reduction)}"
                f"={value}"
            )
            build.setdefault('variables', []).append({
                'intrinsic':  'integer',
                'attributes': ['public', 'parameter'],
                'variables':  [name],
            })


def Tree_Node_Map_Void(build):
    """Generate `treeNodeMapVoid`.

    Mirrors `Tree_Node_Map_Void`.  Iterates each active component class
    and applies `mapFunction(self%component<Class>(i))` to every member.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeMapVoid',
        'description': "Map a void function over components of the node.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'Node_Component_Null_Void0_InOut',
                'attributes': ['pointer'],
                'variables':  ['mapFunction'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"do i=1,size(self%component{cap})\n"
            f"  call mapFunction(self%component{cap}(i))\n"
            f"end do\n"
        )
    function['content'] = content
    _bind(build, 'mapVoid', function)


def Tree_Node_Map_Double0(build):
    """Generate `treeNodeMapDouble0`.  Mirrors `Tree_Node_Map_Double0`."""
    _build_reducing_map(
        build,
        function_name='treeNodeMapDouble0',
        method_name='mapDouble0',
        result_var='treeNodeMapDouble0',
        return_type='double precision',
        identity=_DOUBLE_IDENTITY,
        operator=_DOUBLE_OPERATOR,
        result_size=False,
        component_value_decls=[
            {
                'intrinsic': 'double precision',
                'variables': ['componentValue'],
            },
        ],
        map_arg_suffix='',
        generic_select=(
            "select case (reduction)\n"
            "case (reductionSummation)\n"
            "  treeNodeMapDouble0=0.0d0\n"
            "case (reductionProduct  )\n"
            "  treeNodeMapDouble0=1.0d0\n"
            "case default\n"
            "  treeNodeMapDouble0=1.0d0\n"
            "  call Error_Report('unknown reduction'//{introspection:location})\n"
            "end select\n"
        ),
        generic_per_class=(
            "  select case (reduction)\n"
            "  case (reductionSummation)\n"
            "    treeNodeMapDouble0=treeNodeMapDouble0+componentValue\n"
            "  case (reductionProduct  )\n"
            "    treeNodeMapDouble0=treeNodeMapDouble0*componentValue\n"
            "  end select\n"
        ),
        proc_type='Node_Component_Null_Double0_InOut',
        extra_in_args=[],
    )


def Tree_Node_Map_Double1(build):
    """Generate `treeNodeMapDouble1`.  Mirrors `Tree_Node_Map_Double1`."""
    _build_reducing_map(
        build,
        function_name='treeNodeMapDouble1',
        method_name='mapDouble1',
        result_var='double1Result',
        return_type=(
            'double precision, dimension(resultSize) => double1Result'
        ),
        identity=_DOUBLE_IDENTITY,
        operator=_DOUBLE_OPERATOR,
        result_size=True,
        component_value_decls=[
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(resultSize)'],
                'variables':  ['componentValue'],
            },
        ],
        map_arg_suffix=',resultSize',
        generic_select=(
            "select case (reduction)\n"
            "case (reductionSummation)\n"
            "  double1Result=0.0d0\n"
            "case (reductionProduct  )\n"
            "  double1Result=1.0d0\n"
            "case default\n"
            "  double1Result=1.0d0\n"
            "  call Error_Report('unknown reduction'//{introspection:location})\n"
            "end select\n"
        ),
        generic_per_class=(
            "  select case (reduction)\n"
            "  case (reductionSummation)\n"
            "    double1Result=double1Result+componentValue\n"
            "  case (reductionProduct  )\n"
            "    double1Result=double1Result*componentValue\n"
            "  end select\n"
        ),
        proc_type='Node_Component_Null_Double1_InOut',
        extra_in_args=['resultSize'],
    )


def Tree_Node_Map_TensorR2D3(build):
    """Generate `treeNodeMapTensorR2D3`.  Mirrors `Tree_Node_Map_TensorR2D3`."""
    _build_reducing_map(
        build,
        function_name='treeNodeMapTensorR2D3',
        method_name='mapTensorR2D3',
        result_var='tensorResult',
        return_type='type(tensorRank2Dimension3Symmetric) => tensorResult',
        identity=_TENSOR_IDENTITY,
        operator=_TENSOR_OPERATOR,
        result_size=False,
        component_value_decls=[
            {
                'intrinsic':  'type(tensorRank2Dimension3Symmetric)',
                'variables':  ['componentValue'],
            },
        ],
        map_arg_suffix='',
        generic_select=(
            "select case (reduction)\n"
            "case (reductionSummation)\n"
            "  tensorResult=tensorNullR2D3Sym\n"
            "case default\n"
            "  tensorResult=tensorNullR2D3Sym\n"
            "  call Error_Report('unknown reduction'//{introspection:location})\n"
            "end select\n"
        ),
        generic_per_class=(
            "  select case (reduction)\n"
            "  case (reductionSummation)\n"
            "    tensorResult=tensorResult+componentValue\n"
            "  end select\n"
        ),
        proc_type='Node_Component_Null_TensorR2D3_InOut',
        extra_in_args=[],
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_reducing_map(build, *, function_name, method_name, result_var,
                        return_type, identity, operator, result_size,
                        component_value_decls, map_arg_suffix,
                        generic_select, generic_per_class,
                        proc_type, extra_in_args):
    """Common skeleton for `Tree_Node_Map_Double0` / `Double1` /
    `TensorR2D3`.

    Each emits an optimized fast-path for every (mappable bound
    function, reduction) pair on the base `nodeComponent` type, then
    falls back to a generic body that switches on the run-time
    `reduction` argument.
    """
    in_args = ['reduction'] + list(extra_in_args)
    variables = [
        {
            'intrinsic':  'class',
            'type':       'treeNode',
            'attributes': ['intent(inout)'],
            'variables':  ['self'],
        },
        {
            'intrinsic':  'procedure',
            'type':       proc_type,
            'attributes': ['pointer'],
            'variables':  ['mapFunction'],
        },
    ]
    if extra_in_args:
        variables.append({
            'intrinsic':  'integer',
            'attributes': ['intent(in   )'],
            'variables':  ['resultSize', 'reduction'],
        })
    else:
        variables.append({
            'intrinsic':  'integer',
            'attributes': ['intent(in   )'],
            'variables':  ['reduction'],
        })
    variables.append({
        'intrinsic':  'integer',
        'attributes': ['intent(in   )', 'optional'],
        'variables':  ['optimizeFor'],
    })
    variables.extend(component_value_decls)
    variables.append({
        'intrinsic':  'integer',
        'variables':  ['i'],
    })

    base_type = build.get('types', {}).get('nodeComponent') or {}
    has_optimizations = _has_mappable_bound_functions(base_type)

    content = ''
    first_optimization = True
    for bf in base_type.get('boundFunctions') or []:
        if 'mappable' not in bf:
            continue
        for reduction in bf['mappable'].split(':'):
            if reduction not in identity:
                # Reduction not supported by this map kind (e.g. tensor
                # has no `product`). Match Perl: KeyError-on-undef,
                # which is benign when the table simply omits an entry.
                continue
            else_ = '' if first_optimization else 'else '
            content += (
                f"{else_}if (present(optimizeFor)"
                f".and.optimizeFor == optimizeFor"
                f"{_ucfirst(bf['name'])}{_ucfirst(reduction)}) then\n"
                f"    if (reduction /= reduction{_ucfirst(reduction)}) "
                f"call Error_Report('reduction mismatch'"
                f"//{{introspection:location}})\n"
                f"    {result_var}={identity[reduction]}\n"
            )
            for class_dict in _active_classes(build):
                if not _class_overrides(class_dict, bf['name']):
                    continue
                cap = _ucfirst(class_dict['name'])
                content += (
                    f"do i=1,size(self%component{cap})\n"
                    f"  {result_var}={result_var}{operator[reduction]}"
                    f"mapFunction(self%component{cap}(i){map_arg_suffix})\n"
                    f"end do\n"
                )
            first_optimization = False

    if has_optimizations:
        content += "else\n"
    content += generic_select
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"do i=1,size(self%component{cap})\n"
            f"  componentValue=mapFunction(self%component{cap}(i)"
            f"{map_arg_suffix})\n"
            f"{generic_per_class}"
            f"end do\n"
        )
    if has_optimizations:
        content += "end if\n"

    function = {
        'type':        return_type,
        'name':        function_name,
        'description': (
            "Map a rank-0, double function over components of the node."
            if function_name == 'treeNodeMapDouble1'
            else (
                "Map a rank-2, dimension-3 tensor function over components "
                "of the node."
                if function_name == 'treeNodeMapTensorR2D3'
                else (
                    "Map a rank-0, double function over components of "
                    "the node."
                )
            )
        ),
        'modules':     ['Error'],
        'variables':   variables,
        'content':     content,
    }
    _bind(build, method_name, function)


def _has_mappable_bound_functions(base_type):
    return any(
        'mappable' in bf
        for bf in base_type.get('boundFunctions') or []
    )


def _active_classes(build):
    """Yield each component-class dict whose name is in
    `componentClassListActive`.  Mirrors the Perl pattern of
    `foreach class ( hashList(componentClasses) ) { next unless grep ... }`.
    """
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        if class_dict['name'] in active:
            yield class_dict


def _class_overrides(class_dict, method_name):
    """Return True if any member of `class_dict` carries a binding whose
    `method` field matches `method_name`.

    Mirrors the Perl `grep {…} map {@{$_->{'bindings'}->{'binding'}}} @members`.
    """
    for member in class_dict.get('members') or []:
        bindings = (member.get('bindings') or {}).get('binding') or []
        if not isinstance(bindings, list):
            bindings = [bindings]
        for binding in bindings:
            if isinstance(binding, dict) and binding.get('method') == method_name:
                return True
    return False


def _bind(build, method_name, function):
    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Map.pm:23-29.
# ---------------------------------------------------------------------------

register('treeNodeMap', 'functions', Tree_Node_Map_Optimizations)
register('treeNodeMap', 'functions', Tree_Node_Map_Void)
register('treeNodeMap', 'functions', Tree_Node_Map_Double0)
register('treeNodeMap', 'functions', Tree_Node_Map_Double1)
register('treeNodeMap', 'functions', Tree_Node_Map_TensorR2D3)
