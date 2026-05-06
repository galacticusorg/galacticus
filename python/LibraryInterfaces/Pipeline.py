# LibraryInterfaces.Pipeline — four pipeline stages that enrich ArgSpec objects.
# Andrew Benson (ported to Python with assistance from Claude 2026)
#
# Each stage accepts a list of ArgSpec (or raw dicts for assign_c_types) and
# returns the enriched list.  Stages must be applied in order:
#
#   1. assign_c_types(raw_list, lib_function_classes)
#      Converts raw Fortran-declaration dicts to ArgSpec objects and sets
#      ctype / fort_type / is_optional / is_function_class.  Inserts _ID
#      companion args for functionClass parameters.
#
#   2. assign_c_attributes(arg_list)
#      Sets fort_attributes, pass_by, ctype_pointer based on the type and
#      direction of each argument.
#
#   3. build_python_reassignments(arg_list)
#      Populates py_pass_as and py_reassignment for functionClass arguments.
#
#   4. build_fortran_reassignments(arg_list, func_class, implementation,
#                                  extensions, module_uses_impls,
#                                  lib_function_classes=None)
#      Populates fort_pass_as, fort_reassignment, fort_declarations,
#      fort_modules, and fort_iso_c_symbols for all arguments that need
#      cross-language type conversion.

import re

from List.ExtraUtils import as_array
from LibraryInterfaces.ArgSpec import ArgSpec

__all__ = [
    'assign_c_types', 'assign_c_attributes',
    'build_python_reassignments', 'build_fortran_reassignments',
]


# Shared derived types whose defining module is NOT the functionClass's own
# module.  Without these explicit overrides, build_fortran_reassignments'
# fall-back imports the type from func_class['module'] — which compiles for
# types defined alongside the class but emits broken `use ::` lines for
# shared types like treeNode, mergerTree, multiCounter, etc.
#
# This table is consulted by every site that has to decide which module a
# referenced type lives in: the constructor-arg fall-back below, the
# enumeration-arg lookup (where it short-circuits the moduleUses walk),
# and `_find_enum_module` in libraryInterfaces.py for method return types.
# An explicit entry always wins over the smart-lookup fall-backs.
#
# Add to this table when a new shared type starts appearing as a method or
# constructor argument in libraryClasses.xml.
_SHARED_TYPE_MODULES = {
    'treeNode'                            : 'Galacticus_Nodes',
    'mergerTree'                          : 'Galacticus_Nodes',
    'universe'                            : 'Galacticus_Nodes',
    'multiCounter'                        : 'Multi_Counters',
    'history'                             : 'Histories',
    'keplerOrbit'                         : 'Kepler_Orbits',
    'hdf5Object'                          : 'IO_HDF5',
    'abundances'                          : 'Abundances_Structure',
    'chemicalAbundances'                  : 'Chemical_Abundances_Structure',
    'multiExtractorList'                  : 'Node_Property_Extractors',
    'enumerationFrameType'                : 'Stellar_Luminosities_Structure',
    'enumerationDestroyStubsType'         : 'Merger_Tree_Build_Controllers',
    'enumerationComponentTypeType'        : 'Galactic_Structure_Options',
    'enumerationCoolingFromType'          : 'Cooling_Options',
    'enumerationIntervalTypeType'         : 'Nodes_Operators',
    'enumerationParticulateKernelType'    : 'Merger_Tree_Operators',
    'enumerationRandomSampleCountTypeType': 'Tasks',
    'enumerationRelativeToType'           : 'Nodes_Operators',
    'enumerationSelectionType'            : 'Merger_Tree_Operators',
}


# Match a fixed-size 1D dimension attribute (e.g. dimension(3)).  Used to
# distinguish a fixed-size array — whose length is known at codegen time
# and so needs no count companion — from a deferred-shape dimension(:),
# which does.
_DIM_FIXED_RX = re.compile(r'^dimension\s*\(\s*(\d+)\s*\)$')


def assign_c_types(argument_list, lib_function_classes):
    """Assign appropriate C types for each argument.

    Mirrors Perl assignCTypes().  Accepts a list of raw Fortran-declaration
    dicts and returns a new list of :class:`ArgSpec` objects.  Processes in
    reverse so that the ``_ID`` companion argument for functionClass parameters
    can be inserted immediately after its parent without disturbing the rest
    of the list.
    """
    new_list = []
    for raw in reversed(argument_list):
        arg = ArgSpec.from_raw(raw) if isinstance(raw, dict) else raw

        # Initialise presence flags (mirrors the original explicit reset).
        arg.fort_is_present       = True
        arg.py_is_present         = True
        arg.galacticus_is_present = True
        arg.is_function_class     = False
        arg.is_optional           = bool('optional' in arg.attributes)

        intrinsic     = arg.intrinsic
        type_spec_val = arg.type_spec

        if intrinsic == 'double precision':
            arg.ctype     = 'c_double'
            arg.fort_type = 'real(c_double)'
        elif intrinsic == 'integer':
            # Default kind maps to c_int; explicit C-interop kinds (c_long,
            # c_size_t) pass through with matching ctypes wrappers so that
            # 64-bit values aren't silently truncated.  `kind_int8` is the
            # Galacticus alias for `selected_int_kind(18)` (a 64-bit
            # integer) — without this branch its arrays would be emitted
            # as `integer(c_int)` and mismatch the inner method's
            # `integer(kind_int8)` signature, breaking the build.
            if type_spec_val in ('c_long', 'kind_int8'):
                arg.ctype     = 'c_long'
                arg.fort_type = 'integer(c_long)'
            elif type_spec_val == 'c_size_t':
                arg.ctype     = 'c_size_t'
                arg.fort_type = 'integer(c_size_t)'
            else:
                arg.ctype     = 'c_int'
                arg.fort_type = 'integer(c_int)'
        elif intrinsic == 'logical':
            arg.ctype     = 'c_bool'
            arg.fort_type = 'logical(c_bool)'
        elif intrinsic == 'character':
            arg.ctype     = 'c_char_p'
            arg.fort_type = 'character(c_char)'
        elif intrinsic == 'type':
            if type_spec_val == 'varying_string':
                arg.ctype     = 'c_char_p'
                arg.fort_type = 'character(c_char)'
            elif re.match(r'^enumeration[a-z0-9_]+type$', type_spec_val, re.IGNORECASE):
                # Enumeration types map to C int.
                arg.ctype     = 'c_int'
                arg.fort_type = 'integer(c_int)'
            else:
                arg.ctype     = 'c_void_p'
                arg.fort_type = 'type(c_ptr)'
        elif intrinsic == 'class':
            arg.ctype     = 'c_void_p'
            arg.fort_type = 'type(c_ptr)'
            # Check whether this is a functionClass argument.
            if type_spec_val.endswith('Class'):
                class_key = type_spec_val[:-5]  # strip trailing 'Class'
                if class_key in lib_function_classes:
                    arg.is_function_class = True
                    # 'self' is dispatched via the method binding, not passed directly.
                    if arg.name == 'self':
                        arg.galacticus_is_present = False
                    # Insert a companion _ID argument (carries the concrete class ID).
                    arg_id = ArgSpec(
                        name                  = arg.name + '_ID',
                        intrinsic             = 'integer',
                        attributes            = ['intent(in)'],
                        ctype                 = 'c_int',
                        fort_type             = 'integer(c_int)',
                        fort_is_present       = True,
                        py_is_present         = False,
                        galacticus_is_present = False,
                        is_optional           = False,
                        is_function_class     = False,
                    )
                    if arg.is_optional:
                        arg_id.attributes.append('optional')
                        arg_id.is_optional = True
                        arg_id.py_present  = arg.name
                    # Insert _ID before the current front, then arg before that.
                    new_list.insert(0, arg_id)

        # Detect 1D numeric arrays (deferred-shape or fixed-size) and set
        # is_array / array_size so the rest of the pipeline can recognise
        # them.  Deferred-shape gets a hidden integer(c_size_t) count
        # companion immediately after it — same trick as the _ID
        # companion for class(FooClass) args; the Python wrapper computes
        # the value from the input numpy array's .size, and the inner
        # Galacticus method receives an `arr(1:arr_count)` slice (see
        # build_fortran_reassignments).  Fixed-size needs no companion:
        # the length is in the dimension spec itself.
        if intrinsic in ('double precision', 'integer'):
            if 'dimension(:)' in arg.attributes:
                arg.is_array   = True
                arg.array_size = None
                count_arg = ArgSpec(
                    name                  = arg.name + '_count',
                    intrinsic             = 'integer',
                    type_spec             = 'c_size_t',
                    attributes            = ['intent(in)'],
                    ctype                 = 'c_size_t',
                    fort_type             = 'integer(c_size_t)',
                    fort_is_present       = True,
                    py_is_present         = False,
                    galacticus_is_present = False,
                    is_optional           = False,
                    is_function_class     = False,
                )
                new_list.insert(0, count_arg)
            else:
                for a in arg.attributes:
                    m = _DIM_FIXED_RX.match(a)
                    if m:
                        arg.is_array   = True
                        arg.array_size = int(m.group(1))
                        break

        new_list.insert(0, arg)

    return new_list


def assign_c_attributes(argument_list):
    """Assign C attributes to arguments."""
    for arg in argument_list:
        arg.fort_attributes = []

        if arg.is_optional:
            arg.fort_attributes.append('optional')

        # dimension(:) (assumed-shape) isn't allowed in bind(c); rewrite to
        # dimension(*) (assumed-size) so the C side sees a plain pointer.
        # The companion <name>_count carries the runtime length, and the
        # inner Galacticus call slices arr(1:arr_count) to recover an
        # ordinary dimension(:) section — see build_fortran_reassignments.
        attr_filters = []
        for a in arg.attributes:
            if a == 'dimension(:)':
                attr_filters.append('dimension(*)')
            elif a.startswith('dimension') or a == 'allocatable':
                attr_filters.append(a)
        arg.fort_attributes.extend(attr_filters)

        if arg.ctype == 'c_char_p':
            arg.fort_attributes.append('dimension(*)')

        # Determine pass-by method
        is_ptr_type   = arg.ctype.endswith('_p')
        is_intent_in  = any('intent(in)' in a for a in arg.attributes)
        is_non_scalar = any(a.startswith('dimension') for a in arg.fort_attributes)

        if (is_ptr_type or is_intent_in) and not arg.is_optional and not is_non_scalar:
            arg.pass_by = 'value'
        else:
            arg.pass_by = 'reference'

        if arg.pass_by == 'value':
            arg.fort_attributes.append('value')

        arg.ctype_pointer = (arg.pass_by == 'reference' and arg.ctype != 'c_char_p')

    return argument_list


def build_python_reassignments(argument_list):
    """Set py_pass_as and py_reassignment for functionClass args.

    Mirrors Perl buildPythonReassignments().  Processes in reverse so that when
    a functionClass arg is encountered, its _ID companion is already sitting at
    the front of new_list (it was the immediately preceding arg in forward
    order, so the last one pushed in reverse order).

    Non-optional:  py_pass_as = 'name._glcObj' / 'name._classID'
    Optional:      py_pass_as = 'name_glcObj'  / 'name_classID' plus a
                   py_reassignment block that extracts the values or sets them
                   to None when the argument is absent.
    """
    new_list = []
    for arg in reversed(argument_list):
        if arg.is_function_class:
            arg_id = new_list.pop(0)          # shift _ID off front of new list
            name = arg.name
            if arg.is_optional:
                arg.py_pass_as       = name + '_glcObj'
                arg_id.py_pass_as    = name + '_classID'
                arg.py_reassignment  = (
                    f'    if {name}:\n'
                    f'        {name}_glcObj ={name}._glcObj\n'
                    f'        {name}_classID={name}._classID\n'
                    f'    else:\n'
                    f'        {name}_glcObj =None\n'
                    f'        {name}_classID=None\n'
                )
            else:
                arg.py_pass_as    = name + '._glcObj'
                arg_id.py_pass_as = name + '._classID'
            new_list.insert(0, arg_id)        # unshift _ID back
        elif arg.is_array:
            # 1D numeric array: convert the user's input (numpy array, list,
            # tuple) to a contiguous ndarray of the right dtype, then pass
            # its data pointer to the bind(c) function.  Deferred-shape
            # also passes a count: same pop/unshift trick as the _ID case
            # above to grab the count companion that's sitting at the
            # front of new_list.  Fixed-size validates the length so a
            # mismatched input raises a clear ValueError instead of
            # silently corrupting Fortran-side memory.
            np_dtype  = _ARRAY_NUMPY_DTYPE.get(arg.ctype, 'float64')
            if arg.array_size is None:
                count_arg = new_list.pop(0)
                arg.py_reassignment = (
                    f'    {arg.name} = np.ascontiguousarray({arg.name},'
                    f' dtype=np.{np_dtype})\n'
                )
                count_arg.py_pass_as = f'c_size_t({arg.name}.size)'
                new_list.insert(0, count_arg)
            else:
                arg.py_reassignment = (
                    f'    {arg.name} = np.ascontiguousarray({arg.name},'
                    f' dtype=np.{np_dtype})\n'
                    f'    if {arg.name}.size != {arg.array_size}:\n'
                    f'        raise ValueError('
                    f'f"{arg.name} expects {arg.array_size} elements, got '
                    f'{{{arg.name}.size}}")\n'
                )
            arg.py_pass_as = (
                f'{arg.name}.ctypes.data_as(POINTER({arg.ctype}))'
            )
        new_list.insert(0, arg)               # unshift current arg
    return new_list


# Numpy dtype string used when converting a Python input to a contiguous
# array for an array-typed argument.  Keep in sync with the ctypes mapping
# in assign_c_types.
_ARRAY_NUMPY_DTYPE = {
    'c_double': 'float64',
    'c_int'   : 'int32',
    'c_long'  : 'int64',
    'c_size_t': 'uint64',
}


def build_fortran_reassignments(argument_list, func_class, implementation,
                                extensions, module_uses_impls,
                                lib_function_classes=None):
    """Generate Fortran reassignments for cross-language type conversions.

    Mirrors Perl buildFortranReassignments().  Processes in reverse (same
    pop/unshift skeleton as the other builders); no new args are inserted so
    the order is unchanged.

    Cases handled:
      logical        — c_bool   → logical via logical() cast
      character      — c_char_p → varying_string via char(String_C_to_Fortran())
      varying_string — c_char_p → varying_string via String_C_to_Fortran()
      enumeration    — c_int    → type(enumXxx) via %ID assignment; module located
                                  by walking implementation/extension/functionClass uses
      treeNode       — c_ptr    → type(treeNode) via c_f_pointer
      other types    — c_ptr    → type(X) via c_f_pointer
      isFunctionClass — c_ptr   → class(XClass) via XGetPtr(ptr, ID)
      class(*)       — c_ptr    → concrete type via c_f_pointer (type from libraryClasses config)
    """
    if lib_function_classes is None:
        lib_function_classes = {}

    new_list = []
    for arg in reversed(argument_list):
        intrinsic     = arg.intrinsic
        type_spec_val = arg.type_spec
        name          = arg.name
        is_optional   = arg.is_optional
        opt_prefix    = f'if (present({name})) ' if is_optional else ''

        if arg.is_array:
            # 1D numeric array.  Deferred-shape needs slicing — the bind(c)
            # signature has `dimension(*)` (assumed-size) plus a separate
            # <name>_count, so we hand the inner method a section
            # `name(1:name_count)` to recover a proper rank-1 descriptor of
            # the right length.  Fixed-size doesn't need slicing: the
            # bind(c) declaration is already explicit-shape `dimension(N)`,
            # which the inner method accepts directly.
            if arg.array_size is None:
                arg.fort_pass_as = f'{name}(1:{name}_count)'
        elif intrinsic == 'logical':
            # c_bool must be recast to a plain Fortran logical.
            arg.fort_reassignment = f'{opt_prefix}{name}_=logical({name})\n'
            arg.fort_declarations = f'logical :: {name}_\n'
            arg.fort_pass_as      = name + '_'

        elif intrinsic == 'character':
            # c_char_p -> character via String_C_to_Fortran then char()
            arg.fort_modules.setdefault('String_Handling',    {})['String_C_to_Fortran'] = 1
            arg.fort_modules.setdefault('ISO_Varying_String', {})['char']                = 1
            arg.fort_pass_as = f'char(String_C_to_Fortran({name}))'

        elif intrinsic == 'type':
            if type_spec_val == 'varying_string':
                # c_char_p -> varying_string via String_C_to_Fortran
                arg.fort_modules.setdefault('String_Handling', {})['String_C_to_Fortran'] = 1
                arg.fort_pass_as = f'String_C_to_Fortran({name})'

            elif re.match(r'^enumeration[a-z0-9_]+type$', type_spec_val, re.IGNORECASE):
                # c_int -> type(enumXxx) via %ID assignment.
                arg.fort_declarations = f'type({type_spec_val}) :: {name}_\n'
                arg.fort_pass_as      = name + '_'
                arg.fort_reassignment = f'{opt_prefix}{name}_%ID={name}\n'
                # Locate the module that imports this enumeration type.
                # 0. an explicit override in _SHARED_TYPE_MODULES wins outright
                #    — used when the impl file's own moduleUses don't carry the
                #    type explicitly (e.g. when the enum is defined in the same
                #    module the impl is included into) so the walk below would
                #    otherwise miss it.
                import_module = _SHARED_TYPE_MODULES.get(type_spec_val)
                # 1. walk implementation's module uses, following the extends chain.
                if implementation:
                    cls = implementation['name']
                    while cls and not import_module:
                        for use_block in module_uses_impls.get(cls, []):
                            for mod_name, mod_data in use_block.items():
                                if (isinstance(mod_data, dict)
                                        and type_spec_val in mod_data.get('only', {})):
                                    import_module = mod_name
                                    break
                            if import_module:
                                break
                        cls = extensions.get(cls)
                # 2. fall back to the functionClass file's own module uses.
                if not import_module:
                    for use_block in func_class.get('moduleUses', []):
                        for mod_name, mod_data in use_block.items():
                            if (isinstance(mod_data, dict)
                                    and type_spec_val in mod_data.get('only', {})):
                                import_module = mod_name
                                break
                        if import_module:
                            break
                # 3. last resort: the functionClass's own module.
                if not import_module:
                    import_module = func_class.get('module')
                if import_module:
                    arg.fort_modules.setdefault(import_module, {})[type_spec_val] = 1

            elif type_spec_val in _SHARED_TYPE_MODULES:
                # c_ptr -> type(X) via c_f_pointer for derived types that
                # live in a known module other than the functionClass's own
                # (see _SHARED_TYPE_MODULES).  Includes the proper
                # null-on-absent handling for optional args.
                arg.fort_declarations = f'type({type_spec_val}), pointer :: {name}_\n'
                arg.fort_pass_as      = name + '_'
                reassign = f'call c_f_pointer({name},{name}_)\n'
                if is_optional:
                    reassign = (f'if (present({name})) then\n '
                                f'{reassign}else\n {name}_ => null()\nend if\n')
                arg.fort_reassignment   = reassign
                arg.fort_iso_c_symbols  = ['c_f_pointer']
                shared_mod = _SHARED_TYPE_MODULES[type_spec_val]
                arg.fort_modules.setdefault(shared_mod, {})[type_spec_val] = 1

            else:
                # Other derived types: c_ptr -> type(X) via c_f_pointer.
                arg.fort_declarations  = f'type({type_spec_val}), pointer :: {name}_\n'
                arg.fort_pass_as       = name + '_'
                arg.fort_reassignment  = f'{opt_prefix}call c_f_pointer({name},{name}_)\n'
                arg.fort_iso_c_symbols = ['c_f_pointer']
                if type_spec_val == 'inputParameters':
                    arg.fort_modules.setdefault('Input_Parameters', {})['inputParameters'] = 1
                else:
                    mod = func_class.get('module', '')
                    if mod:
                        arg.fort_modules.setdefault(mod, {})[type_spec_val] = 1

        elif arg.is_function_class:
            # class(FooClass) -> class(FooClass) pointer via FooGetPtr(ptr, ID).
            class_key = type_spec_val[:-5] if type_spec_val.endswith('Class') else type_spec_val
            mod_name  = lib_function_classes.get(class_key, {}).get('module', '')
            arg.fort_declarations  = f'class({type_spec_val}), pointer :: {name}_\n'
            arg.fort_pass_as       = name + '_'
            arg.fort_reassignment  = (
                f'{opt_prefix}{name}_ => {class_key}GetPtr({name},{name}_ID)\n'
            )
            arg.fort_function_class = class_key
            if mod_name:
                arg.fort_modules.setdefault(mod_name, {})[type_spec_val] = 1

        elif intrinsic == 'class' and type_spec_val == '*':
            # Unlimited polymorphic: look up the concrete type in libraryClasses config.
            impl_name = implementation['name'] if implementation else None
            fc_name   = func_class.get('name', '')
            impl_info = (lib_function_classes.get(fc_name, {}).get(impl_name, {})
                         if impl_name else {})
            ctor_args = as_array(impl_info.get('constructor', {}).get('argument', []))
            concrete  = [a for a in ctor_args
                         if isinstance(a, dict) and a.get('name') == name]
            if concrete:
                ct      = concrete[0]
                ct_type = ct.get('type', '')
                ct_mod  = ct.get('module', '')
                arg.fort_declarations  = f'type({ct_type}), pointer :: {name}_\n'
                arg.fort_pass_as       = name + '_'
                arg.fort_reassignment  = f'{opt_prefix}call c_f_pointer({name},{name}_)\n'
                arg.fort_iso_c_symbols = ['c_f_pointer']
                if ct_mod and ct_type:
                    arg.fort_modules.setdefault(ct_mod, {})[ct_type] = 1

        new_list.insert(0, arg)

    return new_list
