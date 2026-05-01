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
import sys
import os

from List.ExtraUtils import as_array
from LibraryInterfaces.ArgSpec import ArgSpec


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

        new_list.insert(0, arg)

    return new_list


def assign_c_attributes(argument_list):
    """Assign C attributes to arguments."""
    for arg in argument_list:
        arg.fort_attributes = []

        if arg.is_optional:
            arg.fort_attributes.append('optional')

        attr_filters = [a for a in arg.attributes
                       if a.startswith('dimension') or a == 'allocatable']
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
        new_list.insert(0, arg)               # unshift current arg
    return new_list


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

        if intrinsic == 'logical':
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
                # Locate the module that imports this enumeration type:
                # 1. walk implementation's module uses, following the extends chain.
                import_module = None
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

            elif type_spec_val == 'treeNode':
                # c_ptr -> type(treeNode) via c_f_pointer.
                arg.fort_declarations = f'type({type_spec_val}), pointer :: {name}_\n'
                arg.fort_pass_as      = name + '_'
                reassign = f'call c_f_pointer({name},{name}_)\n'
                if is_optional:
                    reassign = (f'if (present({name})) then\n '
                                f'{reassign}else\n {name}_ => null()\nend if\n')
                arg.fort_reassignment   = reassign
                arg.fort_iso_c_symbols  = ['c_f_pointer']
                arg.fort_modules.setdefault('Galacticus_Nodes', {})['treeNode'] = 1

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
