"""LibraryInterfaces.Emitters — ten code-string emitter functions.

Andrew Benson (ported to Python with assistance from Claude 2026)

All functions accept a list of fully-enriched ArgSpec objects (produced by
the four Pipeline stages) and return a code string or a list of strings.
"""

import keyword
import re
from itertools import chain, combinations

__all__ = [
    'ctypes_arg_types', 'fortran_arg_list', 'fortran_declarations',
    'fortran_reassignments', 'fortran_module_uses', 'fortran_call_code',
    'iso_c_binding_import',
    'python_arg_list', 'python_reassignments', 'python_call_code',
    'python_safe_name',
]


def python_safe_name(name):
    """Return *name* escaped so it can be used as a Python identifier.

    Galacticus argument and method names (`lambda`, `yield`, `class`, ...)
    can collide with Python reserved words, which would emit a syntax
    error in the generated wrapper module.  Follow PEP 8's
    trailing-underscore convention to escape; the Fortran identifier is
    unaffected because every emitter that touches the Fortran side uses
    `arg.name` (or `fort_pass_as`), not this helper.
    """
    return name + '_' if keyword.iskeyword(name) else name


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _powerset(iterable):
    """Return all subsets of *iterable* in ascending-size order."""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


# ---------------------------------------------------------------------------
# ctypes / shared-library interface emitters
# ---------------------------------------------------------------------------

def ctypes_arg_types(argument_list):
    """Generate ctypes argument type list.

    Returns a list of strings such as ``['c_double', 'POINTER(c_int)', …]``
    suitable for assignment to ``c_lib.funcName.argtypes``.
    """
    types = []
    for arg in argument_list:
        ctype = arg.ctype or 'c_int'
        if arg.ctype_pointer:
            ctype = f'POINTER({ctype})'
        types.append(ctype)
    return types


# ---------------------------------------------------------------------------
# Fortran emitters
# ---------------------------------------------------------------------------

def fortran_arg_list(argument_list):
    """Return the list of argument names for the bind(C) Fortran wrapper."""
    return [arg.name for arg in argument_list if arg.fort_is_present]


def fortran_declarations(argument_list):
    """Generate Fortran declarations.

    Mirrors Perl fortranDeclarations(): emits one declaration per argument
    plus any extra declarations stored in ``fort_declarations``.  Also
    emits an interface block for every distinct functionClass GetPtr function
    that is referenced by an is_function_class argument.
    """
    code = ''
    function_classes = {}   # {className: True} — deduplicates interface blocks

    for arg in argument_list:
        # Skip args that don't appear in the bind(c) signature
        # (fort_is_present=False).  These are typically null-filled
        # constructor overrides whose local-pointer declaration lives
        # entirely in `fort_declarations`; emitting the default
        # `integer(c_int)` line for them would shadow that local and
        # produce conflicting declarations of the same name.
        if not arg.fort_is_present:
            if arg.fort_declarations:
                code += arg.fort_declarations
            if arg.fort_function_class:
                function_classes[arg.fort_function_class] = True
            continue
        attr_str  = (', ' + ', '.join(arg.fort_attributes)) if arg.fort_attributes else ''
        fort_type = arg.fort_type or 'integer(c_int)'
        code += f'  {fort_type}{attr_str} :: {arg.name}\n'
        if arg.fort_declarations:
            code += arg.fort_declarations
        if arg.fort_function_class:
            function_classes[arg.fort_function_class] = True

    # Emit interface blocks so the compiler knows the GetPtr signatures.
    for fc in sorted(function_classes):
        code += (
            f'interface\n'
            f' function {fc}GetPtr(ptr_,classID)\n'
            f'  import c_int, c_ptr, {fc}Class\n'
            f'  class({fc}Class), pointer :: {fc}GetPtr\n'
            f'  type   (c_ptr), intent(in   ) :: ptr_\n'
            f'  integer(c_int), intent(in   ) :: classID\n'
            f' end function {fc}GetPtr\n'
            f'end interface\n'
        )
    return code


def fortran_reassignments(argument_list):
    """Concatenate all Fortran pre-call reassignment statements."""
    return ''.join(arg.fort_reassignment for arg in argument_list)


def fortran_module_uses(argument_list):
    """Generate Fortran module use statements.

    Mirrors Perl fortranModuleUses(): accumulates ``{module: {symbol: 1}}``
    dicts from every arg's ``fort_modules`` field, then emits one ``use``
    line per module with a sorted ``only`` list.
    """
    modules = {}
    for arg in argument_list:
        for mod_name, symbols in arg.fort_modules.items():
            modules.setdefault(mod_name, {}).update(symbols)

    code = ''
    for mod_name in sorted(modules):
        syms = ', '.join(sorted(modules[mod_name]))
        code += f'  use :: {mod_name}, only : {syms}\n'
    return code


def fortran_call_code(argument_list, pre_arguments, post_arguments, continuation):
    """Generate Fortran call code, with optional-argument branching.

    Mirrors Perl fortranCallCode().  When N optional args are present, emits
    2^N if/else-if branches using the ``present()`` intrinsic, plus an
    unconditional else fallback to suppress compiler warnings about unset
    function results.
    """
    def make_call(present_set):
        # present_set=None → all galacticus-present args included;
        # otherwise only include optionals whose name is in present_set.
        args = []
        for arg in argument_list:
            if not arg.galacticus_is_present:
                continue
            if (arg.is_optional and present_set is not None
                    and arg.name not in present_set):
                continue
            pass_name = arg.fort_pass_as if arg.fort_pass_as else arg.name
            args.append(f"{arg.name}={pass_name}")
        sep = ',' + continuation + '\n' + '  ' + continuation + ' '
        return (f"{pre_arguments}{continuation} {sep.join(args)} {continuation}\n"
                f"{post_arguments}")

    optional_names = sorted(
        arg.name for arg in argument_list
        if arg.is_optional and arg.galacticus_is_present
    )

    if not optional_names:
        return make_call(None)

    code = ''
    first = True
    for subset in _powerset(optional_names):
        present_set = set(subset)
        conditions = [
            (f'present({n})' if n in present_set else f'.not.present({n})')
            for n in optional_names
        ]
        prefix = '' if first else 'else '
        code += f'{prefix}if ({" .and. ".join(conditions)}) then\n'
        code += make_call(present_set)
        first = False
    # Unconditional else avoids compiler warnings about unset function results.
    code += 'else\n'
    code += make_call(None)
    code += 'end if\n'
    return code


def iso_c_binding_import(argument_list, *extra_symbols):
    """Generate ISO_C_Binding import statement.

    Mirrors Perl isoCBindingImport(): collects the kind symbol from each
    argument's Fortran type (e.g. ``'c_double'`` from ``'real(c_double)'``)
    plus any extra symbols stored in ``fort_iso_c_symbols`` (e.g.
    ``'c_f_pointer'`` added for pointer-dereference reassignments).
    """
    symbols = set(extra_symbols)
    for arg in argument_list:
        m = re.search(r'\(([a-z_]+)\)', arg.fort_type)
        if m:
            symbols.add(m.group(1))
        for sym in arg.fort_iso_c_symbols:
            symbols.add(sym)
    return f"  use, intrinsic :: ISO_C_Binding, only : {', '.join(sorted(symbols))}\n"


# ---------------------------------------------------------------------------
# Python emitters
# ---------------------------------------------------------------------------

def python_arg_list(argument_list):
    """Generate Python argument list.

    Mirrors Perl pythonArgList(): if the first argument is not named ``'self'``
    (i.e. for constructors) prepend ``'self'`` explicitly.  For methods the
    first argument already is ``'self'`` (py_is_present=True) so the loop
    adds it and no explicit prepend is needed.
    """
    first_name = argument_list[0].name if argument_list else None
    args = [] if first_name == 'self' else ['self']

    first_optional = False
    for arg in argument_list:
        if not arg.py_is_present:
            continue
        name = python_safe_name(arg.name)
        if arg.is_optional or first_optional:
            name += '=None'
            first_optional = True
        args.append(name)
    return args


def python_reassignments(argument_list):
    """Concatenate all Python pre-call reassignment statements."""
    return ''.join(arg.py_reassignment for arg in argument_list)


def python_call_code(argument_list, call):
    """Generate Python call code, with optional-argument branching.

    Mirrors Perl pythonCallCode().  When N optional args are present, emits
    2^N if/elif branches.  Once the first optional arg is encountered, ALL
    subsequent args need explicit ctype wrapping (ctypes cannot infer types
    when some args may be None).  Absent optional args are passed as None.
    """
    def make_py_call(present_set, indent='        '):
        # present_set=None → all present; else set of presence-variable names.
        args = []
        first_optional = False
        for arg in argument_list:
            if not arg.fort_is_present:
                continue
            is_opt = arg.is_optional
            pv     = arg.py_present if arg.py_present else python_safe_name(arg.name)
            pa     = arg.py_pass_as if arg.py_pass_as else python_safe_name(arg.name)
            if is_opt:
                first_optional = True
            if first_optional:
                if is_opt and present_set is not None and pv not in present_set:
                    args.append('None')
                elif arg.is_array:
                    # Array args are always passed as pointers
                    # (`arr.ctypes.data_as(POINTER(...))`) — wrapping them
                    # in `{ctype}(...)` would call e.g. c_double(<pointer>)
                    # which ctypes rejects.  Pass `pa` directly; for
                    # optional arrays `pa` is itself a None-aware
                    # conditional expression set up by
                    # build_python_reassignments.
                    args.append(pa)
                else:
                    ctype = arg.ctype or 'c_void_p'
                    args.append(f'{ctype}({pa})')
            else:
                args.append(pa)
        return f"{indent}{call}({','.join(args)})\n"

    # Collect unique presence-variable names for optional fortran-present args.
    pv_seen = {}
    for arg in argument_list:
        if not arg.is_optional:
            continue
        if not arg.fort_is_present:
            continue
        pv = arg.py_present if arg.py_present else python_safe_name(arg.name)
        pv_seen[pv] = True
    optional_pvs = sorted(pv_seen.keys())

    if not optional_pvs:
        return make_py_call(None, indent='    ')

    code = ''
    first = True
    for subset in _powerset(optional_pvs):
        present_set = set(subset)
        conditions = [
            (f'{pv} is not None' if pv in present_set else f'{pv} is None')
            for pv in optional_pvs
        ]
        prefix = '' if first else 'el'
        code += f"    {prefix}if {' and '.join(conditions)}:\n"
        code += make_py_call(present_set)
        first = False
    return code
