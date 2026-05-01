# LibraryInterfaces.ArgSpec — intermediate representation for a single argument.
# Andrew Benson (ported to Python with assistance from Claude 2026)

from dataclasses import dataclass, field

__all__ = ['ArgSpec']


@dataclass
class ArgSpec:
    """Intermediate representation for a single argument in the cross-language pipeline.

    Instances are built from raw Fortran-declaration dicts by ``assign_c_types()``
    and then progressively enriched by ``assign_c_attributes()``,
    ``build_python_reassignments()``, and ``build_fortran_reassignments()``
    before the emitter functions (``fortran_arg_list``, ``python_call_code``,
    etc.) consume them.

    Field layout mirrors the nested-dict structure it replaces::

        arg['ctypes']['type']          →  arg.ctype
        arg['ctypes']['pointer']       →  arg.ctype_pointer
        arg['fortran']['type']         →  arg.fort_type
        arg['fortran']['isPresent']    →  arg.fort_is_present
        arg['python']['isPresent']     →  arg.py_is_present
        arg['galacticus']['isPresent'] →  arg.galacticus_is_present
        arg['isOptional']              →  arg.is_optional
        arg['isFunctionClass']         →  arg.is_function_class
        arg['passBy']                  →  arg.pass_by
    """

    # -------------------------------------------------------------------------
    # Source-level (from parsed Fortran declaration)
    # -------------------------------------------------------------------------
    name:       str
    intrinsic:  str = ''
    type_spec:  str = ''            # kind/type spec, e.g. 'treeNode', 'varying_string'
    attributes: list = field(default_factory=list)  # e.g. ['intent(in)', 'optional']

    # -------------------------------------------------------------------------
    # Set by assign_c_types
    # -------------------------------------------------------------------------
    is_optional:       bool = False
    is_function_class: bool = False

    # ctypes
    ctype:         str  = ''    # e.g. 'c_double', 'c_void_p', 'c_char_p', 'c_int'
    ctype_pointer: bool = False  # wrap as POINTER(ctype) — set by assign_c_attributes

    # Fortran ABI
    fort_type:           str  = ''    # e.g. 'real(c_double)', 'type(c_ptr)'
    fort_is_present:     bool = True  # include in the bind(C) Fortran argument list
    fort_attributes:     list = field(default_factory=list)  # 'optional', 'value', …
    fort_pass_as:        str  = ''    # expression passed to Galacticus (default: name)
    fort_reassignment:   str  = ''    # Fortran code to run before the call
    fort_declarations:   str  = ''    # extra Fortran local-variable declarations
    fort_iso_c_symbols:  list = field(default_factory=list)  # extra ISO_C_Binding symbols
    fort_modules:        dict = field(default_factory=dict)  # {module: {symbol: 1}}
    fort_function_class: str  = ''    # class name for GetPtr interface blocks

    # Python ctypes
    py_is_present:   bool = True   # include in the Python argument list
    py_pass_as:      str  = ''     # expression to pass to c_lib (default: name)
    py_reassignment: str  = ''     # Python code before the call (optional FC args)
    py_present:      str  = ''     # presence-variable name for optional _ID companions

    # Galacticus constructor call
    galacticus_is_present: bool = True  # include in the Galacticus constructor call

    # -------------------------------------------------------------------------
    # Set by assign_c_attributes
    # -------------------------------------------------------------------------
    pass_by: str = ''    # 'value' or 'reference'

    @classmethod
    def from_raw(cls, d):
        """Build an ArgSpec from a raw Fortran-declaration dict.

        The dict must have a ``'name'`` key; ``'intrinsic'``, ``'type'``, and
        ``'attributes'`` are optional and default to empty values when absent.
        """
        return cls(
            name       = d['name'],
            intrinsic  = d.get('intrinsic') or '',
            type_spec  = d.get('type') or '',
            attributes = list(d.get('attributes') or []),
        )
