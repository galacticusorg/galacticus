"""LibraryInterfaces.ArgSpec — intermediate representation for a single argument.

Andrew Benson (ported to Python with assistance from Claude 2026)
"""

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
    is_array:          bool = False    # 1D numeric array (deferred or fixed shape)
    array_size:        int  = None     # element count for fixed-shape arrays;
                                       # None for deferred (dimension(:))
    array_rank:        int  = 1        # 1 for `dimension(:)` / `dimension(N)`
                                       # / character(len=N) arrays; 2 for
                                       # `dimension(:,:)` 2D numeric.
    array_shape:       tuple = ()      # per-axis element counts for fixed-
                                       # shape multi-rank numeric arrays
                                       # (e.g. (3, 2) for `dimension(3,2)`,
                                       # (2, 2, 3) for `dimension(2,2,3)`).
                                       # Empty for 1D fixed (use array_size)
                                       # and deferred-shape arrays.
    char_len:          int  = 0        # per-element length for fixed-length
                                       # character arrays (`character(len=N),
                                       # dimension(:)`); 0 otherwise.

    # Marker for `type(varying_string), dimension(:)` arrays.  These can't
    # be packed at codegen time (no fixed per-element length), so the
    # Python wrapper computes the maximum encoded length at runtime, pads
    # to that width, and ships it through a flat byte buffer plus two
    # c_size_t companions (count + per-element length).  The Fortran
    # wrapper repacks into `type(varying_string), dimension(:)`.
    varying_string_array: bool = False

    # Marker for `type(<className>List), dimension(:)` arrays — Galacticus's
    # idiom for "array of class(<className>Class)" (Fortran disallows arrays
    # of polymorphic types directly).  Each element holds a `class(...)`
    # pointer to a registered functionClass, so the Python wrapper ships
    # parallel arrays of object pointers and class IDs (one pair per
    # element); the Fortran wrapper rebuilds the list locally via the
    # GetPtr helper, the same machinery the scalar `class(FooClass)` arg
    # path uses.
    polymorphic_list_array:    bool = False
    polymorphic_list_class:    str  = ''   # e.g. 'modelParameter' (registered functionClass)
    polymorphic_list_type:     str  = ''   # e.g. 'modelParameterList' (the wrapper struct)
    polymorphic_list_component: str = ''   # e.g. 'modelParameter_' (component on the wrapper)

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

    # Marker for `class(<intermediate>)` constructor args where
    # <intermediate> is an abstract Fortran type that itself extends a
    # registered functionClass (e.g. `class(massDistributionSpherical)`,
    # whose parent chain reaches `massDistributionClass`).  The wrapper
    # accepts these through the root functionClass's GetPtr and then
    # narrows the polymorphic pointer to the intermediate via a
    # `select type` block emitted by build_fortran_reassignments.  Empty
    # when no narrowing is needed (the plain `class(<base>Class)` case).
    narrowing_type: str = ''

    # Marker for `intent(out), allocatable, dimension(:)` numeric OUTPUT-array
    # arguments — the inner method allocates and fills the array, and its
    # data + size flow *back* to Python (it is not supplied by the caller).
    # Handled like the dynamic-array *return* path: the arg itself is dropped
    # from the bind(c) and Python signatures (fort/py_is_present=False) and
    # replaced by a function-local `save, target` allocatable
    # (`glcOut_<name>_`, named via fort_pass_as) that the inner call fills;
    # the generator appends a `(c_ptr, c_size_t)` intent(out) companion pair
    # to the bind(c) signature and the Python wrapper wraps them into a numpy
    # array (copied, since the save buffer is overwritten on the next call).
    # galacticus_is_present stays True so the inner call still receives it.
    is_output_array:    bool = False
    output_elem_ctype:  str  = ''   # 'c_double' / 'c_int'
    output_elem_fort:   str  = ''   # 'real(c_double)' / 'integer(c_int)'
    output_elem_dtype:  str  = ''   # numpy dtype: 'float64' / 'int32'

    # Marker for a pointer write-back argument: `type(X), pointer,
    # intent(out|inout)` on a METHOD, with X a shared type
    # (_SHARED_TYPE_MODULES — e.g. treeNode).  The bind(c) wrapper takes
    # the handle as `type(c_ptr), intent(inout)` BY REFERENCE, converts a
    # c_associated handle to a local Fortran pointer (null handle →
    # disassociated, the tree-walker start-of-iteration idiom), passes the
    # local to the inner method, and after the call writes c_loc of the
    # (re)pointed target (or c_null_ptr) back through the reference.  The
    # Python caller passes a ctypes.c_void_p, which is updated in place —
    # `while walker.next(node): …` iterates.  Constructor args never take
    # this path (the constructor wrapper has no post-call hook; the
    # classify predicate gates on allow_pointer_writeback).
    is_pointer_writeback: bool = False

    # Marker for a "sized output buffer": an `intent(out)` explicit-shape
    # array whose extents are identifiers naming integer intent(in) args of
    # the same method (`dimension(gridCount,…)`), with numeric or
    # complex(c_double_complex) elements.  The dummy stays in the bind(c)
    # signature (rewritten to `dimension(*)`; the inner call sequence-
    # associates it to the explicit-shape inner dummy); the Python wrapper
    # pre-allocates a flat numpy buffer of the product size (it knows the
    # extents — they're its own parameters), passes its data pointer,
    # drops the arg from the Python signature, and returns the buffer
    # reshaped column-major.  No companions are needed.
    is_output_sized: bool = False
    output_extents:  list = field(default_factory=list)  # extent arg names

    # Marker for an inbound `procedure(<iface>)` callback argument whose
    # abstract interface is registered in
    # Pipeline._CALLBACK_PROCEDURE_INTERFACES.  The bind(c) wrapper receives
    # a C function pointer (`type(c_funptr), value`; ctypes passes a
    # CFUNCTYPE-wrapped Python callable), stores it in a per-method module
    # slot, and hands the module's shim function — which adapts the
    # Galacticus-side arguments and invokes the stored pointer — to the
    # inner method.  All wiring (fort_pass_as/fort_reassignment/
    # fort_modules/py_*) is filled in by the generator, which also emits
    # the storage+shim module; assign_c_types only sets the marker and the
    # c_funptr type mapping.
    is_callback: bool = False

    # Marker for a scalar `intent(out)` numeric/logical companion argument
    # (e.g. `integer, intent(out) :: count` or `double precision,
    # intent(out) :: f`) on an output-array method.  Unlike an output array
    # it stays in the bind(c) signature — as an ordinary by-reference
    # intent(out) scalar (ctype_pointer=True) passed straight to the inner
    # method — but it is dropped from the Python input signature
    # (py_is_present=False) and its filled value is instead appended to the
    # Python return, interleaved with the array outputs in declaration
    # order.  Only recognised on methods that also have an output array (see
    # unsupported_output_array_method); elsewhere scalar intent(out) args
    # keep their existing handling.
    is_output_scalar:   bool = False

    # Marker for constructor args that libraryClasses.xml asks the
    # wrapper to fill with a null pointer rather than expose to Python
    # (`<argument name="..." value="null"/>`).  Used for callback-
    # injection escape-hatch args — typically a `procedure(...), pointer`
    # plus paired `class(*), pointer` slots — that the parameter-driven
    # path of the impl already passes as null.  When set, the arg is
    # dropped from the bind(c) signature (fort_is_present=False) and
    # the Python signature (py_is_present=False); the Fortran wrapper
    # declares a local pointer initialised to null() and passes that to
    # the inner constructor (galacticus_is_present stays True).
    is_null_filled: bool = False

    # Marker for *optional* constructor args that libraryClasses.xml
    # asks the wrapper to drop entirely (`<argument name="..." value="absent"/>`).
    # Used for optional Fortran args whose absence triggers a sensible
    # default in the inner constructor (e.g.
    # `type(vector), dimension(3), optional :: axes` defaulting to the
    # Cartesian basis vectors) but whose declared type isn't otherwise
    # plumbable through the pipeline.  All three is_present flags go
    # False — the arg vanishes from the Python signature, the bind(c)
    # signature, and the inner constructor call.  The optional-arg
    # `make_call(present_set)` mechanism in fortran_call_code already
    # omits args with galacticus_is_present=False, so the inner sees
    # the arg as not-present and falls back to its built-in default.
    is_absent_filled: bool = False

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
