"""LibraryInterfaces.Pipeline — four pipeline stages that enrich ArgSpec objects.

Andrew Benson (ported to Python with assistance from Claude 2026)

Each stage accepts a list of ArgSpec (or raw dicts for assign_c_types) and
returns the enriched list.  Stages must be applied in order:

  1. assign_c_types(raw_list, lib_function_classes)
     Converts raw Fortran-declaration dicts to ArgSpec objects and sets
     ctype / fort_type / is_optional / is_function_class.  Inserts _ID
     companion args for functionClass parameters.

  2. assign_c_attributes(arg_list)
     Sets fort_attributes, pass_by, ctype_pointer based on the type and
     direction of each argument.

  3. build_python_reassignments(arg_list)
     Populates py_pass_as and py_reassignment for functionClass arguments.

  4. build_fortran_reassignments(arg_list, func_class, implementation,
                                 extensions, module_uses_impls,
                                 lib_function_classes=None)
     Populates fort_pass_as, fort_reassignment, fort_declarations,
     fort_modules, and fort_iso_c_symbols for all arguments that need
     cross-language type conversion.
"""

import re

from List.ExtraUtils import as_array
from LibraryInterfaces.ArgSpec import ArgSpec
from LibraryInterfaces.Emitters import python_safe_name
from LibraryInterfaces.Hierarchy import resolve_function_class_base

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
    'treeNode'                                    : 'Galacticus_Nodes',
    'mergerTree'                                  : 'Galacticus_Nodes',
    'universe'                                    : 'Galacticus_Nodes',
    'multiCounter'                                : 'Multi_Counters',
    'history'                                     : 'Histories',
    'keplerOrbit'                                 : 'Kepler_Orbits',
    'hdf5Object'                                  : 'IO_HDF5',
    'abundances'                                  : 'Abundances_Structure',
    'chemicalAbundances'                          : 'Chemical_Abundances_Structure',
    'multiExtractorList'                          : 'Node_Property_Extractors',
    'modelParameterList'                          : 'Model_Parameters',
    'stellarPopulationSpectraPostprocessorList'   : 'Stellar_Population_Spectra_Postprocess',
    'enumerationFrameType'                        : 'Stellar_Luminosities_Structure',
    'enumerationDestroyStubsType'                 : 'Merger_Tree_Build_Controllers',
    'enumerationComponentTypeType'                : 'Galactic_Structure_Options',
    'enumerationCoolingFromType'                  : 'Cooling_Options',
    'enumerationIntervalTypeType'                 : 'Nodes_Operators',
    'enumerationParticulateKernelType'            : 'Merger_Tree_Operators',
    'enumerationRandomSampleCountTypeType'        : 'Tasks',
    'enumerationRelativeToType'                   : 'Nodes_Operators',
    'enumerationSelectionType'                    : 'Merger_Tree_Operators',
    'enumerationOutputAnalysisCovarianceModelType': 'Output_Analyses_Options',
}

# Inbound procedure (callback) arguments the pipeline can marshal.  Keyed by
# the abstract-interface name appearing in `procedure(<name>) :: arg`; each
# entry describes how to bridge a Python callable to the Fortran interface:
#
#   cfunctype  — the ctypes callback factory for the C-side signature.
#   py_adapter — format template (placeholder {fn} = the user's Python
#                callable) producing the function handed to `cfunctype`;
#                adapts raw ctypes pointers to friendly numpy/scalar views.
#   c_iface    — format template ({iface}) for the Fortran abstract
#                interface matching the C-side signature.
#   shim       — format template ({shim}, {iface}, {var}) for a Fortran
#                function matching the *Galacticus* interface; it converts
#                the Fortran-side arguments to C-representable values,
#                retrieves the stored Python callback via c_f_procpointer,
#                and invokes it.
#   shim_uses  — {module: (symbols,…)} needed by the shim (beyond
#                ISO_C_Binding, which is always imported).
#
# The generator emits, per callback-taking method, a small module holding a
# `type(c_funptr)` slot ({var}) plus the shim; the bind(c) wrapper stores
# the incoming funptr in the slot and passes the shim to the inner method.
# Like the save-buffer array returns, the slot is process-global state: the
# callback is valid for the duration of the wrapped call, and concurrent
# calls to the same method from multiple threads would race — matching the
# existing library-interface conventions.
#
# Registry gating (rather than parsing arbitrary abstract interfaces from
# source) keeps this explicit: an entry is a promise that the C signature,
# the adapter, and the shim's argument conversion have been checked by hand.
# Consulted by Classification.classify_arg (membership ⇒ supported), the
# inline accept in assign_c_types below, and the emission code in
# libraryInterfaces.interfaces_methods.
_CALLBACK_PROCEDURE_INTERFACES = {
    # double precision function crossSectionFunction(wavelength)
    #   double precision, intent(in) :: wavelength
    # (source/radiation/_class.F90).  The method dummy is
    # `procedure(...), pointer`, so the wrapper passes the module's
    # procedure-pointer slot aimed at the shim (a pointer actual for a
    # pointer dummy) — the generator keys that off the arg's own
    # `pointer` attribute.
    'crossSectionFunctionTemplate': {
        'cfunctype' : 'CFUNCTYPE(c_double, c_double)',
        'py_adapter': '(lambda _wavelength_: float({fn}(_wavelength_)))',
        'c_iface'   : (
            '  abstract interface\n'
            '     function {iface}(wavelength) bind(c)\n'
            '       import c_double\n'
            '       real(c_double), value :: wavelength\n'
            '       real(c_double)        :: {iface}\n'
            '     end function {iface}\n'
            '  end interface\n'
        ),
        'shim'      : (
            '  double precision function {shim}(wavelength)\n'
            '    double precision, intent(in   ) :: wavelength\n'
            '    procedure({iface}), pointer     :: fn_\n'
            '    call c_f_procpointer({var}, fn_)\n'
            '    {shim} = fn_(real(wavelength, c_double))\n'
            '  end function {shim}\n'
        ),
        # Galacticus-side abstract interface, used to declare the module's
        # procedure-pointer slot when the method dummy is
        # `procedure(...), pointer` (gfortran cannot use a module procedure
        # from the same module's contains part as an interface-name in the
        # specification part, so the slot needs its own abstract interface
        # with characteristics matching the shim).
        'g_iface'   : (
            '  abstract interface\n'
            '     double precision function {giface}(wavelength)\n'
            '       double precision, intent(in   ) :: wavelength\n'
            '     end function {giface}\n'
            '  end interface\n'
        ),
        'shim_uses' : {},
    },
    # double precision function integrand(coordinates)
    #   class(coordinate), intent(in) :: coordinates
    # (source/computational_domains/volume_integrators/_class.F90).  The
    # coordinate is delivered to Python as its 3-element Cartesian array
    # via coordinate%toCartesian().
    'computationalDomainVolumeIntegrand': {
        'cfunctype' : 'CFUNCTYPE(c_double, POINTER(c_double))',
        'py_adapter': ('(lambda _pos_: float({fn}('
                       'np.ctypeslib.as_array(_pos_, (3,)))))'),
        'c_iface'   : (
            '  abstract interface\n'
            '     function {iface}(position) bind(c)\n'
            '       import c_double\n'
            '       real(c_double), dimension(3), intent(in   ) :: position\n'
            '       real(c_double)                              :: {iface}\n'
            '     end function {iface}\n'
            '  end interface\n'
        ),
        'shim'      : (
            '  double precision function {shim}(coordinates)\n'
            '    class(coordinate), intent(in   ) :: coordinates\n'
            '    procedure({iface}), pointer      :: fn_\n'
            '    real(c_double), dimension(3)     :: position_\n'
            '    position_ = coordinates%toCartesian()\n'
            '    call c_f_procpointer({var}, fn_)\n'
            '    {shim} = fn_(position_)\n'
            '  end function {shim}\n'
        ),
        'g_iface'   : (
            '  abstract interface\n'
            '     double precision function {giface}(coordinates)\n'
            '       import :: coordinate\n'
            '       class(coordinate), intent(in   ) :: coordinates\n'
            '     end function {giface}\n'
            '  end interface\n'
        ),
        'shim_uses' : {'Coordinates': ('coordinate',)},
    },
}


def _module_for_symbol(use_blocks, symbol):
    """Return the name of the module whose `use ..., only :` list includes
    `symbol`, or None.

    A parsed moduleUse node maps each module to a *list* of entries (one per
    preprocessor condition set); a bare dict is tolerated for safety.
    """
    for use_block in use_blocks or []:
        for mod_name, mod_data in use_block.items():
            entries = mod_data if isinstance(mod_data, list) else [mod_data]
            for entry in entries:
                if isinstance(entry, dict) and symbol in entry.get('only', {}):
                    return mod_name
    return None


# Match a fixed-size dimension attribute of any rank.  Each comma-
# separated dim is either `N` (default lower bound 1) or `L:U` (explicit
# lower bound).  Used to distinguish a fixed-size array — whose shape is
# known at codegen time and so needs no count companions — from a
# deferred-shape `dimension(:)` / `dimension(:,:)`, which do.  The
# wrapper always declares the bind(c) side with default lower bounds;
# the inner constructor's explicit-lower-bound dummy receives the data
# via Fortran's elementwise copy on entry, so callers don't have to
# care about base offsets.
_DIM_FIXED_RX = re.compile(
    r'^dimension\s*\(\s*'
    r'(?:\d+\s*:\s*)?\d+'                      # first axis
    r'(?:\s*,\s*(?:\d+\s*:\s*)?\d+)*'          # additional axes
    r'\s*\)$'
)
# Pulls one axis out of a fixed-shape attr's inner csv.
_DIM_AXIS_RX = re.compile(
    r'^\s*(?:(\d+)\s*:\s*)?(\d+)\s*$'
)


def _fixed_dim_shape(attr):
    """Return the per-axis element counts encoded by a fixed-shape
    dimension attribute as a tuple, or ``None`` if *attr* doesn't match
    a supported fixed-shape form.  Handles `dimension(N)`, `dimension(L:U)`,
    and multi-rank forms like `dimension(N, M)` / `dimension(N, M, K)`
    (each axis independently default-lower-bound or explicit `L:U`)."""
    if not _DIM_FIXED_RX.match(attr):
        return None
    inner = attr[attr.index('(') + 1 : attr.rindex(')')]
    dims = []
    for part in inner.split(','):
        m = _DIM_AXIS_RX.match(part)
        if not m:
            return None
        lower = int(m.group(1)) if m.group(1) is not None else 1
        upper = int(m.group(2))
        dims.append(upper - lower + 1)
    return tuple(dims)


def _fixed_dim_size(attr):
    """Return the total element count encoded by a fixed-shape 1D
    dimension attribute, or None if *attr* doesn't match.  Provided as
    a thin wrapper over :func:`_fixed_dim_shape` so call-sites that
    only care about total length stay terse."""
    shape = _fixed_dim_shape(attr)
    return None if shape is None else _shape_product(shape)


def _shape_product(shape):
    """Element count for a fixed-shape tuple."""
    n = 1
    for d in shape:
        n *= d
    return n

# Match a `len=N` literal in a character type-spec, used to detect
# fixed-length character arrays (e.g. `character(len=2), dimension(:)`).
# Variable-length forms (`len=*`, `len=:`) deliberately don't match — we
# can't pack those into a fixed-stride byte buffer at the Python boundary.
_CHAR_LEN_RX = re.compile(r'^len\s*=\s*(\d+)$')


def _make_id_array_companion(name):
    """Build a hidden 1D ``integer(c_int), dimension(*)`` companion ArgSpec
    for a ``polymorphic_list_array`` argument — the parallel buffer carrying
    each element's concrete class ID.  Mirrors the scalar functionClass
    arg path's ``_ID`` companion (set up inside ``assign_c_types``), only
    1D and never optional.

    The companion participates in ``assign_c_attributes`` like any other
    deferred-shape numeric input — ``dimension(:)`` is rewritten to
    ``dimension(*)`` and ``ctype_pointer`` becomes True so ctypes
    receives ``POINTER(c_int)`` and the Python side hands over a
    ``(c_int * count)(*ids)`` ctypes array.
    """
    return ArgSpec(
        name                  = name,
        intrinsic             = 'integer',
        type_spec             = '',
        attributes            = ['intent(in)', 'dimension(:)'],
        ctype                 = 'c_int',
        fort_type             = 'integer(c_int)',
        fort_is_present       = True,
        py_is_present         = False,
        galacticus_is_present = False,
        is_optional           = False,
        is_function_class     = False,
    )


def _make_count_companion(name):
    """Build a hidden c_size_t count companion ArgSpec for a deferred-shape
    array argument.  The companion appears in the bind(c) signature and
    is filled in by the Python wrapper from the input array's `.size`
    or `.shape[i]`; it's never visible in the user-facing Python signature
    nor passed to the inner Galacticus call.
    """
    return ArgSpec(
        name                  = name,
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


def assign_c_types(argument_list, lib_function_classes, class_hierarchy=None,
                   constructor_overrides=()):
    """Assign appropriate C types for each argument.

    Mirrors Perl assignCTypes().  Accepts a list of raw Fortran-declaration
    dicts and returns a new list of :class:`ArgSpec` objects.  Processes in
    reverse so that the ``_ID`` companion argument for functionClass parameters
    can be inserted immediately after its parent without disturbing the rest
    of the list.

    *class_hierarchy* is the type→parent map produced by
    :func:`LibraryInterfaces.Hierarchy.build_type_hierarchy`.  When supplied,
    ``class(<intermediate>)`` args whose extends-chain reaches a registered
    ``<base>Class`` get the function-class treatment too, with the
    intermediate type recorded as ``arg.narrowing_type`` for the Fortran
    emitter to ``select type``-narrow at runtime.

    *constructor_overrides* is the impl's ``<constructor><argument …/></constructor>``
    entries from libraryClasses.xml.  An entry with ``value="null"`` flags
    the matching arg as null-filled: it's dropped from both the bind(c)
    and the Python signatures (``fort_is_present`` / ``py_is_present`` set
    False) and the Fortran wrapper declares a local null pointer that the
    inner constructor receives in its place.  See ``is_null_filled`` on
    :class:`ArgSpec` for the rationale.
    """
    if class_hierarchy is None:
        class_hierarchy = {}
    null_filled_names = {
        o.get('name') for o in constructor_overrides
        if isinstance(o, dict) and o.get('value') == 'null' and o.get('name')
    }
    absent_filled_names = {
        o.get('name') for o in constructor_overrides
        if isinstance(o, dict) and o.get('value') == 'absent' and o.get('name')
    }

    # Does this argument list contain a 1D `intent(out), allocatable,
    # dimension(:)` numeric output array?  If so, scalar `intent(out)`
    # numeric/logical companions on the same method are also treated as
    # outputs (dropped from the Python input signature, filled value
    # returned) — see is_output_scalar handling near the end of the loop.
    # The whole-method gate (unsupported_output_array_method) guarantees
    # such methods are void-returning and optional-free, so this only ever
    # fires where the bespoke output-array Python emission handles it.
    def _attrs_intr(r):
        if isinstance(r, dict):
            return r.get('attributes', []), r.get('intrinsic')
        return getattr(r, 'attributes', []), getattr(r, 'intrinsic', '')
    has_output_array = False
    for r in argument_list:
        a_attrs, a_intr = _attrs_intr(r)
        if ('allocatable' in a_attrs
                and ('intent(out)' in a_attrs or 'intent(inout)' in a_attrs)
                and ((a_intr in ('double precision', 'integer')
                      and ('dimension(:)' in a_attrs
                           or 'dimension(:,:)' in a_attrs))
                     or (a_intr == 'logical'
                         and 'dimension(:)' in a_attrs))):
            has_output_array = True
            break

    new_list = []
    for raw in reversed(argument_list):
        arg = ArgSpec.from_raw(raw) if isinstance(raw, dict) else raw

        # Initialise presence flags (mirrors the original explicit reset).
        arg.fort_is_present       = True
        arg.py_is_present         = True
        arg.galacticus_is_present = True

        # Null-fill override: drop the arg from both wrappers and have
        # build_fortran_reassignments emit a local null pointer that the
        # inner constructor receives instead.  Short-circuits the rest
        # of the per-intrinsic dispatch — the ArgSpec carries no ctype
        # because the bind(c) signature won't reference it.
        if arg.name in null_filled_names:
            arg.is_null_filled        = True
            arg.fort_is_present       = False
            arg.py_is_present         = False
            arg.galacticus_is_present = True
            arg.is_optional           = False
            new_list.insert(0, arg)
            continue

        # Absent-fill override: drop the arg from both wrappers AND
        # from the inner call.  Only valid for optional args — the
        # inner constructor must handle the absence via its declared
        # default behaviour (e.g. identity-aligned principal axes for
        # the Gaussian ellipsoid's optional `axes`).  Detecting
        # `is_optional` requires the `optional` attribute on the raw
        # decl, which we've already captured in arg.attributes.
        if arg.name in absent_filled_names:
            if 'optional' not in arg.attributes:
                raise ValueError(
                    f"value='absent' override on non-optional argument "
                    f"'{arg.name}' — only optional args may be dropped")
            arg.is_absent_filled      = True
            arg.fort_is_present       = False
            arg.py_is_present         = False
            arg.galacticus_is_present = False
            arg.is_optional           = False
            new_list.insert(0, arg)
            continue

        # Output-array arg: 1D `allocatable, dimension(:)` numeric, either
        # `intent(out)` or `intent(inout)` (the latter is Galacticus's
        # allocation-reuse idiom — see Classification.is_output_array_arg).
        # The inner method allocates and fills it, and its data + size flow
        # *back* to Python via a `(c_ptr, c_size_t)` companion pair that the
        # generator appends to the bind(c) signature (see interfaces_methods).
        # Short-circuit the per-intrinsic dispatch below, and in particular
        # the deferred-shape *input*-array block that would otherwise treat
        # this like an incoming numpy buffer and add a count companion: the
        # arg is a wrapper-local `save, target` allocatable, not a bind(c) or
        # Python formal, so it carries no ctype. Only galacticus_is_present
        # stays True, so the inner call still receives the local buffer.
        # (Predicate inlined rather than imported from Classification, which
        # already imports from this module.)
        if ('allocatable' in arg.attributes
                and ('intent(out)' in arg.attributes
                     or 'intent(inout)' in arg.attributes)
                and ((arg.intrinsic in ('double precision', 'integer')
                      and ('dimension(:)' in arg.attributes
                           or 'dimension(:,:)' in arg.attributes))
                     or (arg.intrinsic == 'logical'
                         and 'dimension(:)' in arg.attributes))):
            if arg.intrinsic == 'double precision':
                elem_ctype, elem_fort, elem_dtype = (
                    'c_double', 'real(c_double)', 'float64')
            elif arg.intrinsic == 'logical':
                elem_ctype, elem_fort, elem_dtype = (
                    'c_bool', 'logical(c_bool)', 'bool_')
            elif arg.type_spec in ('c_long', 'kind_int8'):
                elem_ctype, elem_fort, elem_dtype = (
                    'c_long', 'integer(c_long)', 'int64')
            elif arg.type_spec == 'c_size_t':
                elem_ctype, elem_fort, elem_dtype = (
                    'c_size_t', 'integer(c_size_t)', 'uint64')
            else:
                elem_ctype, elem_fort, elem_dtype = (
                    'c_int', 'integer(c_int)', 'int32')
            rank  = 2 if 'dimension(:,:)' in arg.attributes else 1
            dims  = ':,:' if rank == 2 else ':'
            local = f'glcOut_{arg.name}_'
            arg.is_output_array       = True
            arg.array_rank            = rank
            arg.output_elem_ctype     = elem_ctype
            arg.output_elem_fort      = elem_fort
            arg.output_elem_dtype     = elem_dtype
            arg.fort_is_present       = False
            arg.py_is_present         = False
            arg.galacticus_is_present = True
            arg.is_optional           = False
            if arg.intrinsic == 'logical':
                # The inner dummy is default-kind `logical`, whose storage
                # width differs from c_bool (4 bytes vs 1 on gfortran), so
                # the export buffer can't be filled directly: the inner
                # method fills a default-kind local, and the generator
                # emits a kind-narrowing copy into the c_bool save buffer
                # after the call (the reverse of the logical *input*
                # path's widening copy).
                inner = f'glcOutInner_{arg.name}_'
                arg.fort_pass_as      = inner
                arg.fort_declarations = (
                    f'  {elem_fort}, dimension({dims}), allocatable, save,'
                    f' target :: {local}\n'
                    f'  logical, dimension({dims}), allocatable'
                    f' :: {inner}\n')
                arg.fort_reassignment = (
                    f'  if (allocated({local})) deallocate({local})\n'
                    f'  if (allocated({inner})) deallocate({inner})\n')
            else:
                arg.fort_pass_as      = local
                arg.fort_declarations = (
                    f'  {elem_fort}, dimension({dims}), allocatable, save,'
                    f' target :: {local}\n')
                # Pre-call hygiene: free any allocation the previous call
                # left (Python has already copied it) so the inner method's
                # allocate-on-assignment starts from a clean deallocated
                # state.
                arg.fort_reassignment = (
                    f'  if (allocated({local})) deallocate({local})\n')
            new_list.insert(0, arg)
            continue

        # Inbound callback arg: `procedure(<iface>) :: fn` with <iface>
        # registered in _CALLBACK_PROCEDURE_INTERFACES.  The wrapper receives
        # a C funptr by value; everything else (module emission, shim
        # pass-through, funptr storage, CFUNCTYPE adaptation) is wired by the
        # generator, which knows the class/method names the module and shim
        # are named after.  `pointer` dummies are accepted too — the
        # generator passes the shim module's procedure-pointer slot instead
        # of the shim itself (keyed off the arg's `pointer` attribute) —
        # EXCEPT with intent(out|inout), where the pointer is an output
        # (classify_arg rejects the surrounding method, so the generator
        # never reaches here for those; the attribute check keeps this
        # predicate self-contained anyway).
        if (arg.intrinsic == 'procedure'
                and arg.type_spec in _CALLBACK_PROCEDURE_INTERFACES
                and 'optional' not in arg.attributes
                and not ('pointer' in arg.attributes
                         and ('intent(out)' in arg.attributes
                              or 'intent(inout)' in arg.attributes))):
            arg.is_callback = True
            arg.ctype       = 'c_void_p'
            arg.fort_type   = 'type(c_funptr)'
            arg.is_optional = False
            new_list.insert(0, arg)
            continue

        arg.is_function_class     = False
        arg.is_optional           = bool('optional' in arg.attributes)

        intrinsic     = arg.intrinsic
        type_spec_val = arg.type_spec

        # Drop optional `integer(omp_lock_kind)` args entirely.  The OpenMP
        # lock kind is platform-dependent (INTEGER(4) on Linux GCC,
        # INTEGER(8) on macOS GCC); mapping it to a fixed C-interop kind
        # picks the wrong size on one platform and produces a "passed
        # INTEGER(4) to INTEGER(8)" type-mismatch when the wrapper calls
        # the inner Galacticus method.  These args are also semantically
        # meaningless to a Python caller — they coordinate threads inside
        # the Fortran code — so simply omitting them and letting the
        # inner method's `optional` default kick in is correct.
        # `_unsupported_arg` rejects the surrounding method outright when
        # the arg is non-optional, since dropping a required arg would
        # silently mis-call the inner method.
        if intrinsic == 'integer' and type_spec_val == 'omp_lock_kind':
            continue

        if intrinsic == 'double precision':
            arg.ctype     = 'c_double'
            arg.fort_type = 'real(c_double)'
        elif intrinsic == 'real':
            # Default-kind `real` is REAL(4) on every platform Galacticus
            # targets; the matching C-interop kind is `c_float`.  Without
            # this branch a scalar `real` arg fell through to the
            # ArgSpec default of `c_int`, and the emitted wrapper
            # declared a 4-byte INTEGER where the inner method expected
            # a 4-byte REAL — broken at compile time (see
            # `posteriorSampleLikelihood::evaluate`'s `timeEvaluate`).
            arg.ctype     = 'c_float'
            arg.fort_type = 'real(c_float)'
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
            char_len_m = _CHAR_LEN_RX.match((type_spec_val or '').strip())
            if char_len_m and 'dimension(:)' in arg.attributes:
                # Fixed-length character array (`character(len=N),
                # dimension(:)`).  bind(c) receives a flat
                # `character(c_char), dimension(*)` byte buffer of total
                # length count*N, plus the count companion inserted
                # below; the wrapper repacks into `character(len=N),
                # dimension(:)` before calling the inner method.  The
                # `_p` suffix would otherwise drag this through the
                # null-terminated-string scalar path in
                # `build_fortran_reassignments`, which is wrong here.
                arg.ctype     = 'c_char'
                arg.fort_type = 'character(c_char)'
                arg.char_len  = int(char_len_m.group(1))
                # Fall through into the array-detection block below; it
                # already handles deferred-shape `dimension(:)` for the
                # numeric case, and we set `is_array` etc. there.
            else:
                # Scalar `character(c_char)` / `character(len=*)` / etc.
                # — null-terminated C string round-trip.
                arg.ctype     = 'c_char_p'
                arg.fort_type = 'character(c_char)'
        elif intrinsic == 'type':
            if type_spec_val == 'varying_string':
                if 'dimension(:)' in arg.attributes:
                    # `type(varying_string), dimension(:)`.  No fixed
                    # per-element width is known at codegen time, so we
                    # ship a flat `character(c_char), dimension(*)` byte
                    # buffer plus a count companion AND a per-element
                    # length companion (both c_size_t, by value).  The
                    # Python wrapper picks the max ASCII-encoded length
                    # at runtime, pads each element to that width, and
                    # passes the contiguous buffer; the Fortran wrapper
                    # repacks into `type(varying_string), dimension(:)`
                    # via element-wise assignment from a per-element
                    # `character(len=:), allocatable` scratch buffer.
                    arg.ctype                = 'c_char'
                    arg.fort_type            = 'character(c_char)'
                    arg.varying_string_array = True
                    # Fall through to the array-detection block below;
                    # it handles the deferred-shape branch and inserts
                    # the (single) count companion.  The second
                    # per-element-length companion is inserted there too
                    # so its position is fixed relative to the count.
                else:
                    arg.ctype     = 'c_char_p'
                    arg.fort_type = 'character(c_char)'
            elif re.match(r'^enumeration[a-z0-9_]+type$', type_spec_val, re.IGNORECASE):
                # Enumeration types map to C int.
                arg.ctype     = 'c_int'
                arg.fort_type = 'integer(c_int)'
            elif (type_spec_val.endswith('List')
                  and type_spec_val[:-4] in lib_function_classes
                  and 'dimension(:)' in arg.attributes
                  and type_spec_val in _SHARED_TYPE_MODULES):
                # `type(<class>List), dimension(:)` — Galacticus's idiom
                # for "array of class(<class>Class)".  We only handle
                # this shape when the List wrapper itself is registered
                # in `_SHARED_TYPE_MODULES`: that's our promise that
                # (a) the type is exported from a known module so the
                # bind(c) host can `use ::` it, and (b) the wrapper has
                # the canonical thin-pointer layout (one component
                # named `<class>_`).  Locally-defined `<class>List`
                # structs that happen to follow the naming convention
                # but carry extra members — e.g.
                # `virialDensityContrastList` in
                # tasks.halo_mass_function.F90, which adds a `label`
                # field — would otherwise silently drop those extra
                # members when the wrapper rebuilds the list.  Falling
                # through to the generic derived-type branch lets
                # `_unsupported_arg` reject the surrounding constructor,
                # which is the correct outcome.
                stem                          = type_spec_val[:-4]
                arg.ctype                     = 'c_void_p'
                arg.fort_type                 = 'type(c_ptr)'
                arg.polymorphic_list_array    = True
                arg.polymorphic_list_class    = stem
                arg.polymorphic_list_type     = type_spec_val
                arg.polymorphic_list_component = stem + '_'
                # Insert the IDs array companion (sits between the
                # parent and the count once the count is added below).
                # Order at the end of this iteration: parent, count, IDs.
                new_list.insert(0, _make_id_array_companion(arg.name + '_IDs'))
                # Fall through to the array-detection block below; it
                # marks is_array and inserts the count companion.
            else:
                arg.ctype     = 'c_void_p'
                arg.fort_type = 'type(c_ptr)'
        elif intrinsic == 'class':
            arg.ctype     = 'c_void_p'
            arg.fort_type = 'type(c_ptr)'
            # Check whether this is a functionClass argument.  Two shapes
            # qualify: the plain `class(<base>Class)` form (stem registered
            # in libraryClasses.xml), or `class(<intermediate>)` where
            # <intermediate>'s extends-chain reaches a registered
            # <base>Class — in which case we route through <base>GetPtr and
            # have build_fortran_reassignments insert a `select type`
            # narrowing to the intermediate.
            class_key      = None
            narrowing_type = ''
            if type_spec_val.endswith('Class') \
                    and type_spec_val[:-5] in lib_function_classes:
                class_key = type_spec_val[:-5]
            elif class_hierarchy:
                base, intermediate = resolve_function_class_base(
                    type_spec_val, class_hierarchy,
                    set(lib_function_classes.keys()))
                if base is not None:
                    class_key      = base
                    narrowing_type = intermediate or ''
            if class_key is not None:
                arg.is_function_class   = True
                arg.narrowing_type      = narrowing_type
                # Record the registered <base> here so
                # build_fortran_reassignments doesn't have to re-derive it
                # (the abstract-intermediate path can't recover it from
                # type_spec alone — only the hierarchy walk knows).
                arg.fort_function_class = class_key
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
                    arg_id.py_present  = python_safe_name(arg.name)
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
        # Numeric arrays go straight through the array path; character
        # arrays piggyback on the deferred-shape branch when their
        # per-element length was resolved above (a fixed-len character
        # array sets arg.char_len > 0; variable-len fell back to the
        # scalar path and is_array stays False).
        is_numeric_array_candidate = intrinsic in ('double precision', 'integer')
        is_logical_array_candidate = intrinsic == 'logical'
        is_char_array_candidate    = intrinsic == 'character' and arg.char_len > 0
        is_vstr_array_candidate    = arg.varying_string_array
        is_list_array_candidate    = arg.polymorphic_list_array
        if (is_numeric_array_candidate or is_logical_array_candidate
                or is_char_array_candidate
                or is_vstr_array_candidate or is_list_array_candidate):
            if 'dimension(:)' in arg.attributes:
                arg.is_array   = True
                arg.array_size = None
                count_arg = _make_count_companion(arg.name + '_count')
                new_list.insert(0, count_arg)
                if is_vstr_array_candidate:
                    # Per-element ASCII byte length (decided at runtime by
                    # the Python wrapper from the input list's max length);
                    # placed AFTER the count companion in bind(c) order, so
                    # we insert it first while iterating in reverse.
                    char_len_arg = _make_count_companion(arg.name + '_charLen')
                    new_list.insert(1, char_len_arg)
            elif is_numeric_array_candidate and 'dimension(:,:)' in arg.attributes:
                # 2D deferred-shape numeric array: two count companions,
                # one per axis.  The wrapper passes a flat C buffer plus
                # the two dimensions; inside, we allocate a
                # `dimension(:,:)` local of the right shape and reshape
                # the flat input into it (column-major Fortran order, so
                # the Python boundary uses np.asfortranarray).
                arg.is_array   = True
                arg.array_size = None
                arg.array_rank = 2
                count_arg_1 = _make_count_companion(arg.name + '_count_1')
                count_arg_2 = _make_count_companion(arg.name + '_count_2')
                # Insert in reverse so final order is arg, count_1, count_2.
                new_list.insert(0, count_arg_2)
                new_list.insert(0, count_arg_1)
            elif is_numeric_array_candidate:
                # Fixed-shape numeric arrays (`dimension(N)`,
                # `dimension(L:U)`, or multi-rank `dimension(N1,N2,...)`)
                # need no count companion; every axis size is in the
                # dimension spec.  For 1D we record `array_size` only;
                # for rank > 1 we also record `array_shape` so the
                # Python wrapper can validate ndim/shape rather than
                # only the total element count.  Fixed-size character
                # arrays aren't supported yet — the only registered
                # class needing them was already covered by the
                # deferred-shape path.
                for a in arg.attributes:
                    shape = _fixed_dim_shape(a)
                    if shape is not None:
                        arg.is_array    = True
                        arg.array_size  = _shape_product(shape)
                        if len(shape) > 1:
                            arg.array_shape = shape
                            arg.array_rank  = len(shape)
                        break

        # Scalar `intent(out)` numeric/logical companion on an output-array
        # method: keep it in the bind(c) signature as an ordinary
        # by-reference intent(out) scalar (assign_c_attributes already picks
        # reference pass-by for it) passed straight to the inner method, but
        # drop it from the Python input signature and mark it so the bespoke
        # output-array Python emission returns its filled value.  Guarded by
        # has_output_array so non-output-array methods keep their existing
        # scalar-arg handling unchanged.  (Predicate inlined — Classification
        # imports from this module — but kept in lock-step with
        # Classification.is_output_scalar_arg.)
        if (has_output_array
                and not arg.is_output_array
                and not arg.is_array
                and arg.intrinsic in ('double precision', 'real',
                                      'integer', 'logical')
                and 'intent(out)' in arg.attributes
                and 'allocatable' not in arg.attributes
                and not any(a.startswith('dimension') for a in arg.attributes)):
            arg.is_output_scalar = True
            arg.py_is_present     = False
            # fort_is_present / galacticus_is_present stay True: the arg is a
            # real by-reference intent(out) dummy the inner method fills.

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
        # `dimension(:,:)` (rank-2 assumed-shape) is also illegal in
        # bind(c); the wrapper takes a rank-1 flat buffer plus two
        # count companions and rebuilds a rank-2 view inside via
        # `reshape`, so the bind(c) signature collapses to
        # `dimension(*)` for both rank-1 and rank-2 deferred shapes.
        attr_filters = []
        for a in arg.attributes:
            if a in ('dimension(:)', 'dimension(:,:)'):
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
            # Galacticus arg names that collide with Python keywords
            # (`lambda`, `class`, ...) need to be escaped before being
            # spliced into the generated Python source — every emitted
            # `if {name}:` / `{name}._glcObj` would otherwise be a
            # SyntaxError in the wrapper module.  Fortran-side identifier
            # is unaffected (Fortran allows `lambda` as a name).
            name = python_safe_name(arg.name)
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
        elif arg.is_array and arg.array_rank == 2 and not arg.array_shape:
            # 2D deferred-shape numeric array.  Convert input to a
            # contiguous Fortran-order numpy array (column-major) so the
            # raw byte layout matches Fortran's `dimension(N, M)`
            # convention; pass the data pointer plus shape[0]/shape[1]
            # as the two count companions, in the same order they sit at
            # the front of new_list (count_1 then count_2 — see
            # assign_c_types).  `array_shape` non-empty signals a
            # fixed-shape array — there are no count companions to pop
            # (assign_c_types didn't insert them), and the fixed-shape
            # branch further down handles those.  For optional args, the conversion is
            # gated on the input not being None — otherwise
            # `np.asarray(None, dtype=...)` produces a 0D scalar that
            # `np.asfortranarray` then promotes to 1D, tripping the
            # ndim check below.  When absent we pass NULL as the data
            # pointer and 0 for both counts; the bind(c) wrapper guards
            # the reshape on `present()` so the values are unused.
            np_dtype = _ARRAY_NUMPY_DTYPE.get(arg.ctype, 'float64')
            safe = python_safe_name(arg.name)
            count_arg_1 = new_list.pop(0)
            count_arg_2 = new_list.pop(0)
            if arg.is_optional:
                arg.py_reassignment = (
                    f'    if {safe} is not None:\n'
                    f'        {safe} = np.asfortranarray('
                    f'np.asarray({safe}, dtype=np.{np_dtype}))\n'
                    f'        if {safe}.ndim != 2:\n'
                    f'            raise ValueError('
                    f'f"{arg.name} expects a 2D array, got '
                    f'{{{safe}.ndim}}D")\n'
                )
                # The count companions' `py_pass_as` deliberately doesn't
                # wrap the result in `c_size_t(...)` — `python_call_code`
                # adds that wrap once the optional region begins, and a
                # second wrap (`c_size_t(c_size_t(...))`) is rejected by
                # ctypes ("'c_ulong' object cannot be interpreted as an
                # integer").  Same logic for the 1D and character branches
                # below.
                count_arg_1.py_pass_as = f'{safe}.shape[0] if {safe} is not None else 0'
                count_arg_2.py_pass_as = f'{safe}.shape[1] if {safe} is not None else 0'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER({arg.ctype})) if {safe} is not None else None'
                )
            else:
                arg.py_reassignment = (
                    f'    {safe} = np.asfortranarray('
                    f'np.asarray({safe}, dtype=np.{np_dtype}))\n'
                    f'    if {safe}.ndim != 2:\n'
                    f'        raise ValueError('
                    f'f"{arg.name} expects a 2D array, got '
                    f'{{{safe}.ndim}}D")\n'
                )
                count_arg_1.py_pass_as = f'c_size_t({safe}.shape[0])'
                count_arg_2.py_pass_as = f'c_size_t({safe}.shape[1])'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER({arg.ctype}))'
                )
            new_list.insert(0, count_arg_2)
            new_list.insert(0, count_arg_1)
        elif arg.is_array and arg.varying_string_array:
            # `type(varying_string), dimension(:)`: encode each element to
            # ASCII, pad to the max encoded length, and ship as a flat
            # contiguous `S{N}` numpy buffer.  N is decided at *runtime*
            # (we can't know it at codegen), so the count companion's
            # `py_pass_as` reads `.size` and the per-element-length
            # companion reads `.itemsize` — numpy's `S{N}` dtype packs
            # each element exactly N bytes wide, which is what the
            # bind(c) buffer expects.  Empty / all-empty lists fall back
            # to N=1 so the dtype is well-formed.  Optional-arg handling
            # mirrors the fixed-len-character-array branch below.
            safe = python_safe_name(arg.name)
            count_arg    = new_list.pop(0)
            char_len_arg = new_list.pop(0)
            encode_block = (
                f"[(s if isinstance(s,(bytes,bytearray)) else str(s).encode('ascii'))"
                f" for s in {safe}]"
            )
            # Each element is `ljust`-padded with **spaces** to the
            # max-length width before being handed to numpy.  numpy's
            # `S{N}` dtype otherwise NUL-pads short bytes-strings,
            # which Fortran `trim()` (the unpacking code uses) doesn't
            # strip — producing varying_string elements with embedded
            # NUL bytes that fail equality against cleanly-encoded
            # `varying_string` scalars.  Pre-padding with spaces means
            # every buffer entry is exactly `_glcLen_` bytes wide AND
            # correctly trim-able on the Fortran side.
            if arg.is_optional:
                arg.py_reassignment = (
                    f"    if {safe} is not None:\n"
                    f"        {safe}_glcStrs_ = {encode_block}\n"
                    f"        {safe}_glcLen_ = max((len(b) for b in {safe}_glcStrs_), default=1) or 1\n"
                    f"        {safe}_glcStrs_ = [b.ljust({safe}_glcLen_, b' ') for b in {safe}_glcStrs_]\n"
                    f"        {safe} = np.ascontiguousarray("
                    f"np.array({safe}_glcStrs_, dtype=f'S{{{safe}_glcLen_}}'))\n"
                )
                count_arg.py_pass_as    = f'{safe}.size if {safe} is not None else 0'
                char_len_arg.py_pass_as = f'{safe}.itemsize if {safe} is not None else 0'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER(c_char)) if {safe} is not None else None'
                )
            else:
                arg.py_reassignment = (
                    f"    {safe}_glcStrs_ = {encode_block}\n"
                    f"    {safe}_glcLen_ = max((len(b) for b in {safe}_glcStrs_), default=1) or 1\n"
                    f"    {safe}_glcStrs_ = [b.ljust({safe}_glcLen_, b' ') for b in {safe}_glcStrs_]\n"
                    f"    {safe} = np.ascontiguousarray("
                    f"np.array({safe}_glcStrs_, dtype=f'S{{{safe}_glcLen_}}'))\n"
                )
                count_arg.py_pass_as    = f'c_size_t({safe}.size)'
                char_len_arg.py_pass_as = f'c_size_t({safe}.itemsize)'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER(c_char))'
                )
            new_list.insert(0, char_len_arg)
            new_list.insert(0, count_arg)
        elif arg.is_array and arg.polymorphic_list_array:
            # `type(<class>List), dimension(:)`: ship parallel ctypes
            # arrays of object pointers + class IDs, one per element,
            # built from the `_glcObj` / `_classID` fields the Python
            # wrapper attaches to every constructor-built functionClass
            # instance.  The bind(c) signature receives the pointer
            # buffer as `type(c_ptr), dimension(*)` (see
            # assign_c_attributes; ctypes side is `POINTER(c_void_p)`)
            # and the IDs companion as `integer(c_int), dimension(*)`.
            # Optional-arg handling matches the other array branches.
            safe         = python_safe_name(arg.name)
            count_arg    = new_list.pop(0)
            ids_arg      = new_list.pop(0)
            if arg.is_optional:
                arg.py_reassignment = (
                    f'    if {safe} is not None:\n'
                    f'        {safe}_glcN_    = len({safe})\n'
                    f'        {safe}_glcPtrs_ = (c_void_p * {safe}_glcN_)('
                    f'*[x._glcObj for x in {safe}])\n'
                    f'        {safe}_glcIDs_  = (c_int    * {safe}_glcN_)('
                    f'*[x._classID for x in {safe}])\n'
                )
                count_arg.py_pass_as = f'{safe}_glcN_ if {safe} is not None else 0'
                ids_arg.py_pass_as   = f'{safe}_glcIDs_ if {safe} is not None else None'
                arg.py_pass_as       = f'{safe}_glcPtrs_ if {safe} is not None else None'
            else:
                arg.py_reassignment = (
                    f'    {safe}_glcN_    = len({safe})\n'
                    f'    {safe}_glcPtrs_ = (c_void_p * {safe}_glcN_)('
                    f'*[x._glcObj for x in {safe}])\n'
                    f'    {safe}_glcIDs_  = (c_int    * {safe}_glcN_)('
                    f'*[x._classID for x in {safe}])\n'
                )
                count_arg.py_pass_as = f'c_size_t({safe}_glcN_)'
                ids_arg.py_pass_as   = f'{safe}_glcIDs_'
                arg.py_pass_as       = f'{safe}_glcPtrs_'
            new_list.insert(0, ids_arg)
            new_list.insert(0, count_arg)
        elif arg.is_array and arg.intrinsic == 'character':
            # Fixed-length character array: build a contiguous count*N
            # byte buffer.  Each input element is `str()`-coerced,
            # right-padded (Fortran convention) and clipped to N chars,
            # then ASCII-encoded.  numpy's `S{N}` dtype gives us a
            # fixed-stride buffer whose `.size` is the element count and
            # whose underlying memory is exactly the N-byte-per-element
            # layout the inner Fortran method expects.  Optional-arg
            # handling matches the 2D-array branch above.
            safe = python_safe_name(arg.name)
            n    = arg.char_len
            count_arg = new_list.pop(0)
            if arg.is_optional:
                arg.py_reassignment = (
                    f"    if {safe} is not None:\n"
                    f"        {safe} = np.ascontiguousarray("
                    f"np.array("
                    f"[(str(s).ljust({n}))[:{n}].encode('ascii') for s in {safe}],"
                    f" dtype='S{n}'))\n"
                )
                # See the 2D branch above — let python_call_code's outer
                # `c_size_t(...)` wrap supply the conversion; double-wrap
                # is rejected by ctypes.
                count_arg.py_pass_as = f'{safe}.size if {safe} is not None else 0'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER(c_char)) if {safe} is not None else None'
                )
            else:
                arg.py_reassignment = (
                    f"    {safe} = np.ascontiguousarray("
                    f"np.array("
                    f"[(str(s).ljust({n}))[:{n}].encode('ascii') for s in {safe}],"
                    f" dtype='S{n}'))\n"
                )
                count_arg.py_pass_as = f'c_size_t({safe}.size)'
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER(c_char))'
                )
            new_list.insert(0, count_arg)
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
            # See the function-class branch above for the rationale —
            # `lambda` / `class` / ... can appear as a Galacticus arg name
            # and must be escaped before splicing into Python source.
            safe = python_safe_name(arg.name)
            if arg.array_size is None:
                count_arg = new_list.pop(0)
                if arg.is_optional:
                    # Optional 1D array: gate on None so that
                    # `np.ascontiguousarray(None, dtype=...)` (which
                    # produces a 0-D scalar in modern numpy and breaks
                    # downstream `.size` / `.ctypes` access) never runs.
                    # The count companion's `py_pass_as` doesn't wrap
                    # in `c_size_t(...)` — `python_call_code` adds that
                    # wrap once the optional region begins; double-wrap
                    # is rejected by ctypes.
                    arg.py_reassignment = (
                        f'    if {safe} is not None:\n'
                        f'        {safe} = np.ascontiguousarray({safe},'
                        f' dtype=np.{np_dtype})\n'
                    )
                    count_arg.py_pass_as = f'{safe}.size if {safe} is not None else 0'
                    arg.py_pass_as = (
                        f'{safe}.ctypes.data_as(POINTER({arg.ctype})) if {safe} is not None else None'
                    )
                else:
                    arg.py_reassignment = (
                        f'    {safe} = np.ascontiguousarray({safe},'
                        f' dtype=np.{np_dtype})\n'
                    )
                    count_arg.py_pass_as = f'c_size_t({safe}.size)'
                    arg.py_pass_as = (
                        f'{safe}.ctypes.data_as(POINTER({arg.ctype}))'
                    )
                new_list.insert(0, count_arg)
            elif arg.array_shape:
                # Fixed-shape multi-rank numeric array (e.g.
                # `dimension(3,2)`, `dimension(2,2,3)`).  Convert to a
                # Fortran-ordered (column-major) contiguous buffer so
                # the raw byte layout matches Fortran's shape spec;
                # validate `.shape` against the expected tuple — a
                # `.size`-only check would accept arbitrary
                # reshape-equivalent inputs.  No count companions are
                # needed: the bind(c) declaration carries the full
                # shape, and the inner constructor's dummy is
                # explicit-shape.
                expected_shape = tuple(arg.array_shape)
                arg.py_reassignment = (
                    f'    {safe} = np.asfortranarray('
                    f'np.asarray({safe}, dtype=np.{np_dtype}))\n'
                    f'    if {safe}.shape != {expected_shape}:\n'
                    f'        raise ValueError('
                    f'f"{arg.name} expects shape {expected_shape}, got '
                    f'{{{safe}.shape}}")\n'
                )
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER({arg.ctype}))'
                )
            else:
                # Fixed-size 1D: `array_size` is the literal length.
                arg.py_reassignment = (
                    f'    {safe} = np.ascontiguousarray({safe},'
                    f' dtype=np.{np_dtype})\n'
                    f'    if {safe}.size != {arg.array_size}:\n'
                    f'        raise ValueError('
                    f'f"{arg.name} expects {arg.array_size} elements, got '
                    f'{{{safe}.size}}")\n'
                )
                arg.py_pass_as = (
                    f'{safe}.ctypes.data_as(POINTER({arg.ctype}))'
                )
        elif arg.ctype == 'c_char_p':
            # Scalar `character(...)` and `type(varying_string)` args both
            # cross the boundary as a `c_char_p` (NUL-terminated C string),
            # which ctypes implements with Python `bytes`.  Python users
            # naturally pass `str` though, and ctypes won't auto-encode —
            # it raises "bytes or integer address expected instead of str
            # instance".  Encode at the boundary; `None` and `bytes` pass
            # through unchanged.
            safe = python_safe_name(arg.name)
            arg.py_reassignment = (
                f'    if isinstance({safe}, str):\n'
                f'        {safe} = {safe}.encode("utf-8")\n'
            )
        new_list.insert(0, arg)               # unshift current arg
    return new_list


# Numpy dtype string used when converting a Python input to a contiguous
# array for an array-typed argument.  Keep in sync with the ctypes mapping
# in assign_c_types.
_ARRAY_NUMPY_DTYPE = {
    'c_double': 'float64',
    'c_float' : 'float32',
    'c_int'   : 'int32',
    'c_long'  : 'int64',
    'c_size_t': 'uint64',
    # `c_bool` is one byte per element on every platform numpy supports;
    # numpy's `bool_` dtype is also one byte, so the in-memory layout
    # matches what the Fortran-side `logical(c_bool), dimension(*)`
    # buffer expects.
    'c_bool'  : 'bool_',
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

        if arg.is_absent_filled:
            # Absent-fill override (`<argument name="..." value="absent"/>`).
            # The arg is invisible from Python, the bind(c) signature,
            # AND the inner constructor call — galacticus_is_present is
            # False so fortran_call_code's `make_call` excludes it.
            # No Fortran declaration, no module import, no reassignment;
            # otherwise the type-handling branches below would emit a
            # spurious `use :: <fc_module>, only : <type>` for a type
            # the wrapper neither declares nor passes (e.g. `vector` /
            # `matrix` for the Gaussian-ellipsoid `axes` / `rotation`).
            new_list.insert(0, arg)
            continue

        if arg.is_output_array or arg.is_callback:
            # Output arrays and callbacks were fully wired in
            # assign_c_types / the generator: fort_pass_as points at the
            # wrapper-local buffer (or shim), fort_declarations declares
            # it, and fort_reassignment carries the pre-call dealloc (or
            # funptr store).  The per-intrinsic branches below must not
            # re-dispatch on the arg's declared type — the *logical*
            # branch in particular would overwrite the wiring with the
            # input path's c_bool→logical widening copy, leaving the
            # wrapper referencing undeclared buffer names.
            new_list.insert(0, arg)
            continue

        if arg.is_null_filled:
            # Null-fill override (`<argument name="..." value="null"/>`).
            # The arg is invisible from both Python and the bind(c)
            # signature; the inner constructor receives a local
            # disassociated pointer whose declared type matches the
            # original Fortran arg.  Supported intrinsics today:
            #   procedure(<intf>), pointer :: name => null()
            #   class(*), pointer :: name => null()
            # `<intf>` (e.g. sphericalAdiabaticGnedin2004Initializor)
            # is imported from the functionClass's own module, which is
            # also where the impl's `abstract interface` block lives.
            if intrinsic == 'procedure':
                arg.fort_declarations = (
                    f'procedure({type_spec_val}), pointer :: {name} => null()\n'
                )
                fc_module = func_class.get('module', '')
                if fc_module and type_spec_val:
                    arg.fort_modules.setdefault(fc_module, {})[type_spec_val] = 1
            elif intrinsic == 'class' and type_spec_val == '*':
                arg.fort_declarations = (
                    f'class(*), pointer :: {name} => null()\n'
                )
            else:
                # Other intrinsics aren't supported under value="null"
                # today; if we ever extend coverage, the per-intrinsic
                # declaration goes here.  Falling through with no
                # declaration would silently break the inner call.
                raise ValueError(
                    f"value='null' override on argument '{name}' of "
                    f"intrinsic '{intrinsic}'/type '{type_spec_val}' — "
                    "only procedure and class(*) args are supported")
            arg.fort_pass_as = name
            new_list.insert(0, arg)
            continue

        if arg.is_array and arg.array_rank == 2 and not arg.array_shape:
            # 2D deferred-shape numeric array.  bind(c) receives the
            # buffer as a flat `dimension(*)` block plus two c_size_t
            # counts.  Allocate a `dimension(:,:)` local of the right
            # shape and reshape the flat slice into it (column-major
            # Fortran order — the Python side already laid the buffer
            # out F-contiguous via np.asfortranarray).  Allocatable
            # sidesteps any concerns about forward-referencing the
            # dummy counts in an automatic-array size spec, and the
            # local is auto-deallocated when the wrapper returns.
            # `array_shape` non-empty signals a fixed-shape array
            # whose bind(c) declaration already carries the full
            # rank-N shape (see the fall-through `elif arg.is_array`
            # branch below) — no reshape is needed.
            #
            # Optional args are gated on `present(...)` because the
            # bind(c) buffer is a NULL pointer when the user passed
            # None on the Python side; slicing it would segfault.
            arg.fort_declarations = (
                f'  {arg.fort_type}, dimension(:,:), allocatable'
                f' :: {name}_F_\n'
            )
            reassign = (
                f'allocate({name}_F_({name}_count_1, {name}_count_2))\n'
                f'{name}_F_ = reshape('
                f'{name}(1:{name}_count_1*{name}_count_2),'
                f' [{name}_count_1, {name}_count_2])\n'
            )
            if is_optional:
                reassign = (
                    f'if (present({name})) then\n'
                    f'   {reassign.replace(chr(10), chr(10) + "   ").rstrip()}\n'
                    f'end if\n'
                )
            arg.fort_reassignment = reassign
            arg.fort_pass_as = f'{name}_F_'
        elif arg.is_array and arg.varying_string_array:
            # `type(varying_string), dimension(:)`.  bind(c) declares the
            # buffer as a flat `character(c_char), dimension(*)` block of
            # `name_count * name_charLen` bytes; the inner Galacticus
            # method wants `type(varying_string), dimension(:)`.  Walk
            # the buffer element-by-element, copy each N-byte slice into
            # a `character(len=:)` scratch local, then assign to the
            # varying_string element (with `trim()` to drop the
            # right-padding spaces the Python wrapper inserted to match
            # numpy's fixed-stride `S{N}` dtype).
            arg.fort_declarations = (
                f'  type(varying_string), dimension(:), allocatable'
                f' :: {name}_F_\n'
                f'  character(len=:), allocatable :: {name}_buf_\n'
                f'  integer(c_size_t) :: {name}_glcI_, {name}_glcJ_\n'
            )
            reassign = (
                f'allocate({name}_F_({name}_count))\n'
                f'allocate(character(len=int({name}_charLen)) :: {name}_buf_)\n'
                f'do {name}_glcI_ = 1, {name}_count\n'
                f'  do {name}_glcJ_ = 1, {name}_charLen\n'
                f'    {name}_buf_({name}_glcJ_:{name}_glcJ_)'
                f' = {name}(({name}_glcI_-1)*{name}_charLen + {name}_glcJ_)\n'
                f'  end do\n'
                f'  {name}_F_({name}_glcI_) = trim({name}_buf_)\n'
                f'end do\n'
            )
            if is_optional:
                reassign = (
                    f'if (present({name})) then\n'
                    f'   {reassign.replace(chr(10), chr(10) + "   ").rstrip()}\n'
                    f'end if\n'
                )
            arg.fort_reassignment = reassign
            arg.fort_pass_as      = f'{name}_F_'
            # ISO_Varying_String supplies both the type and the
            # `assignment(=)` overload from char to varying_string.
            arg.fort_modules.setdefault('ISO_Varying_String', {})['varying_string'] = 1
            arg.fort_modules.setdefault('ISO_Varying_String', {})['assignment(=)']  = 1
        elif arg.is_array and arg.polymorphic_list_array:
            # `type(<class>List), dimension(:)`.  bind(c) gets a flat
            # `type(c_ptr), dimension(*)` buffer plus a parallel
            # `integer(c_int), dimension(*)` IDs buffer plus a count.
            # Build a `type(<class>List), dimension(:), allocatable`
            # local of the right size, then for each element use the
            # registered class's `<class>GetPtr(ptr,id)` helper to
            # recover the polymorphic pointer (same dispatch the scalar
            # `class(<class>Class)` arg path uses) and stash it on the
            # element's `<class>_` component.
            stem      = arg.polymorphic_list_class
            list_type = arg.polymorphic_list_type
            comp      = arg.polymorphic_list_component
            arg.fort_declarations = (
                f'  type({list_type}), dimension(:), allocatable'
                f' :: {name}_F_\n'
                f'  integer(c_size_t) :: {name}_glcI_\n'
            )
            reassign = (
                f'allocate({name}_F_({name}_count))\n'
                f'do {name}_glcI_ = 1, {name}_count\n'
                f'  {name}_F_({name}_glcI_)%{comp} =>'
                f' {stem}GetPtr({name}({name}_glcI_), {name}_IDs({name}_glcI_))\n'
                f'end do\n'
            )
            if is_optional:
                reassign = (
                    f'if (present({name})) then\n'
                    f'   {reassign.replace(chr(10), chr(10) + "   ").rstrip()}\n'
                    f'end if\n'
                )
            arg.fort_reassignment   = reassign
            arg.fort_pass_as        = f'{name}_F_'
            # Re-use the scalar class-arg path's interface-block emission
            # (see fortran_declarations in Emitters.py) — the same
            # `<class>GetPtr` is what we're calling element-by-element.
            arg.fort_function_class = stem
            # The list type's home module is needed for the `use ::`
            # line that imports `<class>List` into the bind(c) wrapper.
            # Either an explicit override in _SHARED_TYPE_MODULES or a
            # walk of moduleUses is required; the override is checked
            # first (same scheme the scalar derived-type branch uses).
            list_mod = _SHARED_TYPE_MODULES.get(list_type)
            if not list_mod:
                list_mod = _module_for_symbol(
                    func_class.get('moduleUses', []), list_type)
            if implementation and not list_mod:
                cls = implementation['name']
                while cls and not list_mod:
                    list_mod = _module_for_symbol(
                        module_uses_impls.get(cls, []), list_type)
                    cls = extensions.get(cls)
            if list_mod:
                arg.fort_modules.setdefault(list_mod, {})[list_type] = 1
            # The bind(c) wrapper also needs `<class>Class` visible in
            # the host scope so the GetPtr interface block (emitted by
            # `fortran_declarations` from `fort_function_class`) can
            # `import` it — without this, gfortran reports "Cannot
            # IMPORT '<class>class' from host scoping unit ... does not
            # exist" because the host only sees the List type.  The
            # Class lives in the registered functionClass's home module
            # (looked up via `lib_function_classes`), which for the
            # current set of polymorphic-list-array users happens to be
            # the same module as the List type, but we don't rely on
            # that — Galacticus has classes whose List wrapper lives in
            # a separate module.
            class_type = stem + 'Class'
            class_mod  = lib_function_classes.get(stem, {}).get('module')
            if class_mod:
                arg.fort_modules.setdefault(class_mod, {})[class_type] = 1
        elif arg.is_array and arg.intrinsic == 'character':
            # Fixed-length character array.  bind(c) declares it as a
            # flat `character(c_char), dimension(*)` byte buffer (length
            # count*N); the inner Galacticus method wants
            # `character(len=N), dimension(:)`.  Repack into an
            # allocatable local of the right type, copying byte-by-byte
            # (an explicit do-loop is bulletproof across compilers, and
            # allocatable sidesteps any concerns about forward-reference
            # to the count dummy in an automatic-array size spec).  The
            # local is auto-deallocated when the wrapper returns.
            n = arg.char_len
            # Both loop indices are c_size_t to match the count
            # companion's kind — the outer loop runs to {name}_count
            # (c_size_t), and arithmetic mixing kinds in the index
            # expression `({name}_glcI_-1)*N + {name}_glcJ_` is cleanest
            # when both are the same kind.
            arg.fort_declarations = (
                f'  character(len={n}), dimension(:), allocatable'
                f' :: {name}_F_\n'
                f'  integer(c_size_t) :: {name}_glcI_, {name}_glcJ_\n'
            )
            reassign = (
                f'allocate({name}_F_({name}_count))\n'
                f'do {name}_glcI_ = 1, {name}_count\n'
                f'  do {name}_glcJ_ = 1, {n}\n'
                f'    {name}_F_({name}_glcI_)({name}_glcJ_:{name}_glcJ_)'
                f' = {name}(({name}_glcI_-1)*{n} + {name}_glcJ_)\n'
                f'  end do\n'
                f'end do\n'
            )
            if is_optional:
                # See the 2D branch above for the rationale: optional
                # arrays may be a NULL pointer at bind(c), and indexing
                # them would segfault.
                reassign = (
                    f'if (present({name})) then\n'
                    f'   {reassign.replace(chr(10), chr(10) + "   ").rstrip()}\n'
                    f'end if\n'
                )
            arg.fort_reassignment = reassign
            arg.fort_pass_as = f'{name}_F_'
        elif arg.is_array and arg.intrinsic == 'logical':
            # 1D deferred-shape logical array.  bind(c) signature has
            # `logical(c_bool), dimension(*)` plus a c_size_t count; the
            # inner method's dummy is the default logical kind (usually
            # 4 bytes, vs c_bool's 1 byte), so we can't just slice the
            # raw buffer.  Allocate a local of the right kind, copy
            # element-wise via the elemental `logical()` cast, and pass
            # the local to the inner.  Allocatable sidesteps any concern
            # about forward-reference to the count dummy in an
            # automatic-array size spec; the local is auto-deallocated
            # on wrapper return.  Optionals gate the conversion on
            # `present(name)` — when absent the bind(c) buffer is a
            # NULL pointer and indexing it would segfault.  The inner
            # call's `present(...)` check is handled at the fortran_call_code
            # level (which omits absent optionals from the call), so
            # there's no need to also gate the local's allocation here.
            arg.fort_declarations = (
                f'logical, dimension(:), allocatable :: {name}_F_\n'
            )
            reassign = (
                f'allocate({name}_F_({name}_count))\n'
                f'{name}_F_ = logical({name}(1:{name}_count))\n'
            )
            if is_optional:
                reassign = (
                    f'if (present({name})) then\n'
                    f'   {reassign.replace(chr(10), chr(10) + "   ").rstrip()}\n'
                    f'end if\n'
                )
            arg.fort_reassignment = reassign
            arg.fort_pass_as      = f'{name}_F_'
        elif arg.is_array:
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
                        import_module = _module_for_symbol(
                            module_uses_impls.get(cls, []), type_spec_val)
                        cls = extensions.get(cls)
                # 2. fall back to the functionClass file's own module uses.
                if not import_module:
                    import_module = _module_for_symbol(
                        func_class.get('moduleUses', []), type_spec_val)
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
            # Two shapes land here.  The simple case is `class(FooClass)`:
            #   class(FooClass), pointer :: name_
            #   name_ => FooGetPtr(name, name_ID)
            #
            # The abstract-intermediate case (arg.narrowing_type set) is
            # `class(<intermediate>)` where <intermediate>'s extends-chain
            # reaches a registered <base>Class.  The inner constructor
            # wants a `class(<intermediate>), pointer`, so we recover the
            # base pointer from <base>GetPtr and narrow with `select type`:
            #   class(<base>Class    ), pointer :: name_base_
            #   class(<intermediate> ), pointer :: name_
            #   name_base_ => <base>GetPtr(name, name_ID)
            #   select type (name_base_)
            #   class is (<intermediate>)
            #      name_ => name_base_
            #   class default
            #      call Error_Report('argument <name> is not of type <intermediate>'//{introspection:location})
            #   end select
            # Error_Report aborts with a clear diagnostic; the introspection
            # tag is expanded later by the SourceTree introspection pass
            # (same convention used by the GetPtr classID-default branch).
            #
            # assign_c_types stashes the registered <base> on
            # arg.fort_function_class (the abstract-intermediate path
            # can't recover it from type_spec alone — the hierarchy walk
            # is the only place that knows).  Fall back to the
            # strip-trailing-Class convention for the plain case so
            # direct emitter callers (e.g. unit tests) keep working
            # without first running through assign_c_types.
            class_key = arg.fort_function_class or (
                type_spec_val[:-5] if type_spec_val.endswith('Class')
                else type_spec_val
            )
            mod_name  = lib_function_classes.get(class_key, {}).get('module', '')
            if arg.narrowing_type:
                arg.fort_declarations = (
                    f'class({class_key}Class), pointer :: {name}_base_\n'
                    f'class({arg.narrowing_type}), pointer :: {name}_\n'
                )
                arg.fort_pass_as = name + '_'
                narrow_body = (
                    f'{name}_base_ => {class_key}GetPtr({name},{name}_ID)\n'
                    f'select type ({name}_base_)\n'
                    f'class is ({arg.narrowing_type})\n'
                    f'   {name}_ => {name}_base_\n'
                    f'class default\n'
                    f"   call Error_Report('argument ''{name}'' is not of "
                    f"type {arg.narrowing_type}'//{{introspection:location}})\n"
                    f'end select\n'
                )
                if arg.is_optional:
                    # Optional: skip the narrowing when the ID companion
                    # is absent, leaving name_ as null() so the inner call's
                    # `present()` check sees the same answer.
                    arg.fort_reassignment = (
                        f'if (present({name}_ID)) then\n {narrow_body}'
                        f'else\n {name}_ => null()\nend if\n'
                    )
                else:
                    arg.fort_reassignment = narrow_body
                arg.fort_modules.setdefault('Error', {})['Error_Report'] = 1
                if mod_name:
                    arg.fort_modules.setdefault(mod_name, {})[f'{class_key}Class'] = 1
                    arg.fort_modules.setdefault(mod_name, {})[arg.narrowing_type] = 1
            else:
                arg.fort_declarations  = f'class({type_spec_val}), pointer :: {name}_\n'
                arg.fort_pass_as       = name + '_'
                arg.fort_reassignment  = (
                    f'{opt_prefix}{name}_ => {class_key}GetPtr({name},{name}_ID)\n'
                )
                if mod_name:
                    arg.fort_modules.setdefault(mod_name, {})[type_spec_val] = 1
            arg.fort_function_class = class_key

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
