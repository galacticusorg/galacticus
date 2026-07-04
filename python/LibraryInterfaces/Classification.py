"""Shared classification of constructor/method arguments and return types
for the library-interface pipeline.

This module is the single owner of the rules deciding which Fortran
argument and return types the libgalacticus wrapper generator can
translate. Two consumers use it:

* ``scripts/build/libraryInterfaces.py`` (the generator) — via
  :func:`unsupported_arg`, which returns a reject reason (or ``None``) for
  each argument, and the return-type recognizer regexes used by its
  ``interfaces_methods`` dispatch.
* ``scripts/build/libraryInterfacesAudit.py`` (the manual audit tool) — via
  :func:`classify_arg` in audit mode, which additionally distinguishes
  "blocked by the pipeline" from "would work if another functionClass were
  registered" (a missing dependency).

Historically the audit carried hand-synced copies of these rules marked
"MUST stay in sync", and they drifted (e.g. the generator's rejection of
non-optional ``integer(omp_lock_kind)`` args was missing from the audit,
which therefore over-reported readiness). Keeping the predicate here makes
that class of drift impossible.
"""

import re

from LibraryInterfaces.Hierarchy import resolve_function_class_base
from LibraryInterfaces.Pipeline import (_SHARED_TYPE_MODULES,
                                        _CALLBACK_PROCEDURE_INTERFACES,
                                        _DIM_IDENT_EXTENTS_RX)

__all__ = [
    'ENUM_RETURN_RX', 'CLASS_RETURN_RX', 'ARRAY_RETURN_RX',
    'DYNAMIC_ARRAY_RETURN_RX', 'DYNAMIC_ARRAY_RETURN_2D_RX',
    'DIM_FIXED_RX', 'CHAR_LEN_RX',
    'RETURN_TYPE_ALIASES', 'SCALAR_RETURN_OK', 'OUTPUT_ARRAY_RETURN_OK',
    'normalize_method_return_type', 'is_internal_constructor_name',
    'classify_arg', 'unsupported_arg',
    'is_output_array_arg', 'is_output_scalar_arg',
    'is_output_sized_array_arg', 'output_sized_extents',
    'unsupported_output_array_method',
]


# ---------------------------------------------------------------------------
# Return-type recognizers (used by the generator's interfaces_methods switch
# and by the audit's method classification).
# ---------------------------------------------------------------------------

# Match `type(enumerationFooType)` — the Galacticus convention for an
# enumeration kind named "foo".
ENUM_RETURN_RX = re.compile(
    r'^type\s*\(\s*(enumeration[a-z0-9_]+type)\s*\)$',
    re.IGNORECASE,
)

# Match `class(fooClass)` — a polymorphic pointer to a registered
# functionClass; routed through the per-class XGetIdAndPtr helper.
CLASS_RETURN_RX = re.compile(
    r'^class\s*\(\s*([a-z][a-zA-Z0-9_]*Class)\s*\)$',
    re.IGNORECASE,
)

# Match a 1D fixed-size array RETURN type, e.g. `double precision,
# dimension(3)`. Captures (intrinsic, size).
ARRAY_RETURN_RX = re.compile(
    r'^(double\s+precision|integer)\s*,\s*dimension\s*\(\s*(\d+)\s*\)\s*$',
    re.IGNORECASE,
)

# Match a 1D *dynamic-size* numeric array RETURN type — anything with a
# `dimension(:)`-like attribute whose extent isn't a literal integer.
# Two source-level shapes share this codegen template:
#
#   * `double precision, allocatable, dimension(:)` — inner method returns
#     an allocatable array; captured into a save-target buffer on
#     assignment, with c_loc + size conveyed back to Python.
#   * `double precision, dimension(self%X)` / `dimension(size(Y))` — a
#     fixed-shape array whose extent is computed at run time. Fortran 2003
#     auto-realloc on assignment to an allocatable LHS handles both shapes
#     identically.
#
# The leading lookahead skips literal-integer extents (the fixed-size case
# ARRAY_RETURN_RX handles). Captures (intrinsic, extent-expression); the
# extent expression is diagnostic only.
DYNAMIC_ARRAY_RETURN_RX = re.compile(
    r'^(double\s+precision|integer)'
    r'(?:\s*,\s*allocatable)?'
    r'\s*,\s*dimension\s*\(\s*'
    r'(?!\s*\d+\s*\))'                    # exclude `dimension(<literal>)`
    r'([^,]+)\s*\)\s*$',                  # forbid top-level commas — 1D
                                          # only — but allow internal parens
                                          # so `dimension(size(X))`,
                                          # `dimension(self%X)` etc. match.
    re.IGNORECASE,
)

# 2D variant of the dynamic/allocatable array return — same save-buffer
# codegen plus a second `c_size_t` size companion so Python can reshape the
# flat byte buffer with Fortran's column-major layout. Captures
# (intrinsic, extent1, extent2); extents are diagnostic only.
DYNAMIC_ARRAY_RETURN_2D_RX = re.compile(
    r'^(double\s+precision|integer)'
    r'(?:\s*,\s*allocatable)?'
    r'\s*,\s*dimension\s*\(\s*'
    r'([^,]+)\s*,\s*([^,]+)'
    r'\s*\)\s*$',
    re.IGNORECASE,
)

# Match a fixed-size dimension(N[,M...]) attribute (integer literals, with
# optional explicit lower bounds `L:U`).
DIM_FIXED_RX = re.compile(
    r'^dimension\s*\(\s*'
    r'(?:\d+\s*:\s*)?\d+'
    r'(?:\s*,\s*(?:\d+\s*:\s*)?\d+)*'
    r'\s*\)$'
)

# Match `len=N` literal in a character type-spec. Variable-length forms
# (`len=*`, `len=:`) are unsupported (no fixed stride at the byte
# boundary).
CHAR_LEN_RX = re.compile(r'^len\s*=\s*(\d+)$')

# Aliases the return-type switch accepts as synonyms. Galacticus source
# uses several spellings for the same C-interop integer kind, and the
# switch only has branches for the unprefixed forms. `kind_int8 =
# SELECTED_INT_KIND(18)` (from Kind_Numbers) resolves to the same width as
# `c_long` on Linux/macOS, and `c_long` is what assign_c_types already
# uses for the arg-side spelling, so the calling convention stays
# consistent across constructor args and method returns.
RETURN_TYPE_ALIASES = {
    'integer(kind=c_size_t)' : 'integer(c_size_t)',
    'integer(kind=c_long)'   : 'integer(c_long)',
    'integer(kind=kind_int8)': 'integer(c_long)',
    'integer(kind_int8)'     : 'integer(c_long)',
}

# Scalar return types the generator handles outright (no per-type plumbing
# beyond ISO_C_Binding kind selection). Used by the audit; must match the
# branches of the generator's interfaces_methods switch (which cannot yet
# be derived mechanically from that 500-line dispatch).
SCALAR_RETURN_OK = frozenset({
    'void', 'double precision', 'integer', 'integer(c_long)',
    'integer(c_size_t)', 'logical', 'type(varying_string)',
})

# Return types permitted on a method that also has output-array arguments:
# the direct-restype scalars — those the bind(c) wrapper returns as a plain
# function result with a matching ctypes restype, no subroutine lowering or
# conversion buffer. Return types that lower the wrapper to a subroutine
# with synthetic out-args (varying_string, class(...), array returns) would
# collide with the output-array companion protocol, so they stay blocked
# alongside output arrays until a method needs the combination.
OUTPUT_ARRAY_RETURN_OK = frozenset({
    'void', 'double precision', 'integer', 'integer(c_long)',
    'integer(c_size_t)', 'logical',
})


def normalize_method_return_type(ret_type):
    """Collapse equivalent return-type spellings to the canonical form the
    generator's return-type switch handles natively. Anything not in
    :data:`RETURN_TYPE_ALIASES` passes through unchanged."""
    return RETURN_TYPE_ALIASES.get(ret_type.strip(), ret_type)


def is_internal_constructor_name(name):
    """Recognise an Internal-suffixed module-procedure name following the
    Galacticus convention ``<short>Constructor[Internal[Suffix]]`` or the
    alternative ``<short>Internal`` form (used by the merger-tree walkers,
    e.g. ``allAndFormationNodesInternal``). Accepts either:

    * a name ending in "internal" (catches both
      ``<short>ConstructorInternal`` and ``<short>Internal``), or
    * a name containing "constructorinternal" (catches the rarer
      ``ConstructorInternalType`` / ``ConstructorInternalDefined``
      variants used to disambiguate multiple internal constructors).

    ``<name>ConstructorParameters`` and similar satisfy neither rule and
    are correctly rejected.
    """
    lower = name.lower()
    return lower.endswith('internal') or 'constructorinternal' in lower


# ---------------------------------------------------------------------------
# Argument classification
# ---------------------------------------------------------------------------

def classify_arg(arg, registered_classes, *, constructor_overrides=(),
                 class_hierarchy=None, known_function_classes=None):
    """Classify one constructor/method argument.

    Returns ``None`` when the pipeline supports the argument;
    ``('blocked', reason)`` when it cannot be translated; or — in audit
    mode only — ``('missing-dep', {stems})`` when the argument would be
    supported if the named functionClasses were also registered in
    libraryClasses.xml.

    `registered_classes` is the mapping (or set) of functionClass stems
    registered in libraryClasses.xml. `known_function_classes` enables
    audit mode: the set of ALL functionClass stems that exist in the source
    tree, used to distinguish an unregistered-but-real functionClass
    dependency (missing-dep) from a type the pipeline can never handle
    (blocked). Generator mode (``known_function_classes=None``) reports
    both as blocked, since the generator must skip the implementation
    either way.

    Rejection rules (generator ground truth):

    * ``complex`` / ``double complex`` — ctypes has no built-in c_complex.
    * ``procedure(...)`` — no ctypes counterpart; without an early reject
      the arg falls through assign_c_types with empty ctype/fort_type and
      the bind(c) wrapper emits a broken declaration.
    * non-optional ``integer(omp_lock_kind)`` — platform-dependent kind
      (INTEGER(4) on Linux GCC vs INTEGER(8) on macOS GCC), no portable
      C-interop kind exists. Optional ones are dropped silently by
      assign_c_types, letting the inner method's default take over.
    * ``class(FooClass)`` whose stem isn't registered (generator mode) —
      without registration assign_c_types can't set is_function_class, so
      no ``_ID`` companion or GetPtr reassignment is generated.
    * ``class(SomethingElse)`` that is not a functionClass and whose
      extends-chain doesn't reach a registered base.
    * ``class(*)`` without a libraryClasses.xml ``<argument>`` override.
    * Array shapes outside the supported set (see the dimension block).

    `constructor_overrides` entries from libraryClasses.xml modify these
    rules: ``value='null'`` accepts procedure / class(*) args by passing a
    local null pointer; ``value='absent'`` drops an *optional* arg
    entirely; a plain ``<argument name=… type=… module=…/>`` hint resolves
    ``class(*)``.
    """
    audit_mode  = known_function_classes is not None
    missing     = set()
    intrinsic   = arg.get('intrinsic')
    # A `<argument name="..." value="null"/>` override tells the wrapper to
    # drop the arg and pass a local null pointer to the inner constructor —
    # sufficient for the procedure-pointer + class(*) callback-injection
    # escape hatch (see assign_c_types / build_fortran_reassignments).
    if any(isinstance(o, dict)
           and o.get('name') == arg.get('name')
           and o.get('value') == 'null'
           for o in constructor_overrides):
        type_spec = (arg.get('type') or '').strip()
        if intrinsic == 'procedure' \
                or (intrinsic == 'class' and type_spec == '*'):
            return None
        return ('blocked',
                f"value='null' override on {intrinsic}({type_spec}) — "
                f"only procedure and class(*) args are supported")
    # A `<argument name="..." value="absent"/>` override drops an
    # *optional* arg entirely — from Python, from the bind(c) signature,
    # and from the inner constructor call.
    if any(isinstance(o, dict)
           and o.get('name') == arg.get('name')
           and o.get('value') == 'absent'
           for o in constructor_overrides):
        if 'optional' in arg.get('attributes', []):
            return None
        return ('blocked',
                f"value='absent' override on non-optional argument — only "
                f"optional args may be dropped")
    # Sized output buffers (`intent(out), dimension(<argName>,…)`) are
    # accepted — including complex(c_double_complex) elements — BEFORE the
    # blanket complex reject below; see is_output_sized_array_arg.  Their
    # whole-method requirements (extents name integer intent(in) args,
    # etc.) are enforced by unsupported_output_array_method.
    if is_output_sized_array_arg(arg):
        return None
    if intrinsic in ('complex', 'double complex'):
        return ('blocked', f"{intrinsic}({arg.get('type','')})")
    if intrinsic == 'procedure':
        type_spec = (arg.get('type') or '').strip()
        attrs     = arg.get('attributes', [])
        # `procedure(...), pointer, intent(out|inout)` — the method hands a
        # *Fortran procedure pointer* back to the caller.  There is nothing
        # a Python caller could do with one, so this is unsupportable in
        # principle (not merely unimplemented); the audit buckets this
        # message as out-of-scope so it doesn't sit in the actionable
        # worklist.  (These methods — timeEvolveTo, evolve,
        # differentialEvolution — also all take treeNode/mergerTree args.)
        if 'pointer' in attrs and ('intent(out)' in attrs
                                   or 'intent(inout)' in attrs):
            return ('blocked',
                    f"procedure({type_spec}) pointer output — a Fortran "
                    f"procedure pointer returned to the caller cannot be "
                    f"exposed to Python")
        # Inbound callback with a registered, hand-checked marshalling
        # recipe (see Pipeline._CALLBACK_PROCEDURE_INTERFACES): the wrapper
        # accepts a Python callable via CFUNCTYPE → c_funptr → a Fortran
        # shim adapting the Galacticus-side interface.  Both plain and
        # `pointer` dummies qualify — for a pointer dummy the wrapper
        # passes the shim module's procedure-pointer slot aimed at the shim
        # (a pointer actual for a pointer dummy; if the callee repoints it,
        # the slot is re-aimed on the next call).  Pointer dummies with
        # intent(out|inout) were already rejected above as outputs;
        # optional callbacks are not yet supported.
        if (type_spec in _CALLBACK_PROCEDURE_INTERFACES
                and 'optional' not in attrs):
            return None
        return ('blocked',
                f"procedure({type_spec}) — procedure-pointer args "
                f"are not supported")
    # Object-pointer dummies with intent(out|inout): the callee may (re)point
    # the pointer, and the wrapper's local-pointer passing loses that
    # repointing silently — e.g. a tree-walker's `next(node)` would return
    # success but never advance the caller's node.  Two flavours:
    #
    # * `class(...), pointer` outputs (only `class(*)` exists in the tree:
    #   timeEvolveTo's taskSelf) — an unlimited-polymorphic pointer handed
    #   back to the caller is unsupportable in principle; the "pointer
    #   output" wording routes it to the audit's out-of-scope bucket.
    # * `type(X), pointer` in/outs (mergerTreeWalker.next, nodeEvolver's
    #   promote, buildController.control, …) — supportable via a pointer
    #   write-back protocol (bind(c) takes the c_ptr by reference and
    #   writes back c_loc of the repointed target); until that exists the
    #   wrapper is silently broken, so reject.  The wording deliberately
    #   avoids the literal "type(" so the audit's internal-derived-type
    #   rule doesn't misfile this *actionable* blocker as deferred.
    pointer_attrs = arg.get('attributes', [])
    if intrinsic in ('class', 'type') and 'pointer' in pointer_attrs \
            and ('intent(out)' in pointer_attrs
                 or 'intent(inout)' in pointer_attrs):
        type_spec = (arg.get('type') or '').strip()
        if intrinsic == 'class':
            return ('blocked',
                    f"class({type_spec}) pointer output — a Fortran object "
                    f"pointer returned to the caller cannot be exposed to "
                    f"Python")
        return ('blocked',
                f"pointer dummy of derived {intrinsic} '{type_spec}' — "
                f"repointing by the method would be silently lost by the "
                f"wrapper (needs a pointer write-back protocol)")
    if intrinsic == 'integer' \
            and (arg.get('type') or '').strip() == 'omp_lock_kind' \
            and 'optional' not in arg.get('attributes', []):
        return ('blocked',
                "integer(omp_lock_kind) non-optional argument — kind is "
                "platform-dependent (INTEGER(4) on Linux vs INTEGER(8) on "
                "macOS), no portable C-interop kind exists")
    if intrinsic == 'class':
        type_spec = (arg.get('type') or '').strip()
        if type_spec == '*':
            has_override = any(
                isinstance(o, dict) and o.get('name') == arg.get('name')
                for o in constructor_overrides
            )
            if not has_override:
                return ('blocked',
                        "class(*) without a libraryClasses.xml override "
                        "(<constructor><argument name=… type=… module=…/>)")
        elif type_spec.endswith('Class'):
            stem = type_spec[:-5]
            if stem not in registered_classes:
                if audit_mode and stem in known_function_classes:
                    missing.add(stem)
                else:
                    return ('blocked',
                            f"class({type_spec}) — '{stem}' is not a "
                            f"registered functionClass in libraryClasses.xml")
        else:
            # Abstract intermediate: a class whose extends-chain reaches a
            # <base>Class is acceptable — the pipeline routes through
            # <base>GetPtr and narrows to <intermediate> at run time. In
            # audit mode resolve against every functionClass that exists
            # (an unregistered base is a missing dependency); in generator
            # mode only a registered base will do.
            base = None
            if class_hierarchy:
                resolve_set = (set(known_function_classes) if audit_mode
                               else set(registered_classes))
                base, _ = resolve_function_class_base(
                    type_spec, class_hierarchy, resolve_set)
            if base is None:
                return ('blocked',
                        f"class({type_spec}) — only registered "
                        f"functionClasses and class(*) (with override) are "
                        f"supported")
            if audit_mode and base not in registered_classes:
                missing.add(base)
    # Dimension check — same predicate for constructor and method args.
    attrs = arg.get('attributes', [])
    for attr in attrs:
        if attr.startswith('dimension'):
            # 1D numeric arrays are plumbed through with a numpy ndarray
            # conversion at the Python boundary; deferred-shape gets a
            # count companion, fixed-size uses the literal length from the
            # dimension spec (see assign_c_types /
            # build_python_reassignments in Pipeline.py).
            is_supported_dim = (
                intrinsic in ('double precision', 'integer')
                and (attr == 'dimension(:)'
                     or attr == 'dimension(:,:)'
                     or DIM_FIXED_RX.match(attr))
            )
            # 1D deferred-shape logical: same plumbing as numeric 1D
            # arrays, except build_fortran_reassignments emits a
            # kind-narrowing copy from `logical(c_bool)` to the default
            # `logical` kind the inner method's dummy declares.
            is_supported_logical = (
                intrinsic == 'logical' and attr == 'dimension(:)'
            )
            # Fixed-length character arrays (`character(len=N),
            # dimension(:)`): the bind(c) boundary receives a contiguous
            # count*N byte buffer plus a count companion.
            is_supported_char_array = (
                intrinsic == 'character'
                and attr == 'dimension(:)'
                and CHAR_LEN_RX.match((arg.get('type') or '').strip())
            )
            # `type(varying_string), dimension(:)`: shipped via a flat
            # contiguous byte buffer plus count + per-element-length
            # companions.
            type_spec = (arg.get('type') or '').strip()
            is_supported_vstring_array = (
                intrinsic == 'type'
                and type_spec == 'varying_string'
                and attr == 'dimension(:)'
            )
            # `type(<class>List), dimension(:)`: an array of polymorphic
            # pointers — supported when `<class>` is a registered
            # functionClass AND the wrapper type is registered in
            # _SHARED_TYPE_MODULES (Pipeline.py's gatekeeper); other
            # locally-defined `<class>List` structs sharing the naming
            # convention fall through to the generic rejection. In audit
            # mode an unregistered-but-real stem is a missing dependency.
            list_stem = type_spec[:-4] if type_spec.endswith('List') else None
            is_supported_list_array = (
                intrinsic == 'type'
                and attr == 'dimension(:)'
                and list_stem is not None
                and type_spec in _SHARED_TYPE_MODULES
                and (list_stem in registered_classes
                     or (audit_mode and list_stem in known_function_classes))
            )
            if (is_supported_list_array and audit_mode
                    and list_stem not in registered_classes):
                missing.add(list_stem)
            if (is_supported_dim or is_supported_logical
                    or is_supported_char_array
                    or is_supported_vstring_array or is_supported_list_array):
                if 'allocatable' in attrs:
                    # 1D `intent(out), allocatable, dimension(:)` numeric
                    # args are supported as OUTPUT arrays: the inner method
                    # allocates and fills them, and the data + size flow
                    # back to Python (see ArgSpec.is_output_array). The
                    # surrounding method must additionally satisfy
                    # `unsupported_output_array_method` — a whole-method
                    # property (void return, no optionals, other args are
                    # inputs) that the generator and audit enforce, not a
                    # per-arg one. Every other allocatable shape
                    # (intent(in)/inout, non-numeric, multi-dim) stays
                    # blocked.
                    if is_output_array_arg(arg):
                        continue
                    return ('blocked',
                            '1D allocatable array argument'
                            ' (output arrays not yet supported)')
                # `intent(inout)` / `intent(out)` non-allocatable arrays
                # are accepted via the in-place-mutable-buffer path: the
                # inner Galacticus method mutates the caller's contiguous
                # numpy buffer in place.
                continue
            # Arrays of polymorphic objects (`class(...), dimension(:)`)
            # cannot be assembled from the per-object pointers a Python
            # caller holds: a polymorphic array has ONE dynamic type, and
            # intrinsic assignment into polymorphic array elements is not
            # permitted.  Distinct message so the audit defers these
            # rather than listing them as actionable.
            if intrinsic == 'class':
                return ('blocked',
                        f"class({type_spec}) array argument — arrays of "
                        f"polymorphic objects cannot be assembled from "
                        f"Python-held object pointers")
            # Name the element type in the generic reject so the audit's
            # bucketing can tell `type(nBodyData)` (an internal derived
            # type — deferred) apart from numeric shapes that are merely
            # not-yet-supported (actionable).
            subject = (f"{intrinsic}({type_spec})" if type_spec
                       else f"{intrinsic}")
            return ('blocked',
                    f'dimensioned argument of {subject}'
                    ' (only 1D deferred-shape or fixed-size numeric input,'
                    ' 2D deferred-shape numeric input, or 1D deferred-shape'
                    ' fixed-length character arrays, are supported)')
    if missing:
        return ('missing-dep', missing)
    return None


def unsupported_arg(arg, lib_function_classes, *,
                    constructor_overrides=(), class_hierarchy=None):
    """Generator-mode wrapper around :func:`classify_arg`: return a
    human-readable reason if ``arg`` has a type the pipeline can't
    translate, otherwise ``None``."""
    verdict = classify_arg(
        arg, lib_function_classes,
        constructor_overrides=constructor_overrides,
        class_hierarchy=class_hierarchy,
    )
    if verdict is None:
        return None
    return verdict[1]


# ---------------------------------------------------------------------------
# Output-array arguments (`intent(out), allocatable, dimension(:)`)
# ---------------------------------------------------------------------------

def is_output_array_arg(arg):
    """True if *arg* is a 1D ``allocatable, dimension(:)`` numeric
    (``double precision`` / ``integer``) OUTPUT-array argument — either
    ``intent(out)`` or ``intent(inout)`` — the shape the generator lowers to
    a ``save, target`` buffer with a ``(c_ptr, c_size_t)`` intent(out)
    companion pair (see :attr:`ArgSpec.is_output_array`).

    ``intent(inout)`` is accepted because the Galacticus idiom for such args
    is allocation-*reuse*, not input: the inner method checks ``allocated()``
    / ``size()`` and reallocates as needed, then overwrites the contents
    (e.g. ``computationalDomain.indicesFromPosition``), so it never reads an
    incoming value. The wrapper deallocates its ``save`` buffer before each
    call, so the inner always sees an unallocated array and allocates a fresh
    one — identical to the ``intent(out)`` path. (A hypothetical method that
    genuinely read an ``intent(inout)`` allocatable before allocating it
    could not be driven output-only, but that is not the convention and such
    a method would already be ill-defined for an unallocated actual.)

    An ``optional`` output array (e.g. ``timeLightconeCrossing``'s
    ``timesCrossing``) also qualifies: the wrapper passes its buffer
    unconditionally, so the inner method always sees the dummy as present
    and always fills it — the Python caller unconditionally receives the
    array, which is the only sensible contract for a Python API (there is
    no way to "omit" a return value).

    Shapes: 1D and 2D numeric (``dimension(:)`` / ``dimension(:,:)``; 2D
    gets a second size companion and a column-major reshape on the Python
    side), and 1D ``logical`` (the wrapper copies the inner method's
    default-kind result into a ``logical(c_bool)`` export buffer, mirroring
    the kind-narrowing the logical *input* path already does).
    """
    attrs = arg.get('attributes', [])
    if 'allocatable' not in attrs:
        return False
    if not ('intent(out)' in attrs or 'intent(inout)' in attrs):
        return False
    if arg.get('intrinsic') in ('double precision', 'integer'):
        return 'dimension(:)' in attrs or 'dimension(:,:)' in attrs
    if arg.get('intrinsic') == 'logical':
        return 'dimension(:)' in attrs
    return False


def is_output_scalar_arg(arg):
    """True if *arg* is a scalar ``intent(out)`` numeric/logical companion —
    an ``intent(out)`` argument that is not allocatable and has no
    ``dimension`` (e.g. ``integer, intent(out) :: count`` or ``double
    precision, intent(out) :: f``).

    These accompany output arrays on methods like
    ``accretionDiskSpectra.wavelengths`` (a count beside the array) or
    ``variogram.modelFDF`` (a scalar residual beside the gradient array).
    The generator passes them straight through as by-reference intent(out)
    scalars and appends their filled values to the Python return. Derived
    types, ``class(...)``, and ``character`` scalars are excluded — those
    need their own handling.
    """
    attrs = arg.get('attributes', [])
    if 'intent(out)' not in attrs:
        return False
    if 'allocatable' in attrs:
        return False
    if any(a.startswith('dimension') for a in attrs):
        return False
    return arg.get('intrinsic') in (
        'double precision', 'real', 'integer', 'logical')


def is_output_sized_array_arg(arg):
    """True if *arg* is an ``intent(out)`` explicit-shape array whose every
    extent is a plain identifier (``dimension(gridCount,…)``) — a "sized
    output buffer".  Because the extents are other (integer, intent(in))
    arguments of the same method, the Python wrapper can pre-allocate a
    buffer of the right size, pass it in place (no companions needed), and
    return it reshaped column-major.  Elements may be numeric or
    ``complex(c_double_complex)`` (bit-compatible with numpy complex128;
    the canonical case is surveyGeometry.windowFunctions' FFT grids).

    Extent-name validity (each names an integer ``intent(in)`` arg) is a
    whole-method property checked by
    :func:`unsupported_output_array_method`.
    """
    attrs = arg.get('attributes', [])
    if 'intent(out)' not in attrs:
        return False
    if 'allocatable' in attrs or 'pointer' in attrs:
        return False
    intr = arg.get('intrinsic')
    ts   = (arg.get('type') or '').strip()
    if not (intr in ('double precision', 'integer')
            or (intr == 'complex' and ts == 'c_double_complex')):
        return False
    return any(_DIM_IDENT_EXTENTS_RX.match(a) for a in attrs
               if a.startswith('dimension'))


def output_sized_extents(arg):
    """Return the list of extent identifiers of a sized output buffer's
    ``dimension(...)`` attribute (see :func:`is_output_sized_array_arg`)."""
    for a in arg.get('attributes', []):
        m = _DIM_IDENT_EXTENTS_RX.match(a) if a.startswith('dimension') \
            else None
        if m:
            return [e.strip() for e in m.group(1).split(',')]
    return []


def unsupported_output_array_method(args, return_type):
    """Whole-method gate for the output-array feature.

    Returns a human-readable reason when *args* contains at least one
    output-array argument (see :func:`is_output_array_arg`) but the
    surrounding method has a feature the current output-array codegen can't
    handle; otherwise ``None`` — including when there are no output arrays at
    all, which is the normal (non-output-array) path.

    Restriction: an output-array method must return a direct-restype scalar
    (:data:`OUTPUT_ARRAY_RETURN_OK` — ``void`` or a plain numeric/logical
    the wrapper returns as an ordinary bind(c) function result), have no
    optional arguments, and every non-``self`` argument must be one of:

    * a supported input,
    * an output array or scalar ``intent(out)`` companion (returned),
    * an ``intent(out|inout)`` *dimensioned non-allocatable* array — the
      in-place mutable-buffer path: the Python caller supplies a numpy
      buffer the method fills, exactly as on non-output-array methods,
    * an ``intent(inout)`` ``class``/``type`` argument — pointer-passed;
      mutation happens on the heap object the caller already holds.

    What stays blocked: return types that lower the wrapper to a subroutine
    (varying_string / class / array returns — they collide with the
    companion protocol), optional args, and scalar ``intent(inout)``
    numeric/logical or ``intent(out|inout)`` character args (outputs the
    wrapper would silently drop).

    Shared by the generator (:mod:`libraryInterfaces`) and the audit
    (:mod:`libraryInterfacesAudit`) so their verdicts can't drift.

    *args* are the raw declaration dicts for the method's arguments
    EXCLUDING ``self``. *return_type* is the method's normalized return-type
    string (``'void'`` for subroutines).
    """
    if not any(is_output_array_arg(a) or is_output_sized_array_arg(a)
               for a in args):
        return None
    # Sized output buffers: every extent identifier must name an integer
    # intent(in) argument of the same method — that's what lets the Python
    # wrapper compute the allocation size before the call.  Compared
    # case-insensitively: Fortran names are case-insensitive and the
    # declaration parser lowercases attribute text (extents) while
    # variableNames preserve source case.
    int_inputs = {(a.get('name') or '').lower() for a in args
                  if a.get('intrinsic') == 'integer'
                  and 'intent(in)' in a.get('attributes', [])}
    for a in args:
        if not is_output_sized_array_arg(a):
            continue
        for extent in output_sized_extents(a):
            if extent.lower() not in int_inputs:
                return (f"sized output buffer '{a.get('name', '?')}' has "
                        f"extent '{extent}' that is not an integer "
                        f"intent(in) argument of the method")
    ret = normalize_method_return_type(return_type or 'void')
    # Dynamic-array returns (1D/2D allocatable or runtime-extent) lower the
    # wrapper to a subroutine whose own (c_ptr, size…) companions are
    # appended BEFORE the per-argument output-array companions, so the two
    # protocols compose mechanically; the Python wrapper leads the returned
    # tuple with the reshaped result array.
    if (ret not in OUTPUT_ARRAY_RETURN_OK
            and not DYNAMIC_ARRAY_RETURN_RX.match(ret)
            and not DYNAMIC_ARRAY_RETURN_2D_RX.match(ret)):
        return (f"output-array method with non-scalar return type "
                f"({return_type}) — only void, direct-scalar, or "
                f"dynamic-array returns may accompany output arrays")
    for a in args:
        if (is_output_array_arg(a) or is_output_scalar_arg(a)
                or is_output_sized_array_arg(a)):
            continue
        attrs = a.get('attributes', [])
        if 'optional' in attrs:
            return (f"output-array method with optional argument "
                    f"'{a.get('name', '?')}' — not yet supported")
        if 'intent(out)' in attrs or 'intent(inout)' in attrs:
            # Dimensioned non-allocatable → the in-place mutable-buffer
            # path (allocatable dimensioned args were caught by
            # is_output_array_arg or, for non-numeric elements, are
            # rejected per-arg by classify_arg before this gate runs).
            if any(at.startswith('dimension') for at in attrs):
                continue
            # class/type scalars are pointer-passed; in-place mutation is
            # visible to the caller through the pointer it already holds.
            if a.get('intrinsic') in ('class', 'type'):
                continue
            return (f"output-array method with non-input argument "
                    f"'{a.get('name', '?')}' — only supported inputs, "
                    f"in-place buffers, and scalar intent(out) companions "
                    f"may accompany output arrays")
    return None
