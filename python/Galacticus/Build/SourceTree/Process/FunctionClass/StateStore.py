"""State store/restore code generation helpers for the functionClass
pipeline.

Provides `state_store_explicit_function`,
`generate_allocatable_state_store_code`, and `state_store_variables`.
The stateStorables data is taken as an explicit parameter.

Andrew Benson (ported to Python 2026)
"""

import re


from List.ExtraUtils                                         import as_array
from Galacticus.Build.StateStorables                         import (
    function_class_names    as _shared_function_class_names,
    function_class_instances as _shared_function_class_instances,
    is_self_storing_type, is_unstored_type,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    strip_variable_name, declaration_rank, variable_rank,
)


# ---------------------------------------------------------------------------
# stateStorables shape helpers
# ---------------------------------------------------------------------------

def _function_class_name_set(state_storables):
    return _shared_function_class_names(state_storables)


def _state_storables_by_type(state_storables):
    """Return `{type: entry}` for every entry in `stateStorables/stateStorables`."""
    raw = (state_storables or {}).get('stateStorables') or []
    if isinstance(raw, dict):
        raw = [raw]
    return {e.get('type'): e for e in raw
            if isinstance(e, dict) and 'type' in e}


def _function_class_instances(state_storables):
    return _shared_function_class_instances(state_storables)


# ---------------------------------------------------------------------------
# Explicit store/restore function support
# ---------------------------------------------------------------------------

def state_store_explicit_function(non_abstract_class):
    """Generate store/restore code + required modules for classes that
    declare a `<stateStore>` block with an explicit `store` / `restore`
    function in their directive.

    Returns `(input_code, output_code, modules_dict)`.
    """
    input_code  = ''
    output_code = ''
    modules     = {}

    ss = (non_abstract_class or {}).get('stateStore', {})
    inner = ss.get('stateStore') if isinstance(ss, dict) else None
    if not isinstance(inner, dict):
        return input_code, output_code, modules
    if 'variables' not in inner:
        return input_code, output_code, modules

    for explicit in inner['variables'].split():
        if 'store' in inner:
            output_code += (
                f"if (associated(self%{explicit})) then\n"
                f" write (stateFile) .true.\n"
                f" call {inner['store']}(self%{explicit},stateFile,"
                f"gslStateFile,stateOperationID)\n"
                f"else\n"
                f" write (stateFile) .false.\n"
                f"end if\n"
            )
        if 'restore' in inner:
            input_code += (
                f"read (stateFile) wasAssociated\n"
                f"if (wasAssociated) then\n"
                f" call {inner['restore']}(self%{explicit},stateFile,"
                f"gslStateFile,stateOperationID)\n"
                f"else\n"
                f" nullify(self%{explicit})\n"
                f"end if\n"
            )
        if 'module' in inner:
            modules[inner['module']] = True

    return input_code, output_code, modules


# ---------------------------------------------------------------------------
# Allocatable store/restore (shared between functionClass and intrinsic paths)
# ---------------------------------------------------------------------------

def bounds_expression(variable_name, rank):
    """Return the expression written to the state file to record the bounds
    of the allocatable array `self%<variable_name>` of the given rank: its
    lower bounds followed by its upper bounds, as `c_size_t` integers."""
    return (f"lbound(self%{variable_name},kind=c_size_t),"
            f"ubound(self%{variable_name},kind=c_size_t)")


def bounds_specification(rank):
    """Return the bounds specification (`l1:u1,l2:u2,...`) with which to
    reallocate an array of the given rank from the `storedBounds` array read
    back from a record written by `bounds_expression`."""
    return ",".join(f"storedBounds({i}):storedBounds({rank + i})"
                    for i in range(1, rank + 1))


def generate_allocatable_state_store_code(state_store, state_stores,
                                          variable_name, rank,
                                          write_expr, read_expr):
    """Mutate `state_store` (outputCode / inputCode) and `state_stores`
    (allocatables / dimensionals / stateFileUsed / labelUsed flags) with
    the store/restore code for one allocatable member variable.
    """
    state_store.setdefault('outputCode', '')
    state_store.setdefault('inputCode',  '')
    state_stores['allocatablesFound'] = True
    state_stores['dimensionalsFound'] = True
    state_stores['stateFileUsed']     = True
    state_stores['labelUsed']         = True

    state_store['outputCode'] += (
        f" if (allocated(self%{variable_name})) then\n"
        "  if (displayVerbosity() >= verbosityLevelWorking) then\n"
        # <workaround type="gfortran" PR="94446" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94446"/>
        # Using sizeof() on a treeNode object causes a bogus "type mismatch" error, so the size is
        # reported as zero here and below.
        "   write (label,'(i16)') 0\n"
        f"   call displayMessage('storing \"{variable_name}\" with size '"
        "//trim(adjustl(label))//' bytes')\n"
        "  end if\n"
        "  write (stateFile) .true.\n"
        # Store the bounds, not just the shape: arrays declared with a
        # non-unit lower bound (e.g. `dimension(0:n)`) must be restored with
        # that same lower bound, or every index into them is shifted (#1786).
        f"  write (stateFile) {bounds_expression(variable_name, rank)}\n"
        f"  write (stateFile) {write_expr}\n"
        " else\n"
        "  write (stateFile) .false.\n"
        " end if\n"
    )
    state_store['inputCode'] += (
        " read (stateFile) wasAllocated\n"
        f" if (allocated(self%{variable_name})) deallocate(self%{variable_name})\n"
        " if (wasAllocated) then\n"
        f"  call displayMessage('restoring \"{variable_name}\"',"
        "verbosity=verbosityLevelWorking)\n"
        f"  allocate(storedBounds({2 * rank}))\n"
        "  read (stateFile) storedBounds\n"
        f"  allocate(self%{variable_name}({bounds_specification(rank)}))\n"
        "  deallocate(storedBounds)\n"
        f"  read (stateFile) {read_expr}\n"
        " end if\n"
    )


# ---------------------------------------------------------------------------
# Main declaration-driven state-store/restore emitter
# ---------------------------------------------------------------------------

def state_store_variables(state_stores, state_store, class_record,
                          declarations, explicit_names_found,
                          state_storables):
    """Walk `declarations` and emit store/restore code for every member.

    `state_stores`
    and `state_store` are the shared accumulators populated by
    FunctionClass's main driver.  Keys touched:
      - `state_stores`: `allocatablesFound`, `dimensionalsFound`,
        `stateFileUsed`, `gslStateFileUsed`, `labelUsed`, `rankMaximum`,
        `stateRestoreModules`.
      - `state_store`: `outputCode`, `inputCode`, `staticVariables`,
        `excludes`, `hasCustomStateStore`, `hasCustomStateRestore`, and
        `unstored` (a list of `(owner, member, type, source)` for each owned
        member of a type that is not state-storable, where `owner` and
        `source` are taken from `state_store['owner']` and
        `state_store['ownerSource']`).
    `explicit_names_found` is a list appended to when a member is named
    in `class.stateStorable.functionClass.variables`.
    """
    # Fortran is case-insensitive, so match type names case-insensitively
    # (a member may be declared as, e.g., `type(table1dGeneric)` while its type
    # is registered as `table1DGeneric`).
    fc_names       = {n.lower() for n in _function_class_name_set(state_storables)}
    fc_instances   = {n.lower() for n in _function_class_instances(state_storables)}
    ss_by_type     = {t.lower(): e for t, e in _state_storables_by_type(state_storables).items() if t}

    excludes = [e.lower() for e in state_store.get('excludes') or []]

    state_store.setdefault('outputCode', '')
    state_store.setdefault('inputCode',  '')
    state_store.setdefault('staticVariables', [])
    state_stores.setdefault('rankMaximum', 0)
    state_stores.setdefault('stateRestoreModules', {})

    for declaration in as_array(declarations):
        if not isinstance(declaration, dict):
            continue
        intrinsic  = declaration.get('intrinsic') or ''
        attributes = declaration.get('attributes') or []
        variables  = declaration.get('variables') or []
        type_raw   = declaration.get('type') or ''
        type_stripped = re.sub(r'\s', '', type_raw)
        has_pointer    = any(a == 'pointer' for a in attributes)
        is_allocatable = any(a == 'allocatable' for a in attributes)

        # ---- Skip type-bound procedures, finalizers, and generic operator
        # bindings.  None declare data members, and `generic :: assignment(=)
        # => …` parses as a "variable" named `assignment(=)=>…` that would
        # otherwise be emitted as `self%assignment(=)…` — invalid Fortran.
        if intrinsic in ('procedure', 'final', 'generic'):
            _maybe_flag_custom_hooks(state_store, declaration)
            continue

        # ---- class(...) / type(...) members ----
        if intrinsic in ('class', 'type'):
            # Pointer to a functionClass.
            if (intrinsic == 'class' and has_pointer
                    and type_stripped.lower() in fc_names):
                for v in variables:
                    name = strip_variable_name(v)
                    if name.lower() in excludes:
                        continue
                    state_stores['labelUsed']        = True
                    state_stores['stateFileUsed']    = True
                    state_stores['gslStateFileUsed'] = True
                    state_store['outputCode'] += (
                        " if (displayVerbosity() >= verbosityLevelWorking) then\n"
                        f"  select type (c__ => self%{name})\n"
                        f"  class is ({type_raw})\n"
                        # gfortran PR 94446 workaround.
                        "   write (label,'(i16)') 0\n"
                        "  end select\n"
                        f"  call displayMessage('storing \"{name}\" with size '"
                        "//trim(adjustl(label))//' bytes')\n"
                        " end if\n"
                        f" call self%{name}%stateStore  (stateFile,gslStateFile,stateOperationID)\n"
                    )
                    state_store['inputCode'] += (
                        f" call displayMessage('restoring \"{name}\"',"
                        "verbosity=verbosityLevelWorking)\n"
                        f" call self%{name}%stateRestore(stateFile,gslStateFile,stateOperationID)\n"
                    )
                _maybe_flag_custom_hooks(state_store, declaration)
                continue

            # Enumeration.
            if (intrinsic == 'type'
                    and re.match(r'^enumeration[a-z0-9_]+type', type_raw,
                                 re.IGNORECASE)):
                if is_allocatable:
                    rank = declaration_rank(declaration)
                    for variable_name in variables:
                        if variable_name.lower() in excludes:
                            continue
                        generate_allocatable_state_store_code(
                            state_store, state_stores, variable_name, rank,
                            f"self%{variable_name}%ID",
                            f"self%{variable_name}%ID",
                        )
                else:
                    state_store['outputCode'] += (
                        " if (displayVerbosity() >= verbosityLevelWorking) then\n"
                    )
                    names = declaration.get('variableNames') or []
                    for n in names:
                        # gfortran PR 94446 workaround.
                        state_store['outputCode'] += (
                            "   write (label,'(i16)') 0\n"
                            f"  call displayMessage('storing \"{n}\" with size '"
                            "//trim(adjustl(label))//' bytes')\n"
                        )
                    state_store['outputCode'] += (
                        " end if\n"
                        "  write (stateFile) "
                        + ",".join(f"self%{n}%ID" for n in names)
                        + "\n"
                    )
                    state_store['inputCode'] += (
                        "  read  (stateFile) "
                        + ",".join(f"self%{n}%ID" for n in names)
                        + "\n"
                    )
                _maybe_flag_custom_hooks(state_store, declaration)
                continue

            # Explicitly stateStorable or functionClass instance, or a type
            # that stores itself (via `stateStore(stateFile)` and
            # `stateRestore(stateFile)` procedures).
            is_self_storing = (intrinsic == 'type'
                               and is_self_storing_type(type_stripped))
            if (type_stripped.lower() in ss_by_type
                    or type_stripped.lower() in fc_instances
                    or is_self_storing):
                explicits = _explicit_function_class_variables(class_record)
                is_function_class = type_stripped.lower() in fc_instances
                for v in variables:
                    name = strip_variable_name(v)
                    if name.lower() in excludes:
                        continue
                    is_explicit = any(
                        e.lower() == name.lower() for e in explicits)
                    if has_pointer and not is_explicit:
                        continue
                    if is_explicit and explicit_names_found is not None:
                        explicit_names_found.append(name.lower())

                    rank = variable_rank(declaration, v)
                    if rank > state_stores['rankMaximum']:
                        state_stores['rankMaximum'] = rank

                    if is_allocatable:
                        state_stores['allocatablesFound'] = True
                        if rank > 0:
                            state_stores['dimensionalsFound'] = True
                        state_store['outputCode'] += (
                            f" if (allocated(self%{name})) then\n"
                            "  write (stateFile) .true.\n"
                        )
                        if rank > 0:
                            state_store['outputCode'] += (
                                f"  write (stateFile) {bounds_expression(name, rank)}\n"
                            )
                        state_store['inputCode'] += (
                            " read (stateFile) wasAllocated\n"
                            f" if (allocated(self%{name})) deallocate(self%{name})\n"
                            " if (wasAllocated) then\n"
                        )
                        if rank > 0:
                            state_store['inputCode'] += (
                                f"  allocate(storedBounds({2 * rank}))\n"
                                "  read (stateFile) storedBounds\n"
                            )
                        if intrinsic == 'class':
                            storable   = ss_by_type.get(type_stripped.lower(), {}) or {}
                            function_name = (
                                f"{type_stripped}ClassRestore"
                                + (f"{rank}D" if rank > 0 else '')
                            )
                            module = storable.get('module') or ''
                            state_stores['stateRestoreModules'][
                                f"{module},only:{function_name}"
                            ] = True
                            state_store['inputCode'] += (
                                f"  call {function_name}(self%{name},stateFile"
                                + (',storedBounds' if rank > 0 else '')
                                + ")\n"
                            )
                        else:
                            if rank > 0:
                                state_store['inputCode'] += (
                                    f"  allocate(self%{name}({bounds_specification(rank)}))\n"
                                )
                            else:
                                state_store['inputCode'] += (
                                    f"  allocate(self%{name})\n"
                                )
                        if rank > 0:
                            state_store['inputCode'] += "  deallocate(storedBounds)\n"

                    for i in range(1, rank + 1):
                        state_store['outputCode'] += (
                            (' ' * i)
                            + f"do i{i}=lbound(self%{name},dim={i}),"
                            + f"ubound(self%{name},dim={i})\n"
                        )
                        state_store['inputCode'] += (
                            (' ' * i)
                            + f"do i{i}=lbound(self%{name},dim={i}),"
                            + f"ubound(self%{name},dim={i})\n"
                        )
                    array_element = (
                        '(' + ','.join(f"i{i}" for i in range(1, rank + 1)) + ')'
                        if rank > 0 else ''
                    )
                    state_stores['labelUsed'] = True
                    state_store['outputCode'] += (
                        " if (displayVerbosity() >= verbosityLevelWorking) then\n"
                    )
                    if intrinsic == 'class':
                        state_store['outputCode'] += (
                            f"  select type (c__ => self%{name}{array_element})\n"
                            f"  class is ({type_raw})\n"
                            "   write (label,'(i16)') 0\n"
                            "  end select\n"
                        )
                    else:
                        state_store['outputCode'] += (
                            "   write (label,'(i16)') 0\n"
                        )
                    state_store['outputCode'] += (
                        f"  call displayMessage('storing \"{name}{array_element}\""
                        " with size '//trim(adjustl(label))//' bytes')\n"
                        " end if\n"
                    )
                    state_store['inputCode'] += (
                        f" call displayMessage('restoring \"{name}{array_element}\"',"
                        "verbosity=verbosityLevelWorking)\n"
                    )
                    if is_self_storing:
                        store_arguments   = "stateFile"
                        restore_arguments = "stateFile"
                    else:
                        store_arguments   = (
                            "stateFile,gslStateFile"
                            + (',stateOperationID' if is_function_class
                               else ',storeIdentifier='
                               + ('.true.' if intrinsic == 'class'
                                  else '.false.'))
                        )
                        restore_arguments = (
                            "stateFile,gslStateFile"
                            + (',stateOperationID' if is_function_class
                               else '')
                        )
                    state_store['inputCode'] += (
                        (' ' * rank)
                        + f" call self%{name}{array_element}"
                        + f"%stateRestore({restore_arguments})\n"
                    )
                    state_store['outputCode'] += (
                        (' ' * rank)
                        + f" call self%{name}{array_element}"
                        + f"%stateStore  ({store_arguments})\n"
                    )
                    for i in range(1, rank + 1):
                        state_store['outputCode'] += (
                            (' ' * (rank + 1 - i)) + "end do\n"
                        )
                        state_store['inputCode'] += (
                            (' ' * (rank + 1 - i)) + "end do\n"
                        )
                    if is_allocatable:
                        state_store['inputCode'] += " end if\n"
                        state_store['outputCode'] += (
                            " else\n  write (stateFile) .false.\n end if\n"
                        )
                    state_stores['stateFileUsed']    = True
                    if not is_self_storing:
                        state_stores['gslStateFileUsed'] = True
                _maybe_flag_custom_hooks(state_store, declaration)
                continue

            # Unrecognized derived type — deliberately emit nothing (the
            # `class`/`type` arm has no fall-through into the intrinsic
            # branch).  Without this explicit continue our walk dropped into
            # the `is_allocatable` intrinsic path below and emitted, e.g.,
            # `write (stateFile) self%postprocessors` for a
            # `type(stellarPopulationSpectraPostprocessorList), allocatable,
            # dimension(:)` member whose element type carries pointer
            # components — invalid Fortran I/O.
            #
            # Record each such member that the object owns (i.e. is not a
            # pointer — pointers reference objects owned elsewhere), so that
            # the caller can warn that it is not stored (issue #1786). Handles
            # (locks, open files, etc.) can never be stored, so are skipped
            # silently.
            if not has_pointer and not is_unstored_type(type_stripped):
                names_original = declaration.get('variableNames') or []
                for i, v in enumerate(variables):
                    name = strip_variable_name(v)
                    if name.lower() in excludes:
                        continue
                    if (len(names_original) == len(variables)
                            and names_original[i].lower() == name.lower()):
                        name = names_original[i]
                    state_store.setdefault('unstored', []).append(
                        (state_store.get('owner'), name, type_stripped,
                         state_store.get('ownerSource')))
            _maybe_flag_custom_hooks(state_store, declaration)
            continue

        # ---- Intrinsic types ----
        if has_pointer:
            # Pointers to intrinsics are not handled.
            _maybe_flag_custom_hooks(state_store, declaration)
            continue
        if re.match(r'^\s*omp_lock_kind\s*', type_raw):
            _maybe_flag_custom_hooks(state_store, declaration)
            continue
        if is_allocatable:
            rank = declaration_rank(declaration)
            for variable_name in variables:
                if variable_name.lower() in excludes:
                    continue
                generate_allocatable_state_store_code(
                    state_store, state_stores, variable_name, rank,
                    f"self%{variable_name}", f"self%{variable_name}",
                )
        else:
            # Static-size intrinsic.  Honor any `restoreTo` overrides from
            # the class's `<stateStorable>` block.
            restore_tos = []
            if class_record and isinstance(class_record.get('stateStorable'), dict):
                restore_tos = list(as_array(
                    class_record['stateStorable'].get('restoreTo')))
            for v in variables:
                name = strip_variable_name(v)
                if name.lower() in excludes:
                    continue
                store = True
                for rt in restore_tos:
                    if not isinstance(rt, dict):
                        continue
                    target_vars = [
                        x.strip()
                        for x in re.split(r'\s*,\s*', rt.get('variables', ''))
                        if x.strip()
                    ]
                    if any(t.lower() == name.lower() for t in target_vars):
                        store = False
                        state_store['inputCode'] += (
                            f" self%{name}={rt.get('state', '')}\n"
                        )
                if store:
                    state_store['staticVariables'].append(name)

        _maybe_flag_custom_hooks(state_store, declaration)


def _explicit_function_class_variables(class_record):
    """Return the comma-split variables list from
    `class.stateStorable.functionClass.variables`, or [].
    """
    if not class_record:
        return []
    ss = class_record.get('stateStorable') or {}
    fc = ss.get('functionClass') if isinstance(ss, dict) else None
    if not isinstance(fc, dict) or 'variables' not in fc:
        return []
    return [v.strip() for v in re.split(r'\s*,\s*', fc['variables']) if v.strip()]


def _maybe_flag_custom_hooks(state_store, declaration):
    """Set `hasCustomStateStore` / `hasCustomStateRestore` on `state_store`
    when the declaration is a `procedure :: stateStore => …` or
    `stateRestore => …` binding.
    """
    if (declaration.get('intrinsic') != 'procedure'
            or not declaration.get('variables')):
        return
    first = declaration['variables'][0]
    if re.match(r'^stateStore=>',   first):
        state_store['hasCustomStateStore']   = True
    if re.match(r'^stateRestore=>', first):
        state_store['hasCustomStateRestore'] = True
