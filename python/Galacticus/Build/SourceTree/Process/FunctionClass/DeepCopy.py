# Deep-copy code generation helpers for the functionClass pipeline.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/FunctionClass/DeepCopy.pm
# — all three exported functions (`deep_copy_copied_self_block`,
# `generate_assignment_allocatable_code`, `deep_copy_declarations`).  Perl
# aliases `$stateStorables` / `$deepCopyActions` from the parent package
# at compile time; we take both as explicit parameters instead — cleaner
# Python and decouples the helper from a global.

import os
import re
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from List.ExtraUtils                                         import as_array
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    strip_variable_name, declaration_rank,
)
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


# ---------------------------------------------------------------------------
# stateStorables shape helpers — bridge our xml_to_dict output to the
# Perl `$stateStorables->{functionClasses}` keyAttr-keyed shape.
# ---------------------------------------------------------------------------

def _function_class_name_set(state_storables):
    """Return the set of functionClass name keys."""
    fc = (state_storables or {}).get('functionClasses') or {}
    if not isinstance(fc, dict):
        return set()
    entries = fc.get('functionClass')
    if entries is None:
        return set(fc.keys())
    if isinstance(entries, dict):
        entries = [entries]
    return {e.get('name') for e in entries
            if isinstance(e, dict) and 'name' in e}


def _function_class_instances(state_storables):
    """Return the list of functionClassInstances names."""
    raw = (state_storables or {}).get('functionClassInstances') or []
    if isinstance(raw, str):
        return [raw] if raw else []
    if isinstance(raw, dict):
        name = raw.get('content') or raw.get('name')
        return [name] if name else []
    out = []
    for item in raw:
        if isinstance(item, str):
            out.append(item)
        elif isinstance(item, dict):
            name = item.get('content') or item.get('name')
            if name:
                out.append(name)
    return out


# ---------------------------------------------------------------------------
# Small public helpers
# ---------------------------------------------------------------------------

def deep_copy_copied_self_block(deep_copy, name, declaration,
                                non_abstract_class, indent=''):
    """Append the `copiedSelf` select-type deep-copy block to
    `deep_copy['assignments']`.

    Mirrors deepCopyCopiedSelfBlock() at DeepCopy.pm:20-37.
    """
    node = non_abstract_class.get('node') or {}
    loc_expr = location(node, node.get('line', 0))
    intrinsic = declaration.get('intrinsic') or ''
    type_name = declaration.get('type') or ''
    deep_copy.setdefault('assignments', '')
    deep_copy['assignments'] += (
        f"{indent}if (associated(self%{name}%copiedSelf)) then\n"
        f"{indent} select type(s => self%{name}%copiedSelf)\n"
        f"{indent} {intrinsic} is ({type_name})\n"
        f"{indent}  destination%{name} => s\n"
        f"{indent} class default\n"
        f"{indent}  call Error_Report('copiedSelf has incorrect type'//{loc_expr})\n"
        f"{indent} end select\n"
        f"{indent} call self%{name}%copiedSelf%referenceCountIncrement()\n"
        f"{indent}else\n"
        f"{indent} allocate(destination%{name},mold=self%{name})\n"
    )


def generate_assignment_allocatable_code(assignment, declaration, name,
                                         allocated, rank_maximum_ref):
    """Generate code that assigns an allocatable/association-forced member.

    `assignment['code']` is mutated in place.  `rank_maximum_ref` is a
    mutable single-element list used to signal the caller the maximum rank
    seen across calls (Perl passes a scalar by reference; Python idiomatic
    equivalent is a 1-element list).

    Mirrors generateAssignmentAllocatableCode() at DeepCopy.pm:39-75.
    """
    assignment.setdefault('code', '')
    assignment['code'] += (
        f"    if ({allocated}(self%{name})) deallocate(self%{name})\n"
        f"    if ({allocated}(from%{name})) then\n"
    )
    rank = declaration_rank(declaration)
    bounds = (
        '(' + ','.join(
            f"lbound(from%{name},dim={i}):ubound(from%{name},dim={i})"
            for i in range(1, rank + 1)
        ) + ')'
        if rank > 0 else ''
    )
    assignment['code'] += (
        f"      allocate(self%{name}{bounds},mold=from%{name})\n"
    )
    intrinsic = declaration.get('intrinsic')
    if rank > 0 and intrinsic == 'type':
        # gfortran PR 46897 / 57696 workaround — use explicit element-wise
        # assignment rather than a single slice assignment.
        for i in range(1, rank + 1):
            assignment['code'] += (
                f"      do i{i}__=lbound(from%{name},dim={i}),"
                f"ubound(from%{name},dim={i})\n"
            )
        indices = ','.join(f"i{i}__" for i in range(1, rank + 1))
        assignment['code'] += (
            f"      self%{name}({indices})=from%{name}({indices})\n"
        )
        for _ in range(rank):
            assignment['code'] += "      end do\n"
        if rank > rank_maximum_ref[0]:
            rank_maximum_ref[0] = rank
    else:
        assignment['code'] += f"      self%{name}=from%{name}\n"
    assignment['code'] += "    end if\n"


# ---------------------------------------------------------------------------
# Main declaration-driven deep-copy emitter
# ---------------------------------------------------------------------------

def deep_copy_declarations(class_record, non_abstract_class, node,
                           declarations, ignore, line_number,
                           deep_copy, found_deep_copy_names,
                           state_storables, deep_copy_actions):
    """Process one node's declarations for deep copy.

    Mirrors deepCopyDeclarations() at DeepCopy.pm:77-257.  `deep_copy` is
    the shared accumulator dict populated by FunctionClass's main driver
    (see the Perl parent package).  Keys we touch: `assignments`,
    `resetCode`, `finalizeCode`, `needReferenceCount`, `rankMaximum`,
    `modules`, `resetModules`, `finalizeModules`.
    `found_deep_copy_names` is appended to when a member named in
    `class.deepCopy.functionClass.variables` is matched.
    """
    fc_names      = _function_class_name_set(state_storables)
    fc_instances  = set(_function_class_instances(state_storables))
    actions_set   = _deep_copy_actions_type_set(deep_copy_actions)
    ignore_lc     = {i.lower() for i in (ignore or [])}

    deep_copy.setdefault('assignments',  '')
    deep_copy.setdefault('resetCode',    '')
    deep_copy.setdefault('finalizeCode', '')
    deep_copy.setdefault('modules',         {})
    deep_copy.setdefault('resetModules',    {})
    deep_copy.setdefault('finalizeModules', {})
    deep_copy.setdefault('rankMaximum', 0)

    for declaration in as_array(declarations):
        if not isinstance(declaration, dict):
            continue
        intrinsic  = declaration.get('intrinsic') or ''
        type_raw   = declaration.get('type') or ''
        # Perl strips only leading/trailing whitespace when intrinsic is
        # class or type; we replicate for the same cases.
        type_stripped = (type_raw.strip()
                         if intrinsic in ('class', 'type')
                         else type_raw)
        attributes = declaration.get('attributes') or []
        variables  = declaration.get('variables') or []
        has_pointer     = any(a == 'pointer'     for a in attributes)
        is_allocatable  = any(a == 'allocatable' for a in attributes)

        # ---- 1. Pointer-to-functionClass objects ----
        if (intrinsic == 'class'
                and (type_stripped in fc_names
                     or type_stripped in fc_instances)
                and has_pointer):
            for obj in variables:
                name = strip_variable_name(obj)
                if name.lower() in ignore_lc:
                    continue
                deep_copy['resetCode'] += (
                    f"if (associated(self%{name})) "
                    f"call self%{name}%deepCopyReset   ()\n"
                )
                deep_copy['finalizeCode'] += (
                    f"if (associated(self%{name})) "
                    f"call self%{name}%deepCopyFinalize()\n"
                )
                deep_copy['assignments'] += (
                    f"nullify(destination%{name})\n"
                    f"if (associated(self%{name})) then\n"
                )
                deep_copy['needReferenceCount'] = 1
                deep_copy['assignments'] += (
                    f" referenceCount__=self%{name}%referenceCountDecrement()\n"
                )
                deep_copy_copied_self_block(
                    deep_copy, name, declaration, non_abstract_class, " ")
                deep_copy['assignments'] += (
                    f"  call self%{name}%deepCopy(destination%{name})\n"
                    f"  self%{name}%copiedSelf => destination%{name}\n"
                    f"  call destination%{name}%autoHook()\n"
                    f" end if\n"
                    f"end if\n"
                )

        # ---- 2. Objects with explicit deepCopyActions bindings ----
        if (intrinsic in ('class', 'type')
                and type_stripped in actions_set):
            rank = declaration_rank(declaration)
            if rank > deep_copy['rankMaximum']:
                deep_copy['rankMaximum'] = rank
            for variable_name in declaration.get('variableNames') or []:
                if is_allocatable:
                    deep_copy['assignments'] += (
                        f"if (allocated(self%{variable_name})) then\n"
                    )
                if has_pointer:
                    deep_copy['assignments'] += (
                        f"if (associated(self%{variable_name})) then\n"
                    )
                for i in range(1, rank + 1):
                    deep_copy['assignments'] += (
                        (' ' * i)
                        + f"do i{i}=lbound(self%{variable_name},dim={i}),"
                        + f"ubound(self%{variable_name},dim={i})\n"
                    )
                array_element = (
                    '(' + ','.join(f"i{i}" for i in range(1, rank + 1)) + ')'
                    if rank > 0 else ''
                )
                deep_copy['assignments'] += (
                    (' ' * rank)
                    + f"call destination%{variable_name}{array_element}"
                    + "%deepCopyActions()\n"
                )
                for i in range(1, rank + 1):
                    deep_copy['assignments'] += (
                        (' ' * (rank + 1 - i)) + "end do\n"
                    )
                if is_allocatable or has_pointer:
                    deep_copy['assignments'] += "end if\n"

        # ---- 3. HDF5 objects ----
        if intrinsic == 'type' and re.match(
                r'^\s*hdf5object\s*$', type_raw, re.IGNORECASE):
            deep_copy['modules']['HDF5_Access'] = True
            deep_copy['assignments'] += "!$ call hdf5Access%set  ()\n"
            for v in variables:
                deep_copy['assignments'] += (
                    f"call self%{v}%deepCopy(destination%{v})\n"
                )
            deep_copy['assignments'] += "!$ call hdf5Access%unset()\n"

        # ---- 4. Non-(class,pointer) functionClass objects from explicit list ----
        explicit_fc = _variables_attr(
            (class_record or {}).get('deepCopy', {}).get('functionClass'))
        if explicit_fc:
            for obj in variables:
                name = strip_variable_name(obj)
                if any(v.lower() == name.lower() for v in explicit_fc):
                    if found_deep_copy_names is not None:
                        found_deep_copy_names.append(name)
                    if has_pointer:
                        deep_copy['assignments'] += (
                            f"nullify(destination%{name})\n"
                            f"if (associated(self%{name})) then\n"
                        )
                        deep_copy['resetCode'] += (
                            f"if (associated(self%{name})) then\n"
                        )
                        deep_copy['finalizeCode'] += (
                            f"if (associated(self%{name})) then\n"
                        )
                        deep_copy['needReferenceCount'] = 1
                        deep_copy['assignments'] += (
                            f"referenceCount__=self%{name}"
                            f"%referenceCountDecrement()\n"
                        )
                        deep_copy_copied_self_block(
                            deep_copy, name, declaration, non_abstract_class)
                    deep_copy['resetCode'] += (
                        f"call self%{name}%deepCopyReset   ()\n"
                    )
                    deep_copy['finalizeCode'] += (
                        f"call self%{name}%deepCopyFinalize()\n"
                    )
                    deep_copy['assignments'] += (
                        f"call self%{name}%deepCopy(destination%{name})\n"
                        f"self%{name}%copiedSelf => destination%{name}\n"
                        f"call destination%{name}%autoHook()\n"
                    )
                    if has_pointer:
                        deep_copy['assignments']  += "end if\n"
                        deep_copy['resetCode']    += "end if\n"
                        deep_copy['finalizeCode'] += "end if\n"

        # ---- 5. Increments (atomic or not) ----
        increment_block = (class_record or {}).get('deepCopy', {}).get('increment')
        if increment_block:
            inc_vars = _variables_attr(increment_block)
            atomic   = increment_block.get('atomic') == 'yes' \
                       if isinstance(increment_block, dict) else False
            increments = [{
                'variable': v,
                'host':     re.sub(r'^([^%]+)%.+', r'\1', v),
            } for v in inc_vars]
            for obj in variables:
                name = strip_variable_name(obj)
                for inc in increments:
                    if inc['host'].lower() == name.lower():
                        if atomic:
                            deep_copy['assignments'] += "!$omp atomic\n"
                        deep_copy['assignments'] += (
                            f"destination%{inc['variable']}="
                            f"destination%{inc['variable']}+1\n"
                        )

        # ---- 6. setTo blocks ----
        set_to_block = (class_record or {}).get('deepCopy', {}).get('setTo')
        if set_to_block:
            set_vars = _variables_attr(set_to_block)
            state    = (set_to_block.get('value') if isinstance(set_to_block, dict)
                        else '')
            set_tos  = [{
                'variable': v,
                'host':     re.sub(r'^([^%]+)%.+', r'\1', v),
            } for v in set_vars]
            for obj in variables:
                name = strip_variable_name(obj)
                for st in set_tos:
                    if st['host'].lower() == name.lower():
                        deep_copy['assignments'] += (
                            f"destination%{st['variable']}={state}\n"
                        )

        # ---- 7. Explicit deepCopy blocks ----
        dc_block = (class_record or {}).get('deepCopy', {}).get('deepCopy')
        if dc_block and isinstance(dc_block, dict):
            dc_vars = _variables_attr(dc_block)
            for obj in declaration.get('variableNames') or []:
                for target in dc_vars:
                    if obj.lower() != target.lower():
                        continue
                    deep_copy['assignments'] += (
                        f"nullify(destination%{obj})\n"
                        f"allocate(destination%{obj},mold=self%{obj})\n"
                    )
                    copy_fn = dc_block.get('copy')
                    module  = dc_block.get('module')
                    reset   = dc_block.get('reset')
                    finalize = dc_block.get('finalize')
                    if copy_fn is not None:
                        if module:
                            deep_copy['modules'][module] = True
                        deep_copy['assignments'] += (
                            f"if (associated(self%{obj})) call {copy_fn}"
                            f"(self%{obj},destination%{obj})\n"
                        )
                    else:
                        deep_copy['assignments'] += (
                            f"if (associated(self%{obj})) "
                            f"call self%{obj}%deepCopy(destination%{obj})\n"
                        )
                    if reset is not None:
                        if module:
                            deep_copy['resetModules'][module] = True
                        deep_copy['resetCode'] += (
                            f"if (associated(self%{obj})) "
                            f"call {reset}(self%{obj})\n"
                        )
                    if finalize is not None:
                        if module:
                            deep_copy['finalizeModules'][module] = True
                        deep_copy['finalizeCode'] += (
                            f"if (associated(self%{obj})) "
                            f"call {finalize}(self%{obj})\n"
                        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _variables_attr(block):
    """Extract the `variables` attribute of a directive block, split on
    commas.  Returns [] if the block is missing or malformed.
    """
    if not isinstance(block, dict):
        return []
    vars_text = block.get('variables') or ''
    return [v.strip() for v in re.split(r'\s*,\s*', vars_text) if v.strip()]


def _deep_copy_actions_type_set(deep_copy_actions):
    """Return the set of types named in `deepCopyActions/deepCopyActions` —
    those are classes with an explicit deepCopyActions hook.
    """
    raw = (deep_copy_actions or {}).get('deepCopyActions') or []
    if isinstance(raw, dict):
        raw = [raw]
    return {
        e.get('type') for e in raw
        if isinstance(e, dict) and 'type' in e
    }
