# Processes `stateStorable` directives: for each directive naming a base
# class, synthesizes `XStateStore(self,…)` and `XStateRestore(self,…)`
# subroutines that walk the class hierarchy, emitting per-variable
# write/read of every non-excluded, state-storable member; also synthesizes
# the `XClassRestore` / `XClassRestore1D` dispatchers that rebuild the
# correct dynamic type on read.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/StateStorable.pm

import os
import re
import sys
import xml.etree.ElementTree as ET


from List.ExtraUtils                     import as_array
from XML.Utils                           import xml_to_dict
from Galacticus.Build.SourceTree         import (
    walk_tree, parse_code, children, set_visibility, insert_post_contains,
)
from Galacticus.Build.SourceTree.Process import register_process, process_tree


_TYPE_OPENER_RE = re.compile(
    r'^\s*type\s*'
    r'((?:,\s*(?:abstract|public|private|extends\s*\([a-zA-Z0-9_]+\))\s*)*)'
    r'(?:::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)
_EXTENDS_ATTR_RE = re.compile(r'extends\(([a-zA-Z0-9_]+)\)', re.IGNORECASE)


# ---------------------------------------------------------------------------
# Shared parsers / lookups
# ---------------------------------------------------------------------------

def _parse_type_opener(opener):
    """Return `(name, extends_or_None, abstract_bool)` for a `type … :: name`
    opener line, or None for openers carrying unexpanded `{Type¦…}` generic
    placeholders (StateStorable doesn't support generics, so those types are
    silently skipped — DeepCopyActions makes the same carve-out at
    DeepCopyActions.pm:62).  Raises only on genuinely unparseable openers.
    """
    if '{' in opener:
        return None
    m = _TYPE_OPENER_RE.match(opener)
    if not m:
        raise RuntimeError(
            "process_state_storable: unable to parse type definition opener")
    attrs_text = m.group(1).strip().strip(',').strip()
    type_name  = m.group(2)
    extends    = None
    abstract   = False
    for attr in re.split(r'\s*,\s*', attrs_text) if attrs_text else []:
        if attr == 'abstract':
            abstract = True
        em = _EXTENDS_ATTR_RE.match(attr)
        if em:
            extends = em.group(1)
    return type_name, extends, abstract


def _load_state_storables_xml():
    """Read `$BUILDPATH/stateStorables.xml` and return its parsed dict.
    Missing file or env var → empty dict (callers treat that as "no registered
    stateStorables", which is the fail-safe path).
    """
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        return {}
    path = os.path.join(build_path, 'stateStorables.xml')
    if not os.path.exists(path):
        return {}
    return xml_to_dict(ET.parse(path).getroot())


def _storable_types(state_storables):
    """Return the set of registered stateStorable type names.  Matches Perl's
    `grep {$_->{'type'} eq $type} @{$stateStorables->{'stateStorables'}}`.
    """
    entries = (state_storables or {}).get('stateStorables')
    return {e.get('type') for e in as_array(entries) if isinstance(e, dict)}


def _descends_from(class_name, target, classes):
    """True if class_name is target or any of its transitive parents is target."""
    cursor = class_name
    while cursor is not None:
        if cursor == target:
            return True
        cursor = classes.get(cursor, {}).get('extends')
    return False


# ---------------------------------------------------------------------------
# Declaration-level helpers
# ---------------------------------------------------------------------------

def _attribute_rank(attributes):
    """Return the rank extracted from `dimension(:,…)` attributes, 0 if scalar."""
    dim_inner = ''.join(
        (m.group(1) if m else '')
        for m in (re.match(r'^dimension\s*\(([:,]+)\)', a) for a in attributes)
    )
    return dim_inner.count(':')


def _strip_init(var):
    """Strip any `=…` initializer from a variable token and return the name."""
    return re.sub(r'\s*=.*$', '', var)


# ---------------------------------------------------------------------------
# Per-declaration code generators
# ---------------------------------------------------------------------------

def _process_derived_declaration(declaration, scope_directive, exclude, storable_types):
    """Generate store/restore text for a `class(…)` or `type(…)` member.

    Mirrors StateStorable.pm:205-303.  Returns a dict of fragments plus flags
    needed to emit local-variable declarations in the outer subroutine.
    """
    fragments = {
        'output':                '',
        'input':                 '',
        'class_functions_used':  False,
        'stored_shape_required': False,
        'was_allocated_required': False,
        'label_required':        False,
        'rank_seen':             0,
    }

    type_name = re.sub(r'\s', '', declaration.get('type') or '')
    attributes = declaration.get('attributes') or []
    pointer_store = (scope_directive is not None
                     and 'pointerStore' in scope_directive)
    pointer_store_names = []
    if pointer_store:
        ps = scope_directive['pointerStore']
        if isinstance(ps, dict):
            pointer_store_names = (ps.get('variables') or '').split()

    if type_name not in storable_types:
        return fragments
    has_pointer = any(a == 'pointer' for a in attributes)
    if has_pointer and not pointer_store:
        return fragments

    rank       = _attribute_rank(attributes)
    fragments['rank_seen'] = rank
    allocatable = any(a in ('allocatable', 'pointer') for a in attributes)

    variables = declaration.get('variables') or []
    for v in variables:
        name = _strip_init(v)
        if any(x.lower() == name.lower() for x in exclude):
            continue
        if has_pointer and not (pointer_store and
                                any(n.lower() == name.lower()
                                    for n in pointer_store_names)):
            continue
        fragments['label_required'] = True

        out = ''
        inp = ''
        if allocatable:
            assoc = 'associated' if has_pointer else 'allocated'
            out += f"  if ({assoc}(self%{name})) then\n"
            out += "   write (stateFile) .true.\n"
            if rank > 0:
                out += (
                    f"   write (stateFile) shape(self%{name},kind=c_size_t)\n"
                )

        out += " if (displayVerbosity() >= verbosityLevelWorking) then\n"
        if declaration.get('intrinsic') == 'class':
            out += f"  select type (c__ => self%{name})\n"
            out += f"  class is ({declaration['type']})\n"
            out += "   write (label,'(i16)') sizeof(c__)\n"
            out += "  end select\n"
        else:
            out += f"   write (label,'(i16)') sizeof(self%{name})\n"
        out += (
            f"  call displayMessage('storing \"{name}\" with size '"
            "//trim(adjustl(label))//' bytes')\n"
        )
        out += " end if\n"

        if rank > 0:
            for i in range(1, rank + 1):
                out += (" " * i) + (
                    f"do i{i}=1,size(self%{name},dim={i})\n"
                )
        idx_list = (
            '(' + ','.join(f"i{i}" for i in range(1, rank + 1)) + ')'
            if rank > 0 else ''
        )
        store_id_flag = (
            '.true.' if declaration.get('intrinsic') == 'class' else '.false.'
        )
        out += (
            f" call self%{name}{idx_list}%stateStore  "
            f"(stateFile,gslStateFile,storeIdentifier={store_id_flag})\n"
        )
        if rank > 0:
            for i in range(rank, 0, -1):
                out += (" " * i) + "end do\n"
        if allocatable:
            out += "  else\n"
            out += "   write (stateFile) .false.\n"
            out += "  end if\n"

        if allocatable:
            fragments['was_allocated_required'] = True
            inp += " read (stateFile) wasAllocated\n"
            assoc = 'associated' if has_pointer else 'allocated'
            inp += (
                f" if ({assoc}(self%{name})) deallocate(self%{name})\n"
            )
            inp += " if (wasAllocated) then\n"
        inp += (
            f"  call displayMessage('restoring \"{name}\"',"
            "verbosity=verbosityLevelWorking)\n"
        )
        if allocatable:
            if rank > 0:
                fragments['stored_shape_required'] = True
                inp += f"  allocate(storedShape({rank}))\n"
                inp += "  read (stateFile) storedShape\n"
                inp += (
                    f"  allocate(self%{name}("
                    + ",".join(f"storedShape({i})" for i in range(1, rank + 1))
                    + "))\n"
                )
                inp += "  deallocate(storedShape)\n"
            else:
                inp += f"  allocate(self%{name})\n"
        if declaration.get('intrinsic') == 'class':
            inp += f" call {type_name}ClassRestore(self%{name},stateFile)\n"
            fragments['class_functions_used'] = True
        if rank > 0:
            for i in range(1, rank + 1):
                inp += (" " * i) + (
                    f"do i{i}=1,size(self%{name},dim={i})\n"
                )
        inp += (
            f" call self%{name}{idx_list}%stateRestore(stateFile,gslStateFile)\n"
        )
        if rank > 0:
            for i in range(rank, 0, -1):
                inp += (" " * i) + "end do\n"
        if allocatable:
            inp += " end if\n"

        fragments['output'] += out
        fragments['input']  += inp

    return fragments


def _process_allocatable_intrinsic(declaration, exclude):
    """Store/restore for an allocatable intrinsic-type member.

    Mirrors StateStorable.pm:310-347.
    """
    fragments = {
        'output': '',
        'input':  '',
        'stored_shape_required':  False,
        'was_allocated_required': False,
        'label_required':         False,
        'rank_seen':              0,
    }
    attributes = declaration.get('attributes') or []
    rank = _attribute_rank(attributes)
    # NB: deliberately do NOT propagate `rank` into `fragments['rank_seen']`.
    # The output/input bodies below use Fortran whole-array I/O
    # (`write (stateFile) self%name`) and never emit a `do iN=...` loop, so
    # the caller's `rank_maximum` (which gates the `integer(c_size_t) :: i1
    # … iN` declaration) must not be bumped on our account — otherwise the
    # compiler emits "Unused variable" warnings for `i2`, `i3`, … in the
    # generated stateStore subroutine when the only multi-rank field in
    # the type is an allocatable intrinsic.
    for name in declaration.get('variables', []):
        if any(x.lower() == name.lower() for x in exclude):
            continue
        if rank > 0:
            fragments['stored_shape_required'] = True
        fragments['was_allocated_required'] = True
        fragments['label_required']         = True

        out  = f"  if (allocated(self%{name})) then\n"
        out += "   if (displayVerbosity() >= verbosityLevelWorking) then\n"
        out += f"    write (label,'(i16)') sizeof(self%{name})\n"
        out += (
            f"    call displayMessage('storing \"{name}\" with size '"
            "//trim(adjustl(label))//' bytes')\n"
        )
        out += "   end if\n"
        out += "   write (stateFile) .true.\n"
        if rank > 0:
            out += f"   write (stateFile) shape(self%{name},kind=c_size_t)\n"
        out += f"   write (stateFile) self%{name}\n"
        out += "  else\n"
        out += "   write (stateFile) .false.\n"
        out += "  end if\n"

        inp  = " read (stateFile) wasAllocated\n"
        inp += f" if (allocated(self%{name})) deallocate(self%{name})\n"
        inp += " if (wasAllocated) then\n"
        inp += (
            f"  call displayMessage('restoring \"{name}\"',"
            "verbosity=verbosityLevelWorking)\n"
        )
        if rank > 0:
            inp += f"  allocate(storedShape({rank}))\n"
            inp += "  read (stateFile) storedShape\n"
            inp += (
                f"  allocate(self%{name}("
                + ",".join(f"storedShape({i})" for i in range(1, rank + 1))
                + "))\n"
            )
            inp += "  deallocate(storedShape)\n"
        else:
            inp += f"  allocate(self%{name})\n"
        inp += f"  read (stateFile) self%{name}\n"
        inp += " end if\n"

        fragments['output'] += out
        fragments['input']  += inp
    return fragments


# ---------------------------------------------------------------------------
# Class-restore dispatcher generator
# ---------------------------------------------------------------------------

def _emit_class_restore(parent_class, target_class, classes, class_identifiers):
    """Emit `XClassRestore` + `XClassRestore1D` dispatch subroutines.

    Mirrors StateStorable.pm:462-523.  Each subroutine reads a stored
    classIdentifier and reallocates `self` to the appropriate concrete type.
    """
    text = ''
    for rank in (0, 1):
        suffix       = f"{rank}D" if rank > 0 else ''
        stored_shape = ',storedShape' if rank > 0 else ''
        dimensions   = (
            ', dimension(' + ','.join(':' for _ in range(rank)) + ')'
            if rank > 0 else ''
        )
        code  = f"subroutine {parent_class}ClassRestore{suffix}(self,stateFile{stored_shape})\n"
        code += " !!{\n"
        code += " Restore the class of this object from file.\n"
        code += " !!}\n"
        code += " use            :: Error        , only : Error_Report\n"
        code += " use, intrinsic :: ISO_C_Binding, only : c_size_t\n"
        code += " implicit none\n"
        code += (
            f" class  ({parent_class}), intent(inout), allocatable{dimensions} :: self\n"
        )
        code += " integer                    , intent(in   )               :: stateFile\n"
        code += " integer                                                  :: classIdentifier\n"
        if rank > 0:
            code += (
                "integer(c_size_t), intent(in   ), dimension("
                + ",".join(':' for _ in range(rank))
                + ") :: storedShape\n"
            )
        code += " read (stateFile) classIdentifier\n"
        code += " select case (classIdentifier)\n"
        for child_class in sorted(class_identifiers.keys()):
            # Include only children that extend `parent_class`.
            if not _descends_from(child_class, parent_class, classes):
                continue
            code += f" case ({class_identifiers[child_class]})\n"
            code += "  if (allocated(self)) deallocate(self)\n"
            shape_args = (
                '(' + ','.join(f"storedShape({i})" for i in range(1, rank + 1)) + ')'
                if rank > 0 else ''
            )
            code += f"  allocate({child_class} :: self{shape_args})\n"
        code += " case default\n"
        code += "  call Error_Report('serialized object is of incorrect class')\n"
        code += " end select\n"
        code += " return\n"
        code += f"end subroutine {parent_class}ClassRestore{suffix}\n"
        text += code + "\n"
    return text


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------

_TYPE_BINDING = (
    "    !![\n"
    "    <methods>\n"
    "     <method method=\"stateStore\"   description=\"Store the state of this object to file.\"    />\n"
    "     <method method=\"stateRestore\" description=\"Restore the state of this object from file.\"/>\n"
    "    </methods>\n"
    "    !!]\n"
    "    procedure :: stateStore   => {class_name}StateStore\n"
    "    procedure :: stateRestore => {class_name}StateRestore\n"
)


def _insert_parsed(parent, source_text, run_process_tree=False):
    sub_tree = parse_code(source_text, name='StateStorable')
    if run_process_tree:
        process_tree(sub_tree)
    kids = children(sub_tree)
    for k in kids:
        k['parent'] = None
    insert_post_contains(parent, kids)


def process_state_storable(tree, options):
    """Mirrors Process_StateStorable() from StateStorable.pm."""
    directive_nodes = []
    classes         = {}
    module_node     = None

    # Pass 1: collect directives + type nodes in a single walk.
    for node in walk_tree(tree):
        ntype = node.get('type')
        if ntype == 'stateStorable':
            directive = node.setdefault('directive', {})
            if directive.get('processed'):
                continue
            parent = node.get('parent')
            if parent is not None and parent.get('type') == 'module':
                module_node = parent
            directive['processed'] = True
            directive_nodes.append(node)
            continue
        if ntype == 'type':
            parsed = _parse_type_opener(node.get('opener', ''))
            if parsed is None:
                # Generic-templated opener (e.g. `type :: {Type¦label}Foo`) —
                # ignore it; the expanded copies will be re-processed after
                # generics runs.
                continue
            name, extends, abstract = parsed
            classes[name] = {
                'node':     node,
                'extends':  extends,
                'abstract': abstract,
            }

    if not directive_nodes:
        return

    state_storables = _load_state_storables_xml()
    storable_types  = _storable_types(state_storables)

    # Cull classes not in any directive's hierarchy.  Matches Perl:78-92.
    targets = [n['directive']['class'] for n in directive_nodes]
    for cname in list(classes.keys()):
        matched = False
        cursor  = cname
        while cursor is not None:
            if cursor in targets:
                matched = True
                break
            cursor = classes.get(cursor, {}).get('extends')
        if not matched:
            del classes[cname]

    class_functions_required = module_node is not None

    for directive_node in directive_nodes:
        directive = directive_node['directive']
        target    = directive['class']
        if target not in classes:
            raise RuntimeError(
                f"process_state_storable: class '{target}' not found")

        function_code, class_fn_used, restore_funcs = _emit_store_restore(
            directive, classes, storable_types)

        # Type-binding on the base class.
        _insert_parsed(
            classes[target]['node'],
            _TYPE_BINDING.format(class_name=target),
        )

        if (class_functions_required or class_fn_used) and module_node is not None:
            function_code += restore_funcs
            for name in _class_restore_visibility_names(
                    directive, classes):
                set_visibility(module_node, name, 'public')

        _insert_parsed(
            directive_node['parent'], function_code, run_process_tree=True)


def _class_restore_visibility_names(directive, classes):
    """Return every `XClassRestore{,1D}` symbol that should be made public for
    the directive's class hierarchy.  Mirrors the SetVisibility loop inside
    the classFunctionsRequired branch at StateStorable.pm:519-521.
    """
    names = []
    target = directive['class']
    for parent_class in sorted(classes.keys()):
        if not _descends_from(parent_class, target, classes):
            continue
        for rank in (0, 1):
            suffix = f"{rank}D" if rank > 0 else ''
            names.append(f"{parent_class}ClassRestore{suffix}")
    return names


def _emit_store_restore(directive, classes, storable_types):
    """Render both the `XStateStore` and `XStateRestore` subroutines (joined)
    plus the class-restore dispatchers.

    Returns (function_code_text, class_functions_used, class_restore_text).
    """
    target = directive['class']
    class_name_ = target

    # Openers (grow if label/storedShape/wasAllocated/ranks are needed).
    output_opener = (
        f"subroutine {target}StateStore(self,stateFile,gslStateFile,storeIdentifier)\n"
        " !!{\n"
        " Store the state of this object to file.\n"
        " !!}\n"
        " use, intrinsic :: ISO_C_Binding     , only : c_size_t, c_ptr\n"
        " use            :: ISO_Varying_String\n"
        " use            :: Display\n"
        " implicit none\n"
        f" class    ({target}), intent(inout)              :: self\n"
        " integer                , intent(in   )              :: stateFile\n"
        " type     (c_ptr       ), intent(in   )              :: gslStateFile\n"
        " logical                , intent(in   ), optional    :: storeIdentifier\n"
    )
    input_opener = (
        f"subroutine {target}StateRestore(self,stateFile,gslStateFile)\n"
        " !!{\n"
        " Store the state of this object to file.\n"
        " !!}\n"
        " use, intrinsic :: ISO_C_Binding     , only : c_size_t, c_ptr\n"
        " use            :: ISO_Varying_String\n"
        " use            :: Display\n"
        " implicit none\n"
        f" class  ({target}), intent(inout)               :: self\n"
        " integer              , intent(in   )               :: stateFile\n"
        " type   (c_ptr       ), intent(in   )               :: gslStateFile\n"
    )
    output_closer = (
        " end select\n"
        " call displayUnindent('done',verbosity=verbosityLevelWorking)\n"
        " return\n"
        f"end subroutine {target}StateStore\n"
    )
    input_closer = (
        " end select\n"
        " call displayUnindent('done',verbosity=verbosityLevelWorking)\n"
        " return\n"
        f"end subroutine {target}StateRestore\n"
    )

    output_body = (
        f" call displayIndent('storing state for \"{target}\"',verbosity=verbosityLevelWorking)\n"
        " select type (self)\n"
    )
    input_body = (
        f" call displayIndent('restoring state for \"{target}\"',verbosity=verbosityLevelWorking)\n"
        " select type (self)\n"
    )

    stored_shape_required  = False
    was_allocated_required = False
    label_required         = False
    class_functions_used   = False
    rank_maximum           = 0

    class_identifier = -1
    class_identifiers = {}

    for cname in sorted(classes.keys()):
        if not _descends_from(cname, target, classes) or classes[cname]['abstract']:
            continue

        class_identifier += 1
        class_identifiers[cname] = class_identifier
        output_body += f" type is ({cname})\n"
        input_body  += f" type is ({cname})\n"
        output_body += (
            f"  if (present(storeIdentifier).and.storeIdentifier) "
            f"write (stateFile) {class_identifier}\n"
        )

        static_variables = []
        method_calls = []

        cursor = cname
        while cursor is not None:
            scope_directive = directive.get(cursor) if isinstance(
                directive.get(cursor), dict) else None
            exclude = []
            if (scope_directive is not None
                    and isinstance(scope_directive.get('exclude'), dict)
                    and 'variables' in scope_directive['exclude']):
                exclude = [
                    v.strip()
                    for v in re.split(
                        r'\s*,\s*', scope_directive['exclude']['variables'])
                    if v.strip()
                ]

            class_node = classes[cursor]['node']
            child = class_node.get('firstChild')
            while child is not None:
                if child.get('type') == 'declaration':
                    for declaration in child.get('declarations', []):
                        intr = declaration.get('intrinsic')
                        # Skip type-bound markers (`procedure`/`generic`/
                        # `final`).  None of these declare data members;
                        # treating them as data emits invalid Fortran like
                        # `sizeof(self%assignment(=)=>integratorAssign)` for
                        # generic-operator bindings.
                        if intr in ('procedure', 'generic', 'final'):
                            continue
                        if intr in ('class', 'type'):
                            frag = _process_derived_declaration(
                                declaration, scope_directive, exclude,
                                storable_types)
                            output_body += frag['output']
                            input_body  += frag['input']
                            class_functions_used |= frag['class_functions_used']
                            stored_shape_required |= frag['stored_shape_required']
                            was_allocated_required |= frag['was_allocated_required']
                            label_required |= frag['label_required']
                            rank_maximum = max(rank_maximum, frag['rank_seen'])
                        else:
                            attributes = declaration.get('attributes') or []
                            type_text  = declaration.get('type') or ''
                            if any(a == 'pointer' for a in attributes):
                                continue
                            if re.match(r'^\s*omp_lock_kind\s*', type_text):
                                continue
                            if any(a == 'allocatable' for a in attributes):
                                frag = _process_allocatable_intrinsic(
                                    declaration, exclude)
                                output_body += frag['output']
                                input_body  += frag['input']
                                stored_shape_required |= frag['stored_shape_required']
                                was_allocated_required |= frag['was_allocated_required']
                                label_required |= frag['label_required']
                                rank_maximum = max(rank_maximum, frag['rank_seen'])
                            else:
                                # Statically-sized variables — batched later.
                                for v in declaration.get('variables', []):
                                    name = _strip_init(v)
                                    if any(x.lower() == name.lower() for x in exclude):
                                        continue
                                    restore_to = (
                                        scope_directive.get('restoreTo')
                                        if scope_directive is not None else None
                                    )
                                    skip_store = False
                                    for rt in as_array(restore_to):
                                        if not isinstance(rt, dict):
                                            continue
                                        target_vars = [
                                            x.strip() for x in re.split(
                                                r'\s*,\s*', rt.get('variables', ''))
                                            if x.strip()
                                        ]
                                        if any(x.lower() == name.lower() for x in target_vars):
                                            skip_store = True
                                            input_body += (
                                                f" self%{name}={rt.get('state', '')}\n"
                                            )
                                    if not skip_store:
                                        static_variables.append(name)
                child = child.get('sibling')

            # methodCalls from this level.
            if scope_directive is not None and 'methodCall' in scope_directive:
                for mc in as_array(scope_directive.get('methodCall')):
                    if not isinstance(mc, dict):
                        continue
                    args = mc.get('arguments', '')
                    method_calls.append(
                        f"  call self%{mc.get('method', '')}({args})"
                    )
            cursor = classes[cursor].get('extends')

        for name in static_variables:
            label_required = True
            output_body += " if (displayVerbosity() >= verbosityLevelWorking) then\n"
            output_body += f"  write (label,'(i16)') sizeof(self%{name})\n"
            output_body += (
                f"  call displayMessage('storing \"{name}\" with size '"
                "//trim(adjustl(label))//' bytes')\n"
            )
            output_body += " end if\n"
        for name in static_variables:
            input_body += (
                f" call displayMessage('restoring \"{name}\"',"
                "verbosity=verbosityLevelWorking)\n"
            )
        if static_variables:
            output_body += (
                "  write (stateFile) "
                + ", &\n  & ".join(f"self%{n}" for n in static_variables)
                + "\n"
            )
            input_body += (
                "  read  (stateFile) "
                + ", &\n  & ".join(f"self%{n}" for n in static_variables)
                + "\n"
            )
        if method_calls:
            input_body += "\n".join(method_calls) + "\n"

    # Grow openers with the locals we need.
    if label_required:
        output_opener += " character(len=16      )                             :: label\n"
    if stored_shape_required:
        input_opener += (
            " integer(c_size_t    ), allocatable  , dimension(:) :: storedShape\n"
        )
    if was_allocated_required:
        input_opener += (
            " logical                                            :: wasAllocated\n"
        )
    if rank_maximum > 0:
        idx_list = ", ".join(f"i{i}" for i in range(1, rank_maximum + 1))
        output_opener += (
            f" integer  (c_size_t    )                             :: {idx_list}\n"
        )
        input_opener += (
            f" integer  (c_size_t    )                             :: {idx_list}\n"
        )

    # Unused-variable annotations.  gslStateFile is never consumed directly
    # — the generated body only passes it along to inner stateStore calls.
    output_unused = ['gslStateFile']
    input_unused  = ['gslStateFile']
    if not label_required:
        output_unused.append('label')
    output_unused_line = (
        " !$GLC attributes unused :: " + ", ".join(output_unused) + "\n"
    )
    input_unused_line = (
        " !$GLC attributes unused :: " + ", ".join(input_unused) + "\n"
    )

    function_code = (
        output_opener
        + output_unused_line
        + output_body
        + output_closer
        + "\n"
        + input_opener
        + input_unused_line
        + input_body
        + input_closer
        + "\n"
    )

    # Class-restore dispatchers for every ancestor/descendant in the target
    # hierarchy.  Assembled outside the main function_code so the caller can
    # decide whether to append them (they are only useful when the generated
    # code is inside a module so the `SetVisibility` call is meaningful).
    restore_text = ''
    if class_identifiers:
        for parent_class in sorted(classes.keys()):
            if not _descends_from(parent_class, target, classes):
                continue
            restore_text += _emit_class_restore(
                parent_class, target, classes, class_identifiers)

    return function_code, class_functions_used, restore_text


register_process('stateStorable', process_state_storable, before=['generics'])
