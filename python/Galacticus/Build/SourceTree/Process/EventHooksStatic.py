# Processes `eventHookStatic` directives: at each call site, synthesizes the
# `call X(callWith)` ladder over every function that has opted into the
# event via the named directive, and imports each function from its
# providing module.  Also marks every hooked function public in its
# enclosing module.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/EventHooksStatic.pm

import os
import pickle
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from List.ExtraUtils                              import as_array
from XML.Utils                                    import xml_to_dict
from Galacticus.Build.Directives                  import extract_directives
from Galacticus.Build.Dependencies                import dependency_sort
from Galacticus.Build.SourceTree                  import (
    walk_tree, parse_file, insert_after_node, set_visibility,
)
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


_EVENT_HOOK_NAMES = None
_DIRECTIVE_LOCATIONS = None
_STATE_STORABLES = None


def _load_xml_once():
    global _DIRECTIVE_LOCATIONS, _STATE_STORABLES
    if _DIRECTIVE_LOCATIONS is not None and _STATE_STORABLES is not None:
        return _DIRECTIVE_LOCATIONS, _STATE_STORABLES
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_event_hooks_static: BUILDPATH is not set")
    if _DIRECTIVE_LOCATIONS is None:
        path = os.path.join(build_path, 'directiveLocations.xml')
        _DIRECTIVE_LOCATIONS = xml_to_dict(ET.parse(path).getroot())
    if _STATE_STORABLES is None:
        path = os.path.join(build_path, 'stateStorables.xml')
        _STATE_STORABLES = xml_to_dict(ET.parse(path).getroot())
    return _DIRECTIVE_LOCATIONS, _STATE_STORABLES


def _load_event_hook_names(directive_locations):
    """Return the set of all `<eventHookStatic>` names discovered across every
    file listed under `directiveLocations/eventHookStatic/file`.

    Uses a `$BUILDPATH/eventHooksStatic.blob` pickle cache between runs,
    mirroring the Perl `Storable`-based cache at EventHooksStatic.pm:33-39.
    """
    global _EVENT_HOOK_NAMES
    if _EVENT_HOOK_NAMES is not None:
        return _EVENT_HOOK_NAMES

    build_path = os.environ.get('BUILDPATH')
    blob_path = (os.path.join(build_path, 'eventHooksStatic.blob')
                 if build_path else None)
    if blob_path and os.path.exists(blob_path):
        with open(blob_path, 'rb') as fh:
            _EVENT_HOOK_NAMES = list(pickle.load(fh))
        return _EVENT_HOOK_NAMES

    block = (directive_locations.get('eventHookStatic') or {})
    files = list(as_array(block.get('file')))
    names = []
    for f in files:
        for d in extract_directives(f, 'eventHookStatic'):
            if 'name' in d:
                names.append(d['name'])
    _EVENT_HOOK_NAMES = names
    if blob_path:
        try:
            with open(blob_path, 'wb') as fh:
                pickle.dump(names, fh)
        except OSError:
            pass
    return _EVENT_HOOK_NAMES


def _function_class_entry_for(type_name, state_storables):
    """Return the functionClass entry (dict) for `<type>+'Class'` or None.

    Matches the lookup `$stateStorables->{'functionClasses'}{$type.'Class'}`
    from EventHooksStatic.pm:84 / :90 / :157.  Our xml_to_dict produces a
    list of entries where the Perl KeyAttr grouping would have produced a
    name-keyed dict; bridge accordingly.
    """
    fc = (state_storables or {}).get('functionClasses') or {}
    target = type_name + 'Class'
    if not isinstance(fc, dict):
        return None
    entries = fc.get('functionClass')
    if entries is None:
        val = fc.get(target)
        return val if isinstance(val, dict) else None
    if isinstance(entries, dict):
        entries = [entries]
    for e in entries:
        if isinstance(e, dict) and e.get('name') == target:
            return e
    return None


def _resolve_hooked_function(source_file, target_directive_name, state_storables):
    """Walk one source file to find the module name + hooked-function metadata.

    Mirrors the inner loop at EventHooksStatic.pm:58-102.  Returns a dict
    shaped like `{function, module, after, before}`, or raises on missing
    info (the Perl code dies in the same spot).
    """
    tree = parse_file(source_file)
    module_name    = None
    function_name  = None
    function_class = None
    after          = None
    before         = None
    use_global     = None

    for n in walk_tree(tree):
        ntype = n.get('type')
        if ntype == 'module':
            module_name = n.get('name')
        if ntype == target_directive_name:
            d = n.get('directive') or {}
            function_name = d.get('function')
            if 'after'     in d: after      = d['after']
            if 'before'    in d: before     = d['before']
            if 'useGlobal' in d: use_global = d['useGlobal']
        # functionClass instance lookup.
        if _function_class_entry_for(ntype or '', state_storables) is not None:
            function_class = (ntype or '') + 'Class'

    if module_name is None and function_class is not None:
        fc_entry = _function_class_entry_for(
            function_class[:-len('Class')], state_storables)
        if fc_entry is not None:
            module_name = fc_entry.get('module')

    if module_name is None:
        raise RuntimeError(
            "process_event_hooks_static: unable to find module containing "
            "hooked function")
    if function_name is None:
        raise RuntimeError(
            "process_event_hooks_static: unable to find name of hooked function")

    if use_global == 'yes':
        function_name += '_'
        module_name    = 'Functions_Global'

    return {
        'function': function_name,
        'module':   module_name,
        'after':    after,
        'before':   before,
    }


def _dependency_sorted(hooked_functions):
    """Return `hooked_functions` reordered per `after`/`before` constraints.

    Matches EventHooksStatic.pm:103-119: first sort alphabetically by
    `function`, then if any entry carries after/before, topo-sort by those.
    """
    hooked_functions = sorted(hooked_functions, key=lambda h: h['function'])
    has_constraints = any(
        h.get('after') is not None or h.get('before') is not None
        for h in hooked_functions
    )
    if not has_constraints:
        return hooked_functions

    tasks = {}
    for h in hooked_functions:
        task = {}
        if h.get('after')  is not None: task['after']  = h['after']
        if h.get('before') is not None: task['before'] = h['before']
        tasks[h['function']] = task
    order = dependency_sort(tasks)

    by_name = {h['function']: h for h in hooked_functions}
    return [by_name[name] for name in order if name in by_name]


def process_event_hooks_static(tree, options):
    """Mirrors Process_EventHooksStatic() from EventHooksStatic.pm."""
    directive_locations, state_storables = _load_xml_once()
    event_hook_names = _load_event_hook_names(directive_locations)

    module_node        = None
    is_function_class  = False
    hooked_functions_in_file = []

    for node in walk_tree(tree):
        ntype = node.get('type', '')

        if ntype == 'module':
            module_node = node

        if ntype == 'eventHookStatic':
            directive = node.get('directive') or {}
            if directive.get('processed'):
                continue
            directive['processed'] = True

            name = directive['name']
            loc_block = directive_locations.get(name) or {}
            locations = list(as_array(loc_block.get('file')))

            hooked = []
            for loc in locations:
                hooked.append(_resolve_hooked_function(
                    loc, name, state_storables))
            hooked = _dependency_sorted(hooked)

            call_with = directive.get('callWith', '')
            on_return = directive.get('onReturn', '')
            call_lines = []
            for hf in hooked:
                call_lines.append(f"call {hf['function']}({call_with})")
                if on_return:
                    call_lines.append(on_return)
            content = "\n".join(call_lines) + ("\n" if call_lines else '')

            insert_after_node(node, [{
                'type':       'code',
                'content':    content,
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':
                    'Galacticus.Build.SourceTree.Process.EventHooksStatic'
                    '.process_event_hooks_static()',
                'line':       1,
            }])

            module_uses = {
                hf['module']: {'intrinsic': False,
                               'only': {hf['function']: True}}
                for hf in hooked
            }
            if module_uses:
                add_uses(node['parent'], {
                    'moduleUse':   module_uses,
                    'moduleOrder': sorted(module_uses.keys()),
                })

        # Collect function names declared via any registered static hook.
        if ntype and ntype in event_hook_names:
            directive = node.get('directive') or {}
            if not directive.get('processed'):
                directive['processed'] = True
                fn = directive.get('function')
                if fn:
                    hooked_functions_in_file.append(fn)

        # Track whether this file hosts a functionClass instance.
        if _function_class_entry_for(ntype, state_storables) is not None:
            is_function_class = True

    if hooked_functions_in_file:
        if module_node is not None:
            for fn in hooked_functions_in_file:
                set_visibility(module_node, fn, 'public')
        elif not is_function_class:
            label  = ('functions' if len(hooked_functions_in_file) > 1
                      else 'function')
            verb   = ('are'       if len(hooked_functions_in_file) > 1 else 'is')
            names  = ", ".join(f"'{n}'" for n in hooked_functions_in_file)
            raise RuntimeError(
                f"process_event_hooks_static: hooked {label} {names} "
                f"{verb} not in a module")


register_process('eventHooksStatic', process_event_hooks_static)
