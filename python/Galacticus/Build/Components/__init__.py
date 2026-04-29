# Driver for the components-build pipeline.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components.pm: the top-level `component`
# handler invoked by scripts/build/buildCode.py.  The driver itself is
# small — it owns three things:
#
#   1. `validate` — XSD-validate the directive body against
#      schema/componentSchema.xsd.
#   2. `parse`    — store each `<component>` directive's parsed body in
#      `build['components']`, keyed by `<class><name>`.
#   3. `generate` — run the phased hook pipeline (preValidate, default,
#      gather, scatter, postValidate, content, types, interfaces,
#      functions), then serialise the accumulated derived types,
#      interfaces, module-scope variables, and functions into the
#      Fortran source held on `build['content']`.
#
# Sub-modules under Galacticus.Build.Components.* register their own
# functions on `Utils.component_utils` at import time; the driver simply
# walks that registry.

import os
import re
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from build.file_changes        import update as file_changes_update
from Galacticus.Build          import Hooks
from Galacticus.Build.Components import Utils
from Galacticus.Build.Components.CodeGeneration import (
    function_arguments,
    importables,
)

# Import sister modules eagerly so they register their own hooks at the
# same time the framework itself registers `component` with `Hooks`.
# Mirrors the chain of `use Galacticus::Build::Components::*;` lines in
# perl/Galacticus/Build/Components.pm.
from Galacticus.Build.Components import (  # noqa: F401
    DataTypes,
    NullFunctions,
    Attributes,
    Components as _Components,
    Classes,
    TreeNodes,
)
from Galacticus.Build.Components.Classes   import Utils as _ClassesUtils    # noqa: F401
from Galacticus.Build.Components.TreeNodes import State          as _TreeNodesState           # noqa: F401
from Galacticus.Build.Components.TreeNodes import Classes        as _TreeNodesClasses         # noqa: F401
from Galacticus.Build.Components.TreeNodes import CreateDestroy  as _TreeNodesCreateDestroy   # noqa: F401
from Galacticus.Build.Components.TreeNodes import Map            as _TreeNodesMap             # noqa: F401


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

def _register():
    """Register the component handler with the global hooks registry.

    Imported sister modules add themselves to `Utils.component_utils` at
    their own import time; the registration here is the build's `validate`
    / `parse` / `generate` entry points.
    """
    Hooks.module_hooks['component'] = {
        'validate': validate,
        'parse':    parse_directive,
        'generate': generate_output,
    }


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

def validate(document_string, file_name):
    """Validate `document_string` against schema/componentSchema.xsd.

    Mirrors `Components_Validate`.  Uses `lxml.etree.XMLSchema`; if `lxml`
    is unavailable we silently skip validation (matching the same fallback
    pattern as scripts/aux/validateParameters.py).
    """
    try:
        from lxml import etree as lxml_etree
    except ImportError:
        return
    schema_path = os.path.join(
        os.environ['GALACTICUS_EXEC_PATH'], 'schema', 'componentSchema.xsd',
    )
    schema_doc = lxml_etree.parse(schema_path)
    schema     = lxml_etree.XMLSchema(schema_doc)
    document   = lxml_etree.fromstring(document_string.encode('utf-8'))
    try:
        schema.assertValid(document)
    except lxml_etree.DocumentInvalid as exc:
        sys.exit(
            f"Galacticus::Build::Components::Components_Validate(): "
            f"validation failed in file {file_name}:\n{exc}"
        )


# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------

def parse_directive(build):
    """Store the current `<component>` directive on `build['components']`.

    Keys are `<class><name>` with each side `ucfirst`-ed, matching Perl's
    `ucfirst($class).ucfirst($name)`.  Duplicate IDs are fatal (Perl
    `die`s in the same spot).
    """
    document = build.get('currentDocument')
    if document is None:
        sys.exit("Galacticus::Build::Components::Components_Parse_Directive: "
                 "no currentDocument present")
    if 'name' not in document:
        sys.exit("Galacticus::Build::Components::Components_Parse_Directive: "
                 "no name present")
    if 'class' not in document:
        sys.exit("Galacticus::Build::Components::Components_Parse_Directive: "
                 "no class present")

    component_id = _ucfirst(document['class']) + _ucfirst(document['name'])
    components   = build.setdefault('components', {})
    if component_id in components:
        sys.exit(
            "Galacticus::Build::Components::Components_Parse_Directive: "
            f"multiple components with ID '{component_id}'"
        )
    components[component_id] = document


# ---------------------------------------------------------------------------
# Generate
# ---------------------------------------------------------------------------

# Phase order matches Perl Components_Generate_Output (Components.pm:122).
_PHASES = (
    'preValidate',
    'default',
    'gather',
    'scatter',
    'postValidate',
    'content',
    'types',
    'interfaces',
    'functions',
)


def generate_output(build):
    """Run the phased hook pipeline and serialise the result.

    Each phase calls every hook registered under that phase, in the
    sub-module-registration order returned by `sorted()` over the owner
    keys (matching Perl `sortedKeys`).  Output is appended to
    `build['content']`, which the buildCode.py driver writes to disk.
    """
    build.setdefault('content',   '')
    build.setdefault('types',     {})
    build.setdefault('functions', [])
    build.setdefault('variables', [])
    build.setdefault('interfaces', {})

    hooks = sorted(Utils.component_utils.items(), key=lambda kv: kv[0])

    print("--> Phase:")
    for phase in _PHASES:
        print(f"   --> {_ucfirst(phase)}...")
        for owner, owner_hooks in hooks:
            functions = owner_hooks.get(phase) or []
            if isinstance(functions, list):
                fns = functions
            else:
                fns = [functions]
            for fn in fns:
                marker = (
                    f" {{{getattr(fn, '__name__', '<fn>')}}}"
                    if len(fns) > 1 else ''
                )
                print(f"      --> {owner}{marker}")
                fn(build)

    derived_types_serialize(build)
    interfaces_serialize   (build)

    if build['variables']:
        build['content'] += format_variable_definitions(build['variables']) + '\n'

    build['content'] += "contains\n\n"

    functions_serialize(build)

    # Include statements for every per-component "functions" file.
    include_dependencies = [
        c['functions']
        for c in build.get('components', {}).values()
        if 'functions' in c
    ]
    build['content'] += '\n'.join(
        f'  include "{p}"\n' for p in include_dependencies
    ) + '\n'

    # Per-build Makefile dependency snippet.
    build_path  = os.environ['BUILDPATH']
    makefile    = os.path.join(build_path, 'Makefile_Component_Includes')
    tmp_path    = makefile + '.tmp'
    with open(tmp_path, 'w') as fh:
        if include_dependencies:
            target = (
                f"{build_path}/objects.nodes.o "
                f"{build_path}/objects.nodes.p.F90:"
                + ''.join(f' {build_path}/{p}' for p in include_dependencies)
            )
            fh.write(target)
    file_changes_update(makefile, tmp_path)


# ---------------------------------------------------------------------------
# Serialisers
# ---------------------------------------------------------------------------

def derived_types_serialize(build):
    """Topo-sort `build['types']` by `extends` + non-pointer member types
    and append each `type … end type` block to `build['content']`.
    """
    from Sort.Topo import sort as topo_sort

    types        = build.get('types') or {}
    dependencies = {}
    for type_name, type_def in types.items():
        # A type depends on its parent.
        parent = type_def.get('extends')
        if parent:
            dependencies.setdefault(type_name, []).append(parent)
        # …and on any member that is itself one of our generated types,
        # unless that member is held by pointer (in which case the parent
        # type only needs the forward declaration that comes for free).
        for member in type_def.get('dataContent') or []:
            intrinsic = member.get('intrinsic')
            if intrinsic not in ('type', 'class'):
                continue
            mtype = member.get('type')
            if mtype not in types:
                continue
            attrs = member.get('attributes') or []
            if 'pointer' in attrs:
                continue
            dependencies.setdefault(type_name, []).append(mtype)

    type_order = topo_sort(sorted(types.keys()), dependencies)

    for name in type_order:
        type_def = types[name]
        line     = "  type"
        if type_def.get('isPublic'):
            line += ", public"
        if type_def.get('extends'):
            line += f", extends({type_def['extends']})"
        line += f" :: {type_def['name']}\n"
        build['content'] += line

        comment = type_def.get('comment')
        if comment:
            build['content'] += "  !!{\n"
            build['content'] += f"  {comment}\n"
            build['content'] += "  !!}\n"

        build['content'] += "    private\n"

        if 'dataContent' in type_def:
            build['content'] += format_variable_definitions(type_def['dataContent'])

        if type_def.get('boundFunctions'):
            build['content'] += "   contains\n"
            build['content'] += bound_function_table(
                type_def['name'], type_def['boundFunctions']
            )

        build['content'] += f"  end type {type_def['name']}\n\n"


def interfaces_serialize(build):
    """Append every abstract interface in `build['interfaces']` to
    `build['content']`.

    The Perl version uses `Text::Template`'s embedded-Perl form to render
    each interface; we explicitly compute the parts we need (subroutine vs.
    function, return-name conventions, importables) and string-build the
    block, since lifting Perl expression evaluation into Python isn't
    worth the complexity.
    """
    print("   --> Serialize interfaces...")
    for interface in (build.get('interfaces') or {}).values():
        print(f"      ---> {interface['name']}")
        is_void    = interface['intrinsic'] == 'void'
        result_yes = interface.get('result') == 'yes'
        if is_void:
            head = "subroutine"
        elif result_yes:
            head = " function"
        else:
            head = f"{interface['intrinsic']} function"
        end_kw = "subroutine" if is_void else "function"

        args = function_arguments(interface.get('data') or [])
        imp  = importables       (interface.get('data') or [])

        block  = f"! {interface.get('comment','')}\n"
        block += "abstract interface\n"
        block += f"  {head} {interface['name']}({','.join(args)})\n"
        if imp:
            block += f"    import {', '.join(imp)}\n"
        block += format_variable_definitions(
            interface.get('data') or [], indent=4,
        )
        block += f"  end {end_kw} {interface['name']}\n"
        block += "end interface\n"
        build['content'] += block + '\n'


def functions_serialize(build):
    """Append every entry on `build['functions']` (and every bound function
    descriptor referenced from `build['types']`) to `build['content']`.
    """
    print("   --> Serialize functions...")

    seen = []
    for fn in build.get('functions') or []:
        seen.append(fn)
    for type_def in (build.get('types') or {}).values():
        for binding in type_def.get('boundFunctions') or []:
            if 'descriptor' in binding:
                seen.append(binding['descriptor'])

    for function in seen:
        print(f"      --> {function['name']}")

        return_name = None
        # `type => name`-style return signature.  Lift the parsed parts
        # explicitly so we can emit the matching declaration.
        m = re.match(
            r'^([a-zA-Z0-9_(),:=\s]+)\s+=>\s+([a-zA-Z0-9_]+)',
            function['type'],
        )
        if m:
            return_descriptor = m.group(1).strip()
            return_name       = m.group(2)
            result            = f"result({return_name})"
            type_text         = ""
            form              = "function"

            attr_match = re.match(
                r'([a-zA-Z0-9_=()\s]+)\s*,\s*([a-zA-Z0-9_,():=\s]+)',
                return_descriptor,
            )
            if attr_match:
                return_type = attr_match.group(1).strip()
                attributes  = [a.strip() for a in attr_match.group(2).split(',')]
            else:
                return_type = return_descriptor
                attributes  = []

            declaration = {
                'attributes': attributes,
                'variables':  [return_name],
            }
            tm = re.match(
                r'^(type|class)\s*\(\s*([a-zA-Z0-9_]+)\s*\)',
                return_type,
            )
            if tm:
                declaration['intrinsic'] = tm.group(1)
                declaration['type'     ] = tm.group(2)
            else:
                declaration['intrinsic'] = return_type
            function.setdefault('variables', []).append(declaration)
        else:
            form      = "subroutine" if function['type'] == 'void' else "function"
            type_text = ""           if function['type'] == 'void' else function['type']
            result    = ""

        function_attributes = []
        if function.get('recursive'):
            function_attributes.append('recursive')

        arguments = Utils.argument_list(*(function.get('variables') or []))

        opener = (
            f"{' '.join(function_attributes)} {type_text} {form} "
            f"{function['name']}({','.join(arguments)}) {result}\n"
        )
        build['content'] += opener
        build['content'] += "   !!{\n"
        build['content'] += f"   {function.get('description', '')}\n"
        build['content'] += "   !!}\n"

        intrinsic_modules = ('iso_c_binding',)
        for module in function.get('modules') or []:
            intrinsic = ", intrinsic" if module.lower() in intrinsic_modules else ""
            build['content'] += f"   use{intrinsic} :: {module}\n"

        build['content'] += "   implicit none\n"
        if function.get('variables'):
            build['content'] += format_variable_definitions(function['variables'])

        if function.get('content'):
            build['content'] += function['content']

        build['content'] += "   return\n"
        build['content'] += f"end {form} {function['name']}\n\n"


def bound_function_table(object_name, bindings):
    """Render a list of type-bound function bindings as a Fortran block.

    Mirrors Perl `boundFunctionTable` (Components.pm:153-238) including its
    optional `<methods>…</methods>` description preamble.  We don't try to
    reproduce Text::Table's column-aligned layout — the Fortran compiler
    accepts any reasonable spacing.
    """
    enriched = []
    for b in bindings:
        # Match Perl's `cmp` semantics: descriptor name takes precedence;
        # otherwise stringify `function` (which may be a list for `generic`
        # bindings).  For lists we join with commas so the sort is at least
        # deterministic — Perl stringified array refs to `ARRAY(0xHEX)`,
        # making the order effectively random there.
        if 'descriptor' in b:
            function_name = b['descriptor']['name']
        else:
            target = b.get('function')
            if isinstance(target, list):
                function_name = ','.join(target)
            else:
                function_name = target or ''
        enriched.append((function_name, b))
    enriched.sort(key=lambda pair: pair[0])

    rows = []
    for _, binding in enriched:
        bind_type = binding.get('type', '')
        pass_text = (', ' + binding['pass']) if 'pass' in binding else ''
        name      = binding.get('name', '')
        connector = ' => ' if 'name' in binding else ''

        target = binding.get('function')
        if isinstance(target, list):
            sorted_targets = sorted(target)
            for i, fn_name in enumerate(sorted_targets):
                suffix = ', &' if i < len(sorted_targets) - 1 else ''
                if i == 0:
                    rows.append(
                        f"     {bind_type}{pass_text} :: {name}{connector}{fn_name}{suffix}"
                    )
                else:
                    rows.append(f"     &{fn_name}{suffix}")
        elif 'descriptor' in binding:
            rows.append(
                f"     {bind_type}{pass_text} :: {name}{connector}"
                f"{binding['descriptor']['name']}"
            )
        else:
            rows.append(
                f"     {bind_type}{pass_text} :: {name}{connector}{target}"
            )

    description       = ""
    method_count      = 0
    for _, binding in enriched:
        text = None
        if 'descriptor' in binding:
            text = binding['descriptor'].get('description')
        elif 'description' in binding:
            text = binding['description']
        if text is None:
            continue
        method_count += 1
        if 'descriptor' in binding and 'methodName' in binding['descriptor']:
            method_name = binding['descriptor']['methodName']
        else:
            method_name = binding.get('name', '')
        description += (
            f"      <method method=\"{method_name}\" "
            f"description=\"{text}\"/>\n"
        )

    out = ""
    if method_count >= 1:
        out += "     !![\n     <methods>\n" + description \
             + "     </methods>\n     !!]\n\n"
    out += '\n'.join(rows)
    if rows:
        out += '\n'
    return out


# ---------------------------------------------------------------------------
# Variable-definition formatter
#
# Stand-in for Perl `Fortran::Utils::Format_Variable_Definitions` — the
# alignment-rich Text::Table output that the Perl tooling builds.  The
# downstream Fortran compiler does not care about column alignment, so we
# produce one declaration per logical variable group.  This is sufficient
# for the components build whose only consumer is the Fortran source
# preprocessor.
# ---------------------------------------------------------------------------

def format_variable_definitions(declarations, indent=2):
    """Return Fortran declaration lines for every dict in `declarations`.

    Honours the keys actually used by the components pipeline: `intrinsic`,
    `type`, `attributes`, `variables`, plus an optional `comment` field for
    a trailing `! …` comment.
    """
    lines = []
    pad   = ' ' * indent
    for decl in declarations or []:
        if not isinstance(decl, dict):
            continue
        intrinsic = decl.get('intrinsic') or ''
        type_text = decl.get('type')
        attrs     = decl.get('attributes') or []
        variables = decl.get('variables')  or []

        line = pad + intrinsic
        if type_text is not None:
            wrapped = (
                type_text
                if (type_text.startswith('(') and type_text.endswith(')'))
                else f'({type_text})'
            )
            line += wrapped
        for attr in attrs:
            line += f', {attr}'
        line += ' :: ' + ', '.join(variables)
        if decl.get('comment'):
            line += f' ! {decl["comment"]}'
        lines.append(line + '\n')
    return ''.join(lines)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# Register handlers eagerly so `import Galacticus.Build.Components` is
# enough to wire the dispatcher up — matches Perl's load-time
# `%Galacticus::Build::Hooks::moduleHooks = (..., component => {...})`.
_register()
