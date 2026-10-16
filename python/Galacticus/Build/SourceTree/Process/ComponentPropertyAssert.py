"""Processes `componentPropertyAssert` directives: expands a compact
statement of which properties of a component class must be gettable /
settable / evolvable into the `if (.not. … ) call Error_Report( … //
Component_List( … ) … )` idiom, and wires in the module uses it needs.

The directive replaces boilerplate of the form

    if     (                                                                 &
         &  .not.                                                            &
         &       (                                                           &
         &        defaultSpheroidComponent%    massStellarIsGettable().and.  &
         &        defaultSpheroidComponent%angularMomentumIsGettable()       &
         &  )                                                                &
         & ) call Error_Report                                               &
         &        (                                                          &
         &         '…'//                                                     &
         &         Component_List(                                           &
         &                        'spheroid'                              ,  &
         &                        defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection. &
         &                        defaultSpheroidComponent%angularMomentumAttributeMatch(requireGettable=.true.)               &
         &                       )                                        // &
         &         {introspection:location}                                  &
         &        )

with

    !![
    <componentPropertyAssert class="spheroid" properties="massStellar angularMomentum" require="gettable"/>
    !!]

Three modes are supported:

  * default            — generate both the test and the error report;
  * `assignTo="…"`     — generate only the test, assigning its result to the
                         named logical variable (for classes that cache the
                         result at construction and report later);
  * `condition="…"`    — generate only the error report, guarded by the
                         supplied expression (the partner of `assignTo`).

Andrew Benson (2026; with assistance from Claude Code)
"""

import re


from Galacticus.Build.SourceTree                             import (
    walk_tree, insert_after_node,
)
from Galacticus.Build.SourceTree.Process                     import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses            import add_uses
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


_SOURCE_TAG = ('Galacticus.Build.SourceTree.Process.ComponentPropertyAssert'
               '.process_component_property_assert()')

# Requirement name → the adjective used in both the `<property>Is<Adjective>()`
# accessor and the `require<Adjective>=` argument of `<property>AttributeMatch`.
_REQUIREMENT_ADJECTIVE = {
    'gettable':  'Gettable',
    'settable':  'Settable',
    'evolvable': 'Evolvable',
}

# Canonical ordering of requirements in generated code — chosen so that output
# is independent of the order in which the directive lists them.
_REQUIREMENT_ORDER = ('gettable', 'settable', 'evolvable')

_NAME_PATTERN = re.compile(r'^[A-Za-z][A-Za-z0-9_]*$')


def _split_list(text, attribute, node):
    """Split a whitespace- and/or comma-separated directive attribute."""
    items = [item for item in re.split(r'[\s,]+', (text or '').strip()) if item]
    if not items:
        raise RuntimeError(
            f"componentPropertyAssert: attribute `{attribute}` is empty"
            f" [{node.get('source', 'unknown')}:{node.get('line', 0)}]")
    return items


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


def _quote(text):
    """Render `text` as a Fortran character literal, doubling any apostrophe."""
    return "'" + text.replace("'", "''") + "'"


def _english_list(items):
    """Join `items` as an English list: `a`, `a and b`, `a, b, and c`."""
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return items[0] + ' and ' + items[1]
    return ', '.join(items[:-1]) + ', and ' + items[-1]


def _default_message(class_name, properties, requirements):
    """Compose the diagnostic used when the directive supplies no `message`."""
    adjectives = _english_list(requirements)
    quoted     = _english_list(['"' + p + '"' for p in properties])
    if len(properties) == 1:
        subject = f'a {adjectives} {quoted} property'
    else:
        subject = f'{adjectives} {quoted} properties'
    return (f'the "{class_name}" component must provide {subject}'
            f' for this functionality to be available.')


def _test_expression(component, properties, requirements):
    """Build the `.and.`-joined accessor test, one term per property per
    requirement."""
    terms = []
    for prop in properties:
        for requirement in requirements:
            adjective = _REQUIREMENT_ADJECTIVE[requirement]
            terms.append(f"{component}%{prop}Is{adjective}()")
    return terms


def _match_expression(component, properties, requirements):
    """Build the `.intersection.`-joined `AttributeMatch` list, one term per
    property."""
    arguments = ','.join(
        f"require{_REQUIREMENT_ADJECTIVE[requirement]}=.true."
        for requirement in requirements
    )
    return [
        f"{component}%{prop}AttributeMatch({arguments})"
        for prop in properties
    ]


def _code_node(content, line):
    return {
        'type':       'code',
        'content':    content,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     _SOURCE_TAG,
        'line':       line,
    }


def _statement(first, rest):
    """Render a single Fortran free-form statement.

    `first` is the opening physical line; each entry of `rest` becomes a
    continuation line carrying a leading `&`.
    """
    if not rest:
        return "   " + first + "\n"
    text = "   " + first + " &\n"
    for i, line in enumerate(rest):
        text += "    & " + line + (" &\n" if i < len(rest) - 1 else "\n")
    return text


def process_component_property_assert(tree, options):
    """Process `componentPropertyAssert` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'componentPropertyAssert':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        where = f"{node.get('source', 'unknown')}:{node.get('line', 0)}"

        class_name = (directive.get('class') or '').strip()
        if not _NAME_PATTERN.match(class_name):
            raise RuntimeError(
                f"componentPropertyAssert: `class` must be a component class"
                f" name, got '{class_name}' [{where}]")

        properties = _split_list(directive.get('properties'), 'properties', node)
        for prop in properties:
            if not _NAME_PATTERN.match(prop):
                raise RuntimeError(
                    f"componentPropertyAssert: '{prop}' is not a valid property"
                    f" name [{where}]")

        requested = _split_list(directive.get('require'), 'require', node)
        unknown   = [r for r in requested if r not in _REQUIREMENT_ADJECTIVE]
        if unknown:
            raise RuntimeError(
                f"componentPropertyAssert: unknown requirement(s)"
                f" {', '.join(sorted(unknown))} — expected any of"
                f" {', '.join(_REQUIREMENT_ORDER)} [{where}]")
        # Deduplicate while imposing the canonical ordering.
        requirements = [r for r in _REQUIREMENT_ORDER if r in requested]

        assign_to = (directive.get('assignTo') or '').strip()
        condition = (directive.get('condition') or '').strip()
        if assign_to and condition:
            raise RuntimeError(
                "componentPropertyAssert: `assignTo` and `condition` are"
                f" mutually exclusive [{where}]")

        component = 'default' + _ucfirst(class_name) + 'Component'
        parent    = node['parent']
        line      = node.get('line', 0)

        # Every mode that mentions the component object needs it in scope.
        add_uses(parent, {
            'moduleUse':   {'Galacticus_Nodes': {'intrinsic': False,
                                                 'only': {component: True}}},
            'moduleOrder': ['Galacticus_Nodes'],
        })

        if assign_to:
            # Test-only mode: cache the result here, report it elsewhere.
            terms = _test_expression(component, properties, requirements)
            code  = (f"   ! Auto-generated test of properties of the"
                     f" \"{class_name}\" component.\n")
            code += _statement(f"{assign_to}={terms[0]}",
                               [f".and.{term}" for term in terms[1:]])
            insert_after_node(node, [_code_node(code, line)])
            continue

        # Both remaining modes emit the error report.
        message = directive.get('message') or _default_message(
            class_name, properties, requirements)
        matches = _match_expression(component, properties, requirements)
        test    = (condition if condition
                   else '.and.'.join(
                       _test_expression(component, properties, requirements)))

        report  = [
            "call Error_Report(",
            f"                  {_quote(message)}",
            f"                  //Component_List({_quote(class_name)},",
            f"                                   {matches[0]}",
        ]
        report += [f"                                   .intersection.{match}"
                   for match in matches[1:]]
        report += [
            "                                  )",
            f"                  //{location(node, line)}",
            "                 )",
        ]

        code  = (f"   ! Auto-generated assertion on properties of the"
                 f" \"{class_name}\" component.\n")
        code += _statement(f"if (.not.({test}))", report)

        insert_after_node(node, [_code_node(code, line)])

        module_uses = {
            'Error': {'intrinsic': False,
                      'only': {'Error_Report': True, 'Component_List': True}},
        }
        if len(matches) > 1:
            module_uses['Array_Utilities'] = {
                'intrinsic': False,
                'only':      {'operator(.intersection.)': True},
            }
        add_uses(parent, {
            'moduleUse':   module_uses,
            'moduleOrder': sorted(module_uses.keys()),
        })


register_process('componentPropertyAssert', process_component_property_assert)
