"""Processes `modelFact` and `modelRequire` directives: the declarative
interface to the model-scoped fact registry held by the root parameter
document (see `inputParameters%modelFact*` in
`source/utility/input_parameters.F90`).

A *fact* is a string token (e.g. `satellite:tidalHeatingNormalized:evolved`)
which an object establishes for the whole model simply by being built. An
object whose results are only correct if some *other* object (which it neither
owns nor builds) is present declares a *requirement* on that fact. Requirements
are recorded when the consumer is built and checked at explicit validation
points (`inputParameters%modelFactsValidate()`), so the check is independent of
the order in which producer and consumer are built.

    !![
    <modelFact name="satellite:tidalHeatingNormalized:evolved"/>
    !!]

expands to

    call parameters%modelFactProvide('satellite:tidalHeatingNormalized:evolved')

and

    !![
    <modelRequire name="satellite:tidalHeatingNormalized:evolved"
                  reason="heating of disks requires that the tidal heating of satellites is evolved"/>
    !!]

expands to

    call parameters%modelFactRequire('satellite:tidalHeatingNormalized:evolved','<reason>','<help>',<location>)

where `<help>` is generated at build time from every `modelFact` directive in
the source tree, naming the class and implementation which provides the fact
and the parameter-file XML which selects it. A `modelRequire` naming a fact
which no `modelFact` directive provides is a build-time error, so the token
vocabulary is checked by the build.

Both directives accept `source="<variable>"` to name the `inputParameters`
object to register with (default `parameters`). Directives should be placed in
parameter constructors, which always have the parameters to hand.

Andrew Benson (2026; with assistance from Claude Code)
"""

import os
import re


from List.ExtraUtils                                         import as_array
from Galacticus.Build                                        import Catalogs
from Galacticus.Build.Directives                             import extract_directives
from Galacticus.Build.SourceTree                             import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process                     import register_process
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location
from Galacticus.Build.SourceTree.Process.MetaPropertyDatabase import (
    _function_class_names, _locate_function_class,
)


_SOURCE_TAG    = 'Galacticus.Build.SourceTree.Process.ModelFacts'
_TOKEN_PATTERN = re.compile(r'^[A-Za-z][A-Za-z0-9_]*(:[A-Za-z][A-Za-z0-9_]*)*$')
_NAME_PATTERN  = re.compile(r'^[A-Za-z][A-Za-z0-9_%]*$')

# Cache of the provider table, keyed by build path. The table is a pure
# function of the source tree, so it is computed once per preprocessor run.
_PROVIDERS_CACHE = {}


def _quote(text):
    """Render `text` as a Fortran character literal, doubling any apostrophe."""
    return "'" + text.replace("'", "''") + "'"


def _token(directive, kind, node):
    """Return the validated fact token named by a directive."""
    token = (directive.get('name') or '').strip()
    if not _TOKEN_PATTERN.match(token):
        raise RuntimeError(
            f"{kind}: fact name '{token}' is not a valid token (colon-separated"
            " identifiers, e.g. `satellite:tidalHeatingNormalized:evolved`)"
            f" [{node.get('source', 'unknown')}:{node.get('line', 0)}]")
    return token


def _source(directive, kind, node):
    """Return the name of the `inputParameters` variable to register with."""
    source = (directive.get('source') or 'parameters').strip()
    if not _NAME_PATTERN.match(source):
        raise RuntimeError(
            f"{kind}: `source` attribute '{source}' is not a variable name"
            f" [{node.get('source', 'unknown')}:{node.get('line', 0)}]")
    return source


def provider_table(build_path):
    """Return `{token: [(functionClass, implementation, fileName), ...]}` for
    every `modelFact` directive in the source tree.

    `functionClass` and `implementation` are `None` for a provider which is
    not a functionClass implementation (e.g. a node component).
    """
    if build_path in _PROVIDERS_CACHE:
        return _PROVIDERS_CACHE[build_path]
    locations_path = os.path.join(build_path, 'directiveLocations.xml')
    storables_path = os.path.join(build_path, 'stateStorables.xml')
    for path in (locations_path, storables_path):
        if not os.path.exists(path):
            raise RuntimeError(f"modelRequire: {path} not found")
    directive_locations = Catalogs.load(locations_path)
    class_names         = _function_class_names(Catalogs.load(storables_path))
    block               = directive_locations.get('modelFact') or {}
    providers           = {}
    for file_name in sorted(set(as_array(block.get('file')))):
        if not os.path.exists(file_name):
            # Stale directiveLocations.xml entry for a deleted file.
            continue
        facts = extract_directives(file_name, 'modelFact')
        if not facts:
            continue
        function_class, implementation = _locate_function_class(file_name, class_names)
        for fact in facts:
            token = (fact.get('name') or '').strip()
            entry = (function_class, implementation, file_name)
            if entry not in providers.setdefault(token, []):
                providers[token].append(entry)
    _PROVIDERS_CACHE[build_path] = providers
    return providers


def help_text(token, providers):
    """Compose the HELP text for an unsatisfied requirement on `token`."""
    hints = []
    for function_class, implementation, file_name in providers:
        if function_class:
            hints.append(
                f"the '{implementation}' implementation of the '{function_class}'"
                f" class (add <{function_class} value=\"{implementation}\"/> to"
                f" your parameter file)")
        else:
            hints.append(f"code in {os.path.basename(file_name)}")
    if len(hints) == 1:
        provided_by = hints[0]
    else:
        provided_by = "one of: " + "; ".join(hints)
    return f"the fact '{token}' is established by {provided_by}"


def _providers(options):
    """Return the provider table, from `options` (tests) or the build tree."""
    if 'modelFactProviders' in options:
        return options['modelFactProviders']
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError("modelRequire: BUILDPATH is not set")
    return provider_table(build_path)


def _insert(node, code, function):
    insert_after_node(node, [{
        'type':       'code',
        'content':    code,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     f"{_SOURCE_TAG}.{function}()",
        'line':       1,
    }])


def process_model_fact(tree, options):
    """Process `modelFact` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'modelFact':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True
        token  = _token (directive, 'modelFact', node)
        source = _source(directive, 'modelFact', node)
        _insert(node, f"call {source}%modelFactProvide({_quote(token)})\n", 'process_model_fact')


def process_model_require(tree, options):
    """Process `modelRequire` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'modelRequire':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True
        token  = _token (directive, 'modelRequire', node)
        source = _source(directive, 'modelRequire', node)
        reason = (directive.get('reason') or '').strip()
        if not reason:
            raise RuntimeError(
                "modelRequire: a `reason` attribute is required"
                f" [{node.get('source', 'unknown')}:{node.get('line', 0)}]")
        providers = _providers(options).get(token)
        if not providers:
            raise RuntimeError(
                f"modelRequire: no `modelFact` directive anywhere in the source"
                f" provides the fact '{token}' - check the spelling, or add a"
                " `modelFact` directive to the class which establishes it"
                f" [{node.get('source', 'unknown')}:{node.get('line', 0)}]")
        code = (
            f"call {source}%modelFactRequire("
            f"{_quote(token)},"
            f"{_quote(reason)},"
            f"{_quote(help_text(token, providers))},"
            f"{location(node, node.get('line', 0))})\n"
        )
        _insert(node, code, 'process_model_require')


register_process('modelFact',    process_model_fact   )
register_process('modelRequire', process_model_require)
