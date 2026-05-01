# Processes `generic` directives: duplicates the following sibling subtree
# once per `<instance>`, substituting generic placeholders of the form
# `{identifier¦modifier}` or `{identifier¦modifier¦regex_from}` with the
# instance's attribute values.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/Generics.pm

import copy
import io
import os
import re
import sys


from List.ExtraUtils                                  import as_array
from Galacticus.Build.SourceTree                      import (
    walk_tree, replace_node, parse_code, serialize,
)
from Galacticus.Build.SourceTree.Process              import register_process, process_tree

# The broken-bar character used as field separator in generic placeholders.
# Defined as a constant to make the intent obvious at call sites.
_SEP = '¦'   # ¦


# ---------------------------------------------------------------------------
# Substitution helpers
# ---------------------------------------------------------------------------

def _dollar_to_backref(template):
    r"""Convert Perl-style `$1`…`$9` backrefs in a `to` string to Python `\1`…`\9`.

    The Perl code performs `eval qq{"$to"}` which interpolates scalar variables,
    and in practice the only variables used are capture references like `$1`.
    Translate them to Python regex backreferences for use with `re.sub`.
    Any literal `\` is escaped first so it isn't interpreted by `re.sub`.
    """
    # Escape existing backslashes first so they survive re.sub.
    out = template.replace('\\', '\\\\')
    # $1..$9  →  \1..\9
    out = re.sub(r'\$([1-9])', r'\\\1', out)
    return out


def _replace_generic(line, identifier, instance, modifier):
    """Substitute `{identifier¦modifier}` (and bare `{identifier}` when
    modifier is absent) or `{identifier¦modifier¦FROM}` with the instance's
    attribute value.

    Mirrors Perl ReplaceGeneric() at Generics.pm:150-167.
    """
    value = instance[modifier]
    m_regex = re.match(r'^regEx' + _SEP + r'(.*)' + _SEP + r'(.*)' + _SEP, value)
    if m_regex:
        # regEx-form modifier: the instance declares FROM and TO patterns.
        # Build the outer regex around the identifier¦modifier tag and apply
        # the instance's FROM/TO to inline matches.
        from_pat = m_regex.group(1)
        to_pat   = _dollar_to_backref(m_regex.group(2))
        pattern  = re.compile(
            r'\{' + re.escape(identifier) + _SEP + re.escape(modifier) + _SEP
            + from_pat + r'\}', re.DOTALL)
        return pattern.sub(to_pat, line)
    # Plain substitution: `{id¦mod}` (or `{id}`) → instance[modifier]
    pattern = re.compile(
        r'\{' + re.escape(identifier)
        + r'(?:' + _SEP + re.escape(modifier) + r')?\}')
    return pattern.sub(lambda _m: str(value), line)


def _replace_generic_conditional(line, identifier, instance):
    """Substitute `{identifier¦match¦REGEX¦MATCH_TEXT¦NOMATCH_TEXT}` with
    either MATCH_TEXT or NOMATCH_TEXT depending on whether REGEX matches the
    instance's `label`.

    Mirrors Perl ReplaceGenericConditional() at Generics.pm:169-187.
    """
    pattern = re.compile(
        r'\{' + re.escape(identifier) + _SEP + r'match' + _SEP
        + r'([^' + _SEP + r']*)' + _SEP
        + r'([^' + _SEP + r']*)' + _SEP
        + r'([^' + _SEP + r']*)\}')
    while True:
        m = pattern.search(line)
        if not m:
            break
        instance_regex = re.compile(m.group(1))
        match_text     = m.group(2)
        no_match_text  = m.group(3)
        replacement = match_text if instance_regex.search(instance.get('label', '')) else no_match_text
        line = line[:m.start()] + replacement + line[m.end():]
    return line


def _substitute_all(text, identifier, instance):
    """Apply `_replace_generic` for every modifier in instance plus the final
    `_replace_generic_conditional` pass.  Matches the two-step loop in
    Generics.pm:53-57 and :65-67.
    """
    for modifier in sorted(instance.keys()):
        text = _replace_generic(text, identifier, instance, modifier)
    text = _replace_generic_conditional(text, identifier, instance)
    return text


def _substitute_content_lines(content, identifier, instance):
    """Apply `_substitute_all` line-by-line to content (preserving line boundaries)."""
    out = []
    for line in content.splitlines(keepends=True):
        out.append(_substitute_all(line, identifier, instance))
    return ''.join(out)


# ---------------------------------------------------------------------------
# Tree-walking helpers
# ---------------------------------------------------------------------------

def _stack_it(root, initial_depth=-1):
    """Return a pre-order list of (node, depth) pairs rooted at `root`.

    Mirrors Perl StackIt().  The resulting list, when popped, yields a
    reverse-pre-order traversal (children processed before their parent),
    which is what Generics needs so that replacements at a given depth don't
    invalidate earlier stack entries.
    """
    result = []

    def _recurse(node, depth):
        depth += 1
        result.append((node, depth))
        child = node.get('firstChild')
        while child is not None:
            _recurse(child, depth)
            child = child.get('sibling')

    _recurse(root, initial_depth)
    return result


def _reparse_declaration(code_node):
    """If `code_node` is the raw-text child of a `declaration` node, re-run
    the declarations parser over its content so the structured `declarations`
    list reflects whatever the generic expansion produced.  Mirrors the
    `$copyNode->{'parent'}->{'type'} eq "declaration"` branch at
    Generics.pm:73-79 / :132-138.
    """
    parent = code_node.get('parent')
    if parent is None or parent.get('type') != 'declaration':
        return
    node_copy = copy.deepcopy(code_node)
    node_copy['parent']     = None
    node_copy['sibling']    = None
    tree_tmp = {
        'type':       'null',
        'firstChild': node_copy,
        'sibling':    None,
        'parent':     None,
        'source':     code_node.get('source', 'unknown'),
        'line':       code_node.get('line', 0),
    }
    node_copy['parent'] = tree_tmp
    # Delegate to the (already-ported) declarations parse pass.
    from Galacticus.Build.SourceTree import _pass_declarations, children
    _pass_declarations(tree_tmp)
    # Replace the original declaration parent with EVERY child the parse
    # pass produced — `_pass_declarations` may split the input into a
    # sequence (declaration + code + declaration + …) when preprocessor
    # directives like `#ifdef` separate runs of declarations.  Passing only
    # the first child to `replace_node` (which overwrites the new node's
    # sibling pointer) silently dropped any preprocessor-wrapped
    # declarations such as the `commitHash` array added by
    # ParameterMigration.
    replace_node(parent, children(tree_tmp))


def _opener_matches_generic(node, generic_re):
    opener = node.get('opener')
    return opener is not None and generic_re.search(opener) is not None


def _has_generic_ancestor(node, generic_re):
    p = node.get('parent')
    while p is not None:
        if _opener_matches_generic(p, generic_re):
            return True
        p = p.get('parent')
    return False


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------

def process_generics(tree, options):
    """Mirrors Process_Generics() from Generics.pm."""
    # We iterate by materializing the walk up-front because the walk mutates
    # the tree underneath us (Perl uses Walk_Tree which tolerates this; our
    # walk_tree is a generator that wouldn't).
    for node in list(walk_tree(tree)):
        if node.get('type') != 'generic':
            continue
        directive = node.setdefault('directive', {})
        directive['processed'] = True
        identifier = directive['identifier']
        generic_re = re.compile(
            r'\{' + re.escape(identifier) + _SEP + r'.*\}')
        instances = list(as_array(directive.get('instance')))

        sibling = node.get('sibling')
        tree_name = tree.get('name', '<generic>')
        while sibling is not None:
            # Walk the sibling's subtree top-down (parent before children) so
            # that an *outer* opener with a generic placeholder is expanded as
            # a single unit — its descendants are cloned once-per-instance via
            # the per-clone substitution walk inside `_expand_subtree`, and
            # the descendants must NOT be expanded a second time afterwards.
            #
            # The bottom-up stack pop order would expand inner placeholders
            # first (e.g. the `module procedure {Type¦…}` *inside* an
            # `interface {Type¦…}`), and then the outer expansion would clone
            # the interface — now containing 13 already-expanded modprocs —
            # 13 times, giving every interface a copy of every constructor.
            _expand_top_down(sibling, generic_re, identifier, instances,
                             tree_name)
            sibling = sibling.get('sibling')


def _expand_top_down(node, generic_re, identifier, instances, tree_name):
    """Walk `node` and its subtree top-down, expanding every opener that
    matches `generic_re`.  When an opener matches, the entire subtree at that
    point is cloned once per instance; descendants of the matched node are
    NOT visited a second time (they were already cloned-and-substituted as
    part of the parent's `_expand_subtree`).

    Code-content nodes (no `opener`) are visited unconditionally so inline
    `{Type¦…}` placeholders inside raw Fortran lines still get substituted.
    """
    if _opener_matches_generic(node, generic_re):
        copies = _expand_subtree(node, identifier, instances, tree_name)
        replace_node(node, copies)
        return

    if 'content' in node:
        if not _has_generic_ancestor(node, generic_re):
            _expand_content_lines(node, identifier, instances, generic_re)

    # Recurse into children.  Materialise the child list first because each
    # child may rewrite the parent's child chain via replace_node above.
    children_list = []
    child = node.get('firstChild')
    while child is not None:
        children_list.append(child)
        child = child.get('sibling')
    for c in children_list:
        _expand_top_down(c, generic_re, identifier, instances, tree_name)


def _strip_processed_directive_content(root):
    """Blank the firstChild content of every already-processed directive
    in the subtree at `root`, EXCEPT those whose type is in the
    `is_non_processed_type` exemption list (`<methods>`, `<workaround>`,
    …) — those carry data rather than generating code, so their content
    must survive into the cloned subtree for the inner `process_tree`'s
    classDocumentation pass to find them.

    See `_expand_subtree` for why we strip in the first place.
    """
    from Galacticus.Build.SourceTree.Process.NonProcessed import (
        is_non_processed_type,
    )
    for cn in walk_tree(root):
        directive = cn.get('directive')
        if not (directive and directive.get('processed')):
            continue
        if is_non_processed_type(cn.get('type')):
            continue
        first = cn.get('firstChild')
        if first is not None:
            first['content'] = ''


def _expand_subtree(sub_node, identifier, instances, tree_name):
    """Clone sub_node per instance, apply substitutions, serialize each clone,
    and re-parse it so the resulting subtree has fresh directive / declaration
    structure.  Matches the inner block at Generics.pm:40-90.
    """
    copies = []
    for instance in instances:
        copied = copy.deepcopy(sub_node)
        copied['parent']  = None
        copied['sibling'] = None

        # Walk the cloned subtree, substituting in every opener/closer/name
        # field and in every code-content field.
        for cn in walk_tree(copied):
            for field in ('opener', 'closer', 'name'):
                if field in cn and cn[field] is not None:
                    cn[field] = _substitute_all(cn[field], identifier, instance)
            if 'content' in cn and cn['content'] is not None:
                cn['content'] = _substitute_content_lines(
                    cn['content'], identifier, instance)
                _reparse_declaration(cn)

        # Strip the `!![…!!]` XML markers from every already-processed
        # directive — see `_strip_processed_directive_content` for the
        # exemption rules and rationale.
        _strip_processed_directive_content(copied)

        # Re-parse the serialized copy so generic-expanded declarations /
        # directives are re-discovered, then run the full process pipeline on
        # the isolated subtree so any directive types that were already
        # handled in the outer tree (e.g. `optionalArgument`, `inputParameter`)
        # also get processed in the new copy.  Hooks are idempotent — they
        # check `directive.get('processed')` and skip — so re-running them is
        # safe.  Mirrors `_insert_parsed(..., run_process_tree=True)`.
        # `instrument=False`: the original tree was already source-introspection-
        # instrumented when first parsed, so its serialised form already
        # contains `{introspection:location:NNN}` tags.  Re-instrumenting
        # would re-tag those, replacing baked-in line numbers with the line
        # number of the position they happen to land on in the synthesised
        # text — matches Perl ParseCode's `instrument => 0` option.
        reparsed = parse_code(serialize(copied), name=tree_name, instrument=False)
        process_tree(reparsed)
        copies.append(reparsed)
    return copies


def _expand_content_lines(code_node, identifier, instances, generic_re):
    """For each line in code_node['content'] that contains a generic
    placeholder, emit one substituted copy per instance; other lines are
    copied verbatim.  Mirrors Generics.pm:102-139.
    """
    out = []
    for line in code_node['content'].splitlines(keepends=True):
        if generic_re.search(line):
            for instance in instances:
                out.append(_substitute_all(line, identifier, instance))
        else:
            out.append(line)
    code_node['content'] = ''.join(out)
    _reparse_declaration(code_node)


register_process('generics', process_generics)
