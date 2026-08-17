"""Contains a Python module which implements parsing of module uses in the Galacticus preprocessor system.

Andrew Benson (ported to Python 2026)
"""

import re
import io
import copy


from Galacticus.Build.FortranUtils import get_fortran_line

# Formatting rules for emitted `use` statements.  `SYMBOLS_PER_ROW_MAX` is the
# usual number of `only` symbols per row; `LINE_LENGTH_MAX` is the column
# beyond which `update_uses()` drops to fewer symbols per row.  132 is the
# Fortran free-form standard limit — not enforced here (the build passes
# `-ffree-line-length-none`) but kept as the house style.
SYMBOLS_PER_ROW_MAX = 4
LINE_LENGTH_MAX     = 132

# NB: the `(?![a-zA-Z0-9_])` after `use` is essential — without it, an
# assignment to a variable whose name merely *starts* with "use" (e.g.
# `useCache=lastCache`) parses as `use Cache, only : =lastCache`, fabricating
# a bogus module-use node in the middle of a procedure body (issue #1030).
# The lookahead requires a non-identifier character (whitespace, `,`, or `::`)
# to follow `use`, while still accepting the `use::module` spelling that was
# historically rejected (issue #1030).
_MODULE_USE_RE = re.compile(
    r'^\s*(!\$)?\s*use(?![a-zA-Z0-9_])\s*(,\s*(intrinsic))?\s*(::)?\s*([a-zA-Z0-9_]+)'
    r'\s*(,\s*only\s*:)?\s*([a-zA-Z0-9_()/=*\-+.,\s]+)?\s*$',
    re.IGNORECASE,
)


# Preprocessor directives which `update_uses` is able to *rebuild* from an
# entry's `conditions` list — and therefore the only ones that may be absorbed
# into a moduleUse block.  Each must match the whole line: anything trailing
# (`#endif // USEMPI`) would be lost when the block is rebuilt.  Macro names
# are upper case as elsewhere in the build system; a name outside that set
# (`#ifndef __aarch64__`) is left as ordinary code rather than half-matched.
_PP_IFDEF_RE  = re.compile(r'^#ifdef\s+([A-Z0-9_]+)\s*$' )
_PP_IFNDEF_RE = re.compile(r'^#ifndef\s+([A-Z0-9_]+)\s*$')
_PP_ELSE_RE   = re.compile(r'^#else\s*$'                 )
_PP_ENDIF_RE  = re.compile(r'^#endif\s*$'                )
# Openers and closers we cannot rebuild, but which must still be *paired* so
# that the region structure is tracked correctly.
_PP_IF_ANY_RE    = re.compile(r'^#if(def|ndef)?\b')
_PP_ELSE_ANY_RE  = re.compile(r'^#(else|elif)\b' )
_PP_ENDIF_ANY_RE = re.compile(r'^#endif\b'       )


def _pp_directive(processed_line):
    """Classify a preprocessor line.

    Returns `(kind, name, invert)` where `kind` is one of:

      `'open'`       — `#ifdef NAME` / `#ifndef NAME`, rebuildable as a
                       condition on the `use` statements it guards;
      `'openOther'`  — any other conditional opener (`#if …`), which opens a
                       region but can never be rebuilt;
      `'else'`       — a bare `#else`, which inverts the enclosing condition;
      `'elseOther'`  — `#elif …`, or an `#else` with trailing text;
      `'close'`      — a bare `#endif`;
      `'closeOther'` — an `#endif` with trailing text;
      `'other'`      — a non-conditional directive (`#include`, `#define`, …).

    `name`/`invert` are meaningful only for `kind == 'open'`.
    """
    m = _PP_IFDEF_RE.match(processed_line)
    if m:
        return 'open', m.group(1), False
    m = _PP_IFNDEF_RE.match(processed_line)
    if m:
        return 'open', m.group(1), True
    if _PP_IF_ANY_RE.match(processed_line):
        return 'openOther', None, False
    if _PP_ELSE_RE.match(processed_line):
        return 'else', None, False
    if _PP_ELSE_ANY_RE.match(processed_line):
        return 'elseOther', None, False
    if _PP_ENDIF_RE.match(processed_line):
        return 'close', None, False
    if _PP_ENDIF_ANY_RE.match(processed_line):
        return 'closeOther', None, False
    return 'other', None, False


def _classify_preprocessor_directives(content):
    """Decide which preprocessor directives in a code node may be *claimed*.

    A directive is "claimed" when `parse_module_uses` absorbs it into a
    moduleUse block and records it in the `conditions` of the `use` statements
    it guards — from which `update_uses` rebuilds a matching
    `#ifdef … #endif` wrapper.  That is faithful only when the guarded region
    contains **nothing but `use` statements and other claimed directives**.

    Where the region also covers the code *after* the `use` statements, the
    rebuilt block closes its own wrapper early while the original `#endif`
    survives in the following code node — one `#endif` too many, which fails
    to compile (issue #1385).  The parsed structure round-trips perfectly in
    that case, so nothing downstream can catch it; the region has to be
    rejected here.

    A region is therefore claimable only if all of the following hold:

    * it opens and closes within this code node — a guard that is still open
      when the node ends cannot be closed by the rebuilt block;
    * every line inside it is a `use` statement or a claimable directive —
      any other line (plain code, a comment, a blank line) closes the
      moduleUse block while the guard is still open;
    * its own directives are rebuildable — `#if`, `#elif`, `#include`,
      `#define` and friends are not, and are emitted verbatim as code, which
      likewise splits a block that a surrounding guard was claimed for.

    The last point propagates outwards: an unclaimable region nested inside
    another makes the outer one unclaimable too, because its directives are
    emitted as ordinary code *inside* the outer region.

    Returns `(claimable, enclosing)`, both keyed by logical-line index within
    `content`: `claimable[i]` is `True` for a directive line that may be
    claimed, and `enclosing[i]` is the list of unclaimable regions enclosing
    line `i` (used to tell whether a moduleUse block sits inside a guard that
    is left in place — see `add_uses`).
    """
    claimable = {}
    enclosing = {}
    regions   = []   # open regions, outermost first
    index     = 0

    if '#' not in content:
        # Nothing to classify — skip a second pass over the whole node.
        return claimable, enclosing

    fh = io.StringIO(content)
    while True:
        raw_line, processed_line, _ = get_fortran_line(fh)
        if not raw_line and not processed_line:
            break
        enclosing[index] = list(regions)

        if processed_line.startswith('#'):
            kind, _name, _invert = _pp_directive(processed_line)
            if kind in ('open', 'openOther'):
                regions.append({'lines': [index], 'claimable': kind == 'open'})
            elif kind in ('else', 'elseOther'):
                if regions:
                    regions[-1]['lines'].append(index)
                    if kind == 'elseOther':
                        regions[-1]['claimable'] = False
                else:
                    claimable[index] = False
            elif kind in ('close', 'closeOther'):
                if regions:
                    region = regions.pop()
                    region['lines'].append(index)
                    verdict = region['claimable'] and kind == 'close'
                    for i in region['lines']:
                        claimable[i] = verdict
                    if not verdict:
                        # The region's own directives become ordinary code
                        # inside every region enclosing it.
                        for outer in regions:
                            outer['claimable'] = False
                else:
                    claimable[index] = False
            else:
                # A directive we cannot rebuild — leave it in place as code.
                claimable[index] = False
                for outer in regions:
                    outer['claimable'] = False
        elif not _MODULE_USE_RE.match(processed_line):
            # Plain code (or a comment, or a blank line) — it will close any
            # open moduleUse block, so no enclosing guard can be claimed.
            for outer in regions:
                outer['claimable'] = False

        index += 1

    # Regions left open at the end of the node never close within it.
    for region in regions:
        region['claimable'] = False
        for i in region['lines']:
            claimable[i] = False

    for index, open_regions in enclosing.items():
        enclosing[index] = [r for r in open_regions if not r['claimable']]

    return claimable, enclosing


def _as_entry_list(value):
    """Normalize a `moduleUse[module]` value to a list of entry dicts.

    Parsed `moduleUse` nodes store a *list* of entries per module — one per
    distinct preprocessor condition set — so that the same module `use`d
    under different `#ifdef` guards within one unit is preserved and each
    occurrence re-emitted under its own guard (issue #1030).  Callers of
    `add_uses` may still pass a bare entry dict for the common unconditional
    case; treat that as a single-element list.
    """
    if isinstance(value, list):
        return value
    return [value]


def _condition_key(entry):
    """A hashable identity for an entry's preprocessor condition set."""
    return tuple((c['name'], c['invert']) for c in entry.get('conditions', []))


def _insert_module_ordered(module_order, module_uses, mod_name, is_intrinsic):
    """Record `mod_name` in `module_order`, intrinsic modules first.

    Intrinsic modules (`use, intrinsic :: ISO_C_Binding`) are hoisted ahead of
    the rest, but keep their order of appearance *among themselves* — they are
    inserted after any intrinsic modules already recorded, not at index zero.

    Inserting at index zero instead would reverse the intrinsic modules on
    every pass, so a block naming two of them oscillated between two layouts
    and reformatting never reached a fixed point.
    """
    if not is_intrinsic:
        module_order.append(mod_name)
        return
    index = 0
    for name in module_order:
        if not any(entry.get('intrinsic')
                   for entry in _as_entry_list(module_uses.get(name, []))):
            break
        index += 1
    module_order.insert(index, mod_name)


def parse_module_uses(tree):
    """Walk the tree parsing Fortran 'use' statements into moduleUse nodes.

    Called by _pass_module_uses() in Galacticus.Build.SourceTree.
    """
    from Galacticus.Build.SourceTree import walk_tree, replace_node, _make_code_node

    nodes_to_replace = []

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        source  = node.get('source', 'unknown')
        line_no = node.get('line', 0)

        new_nodes       = []
        raw_code_buf    = []
        raw_use_buf     = []
        raw_pp_buf      = []   # preprocessor lines between/before use lines
        module_uses     = {}
        module_order    = []
        pp_stack        = []   # list of {'name': str, 'invert': bool}

        raw_code_line   = line_no
        raw_use_line    = line_no
        pp_first_line   = line_no   # first line of the buffered pp-directive run
        current_line    = line_no
        line_index      = 0         # logical-line index within this node
        use_line_index  = 0         # …of the first line of the open use block

        # Which preprocessor directives in this node may be absorbed into a
        # moduleUse block and re-emitted from its `conditions` (issue #1385).
        claimable, enclosing = _classify_preprocessor_directives(content)

        fh = io.StringIO(content)

        def _flush_code():
            nonlocal raw_code_buf, raw_code_line, raw_use_line
            # NB: deliberately leaves `raw_pp_buf` alone — preprocessor lines
            # immediately preceding a `use` belong with the moduleUse block
            # (they record the directive's condition in `entry['conditions']`
            # and `_build_module_uses` re-emits a matching `#ifdef … #endif`
            # wrapper).  Absorbing them into the code node here would cause
            # the wrapper to be emitted twice — once by the surrounding code
            # node, once by the rebuild — leaving the output with one extra
            # `#ifdef` and an unmatched preprocessor block.  The plain-code
            # `else` branch and the EOF handler still absorb any pp lines
            # that were *not* followed by a `use`.
            if raw_code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(raw_code_buf), source, raw_code_line))
                raw_code_buf.clear()
                raw_code_line = current_line
                raw_use_line  = current_line

        def _flush_uses():
            nonlocal module_uses, module_order, raw_use_buf, raw_code_line, raw_use_line
            if not module_uses:
                return
            mu_node = {
                'type':        'moduleUse',
                'moduleUse':   module_uses,
                'moduleOrder': module_order,
                'parent':      None,
                'firstChild':  None,
                'sibling':     None,
                'source':      source,
                'line':        raw_use_line,
            }
            # Record whether the block sits inside a preprocessor guard that
            # is left in place as ordinary code.  `add_uses` avoids merging
            # into such a block: an unconditional import added there would be
            # compiled only when that guard happens to be satisfied.
            if enclosing.get(use_line_index):
                mu_node['inUnclaimedGuard'] = True
            child = _make_code_node(''.join(raw_use_buf), source, raw_use_line)
            child['parent'] = mu_node
            mu_node['firstChild'] = child
            new_nodes.append(mu_node)
            module_uses   = {}
            module_order  = []
            raw_use_buf   = []
            raw_code_line = current_line
            raw_use_line  = current_line

        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                # EOF — flush whatever is pending
                if raw_pp_buf:
                    if not raw_code_buf:
                        raw_code_line = pp_first_line
                    raw_code_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                _flush_code()
                _flush_uses()
                break

            n_newlines  = raw_line.count('\n')
            line_after  = current_line + n_newlines

            m = _MODULE_USE_RE.match(processed_line)
            is_use = bool(m)

            if is_use:
                # Flush any preceding plain code first.
                _flush_code()
                # A new use block starts on the first absorbed preprocessor
                # line (if any) — the block's emitted content begins there,
                # so its `.lmap` anchor must too.
                if not raw_use_buf:
                    raw_use_line   = pp_first_line if raw_pp_buf else current_line
                    use_line_index = line_index
                # Absorb any buffered preprocessor lines into the use block.
                if raw_pp_buf:
                    raw_use_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                raw_use_buf.append(raw_line)

                is_openmp    = bool(m.group(1))
                is_intrinsic = bool(m.group(3))
                mod_name     = m.group(5)
                only_text    = m.group(7)

                # Entries are keyed by module name *and* preprocessor
                # condition set: the same module may be `use`d under
                # different `#ifdef` guards within one unit, and each such
                # occurrence must be preserved and re-emitted under its own
                # guard (issue #1030).
                conditions = copy.deepcopy(pp_stack) if pp_stack else []
                cond_key   = tuple((c['name'], c['invert']) for c in conditions)

                entries = module_uses.setdefault(mod_name, [])
                if not entries:
                    # First occurrence of this module name fixes its order.
                    _insert_module_ordered(
                        module_order, module_uses, mod_name, is_intrinsic)

                entry = next(
                    (e for e in entries if _condition_key(e) == cond_key), None)
                if entry is None:
                    entry = {'openMP': is_openmp, 'intrinsic': is_intrinsic}
                    if conditions:
                        entry['conditions'] = conditions
                    entries.append(entry)
                else:
                    entry['openMP']    = is_openmp
                    entry['intrinsic'] = is_intrinsic

                if only_text and 'all' not in entry:
                    only_text = only_text.strip()
                    only_dict = entry.setdefault('only', {})
                    for sym in re.split(r'\s*,\s*', only_text):
                        sym = re.sub(r'\s', '', sym)
                        if sym:
                            only_dict[sym] = True
                else:
                    entry['all'] = True
                    entry.pop('only', None)

            elif processed_line.startswith('#') and claimable.get(line_index):
                # Preprocessor directive guarding nothing but `use` statements
                # — absorb it into the use block and record it as a condition,
                # from which `update_uses` rebuilds the wrapper.  Directives
                # the pre-scan did not clear fall through to the plain-code
                # branch below and are emitted where they stand (issue #1385).
                if module_uses:
                    # Inside a use block — keep it in the use buffer.
                    is_use = True
                    raw_use_buf.append(raw_line)
                else:
                    if not raw_pp_buf:
                        pp_first_line = current_line
                    raw_pp_buf.append(raw_line)

                # Update the preprocessor stack.
                kind, name, invert = _pp_directive(processed_line)
                if kind == 'open':
                    pp_stack.append({'name': name, 'invert': invert})
                elif kind == 'close' and pp_stack:
                    pp_stack.pop()
                elif kind == 'else' and pp_stack:
                    top = pp_stack.pop()
                    top['invert'] = not top['invert']
                    pp_stack.append(top)

            else:
                # Plain code line — close any open use block first.
                _flush_uses()
                if not raw_code_buf:
                    raw_code_line = pp_first_line if raw_pp_buf else current_line
                if raw_pp_buf:
                    raw_code_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                raw_code_buf.append(raw_line)

            current_line = line_after
            line_index  += 1

        # Only replace if the node actually changed.
        single_code = (
            len(new_nodes) == 1
            and new_nodes[0].get('type') == 'code'
            and new_nodes[0].get('content') == content
        )
        if not single_code and new_nodes:
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)


def update_uses(uses_node):
    """Regenerate formatted Fortran 'use' statements from structured node data.

    Rewrites uses_node['firstChild']['content'] in place.
    """
    from Galacticus.Build.FortranUtils import get_fortran_line

    # Determine the block's indentation from the existing raw content: the
    # smallest indent of any `use` statement in it.
    #
    # The minimum is what keeps this function idempotent.  The emitter below
    # starts every statement at the block's indent, so re-reading its own
    # output recovers that indent exactly — but a block laid out by hand (or
    # by an earlier convention, which padded unconditional `use` statements
    # three columns to the right so that `use` lined up with the `use` after
    # a `!$ ` sentinel) may indent some of its statements further.  Taking the
    # minimum recovers the base indent in either case.
    #
    # Only `use` lines are considered: any `#ifdef` guards in the block start
    # at column zero and would otherwise drag the indent to nothing.
    raw_content = uses_node['firstChild']['content'] if uses_node.get('firstChild') else ""
    indents = []
    fh = io.StringIO(raw_content)
    while True:
        raw_line, processed_line, _ = get_fortran_line(fh)
        if not raw_line and not processed_line:
            break
        if _MODULE_USE_RE.match(processed_line):
            indents.append(re.match(r'^(\s*)', raw_line).group(1))
    indent = min(indents, key=len) if indents else ""

    module_use  = uses_node.get('moduleUse', {})
    module_order = uses_node.get('moduleOrder', list(module_use.keys()))

    openmp_any = any(
        e.get('openMP')
        for v in module_use.values() for e in _as_entry_list(v))
    intrinsic_any = any(
        e.get('intrinsic')
        for v in module_use.values() for e in _as_entry_list(v))

    name_len_max = max((len(n) for n in module_use), default=0)

    def _column_widths(symbols_per_row):
        """Max width of each symbol column, over every `use` in the block.

        Widths are shared across the whole block so that symbols line up
        between `use` statements, not merely within one.
        """
        col_max = [0] * symbols_per_row
        for v in module_use.values():
            for entry in _as_entry_list(v):
                for i, sym in enumerate(sorted(entry.get('only', {}).keys())):
                    j = i % symbols_per_row
                    if len(sym) > col_max[j]:
                        col_max[j] = len(sym)
        return col_max

    def _emit_entry(mod_name, entry, symbols_per_row, col_max):
        text = ""

        # Emit preprocessor conditions.
        for cond in entry.get('conditions', []):
            directive = "#ifndef" if cond['invert'] else "#ifdef"
            text += directive + " " + cond['name'] + "\n"

        use_line = indent

        # The `!$ ` sentinel is a prefix, so an OpenMP-conditional `use`
        # starts three columns to the right of an unconditional one.  The
        # padding that brings the `::` back into line goes *after* the module
        # attributes rather than before `use`: every statement then begins at
        # the block's own indent, with `use` itself aligned down the block
        # and the sentinel sitting out to its left where it reads as the
        # annotation it is.
        if entry.get('openMP'):
            use_line += "!$ "

        use_line += "use"

        if entry.get('intrinsic'):
            use_line += ", intrinsic"
        elif intrinsic_any:
            use_line += "           "

        if openmp_any and not entry.get('openMP'):
            use_line += "   "

        use_line += " :: " + mod_name

        only = entry.get('only', {})
        symbol_count = len(only)
        if symbol_count > 0:
            use_line += " " * (name_len_max - len(mod_name)) + ", only : "
            offset_len = len(use_line)
            symbols = sorted(only.keys())
            remaining = symbol_count
            for i, sym in enumerate(symbols):
                j = i % symbols_per_row
                use_line += sym
                remaining -= 1
                if remaining > 0:
                    use_line += " " * (col_max[j] - len(sym))
                    use_line += ", "
                    if j == symbols_per_row - 1:
                        use_line += "&\n"
                        # Continuations of an OpenMP-conditional `use` need
                        # the sentinel too: without it the continuation is
                        # compiled in serial builds, where the statement it
                        # continues is only a comment.
                        if entry.get('openMP'):
                            cont_prefix = indent + "!$ "
                        else:
                            cont_prefix = indent
                        use_line += cont_prefix + "&" + " " * (offset_len - len(cont_prefix) - 1)

        use_line += "\n"
        text += use_line

        # Close preprocessor conditions.
        for _ in entry.get('conditions', []):
            text += "#endif\n"
        return text

    def _emit_block(symbols_per_row):
        col_max = _column_widths(symbols_per_row)
        text = ""
        for mod_name in module_order:
            for entry in _as_entry_list(module_use.get(mod_name, [])):
                text += _emit_entry(mod_name, entry, symbols_per_row, col_max)
        return text

    # Use as many symbols per row as fit: SYMBOLS_PER_ROW_MAX normally, fewer
    # when the resulting rows would run past LINE_LENGTH_MAX.  The count is
    # chosen once for the whole block rather than per statement, because the
    # column widths are shared across the block — varying it per statement
    # would destroy the alignment the widths exist to produce.
    #
    # A block whose symbols are too long to fit even one per row still emits
    # one per row and overflows: there is nowhere to break a single symbol,
    # and Galacticus compiles with `-ffree-line-length-none`, so the limit is
    # a style rule rather than a constraint.
    for symbols_per_row in range(SYMBOLS_PER_ROW_MAX, 0, -1):
        content = _emit_block(symbols_per_row)
        if max((len(line) for line in content.splitlines()),
               default=0) <= LINE_LENGTH_MAX:
            break

    if uses_node.get('firstChild') is None:
        uses_node['firstChild'] = {
            'type':       'code',
            'content':    content,
            'parent':     uses_node,
            'sibling':    None,
            'firstChild': None,
            'source':     uses_node.get('source', 'unknown'),
            'line':       uses_node.get('line', 0),
        }
    else:
        uses_node['firstChild']['content'] = content


def add_uses(node, module_uses_node):
    """Add module uses to an existing node, creating a moduleUse child if needed.

    Parameters
    ----------
    node : dict
        An AST node whose children will receive the new module uses.
    module_uses_node : dict
        A moduleUse-structured dict with 'moduleUse' and 'moduleOrder' keys.
    """
    from Galacticus.Build.SourceTree import insert_before_node

    # Find an existing moduleUse child.  Blocks sitting inside a preprocessor
    # guard that the parser left in place are skipped: adding an import there
    # would compile it only when that guard is satisfied.  If every block is
    # guarded, a fresh, unguarded one is created ahead of them all.
    uses_node = None
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'moduleUse' and not child.get('inUnclaimedGuard'):
            uses_node = child
        child = child.get('sibling')

    if uses_node is None:
        # Create a new, empty moduleUse node.
        uses_node = {
            'type':        'moduleUse',
            'moduleUse':   {},
            'moduleOrder': [],
            'parent':      None,
            'sibling':     None,
            'source':      'Galacticus.Build.SourceTree.Parse.ModuleUses.add_uses()',
            'line':        1,
        }
        uses_node['firstChild'] = {
            'type':       'code',
            'content':    '',
            'parent':     uses_node,
            'sibling':    None,
            'firstChild': None,
            'source':     uses_node['source'],
            'line':       1,
        }
        first_child = node.get('firstChild')
        if first_child is None:
            uses_node['parent'] = node
            uses_node['sibling'] = None
            node['firstChild'] = uses_node
        else:
            insert_before_node(first_child, [uses_node])

    # Merge module uses.
    new_module_use   = module_uses_node.get('moduleUse', {})
    new_module_order = module_uses_node.get('moduleOrder', sorted(new_module_use.keys()))

    for mod_name in new_module_order:
        for new_entry in _as_entry_list(new_module_use[mod_name]):
            cond_key = _condition_key(new_entry)
            entries  = uses_node['moduleUse'].setdefault(mod_name, [])

            if mod_name not in uses_node['moduleOrder']:
                _insert_module_ordered(
                    uses_node['moduleOrder'], uses_node['moduleUse'],
                    mod_name, new_entry.get('intrinsic'))

            # Merge into the entry with the *same* condition set.  An
            # unconditional `use` thus joins (or creates) the unconditional
            # entry and leaves any `#ifdef`-guarded entry of the same module
            # intact — and vice versa.  This is what fixes the condition
            # conflation of issue #1030: with a single per-module entry, a
            # later occurrence's guard used to overwrite (or be merged
            # under) the earlier one's, moving symbols into the wrong block.
            existing = next(
                (e for e in entries if _condition_key(e) == cond_key), None)
            if existing is None:
                existing = {}
                if 'conditions' in new_entry:
                    existing['conditions'] = copy.deepcopy(new_entry['conditions'])
                entries.append(existing)
            existing['openMP']    = new_entry.get('openMP', False)
            existing['intrinsic'] = new_entry.get('intrinsic', False)

            if 'all' not in existing:
                if new_entry.get('all'):
                    existing['all'] = True
                    existing.pop('only', None)
                else:
                    only = existing.setdefault('only', {})
                    for sym in new_entry.get('only', {}):
                        only[sym] = True

    update_uses(uses_node)
