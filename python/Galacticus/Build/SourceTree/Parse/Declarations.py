"""Provides parsing of Fortran variable declaration lines.

Andrew Benson (ported to Python 2026)
"""

import re
import copy

from Fortran.Utils import INTRINSIC_DECLARATIONS, extract_variables
from Galacticus.Build.FortranUtils import extract_bracketed


def parse_declaration(line):
    """Parse a single Fortran variable declaration line.

    Parameters
    ----------
    line : str
        A processed (continuation-joined, comment-stripped) Fortran line.

    Returns
    -------
    dict or None
        None if the line is not a recognised variable declaration.
        Otherwise a dict with keys:
          'intrinsic'     : str   — e.g. 'integer', 'real', 'double precision', …
          'type'          : str or None — kind/type spec (parentheses stripped)
          'openMP'        : bool
          'attributes'    : list of str — e.g. ['intent(in)', 'allocatable']
          'attributesOriginal' : list of str — as above, original case
          'variables'     : list of str — lowercase names, qualifiers preserved
          'variableNames' : list of str — original-case names, no qualifiers
          'variablesOriginal' : list of str — original-case, qualifiers preserved
    """
    # Use a unified approach: match declaration pattern, then split manually
    # This avoids complex capture group indexing issues.
    # Match: [!$] intrinsic-type [attributes/type-spec] :: variables
    decl_pattern = re.compile(
        r'^\s*(!\$)?\s*'  # Optional OpenMP marker
        r'(integer|real|double\s+precision|double\s+complex|logical|character|type|class|procedure|generic|final|complex)'  # intrinsic
        r'(.*)$',  # Everything else
        re.IGNORECASE
    )
    m = decl_pattern.match(line)
    if not m or '::' not in line:
        return None

    # Now split on the :: delimiter
    before_decl, variables_raw = line.split('::', 1)
    if not variables_raw.strip():
        return None

    openmp_marker = m.group(1)
    intrinsic_raw = m.group(2).lower()
    # Normalize whitespace in intrinsic (e.g. double   precision -> double precision)
    intrinsic_raw = re.sub(r'\s+', ' ', intrinsic_raw)

    # Extract type/attributes from before :: delimiter
    type_attrs_raw = before_decl.split(m.group(2), 1)[1].strip() if m.group(2) else ''
    variables_raw = variables_raw.strip()

    # Normalize intrinsic name
    intrinsic_map = {
        'double precision': 'double precision',
        'double complex': 'double complex',
    }
    intrinsic = intrinsic_map.get(intrinsic_raw, intrinsic_raw)

    # Parse type and attributes
    type_val = None
    attributes = []
    attributes_original = []
    if type_attrs_raw:
        # Carve a leading balanced parenthesised group off `type_attrs_raw` if
        # one is present.  Uses `extract_bracketed` so that nested parens
        # — e.g. `character(len=len(tagName))` — are matched correctly; the
        # earlier `\(\s*([^)]+)\s*\)` regex stopped at the first `)` and
        # produced `len=len(tagName` as the type, dropping the outer paren
        # and breaking the round-trip.
        leading_parens = None
        rest_after_parens = type_attrs_raw.lstrip()
        if rest_after_parens.startswith('('):
            extracted, remainder, _ = extract_bracketed(
                rest_after_parens, brackets="()")
            if extracted is not None:
                leading_parens   = extracted          # includes the outer ()
                rest_after_parens = remainder.lstrip()

        # The remainder (after the optional leading parens) is the attributes
        # list.  It still uses a single comma-separated form, so split on the
        # leading comma if present.
        first_part = leading_parens or ''
        rest       = rest_after_parens.lstrip(', ').strip() \
            if rest_after_parens.startswith(',') else rest_after_parens.strip()

        # Decide whether the leading parens should be consumed as the
        # type-spec.  If `first_part` *starts* with `(` (i.e. the parens are
        # the very first non-space content after the intrinsic), they are
        # always the type-spec — `integer(c_size_t) :: …`,
        # `type(varying_string) :: …`, `procedure(template), nopass :: …`,
        # `character(len=len(tag)) :: …`, etc.
        consumed_as_type = False
        if first_part.startswith('(') and first_part.endswith(')'):
            # Strip the outer parens to get the type-spec body.
            type_val = first_part[1:-1].strip()
            consumed_as_type = True

        # Attributes are extracted twice for the same reason as the variables
        # below: the lower-cased form is what the code generators match
        # against, but a formatter must write `dimension(nBins)` back with its
        # original case rather than `dimension(nbins)`.
        if rest:
            attributes          = extract_variables(rest, keep_qualifiers=True)
            attributes_original = extract_variables(rest, keep_qualifiers=True,
                                                    lower_case=False)
        elif not consumed_as_type and first_part:
            # `first_part` carries parens that aren't a type-spec (rare —
            # `dimension(:)` written without a leading comma).
            attributes          = extract_variables(first_part, keep_qualifiers=True)
            attributes_original = extract_variables(first_part, keep_qualifiers=True,
                                                    lower_case=False)
        else:
            attributes          = []
            attributes_original = []

    variables = extract_variables(variables_raw, keep_qualifiers=True, lower_case=True)
    variable_names = extract_variables(variables_raw, keep_qualifiers=False, lower_case=False)
    # `variables` is lower-cased (the code generators match against it
    # case-insensitively) and `variableNames` drops qualifiers.  A formatter
    # needs both preserved at once — it must write `massStellar` back as
    # `massStellar`, initializer and dimension qualifiers included — so keep a
    # third rendering rather than perturbing either existing one.
    variables_original = extract_variables(
        variables_raw, keep_qualifiers=True, lower_case=False)

    return {
        'intrinsic':          intrinsic,
        'type':               type_val,
        'openMP':             bool(openmp_marker),
        'attributes':         attributes,
        'attributesOriginal': attributes_original,
        'variables':          variables,
        'variableNames':      variable_names,
        'variablesOriginal':  variables_original,
    }


def build_declarations(node):
    """Build Fortran declaration text from a declaration node's structured data.

    Rewrites node['firstChild']['content'] in place.
    """
    content = "implicit none\n" if node.get('implicitNone') else ""
    for declaration in node.get('declarations', []):
        line = ""
        if declaration.get('preprocessor'):
            line += "#ifdef " + declaration['preprocessor'] + "\n"
        line += "  "
        if declaration.get('openMP'):
            line += "!$ "
        line += declaration['intrinsic']
        type_val = declaration.get('type')
        if type_val is not None:
            has_parens = type_val.startswith('(') and type_val.endswith(')')
            if has_parens:
                line += type_val
            else:
                line += "(" + type_val + ")"
        attrs = declaration.get('attributes') or []
        if attrs:
            line += ", " + ", ".join(attrs)
        line += " :: " + ", ".join(declaration.get('variables', [])) + "\n"
        if declaration.get('threadprivate'):
            names = [re.sub(r'([a-zA-Z0-9_]+).*', r'\1', v)
                     for v in declaration.get('variables', [])]
            line += " !$omp threadprivate(" + ",".join(names) + ")\n"
        if declaration.get('preprocessor'):
            line += "#endif\n"
        content += line
    if node.get('firstChild') is None:
        node['firstChild'] = {
            'type':       'code',
            'content':    content,
            'parent':     node,
            'sibling':    None,
            'firstChild': None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }
    else:
        node['firstChild']['content'] = content


def add_declarations(node, declarations):
    """Add declarations to an existing node, creating a declaration child if needed.

    The new declaration node is inserted after any existing moduleUse child,
    otherwise before the first child.
    """
    from Galacticus.Build.SourceTree import insert_before_node, insert_after_node

    declarations_node = None
    uses_node = None
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'declaration' and declarations_node is None:
            declarations_node = child
        if child.get('type') == 'moduleUse' and uses_node is None:
            uses_node = child
        child = child.get('sibling')

    if declarations_node is None:
        declarations_node = {
            'type':         'declaration',
            'declarations': [],
            'parent':       None,
            'firstChild':   None,
            'sibling':      None,
            'source':       node.get('source', 'unknown'),
            'line':         node.get('line', 0),
        }
        declarations_node['firstChild'] = {
            'type':       'code',
            'content':    '',
            'parent':     declarations_node,
            'sibling':    None,
            'firstChild': None,
            'source':     declarations_node['source'],
            'line':       declarations_node['line'],
        }
        if uses_node is not None:
            insert_after_node(uses_node, [declarations_node])
        else:
            first_child = node.get('firstChild')
            if first_child is None:
                declarations_node['parent'] = node
                node['firstChild'] = declarations_node
            else:
                insert_before_node(first_child, [declarations_node])

    declarations_node.setdefault('declarations', []).extend(declarations)
    build_declarations(declarations_node)


def add_attributes(node, variable_name, attributes):
    """Add attributes to the declaration of a named variable.

    If the variable shares a declaration line with others, that declaration is
    split so the target variable can receive attributes independently.
    """
    declarations_found = None
    declaration_found  = None
    child = node.get('firstChild')
    target = variable_name.lower()
    while child:
        if child.get('type') == 'declaration':
            declarations_found = child
            for declaration in child.get('declarations', []):
                if target in declaration.get('variables', []):
                    declaration_found = declaration
                    break
            if declaration_found is not None:
                break
        child = child.get('sibling')

    if declarations_found is None:
        raise RuntimeError(
            f'add_attributes: no declarations present in '
            f'{node.get("type")} "{node.get("name")}"')
    if declaration_found is None:
        raise RuntimeError(
            f'add_attributes: variable declaration [{variable_name}] not found in '
            f'{node.get("type")} "{node.get("name")}"')

    variables = declaration_found.get('variables', [])
    if len(variables) > 1:
        declaration_copy = copy.deepcopy(declaration_found)
        index = declaration_copy['variables'].index(target)
        declaration_found['variables'    ] = [declaration_copy['variables'    ][index]]
        declaration_found['variableNames'] = [declaration_copy['variableNames'][index]]
        declaration_copy['variables'    ].pop(index)
        declaration_copy['variableNames'].pop(index)
        declarations_found.setdefault('declarations', []).append(declaration_copy)

    declaration_found.setdefault('attributes', [])
    for attribute in attributes:
        if attribute not in declaration_found['attributes']:
            declaration_found['attributes'].append(attribute)

    build_declarations(declarations_found)


def get_declaration(node, variable_name):
    """Return a descriptor of the declaration for a named variable.

    The returned dict is a deep copy with 'variables' reduced to the single
    target variable (original case).
    Raises RuntimeError if no declarations are present or the variable is not found.
    """
    declarations_found = False
    declaration_found  = None
    target = variable_name.lower()
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'declaration':
            declarations_found = True
            for declaration in child.get('declarations', []):
                # Strip any initializer value before comparing.
                bare_names = []
                for v in declaration.get('variables', []):
                    m_eq = re.match(r'^([^=]+)\s*=', v)
                    bare_names.append(m_eq.group(1) if m_eq else v)
                if target in bare_names:
                    declaration_found = copy.deepcopy(declaration)
                    declaration_found['variables'] = [variable_name]
                    break
            if declaration_found is not None:
                break
        child = child.get('sibling')

    if not declarations_found:
        raise RuntimeError(
            'get_declaration: no declarations present')
    if declaration_found is None:
        raise RuntimeError(
            f'get_declaration: variable declaration for "{variable_name}" not found '
            f'in node "{node.get("opener")}"')
    return declaration_found


def declaration_exists(node, variable_name):
    """Return True if the named variable has a declaration in node.

    Case-insensitive.  Strips any `=…` initializer off each stored variable name before
    comparing — declarations are stored as raw `name=value` tokens (e.g.
    `warnObjectBuilder0__=.false.`) but callers query by bare name.  Also
    falls back to `variableNames` (the parser's already-stripped list)
    when present.
    """
    target = variable_name.lower()
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'declaration':
            for declaration in child.get('declarations', []):
                for v in declaration.get('variableNames') or []:
                    if v.lower() == target:
                        return True
                for v in declaration.get('variables', []):
                    bare = re.match(r'^([^=\s]+)', v)
                    if bare and bare.group(1).lower() == target:
                        return True
        child = child.get('sibling')
    return False


# ---------------------------------------------------------------------------
# Formatting
#
# `build_declarations()` above serializes declarations for the *code
# generators*: fixed indent, no alignment, minimal spacing.  The functions
# below instead reproduce the hand-maintained house style described in the
# developer guide, for use by `scripts/aux/formatDeclarations.py`.  They are
# deliberately separate so that changing the formatting rules cannot perturb
# generated code.
# ---------------------------------------------------------------------------

# Line width beyond which the emitter packs fewer variables per row, and the
# usual maximum number of variables on one row.
LINE_LENGTH_MAX      = 132
VARIABLES_PER_ROW_MAX = 2

# Canonical ordering of declaration attributes.
#
# Fortran attaches no meaning to attribute order, but fixing one turns "the
# same attribute should sit in the same column" from a multiple-sequence
# alignment problem into a per-position maximum width.  This order was derived
# from pairwise precedence counts over the whole source tree — for every pair
# of attributes, which one appears first more often — and reorders about 17% of
# the multi-attribute declarations in the tree (roughly 3% of all
# declarations).  Attributes not listed sort last, alphabetically.
ATTRIBUTE_ORDER = (
    'intent',
    'pointer',
    'allocatable',
    'dimension',
    'target',
    'optional',
    'parameter',
    'save',
    'value',
    'public',
    'private',
)

_ATTRIBUTE_RANK = {name: index for index, name in enumerate(ATTRIBUTE_ORDER)}

# Type-bound declarations (`procedure :: name => target`) carry no attribute or
# variable columns and are laid out on their own axis.
_TYPE_BOUND_INTRINSICS = frozenset(('procedure', 'generic', 'final'))

# The `intent` field is padded to the width of the longest spelling, `inout`,
# so that the three forms occupy one column: `in` is left-aligned within it and
# `out` right-aligned, giving `intent(in   )`, `intent(inout)`, `intent(  out)`.
_INTENT_FIELD = {'in': 'in   ', 'out': '  out', 'inout': 'inout'}


def _attribute_name(attribute):
    """The bare name of an attribute, without any parenthesized argument."""
    return re.sub(r'\(.*', '', attribute).strip().lower()


def _normalize_attribute(attribute):
    """Return `attribute` in canonical spelling."""
    name = _attribute_name(attribute)
    if name == 'intent':
        match = re.match(r'^\s*intent\s*\(\s*([a-zA-Z]*)\s*\)\s*$', attribute,
                         re.IGNORECASE)
        if match:
            field = _INTENT_FIELD.get(match.group(1).lower())
            if field is not None:
                return 'intent(' + field + ')'
    return attribute


def _ordered_attributes(declaration):
    """Attributes of `declaration`, canonically ordered and spelled.

    The sort is stable on the canonical rank, so attributes sharing a rank
    (only those outside `ATTRIBUTE_ORDER`, which all rank last) keep their
    source order relative to one another apart from the alphabetical tie-break.
    """
    source = (declaration.get('attributesOriginal')
              or declaration.get('attributes')
              or [])
    attributes = [_normalize_attribute(a) for a in source]
    return sorted(
        attributes,
        key=lambda a: (_ATTRIBUTE_RANK.get(_attribute_name(a), len(ATTRIBUTE_ORDER)),
                       _attribute_name(a)),
    )


def _split_initializer(variable):
    """Split `variable` into (name, operator, value).

    `operator` is `'=>'`, `'='` or `''`.  The scan ignores `=` characters
    inside brackets, so that a dimension qualifier such as `x(kind=c_int)` is
    not mistaken for an initializer.
    """
    depth = 0
    for index, character in enumerate(variable):
        if character in '([':
            depth += 1
        elif character in ')]':
            depth -= 1
        elif depth == 0 and character == '=':
            if variable[index + 1:index + 2] == '>':
                return variable[:index], '=>', variable[index + 2:]
            return variable[:index], '=', variable[index + 1:]
    return variable, '', ''


def _paren_body(declaration):
    """The body of a declaration's parenthetical, without its brackets.

    `parse_declaration` strips the outer brackets from most type-specs but
    leaves them on some, so normalize both spellings here.
    """
    type_value = declaration.get('type')
    if type_value is None:
        return ''
    if type_value.startswith('(') and type_value.endswith(')'):
        return type_value[1:-1].strip()
    return type_value.strip()


def _variables_of(declaration):
    """The declaration's variables, original case and qualifiers preserved.

    Falls back to the lower-cased `variables` for declaration dicts built by
    the code generators, which do not carry `variablesOriginal`.
    """
    return (declaration.get('variablesOriginal')
            or declaration.get('variables')
            or [])


def _operator_text(operator):
    """The operator as written, including any spaces around it.

    `=>` is spaced on both sides (`x => null()`) while `=` is written tight
    against both name and value (`x=0.0d0`).  The spaces belong to the operator
    rather than to the name column, so that the widest name in a column keeps
    its separator too — `nameThatIsLongest => target`, not `nameThatIsLongest=>`.
    """
    return ' => ' if operator == '=>' else operator


def _variable_fields(declaration, variable):
    """Split one variable into its (name, operator, value) columns.

    A `final` binding names no target, and the house style places it in the
    *target* column rather than the name column so that it lines up with the
    targets of the `procedure` bindings beside it — `final :: fooDestructor`
    sits under `procedure :: bar => fooBar`'s `fooBar`.
    """
    if declaration['intrinsic'] == 'final':
        return '', '', variable
    return _split_initializer(variable)


def _variables_per_row(declaration, per_row):
    """How many variables share one row for this declaration.

    Type-bound declarations (`procedure :: name => target`) take one binding
    per row regardless of the group's setting; everything else takes `per_row`.
    """
    if declaration['intrinsic'] in _TYPE_BOUND_INTRINSICS:
        return 1
    return per_row


class _Layout:
    """Column widths shared by every declaration in one alignment group.

    Attribute columns are *positional*, not per-attribute: column 0 holds
    whichever attribute comes first on each row, so a row with `pointer` and a
    row with `intent(inout)` share one column sized to the wider of the two.
    Canonical ordering (`ATTRIBUTE_ORDER`) is what makes that produce a
    coherent table — without it, position 0 would hold an arbitrary attribute
    on each row.
    """

    def __init__(self, declarations, per_row):
        self.intrinsic = max((len(d['intrinsic']) for d in declarations), default=0)

        parenthesized = [d for d in declarations if d.get('type') is not None]
        self.paren = max((len(_paren_body(d)) for d in parenthesized), default=0)
        # Width of the whole parenthetical field, brackets included.  Zero when
        # no declaration in the group has one, so that a group of plain
        # `double precision` declarations carries no dead column.
        self.paren_field = self.paren + 2 if parenthesized else 0

        attribute_lists = [_ordered_attributes(d) for d in declarations]
        self.attributes = [
            max(a[column] and len(a[column]) or 0
                for a in attribute_lists if len(a) > column)
            for column in range(max((len(a) for a in attribute_lists), default=0))
        ]

        # Variable columns.  Name, operator and value are sized independently
        # so that `=` and `=>` line up down each column.
        self.names     = []
        self.operators = []
        self.values    = []
        for declaration in declarations:
            columns = _variables_per_row(declaration, per_row)
            for index, variable in enumerate(_variables_of(declaration)):
                column = index % columns
                while len(self.names) <= column:
                    self.names.append(0)
                    self.operators.append(0)
                    self.values.append(0)
                name, operator, value = _variable_fields(declaration, variable)
                self.names[column]     = max(self.names[column], len(name))
                self.operators[column] = max(self.operators[column],
                                             len(_operator_text(operator)))
                self.values[column]    = max(self.values[column], len(value))

    def prefix_width(self, indent, openmp_any):
        """The column at which the variable list starts."""
        width = len(indent) + (3 if openmp_any else 0)
        width += self.intrinsic + self.paren_field
        for column_width in self.attributes:
            width += 2 + column_width
        return width + 4          # " :: "


def _emit_variable(declaration, variable, layout, column, pad):
    """Render one variable into its column.

    `pad` is False for the last variable on the last row, so that the line ends
    ragged rather than carrying trailing whitespace.  The *name* is still
    padded whenever the column holds any initializer, since otherwise the `=`
    and `=>` of a column would not line up.
    """
    name, operator, value = _variable_fields(declaration, variable)
    # Anything after the name — an operator column shared with other rows, or
    # this row's own value — means the name must be padded out to its column so
    # that what follows lines up.  A `final` binding has *only* a value, so the
    # value must be emitted whether or not the column carries an operator.
    trailing = bool(layout.operators[column]) or bool(value)
    text = name.ljust(layout.names[column]) if (pad or trailing) else name
    if trailing:
        text += _operator_text(operator).ljust(layout.operators[column])
        text += value.ljust(layout.values[column]) if pad else value
    return text


def _emit_declaration(declaration, layout, indent, openmp_any, per_row):
    """Render one declaration statement, with continuation lines as needed."""
    prefix = indent
    if declaration.get('openMP'):
        prefix += '!$ '
    elif openmp_any:
        prefix += '   '

    head = prefix + declaration['intrinsic'].ljust(layout.intrinsic)

    if layout.paren_field:
        if declaration.get('type') is not None:
            head += '(' + _paren_body(declaration).ljust(layout.paren) + ')'
        else:
            head += ' ' * layout.paren_field

    attributes = _ordered_attributes(declaration)
    for column, width in enumerate(layout.attributes):
        if column < len(attributes):
            head += ', ' + attributes[column].ljust(width)
        else:
            head += ' ' * (2 + width)
    head += ' :: '

    variables = _variables_of(declaration)
    if not variables:
        # No variables to lay out — emit the head alone rather than dropping
        # the statement.  `parse_declaration` should not produce this, but a
        # declaration dict synthesized by a code generator can.
        return head.rstrip() + '\n'
    columns = _variables_per_row(declaration, per_row)
    # Continuation lines carry the `&` five columns in from the block indent,
    # then run out to the variable column.
    continuation = (indent + ' ' * 5 + '&').ljust(
        layout.prefix_width(indent, openmp_any))

    lines = []
    for start in range(0, len(variables), columns):
        row      = variables[start:start + columns]
        last_row = start + columns >= len(variables)
        text     = head if start == 0 else continuation
        for column, variable in enumerate(row):
            is_last = last_row and column == len(row) - 1
            text += _emit_variable(
                declaration, variable, layout, column, pad=not is_last)
            if not is_last:
                text += ', '
        lines.append(text.rstrip() if last_row else text + '&')
    return ''.join(line + '\n' for line in lines)


def _choose_variables_per_row(declarations, indent, openmp_any):
    """Largest variables-per-row (up to the maximum) whose rows fit the ruler.

    Chosen once for the whole group, because the variable column widths are
    shared across it — varying it per statement would break the alignment those
    widths exist to produce.  A group that cannot fit even one variable per row
    emits one per row and overflows; there is nowhere else to break.
    """
    for per_row in range(VARIABLES_PER_ROW_MAX, 0, -1):
        layout = _Layout(declarations, per_row)
        widest = 0
        for declaration in declarations:
            text = _emit_declaration(
                declaration, layout, indent, openmp_any, per_row)
            widest = max(widest,
                         max((len(line) for line in text.splitlines()),
                             default=0))
        if widest <= LINE_LENGTH_MAX:
            return per_row
    return 1


def _minimum_width(declaration, indent):
    """Width of this declaration at its narrowest — every column unpadded."""
    width = len(indent) + len(declaration['intrinsic'])
    if declaration.get('type') is not None:
        width += len(_paren_body(declaration)) + 2
    for attribute in _ordered_attributes(declaration):
        width += 2 + len(attribute)
    width += 4                                        # " :: "
    longest = 0
    for variable in _variables_of(declaration):
        name, operator, value = _variable_fields(declaration, variable)
        longest = max(longest,
                      len(name) + len(_operator_text(operator)) + len(value))
    return width + longest


def _is_unformattable(declaration, indent):
    """True if no layout of this declaration can fit within the ruler.

    A single variable carrying a large initializer — a `reshape([...])` data
    table, say — cannot be broken by a formatter that only breaks *between*
    variables.  `get_fortran_line` has already joined away the continuation
    lines its author used, so re-emitting it from the parsed fields would
    collapse a carefully wrapped table onto one enormous line.  Such
    declarations are written back verbatim instead, and are kept out of the
    column-width calculation so that they cannot pad their neighbours out to
    their own width.
    """
    return (declaration.get('rawText') is not None
            and _minimum_width(declaration, indent) > LINE_LENGTH_MAX)


def _is_type_definition(declaration):
    """True if this is a derived-type *definition*, not a variable declaration.

    `type, bind(c) :: unitType` opens a derived type; `type(unitType) :: x`
    declares a variable of one.  The unit parser normally claims the former, so
    it never reaches a declaration node — but it recognizes the statement by
    pattern, and a few in the tree carry enough internal padding
    (`type, extends(hdf5Group             ) :: hdf5File`) to slip past it and
    land here instead.

    Reformatting such a statement would tidy away exactly the padding that hid
    it, so on the next parse the unit parser *would* claim it and the shape of
    the tree would change underneath the formatter.  There are four in the tree;
    leaving them untouched is simpler and safer than teaching the emitter to
    reproduce their disguise.
    """
    return (declaration.get('intrinsic') in ('type', 'class')
            and declaration.get('type') is None)


def _is_separator(node):
    """True if `node` is code that may sit *inside* an alignment group.

    `_pass_declarations` ends a declaration block at any line it does not
    recognize as a declaration, so a comment or `#ifdef` between two runs of
    declarations splits them into separate nodes.  Hand-formatted Galacticus
    source aligns straight through both, so the formatter rejoins runs
    separated only by comments, preprocessor directives and blank lines.
    Anything else — a real statement — genuinely ends the group.
    """
    if node.get('type') != 'code':
        return False
    for line in node.get('content', '').splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith(('!', '#')):
            return False
    return True


def declaration_groups(node):
    """Yield runs of sibling declaration nodes that share one alignment.

    Each group is a list of `declaration` nodes, in source order, which may be
    separated in the sibling chain by comment/preprocessor/blank code nodes.
    """
    from Galacticus.Build.SourceTree import children

    group     = []
    pending   = []          # separators seen since the last declaration node
    for child in children(node):
        if child.get('type') == 'declaration':
            group.append(child)
            pending = []
        elif group and _is_separator(child):
            pending.append(child)
        else:
            if group:
                yield group
            group   = []
            pending = []
    if group:
        yield group


def format_declarations(group):
    """Re-emit a group of declaration nodes in the house style.

    Rewrites each node's `firstChild['content']` in place.  Column widths are
    computed across the whole group so that declarations separated by a comment
    or a `#ifdef` still line up with one another.
    """
    declarations = [d for node in group for d in node.get('declarations', [])]
    if not declarations:
        return
    if any(_is_type_definition(d) for d in declarations):
        # Leave the whole group alone — see `_is_type_definition`.
        return

    # The block indent is the smallest indent of any declaration line in the
    # group.  As with `use` blocks, the minimum is what makes this idempotent:
    # the emitter indents non-OpenMP declarations three further columns when
    # the group also contains `!$` ones, and reading the indent off the first
    # line would pick that padding up again on every pass.
    indents = []
    for node in group:
        for line in node['firstChild']['content'].splitlines():
            stripped = line.lstrip()
            if not stripped or stripped.startswith('#'):
                continue
            # `!$` is an OpenMP sentinel, not a comment: such a line *is* a
            # declaration and sits at the base indent, while its unconditional
            # neighbours are padded three columns past it.  Skipping it here
            # would leave only the padded lines to measure, and the block would
            # march three columns right on every pass.
            if stripped.startswith('!') and not stripped.startswith('!$'):
                continue
            indents.append(re.match(r'^(\s*)', line).group(1))
    indent = min(indents, key=len) if indents else '  '

    openmp_any = any(d.get('openMP') for d in declarations)

    # Declarations too wide for any layout are written back as they were, and
    # take no part in sizing the columns.
    laid_out = [d for d in declarations if not _is_unformattable(d, indent)]
    if not laid_out:
        return
    per_row = _choose_variables_per_row(laid_out, indent, openmp_any)
    layout  = _Layout(laid_out, per_row)

    for node in group:
        content = indent + 'implicit none\n' if node.get('implicitNone') else ''
        for declaration in node.get('declarations', []):
            if _is_unformattable(declaration, indent):
                content += declaration['rawText']
            else:
                content += _emit_declaration(
                    declaration, layout, indent, openmp_any, per_row)
        node['firstChild']['content'] = content


def format_declarations_tree(tree):
    """Apply `format_declarations` to every alignment group in `tree`."""
    from Galacticus.Build.SourceTree import walk_tree

    for node in list(walk_tree(tree)):
        for group in declaration_groups(node):
            format_declarations(group)
