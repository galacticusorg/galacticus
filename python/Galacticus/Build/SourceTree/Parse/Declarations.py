# Provides parsing of Fortran variable declaration lines.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the parseDeclaration() function from
# perl/Galacticus/Build/SourceTree/Parse/Declarations.pm

import re
import copy

from Fortran.Utils import INTRINSIC_DECLARATIONS, extract_variables
from build.fortran_utils import extract_bracketed


def parse_declaration(line):
    """Parse a single Fortran variable declaration line.

    Mirrors Perl Galacticus::Build::SourceTree::Parse::Declarations::parseDeclaration().

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
          'variables'     : list of str — lowercase names, qualifiers preserved
          'variableNames' : list of str — original-case names, no qualifiers
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

        if rest:
            attributes = extract_variables(rest, keep_qualifiers=True)
        elif not consumed_as_type and first_part:
            # `first_part` carries parens that aren't a type-spec (rare —
            # `dimension(:)` written without a leading comma).
            attributes = extract_variables(first_part, keep_qualifiers=True)
        else:
            attributes = []

    variables = extract_variables(variables_raw, keep_qualifiers=True, lower_case=True)
    variable_names = extract_variables(variables_raw, keep_qualifiers=False, lower_case=False)

    return {
        'intrinsic':     intrinsic,
        'type':          type_val,
        'openMP':        bool(openmp_marker),
        'attributes':    attributes,
        'variables':     variables,
        'variableNames': variable_names,
    }


def build_declarations(node):
    """Build Fortran declaration text from a declaration node's structured data.

    Mirrors BuildDeclarations() from Parse/Declarations.pm.  Rewrites
    node['firstChild']['content'] in place.
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

    Mirrors AddDeclarations() from Parse/Declarations.pm.  The new declaration node
    is inserted after any existing moduleUse child, otherwise before the first child.
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

    Mirrors AddAttributes() from Parse/Declarations.pm.  If the variable shares a
    declaration line with others, that declaration is split so the target variable
    can receive attributes independently.
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

    Mirrors GetDeclaration() from Parse/Declarations.pm.  The returned dict is a
    deep copy with 'variables' reduced to the single target variable (original case).
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
                # Strip any initializer value before comparing (mirrors Perl
                # regex `m/^([^=]+)\s*=/ ? $1 : $_`).
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

    Mirrors DeclarationExists() from Parse/Declarations.pm.  Case-insensitive.
    Strips any `=…` initializer off each stored variable name before
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
