# Provides parsing of Fortran variable declaration lines.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the parseDeclaration() function from
# perl/Galacticus/Build/SourceTree/Parse/Declarations.pm

import re
import sys
import os
import copy
sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Fortran.Utils import INTRINSIC_DECLARATIONS, extract_variables


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
        # Look for a type specification: intrinsic(something)
        # But avoid matching dimension(:) or other attributes that follow a comma
        parts = type_attrs_raw.split(',', 1)  # Split at first comma
        first_part = parts[0].strip()
        rest = parts[1].strip() if len(parts) > 1 else ''

        # Check if first part has parentheses (type specification)
        m_type = re.search(r'\(\s*([^)]+)\s*\)', first_part)
        consumed_as_type = False
        if m_type and ('kind' in first_part.lower() or 'len' in first_part.lower() or
                       intrinsic in ['type', 'class']):
            type_val = m_type.group(1).strip()
            consumed_as_type = True

        # Attributes come after the comma; if there is no comma, anything in
        # `first_part` is an attribute *unless* we already consumed it as the
        # type-spec parens.
        if rest:
            attributes = extract_variables(rest, keep_qualifiers=True)
        elif not consumed_as_type and m_type is not None:
            # `first_part` carries parens that aren't a type-spec — must be
            # attributes (e.g. `dimension(:)` written without a leading comma).
            attributes = extract_variables(first_part, keep_qualifiers=True)
        elif not m_type:
            # No parens in `first_part` and no comma — nothing to extract.
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
    """
    target = variable_name.lower()
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'declaration':
            for declaration in child.get('declarations', []):
                for v in declaration.get('variables', []):
                    if v.lower() == target:
                        return True
        child = child.get('sibling')
    return False
