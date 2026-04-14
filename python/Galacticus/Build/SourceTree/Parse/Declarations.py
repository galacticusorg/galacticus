# Provides parsing of Fortran variable declaration lines.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the parseDeclaration() function from
# perl/Galacticus/Build/SourceTree/Parse/Declarations.pm

import re
import sys
import os
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
        if m_type and ('kind' in first_part.lower() or 'len' in first_part.lower() or
                       intrinsic in ['type', 'class']):
            type_val = m_type.group(1).strip()

        # Attributes come after comma or are everything else
        if rest:
            attributes = extract_variables(rest, keep_qualifiers=True)
        elif not m_type or not ('kind' in first_part.lower() or 'len' in first_part.lower()):
            # No comma but has parens - might be attributes without type
            attributes = extract_variables(first_part, keep_qualifiers=True)

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
