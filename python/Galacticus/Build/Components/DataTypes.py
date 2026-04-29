# Helpers for translating component-XML data descriptors into the matching
# Fortran type / attribute strings.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/DataTypes.pm.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils import intrinsic_types


# Minimal LaTeX special-character escape — covers the cases that the Perl
# `LaTeX::Encode` module produces for type strings ("integer", type names,
# etc.).  Matches the subset already used in
# python/Galacticus/Build/SourceTree/Process/FunctionClass/__init__.py.
_LATEX_ESCAPES = {
    '\\': r'\textbackslash{}',
    '_':  r'\_',
    '&':  r'\&',
    '%':  r'\%',
    '$':  r'\$',
    '#':  r'\#',
    '{':  r'\{',
    '}':  r'\}',
    '~':  r'\textasciitilde{}',
    '^':  r'\textasciicircum{}',
}


def _latex_encode(text):
    return ''.join(_LATEX_ESCAPES.get(c, c) for c in str(text))


def data_object_primitive_name(data_object, *, match_only=False):
    """Return `(name, type, attribute_list)` for `data_object`.

    Mirrors Perl `dataObjectPrimitiveName`:

    * `name`            — the Fortran type spec, e.g. `"double precision"`
                          or `"type(massDistribution)"`.
    * `type`            — a CamelCase label derived from
                          `data_object['type']`, suitable for embedding
                          inside generated identifiers.
    * `attribute_list`  — a leading-comma string of attributes
                          (`", dimension(:), allocatable"` etc.) or `""`.

    `match_only=True` suppresses the `allocatable` attribute on
    rank > 0 — used when the descriptor is being matched against an
    existing variable, not declared.
    """
    if 'type' not in data_object:
        raise RuntimeError(
            "DataTypes::dataObjectPrimitiveName: no 'type' specifier present"
        )
    has_rank  = 'rank'  in data_object
    has_shape = 'shape' in data_object
    if not (has_rank or has_shape):
        raise RuntimeError(
            "DataTypes::dataObjectPrimitiveName: no 'rank' or 'shape' specifier present"
        )
    if has_rank and has_shape:
        raise RuntimeError(
            "DataTypes::dataObjectPrimitiveName: can not have both 'rank' and 'shape' specifiers present"
        )

    type_label = data_object['type']
    if type_label in intrinsic_types:
        name = intrinsic_types[type_label]
    else:
        name = f"type({type_label})"

    type_camel = ''.join(_ucfirst(part) for part in type_label.split())

    attributes = []
    if has_rank and int(data_object['rank']) > 0:
        rank = int(data_object['rank'])
        attributes.append('dimension(' + ','.join([':'] * rank) + ')')
        if not match_only:
            attributes.append('allocatable')
    if has_shape:
        attributes.append(f"dimension({data_object['shape']})")

    attribute_list = (', ' + ', '.join(attributes)) if attributes else ''
    return name, type_camel, attribute_list


def data_object_doc_name(data_object):
    """Return the LaTeX-encoded display name used in component documentation.

    Mirrors Perl `dataObjectDocName`.  Always returns the `\\textcolor{red}
    {\\textless …\\textgreater}` form, with `\\textless` and `\\textgreater`
    typeset literally so the LaTeX renderer sees them as `<` and `>`.
    """
    if 'type' not in data_object:
        raise RuntimeError(
            "DataTypes::dataObjectDocName: no 'type' specifier present"
        )
    rank = int(data_object.get('rank') or 0)
    type_label = data_object['type']
    if type_label in intrinsic_types:
        body = _latex_encode(intrinsic_types[type_label])
    else:
        body = f"type({_latex_encode(type_label)})"
    if rank > 0:
        body += "[" + ",".join([":"] * rank) + "]"
    return f"\\textcolor{{red}}{{\\textless {body}\\textgreater}}"


def data_object_name(data_object):
    """Return the canonical `nodeData…` identifier for a property.

    Mirrors Perl `dataObjectName`.  The generated name carries the
    intrinsic-type-derived suffix (`Integer`, `DoublePrecision`, …) plus
    `Scalar` / `<rank>D` and an optional `Evolvable` marker.
    """
    if 'type' not in data_object:
        raise RuntimeError(
            "DataTypes::dataObjectName: no 'type' specifier present"
        )
    name = "nodeData"
    type_label = data_object['type']
    if type_label in intrinsic_types:
        name += ''.join(_ucfirst(part) for part in intrinsic_types[type_label].split())
    else:
        name += _ucfirst(type_label)
    # Perl: `exists($rank)` always takes the "Scalar" branch; the "<rank>D"
    # branch is unreachable in practice (it reads $rank under elsif !exists).
    # Mirrored faithfully — callers only ever pass descriptors with rank set.
    if 'rank' in data_object:
        name += "Scalar"
    elif type_label != "void":
        name += f"{data_object.get('rank', '0')}D"
    if data_object.get('isEvolvable'):
        name += "Evolvable"
    return name.replace(' ', '')


def data_object_definition(data_object, *, match_only=False):
    """Return `(declaration_dict, label)` describing a property's variable.

    Mirrors Perl `dataObjectDefinition`.  `declaration_dict` is suitable
    for handing to `format_variable_definitions` (carrying `intrinsic`,
    `type` if applicable, and `attributes`).  `label` is the
    CamelCase-ified type used elsewhere in the build pipeline.
    """
    if 'type' not in data_object:
        raise RuntimeError(
            "dataObjectDefinition: no 'type' specifier present"
        )
    if 'rank' not in data_object:
        raise RuntimeError(
            "dataObjectDefinition: no 'rank' specifier present"
        )

    type_label = data_object['type']
    declaration = {}
    if type_label in intrinsic_types:
        declaration['intrinsic'] = intrinsic_types[type_label]
    else:
        declaration['intrinsic'] = 'type'
        declaration['type'     ] = type_label

    label = _ucfirst(type_label)

    rank = int(data_object['rank'])
    attributes = []
    if rank > 0:
        attributes.append('dimension(' + ','.join([':'] * rank) + ')')
        if not match_only:
            attributes.append('allocatable')
    if attributes:
        declaration['attributes'] = attributes

    return declaration, label


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
