#!/usr/bin/env python3
"""Generate an XSD for Galacticus parameter FILES from the parameter catalog.

The schema is intended for editor assistance (e.g. the VSCode XML extension):
it declares each ``functionClass`` selector element with its ``value`` attribute
restricted to the valid implementation labels, and each *unambiguous*
enumeration-valued parameter with its ``value`` restricted to the allowed
labels.  All element content is otherwise lax (``processContents="lax"``), so:

* any nesting is permitted, and unknown elements (global/meta parameters,
  components, custom parameters) are allowed;
* known elements (selectors, enum parameters) are validated *wherever* they
  appear, because lax validation matches them against their global declarations.

This is a pragmatic XSD 1.0 schema (works in any XSD-aware editor).  Precise,
per-implementation validation (a parameter must be accepted by the *selected*
implementation) is provided separately by scripts/build/parameterValidate.py.

Usage:
    parameterSchema.py <sourceDirectory> [<outputFile>]

`<outputFile>` defaults to `<sourceDirectory>/schema/parameters.xsd`.

Andrew Benson (2026).
"""

import collections
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, 'python'))

from Galacticus.Parameters.catalog import build_catalog


# Root elements of non-parameter files that live in the same directories as
# parameter files (so an editor association scoped to those directories does not
# spuriously flag them). Their content is not validated.
_TOLERATED_ROOTS = ('changes', 'mergerTrees', 'parameterGrid', 'tree')


def _escape(text):
    return (str(text).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;'))


def _value_enum_attribute(values):
    lines = ['    <xs:attribute name="value" use="optional">',
             '     <xs:simpleType>',
             '      <xs:restriction base="xs:string">']
    for value in values:
        lines.append(f'       <xs:enumeration value="{_escape(value)}"/>')
    lines += ['      </xs:restriction>',
              '     </xs:simpleType>',
              '    </xs:attribute>']
    return lines


def _value_constrained_element(name, values):
    """A global element whose `value` is restricted to `values`, with lax
    children and arbitrary other attributes (id, idRef, ignoreWarnings, ...).
    `mixed` content tolerates the `<name>value</name>` text form."""
    lines = [f'  <xs:element name="{_escape(name)}">',
             '   <xs:complexType mixed="true">',
             '    <xs:sequence>',
             '     <xs:any minOccurs="0" maxOccurs="unbounded" processContents="lax"/>',
             '    </xs:sequence>']
    lines += ['   ' + ln for ln in _value_enum_attribute(values)]
    lines += ['    <xs:anyAttribute processContents="lax"/>',
              '   </xs:complexType>',
              '  </xs:element>']
    return lines


def _generic_element(name):
    """A global element with an unconstrained `value` -- used for names that are
    both a functionClass base and a scalar parameter (so the value may be either
    an implementation label or an ordinary value)."""
    return [f'  <xs:element name="{_escape(name)}">',
            '   <xs:complexType mixed="true">',
            '    <xs:sequence>',
            '     <xs:any minOccurs="0" maxOccurs="unbounded" processContents="lax"/>',
            '    </xs:sequence>',
            '    <xs:attribute name="value" type="xs:string" use="optional"/>',
            '    <xs:anyAttribute processContents="lax"/>',
            '   </xs:complexType>',
            '  </xs:element>']


def _enum_parameter_elements(catalog):
    """Return ``{paramName: [allowed values]}`` for parameters that are
    unambiguously enumeration-valued: every occurrence of the name is linked to
    the *same* single enumeration, the name is not a functionClass base, and the
    enumeration's labels are known."""
    enumerations = catalog.get('enumerations', {})
    base_names = set(catalog.get('functionClasses', {}))
    name_enums = collections.defaultdict(set)
    name_has_free = collections.defaultdict(bool)
    for impl in catalog.get('implementations', {}).values():
        for parameter in impl.get('parameters', []):
            name = parameter.get('name')
            if not name:
                continue
            if parameter.get('enumeration'):
                name_enums[name].add(parameter['enumeration'])
            else:
                name_has_free[name] = True
    result = {}
    for name, enums in name_enums.items():
        if name in base_names or name_has_free.get(name) or len(enums) != 1:
            continue
        values = enumerations.get(next(iter(enums)))
        if values:
            result[name] = sorted(set(values))
    return result


def build_schema(catalog):
    function_classes = catalog.get('functionClasses', {})
    enum_parameters = _enum_parameter_elements(catalog)

    out = ['<?xml version="1.0"?>',
           '<!-- Galacticus parameter-file schema. GENERATED from the parameter',
           '     catalog by scripts/build/parameterSchema.py (do not edit by hand).',
           '     For editor assistance only; precise validation is done by',
           '     scripts/build/parameterValidate.py. -->',
           '<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">',
           '']

    # Names that are both a functionClass base and a scalar parameter: their
    # value may be an implementation label OR an ordinary value, so leave it
    # unconstrained to avoid false positives.
    parameter_names = {p.get('name')
                       for impl in catalog.get('implementations', {}).values()
                       for p in impl.get('parameters', []) if p.get('name')}

    out.append('  <!-- functionClass selectors: value restricted to valid implementations. -->')
    for base in sorted(function_classes):
        # Sort labels so the emitted order is stable regardless of catalog
        # ordering (keeps regeneration diffs minimal).
        labels = sorted(function_classes[base].get('implementations', []))
        if not labels:
            continue
        if base in parameter_names:
            out += _generic_element(base)
        else:
            out += _value_constrained_element(base, labels)
    out.append('')

    out.append('  <!-- Enumeration-valued parameters: value restricted to allowed labels. -->')
    for name in sorted(enum_parameters):
        out += _value_constrained_element(name, enum_parameters[name])
    out.append('')

    out += ['  <!-- Root element. Any parameters are permitted; known elements above',
            '       are validated wherever they appear. -->',
            '  <xs:element name="parameters">',
            '   <xs:complexType mixed="true">',
            '    <xs:sequence>',
            '     <xs:any minOccurs="0" maxOccurs="unbounded" processContents="lax"/>',
            '    </xs:sequence>',
            '    <xs:anyAttribute processContents="lax"/>',
            '   </xs:complexType>',
            '  </xs:element>',
            '']

    out.append('  <!-- Non-parameter file roots that share parameter directories')
    out.append('       (merger trees, changes files, parameter grids): declared so that a')
    out.append('       directory-scoped editor association does not flag them; their content')
    out.append('       is not validated. -->')
    for root_name in _TOLERATED_ROOTS:
        out += _generic_element(root_name)

    out += ['', '</xs:schema>', '']
    return '\n'.join(out)


def main(argv):
    if len(argv) not in (2, 3):
        print("Usage: parameterSchema.py <sourceDirectory> [<outputFile>]",
              file=sys.stderr)
        return 1
    source_directory = argv[1]
    output_path = argv[2] if len(argv) == 3 else os.path.join(
        source_directory, 'schema', 'parameters.xsd')
    os.environ.setdefault('GALACTICUS_EXEC_PATH', source_directory)
    catalog = build_catalog(os.path.join(source_directory, 'source'))
    with open(output_path, 'w') as fh:
        fh.write(build_schema(catalog))
    n_fc = sum(1 for b in catalog['functionClasses'].values()
               if b.get('implementations'))
    print(f"wrote {output_path} ({n_fc} functionClass selectors, "
          f"{len(_enum_parameter_elements(catalog))} enumeration parameters)",
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
