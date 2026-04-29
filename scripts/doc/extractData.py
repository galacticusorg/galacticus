#!/usr/bin/env python3
"""Scan Fortran90 source code and extract class/method data from .classes.xml files.

Andrew Benson 12-Mar-2010 (original Perl); Python port 2026.

Usage: extractData.py <sourceDir> <outputRoot>
"""

import os
import re
import sys
import xml.etree.ElementTree as ET

_exec_path = os.environ.get('GALACTICUS_EXEC_PATH', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.join(_exec_path, 'python'))
from latex_utils import latex_encode  # noqa: E402
from list_utils  import as_array      # noqa: E402


# ---------------------------------------------------------------------------
# XML::Simple compatibility
# ---------------------------------------------------------------------------

def _elem_to_python(elem):
    """Recursively convert an ElementTree element to a Python dict/str.

    Replicates XML::Simple's XMLin default behaviour (no ForceArray):
      - Element with no children and only text  → str
      - Element with children                   → dict
      - Multiple sibling elements with same tag → list
    """
    children = list(elem)
    if not children:
        return (elem.text or '').strip()

    result = {}
    tag_counts = {}
    for child in children:
        tag_counts[child.tag] = tag_counts.get(child.tag, 0) + 1

    for child in children:
        val = _elem_to_python(child)
        if tag_counts[child.tag] == 1:
            result[child.tag] = val
        else:
            result.setdefault(child.tag, [])
            result[child.tag].append(val)

    return result


def xml_simple_in(path):
    """Parse an XML file written by XML::Simple XMLout(…, NoAttr=>1).

    Returns a dict keyed by the direct children of the root element.
    """
    tree = ET.parse(path)
    root = tree.getroot()
    result = {}
    for child in root:
        result[child.tag] = _elem_to_python(child)
    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

FUNCTION_CLASS_EXCLUDES = {
    'allowedParameters', 'autoHook', 'deepCopy', 'deepCopyFinalize',
    'deepCopyReset', 'descriptor', 'destructor', 'hashedDescriptor',
    'objectType', 'stateRestore', 'stateStore', 'assignment(=)',
}
_FUNCTION_CLASS_EXCLUDES_LC = {n.lower() for n in FUNCTION_CLASS_EXCLUDES}


def declaration_builder(type_val, variables=True):
    """Build a \\mono{...} LaTeX string for a Fortran type declaration.

    Port of Perl declarationBuilder().  type_val is either a dict (complex type)
    or one of the literal strings 'subroutine' / 'void' for the void-return
    case.  Both string forms appear in the wild — `_process_function` emits
    `'subroutine'` (the AST node type) for type-bound subroutine bindings,
    and FunctionClass auto-generated method stubs (autoHook, descriptor,
    deepCopy, …) carry `'void'` straight from the directive.
    """
    if isinstance(type_val, dict):
        inner = type_val.get('intrinsic', '')
        t = type_val.get('type')
        if t and not isinstance(t, dict):
            inner += '(' + latex_encode(t) + ')'
        attrs = as_array(type_val.get('attributes'))
        if attrs:
            inner += ', ' + ', '.join(latex_encode(a) for a in attrs)
        if variables:
            vs = as_array(type_val.get('variables'))
            if vs:
                inner += ' :: ' + ', '.join(latex_encode(v) for v in vs)
    elif type_val in ('subroutine', 'void'):
        inner = 'void'
    else:
        raise ValueError(f'declaration_builder: unknown type {type_val!r}')
    return f'\\mono{{{inner}}}'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print('Usage: extractData.py <sourceDir> <outputRoot>', file=sys.stderr)
        sys.exit(1)

    # source_dir = sys.argv[1]   # not used directly, kept for compatibility
    output_root = sys.argv[2]
    build_path  = os.environ.get('BUILDPATH')
    if not build_path:
        print('Error: BUILDPATH environment variable is not set.', file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Parse all *.classes.xml files from the build directory.
    # ------------------------------------------------------------------
    classes = {}
    for fname in os.listdir(build_path):
        if not fname.endswith('.classes.xml'):
            continue
        path = os.path.join(build_path, fname)
        try:
            file_classes = xml_simple_in(path)
        except ET.ParseError as exc:
            print(f'Warning: failed to parse {path}: {exc}', file=sys.stderr)
            continue
        for class_name, class_data in file_classes.items():
            classes[class_name] = class_data

    # ------------------------------------------------------------------
    # Build child/parent relationships; detect functionClass descendants;
    # remove missing methods found in ancestors.
    # ------------------------------------------------------------------
    for class_name in sorted(classes):
        classes[class_name]['isFunctionClass'] = False
        parent_name = classes[class_name].get('extends')
        if not parent_name:
            continue

        # Add this class to the parent's children list.
        classes.setdefault(parent_name, {})
        classes[parent_name].setdefault('children', [])
        classes[parent_name]['children'].append(class_name)

        # Detect functionClass descendants.
        ancestor = parent_name
        while ancestor:
            if ancestor == 'functionClass':
                classes[class_name]['isFunctionClass'] = True
                break
            ancestor = classes.get(ancestor, {}).get('extends')

        # Remove missing methods that are actually defined in an ancestor.
        missing_raw = classes[class_name].get('missingMethods', '')
        if isinstance(missing_raw, dict):
            # XML::Simple wraps an element that has no text as a dict; treat as empty.
            continue
        missing_methods = missing_raw.split()
        if not missing_methods:
            continue

        remaining = []
        for method in missing_methods:
            # Auto-generated functionClass methods (autoHook, descriptor,
            # deepCopy, stateStore, …) are inherent to every functionClass —
            # they are emitted by the build's FunctionClass processor, not
            # described anywhere in the source.  Treat them as resolved so
            # they don't surface as "missing method descriptions" warnings
            # for every concrete child class.
            if (classes[class_name]['isFunctionClass']
                    and method.lower() in _FUNCTION_CLASS_EXCLUDES_LC):
                continue
            found = False
            ancestor = parent_name
            while ancestor:
                anc_data = classes.get(ancestor, {})
                # Check descriptions list.
                for desc in as_array(anc_data.get('descriptions')):
                    if isinstance(desc, dict) and desc.get('method', '').lower() == method.lower():
                        found = True
                        break
                if found:
                    break
                # Check genericUses.
                generic_uses = anc_data.get('genericUses', '')
                if not isinstance(generic_uses, dict):
                    if method.lower() in [g.lower() for g in generic_uses.split()]:
                        found = True
                        break
                ancestor = classes.get(ancestor, {}).get('extends')
            if not found:
                remaining.append(method)
        classes[class_name]['missingMethods'] = ' '.join(remaining)

    # Warn about any unresolved missing methods.
    for class_name in sorted(classes):
        missing_raw = classes[class_name].get('missingMethods', '')
        if isinstance(missing_raw, dict):
            continue
        missing_methods = missing_raw.split()
        if missing_methods:
            print(f"Warning: missing method descriptions in class '{class_name}':")
            for m in missing_methods:
                print(f'\t{m}')

    # ------------------------------------------------------------------
    # Write the methods LaTeX file.
    # ------------------------------------------------------------------
    output_path = output_root + 'Methods.tex'
    with open(output_path, 'w', encoding='utf-8') as fh:
        for class_name in sorted(classes):
            cls = classes[class_name]
            name = cls.get('name', class_name)

            fh.write(
                f'\\subsection{{\\large \\mono{{{latex_encode(name)}}}}}'
                f'\\label{{class:{name}}}\\hyperdef{{class}}{{{name}}}{{}}\n\n'
            )

            if cls.get('isFunctionClass'):
                fh.write(f'\\emph{{Physics model:}} \\refPhysics{{{name}}}\n\n')

            if 'extends' in cls:
                fh.write(f'\\noindent\\emph{{Parent class:}} \\refClass{{{cls["extends"]}}}\n\n')

            if 'children' in cls:
                sorted_children = sorted(cls['children'])
                fh.write('\\noindent\\emph{Child classes:}\n\n\\begin{tabular}{ll}\n')
                for i in range(0, len(sorted_children), 2):
                    fh.write(f'\\refClass{{{sorted_children[i]}}}')
                    if i + 1 < len(sorted_children):
                        fh.write(f' & \\refClass{{{sorted_children[i + 1]}}}')
                    fh.write('\\\\\n')
                fh.write('\\end{tabular}\n\n')

            descriptions = as_array(cls.get('descriptions'))
            if descriptions or cls.get('isFunctionClass'):
                fh.write('\\begin{description}\n')
                if cls.get('isFunctionClass'):
                    fh.write(
                        '\\item[] All standard \\refClass{functionClass} methods'
                        ' (see \\S\\ref{sec:functionClassAll})\n'
                    )

                if descriptions:
                    methods = []
                    for method in descriptions:
                        if not isinstance(method, dict):
                            continue
                        if cls.get('isFunctionClass') and method.get('method') in FUNCTION_CLASS_EXCLUDES:
                            continue
                        methods.append(method)

                    for method in methods:
                        method_name = method.get('method', '')
                        if 'type' not in method:
                            print(
                                f"Warning: missing function type for method '{method_name}'"
                                f" of class '{class_name}'"
                            )
                            continue
                        # Escape underscores in the method label (but not already-escaped ones).
                        method_label = re.sub(r'([^\\])_', r'\1\\_', method_name)
                        description  = method.get('description', '').rstrip('\n')
                        fh.write(f'\\item[]\\mono{{{method_label}}} {description}\n')
                        fh.write('\\begin{itemize}\n')
                        fh.write(f'\\item Return type: {declaration_builder(method["type"], variables=False)}\n')

                        arg_lists = as_array(method.get('argumentList'))
                        arg_groups = as_array(method.get('arguments'))
                        for i, arg_list in enumerate(arg_lists):
                            fh.write('\\item Interface: ')
                            if isinstance(arg_list, dict) or not arg_list:
                                fh.write('\\mono{()}\\\\\n')
                            else:
                                fh.write(f'\\mono{{({latex_encode(arg_list)})}}\\\\\n')
                                if i < len(arg_groups) and isinstance(arg_groups[i], dict):
                                    for arg in as_array(arg_groups[i].get('argument')):
                                        fh.write(declaration_builder(arg) + '\\\\\n')
                        fh.write('\\end{itemize}\n')

                fh.write('\\end{description}\n')


if __name__ == '__main__':
    main()
