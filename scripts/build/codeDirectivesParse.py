#!/usr/bin/env python3
# Scans source code for "!![...!!]" directives and generates a Makefile.
# Python port of scripts/build/codeDirectivesParse.pl.
# Andrew Benson (ported to Python 2026)

import os
import pickle
import re
import sys
import filecmp
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.Directives import extract_directives


def update_if_changed(old_path, new_path):
    """Replace old_path with new_path only if content differs; delete new_path
    otherwise. Mirrors Perl File::Changes::Update()."""
    if not os.path.exists(old_path):
        os.replace(new_path, old_path)
        return
    if filecmp.cmp(old_path, new_path, shallow=False):
        os.remove(new_path)
    else:
        os.replace(new_path, old_path)


def _indent(element, level=0):
    pad = '\n' + '  ' * level
    if len(element):
        if not element.text or not element.text.strip():
            element.text = pad + '  '
        for child in element:
            _indent(child, level + 1)
        last = list(element)[-1]
        if not last.tail or not last.tail.strip():
            last.tail = pad
    if level and (not element.tail or not element.tail.strip()):
        element.tail = pad


def _build_element(tag, value):
    element = ET.Element(tag)
    if isinstance(value, dict):
        for key, sub in value.items():
            if isinstance(sub, list):
                for item in sub:
                    element.append(_build_element(key, item))
            else:
                element.append(_build_element(key, sub))
    elif isinstance(value, list):
        for item in value:
            element.append(_build_element(tag, item))
    elif value is None:
        pass
    else:
        element.text = str(value)
    return element


def dict_to_xml(root_name, value):
    """Serialize a dict to an XML string, mirroring XML::Simple XMLout with
    NoAttr=>1 semantics (every field becomes a child element, never an
    attribute). Output is pretty-printed with two-space indentation."""
    root = _build_element(root_name, value)
    _indent(root)
    return ET.tostring(root, encoding='unicode') + '\n'


def add_implicit_directives(directive, file_identifier, file_name,
                            file_path_name, directives_per_file):
    """Append implicit directives required by certain directive kinds
    (currently functionClass: stateful store/retrieve and destroy tasks)."""
    is_function_class = directive.get('rootElementType') == 'functionClass'
    implicit = {
        'stateful': {
            'always': is_function_class,
            'tasks': ['galacticusStateRetrieveTask', 'galacticusStateStoreTask'],
        },
        'functionClassDestroy': {
            'always': is_function_class and 'functionClassDestroy' not in directive,
            'tasks': ['functionClassDestroyTask'],
        },
    }
    non_include = directives_per_file[file_identifier].setdefault('nonIncludeDirectives', {})
    for key in sorted(implicit):
        spec = implicit[key]
        if directive.get(key) == 'yes' or spec['always']:
            for task in spec['tasks']:
                entry = non_include.setdefault(task, {})
                entry.setdefault('files', []).append(file_path_name)
                entry.setdefault('dependency', []).append(file_name)


def main(install_directory):
    build_path = os.environ['BUILDPATH']
    blob_path = os.path.join(build_path, 'codeDirectives.blob')
    directives_per_file = {}
    update_time = None
    if os.path.exists(blob_path):
        try:
            with open(blob_path, 'rb') as handle:
                directives_per_file = pickle.load(handle)
            update_time = os.path.getmtime(blob_path)
        except Exception:
            directives_per_file = {}
            update_time = None

    source_directory = os.path.join(install_directory, 'source')
    source_file_pattern = re.compile(r'\.(f|f90|c|cpp|h)$', re.IGNORECASE)
    source_file_names = sorted(
        f for f in os.listdir(source_directory)
        if source_file_pattern.search(f) and not f.startswith('.#')
    )

    def identifier_for(path):
        ident = path.replace('/', '_')
        ident = re.sub(r'^\._?', '', ident)
        return ident

    file_identifiers = set(
        identifier_for(os.path.join(source_directory, f)) for f in source_file_names
    )
    force_rescan = any(ident not in directives_per_file for ident in file_identifiers)
    if not force_rescan:
        force_rescan = any(ident not in file_identifiers for ident in directives_per_file)

    include_pattern = re.compile(r"^\s*include\s*['\"]([^'\"]+)['\"]\s*$")
    include_ref_pattern = re.compile(r"""^\s*\#??include\s*["'<](.+)["'>]""", re.IGNORECASE | re.MULTILINE)

    for file_name in source_file_names:
        full_path = os.path.join(source_directory, file_name)
        file_identifier = identifier_for(full_path)
        rescan = True
        if update_time is not None and file_identifier in directives_per_file:
            files_cached = directives_per_file[file_identifier].get('files', [])
            rescan = False
            for f in files_cached:
                try:
                    if os.stat(f).st_mtime > update_time:
                        rescan = True
                        break
                except OSError:
                    pass
        if not (rescan or force_rescan):
            continue
        directives_per_file.pop(file_identifier, None)
        entry = directives_per_file.setdefault(file_identifier, {})
        entry['files'] = [full_path]
        stack = [full_path]
        while stack:
            file_path_name = stack.pop()
            try:
                with open(file_path_name, 'r', encoding='utf-8', errors='replace') as handle:
                    lines = handle.readlines()
            except OSError as err:
                raise SystemExit(
                    "codeDirectivesParse.py: can not open input file: " + str(err)
                )
            for line in lines:
                m = include_pattern.match(line)
                if not m:
                    continue
                included = os.path.join(source_directory, m.group(1))
                included = re.sub(r'\.inc$', '.Inc', included)
                if os.path.exists(included):
                    stack.append(included)
                    entry['files'].append(included)
            for directive in extract_directives(file_path_name, '*', set_root_element_type=True):
                root_type = directive.get('rootElementType')
                if root_type == 'include':
                    directive['source'] = file_path_name
                    content_match = include_ref_pattern.search(directive.get('content', ''))
                    if content_match:
                        included_name = content_match.group(1)
                        directive['fileName'] = os.path.join(
                            build_path, re.sub(r'\.inc$', '.Inc', included_name)
                        )
                    directive.pop('content', None)
                    directive_key = (
                        directive.get('name', directive.get('directive', ''))
                        + '.' + directive.get('type', '')
                    )
                    entry.setdefault('includeDirectives', {})[directive_key] = {
                        'source': file_path_name,
                        'fileName': directive.get('fileName'),
                        'xml': dict_to_xml(root_type, directive),
                    }
                else:
                    non_inc = entry.setdefault('nonIncludeDirectives', {})
                    non_inc.setdefault(root_type, {}).setdefault('files', []).append(full_path)
                    if root_type == 'functionClass':
                        preprocessed = re.sub(
                            r'\.F90$', '.p.F90',
                            os.path.join(build_path, file_name)
                        )
                        entry.setdefault('functionClasses', {})[directive['name']] = preprocessed
                        add_implicit_directives(
                            directive, file_identifier, preprocessed, preprocessed,
                            directives_per_file,
                        )

    include_directives = {}
    non_include_directives = {}
    function_classes = {}
    for ident, record in directives_per_file.items():
        for key, value in record.get('includeDirectives', {}).items():
            include_directives[key] = value
        for directive, details in record.get('nonIncludeDirectives', {}).items():
            target = non_include_directives.setdefault(directive, {})
            for field in ('files', 'dependency'):
                if field in details:
                    target.setdefault(field, []).extend(details[field])
        for key, value in record.get('functionClasses', {}).items():
            function_classes[key] = value

    for directive, details in non_include_directives.items():
        for field in ('files', 'dependency'):
            if field in details:
                details[field] = sorted(set(details[field]))

    output_directives = {
        directive: {'file': list(details['files'])}
        for directive, details in non_include_directives.items()
        if 'files' in details
    }
    locations_xml_path = os.path.join(build_path, 'directiveLocations.xml')
    tmp_locations_path = locations_xml_path + '.tmp'
    with open(tmp_locations_path, 'w', encoding='utf-8') as handle:
        handle.write(dict_to_xml('directives', output_directives))
    update_if_changed(locations_xml_path, tmp_locations_path)

    makefile_path = os.path.join(build_path, 'Makefile_Directives')
    with open(makefile_path, 'w', encoding='utf-8') as makefile:
        for directive in sorted(include_directives):
            include_entry = include_directives[directive]
            file_name = re.sub(r'\.inc$', '.Inc', include_entry['fileName'])
            extra_dependencies = []
            function_match = re.match(r'^([a-zA-Z0-9_]+)\.function$', directive)
            if function_match:
                base = function_match.group(1)
                files = non_include_directives.get(base, {}).get('files', [])
                extra_dependencies.extend(sorted(files))
            pair_match = re.match(r'^([a-zA-Z0-9_]+)\.(moduleUse|functionCall)$', directive)
            if pair_match:
                base = pair_match.group(1)
                deps = non_include_directives.get(base, {}).get('dependency')
                if deps:
                    extra_dependencies.extend(deps)
            makefile.write(
                file_name + '.up: ' + os.path.join(build_path, directive + '.xml')
                + ' ' + ' '.join(extra_dependencies)
                + ' $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml\n'
            )
            makefile.write(
                '\t./scripts/build/buildCode.pl ' + install_directory
                + ' ' + os.path.join(build_path, directive + '.xml') + '\n'
            )
            makefile.write(file_name + ': ' + file_name + '.up\n')
            makefile.write('\n')
            directive_xml_path = os.path.join(build_path, directive + '.xml')
            tmp_directive_path = directive_xml_path + '.tmp'
            with open(tmp_directive_path, 'w', encoding='utf-8') as directive_file:
                directive_file.write(include_entry['xml'])
            update_if_changed(directive_xml_path, tmp_directive_path)

        for directive_name in sorted(function_classes):
            dep_files = sorted(non_include_directives.get(directive_name, {}).get('files', []))
            makefile.write(
                function_classes[directive_name] + '.up: '
                + ' '.join(dep_files) + '\n\n'
            )

        use_deps_path = os.path.join(build_path, 'Makefile_Use_Dependencies')
        include_file_names = sorted(
            re.sub(r'\.inc$', '.Inc', include_directives[key]['fileName'])
            for key in include_directives
        )
        makefile.write(use_deps_path + ': ' + ' '.join(include_file_names) + '\n\n')

        makefile.write('-include ' + os.path.join(build_path, 'Makefile_Component_Includes') + '\n')
        makefile.write(
            os.path.join(build_path, 'Makefile_Component_Includes') + ': '
            + os.path.join(build_path, 'objects.nodes.components.Inc') + '\n\n'
        )

    tmp_blob_path = blob_path + '.tmp'
    with open(tmp_blob_path, 'wb') as handle:
        pickle.dump(directives_per_file, handle)
    update_if_changed(blob_path, tmp_blob_path)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: codeDirectivesParse.py <installDirectory>')
    main(sys.argv[1])
