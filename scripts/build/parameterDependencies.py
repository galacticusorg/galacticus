#!/usr/bin/env python3
# Scan source files for input parameter definitions for a given executable.
# Python port of scripts/build/parameterDependencies.pl.
# Andrew Benson (ported to Python 2026)

import os
import pickle
import re
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.Directives import extract_directives


INCLUDE_FILES_EXCLUDED = {'fftw3.f03'}
INCLUDE_LINE = re.compile(r"^\s*include\s*['\"]([^'\"]+)['\"]")


def _identifier_for(path):
    ident = path.replace('/', '_')
    return re.sub(r'^\._?', '', ident)


def _read_lines(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as handle:
            return handle.read().splitlines()
    except OSError:
        return []


def main(source_directory, target_name):
    build_path = os.environ['BUILDPATH']

    # Read the list of object files that make up this target.
    dependency_file_name = os.path.join(
        build_path, re.sub(r'\.(exe|o)$', '.d', target_name)
    )
    object_files = set()
    build_prefix = build_path.rstrip('/') + '/'
    for line in _read_lines(dependency_file_name):
        if line.startswith(build_prefix) and line.endswith('.o'):
            object_files.add(line[len(build_prefix):])

    blob_file_name = re.sub(r'\.(exe|o)$', '.blob', target_name)
    blob_path = os.path.join(build_path, blob_file_name)
    parameters_per_file = {}
    update_time = None
    if os.path.exists(blob_path):
        try:
            with open(blob_path, 'rb') as handle:
                parameters_per_file = pickle.load(handle)
            update_time = os.path.getmtime(blob_path)
        except Exception:
            parameters_per_file = {}
            update_time = None

    src_subdir = os.path.join(source_directory, 'source')
    for file_name in sorted(os.listdir(src_subdir)):
        if file_name.startswith('.#'):
            continue
        if not re.search(r'\.(F90|cpp)$', file_name):
            continue
        object_file_name = re.sub(r'\.(F90|cpp)$', '.o', file_name)
        if object_file_name not in object_files:
            continue

        is_fortran = file_name.endswith('.F90')
        preprocessed = os.path.join(build_path, file_name[:-len('.F90')] + '.p.F90') if is_fortran else None
        raw_path = os.path.join(src_subdir, file_name)
        start_file = preprocessed if is_fortran and preprocessed and os.path.exists(preprocessed) else raw_path
        file_stack = [start_file]

        file_identifier = _identifier_for(start_file)
        rescan = True
        if update_time is not None and file_identifier in parameters_per_file:
            rescan = False
            for f in parameters_per_file[file_identifier].get('files', []):
                try:
                    if os.stat(f).st_mtime > update_time:
                        rescan = True
                        break
                except OSError:
                    pass
        if not rescan:
            continue

        parameters_per_file.pop(file_identifier, None)
        entry = parameters_per_file.setdefault(file_identifier, {})
        entry['files'] = [start_file]
        entry.setdefault('parameter', [])

        # Pick up parameters from the .p parameter file written by the
        # FunctionClass source-tree processor, if present.
        if is_fortran:
            parameter_file = os.path.join(
                build_path, file_name[:-len('.F90')] + '.p'
            )
            if os.path.exists(parameter_file):
                entry['parameter'].extend(_read_lines(parameter_file))

        while file_stack:
            current = file_stack.pop(0)
            for line in _read_lines(current):
                m = INCLUDE_LINE.match(line)
                if not m:
                    continue
                included_name = m.group(1)
                if included_name in INCLUDE_FILES_EXCLUDED:
                    continue
                included_path = os.path.join(build_path, included_name)
                file_stack.append(included_path)
                entry['files'].append(included_path)

            for directive in extract_directives(current, 'inputParameter'):
                name = directive.get('name')
                if name is not None:
                    entry['parameter'].append(name)

            for directive in extract_directives(current, 'objectBuilder'):
                parameter_name = directive.get('parameterName')
                class_name = directive.get('class')
                if parameter_name is not None and parameter_name != class_name:
                    entry['parameter'].append(parameter_name)

    all_parameters = []
    for record in parameters_per_file.values():
        all_parameters.extend(record.get('parameter', []))
    unique_parameters = sorted(set(all_parameters))

    output_file_name = os.path.join(
        build_path, re.sub(r'\.(exe|o)$', '.parameters.F90', target_name)
    )
    with open(output_file_name, 'w', encoding='utf-8') as output:
        output.write('subroutine knownParameterNames(names)\n')
        output.write('  use ISO_Varying_String\n')
        output.write('  implicit none\n')
        output.write('  type(varying_string), dimension(:), allocatable, intent(inout) :: names \n')
        output.write('  allocate(names(' + str(len(unique_parameters)) + '))\n')
        for index, parameter in enumerate(unique_parameters, start=1):
            output.write("  names(" + str(index) + ")='" + parameter + "'\n")
        output.write('end subroutine knownParameterNames\n')

    with open(blob_path, 'wb') as handle:
        pickle.dump(parameters_per_file, handle)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('Usage: parameterDependencies.py <sourceDirectory> <target>')
    main(sys.argv[1], sys.argv[2])
