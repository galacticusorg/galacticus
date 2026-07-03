#!/usr/bin/env python3
import os
import re
import subprocess
import sys
import shlex

# Make a static-linked copy of an executable by replacing dynamic library links
# with their static equivalents. macOS-specific (uses otool).
# Andrew Benson (ported to Python 2026)

arguments = list(sys.argv[1:])

# Find the output executable name from the -o flag.
executable = None
for i, arg in enumerate(arguments):
    if arg == '-o' and i + 1 < len(arguments):
        executable = arguments[i + 1]
        break

if executable is None:
    print("Error: unable to determine executable", file=sys.stderr)
    sys.exit(1)

print(f"Executable is '{executable}'")

# Shell operators that begin the plumbing appended to the link recipe in the
# build log this line was grepped from: output redirection, pipes, command
# terminators, and a trailing line-continuation backslash. None of it is part
# of the link command this script must re-run (it reconstructs its own
# invocation and lets its output flow to the console), so parsing stops at the
# first such token and the remainder is discarded. This handles both the older
# recipe shape ('... 2>&1 | ./scripts/build/postprocessLinker.py') and the
# current one ('... > <stem>.link.diag 2>&1; \'), which redirects diagnostics
# to a file for exit-status capture.
_SHELL_PLUMBING = {'2>&1', '&>', '1>', '2>', '|', '||', '&&', ';', '&', '\\'}

# Assemble the link command, expanding backtick-quoted sub-expressions.
compile_parts = []
i = 0
while i < len(arguments):
    arg = arguments[i]
    if arg in _SHELL_PLUMBING or arg.startswith('>') or arg.startswith('<'):
        break
    # Expand backtick-quoted shell expressions.
    if arg.startswith('`'):
        to_expand = arg
        while not to_expand.endswith('`'):
            if i + 1 >= len(arguments):
                print(
                    "Error: unterminated backtick-quoted expression in arguments",
                    file=sys.stderr,
                )
                sys.exit(1)
            i += 1
            to_expand += ' ' + arguments[i]
        to_expand = to_expand.strip('`')
        expanded = subprocess.run(
            to_expand, shell=True, capture_output=True, text=True
        ).stdout.replace('\n', ' ')
        compile_parts.append(expanded.strip())
    else:
        compile_parts.append(arg)
    i += 1

compile_command = ' '.join(compile_parts)

# Inspect the linked dynamic libraries.
is_gcc      = False
is_gfortran = False
is_gpp      = False

# Compiler runtime libraries that the compiler driver links implicitly and for which there is no
# dedicated '-static-lib*' option (unlike libgfortran, libgcc, and libstdc++). When these are linked
# statically (by naming their '.a' archive) any coexisting dynamic library must be temporarily moved
# aside so the linker uses the static archive rather than preferring the dynamic one.
hide_dylib_libraries = ('quadmath', 'gomp')

# Full paths of dynamic libraries to temporarily move aside before relinking. A set, so that each is
# recorded (and therefore moved) at most once even if encountered repeatedly.
dylibs_to_hide = set()

def register_dylibs_to_hide(static_lib_path, library_name):
    """Record the dynamic versions of `library_name` that sit alongside its static archive, so they
    can be temporarily moved aside to force the linker to use the static archive."""
    directory = os.path.dirname(static_lib_path)
    if not directory or not os.path.isdir(directory):
        return
    pattern = re.compile(r'^lib' + re.escape(library_name) + r'(\.\d+)*\.dylib$')
    for fname in os.listdir(directory):
        if pattern.match(fname):
            dylibs_to_hide.add(os.path.join(directory, fname))

try:
    otool_out = subprocess.run(
        ['otool', '-L', executable], capture_output=True, text=True
    ).stdout
except FileNotFoundError:
    print("Error: otool not found (this script is macOS-specific)", file=sys.stderr)
    sys.exit(1)

print("Parsing otool output:")
for line in otool_out.splitlines():
    columns = line.split()
    if not columns:
        continue
    lib_path = columns[0]    
    if not (re.search(r'\.dylib$', lib_path) or re.search(r'\.so[0-9.]+$', lib_path)):
        continue

    dynamic_name = lib_path
    m = re.search(r'^.*/lib([a-zA-Z0-9_\-\+]+)\.', dynamic_name)
    library_name          = m.group(1) if m else dynamic_name
    library_name_original = library_name

    if re.search(r'\.dylib$', dynamic_name):
        static_name = re.sub(r'(\.\d+)?\.dylib$', '.a', dynamic_name)
    else:
        static_name = re.sub(r'\.so[0-9.]+$', '.a', dynamic_name)

    if library_name == 'qhull_r':
        library_name = 'qhullstatic_r'

    print(f"Line: '{line}'")
    print(f"   Original path was '{dynamic_name}'")
    print(f"   Looking for static library for '{library_name}' (originally '{library_name_original}')")
    print(f"   Static name is '{static_name}'")

    if re.match(r'^gcc', library_name):
        is_gcc = True
        print(" -> Can use gcc compiler option")
    elif re.match(r'^gfortran', library_name):
        is_gfortran = True
        print(" -> Can use gfortran compiler option")
    elif re.match(r'^stdc\+\+', library_name):
        is_gpp = True
        print(" -> Can use g++ compiler option")
    elif os.path.exists(static_name):
        print(f" -> Found static library at '{static_name}'")
        escaped = re.escape(library_name_original)
        if re.search(r'-l' + escaped + r'(\s|$)', compile_command):
            compile_command = re.sub(r'-l' + escaped + r'(\s|$)', static_name + ' ', compile_command)
        else:
            compile_command += ' ' + static_name
            if library_name in hide_dylib_libraries:
                register_dylibs_to_hide(static_name, library_name)
    else:
        found      = False
        candidates = []
        # First, ask the compiler (the first token of the link command) where its own copy of the
        # static library lives. This locates runtime libraries such as libgomp.a and libquadmath.a
        # that reside in the compiler's installation directory (e.g. /opt/gcc-16/lib/...) rather
        # than in a standard system location, without having to hard-code that path.
        if compile_parts:
            try:
                printed = subprocess.run(
                    [compile_parts[0], f"-print-file-name=lib{library_name}.a"],
                    capture_output=True, text=True
                ).stdout.strip()
            except (FileNotFoundError, OSError):
                printed = ''
            # '-print-file-name' echoes back the bare library name if it is not found, so only
            # accept the result if it is an absolute path.
            if printed and os.path.isabs(printed):
                candidates.append(printed)
        # Next, search standard system locations and any directories on LD_LIBRARY_PATH.
        locations = ['/usr/local/lib']
        ld_path   = os.environ.get('LD_LIBRARY_PATH', '')
        if ld_path:
            locations.extend(ld_path.split(':'))
        candidates.extend(os.path.join(loc, f"lib{library_name}.a") for loc in locations)
        for candidate in candidates:
            if os.path.exists(candidate):
                print(f" -> Found static library at '{candidate}'")
                escaped = re.escape(library_name_original)
                if re.search(r'-l' + escaped + '( |$)', compile_command):
                    compile_command = re.sub(r'-l' + escaped + '( |$)', candidate + ' ', compile_command)
                else:
                    compile_command += ' ' + candidate
                    if library_name in hide_dylib_libraries:
                        register_dylibs_to_hide(candidate, library_name)
                found = True
                break
        if not found:
            print(" -> No static library found")

# Add static compiler flags.
if is_gfortran:
    compile_command += ' -static-libgfortran'
if is_gcc:
    compile_command += ' -static-libgcc'
if is_gpp:
    compile_command += ' -static-libstdc++'

# Temporarily move aside any dynamic libraries that would otherwise be preferred over the static
# archives being linked, forcing the linker to use the static versions. They are restored below.
if dylibs_to_hide:
    print("Must move dylibs temporarily (requires sudo):")
    for dylib in sorted(dylibs_to_hide):
        print(f"   {dylib} -> {dylib}~")
    mv_cmd = "sudo -- sh -c '" + '; '.join(
        "mv " + shlex.quote(dylib) + " " + shlex.quote(dylib + "~") for dylib in sorted(dylibs_to_hide)
    ) + "'"
    subprocess.run(mv_cmd, shell=True)

print(f"Relinking with: {compile_command}")
status = subprocess.run(compile_command, shell=True)

# Restore any temporarily moved dynamic libraries.
if dylibs_to_hide:
    print("Must restore temporarily moved dylibs (requires sudo):")
    mv_cmd = "sudo -- sh -c '" + '; '.join(
        "mv " + shlex.quote(dylib + "~") + " " + shlex.quote(dylib) for dylib in sorted(dylibs_to_hide)
    ) + "'"
    subprocess.run(mv_cmd, shell=True)

# Exit with status.
sys.exit(status.returncode)
