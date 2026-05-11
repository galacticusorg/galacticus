"""Processes `constant` directives: emits a Fortran `parameter` declaration
whose value is either given inline, extracted from a GSL symbol by
compiling-and-running a tiny C program, or similarly extracted from a
kernel header.  In documentation builds, also writes a per-file
`<constants>` XML manifest.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/SourceTree/Process/Constants.pm
"""

import os
import re
import tempfile
import subprocess
import xml.etree.ElementTree as ET


from Galacticus.Build.SourceTree         import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process import register_process


def _compile_and_run(c_source):
    """Compile `c_source` with $CCOMPILER $CFLAGS and return its stdout.

    Used for both the GSL and kernel extraction paths (Constants.pm:56-60,
    Constants.pm:127-132).  Raises on compilation failure.
    """
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError("process_constant: BUILDPATH is not set")
    compiler = os.environ.get('CCOMPILER')
    if not compiler:
        raise RuntimeError("process_constant: CCOMPILER is not set")
    cflags = os.environ.get('CFLAGS', '').split()

    src_fd, src_path = tempfile.mkstemp(prefix='temp', suffix='.c', dir=build_path)
    exe_fd, exe_path = tempfile.mkstemp(prefix='temp', suffix='.x', dir=build_path)
    os.close(src_fd)
    os.close(exe_fd)
    try:
        with open(src_path, 'w') as fh:
            fh.write(c_source)
        compile_res = subprocess.run(
            [compiler, '-o', exe_path, src_path, *cflags],
            capture_output=True,
        )
        if compile_res.returncode != 0:
            raise RuntimeError(
                "process_constant: failed to compile:\n"
                + compile_res.stderr.decode('utf-8', errors='replace'))
        run_res = subprocess.run([exe_path], capture_output=True)
        return run_res.stdout.decode('utf-8', errors='replace')
    finally:
        for p in (src_path, exe_path):
            try:
                os.unlink(p)
            except OSError:
                pass


def _emit_gsl_enum(directive):
    """Build `enum, bind(c) ... end enum` from a GSL-header enum.

    Mirrors Constants.pm:42-70.
    """
    members = [m.strip() for m in re.split(r'\s*,\s*', directive['members']) if m.strip()]
    c_src = (
        "#include <stdio.h>\n"
        "#include <float.h>\n"
        f"#include <gsl/{directive['gslHeader']}.h>\n"
        "int main () {\n"
        f" printf(\"{' '.join(['%i'] * len(members))}\", {directive['members']});\n"
        "}\n"
    )
    output = _compile_and_run(c_src)
    values = output.split()
    definitions = [f"{m}={v}" for m, v in zip(members, values)]
    return (
        "enum, bind(c)\n"
        " enumerator :: " + ", ".join(definitions) + "\n"
        "end enum\n"
    )


def _emit_gsl_scalar(directive, type_str):
    """Build a `type, parameter, public :: var = value` line from a GSL scalar.

    Mirrors Constants.pm:72-110.
    """
    c_src  = "#include <stdio.h>\n"
    c_src += "#include <float.h>\n"
    c_src += f"#include <gsl/{directive['gslHeader']}.h>\n"
    if type_str == 'double precision':
        c_src += (
            "#ifdef LDBL_DECIMAL_DIG\n"
            "  #define OP_LDBL_Digs (LDBL_DECIMAL_DIG)\n"
            "#else\n"
            "  #ifdef DECIMAL_DIG\n"
            "    #define OP_LDBL_Digs (DECIMAL_DIG)\n"
            "  #else\n"
            "    #define OP_LDBL_Digs (LDBL_DIG + 3)\n"
            "  #endif\n"
            "#endif\n"
            "int main () {\n"
            f" printf(\"%.*e\", OP_LDBL_Digs - 1, {directive['gslSymbol']});\n"
            "}\n"
        )
    elif type_str == 'integer':
        c_src += (
            "int main () {\n"
            f" printf(\"%i\", {directive['gslSymbol']});\n"
            "}\n"
        )
    else:
        raise RuntimeError(f"process_constant: unsupported gslSymbol type '{type_str}'")

    value = _compile_and_run(c_src)
    value = value.replace('e', 'd')
    directive['value'] = value
    return f"{type_str}, parameter, public :: {directive['variable']}={value}\n"


def _emit_kernel_scalar(directive, type_str):
    """Build a `parameter` line from a kernel-header integer symbol.

    Mirrors Constants.pm:112-135.
    """
    if type_str != 'integer':
        raise RuntimeError(
            "process_constant: kernelSymbol constants must have type='integer'")
    c_src = (
        "#include <stdio.h>\n"
        f"#include <{directive['kernelHeader']}.h>\n"
        "int main () {\n"
        f" printf(\"%i\", {directive['kernelSymbol']});\n"
        "}\n"
    )
    value = _compile_and_run(c_src)
    directive['value'] = value
    return f"{type_str}, parameter, public :: {directive['variable']}={value}\n"


def _emit_inline(directive, type_str):
    """Build a `parameter` line from a directive-supplied value."""
    return f"{type_str}, parameter, public :: {directive['variable']}={directive['value']}\n"


def process_constant(tree, options):
    """Mirrors Process_Constant() from Constants.pm."""
    constants_orphaned = []   # directives awaiting a module-name tag
    all_constants      = []
    file_name          = None

    for node in walk_tree(tree):
        ntype = node.get('type')
        if ntype == 'file':
            file_name = node.get('name')
            continue
        if ntype == 'module':
            module_name = node.get('name')
            if constants_orphaned:
                for c in constants_orphaned:
                    c['module'] = module_name
                all_constants.extend(constants_orphaned)
                constants_orphaned = []
            continue
        if ntype != 'constant':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue

        type_str = directive.get('type', 'double precision')
        if 'gslSymbol' in directive:
            if type_str == 'enum':
                code = _emit_gsl_enum(directive)
            else:
                code = _emit_gsl_scalar(directive, type_str)
        elif 'kernelSymbol' in directive:
            code = _emit_kernel_scalar(directive, type_str)
        else:
            code = _emit_inline(directive, type_str)

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])
        constants_orphaned.append(directive)
        directive['processed'] = True

    # Tag final file name.  Matches Constants.pm:156-159.
    for c in all_constants:
        c['fileName'] = file_name

    # Documentation mode: dump collected constants to a per-file XML manifest.
    if (all_constants
            and os.environ.get('GALACTICUS_BUILD_DOCS') == 'yes'
            and file_name is not None):
        build_path = os.environ.get('BUILDPATH')
        if build_path:
            out_name = re.sub(r'\.F90$', '.constants.xml', file_name)
            out_path = os.path.join(build_path, out_name)
            root_el = ET.Element('constants')
            for c in all_constants:
                const_el = ET.SubElement(root_el, 'constant')
                for k, v in c.items():
                    if k == 'processed':
                        continue
                    if not isinstance(v, (str, int, float)):
                        continue
                    ET.SubElement(const_el, k).text = str(v)
            ET.ElementTree(root_el).write(out_path, encoding='utf-8', xml_declaration=True)


register_process('constant', process_constant)
