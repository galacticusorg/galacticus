# Processes `dependenciesInitialize` directives: reads the Galacticus
# dependencies manifest and emits initialization calls populating the
# dependency registry at runtime.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/Dependencies.pm

import re
import os


from Galacticus.Build.SourceTree         import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process import register_process


def _read_dependencies():
    """Read `$GALACTICUS_EXEC_PATH/aux/dependencies.yml` into {name: version}.

    Mirrors the reader at Dependencies.pm:23-32.  Any line that does not match
    `<name>: <version>` is a fatal error.
    """
    exec_path = os.environ.get('GALACTICUS_EXEC_PATH')
    if not exec_path:
        raise RuntimeError(
            "process_dependencies: GALACTICUS_EXEC_PATH environment variable is not set")
    path = os.path.join(exec_path, 'aux', 'dependencies.yml')
    deps = {}
    with open(path, 'r') as fh:
        for line in fh:
            m = re.match(r'^(.+):\s+([0-9\.]+)', line)
            if not m:
                raise RuntimeError(
                    f"process_dependencies: cannot parse dependency file line:\n{line}")
            deps[m.group(1)] = m.group(2)
    return deps


def process_dependencies(tree, options):
    """Mirrors Process_Dependencies() from Dependencies.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'dependenciesInitialize':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        deps = _read_dependencies()
        code = "\n".join(
            f"call dependencies_%set(var_str('{name}'),var_str('{deps[name]}'))"
            for name in sorted(deps)
        ) + "\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.Dependencies.process_dependencies()',
            'line':       1,
        }])


register_process('dependencies', process_dependencies)
