# Processes `parameterMigration` directives: emits a `commitHash(i)="..."`
# initialization block from `scripts/aux/migrations.xml` and adds the
# backing character array declaration.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ParameterMigration.pm

import os
import sys
import xml.etree.ElementTree as ET


from List.ExtraUtils                                import as_array
from XML.Utils                                      import xml_to_dict
from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import add_declarations


def _read_migrations():
    """Return the ordered list of `<migration>` commit hashes from
    `$GALACTICUS_EXEC_PATH/scripts/aux/migrations.xml`.
    """
    exec_path = os.environ.get('GALACTICUS_EXEC_PATH')
    if not exec_path:
        raise RuntimeError(
            "process_parameter_migration: GALACTICUS_EXEC_PATH is not set")
    path = os.path.join(exec_path, 'scripts', 'aux', 'migrations.xml')
    root = ET.parse(path).getroot()
    migrations_dict = xml_to_dict(root)
    migrations = list(as_array(migrations_dict.get('migration')))
    return [m.get('commit') for m in migrations]


def process_parameter_migration(tree, options):
    """Mirrors Process_ParameterMigration() from ParameterMigration.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'parameterMigration':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        commit_hashes = _read_migrations()
        code = ''.join(
            f'commitHash({i + 1})="{h}"//c_null_char\n'
            for i, h in enumerate(commit_hashes)
        )

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.ParameterMigration.process_parameter_migration()',
            'line':       1,
        }])

        add_declarations(node['parent'], [{
            'intrinsic':    'character',
            'type':         'len=41',
            'openMP':       False,
            'attributes':   [f'dimension({len(commit_hashes)})'],
            'variables':    ['commitHash'],
            'variableNames': ['commitHash'],
            'preprocessor': 'GIT2AVAIL',
        }])


register_process('parameterMigration', process_parameter_migration)
