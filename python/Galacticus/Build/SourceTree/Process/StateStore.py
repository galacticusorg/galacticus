# Processes `stateStore` and `stateRestore` directives: emits per-variable
# serialize/deserialize calls to/from `stateFile`, and for stateRestore adds
# a `wasAllocated_` tracker plus `use :: Error, only : Error_Report` so the
# generated code can report an inconsistent state.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/StateStore.pm



from Galacticus.Build.SourceTree                          import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process                  import register_process
from Galacticus.Build.SourceTree.Parse.Declarations       import add_declarations
from Galacticus.Build.SourceTree.Parse.ModuleUses         import add_uses
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


def _emit_code_after(node, content):
    insert_after_node(node, [{
        'type':       'code',
        'content':    content,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     node.get('source', 'unknown'),
        'line':       node.get('line', 0),
    }])


def process_state_store(tree, options):
    """Mirrors Process_StateStore() from StateStore.pm, handling BOTH
    `stateStore` and `stateRestore` directive types in a single walk.
    """
    for node in walk_tree(tree):
        ntype = node.get('type')
        if ntype not in ('stateStore', 'stateRestore'):
            continue
        directive = node.setdefault('directive', {})

        if ntype == 'stateStore' and not directive.get('processed'):
            directive['processed'] = True
            variables = (directive.get('variables') or '').split()
            code = ''
            for v in variables:
                code += f"write (stateFile) associated({v})\n"
                code += (
                    f"if (associated({v})) call {v}"
                    "%stateStore(stateFile,gslStateFile,stateOperationID)\n"
                )
            _emit_code_after(node, code)

        elif ntype == 'stateRestore' and not directive.get('processed'):
            directive['processed'] = True
            variables = (directive.get('variables') or '').split()
            loc_expr = location(node, node.get('line', 0))
            code = ''
            for v in variables:
                code += "read (stateFile) wasAllocated_\n"
                code += "if (wasAllocated_) then\n"
                code += (
                    f"   if (.not.associated({v})) call Error_Report("
                    f"\"'{v}' was stored, but is now not allocated\"//"
                    f"{loc_expr})\n"
                )
                code += (
                    f"   call {v}%stateRestore(stateFile,gslStateFile,"
                    "stateOperationID)\n"
                )
                code += "end if\n"
            _emit_code_after(node, code)

            add_declarations(node['parent'], [{
                'intrinsic':     'logical',
                'type':          None,
                'openMP':        False,
                'attributes':    [],
                'variables':     ['wasAllocated_'],
                'variableNames': ['wasAllocated_'],
            }])

            add_uses(node['parent'], {
                'moduleUse': {
                    'Error': {
                        'intrinsic': False,
                        'only':      {'Error_Report': True},
                    },
                },
                'moduleOrder': ['Error'],
            })


register_process('stateStore', process_state_store)
