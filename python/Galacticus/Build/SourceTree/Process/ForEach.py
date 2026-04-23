# Processes `forEach` directives: iterates over every element of a named
# array of arbitrary rank, expanding `{index}` / `{{index}}` / `%index%`
# placeholders in the directive's inline content.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ForEach.pm

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import get_declaration, add_declarations
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses


def _rank_from_declaration(parent, variable):
    """Return the array rank of `variable` declared under `parent`, 0 if scalar.

    Mirrors the `dimension(...)` inspection in ForEach.pm:27-33.
    """
    declaration = get_declaration(parent, variable)
    for attr in declaration.get('attributes') or []:
        if attr.startswith('dimension'):
            return attr.count(',') + 1
    return 0


def process_for_each(tree, options):
    """Mirrors Process_ForEach() from ForEach.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'forEach':
            continue
        directive = node.get('directive') or {}
        if directive.get('processed'):
            continue
        directive['processed'] = True

        variable = directive['variable']
        rank = _rank_from_declaration(node['parent'], variable)

        # Index variables: foreach__1, foreach__2, ..., foreach__rank.
        if rank > 0:
            add_declarations(node['parent'], [{
                'intrinsic':  'integer',
                'type':       None,
                'openMP':     False,
                'variables':  [f'foreach__{i}' for i in range(1, rank + 1)],
                'attributes': [],
            }])

        # Placeholder expansion values.  Mirrors ForEach.pm:51-53.
        if rank == 0:
            indexes_format = "'(a1)'"
            indexes        = "'.'"
            indexer        = ""
        else:
            indexes_format = "indexesFormat__"
            indexes        = ",".join(f"foreach__{i}" for i in range(1, rank + 1))
            indexer        = "(" + ",".join(f"foreach__{i}" for i in range(1, rank + 1)) + ")"

        iterator  = "! Auto-generated iteration over elements\n"
        for i in range(1, rank + 1):
            iterator += f"do foreach__{i}=1,size({variable},dim={i})\n"

        need_format = False
        content = directive.get('content', '') or ''
        for line in content.splitlines(keepends=True):
            if '%index%' in line:
                need_format = True
            line = line.replace('%index%',     indexes_format)
            line = line.replace('{{index}}',   indexes)
            line = line.replace('{index}',     indexer)
            iterator += line

        for _ in range(rank):
            iterator += "end do\n"

        if need_format and rank > 0:
            add_declarations(node['parent'], [
                {
                    'intrinsic':  'character',
                    'type':       f"len=5+{rank}*(7+int(log10(dble(huge(0_c_size_t)))))",
                    'openMP':     False,
                    'variables':  ['indexesFormat__'],
                    'attributes': [],
                },
                {
                    'intrinsic':  'character',
                    'type':       'len=12',
                    'openMP':     False,
                    'variables':  ['indexesFormatMeta_'],
                    'attributes': [],
                },
                {
                    'intrinsic':  'character',
                    'type':       f"len=5+{rank}*(6+int(log10(dble(huge(0_c_size_t)))))",
                    'openMP':     False,
                    'variables':  ['indexesFormatMeta__'],
                    'attributes': [],
                },
            ])
            add_uses(node['parent'], {
                'moduleUse': {
                    'ISO_C_Binding': {
                        'intrinsic': True,
                        'only':      {'c_size_t': True},
                    },
                },
                'moduleOrder': ['ISO_C_Binding'],
            })

            # Build the format-string setup block (ForEach.pm:111-118).
            fmt  = (
                "write (indexesFormatMeta_,'(\"(\"\"i\"\",i\",i2.2,\".\",i2.2,\")\")') "
                "1+int(log10(dble(huge(0_c_size_t)))),"
                "1+int(log10(dble(huge(0_c_size_t)))))\n"
            )
            fmt += "indexesFormat__='(\"[\",'\n"
            for i in range(1, rank + 1):
                fmt += (
                    f"write (indexesFormatMeta__,indexesFormatMeta_) "
                    f"1+int(log10(dble(size({variable},dim={i}))))\n"
                )
                suffix = "//',\",\",'" if i < rank else ""
                fmt += f"indexesFormat__=trim(indexesFormat__)//trim(indexesFormatMeta__){suffix}\n"
            fmt += "indexesFormat__=trim(indexesFormat__)//',\"]\")'\n"

            iterator = fmt + iterator

        iterator = (
            "! Auto-generated iteration over elements\n"
            + iterator
            + "! End auto-generated iteration over elements\n"
        )

        insert_after_node(node, [{
            'type':       'code',
            'content':    iterator,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])


register_process('forEach', process_for_each, before=['generics'])
