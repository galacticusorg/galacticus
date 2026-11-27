#!/usr/bin/env python3
"""Consolidate the procedure-level USE statements of a preprocessed Fortran module or submodule into its specification
part.

Galacticus style places USE statements in each procedure. gfortran reads a module's file again for every USE statement,
and creates a new copy of every symbol that the USEd names reach (for a class module, that is the closure of all of the
types and type-bound procedures reachable from it - often tens of megabytes of module file). The copies persist until the
whole file has been compiled (GCC PR 117474), and the time taken to build the derived types of the file grows as the square of the number
of copies. Consolidating the USE statements reduced the time to compile some files from minutes to seconds, and their
memory use from 9-15 GB to under 1 GB. It also shrinks the module files written by a module, as each procedure's copies
of the symbols were written to the module file.

A name is moved only when doing so cannot change what any name means:
 - every region (the unit's specification part, and each procedure) which mentions the name either imports the same
   entity under that name itself, or is contained in a procedure which does;
 - the import is unconditional (not inside a preprocessor conditional, nor behind an OpenMP "!$" sentinel), lies outside
   any interface body, and comes from a USE statement with an ONLY list;
 - it is not "assignment(=)" or an intrinsic operator, whose visibility can change the meaning of a statement which never
   names it.
Names which are already accessible, as the same entity, from the unit's own specification part, or from the parent of a
submodule, are dropped from the procedures rather than moved, as host association provides them. (gfortran rejects a USE
statement in the specification part of a submodule which imports a generic name, or a name from ISO_C_Binding, which the
parent module already provides - ifx accepts it - so such names cannot be moved there.) In a module without a default
PRIVATE statement, the moved names are declared PRIVATE so that the interface of the module is unchanged.

The physical line count is preserved (moved statements are joined to the unit statement with ";", and lines removed are
left blank), so the line map (".lmap") of the preprocessed file remains valid for the output. If anything about the
structure of the file is not understood, it is written out unchanged.

Usage: consolidateUses.py <input.p.F90> <output.F90> [--moduleDependencies <Makefile_Module_Dependencies>] [--report]
"""
import argparse
import os
import re
import sys

# Andrew Benson (20-November-2026)

# <workaround type="gfortran" PR="117474" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=117474"/>
# This whole script is the workaround: see the module docstring. gfortran re-reads a module, and copies the symbols it
# reaches, for every USE statement (GCC PR 117474 reports the same cost for USE statements in interface bodies).

PREFIXES = r'(?:(?:pure|impure|elemental|recursive|non_recursive|module)\s+)*'
PARENS = r'\((?:[^()]|\([^()]*\))*\)'
TYPESPEC = (r'(?:(?:integer|real|logical|complex|character|double\s*precision|double\s*complex|type\s*' + PARENS
            + r'|class\s*' + PARENS + r')(?:\s*' + PARENS + r'|\s*\*\s*\d+)?\s+)?')
RE_PROC_START = re.compile(r'^' + PREFIXES + TYPESPEC + PREFIXES + r'(subroutine|function)\s+(\w+)')
RE_MODPROC_START = re.compile(r'^module\s+procedure\s+(\w+)\s*$')
RE_PROC_END = re.compile(r'^end\s*(subroutine|function|procedure)\b\s*(\w*)')
RE_IFACE_START = re.compile(r'^(abstract\s+)?interface\b')
RE_IFACE_END = re.compile(r'^end\s*interface\b')
RE_TYPE_START = re.compile(r'^type\s*(,[^:]*)?::\s*\w+|^type\s+(?!is\b)\w+\s*(\(.*\))?$')
RE_TYPE_END = re.compile(r'^end\s*type\b')
RE_UNIT = re.compile(r'^(module|submodule)\b\s*(\(\s*([\w:]+)\s*\))?\s*(\w+)\s*$')
RE_OTHER_UNIT = re.compile(r'^(program|block\s*data)\b')
RE_USE = re.compile(r'^use\b\s*(,\s*(intrinsic|non_intrinsic)\s*)?(::)?\s*(\w+)\s*(,\s*(only\s*:)?\s*(.*))?$')
RE_DECLARATION = re.compile(r'^(integer|real|logical|complex|character|double\s*precision|double\s*complex|type\s*\(|class\s*\('
                            r'|procedure\s*\()')
RE_IDENT = re.compile(r'(?<![%\w])[a-z_]\w*')
RE_DOTOP = re.compile(r'\.[a-z]+\.')
INTRINSIC_DOTOPS = {'.and.', '.or.', '.not.', '.eqv.', '.neqv.', '.eq.', '.ne.', '.lt.', '.le.', '.gt.', '.ge.', '.true.',
                    '.false.'}


class Unsupported(Exception):
    """Raised when the structure of a file is not understood, in which case it is left unchanged."""


def splitCode(line):
    """Return the code part of a physical line (with comments and the contents of strings removed), and whether it ends
    with a continuation mark."""
    out = []
    quote = None
    for character in line:
        if quote:
            if character == quote:
                quote = None
                out.append(character)
            continue
        if character in "'\"":
            quote = character
            out.append(character)
        elif character == '!':
            break
        else:
            out.append(character)
    code = ''.join(out).rstrip()
    return code, code.endswith('&')


def splitSemicolons(text):
    """Split a statement line at ";" separators."""
    pieces = [piece.strip() for piece in text.split(';')]
    return [piece for piece in pieces if piece] or ['']


def statements(lines):
    """Group physical lines into statements, each a tuple (first line, last line, kind, text). Kinds are "code", "sentinel"
    (OpenMP conditional compilation lines, and OpenMP directives, whose clauses may name variables), "cpp", and "blank"
    (which includes comments)."""
    result = []
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        if stripped.startswith('#'):
            if stripped.rstrip().endswith('\\'):
                raise Unsupported('continued preprocessor line')
            result.append((i, i, 'cpp', stripped.strip()))
            i += 1
            continue
        isOpenMP = stripped[:5].lower() == '!$omp'
        isSentinel = stripped[:2] == '!$' and not isOpenMP
        if stripped.startswith('!') and not isSentinel:
            if isOpenMP:
                result.append((i, i, 'sentinel', stripped[5:]))
            else:
                result.append((i, i, 'blank', ''))
            i += 1
            continue
        code, continued = splitCode(stripped[2:] if isSentinel else lines[i])
        first = i
        parts = [code[:-1] if continued else code]
        while continued:
            i += 1
            if i >= len(lines):
                raise Unsupported('unterminated continuation')
            following = lines[i].lstrip()
            if isSentinel:
                if not following.startswith('!$'):
                    raise Unsupported('continuation of an OpenMP conditional compilation line')
                following = following[2:].lstrip()
            elif following.startswith('#'):
                raise Unsupported('preprocessor line inside a continued statement')
            elif following.startswith('!') or not following.strip():
                # Comment lines may appear between continuation lines.
                continue
            code, continued = splitCode(following)
            code = code[1:] if code.startswith('&') else code
            parts.append(code[:-1] if continued else code)
        text = ' '.join(part.strip() for part in parts).strip()
        for piece in splitSemicolons(text):
            result.append((first, i, 'sentinel' if isSentinel else ('code' if piece else 'blank'), piece))
        i += 1
    return result


def parseOnly(text):
    """Split an ONLY list into (local name, original name) pairs, in lower case and without spaces."""
    items, depth, current = [], 0, []
    for character in text:
        if character == '(':
            depth += 1
        elif character == ')':
            depth -= 1
        if character == ',' and depth == 0:
            items.append(''.join(current))
            current = []
        else:
            current.append(character)
    items.append(''.join(current))
    pairs = []
    for item in items:
        item = re.sub(r'\s+', '', item.lower())
        if not item:
            continue
        local, original = item.split('=>', 1) if '=>' in item else (item, item)
        pairs.append((local, original))
    return pairs


def isMovable(local):
    """Defined assignment and intrinsic operators are never moved: their visibility can change the meaning of statements
    which do not name them."""
    if local.startswith('assignment('):
        return False
    if local.startswith('operator('):
        operator = local[len('operator('):-1]
        return bool(re.fullmatch(r'\.[a-z]+\.', operator)) and operator not in INTRINSIC_DOTOPS
    return True


def mentionKey(local):
    """The token by which a use of a local name is recognized in code."""
    return local[len('operator('):-1] if local.startswith('operator(') else local


def mentions(text):
    """The set of names (and defined operators) appearing in a statement."""
    text = text.lower()
    return set(RE_IDENT.findall(text)) | set(RE_DOTOP.findall(text))


def splitDeclaration(low):
    """Split a type declaration statement at its "::" into the type specification and attributes, and the list of
    entities declared. Returns None for other statements."""
    if not RE_DECLARATION.match(low):
        return None
    depth = 0
    for i, character in enumerate(low):
        if character == '(':
            depth += 1
        elif character == ')':
            depth -= 1
        elif depth == 0 and low.startswith('::', i):
            return low[:i], low[i + 2:]
    return None


def declaredNames(entities):
    """The names declared by the entity list of a type declaration statement (the part after its "::")."""
    names, depth, current = set(), 0, []
    for character in entities + ',':
        if character in '([':
            depth += 1
        elif character in ')]':
            depth -= 1
        if character == ',' and depth == 0:
            match = re.match(r'\s*([a-z_]\w*)', ''.join(current))
            if match:
                names.add(match.group(1))
            current = []
        else:
            current.append(character)
    return names


def specificationImports(lines):
    """The imports made unconditionally by the specification part of a module: a set of (module, local, original) triples,
    with (module, '*', '*') for a USE statement without an ONLY list (and without renames)."""
    found = set()
    depth = 0
    for first, last, kind, text in statements(lines):
        low = text.lower()
        if kind == 'cpp':
            if re.match(r'#\s*if', low):
                depth += 1
            elif re.match(r'#\s*endif', low):
                depth -= 1
            continue
        if kind != 'code':
            continue
        if re.match(r'^contains\b', low):
            break
        match = RE_USE.match(low)
        if not match or depth != 0:
            continue
        if match.group(6):
            for local, original in parseOnly(match.group(7)):
                found.add((match.group(4), local, original))
        elif not match.group(5):
            found.add((match.group(4), '*', '*'))
    return found


def findParentSource(parent, moduleDependencies):
    """Find the source from which the parent of a submodule was compiled, using the generated module dependencies
    makefile. If the parent's USE statements were themselves consolidated, that is its consolidated source (".pc.F90"),
    whose specification part includes the names moved there; otherwise it is the preprocessed source (".p.F90")."""
    if not moduleDependencies or not os.path.exists(moduleDependencies) or ':' in parent:
        return None
    target = re.compile(r'^(\S*/)?' + re.escape(parent.lower()) + r'\.mod:\s*(\S+)\.o\s*$')
    with open(moduleDependencies, errors='replace') as dependencies:
        for line in dependencies:
            match = target.match(line)
            if match:
                for suffix in ('.pc.F90', '.p.F90'):
                    if os.path.exists(match.group(2) + suffix):
                        return match.group(2) + suffix
                return None
    return None


def analyze(lines):
    """Parse a file into its unit, its regions (the specification part, and each procedure), and its statements."""
    stmts = statements(lines)
    unit = None
    unitIndex = None
    # For the specification part, names which appear only in the type specification of a declaration (e.g. "type(a)"
    # in "type(a), pointer :: b") are recorded separately: such a name refers to an accessible entity, and is not
    # declared there.
    regions = {0: {'parent': None, 'imports': [], 'mentions': set(), 'typeMentions': set(), 'interface': False,
                   'name': None, 'unknownImports': False, 'declared': set()}}
    stack = []
    interfaceDepth = 0
    typeDepth = 0
    cppDepth = 0
    unitContains = None
    defaultPrivate = False
    for index, (first, last, kind, text) in enumerate(stmts):
        low = text.lower().strip()
        region = stack[-1] if stack else 0
        if kind == 'cpp':
            if re.match(r'#\s*if', low):
                cppDepth += 1
            elif re.match(r'#\s*endif', low):
                cppDepth -= 1
            continue
        if kind == 'blank':
            continue
        low = re.sub(r'^\d+\s+', '', low)  # Remove any statement label.
        if unit is None:
            if RE_OTHER_UNIT.match(low):
                raise Unsupported('not a module or submodule')
            match = RE_UNIT.match(low)
            if match and kind == 'code' and not low.startswith('module procedure'):
                unit = (match.group(1), match.group(3), match.group(4))
                unitIndex = index
                continue
            if kind == 'code':
                raise Unsupported('code before the unit statement')
            continue
        if re.match(r'^end\s*(sub)?module\b', low) or low == 'end':
            if stack or interfaceDepth or typeDepth:
                raise Unsupported('unit ends inside a construct')
            if any(other[2] not in ('blank',) for other in stmts[index + 1:]):
                raise Unsupported('more than one program unit')
            break
        if kind == 'code':
            if RE_IFACE_START.match(low):
                interfaceDepth += 1
            elif RE_IFACE_END.match(low):
                interfaceDepth -= 1
            elif interfaceDepth == 0 and RE_TYPE_START.match(low):
                typeDepth += 1
            elif RE_TYPE_END.match(low):
                typeDepth -= 1
            elif typeDepth == 0 and not stack and interfaceDepth == 0 and re.match(r'^contains\s*$', low):
                unitContains = index
            elif typeDepth == 0 and not stack and interfaceDepth == 0 and low == 'private':
                defaultPrivate = True
            if typeDepth == 0:
                start = RE_PROC_START.match(low)
                moduleProcedure = RE_MODPROC_START.match(low) if interfaceDepth == 0 else None
                if (start or moduleProcedure) and not low.startswith('end'):
                    if unitContains is None and interfaceDepth == 0:
                        raise Unsupported('procedure before CONTAINS')
                    identifier = len(regions)
                    regions[identifier] = {'parent': region, 'imports': [], 'mentions': set(),
                                           'interface': interfaceDepth > 0, 'unknownImports': False, 'declared': set(),
                                           'name': start.group(2) if start else moduleProcedure.group(1)}
                    stack.append(identifier)
                    region = identifier
                else:
                    end = RE_PROC_END.match(low)
                    if end:
                        if not stack:
                            raise Unsupported('unmatched end of procedure')
                        if end.group(2) and end.group(2) != regions[stack[-1]]['name']:
                            raise Unsupported(f'end of procedure "{end.group(2)}" does not match its start')
                        regions[effective(regions, stack[-1])]['mentions'] |= mentions(low)
                        stack.pop()
                        continue
        match = RE_USE.match(low)
        if match:
            entry = {'statement': index, 'module': match.group(4), 'intrinsic': match.group(2) == 'intrinsic',
                     'only': parseOnly(match.group(7)) if match.group(6) else None,
                     'renames': bool(match.group(5)) and not match.group(6),
                     'conditional': cppDepth > 0 or kind == 'sentinel', 'interface': regions[region]['interface']}
            regions[region]['imports'].append(entry)
            if entry['only'] is None:
                regions[effective(regions, region)]['unknownImports'] = True
            continue
        # Interface bodies do not see the names of their host (unless they IMPORT them): their mentions are attributed to
        # the enclosing region, which is conservative.
        target = effective(regions, region)
        declaration = splitDeclaration(low) if kind == 'code' else None
        if declaration and not regions[region]['interface']:
            regions[region]['declared'] |= declaredNames(declaration[1])
        if declaration and target == 0:
            regions[0]['typeMentions'] |= mentions(declaration[0])
            regions[0]['mentions'] |= mentions(declaration[1])
        else:
            regions[target]['mentions'] |= mentions(low)
    if unit is None or unitContains is None:
        raise Unsupported('no module or submodule with contained procedures')
    if stack or interfaceDepth or typeDepth or cppDepth:
        raise Unsupported('unbalanced constructs')
    return stmts, unit, unitIndex, unitContains, regions, defaultPrivate


def effective(regions, region):
    """The region to which the mentions of a region are attributed: interface bodies defer to their host."""
    while regions[region]['interface']:
        region = regions[region]['parent']
    return region


def consolidate(lines, parentImports=None):
    """Consolidate the USE statements of a file. Returns the new lines, the number of USE statements removed, and the
    number of names moved."""
    stmts, unit, unitIndex, unitContains, regions, defaultPrivate = analyze(lines)
    if unit[0] == 'submodule' and parentImports is None:
        raise Unsupported('the parent of the submodule was not found')
    parentImports = parentImports or set()

    def providesItself(region, key):
        module, intrinsic, local, original = key
        return any((imp['module'], imp['intrinsic']) == (module, intrinsic) and (local, original) in imp['only']
                   for imp in regions[region]['imports'] if imp['only'] is not None and not imp['interface'])

    def importsLocalName(region, local):
        return any(any(name == local for name, _ in imp['only'])
                   for imp in regions[region]['imports'] if imp['only'] is not None and not imp['interface'])

    def unitProvides(key):
        module, intrinsic, local, original = key
        for imp in regions[0]['imports']:
            if imp['conditional'] or (imp['module'], imp['intrinsic']) != (module, intrinsic):
                continue
            if imp['only'] is not None and (local, original) in imp['only']:
                return True
            if imp['only'] is None and not imp['renames'] and local == original:
                return True
        return False

    def parentProvides(key):
        module, intrinsic, local, original = key
        return unit[0] == 'submodule' and ((module, local, original) in parentImports
                                           or ((module, '*', '*') in parentImports and local == original))

    def visible(region, key):
        """Would this entity be visible under its local name in the region, if every procedure-level import of it were
        removed, and replaced by one in the specification part (or by the parent's)? Entities imported by the region
        itself, or by an enclosing procedure, are visible; otherwise the name must reach the specification part without
        meeting a different import of the same name."""
        local = key[2]
        while region != 0:
            if providesItself(region, key):
                return True
            if importsLocalName(region, local) or regions[region]['unknownImports']:
                return False
            region = regions[region]['parent']
        # No procedure imports the name, so here it refers to a local entity, an entity of the unit's specification part
        # or its parent, or an intrinsic procedure. Moving the import into the specification part could change that, so
        # it is safe only if the specification part already imports this entity, or the parent provides it (and the
        # specification part does not declare a different entity of the same name).
        if unitProvides(key):
            return True
        if importsLocalName(0, local) or regions[0]['unknownImports'] or mentionKey(local) in regions[0]['mentions']:
            return False
        return parentProvides(key)

    candidates = set()
    for identifier, region in regions.items():
        if identifier == 0 or region['interface']:
            continue
        for imp in region['imports']:
            if imp['conditional'] or imp['only'] is None:
                continue
            for local, original in imp['only']:
                candidates.add((imp['module'], imp['intrinsic'], local, original))
    move, drop = set(), set()
    for key in sorted(candidates):
        if not isMovable(key[2]):
            continue
        # A procedure which both imports a name and declares it is in error. Leave its import in place, so that the
        # compiler still reports that (once the import was moved, the declaration would simply hide it).
        if any(key[2] in region['declared'] and any(
                (imp['module'], imp['intrinsic']) == key[:2] and imp['only'] is not None and (key[2], key[3]) in imp['only']
                for imp in region['imports']) for region in regions.values()):
            continue
        token = mentionKey(key[2])
        if not all(visible(identifier, key) for identifier, region in regions.items()
                   if not region['interface'] and (token in region['mentions'] or token in region.get('typeMentions', ())
                                                   or importsLocalName(identifier, key[2]))):
            continue
        if unitProvides(key) or parentProvides(key):
            drop.add(key)
        else:
            move.add(key)
    removed = move | drop
    if not removed:
        return list(lines), 0, 0
    output = list(lines)
    removedStatements = 0
    for identifier, region in regions.items():
        if identifier == 0 or region['interface']:
            continue
        for imp in region['imports']:
            if imp['conditional'] or imp['only'] is None:
                continue
            keep = [(local, original) for local, original in imp['only']
                    if (imp['module'], imp['intrinsic'], local, original) not in removed]
            if len(keep) == len(imp['only']):
                continue
            first, last = stmts[imp['statement']][0], stmts[imp['statement']][1]
            if any(other[0] <= last and other[1] >= first and other[2] != 'blank'
                   for j, other in enumerate(stmts) if j != imp['statement']):
                raise Unsupported('a USE statement shares a line with another statement')
            indent = re.match(r'^\s*', lines[first]).group(0)
            for j in range(first, last + 1):
                output[j] = ''
            if keep:
                nature = ', intrinsic' if imp['intrinsic'] else ''
                output[first] = (f"{indent}use{nature} :: {imp['module']}, only : "
                                 + ', '.join(local if local == original else f'{local} => {original}'
                                             for local, original in keep))
            else:
                removedStatements += 1
    byModule = {}
    for module, intrinsic, local, original in sorted(move):
        byModule.setdefault((module, intrinsic), []).append(local if local == original else f'{local} => {original}')
    if byModule:
        unitLine = stmts[unitIndex][0]
        if stmts[unitIndex][1] != unitLine:
            raise Unsupported('continued unit statement')
        added = [f"use{', intrinsic' if intrinsic else ''} :: {module}, only : {', '.join(names)}"
                 for (module, intrinsic), names in byModule.items()]
        output[unitLine] = output[unitLine].rstrip() + ' ; ' + ' ; '.join(added)
        if unit[0] == 'module' and not defaultPrivate:
            containsLine = stmts[unitContains][0]
            output[containsLine] = (re.match(r'^\s*', lines[containsLine]).group(0) + 'private :: '
                                    + ', '.join(sorted({local for _, _, local, _ in move})) + ' ; '
                                    + lines[containsLine].lstrip())
    return output, removedStatements, len(move)


def main():
    parser = argparse.ArgumentParser(description='Consolidate procedure-level USE statements of a preprocessed Fortran '
                                                 'module or submodule into its specification part.')
    parser.add_argument('source', help='the preprocessed source file')
    parser.add_argument('output', help='the file to write')
    parser.add_argument('--moduleDependencies', help='the generated module dependencies makefile, used to find the '
                                                     'parent of a submodule')
    parser.add_argument('--report', action='store_true', help='report what was done')
    arguments = parser.parse_args()
    with open(arguments.source, errors='replace') as source:
        lines = source.read().split('\n')
    try:
        try:
            parentImports = None
            match = re.search(r'^\s*submodule\s*\(\s*([\w:]+)\s*\)', '\n'.join(lines), re.IGNORECASE | re.MULTILINE)
            if match:
                parentSource = findParentSource(match.group(1), arguments.moduleDependencies)
                if parentSource:
                    with open(parentSource, errors='replace') as parentFile:
                        parentImports = specificationImports(parentFile.read().split('\n'))
            output, removedStatements, moved = consolidate(lines, parentImports)
            status = f'{removedStatements} USE statements removed, {moved} names moved'
        except Unsupported as error:
            output, status = lines, f'unchanged: {error}'
    except Exception as error:
        # Consolidation is an optimization only: never let a failure in it stop the build.
        print(f'consolidateUses.py: warning: "{arguments.source}" left unchanged: {type(error).__name__}: {error}',
              file=sys.stderr)
        output, status = lines, 'unchanged: error'
    if len(output) != len(lines):
        print(f'consolidateUses.py: warning: "{arguments.source}" left unchanged: line count changed', file=sys.stderr)
        output, status = lines, 'unchanged: error'
    with open(arguments.output, 'w') as destination:
        destination.write('\n'.join(output))
    if arguments.report:
        print(f'{arguments.source}: {status}')


if __name__ == '__main__':
    main()
