#!/usr/bin/env python3
"""Audit which functionClasses can be exposed via the Python interface.

Walks ``source/`` to discover every Galacticus ``<functionClass>`` directive
and every concrete implementation of each one, parses each impl's
``ConstructorInternal`` argument list, and classifies the *class as a whole*
(union of its concrete impls' classifications) into one of three buckets:

  READY              — at least one concrete impl's constructor uses only
                       pipeline-supported arg types and depends only on
                       functionClasses that are already registered in
                       libraryClasses.xml.  Adding ``<className/>`` to
                       libraryClasses.xml today should compile and surface
                       the impl on the Python side.

  MISSING-DEP        — would be ready if N specific other functionClasses
                       were also registered.  The audit reports the unmet
                       dependency edges and runs a closure analysis: starting
                       from the currently-registered set, iteratively pulls
                       in classes whose deps are now satisfied, and reports
                       how many classes the closure ultimately reaches.

  PIPELINE-BLOCKED   — every concrete impl needs a generator feature we
                       don't have yet (output arrays, complex args,
                       procedure pointers, etc.).  Grouped by blocker so
                       you can see which feature would unblock the most
                       classes.

The script is standalone — it doesn't require ``directiveLocations.xml`` to
be built first.  It does need the rest of the python/ tree on PYTHONPATH
(or installed via ``pip install -e .``) so SourceTree.parse_file is
available; the Makefile already exports PYTHONPATH for build commands.

Andrew Benson (drafted with assistance from Claude 2026)
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT  = Path(__file__).resolve().parent.parent.parent
SOURCE_DIR = REPO_ROOT / 'source'
LIB_XML    = SOURCE_DIR / 'libraryClasses.xml'

# Make the python/ tree importable when the script is run directly without
# a `pip install -e .` having happened.
sys.path.insert(0, str(REPO_ROOT / 'python'))

import xml.etree.ElementTree as ET                     # noqa: E402

from Galacticus.Build import SourceTree                # noqa: E402


# Match a fixed-size dimension(N) attribute (N a positive integer literal).
# Same predicate the generator uses; duplicated here so the audit doesn't
# have to import from libraryInterfaces.py (which would also drag in the
# generator's unrelated emitter machinery).
_DIM_FIXED_RX = re.compile(r'^dimension\s*\(\s*(\d+)\s*\)$')

# Match `len=N` in a character type-spec.  Used to recognise fixed-length
# character arrays (`character(len=N), dimension(:)`), which the pipeline
# now plumbs through; variable-length forms (`len=*`, `len=:`) are still
# unsupported because they have no fixed stride at the byte boundary.
_CHAR_LEN_RX = re.compile(r'^len\s*=\s*(\d+)$')

# Recognise an Internal-suffixed module-procedure name following the
# Galacticus convention <short>Constructor[Internal[Suffix]] OR the
# alternative <short>Internal form (used by the merger-tree walkers,
# e.g. allAndFormationNodesInternal).  We accept either:
#   • a name ending in "internal" (catches both `<short>ConstructorInternal`
#     and `<short>Internal`), or
#   • a name containing "constructorinternal" (catches the rarer
#     ConstructorInternalType / ConstructorInternalDefined variants
#     used to disambiguate multiple internal constructors).
# `<name>ConstructorParameters` and similar do not satisfy either rule
# and are correctly rejected.
def _is_internal_constructor_name(name):
    lower = name.lower()
    return lower.endswith('internal') or 'constructorinternal' in lower


# ---------------------------------------------------------------------------
# Discovery: walk source/ and find every functionClass directive and every
# concrete implementation declaration.
# ---------------------------------------------------------------------------

def discover_function_classes(source_dir):
    """Return ``(all_fcs, fc_files)`` where:
      all_fcs   — set of every functionClass name found in source/.
      fc_files  — {fc_name: {set of paths containing concrete <fc_name name="…">}}
    """
    all_fcs  = set()
    fc_files = defaultdict(set)

    paths = sorted(source_dir.glob('*.F90'))

    # Pass 1: collect every functionClass directive's name.
    for path in paths:
        text = path.read_text(errors='replace')
        for m in re.finditer(
                r'<functionClass>(?:[^<]|<(?!/functionClass>))*?'
                r'<name>([^<]+)</name>',
                text, re.DOTALL):
            all_fcs.add(m.group(1))

    # Pass 2: now that we know the fc-name universe, find impl directives.
    # An impl declaration looks like:  <fcName name="impl_name" …>
    impl_directive_rx = re.compile(
        r'<(' + '|'.join(re.escape(c) for c in all_fcs)
        + r')\s+name="([^"]+)"',
    )
    for path in paths:
        text = path.read_text(errors='replace')
        for m in impl_directive_rx.finditer(text):
            fc_files[m.group(1)].add(path)

    return all_fcs, fc_files


def parse_impls_in_file(impl_file, fc_name):
    """Walk *impl_file*'s SourceTree and return a list of dicts describing
    every concrete impl of *fc_name* found in it:

        {'name'       : impl name from the directive
         'is_abstract': bool
         'is_excluded': bool (currently unused — relies on libraryClasses.xml
                              <impl exclude="yes"/> override which we don't
                              read here)
         'internals'  : list of *ConstructorInternal* names found in the
                        impl's generic interface block (>1 ⇒ ambiguous).
         'args'       : list of {name, intrinsic, type, attributes} dicts
                        for the chosen Internal constructor's arguments,
                        or [] if none could be resolved.
        }

    Mirrors the second pass of ``libraryInterfaces._process_implementations``
    (kept as a duplicate rather than an import to avoid coupling the audit
    to the generator's non-classification machinery).
    """
    tree = SourceTree.parse_file(str(impl_file))

    impls            = []
    impl_name        = None
    is_abstract      = False
    name_constructor = None
    internals_seen   = []
    args_constructor = []

    def _flush():
        if impl_name is not None:
            impls.append({
                'name'       : impl_name,
                'is_abstract': is_abstract,
                'internals'  : list(internals_seen),
                'args'       : list(args_constructor),
            })

    for node in SourceTree.walk_tree(tree):
        if node['type'] == fc_name:
            _flush()
            impl_name        = node.get('directive', {}).get('name')
            is_abstract      = (node.get('directive', {})
                                    .get('abstract', 'no') == 'yes')
            name_constructor = None
            internals_seen   = []
            args_constructor = []

        elif (impl_name
              and node['type'] == 'interface'
              and node.get('name', '').lower() == impl_name.lower()):
            child = node.get('firstChild')
            while child:
                if child['type'] == 'moduleProcedure':
                    for n in child.get('names', []):
                        if _is_internal_constructor_name(n):
                            internals_seen.append(n)
                child = child.get('sibling')
            if len(internals_seen) == 1:
                name_constructor = internals_seen[0]

        elif (name_constructor
              and node['type'] == 'function'
              and node.get('name', '').lower() == name_constructor.lower()):
            opener = node.get('opener', '')
            m = re.search(
                r'function\s+' + re.escape(name_constructor) + r'\s*\(([^)]+)\)',
                opener, re.IGNORECASE)
            if m:
                # Strip Fortran line-continuation `&` and any whitespace from
                # each captured argument; multi-line openers otherwise leave
                # `&\n  &` fragments inside the names, breaking name->decl matching.
                args_constructor = [{'name': re.sub(r'[\s&]+', '', a)}
                                    for a in m.group(1).split(',')]
            child = node.get('firstChild')
            while child:
                if child['type'] == 'declaration':
                    for decl in child.get('declarations', []):
                        for var_name in decl.get('variableNames', []):
                            for arg in args_constructor:
                                if arg['name'].lower() == var_name.lower():
                                    arg['intrinsic']  = decl.get('intrinsic')
                                    arg['type']       = decl.get('type')
                                    arg['attributes'] = decl.get('attributes', [])
                child = child.get('sibling')

    _flush()
    return impls


# ---------------------------------------------------------------------------
# Classification: per-impl predicate that distinguishes pipeline-blockers
# from missing-dependency situations.  Roughly mirrors
# libraryInterfaces._unsupported_arg, but with a richer return so the audit
# can tell apart "Foo is a functionClass we just haven't registered" from
# "Foo isn't a functionClass at all".
# ---------------------------------------------------------------------------

def classify_constructor(args, all_fcs, registered):
    """Return ``(missing_deps, pipeline_reasons)``.

    *missing_deps* is the set of functionClass names this constructor depends
    on that aren't currently registered in libraryClasses.xml.  An empty
    set means the constructor's dependencies (if any) are all satisfied.

    *pipeline_reasons* is a list of human-readable strings describing
    arg types the generator can't handle today — empty iff every arg is
    pipeline-supported.
    """
    missing_deps     = set()
    pipeline_reasons = []

    for arg in args:
        name      = arg.get('name', '?')
        intrinsic = arg.get('intrinsic')
        attrs     = list(arg.get('attributes', []))
        type_spec = (arg.get('type') or '').strip()

        if intrinsic in ('complex', 'double complex'):
            pipeline_reasons.append(
                f"complex arg ({name}: {intrinsic}({type_spec}))")
            continue

        if intrinsic == 'class':
            if type_spec == '*':
                # class(*) needs an explicit override in libraryClasses.xml.
                # This audit doesn't look up overrides (it would have to
                # re-parse libraryClasses.xml structure), so flag it as a
                # pipeline-ish blocker — the user knows from existing cases
                # that a hint can resolve it.
                pipeline_reasons.append(
                    f"class(*) without override ({name})")
                continue
            if type_spec.endswith('Class'):
                stem = type_spec[:-5]
                if stem in all_fcs:
                    if stem not in registered:
                        missing_deps.add(stem)
                    continue   # fc — registered or just missing-dep
            # class(SomethingElse) — not a functionClass at all.
            pipeline_reasons.append(
                f"class({type_spec}) — not a functionClass ({name})")
            continue

        # Dimension shape checks.
        dim_attr = next((a for a in attrs if a.startswith('dimension')), None)
        if dim_attr:
            # Fixed-length character arrays at deferred shape are now
            # supported alongside the numeric cases — see the analogous
            # extension in libraryInterfaces._unsupported_arg.
            is_supported_shape = (
                (intrinsic in ('double precision', 'integer')
                 and (dim_attr == 'dimension(:)'
                      or _DIM_FIXED_RX.match(dim_attr)))
                or
                (intrinsic == 'character'
                 and dim_attr == 'dimension(:)'
                 and _CHAR_LEN_RX.match(type_spec))
            )
            if not is_supported_shape:
                pipeline_reasons.append(
                    f"unsupported array shape "
                    f"({name}: {intrinsic} {dim_attr})")
            elif 'allocatable' in attrs:
                pipeline_reasons.append(
                    f"allocatable array ({name})")
            elif 'intent(in)' not in attrs:
                pipeline_reasons.append(
                    f"non-intent(in) array ({name})")
            continue

        if intrinsic == 'procedure':
            pipeline_reasons.append(f"procedure-pointer arg ({name})")
            continue

    return missing_deps, pipeline_reasons


def aggregate_class(impls, all_fcs, registered):
    """Reduce a list of impl dicts into a single class-level verdict.

    Status precedence: ready > missing-dep > pipeline-blocked > no-concrete-impls.
    A class is READY if any concrete impl has empty pipeline_reasons AND
    empty missing_deps.  MISSING-DEP if no impl is ready but at least one
    impl has empty pipeline_reasons (just unmet deps).  Otherwise
    PIPELINE-BLOCKED.
    """
    concrete = [i for i in impls if not i['is_abstract']]
    if not concrete:
        return {'status': 'no-concrete-impls', 'impls': impls}

    classified = []
    for impl in concrete:
        if len(impl['internals']) > 1:
            # Generator skips ambiguous-Internal impls; treat as pipeline-blocked.
            classified.append((impl, set(),
                               [f"ambiguous Internal constructors: "
                                f"{', '.join(impl['internals'])}"]))
            continue
        deps, reasons = classify_constructor(impl['args'], all_fcs, registered)
        classified.append((impl, deps, reasons))

    ready = [(i, d, r) for (i, d, r) in classified if not r and not d]
    if ready:
        return {'status'   : 'ready',
                'impls'    : impls,
                'ready_impl': ready[0][0]['name']}

    missing_dep = [(i, d, r) for (i, d, r) in classified if not r and d]
    if missing_dep:
        union = set()
        for _, d, _ in missing_dep:
            union |= d
        # Pick the representative impl with the *smallest* missing-dep set.
        rep = min(missing_dep, key=lambda x: (len(x[1]), x[0]['name']))
        return {'status'    : 'missing-dep',
                'impls'     : impls,
                'deps'      : union,
                'min_deps'  : rep[1],
                'rep_impl'  : rep[0]['name']}

    # Everything is pipeline-blocked.  Group reasons by category for the report.
    reasons = []
    for impl, _, r in classified:
        for s in r:
            reasons.append((impl['name'], s))
    return {'status' : 'pipeline-blocked',
            'impls'  : impls,
            'reasons': reasons}


# ---------------------------------------------------------------------------
# Closure analysis: starting from the currently-registered set, iteratively
# pull in any class whose unmet deps have been satisfied.  Reports the
# total reach and the order in which classes get pulled in.
# ---------------------------------------------------------------------------

def closure(audit, registered):
    """Return ``(closure_set, levels)``.  *levels* is a list of lists; entry
    *i* is the classes that became reachable in iteration *i*.
    """
    closure_set = set(registered)
    levels      = []
    while True:
        added = []
        for fc, info in audit.items():
            if fc in closure_set:
                continue
            if info['status'] == 'ready':
                added.append(fc)
            elif info['status'] == 'missing-dep' \
                    and info['deps'].issubset(closure_set):
                added.append(fc)
        if not added:
            break
        added.sort()
        levels.append(added)
        closure_set |= set(added)
    return closure_set, levels


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def load_registered_classes(xml_path):
    """Parse libraryClasses.xml and return the set of class tags inside <classes>."""
    if not xml_path.exists():
        return set()
    classes_el = ET.parse(xml_path).getroot().find('classes')
    if classes_el is None:
        return set()
    return {child.tag for child in classes_el}


def fmt_section(title):
    bar = '=' * len(title)
    return f"\n{bar}\n{title}\n{bar}\n"


def main():
    print("Discovering functionClasses…", file=sys.stderr)
    all_fcs, fc_files = discover_function_classes(SOURCE_DIR)
    registered        = load_registered_classes(LIB_XML)
    print(f"  {len(all_fcs)} functionClass(es); "
          f"{len(registered)} currently registered.", file=sys.stderr)

    audit = {}
    for i, fc in enumerate(sorted(all_fcs), 1):
        print(f"\rParsing impls… {i}/{len(all_fcs)} ({fc[:40]:<40})",
              end='', file=sys.stderr)
        impls = []
        for path in sorted(fc_files.get(fc, [])):
            impls.extend(parse_impls_in_file(path, fc))
        audit[fc] = aggregate_class(impls, all_fcs, registered)
    print("\n", file=sys.stderr)

    # Buckets.
    by_status = defaultdict(list)
    for fc, info in audit.items():
        by_status[info['status']].append(fc)

    print(fmt_section("Summary"))
    print(f"  Total functionClasses:    {len(all_fcs)}")
    print(f"  Currently registered:     {len(registered)} "
          f"({100*len(registered)//len(all_fcs)}%)")
    for status in ('ready', 'missing-dep',
                   'pipeline-blocked', 'no-concrete-impls'):
        n = len(by_status[status])
        print(f"  {status:18s}        {n:4d} "
              f"({100*n//len(all_fcs):3d}%)")
    print(f"  (registered classes are also counted in their bucket above —")
    print(f"  they all evaluate as 'ready' so they don't dominate the lists.)")

    # READY (excluding already-registered).
    ready_unregistered = sorted(c for c in by_status['ready']
                                if c not in registered)
    print(fmt_section(
        f"READY ({len(ready_unregistered)}) — register today"))
    if ready_unregistered:
        for fc in ready_unregistered:
            rep = audit[fc].get('ready_impl', '?')
            print(f"  <{fc}/>".ljust(60) + f"# e.g. {rep}")
    else:
        print("  (none — every ready class is already registered)")

    # MISSING-DEP.
    missing = sorted(by_status['missing-dep'])
    print(fmt_section(
        f"MISSING-DEP ({len(missing)}) — would be ready with deps registered"))
    if missing:
        for fc in missing:
            info = audit[fc]
            min_deps = sorted(info['min_deps'])
            print(f"  {fc}".ljust(50) + f"needs: {', '.join(min_deps)}")
    else:
        print("  (none)")

    # Closure analysis.
    closure_set, levels = closure(audit, registered)
    print(fmt_section(
        f"CLOSURE — adding the missing-dep classes pulls in "
        f"{len(closure_set) - len(registered)} more"))
    print(f"  Currently reachable:  {len(registered)}")
    print(f"  After closure:        {len(closure_set)} "
          f"({100*len(closure_set)//len(all_fcs)}% of all functionClasses)")
    if levels:
        print(f"  Levels:")
        for i, lvl in enumerate(levels, 1):
            print(f"    Level {i}:  {len(lvl):3d} classes:  "
                  f"{', '.join(lvl[:6])}"
                  f"{'…' if len(lvl) > 6 else ''}")

    # PIPELINE-BLOCKED.  Group by blocker to highlight which generator
    # feature would unblock the most classes.
    blocked = sorted(by_status['pipeline-blocked'])
    print(fmt_section(
        f"PIPELINE-BLOCKED ({len(blocked)}) — needs generator work"))
    blocker_counts = defaultdict(set)
    for fc in blocked:
        for impl_name, reason in audit[fc]['reasons']:
            # Bucket the reason by its leading category.
            head = re.split(r'\s+\(', reason, maxsplit=1)[0]
            blocker_counts[head].add(fc)
    if blocker_counts:
        print(f"  Top blockers (number of classes affected):")
        for head, classes in sorted(blocker_counts.items(),
                                    key=lambda kv: -len(kv[1])):
            print(f"    {len(classes):3d}  {head}")
    print()
    if blocked:
        print(f"  Per-class reasons (first reason per impl):")
        for fc in blocked:
            reason_set = sorted(set(r for _, r in audit[fc]['reasons']))
            print(f"    {fc}".ljust(50) + f"{reason_set[0]}")

    # NO-CONCRETE-IMPLS — generally these are abstract base classes that
    # only matter as deps; surface them but don't dwell.
    no_impls = sorted(by_status['no-concrete-impls'])
    if no_impls:
        print(fmt_section(
            f"NO-CONCRETE-IMPLS ({len(no_impls)}) — abstract bases only"))
        print(f"  {', '.join(no_impls)}")


if __name__ == '__main__':
    main()
