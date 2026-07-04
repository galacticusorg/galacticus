#!/usr/bin/env python3
"""Audit which functionClasses can be exposed via the Python interface.

Two passes:

1. **Constructors.**  Walk ``source/`` to discover every Galacticus
   ``<functionClass>`` directive and every concrete implementation of
   each one, parse each impl's ``ConstructorInternal`` argument list,
   and classify the *class as a whole* (union of its concrete impls'
   classifications) into one of three buckets:

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

2. **Methods.**  Walk every ``<method>`` directive on every registered
   functionClass, parse its return-type and ``<argument>`` list, and
   classify each method individually.  Blockers are split into
   *in-scope* and *out-of-scope* buckets — the latter are deferred
   categories the team has explicitly chosen not to tackle yet
   (``type(treeNode)``-style internal derived types, non-fc class
   hierarchies); see :func:`_is_out_of_scope_reason` for the predicate.
   Both buckets are reported separately so the in-scope worklist
   doesn't get drowned out by the deferred backlog.

The script is standalone — it doesn't require ``directiveLocations.xml`` to
be built first.  It does need the rest of the python/ tree on PYTHONPATH
(or installed via ``pip install -e .``) so SourceTree.parse_file is
available; the Makefile already exports PYTHONPATH for build commands.

**This is a manual developer tool**: nothing in the Makefile or CI invokes
it — run it by hand when planning library-interface coverage work.  The
classification rules are imported from LibraryInterfaces.Classification —
the same module the generator uses — so the audit's notion of "supported"
tracks the generator by construction.

Andrew Benson (drafted with assistance from Claude 2026)
"""

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
from Galacticus.Build.SourceTree.Parse import Declarations  # noqa: E402
from LibraryInterfaces.Hierarchy import build_type_hierarchy  # noqa: E402


# The classification rules (regexes, return-type tables, and the argument
# predicate) are shared with the generator via
# LibraryInterfaces.Classification — the audit's verdicts follow the
# generator's by construction rather than by hand-synced copies.
from LibraryInterfaces.Classification import (   # noqa: E402
    ENUM_RETURN_RX               as _ENUM_RETURN_RX,
    ARRAY_RETURN_RX              as _ARRAY_RETURN_RX,
    CLASS_RETURN_RX              as _CLASS_RETURN_RX,
    DYNAMIC_ARRAY_RETURN_RX      as _DYNAMIC_ARRAY_RETURN_RX,
    DYNAMIC_ARRAY_RETURN_2D_RX   as _DYNAMIC_ARRAY_RETURN_2D_RX,
    RETURN_TYPE_ALIASES          as _RETURN_TYPE_ALIASES,
    SCALAR_RETURN_OK             as _SCALAR_RETURN_OK,
    classify_arg,
    is_internal_constructor_name as _is_internal_constructor_name,
    unsupported_output_array_method as _unsupported_output_array_method,
)


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

    # rglob: the source tree is hierarchical — a flat glob finds only the
    # handful of top-level files and silently discovers nothing.
    paths = sorted(source_dir.rglob('*.F90'))

    # Pass 1: collect every functionClass directive's name.
    for path in paths:
        text = path.read_text(errors='replace')
        for m in re.finditer(
                r'<functionClass(?:\s[^>]*)?>(?:[^<]|<(?!/functionClass>))*?'
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


def parse_impls_in_file(impl_file, fc_name, internal_selectors=None):
    """Walk *impl_file*'s SourceTree and return a list of dicts describing
    every concrete impl of *fc_name* found in it:

        {'name'       : impl name from the directive
         'is_abstract': bool
         'is_excluded': bool (currently unused — relies on libraryClasses.xml
                              <impl exclude="yes"/> override which we don't
                              read here)
         'internals'  : list of *ConstructorInternal* names found in the
                        impl's generic interface block (>1 ⇒ ambiguous,
                        unless *internal_selectors* names one of them).
         'chosen_internal':
                        the Internal constructor name picked by an
                        `<constructor internal=…/>` hint when there are
                        multiple candidates; None otherwise.
         'args'       : list of {name, intrinsic, type, attributes} dicts
                        for the chosen Internal constructor's arguments,
                        or [] if none could be resolved.
        }

    *internal_selectors* is a ``dict impl_name -> chosen_internal_name``
    built from libraryClasses.xml's
    ``<constructor internal="…"/>`` attributes; mirrors the generator's
    disambiguation hint.

    Mirrors the second pass of ``libraryInterfaces._process_implementations``
    (kept as a duplicate rather than an import to avoid coupling the audit
    to the generator's non-classification machinery).
    """
    if internal_selectors is None:
        internal_selectors = {}
    tree = SourceTree.parse_file(str(impl_file))

    impls            = []
    impl_name        = None
    is_abstract      = False
    name_constructor = None
    internals_seen   = []
    args_constructor = []

    chosen_internal = None

    def _flush():
        if impl_name is not None:
            impls.append({
                'name'            : impl_name,
                'is_abstract'     : is_abstract,
                'internals'       : list(internals_seen),
                'chosen_internal' : chosen_internal,
                'args'            : list(args_constructor),
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
            chosen_internal  = None

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
            elif len(internals_seen) > 1:
                # libraryClasses.xml can break the tie via
                # `<constructor internal="…"/>`; mirrors the generator.
                hinted = internal_selectors.get(impl_name)
                if hinted and hinted in internals_seen:
                    name_constructor = hinted
                    chosen_internal  = hinted

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

def classify_constructor(args, all_fcs, registered, overridden_args=frozenset(),
                         class_hierarchy=None, null_filled_args=frozenset(),
                         absent_filled_args=frozenset()):
    """Return ``(missing_deps, pipeline_reasons)``.

    *missing_deps* is the set of functionClass names this constructor depends
    on that aren't currently registered in libraryClasses.xml.  An empty
    set means the constructor's dependencies (if any) are all satisfied.

    *pipeline_reasons* is a list of human-readable strings describing
    arg types the generator can't handle today — empty iff every arg is
    pipeline-supported.

    *overridden_args* is the set of constructor argument names for which
    libraryClasses.xml supplies a ``<constructor><argument …/></constructor>``
    type hint (bypassing the ``class(*)`` blocker); *null_filled_args* /
    *absent_filled_args* are the names carrying ``value='null'`` /
    ``value='absent'`` overrides.

    Delegates per-argument classification to the generator-shared
    :func:`LibraryInterfaces.Classification.classify_arg` (audit mode), so
    the audit's notion of "supported" is the generator's by construction.
    """
    # Express the audit's name-set override parameters in the override-dict
    # form the shared predicate consumes.
    overrides = (
        [{'name': n, 'value': 'null'}   for n in null_filled_args]
        + [{'name': n, 'value': 'absent'} for n in absent_filled_args]
        + [{'name': n} for n in overridden_args
           if n not in null_filled_args and n not in absent_filled_args]
    )
    missing_deps     = set()
    pipeline_reasons = []
    for arg in args:
        verdict = classify_arg(
            arg, registered,
            constructor_overrides=overrides,
            class_hierarchy=class_hierarchy,
            known_function_classes=set(all_fcs),
        )
        if verdict is None:
            continue
        kind, payload = verdict
        if kind == 'missing-dep':
            missing_deps |= payload
        else:
            pipeline_reasons.append(f"{payload} ({arg.get('name', '?')})")
    return missing_deps, pipeline_reasons


def aggregate_class(impls, all_fcs, registered, overrides, class_hierarchy=None,
                    null_filled=None, absent_filled=None):
    """Reduce a list of impl dicts into a single class-level verdict.

    Status precedence: ready > missing-dep > pipeline-blocked > no-concrete-impls.
    A class is READY if any concrete impl has empty pipeline_reasons AND
    empty missing_deps.  MISSING-DEP if no impl is ready but at least one
    impl has empty pipeline_reasons (just unmet deps).  Otherwise
    PIPELINE-BLOCKED.

    The full per-impl classifications are also returned under
    ``classified`` so the report can surface impl-level drops even when
    the class as a whole evaluates as ready.  Each entry is a tuple
    ``(impl_dict, deps_set, reasons_list)`` for one concrete impl;
    ``deps`` and ``reasons`` are both empty iff the impl is itself ready.
    """
    concrete = [i for i in impls if not i['is_abstract']]
    if not concrete:
        return {'status': 'no-concrete-impls', 'impls': impls,
                'classified': []}

    classified = []
    for impl in concrete:
        if len(impl['internals']) > 1 and not impl.get('chosen_internal'):
            # Generator skips genuinely ambiguous-Internal impls; treat
            # as pipeline-blocked.  When libraryClasses.xml's
            # `<constructor internal="…"/>` selector picked one, the
            # chosen constructor's args were already parsed and we fall
            # through to the normal classify_constructor path.
            classified.append((impl, set(),
                               [f"ambiguous Internal constructors: "
                                f"{', '.join(impl['internals'])}"]))
            continue
        deps, reasons = classify_constructor(
            impl['args'], all_fcs, registered,
            overrides.get(impl['name'], frozenset()),
            class_hierarchy=class_hierarchy,
            null_filled_args=((null_filled or {})
                              .get(impl['name'], frozenset())),
            absent_filled_args=((absent_filled or {})
                                .get(impl['name'], frozenset())))
        classified.append((impl, deps, reasons))

    ready = [(i, d, r) for (i, d, r) in classified if not r and not d]
    if ready:
        return {'status'    : 'ready',
                'impls'     : impls,
                'classified': classified,
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
                'classified': classified,
                'deps'      : union,
                'min_deps'  : rep[1],
                'rep_impl'  : rep[0]['name']}

    # Everything is pipeline-blocked.  Group reasons by category for the report.
    reasons = []
    for impl, _, r in classified:
        for s in r:
            reasons.append((impl['name'], s))
    return {'status'    : 'pipeline-blocked',
            'impls'     : impls,
            'classified': classified,
            'reasons'   : reasons}


# ---------------------------------------------------------------------------
# Methods: discover every `<method>` directive on every functionClass, parse
# its return-type and arg list, and classify each one with the same
# predicate the generator uses.  An "out-of-scope" tag marks reasons that
# the team has explicitly deferred (`type(treeNode)` and other internal
# derived types, non-functionClass polymorphic hierarchies); see
# :func:`_is_out_of_scope_reason` for the predicate.
# ---------------------------------------------------------------------------

# Match a `<method name="…">…</method>` block inside a functionClass
# directive's body.  Bodies contain other XML elements (description, type,
# argument, …); `[\s\S]*?` keeps the inner match minimal without choking on
# the embedded `<…>` markup.
_METHOD_RX   = re.compile(
    r'<method\s+name="(?P<name>[^"]+)"[^>]*>(?P<body>[\s\S]*?)</method>'
)
_M_TYPE_RX   = re.compile(r'<type>([^<]+)</type>')
_M_ARG_RX    = re.compile(r'<argument>([^<]+)</argument>')
# The `(?:\s[^>]*)?` allows the optional attributes real directives carry
# (`<functionClass docformat="rst">`, …); a bare `<functionClass>` form
# does not occur in the source tree, so without it this block matches
# nothing and the whole methods pass silently reports zero methods.
_FC_BLOCK_RX = re.compile(
    r'<functionClass(?:\s[^>]*)?>(?P<body>[\s\S]*?)</functionClass>'
)
_FC_NAME_RX  = re.compile(r'<name>([^<]+)</name>')


def discover_methods(source_dir, all_fcs):
    """Return ``{fc_name: [{'name': …, 'return': str, 'args': [decl, …]}, …]}``
    where each *decl* is the dict produced by
    :func:`Declarations.parse_declaration` (so the per-arg shape matches
    what :func:`classify_constructor` already consumes).

    Methods are extracted by regex from the `<functionClass>…</functionClass>`
    body — the SourceTree pass that the impl walk uses would surface the
    same data, but at ~5× the cost (every functionClass file would have to
    be parsed in full).  The regex form is enough because the directive
    markup is uniformly XML-like and the audit only needs the method
    signature, not the surrounding Fortran.
    """
    methods = defaultdict(list)
    for path in sorted(source_dir.rglob('*.F90')):
        text = path.read_text(errors='replace')
        for fc_match in _FC_BLOCK_RX.finditer(text):
            body    = fc_match.group('body')
            name_m  = _FC_NAME_RX.search(body)
            if not name_m:
                continue
            fc = name_m.group(1).strip()
            if fc not in all_fcs:
                continue
            for mm in _METHOD_RX.finditer(body):
                mbody    = mm.group('body')
                t_m      = _M_TYPE_RX.search(mbody)
                ret_type = (t_m.group(1).strip() if t_m else 'void')
                args     = []
                for am in _M_ARG_RX.finditer(mbody):
                    decl = Declarations.parse_declaration(am.group(1))
                    if not decl:
                        # An <argument> line the Fortran parser refused —
                        # surface it as a synthetic blocker rather than
                        # silently swallowing it (otherwise an obviously
                        # broken method would look ready).
                        args.append({'name': '?',
                                     'intrinsic': '?',
                                     'type': am.group(1).strip(),
                                     'attributes': []})
                        continue
                    for var_name in decl.get('variableNames', []):
                        args.append({
                            'name'     : var_name,
                            'intrinsic': decl.get('intrinsic'),
                            'type'     : decl.get('type'),
                            'attributes': decl.get('attributes', []),
                        })
                methods[fc].append({
                    'name'   : mm.group('name'),
                    'return' : ret_type,
                    'args'   : args,
                })
    return methods


def classify_method_return(ret_type, all_fcs, registered,
                           class_hierarchy=None):
    """Return ``(missing_deps, reasons)`` for a method's return type.

    Mirrors the return-type switch in
    ``libraryInterfaces.interfaces_methods`` — anything not handled there
    surfaces here as a pipeline blocker.  Class-typed returns whose stem is
    a known-but-unregistered functionClass are reported as a missing dep
    instead, so they participate in the closure analysis the same way
    class-typed constructor args do.
    """
    ret = _RETURN_TYPE_ALIASES.get(ret_type.strip(), ret_type.strip())
    if ret in _SCALAR_RETURN_OK:
        return set(), []
    if _ENUM_RETURN_RX.match(ret):
        return set(), []
    if _ARRAY_RETURN_RX.match(ret):
        return set(), []
    if _DYNAMIC_ARRAY_RETURN_RX.match(ret):
        return set(), []
    if _DYNAMIC_ARRAY_RETURN_2D_RX.match(ret):
        return set(), []
    m = _CLASS_RETURN_RX.match(ret)
    if m:
        stem = m.group(1)[:-5]   # strip trailing "Class"
        if stem in all_fcs:
            if stem not in registered:
                return {stem}, []
            return set(), []
        return set(), [
            f"class(non-fc) return ({ret} — '{stem}' is not a registered "
            f"functionClass)"
        ]
    # Tag common sub-categories so the report's `re.split(r'\s+\(', …)`
    # bucketing tells the distinct backlogs apart.  The "kind=" alias path
    # is a one-liner in the generator's switch; the deferred-shape and
    # dynamic-size return-array paths each need a real out-buffer
    # protocol; the scalar `type(<X>)` returns are the deferred internal-
    # derived-type / pointer-return cases.
    if re.match(r'^integer\s*\(\s*kind\s*=', ret, re.IGNORECASE):
        return set(), [f"kind-alias return type ({ret})"]
    if re.search(r'allocatable', ret, re.IGNORECASE) \
            and re.search(r'dimension', ret, re.IGNORECASE):
        return set(), [f"allocatable-array return type ({ret})"]
    if re.search(r'dimension\s*\(\s*[a-z_]', ret, re.IGNORECASE):
        # `dimension(self%…)` / `dimension(size(…))` etc. — extent
        # computed from runtime state rather than a literal integer.
        return set(), [f"dynamic-size array return type ({ret})"]
    if re.match(r'^type\s*\(', ret, re.IGNORECASE):
        # Includes `type(abundances)`, `type(mergerTree), pointer`, … —
        # the internal-derived-type / pointer-return cases that
        # _is_out_of_scope_reason flags as deferred.
        return set(), [f"scalar derived-type return type ({ret})"]
    return set(), [f"unsupported return type ({ret})"]


# Categorise a blocker reason for the in-scope vs deferred split in the
# method report.  "Out-of-scope" tracks the team decision to defer
# `type(treeNode)`, the other internal derived types, and non-fc class
# hierarchies; everything else (kind aliases, deferred-shape return
# arrays, dynamic-size returns, output-array args, procedure-pointer args,
# class(*) in methods, complex args, …) stays in-scope.
def _is_out_of_scope_reason(reason):
    """True if the blocker is one of the categories deferred to a later
    pass.  Used purely for report bucketing; classify_method_return /
    classify_constructor stay agnostic so the unfiltered classification
    is still available to callers."""
    # class(<X>) where X isn't a registered fc (incl. the abstract-
    # intermediate-failed variant emitted by classify_constructor /
    # classify_method_return).  Catches "class(<X>) — not a functionClass
    # …" and "class(<X>) — only registered functionClasses and class(*)
    # (with override) are supported".
    #
    # Anchor on the blocker's *subject* — the leading `class(...)` — rather
    # than searching the whole reason.  The explanatory tail of every
    # concrete-class blocker literally reads "…and class(*) (with override)
    # are supported", so a whole-string search for `class(*)` misfires on
    # that tail and mis-buckets every non-fc class hierarchy as in-scope.
    # The genuinely in-scope `class(*)` blocker reads "class(*) without a
    # libraryClasses.xml override (…)" — its subject is `*` and it carries
    # no em-dash, so the match below correctly leaves it in-scope.
    m_cls = re.match(r'\s*class\s*\(([^)]*)\)\s*(?:—|--)', reason)
    if m_cls and m_cls.group(1).strip() != '*':
        return True
    if 'class(non-fc) return' in reason:
        return True
    # `procedure(...)` / `class(...)` pointers with intent(out|inout) — a
    # Fortran procedure or object pointer handed back to the caller.
    # Unsupportable in principle (a Python caller can do nothing with one),
    # so they belong with the deferred backlog, not the actionable
    # worklist.  Inbound callback blockers ("procedure(...) —
    # procedure-pointer args are not supported") stay in-scope (candidates
    # for Pipeline._CALLBACK_PROCEDURE_INTERFACES), as do `type(X),
    # pointer` in/outs ("pointer dummy of derived type …" — fixable via a
    # pointer write-back protocol; that wording deliberately avoids
    # matching here or the derived-type rule below).
    if ') pointer output' in reason:
        return True
    # `class(<X>), dimension(:)` args — arrays of polymorphic objects
    # cannot be assembled from Python-held per-object pointers (one
    # dynamic type per array; no intrinsic assignment into polymorphic
    # elements), so these are deferred, not actionable.
    if re.match(r'\s*class\s*\(', reason) and 'array argument' in reason:
        return True
    # Scalar `type(<X>)` returns/args where X isn't varying_string or an
    # enumeration*Type — those are the internal-derived-type cases.
    # `unsupported return type (type(<X>))` and `class(<X>List) …` array
    # cases both come through here.
    m = re.search(r'type\s*\(\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\)', reason)
    if m:
        stem = m.group(1)
        if stem == 'varying_string':
            return False
        if re.match(r'enumeration[a-z0-9_]+type$', stem, re.IGNORECASE):
            return False
        return True
    return False


def aggregate_methods(methods_by_fc, all_fcs, registered,
                      class_hierarchy=None):
    """Classify every method on every functionClass.

    Returns ``{fc: [{'method': method_dict, 'deps': set, 'reasons': list}, …]}``.
    Empty `deps` AND `reasons` mean the method is ready; anything else
    means the generator would drop it (or, if only `deps`, that it would
    become reachable once the named functionClasses are registered).
    """
    out = {}
    for fc, methods in methods_by_fc.items():
        rows = []
        for m in methods:
            ret_deps, ret_reasons = classify_method_return(
                m['return'], all_fcs, registered,
                class_hierarchy=class_hierarchy)
            arg_deps, arg_reasons = classify_constructor(
                m['args'], all_fcs, registered,
                overridden_args=frozenset(),
                class_hierarchy=class_hierarchy)
            reasons = list(ret_reasons) + list(arg_reasons)
            # Whole-method gate for output-array args: classify_constructor
            # accepts each `intent(out), allocatable, dimension(:)` arg
            # per-arg, but the generator only emits output-array methods that
            # are void-returning, optional-free, and whose other args are
            # supported inputs.  Mirror that gate so the audit's readiness
            # tracks the generator (shared predicate — can't drift).
            output_block = _unsupported_output_array_method(
                m['args'], m['return'])
            if output_block:
                reasons.append(output_block)
            rows.append({
                'method' : m,
                'deps'   : ret_deps | arg_deps,
                'reasons': reasons,
            })
        out[fc] = rows
    return out


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
    """Parse libraryClasses.xml and return ``(registered, overrides, null_filled,
    absent_filled, internal_selectors)``.

    *registered* is the set of class tags inside ``<classes>``.

    *overrides* is a dict ``impl_name -> set(arg_names)`` listing constructor
    arguments for which the XML supplies a ``<constructor><argument
    name=… type=… module=…/></constructor>`` hint.  These hints resolve
    otherwise-unsupported ``class(*)`` constructor args, so the audit must
    treat the corresponding impls as ready rather than flagging the
    ``class(*) without override`` blocker.

    *null_filled* is a dict ``impl_name -> set(arg_names)`` listing
    constructor arguments tagged with ``value="null"``.  Those args are
    dropped from the wrapper entirely (a local null pointer is passed
    to the inner constructor), so the audit must skip the relevant
    blockers (procedure-pointer / unsupported ``class(*)`` etc.) for
    them.

    *absent_filled* is a dict ``impl_name -> set(arg_names)`` for args
    tagged ``value="absent"``: the wrapper drops them entirely from
    both the Python and the bind(c) signatures *and* from the inner
    call (the inner constructor must handle the absent case via its
    declared default).  Only valid for `optional` args.

    *internal_selectors* is a dict ``impl_name -> chosen_internal_name``
    extracted from ``<constructor internal="…"/>`` attributes.  Used
    to pick one of several Internal-suffixed constructors when an
    impl exposes more than one.
    """
    if not xml_path.exists():
        return set(), {}, {}, {}, {}
    classes_el = ET.parse(xml_path).getroot().find('classes')
    if classes_el is None:
        return set(), {}, {}, {}, {}
    registered         = {child.tag for child in classes_el}
    overrides          = defaultdict(set)
    null_filled        = defaultdict(set)
    absent_filled      = defaultdict(set)
    internal_selectors = {}
    for class_el in classes_el:
        for impl_el in class_el:
            ctor_el = impl_el.find('constructor')
            if ctor_el is None:
                continue
            chosen = ctor_el.get('internal')
            if chosen:
                internal_selectors[impl_el.tag] = chosen
            for arg_el in ctor_el.findall('argument'):
                arg_name = arg_el.get('name')
                if not arg_name:
                    continue
                overrides[impl_el.tag].add(arg_name)
                if arg_el.get('value') == 'null':
                    null_filled[impl_el.tag].add(arg_name)
                if arg_el.get('value') == 'absent':
                    absent_filled[impl_el.tag].add(arg_name)
    return (registered, dict(overrides),
            dict(null_filled), dict(absent_filled),
            internal_selectors)


def fmt_section(title):
    bar = '=' * len(title)
    return f"\n{bar}\n{title}\n{bar}\n"


def main():
    print("Discovering functionClasses…", file=sys.stderr)
    all_fcs, fc_files                                 = discover_function_classes(SOURCE_DIR)
    (registered, overrides,
     null_filled, absent_filled,
     internal_selectors)                              = load_registered_classes(LIB_XML)
    class_hierarchy                                   = build_type_hierarchy(SOURCE_DIR)
    print(f"  {len(all_fcs)} functionClass(es); "
          f"{len(registered)} currently registered.", file=sys.stderr)

    audit = {}
    for i, fc in enumerate(sorted(all_fcs), 1):
        print(f"\rParsing impls… {i}/{len(all_fcs)} ({fc[:40]:<40})",
              end='', file=sys.stderr)
        impls = []
        for path in sorted(fc_files.get(fc, [])):
            impls.extend(parse_impls_in_file(
                path, fc, internal_selectors=internal_selectors))
        audit[fc] = aggregate_class(impls, all_fcs, registered, overrides,
                                    class_hierarchy=class_hierarchy,
                                    null_filled=null_filled,
                                    absent_filled=absent_filled)
    print("\n", file=sys.stderr)

    # Buckets.
    by_status = defaultdict(list)
    for fc, info in audit.items():
        by_status[info['status']].append(fc)

    # Impl-level totals — a class can be READY at the class level (some
    # impl is ready) while still having other impls that the generator
    # drops; the impl-level numbers make those drops visible alongside
    # the class-level totals.  `dropped` counts impls (in registered
    # classes only — drops in unregistered classes don't affect the
    # built library).
    impl_total_concrete   = 0
    impl_ready            = 0
    impl_dropped_in_built = 0
    for fc, info in audit.items():
        for impl, deps, reasons in info.get('classified', []):
            impl_total_concrete += 1
            if not deps and not reasons:
                impl_ready += 1
            elif fc in registered:
                impl_dropped_in_built += 1

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
    print()
    print(f"  Total concrete impls:     {impl_total_concrete}")
    print(f"  Impl-level ready:         {impl_ready} "
          f"({100*impl_ready//impl_total_concrete}%)")
    print(f"  Dropped from BUILT lib:   {impl_dropped_in_built} "
          f"(impls in registered classes the generator skips)")

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

    # DROPPED IMPLS — impls of *registered* classes that the generator
    # nevertheless skips (because their constructor takes a type the
    # pipeline can't translate, or has unmet missing deps even though
    # some sibling impl is ready).  These are the impls that show up
    # as `caution: implementation '<name>' of class '<cls>' has
    # constructor argument …` warnings during a `make` build, and are
    # invisible to the class-level READY/MISSING-DEP/PIPELINE-BLOCKED
    # buckets above — the class itself still surfaces in Python via
    # its other impls, but the dropped ones simply aren't reachable.
    dropped = []
    for fc in sorted(registered):
        info = audit.get(fc, {})
        for impl, deps, reasons in info.get('classified', []):
            if reasons or deps:
                # Pick the most informative line — pipeline reasons
                # dominate; otherwise list missing deps.
                if reasons:
                    summary = reasons[0]
                else:
                    summary = "missing dep(s): " + ', '.join(sorted(deps))
                dropped.append((fc, impl['name'], summary))
    print(fmt_section(
        f"DROPPED IMPLS ({len(dropped)}) — registered classes whose"
        f" non-ready impls are skipped"))
    if dropped:
        # Width chosen so the longest current `<class> :: <impl>` label
        # fits without truncation; reasons line up in the second column.
        label_width = max(
            (len(f"  {fc} :: {impl_name}") for fc, impl_name, _ in dropped),
            default=0,
        )
        for fc, impl_name, summary in dropped:
            label = f"  {fc} :: {impl_name}"
            print(f"{label.ljust(label_width)}  {summary}")
    else:
        print("  (none — every concrete impl of every registered class is built)")

    # NO-CONCRETE-IMPLS — generally these are abstract base classes that
    # only matter as deps; surface them but don't dwell.
    no_impls = sorted(by_status['no-concrete-impls'])
    if no_impls:
        print(fmt_section(
            f"NO-CONCRETE-IMPLS ({len(no_impls)}) — abstract bases only"))
        print(f"  {', '.join(no_impls)}")

    # ------------------------------------------------------------------
    # Methods pass — classify every `<method>` directive on every
    # functionClass and report what the generator would drop.  See
    # :func:`classify_method_return`, :func:`aggregate_methods`, and
    # :func:`_is_out_of_scope_reason` for the predicates.
    # ------------------------------------------------------------------
    print("Discovering methods…", file=sys.stderr)
    methods_by_fc = discover_methods(SOURCE_DIR, all_fcs)
    method_audit  = aggregate_methods(
        methods_by_fc, all_fcs, registered,
        class_hierarchy=class_hierarchy)

    # Aggregates: counts by status, separated into in-scope vs out-of-scope
    # blockers so the report tracks the two backlogs independently.
    n_total           = 0
    n_ready           = 0
    n_blocked_inscope = 0
    n_blocked_oos     = 0
    n_missing_dep     = 0
    in_scope_buckets  = defaultdict(set)   # head -> set((fc, method))
    oos_buckets       = defaultdict(set)
    dropped_in_built  = []                 # rows for registered classes only
    for fc, rows in method_audit.items():
        for row in rows:
            n_total += 1
            mname = row['method']['name']
            if not row['reasons'] and not row['deps']:
                n_ready += 1
                continue
            if row['reasons']:
                # In-scope iff *every* reason is in-scope; mixed cases get
                # tagged in-scope so the in-scope worklist sees them too.
                any_in_scope = any(not _is_out_of_scope_reason(r)
                                   for r in row['reasons'])
                for r in row['reasons']:
                    head = re.split(r'\s+\(', r, maxsplit=1)[0]
                    if _is_out_of_scope_reason(r):
                        oos_buckets[head].add((fc, mname))
                    else:
                        in_scope_buckets[head].add((fc, mname))
                if any_in_scope:
                    n_blocked_inscope += 1
                else:
                    n_blocked_oos += 1
            else:
                # deps only — would be ready once the dep is registered.
                n_missing_dep += 1
            if fc in registered and (row['reasons'] or row['deps']):
                # Choose a single summary line, in-scope reasons preferred.
                summary = None
                for r in row['reasons']:
                    if not _is_out_of_scope_reason(r):
                        summary = r
                        break
                if summary is None:
                    summary = (row['reasons'][0] if row['reasons']
                               else 'missing dep(s): '
                                    + ', '.join(sorted(row['deps'])))
                dropped_in_built.append((fc, mname, summary,
                                         any(_is_out_of_scope_reason(r)
                                             for r in row['reasons'])
                                         and not any(not _is_out_of_scope_reason(r)
                                                     for r in row['reasons'])))

    n_total_classes_with_methods = len(method_audit)
    print(fmt_section("Methods — summary"))
    print(f"  Total method directives:  {n_total} "
          f"across {n_total_classes_with_methods} functionClasses")
    if n_total:
        print(f"  Ready:                    {n_ready} "
              f"({100*n_ready//n_total}%)")
        print(f"  Missing-dep:              {n_missing_dep}")
        print(f"  In-scope blocked:         {n_blocked_inscope}")
        print(f"  Out-of-scope blocked:     {n_blocked_oos} "
              f"(deferred — treeNode / other internal derived types /"
              f" non-fc class hierarchies)")

    # Top in-scope blockers — what the team can actually work on next.
    print(fmt_section(
        f"METHOD BLOCKERS — in-scope ({sum(len(v) for v in in_scope_buckets.values())}"
        f" method occurrences)"))
    if in_scope_buckets:
        for head, occs in sorted(in_scope_buckets.items(),
                                 key=lambda kv: -len(kv[1])):
            n_classes = len({fc for fc, _ in occs})
            print(f"  {len(occs):4d}  ({n_classes:3d} classes)  {head}")
    else:
        print("  (none — every method's in-scope path is unblocked)")

    # Same shape, for the deferred bucket so the size of the future
    # backlog is visible at a glance.
    print(fmt_section(
        f"METHOD BLOCKERS — out-of-scope (deferred)"))
    if oos_buckets:
        for head, occs in sorted(oos_buckets.items(),
                                 key=lambda kv: -len(kv[1])):
            n_classes = len({fc for fc, _ in occs})
            print(f"  {len(occs):4d}  ({n_classes:3d} classes)  {head}")
    else:
        print("  (none)")

    # Dropped methods on registered classes — the actual loss in today's
    # built library, mirroring the DROPPED IMPLS section above.  Separated
    # by in-scope vs out-of-scope so the actionable list isn't drowned by
    # the deferred backlog.
    dropped_in_scope = [r for r in dropped_in_built if not r[3]]
    dropped_oos      = [r for r in dropped_in_built if     r[3]]
    print(fmt_section(
        f"DROPPED METHODS — in-scope ({len(dropped_in_scope)}) on registered classes"))
    if dropped_in_scope:
        label_width = max(
            (len(f"  {fc} :: {m}") for fc, m, _, _ in dropped_in_scope),
            default=0,
        )
        for fc, mname, summary, _ in dropped_in_scope:
            label = f"  {fc} :: {mname}"
            print(f"{label.ljust(label_width)}  {summary}")
    else:
        print("  (none — every in-scope method on a registered class is built)")

    print(fmt_section(
        f"DROPPED METHODS — out-of-scope (deferred, {len(dropped_oos)}) on registered classes"))
    if dropped_oos:
        # Don't enumerate the whole list — too long to be useful.  Just
        # surface the per-class drop count so the report stays scannable.
        by_class = defaultdict(int)
        for fc, _, _, _ in dropped_oos:
            by_class[fc] += 1
        for fc in sorted(by_class, key=lambda k: -by_class[k]):
            print(f"  {fc:50s}  {by_class[fc]:3d} method(s)")
    else:
        print("  (none)")


if __name__ == '__main__':
    main()
