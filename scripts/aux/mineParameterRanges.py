#!/usr/bin/env python3
"""Mine obvious run-time range checks and (optionally) seed `<minimum>`/
`<maximum>` annotations on the corresponding `<inputParameter>` directives.

Heuristic and review-oriented.  It finds single-line error guards of the form

    if (<var> <op> <literal> [.or. <var> <op> <literal> ...]) ... call Error_Report

in an implementation's source file, where `<var>` is a numeric input parameter
of that implementation.  An error condition inverts to a valid bound:

    var <  L  invalid -> minimum L (inclusive)     var >  L  invalid -> maximum L (inclusive)
    var <= L  invalid -> minimum L (exclusive)     var >= L  invalid -> maximum L (exclusive)

Conditions containing `.and.` are skipped (not a simple per-variable bound), as
are clauses comparing to a non-literal (e.g. another variable).

Usage:
    mineParameterRanges.py <repoRoot> [--catalog <path>] [--apply]

Without `--apply` it prints a report of proposed constraints.  With `--apply` it
inserts the annotations into the directive source files (skipping any parameter
whose directive cannot be located unambiguously, or whose bound conflicts with
an existing annotation).  Re-run the catalog build and parameter-file validator
afterwards to confirm.

Andrew Benson (2026).
"""

import collections
import json
import os
import re
import sys

_GUARD_RE = re.compile(r"if\s*\(([^()]*)\)\s*call\s+Error_Report", re.IGNORECASE)
_CLAUSE_RE = re.compile(
    r"([A-Za-z_]\w*)\s*(<=|>=|<|>)\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:d[+-]?\d+)?)",
    re.IGNORECASE)
_INPUT_PARAMETER_OPEN = re.compile(r"<inputParameter\b")
_INPUT_PARAMETER_CLOSE = re.compile(r"</inputParameter>")
_NAME_RE = re.compile(r"<name>\s*([^<\s]+)\s*</name>")


def _invert(op, literal):
    """Map an error-condition comparison to a (kind, value, inclusive) bound."""
    return {
        '<':  ('minimum', literal, True),
        '<=': ('minimum', literal, False),
        '>':  ('maximum', literal, True),
        '>=': ('maximum', literal, False),
    }[op]


def _normalize_literal(literal):
    """Strip a Fortran kind/exponent suffix for the annotation value
    (`0.0d0` -> `0.0`, `1.0d0` -> `1.0`); plain integers/decimals pass through."""
    return re.sub(r'[dD][+-]?\d+$', '', literal)


def mine(catalog):
    """Return ``{(impl_type, paramName): {'minimum': (value, inclusive),
    'maximum': (...), 'evidence': [...]}}`` for the source tree described by
    `catalog`."""
    by_file = collections.defaultdict(list)
    for impl_type, impl in catalog["implementations"].items():
        numeric = {}
        for parameter in impl["parameters"]:
            if parameter["type"] in ("real", "integer") and not parameter.get("constraints"):
                for key in (parameter.get("name"), parameter.get("variable")):
                    if key:
                        numeric[key] = parameter["name"]
        if numeric and impl.get("sourceFile"):
            by_file[impl["sourceFile"]].append((impl_type, numeric))

    proposals = collections.defaultdict(lambda: {"evidence": []})
    for source_file, impls in by_file.items():
        path = os.path.join(REPO, "source", source_file)
        if not os.path.exists(path):
            continue
        raw = re.sub(r"&\s*\n\s*&?", " ", open(path, errors="replace").read())
        for line in raw.split("\n"):
            guard = _GUARD_RE.search(line)
            if not guard or re.search(r"\.and\.", guard.group(1), re.IGNORECASE):
                continue
            for var, op, literal in _CLAUSE_RE.findall(guard.group(1)):
                kind, value, inclusive = _invert(op, literal)
                for impl_type, numeric in impls:
                    if var in numeric:
                        key = (impl_type, numeric[var])
                        proposals[key][kind] = (value, inclusive)
                        proposals[key]["evidence"].append(line.strip())
    return proposals


def _print_report(proposals):
    print(f"proposals: {len(proposals)} parameters across "
          f"{len(set(k[0] for k in proposals))} implementations\n")
    for (impl_type, param), data in sorted(proposals.items()):
        bits = []
        for kind in ("minimum", "maximum"):
            if kind in data:
                value, inclusive = data[kind]
                bits.append(f"{kind}={value}{'' if inclusive else ' (excl)'}")
        print(f"{impl_type}::{param}  ->  {', '.join(bits)}")
        for evidence in dict.fromkeys(data["evidence"]):
            print(f"      | {evidence[:110]}")


def _format_bound(tag, value, inclusive, indent):
    attr = "" if inclusive else ' inclusive="false"'
    return f"{indent}<{tag}{attr}>{_normalize_literal(value)}</{tag}>\n"


def _apply(catalog, proposals):
    """Insert the proposed bounds into the directive source files."""
    # Group by file -> {paramName: {minimum/maximum}} (catalog gives the file).
    source_file_of = {t: impl.get("sourceFile")
                      for t, impl in catalog["implementations"].items()}
    by_file = collections.defaultdict(dict)
    for (impl_type, param), data in proposals.items():
        source_file = source_file_of.get(impl_type)
        if not source_file:
            continue
        bounds = {k: data[k] for k in ("minimum", "maximum") if k in data}
        existing = by_file[source_file].get(param)
        if existing is not None and existing != bounds:
            print(f"  [skip] {source_file}::{param}: conflicting bounds across "
                  f"implementations", file=sys.stderr)
            by_file[source_file][param] = None
            continue
        by_file[source_file][param] = bounds

    applied = skipped = 0
    for source_file, params in by_file.items():
        path = os.path.join(REPO, "source", source_file)
        lines = open(path, errors="replace").read().split("\n")
        # Locate each inputParameter block's name + closing-line index + indent.
        blocks = {}
        in_block = False
        name = indent = None
        for i, line in enumerate(lines):
            if _INPUT_PARAMETER_OPEN.search(line):
                in_block, name, indent = True, None, None
            if in_block:
                match = _NAME_RE.search(line)
                if match and name is None:
                    name = match.group(1)
                    indent = re.match(r"\s*", line).group(0)
            if _INPUT_PARAMETER_CLOSE.search(line) and in_block:
                if name is not None:
                    blocks.setdefault(name, []).append((i, indent))
                in_block = False
        # Queue insertions; skip if the directive is missing or ambiguous.
        insertions = []
        for param, bounds in params.items():
            if not bounds:
                skipped += 1
                continue
            located = blocks.get(param)
            if not located or len(located) != 1:
                print(f"  [skip] {source_file}::{param}: "
                      f"{'no' if not located else 'ambiguous'} inputParameter "
                      f"directive", file=sys.stderr)
                skipped += 1
                continue
            close_index, block_indent = located[0]
            text = ""
            for kind, tag in (("minimum", "minimum"), ("maximum", "maximum")):
                if kind in bounds:
                    value, inclusive = bounds[kind]
                    text += _format_bound(tag, value, inclusive, block_indent)
            insertions.append((close_index, text))
            applied += 1
        for close_index, text in sorted(insertions, reverse=True):
            lines[close_index] = text + lines[close_index]
        if insertions:
            open(path, "w").write("\n".join(lines))
    print(f"\napplied {applied} parameters; skipped {skipped}.", file=sys.stderr)


def main(argv):
    global REPO
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = {a for a in argv[1:] if a.startswith("--")}
    catalog_path = None
    if "--catalog" in argv:
        catalog_path = argv[argv.index("--catalog") + 1]
        args = [a for a in args if a != catalog_path]
    if not args:
        print("Usage: mineParameterRanges.py <repoRoot> [--catalog <path>] [--apply]",
              file=sys.stderr)
        return 2
    REPO = args[0]
    catalog_path = catalog_path or os.path.join(REPO, "work/build/parameters.catalog.json")
    catalog = json.load(open(catalog_path))
    proposals = mine(catalog)
    if "--apply" in flags:
        _apply(catalog, proposals)
    else:
        _print_report(proposals)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
