#!/usr/bin/env python3
# Scan source code for "!![…!!]" XML directives and dispatch the matching
# generator.  Andrew Benson (ported to Python 2026)
#
# The top-level driver that the build system invokes once per
# `<directive>.xml` control file: each include directive's Makefile rule
# (written into Makefile_Directives by codeDirectivesParse.py) expands to
# `./scripts/build/buildCode.py <install_dir> <directive_xml>`.
#
# The driver scans every Fortran source file containing the named
# directive, accumulates each `<directive>` block's parsed XML body together
# with the surrounding module name and source-file metadata, then dispatches
# a `validate → parse → generate` pipeline registered by
# `Galacticus::Build::Hooks` for the build's `type` field.  Currently the
# only registered handler is `component`, which is owned by
# python/Galacticus/Build/Components.

import argparse
import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.FortranUtils                      import get_fortran_line
from Galacticus.Build.FileChanges                       import update as file_changes_update
from Galacticus.Build.ScanCache                         import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
    prune           as _prune_cache,
)
from Fortran.Utils                            import UNIT_OPENERS, UNIT_CLOSERS
from XML.Utils                                import xml_to_dict
from Galacticus.Build                         import Hooks
from Galacticus.Build                         import SourceTree
from Galacticus.Build.SourceTree.Process      import process_tree
import Galacticus.Build.Components            # registers the `component` hook on import
from Galacticus._logging                      import configure_default as _configure_default

# Show INFO-level diagnostic output from the library modules (mirrors the
# verbose `print()`-driven output of the Perl-era driver).
_configure_default()

# Register every source-tree process hook (the full set — see Process/all.py).
import Galacticus.Build.SourceTree.Process.all  # noqa: F401


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_END_MARKER_RE   = re.compile(r'^\s*!!\]')
_START_MARKER_RE = re.compile(r'^\s*!!\[')
_INCLUDE_RE      = re.compile(r"^\s*include\s*['\"]([^'\"]+)['\"]\s*$")
# An XML directive opening tag (single-line or first-line of a multi-line
# block) — `<foo …>` or `<foo …/>`.
_DIRECTIVE_OPEN_RE = re.compile(r'^\s*(<\s*([a-zA-Z]+)\b[^>]*>)\s*$')

# `XMLin($xml, ForceArray => ["data","property","binding"])` — the Perl
# driver always lifts these tags into lists so downstream consumers can
# iterate uniformly.
_FORCE_ARRAY = {'data', 'property', 'binding'}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Scan Galacticus Fortran sources for the directive "
                    "named in <xml-file> and dispatch its registered "
                    "build handler.",
    )
    parser.add_argument('source_directory',
                        help="Top-level source directory to scan.")
    parser.add_argument('xml_file',
                        help="Build-control XML file describing one "
                             "include directive (`type`, `directive`, "
                             "`fileName`, …).")
    args = parser.parse_args(argv)

    source_directory = args.source_directory
    build_path       = os.environ['BUILDPATH']

    # Optional: list of files known to contain each directive.
    locations = _xml_load_optional(
        os.path.join(build_path, 'directiveLocations.xml')
    )
    # Required: function-class registry (used to attribute orphan
    # `<…Class>` directives to the right module).
    storables_path = os.path.join(build_path, 'stateStorables.xml')
    if not os.path.exists(storables_path):
        sys.exit(f"file '{storables_path}' is missing")
    storables = _xml_load(storables_path)

    # The build-control XML.
    build = _xml_load(args.xml_file)
    build['moduleName']      = ""
    build['currentFileName'] = ""
    build['codeType']        = "fortran"

    # Per-file directive cache.  Persists between runs so that files whose
    # mtime is older than the cache get skipped.  Mirrors Perl's `Storable`
    # blob at `<fileName>.blob`.
    blob_path           = build['fileName'] + '.blob'
    directives_per_file, update_time = _load_cache(blob_path)
    have_per_file       = update_time is not None

    # Decide which files to scan.
    file_names_to_scan = _scan_targets(
        build['directive'], locations, source_directory
    )

    # Walk every file, populating the cache.
    for file_name in file_names_to_scan:
        _scan_file(
            file_name, build, source_directory,
            directives_per_file, have_per_file, update_time,
            storables,
        )

    # Drop cache entries for files no longer among the scan targets (deleted,
    # or the directive was removed), or their directives would keep feeding
    # the generated code below forever.
    _prune_cache(directives_per_file,
                 {_file_identifier(f) for f in file_names_to_scan})

    # Persist the cache.
    with open(blob_path, 'wb') as fh:
        pickle.dump(directives_per_file, fh)

    # Dispatch validate / parse for every cached directive, then run the
    # whole-build `generate` step.
    handler = Hooks.module_hooks.get(build['type'])
    if handler is None:
        sys.exit(f"buildCode.py: failed to find a function to parse "
                 f"'{build['type']}' action")

    for entry in directives_per_file.values():
        for directive in entry.get('directives', []):
            build['currentDocument'  ] = directive['directive'      ]
            build['moduleName'       ] = directive['moduleName'     ]
            build['currentFileName'  ] = directive['currentFileName']
            build['codeType'         ] = directive['codeType'       ]
            if 'validate' in handler:
                handler['validate'](directive['xmlCode'], directive['fileName'])
            handler['parse'](build)

    if 'generate' not in handler:
        sys.exit(f"buildCode.py: failed to find a function to generate "
                 f"'{build['type']}' action")
    handler['generate'](build)

    # Emit the result.  `.Inc` outputs are run through the SourceTree
    # processor so that any directives the generator embedded in its own
    # output (e.g. `<allocate>`) are expanded.  Other outputs are dumped
    # verbatim.
    tmp_path = build['fileName'] + '.tmp'
    with open(tmp_path, 'w') as out:
        if build['fileName'].endswith('.Inc'):
            # `annotate=True` emits the `!--> <origLine> <outLine> "<source>"`
            # line-number mapping comments that the rest of the Galacticus
            # build pipeline expects ahead of each node's serialised content.
            # Perl's Serialize() defaults to annotated; the Python port
            # defaults to off so existing in-process callers (Generics,
            # FunctionClass, …) get pure Fortran back, but `buildCode.py`
            # is producing an `.Inc` for downstream stages and needs the
            # markers.
            out.write(SourceTree.serialize(
                process_tree(
                    SourceTree.parse_code(build['content'], build['fileName'])
                ),
                annotate=True,
            ))
        else:
            out.write(build['content'])

    file_changes_update(build['fileName'], tmp_path, prove_update=True)


# ---------------------------------------------------------------------------
# Build-XML loading
# ---------------------------------------------------------------------------

def _xml_load_optional(path):
    """Load `path` as a top-level dict via `xml_to_dict`, or None."""
    return _xml_load(path) if os.path.exists(path) else None


def _xml_load(path):
    """Load `path` and return its root element as a dict.

    Matches Perl `XMLin($file, KeyAttr => [])`: same defaults as
    `xml_to_dict` (no `keyed_tags`, single same-tag children become a dict,
    multiples become a list).
    """
    tree = ET.parse(path)
    return xml_to_dict(tree.getroot())


# ---------------------------------------------------------------------------
# Scan-target selection
# ---------------------------------------------------------------------------

def _scan_targets(directive, locations, source_directory):
    """Return the list of files to scan for a given directive name, from the
    pre-computed `directiveLocations.xml` map.

    The map is required: the Makefile rules that invoke this script are
    themselves written by codeDirectivesParse.py, which generates
    `directiveLocations.xml` in the same run, so the file always exists in a
    build. A missing map means this script was invoked out of sequence.
    """
    if locations is None:
        sys.exit(
            "buildCode.py: directiveLocations.xml not found — run "
            "./scripts/build/codeDirectivesParse.py first (in a build this "
            "is done by the $(BUILDPATH)/directiveCatalogs.stamp rule)."
        )
    if directive not in locations:
        return []
    files = locations[directive].get('file')
    if files is None:
        return []
    return files if isinstance(files, list) else [files]


# ---------------------------------------------------------------------------
# Per-file directive scanner
# ---------------------------------------------------------------------------

def _scan_file(file_name, build, source_directory,
               directives_per_file, have_per_file, update_time, storables):
    """Walk `file_name` (and any `include`d files), collect every matching
    `<{directive}>` block, and stash the result in
    `directives_per_file[file_id]`.

    The work mirrors buildCode.pl's main `foreach my $fileName (@fileNamesToScan)`
    loop: each include file is processed in the module-name context of its
    parent so that a directive appearing inside an inclusion gets the
    correct `moduleName`.  Files whose recorded mtime predates the cache are
    skipped wholesale.
    """
    file_identifier = _file_identifier(file_name)

    # mtime-based skip:  if every file we previously scanned for this id is
    # older than the cache, nothing has changed.
    if have_per_file and file_identifier in directives_per_file:
        recorded = directives_per_file[file_identifier].get('files', [])
        # Perl tests `-M $_ < $updateTime` — i.e. *any* recorded file
        # newer than the cache forces a rescan.
        if not any(
            os.path.exists(p) and os.path.getmtime(p) > update_time
            for p in recorded
        ):
            return
        del directives_per_file[file_identifier]

    directives_per_file[file_identifier] = {
        'files':      [file_name],
        'directives': [],
    }

    build['currentFileName'] = file_name
    build['codeType'       ] = 'c' if re.search(r'\.c(pp)?$', file_name) else 'fortran'

    file_stack = [{'name': file_name, 'position': -1, 'in_xml': False}]
    while file_stack:
        frame = file_stack.pop()
        with open(frame['name'], 'r', errors='replace') as fh:
            if frame['position'] != -1:
                fh.seek(frame['position'])
            include_pushed = _process_until_include_or_eof(
                fh, frame, build, source_directory,
                directives_per_file[file_identifier], storables,
            )
            if include_pushed:
                # Push the parent back, then the include, so the include is
                # processed next.
                include_file, parent_frame = include_pushed
                file_stack.append(parent_frame)
                file_stack.append({
                    'name': include_file,
                    'position': -1,
                    'in_xml': False,
                })
                directives_per_file[file_identifier]['files'].append(include_file)


def _process_until_include_or_eof(fh, frame, build, source_directory,
                                  per_file_entry, storables):
    """Read lines from `fh`, accumulating directive bodies until either
    EOF is reached or an `include` statement requires re-entry.

    On encountering an include we record the current file position in
    `frame`, return `(include_path, frame)`, and let the caller push the
    parent + child frames onto the file stack.  Returning `None` means
    we exhausted the file.
    """
    while True:
        if build['codeType'] == 'fortran':
            raw_line, processed_line, _ = get_fortran_line(fh)
        else:
            # For C/C++ the raw and processed forms are the same single line;
            # `get_fortran_line`'s continuation-joining does not apply. Read
            # ONE line and use it for both — mirrors `_consume_until_close`
            # below and the Perl original. (A previous
            # `(fh.readline(), fh.readline(), '')` here read two lines per
            # iteration, so `raw_line`/`processed_line` were different,
            # consecutive lines and every other line was skipped.)
            line = fh.readline()
            raw_line, processed_line = line, line
        if not raw_line and not processed_line:
            return None

        # Detect end of an XML section *before* reading content for this line.
        if _END_MARKER_RE.match(raw_line):
            frame['in_xml'] = False

        include_file = None

        # Fortran-specific bookkeeping:  module entry/exit and `include`
        # detection.
        if build['codeType'] == 'fortran':
            m = _INCLUDE_RE.match(processed_line)
            if m:
                include_file = re.sub(
                    r'\.inc$', '.Inc',
                    os.path.join(source_directory, 'source', m.group(1)),
                )
            m = UNIT_OPENERS['module']['regex'].search(processed_line)
            if m:
                build['moduleName'] = m.group(
                    UNIT_OPENERS['module']['unit_name'] + 1
                )
            elif UNIT_CLOSERS['module'].search(processed_line):
                build['moduleName'] = ""

        if frame['in_xml']:
            stripped = re.sub(r'^\s*!<\s*', '', raw_line)
            m = _DIRECTIVE_OPEN_RE.match(stripped)
            if m:
                xml_code = m.group(1) + "\n"
                xml_tag  = m.group(2)
                # Read on until the matching close tag, unless this is a
                # self-closing element.
                if not re.search(r'/\s*>', xml_code):
                    xml_code, include_in_block = _consume_until_close(
                        fh, build, xml_tag, frame, source_directory, xml_code,
                    )
                    if include_in_block is not None and include_file is None:
                        include_file = include_in_block

                if xml_tag == build['directive']:
                    directive = _parse_directive(
                        xml_code, frame, build,
                    )
                    if directive is not None:
                        per_file_entry['directives'].append(directive)
                else:
                    # functionClass instance directives report orphan
                    # `<…Class>` tags so we can still attribute them to the
                    # module they belong to during a future scan.
                    fc_key = xml_tag + 'Class'
                    function_classes = (storables or {}).get('functionClasses') or {}
                    if fc_key in function_classes:
                        build['moduleName'] = function_classes[fc_key].get(
                            'module', build['moduleName']
                        )

        # Finally, detect the start of an XML section.
        if _START_MARKER_RE.match(raw_line):
            frame['in_xml'] = True

        if include_file is not None and os.path.exists(include_file):
            frame['position'] = fh.tell()
            return include_file, frame


def _consume_until_close(fh, build, xml_tag, frame, source_directory,
                         xml_code):
    """Read further lines from `fh` until the matching `</xml_tag>` (or
    EOF).  Returns `(xml_code, include_file_or_None)`.
    """
    include_file = None
    close_re = re.compile(r'</' + re.escape(xml_tag) + r'>')

    build_path = os.environ['BUILDPATH']

    while True:
        if build['codeType'] == 'fortran':
            raw, processed, _ = get_fortran_line(fh)
        else:
            line = fh.readline()
            raw, processed = line, line
        if not raw:
            break

        next_line = re.sub(r'^\s*!<\s*', '', raw)

        if _END_MARKER_RE.match(next_line):
            frame['in_xml'] = False

        if build['codeType'] == 'fortran':
            m = _INCLUDE_RE.match(processed)
            if m and include_file is None:
                include_file = re.sub(
                    r'\.inc$', '.p.Inc',
                    os.path.join(source_directory, build_path, m.group(1)),
                )

        next_line = re.sub(r'^\s*', '', next_line)
        xml_code += next_line

        if _START_MARKER_RE.match(next_line):
            frame['in_xml'] = True
        if close_re.search(next_line):
            break

    return xml_code, include_file


def _parse_directive(xml_code, frame, build):
    """Parse one `xml_code` block and return a directive dict.

    Mirrors Perl's `XMLin($xml_code, ForceArray => ["data","property",
    "binding"])`.  Failure is fatal — Perl `die`s in the same spot.
    """
    try:
        elem = ET.fromstring(xml_code)
    except ET.ParseError as exc:
        sys.exit(
            f"buildCode.py: failed in {frame['name']} with message:\n{exc}"
        )

    body = xml_to_dict(elem, force_array=_FORCE_ARRAY)
    if not isinstance(body, dict):
        body = {'content': body}

    return {
        'fileName':        build['currentFileName'],
        'xmlCode':         xml_code,
        'directive':       body,
        'moduleName':      build['moduleName'     ],
        'currentFileName': build['currentFileName'],
        'codeType':        build['codeType'       ],
    }


# ---------------------------------------------------------------------------
# Cache I/O
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    main()
