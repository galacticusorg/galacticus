# Wraps Fortran I/O statements in an `!$omp critical(gfortranInternalIO_)`
# section when $GALACTICUS_FCFLAGS contains `-DTHREADSAFEIO`, to work
# around gfortran PR 92836 (internal-file I/O is not thread-safe there).
# Also rewrites existing `gfortranInternalIO` / `FoX_DOM_Access` critical
# sections to use the same underlying `gfortranInternalIO_` lock name.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ThreadSafeIO.pm

import io
import os
import re
import sys


from build.fortran_utils                          import get_fortran_line
from Galacticus.Build.SourceTree                  import walk_tree
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


_IO_STATEMENTS = ('open', 'close', 'read', 'write')

_CRITICAL_OPEN_RE  = re.compile(
    r'^\s*!\$omp\s+critical\s*\(gfortranInternalIO(_)??\)\s*$', re.IGNORECASE)
_CRITICAL_CLOSE_RE = re.compile(
    r'^\s*!\$omp\s+end\s+critical\s*\(gfortranInternalIO(_)??\)\s*$', re.IGNORECASE)
_FOX_DOM_RE = re.compile(
    r'^\s*!\$omp\s+(end\s+)??critical\s*\(FoX_DOM_Access\)\s*$', re.IGNORECASE)
_CALL_FLUSH_RE = re.compile(
    r'^\s*(!\$)??\s*call\s+flush\s*\(', re.IGNORECASE)


_ACQUIRE_REPORT = (
    "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),"
    "' acquired the ''gfortranInternalIO'' lock'\n"
)
_RELEASE_REPORT = (
    "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),"
    "' released the ''gfortranInternalIO'' lock'\n"
)


def _enclosing_container(node):
    """Walk up to the first `function`/`subroutine`/`moduleProcedure`/
    `module`/`program` ancestor.  Mirrors addUse() at ThreadSafeIO.pm:131-148.
    """
    cursor = node
    while cursor is not None:
        ntype = cursor.get('type')
        if ntype in ('function', 'subroutine', 'moduleProcedure',
                     'module', 'program'):
            return cursor
        cursor = cursor.get('parent')
    return None


def _inject_omp_use(node):
    """Add `use OMP_Lib, only : OMP_Get_Thread_Num` to the enclosing
    container.  Mirrors addUse() body at ThreadSafeIO.pm:149-163.
    """
    container = _enclosing_container(node)
    if container is None:
        return
    add_uses(container, {
        'moduleUse': {
            'OMP_Lib': {
                'openMP':    True,
                'intrinsic': False,
                'only':      {'OMP_Get_Thread_Num': True},
            },
        },
        'moduleOrder': ['OMP_Lib'],
    })


def _is_io_line(processed_line):
    for stmt in _IO_STATEMENTS:
        if re.match(
                r'^\s*(!\$)??\s*' + stmt + r'\s*\(', processed_line):
            return True
    return bool(_CALL_FLUSH_RE.match(processed_line))


def _tree_is_error_module(tree):
    """Return True if any descendant is `module Error` — matches the early
    exit at ThreadSafeIO.pm:34-40 that avoids locking inside the error-
    reporting module (where deadlocks would occur).
    """
    for n in walk_tree(tree):
        if n.get('type') == 'module' and n.get('name') == 'Error':
            return True
    return False


def process_thread_safe_io(tree, options):
    """Mirrors Lock_IO() from ThreadSafeIO.pm."""
    flags = os.environ.get('GALACTICUS_FCFLAGS', '')
    if '-DTHREADSAFEIO' not in flags.split():
        return
    if _tree_is_error_module(tree):
        return

    report = os.environ.get('GALACTICUS_REPORT_THREADSAFEIO') == 'yes'

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        fh = io.StringIO(content)
        new_content = ''
        in_io       = False
        in_critical = False

        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                break
            ignore_line = False

            # Existing gfortranInternalIO critical section — rename the
            # explicit one to the underscore variant, and (when reporting)
            # emit acquire/release messages.
            m_open = _CRITICAL_OPEN_RE.match(processed_line)
            if m_open:
                is_explicit = m_open.group(1) is None
                processed_line = processed_line.replace(
                    '(gfortranInternalIO)', '(gfortranInternalIO_)')
                if report and is_explicit:
                    processed_line = processed_line + _ACQUIRE_REPORT
                    ignore_line = True
                    _inject_omp_use(node)
                in_critical = True

            m_close = _CRITICAL_CLOSE_RE.match(processed_line)
            if m_close:
                is_explicit = m_close.group(1) is None
                processed_line = processed_line.replace(
                    '(gfortranInternalIO)', '(gfortranInternalIO_)')
                if report and is_explicit:
                    processed_line = _RELEASE_REPORT + processed_line
                    ignore_line = True
                    _inject_omp_use(node)
                in_critical = False

            # FoX_DOM_Access critical sections are converted to the same
            # gfortranInternalIO_ lock so they also block internal I/O.
            m_fox = _FOX_DOM_RE.match(raw_line)
            if m_fox:
                raw_line = raw_line.replace(
                    'FoX_DOM_Access', 'gfortranInternalIO_')
                if report:
                    _inject_omp_use(node)
                    if m_fox.group(1):   # "end critical" form
                        raw_line = _RELEASE_REPORT + raw_line
                    else:
                        raw_line = raw_line + _ACQUIRE_REPORT

            is_io = _is_io_line(processed_line)

            if is_io and not in_io and not in_critical and not ignore_line:
                in_io = True
                new_content += "    !$omp critical(gfortranInternalIO_)\n"
                if report:
                    new_content += _ACQUIRE_REPORT
                    _inject_omp_use(node)

            if (not is_io) and in_io and not in_critical and not ignore_line:
                in_io = False
                if report:
                    new_content += _RELEASE_REPORT
                    _inject_omp_use(node)
                new_content += "    !$omp end critical(gfortranInternalIO_)\n"

            new_content += raw_line

        if in_io and not in_critical:
            if report:
                new_content += _RELEASE_REPORT
                _inject_omp_use(node)
            new_content += "    !$omp end critical(gfortranInternalIO_)\n"

        node['content'] = new_content


register_process(
    'threadSafeIO',
    process_thread_safe_io,
    before=['metaPropertyDatabase', 'stateStore', 'functionClass', 'stateStorable'],
)
