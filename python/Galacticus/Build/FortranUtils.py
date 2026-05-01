# Utility functions for reading and parsing Fortran source code.
# Andrew Benson (ported to Python 2026)
"""Utility functions for reading and parsing Fortran source code.

Provides helpers for iterating over logical Fortran lines (handling
continuation lines and stripping comments) and for locating the start of
inline comments, including correct treatment of string literals and
OpenMP sentinels.
"""
from __future__ import annotations

import os
import re

from typing import IO

__all__ = [
    'get_fortran_line', 'extract_bracketed', 'extract_variables',
    'get_matching_lines', 'read_file',
]


def get_fortran_line(file_obj: IO[str]) -> tuple[str, str, str]:
    """Read one logical Fortran line from file_obj, handling & continuation lines.

    Returns a 3-tuple: (raw_line, processed_line, buffered_comments)
      raw_line          — the raw text as it appears in the file (all physical lines joined)
      processed_line    — continuations joined, comments stripped
      buffered_comments — any comment text that was stripped
    """
    processed_full_line = False
    raw_line           = ""
    processed_line     = ""
    buffered_comments  = ""
    first_line         = True

    while not processed_full_line:
        line = file_obj.readline()
        if not line:
            # EOF — return whatever we have
            break

        tmp_line = line

        # Locate the start of a comment in this physical line.
        # Comments begin with '!' except when inside string literals or {braces},
        # and not when '!' is followed by '$' (OpenMP directive).
        if '!' not in tmp_line:
            comment_position = -1
        elif _fast_comment_check(tmp_line):
            # The regex checks for '!' before any quote or brace, not followed by '$'.
            m = re.match(r"^([^'\"{\n]+?)![^\$]", tmp_line)
            if m:
                comment_position = len(m.group(1))
            else:
                comment_position = _find_comment_position(tmp_line)
        else:
            comment_position = _find_comment_position(tmp_line)

        raw_line += line
        # Strip trailing newline from processed_line before appending the next chunk.
        if processed_line.endswith('\n'):
            processed_line = processed_line[:-1]

        if comment_position == -1:
            code_part = line
        else:
            code_part = line[:comment_position] + '\n'
            comment_text = line[comment_position + 1:]
            buffered_comments += comment_text.rstrip('\n')

        # Strip leading continuation marker from continuation lines.
        code_part = re.sub(r'^\s*&\s*', '', code_part)
        code_part = re.sub(r'^\s*!\$\s+&\s*', '', code_part)

        if processed_line:
            processed_line += ' '
        processed_line += code_part

        # If the processed line ends with '&', it continues onto the next line.
        if re.search(r'&\s*$', processed_line):
            processed_line = re.sub(r'\s*&\s*$', '', processed_line)
        elif not first_line and (line.startswith('#') or re.match(r'^\s*![^\$]', line)):
            # Preprocessor directive or plain comment inside continuation — concatenate.
            pass
        else:
            processed_full_line = True

        first_line = False

    return raw_line, processed_line, buffered_comments


def _fast_comment_check(line: str) -> bool:
    """Return True if the line may need the full comment-position scan."""
    return bool(re.match(r"^([^'\"{\n]+?)![^\$]", line))


def _find_comment_position(line: str) -> int:
    """Scan through line character by character to find the comment start position."""
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if c == '!':
            # Check it's not an OpenMP sentinel '!$'
            if i + 1 < n and line[i + 1] == '$':
                i += 2
                continue
            return i
        elif c == '{':
            # Braces can be nested — skip to matching '}'
            depth = 0
            while i < n:
                if line[i] == '{':
                    depth += 1
                elif line[i] == '}':
                    depth -= 1
                    if depth == 0:
                        i += 1
                        break
                i += 1
        elif c in ('"', "'"):
            # Skip quoted string
            quote = c
            i += 1
            while i < n:
                if line[i] == quote:
                    if i + 1 < n and line[i + 1] == quote:
                        i += 2  # escaped quote
                        continue
                    i += 1
                    break
                i += 1
        else:
            i += 1
    return -1


def extract_bracketed(text: str, brackets: str = "()") -> tuple[str | None, str, str | None]:
    """Find and extract the first balanced bracket pair in text.

    Port of Perl Text::Balanced::extract_bracketed.

    Returns (extracted, remainder, prefix) where extracted includes the bracket
    chars themselves, or (None, text, None) if no balanced pair is found.

    brackets: "()" for parentheses, "[]" for square brackets, "()[]" for either
    (whichever opening char appears first in text is used).
    """
    close_map  = {'(': ')', '[': ']'}
    open_chars = [brackets[i] for i in range(0, len(brackets), 2) if brackets[i] in close_map]

    best_pos     = -1
    best_open: str | None = None
    for oc in open_chars:
        pos = text.find(oc)
        if pos != -1 and (best_pos == -1 or pos < best_pos):
            best_pos  = pos
            best_open = oc

    if best_pos == -1 or best_open is None:
        return None, text, None

    close_char = close_map[best_open]
    prefix     = text[:best_pos]
    depth      = 0
    i          = best_pos
    while i < len(text):
        if text[i] == best_open:
            depth += 1
        elif text[i] == close_char:
            depth -= 1
            if depth == 0:
                return text[best_pos:i + 1], text[i + 1:], prefix
        i += 1

    return None, text, None


def extract_variables(variable_list: str, lower_case: bool = True,
                      keep_qualifiers: bool = False, remove_spaces: bool = True) -> list[str]:
    """Given the post-'::' section of a Fortran variable declaration, return a list of names.

    Port of Fortran::Utils::Extract_Variables (perl/Fortran/Utils.pm).
    """
    if variable_list is None:
        return []

    if '::' in variable_list:
        raise ValueError(
            f"extract_variables: variable list '{variable_list}' contains '::'"
            " — regex matching likely failed"
        )

    if lower_case:
        variable_list = variable_list.lower()

    if remove_spaces:
        variable_list = re.sub(r'\s', '', variable_list)
    else:
        variable_list = variable_list.rstrip()

    if keep_qualifiers:
        variable_list = variable_list.replace('*', '%%ASTERISK%%')
    else:
        variable_list = variable_list.replace('*', '')

    # Remove (or encode) text within balanced () and [] pairs.
    iteration = 0
    while '(' in variable_list or '[' in variable_list:
        iteration += 1
        if iteration > 10000:
            raise RuntimeError(
                f"extract_variables: maximum iterations exceeded for input: '{variable_list}'"
            )
        extracted, remainder, prefix = extract_bracketed(variable_list, "()[]")
        if extracted is None:
            break
        # `extract_bracketed` returns prefix as `str` whenever `extracted` is
        # not None — they're correlated; narrow for the type checker.
        assert prefix is not None
        if keep_qualifiers:
            encoded = (extracted
                       .replace('(', '%%OPEN%%').replace(')', '%%CLOSE%%')
                       .replace('[', '%%OPENSQ%%').replace(']', '%%CLOSESQ%%')
                       .replace(',', '%%COMMA%%'))
            variable_list = prefix + encoded + remainder
        else:
            if prefix is None:
                raise ValueError(
                    f'extract_variables: failed to find prefix in "{variable_list}"'
                )
            variable_list = prefix + remainder

    # Remove initialisation / association expressions.
    if not keep_qualifiers:
        variable_list = re.sub(r'=[^,]*(,|$)', lambda m: m.group(1), variable_list)

    variables = [v.strip() for v in variable_list.split(',') if v.strip()]

    if keep_qualifiers:
        result = []
        for v in variables:
            v = (v.replace('%%OPEN%%', '(').replace('%%CLOSE%%', ')')
                  .replace('%%OPENSQ%%', '[').replace('%%CLOSESQ%%', ']')
                  .replace('%%COMMA%%', ',').replace('%%ASTERISK%%', '*'))
            result.append(v)
        return result

    return variables


# ---------------------------------------------------------------------------
# Whole-file readers
# ---------------------------------------------------------------------------

_processed_files_cache = {}


def get_matching_lines(file_name: str, regex: re.Pattern[str]) -> list[dict]:
    """Return every processed Fortran line in `file_name` that matches `regex`.

    Mirrors Perl Fortran::Utils::Get_Matching_Lines.  Each returned entry is
    a dict with keys 'line' (the processed line) and 'submatches' (the list
    of capture groups from the match).  Processed-line lists are cached per
    file across calls, matching the Perl module's `%processedFiles` cache.
    """
    if isinstance(regex, str):
        regex = re.compile(regex)

    if file_name not in _processed_files_cache:
        lines = []
        with open(file_name, 'r', errors='replace') as fh:
            while True:
                raw, processed, _ = get_fortran_line(fh)
                if not raw and not processed:
                    break
                lines.append(processed)
        _processed_files_cache[file_name] = lines

    matches = []
    for processed in _processed_files_cache[file_name]:
        m = regex.search(processed)
        if m:
            matches.append({'line': processed, 'submatches': list(m.groups())})
    return matches


def read_file(file_name: str, *, state: str = 'raw', follow_includes: bool = False,
              include_locations: list[str] | None = None,
              include_files_excluded: set[str] | None = None,
              strip_regex: re.Pattern[str] | str | None = None,
              strip_leading: bool = False, strip_trailing: bool = False,
              strip_empty: bool = False) -> str:
    """Return the (optionally preprocessed) text of `file_name`.

    Mirrors Perl Fortran::Utils::read_file.  Follows `include '…'` statements
    when `follow_includes=True`, searching `include_locations` with `.inc` /
    `.Inc` suffixes.  `state` is one of `'raw'`, `'processed'`, `'comments'`.
    """
    if include_locations is None:
        include_locations = []
    if include_files_excluded is None:
        include_files_excluded = set()
    if strip_regex is not None and isinstance(strip_regex, str):
        strip_regex = re.compile(strip_regex)

    # Perl prepends an empty string to the list of include locations so that
    # the bare `<dir>/<name>` path is tried first.
    all_include_locations = [''] + list(include_locations)

    file_names     = [file_name]
    file_positions = [-1]
    code_buffer    = []

    while file_names:
        with open(file_names[0], 'r', errors='replace') as fh:
            if file_positions[0] != -1:
                fh.seek(file_positions[0])

            include_pushed = False
            while True:
                raw_line, processed_line, buffered_comments = get_fortran_line(fh)
                if not raw_line and not processed_line:
                    break  # EOF

                # Detect `include '<name>'` and recurse.
                if follow_includes:
                    m = re.match(r"^\s*include\s*['\"]([^'\"]+)['\"]\s*$",
                                 processed_line)
                    if m:
                        include_leaf = m.group(1)
                        if include_leaf not in include_files_excluded:
                            # Perl: `$fileNames[0] =~ s/\/[^\/]+$/\//` — strip basename,
                            # leaving the directory with a trailing slash (or the
                            # original path unchanged if there was no slash).
                            cur_dir = re.sub(r'/[^/]+$', '/', file_names[0])
                            if cur_dir == file_names[0]:
                                cur_dir = ''
                            include_file = None
                            for suffix in ('.inc', '.Inc'):
                                for loc in all_include_locations:
                                    candidate = re.sub(
                                        r'\.inc$',
                                        suffix,
                                        cur_dir + loc + '/' + include_leaf,
                                    )
                                    if os.path.exists(candidate):
                                        include_file = candidate
                                        break
                                if include_file is not None:
                                    break
                            if include_file is not None:
                                file_positions[0] = fh.tell()
                                file_names.insert(0, include_file)
                                file_positions.insert(0, -1)
                                include_pushed = True
                                break

                if state == 'raw':
                    line = raw_line
                elif state == 'processed':
                    line = processed_line
                elif state == 'comments':
                    line = buffered_comments
                else:
                    raise ValueError(f"read_file: invalid state '{state}'")

                if strip_regex is not None:
                    line = strip_regex.sub('', line)
                if strip_leading:
                    line = re.sub(r'^\s*', '', line)
                if strip_trailing:
                    line = re.sub(r'\s*$', '', line)
                if line.endswith('\n'):
                    line = line[:-1]

                if not (strip_empty and re.match(r'^\s*$', line)):
                    code_buffer.append(line + '\n')

        if not include_pushed:
            file_names.pop(0)
            file_positions.pop(0)

    return ''.join(code_buffer)
