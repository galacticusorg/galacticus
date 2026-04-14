# Utility functions for reading and parsing Fortran source code.
# Andrew Benson (ported to Python 2026)

def get_fortran_line(file_obj):
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
            import re
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
        import re
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


def _fast_comment_check(line):
    """Return True if the line may need the full comment-position scan."""
    import re
    return bool(re.match(r"^([^'\"{\n]+?)![^\$]", line))


def _find_comment_position(line):
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
