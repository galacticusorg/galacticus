#!/usr/bin/env python3
# Shared spell-checking logic for LaTeX and plain-text fragments.
# Used by spellChecker.py and embeddedAnalyzer.py.
# Andrew Benson (28-February-2023) [Python port]

import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

# Module-level cache so the word list is only loaded once per process.
_spell_words_cache = None


def _load_spell_words():
    """Load custom spell-check words from stateStorables.xml and aux/words.dict."""
    global _spell_words_cache
    if _spell_words_cache is not None:
        return _spell_words_cache

    words = []

    # Load function-class names and instance names from stateStorables.xml.
    storables_path = "work/build/stateStorables.xml"
    if os.path.exists(storables_path):
        try:
            tree = ET.parse(storables_path)
            root = tree.getroot()

            # XML::Simple serialises the functionClasses array-of-hashes as repeated
            # <functionClasses name="..."> elements keyed by their name attribute.
            class_names = [el.get('name') for el in root.findall('functionClasses')
                           if el.get('name')]

            # functionClassInstances are repeated <functionClassInstances> elements
            # whose text content is the instance name.
            instances = [el.text.strip() for el in root.findall('functionClassInstances')
                         if el.text and el.text.strip()]

            # Add raw class names and the prefix-without-"Class" suffix.
            for cn in class_names:
                words.append(cn)
                if cn.endswith('Class'):
                    words.append(cn[:-5])

            # Add instance names.
            words.extend(instances)

            # For each instance, find the longest matching class prefix and add the
            # lowercase-first suffix (e.g. "cosmologyFunctionsSimple" → "simple").
            for inst in instances:
                best_len = 0
                best_suffix = None
                for cn in class_names:
                    prefix = cn[:-5] if cn.endswith('Class') else cn
                    if inst.startswith(prefix) and len(prefix) > best_len:
                        best_len = len(prefix)
                        suffix = inst[len(prefix):]
                        if suffix:
                            best_suffix = suffix[0].lower() + suffix[1:]
                if best_suffix:
                    words.append(best_suffix)
        except ET.ParseError:
            pass

    # Load the project custom word list.
    if os.path.exists("aux/words.dict"):
        with open("aux/words.dict", errors='replace') as f:
            for line in f:
                word = line.rstrip('\n')
                if word:
                    words.append(word)

    _spell_words_cache = words
    return words


def _extract_balanced(text, pos):
    """Return (matched_string, end_pos) for the brace-balanced substring starting
    at pos (which must point at '{').  Returns (None, pos) on failure."""
    if pos >= len(text) or text[pos] != '{':
        return None, pos
    depth = 0
    for i in range(pos, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[pos:i + 1], i + 1
    return None, pos


def _preprocess_line(line, is_latex, spell_words):
    """Return a spell-check-friendly version of one line."""

    # Split camelCase words into space-separated components, unless the word is
    # in the known-words list or is the special token "FoX".
    def _split_camel(m):
        word = m.group(0)
        if word in spell_words:
            return word
        if word == 'FoX':
            return 'fox'
        return re.sub(r'([a-z0-9])([A-Z0-9])', r'\1 \2', word)

    # Match camelCase words not preceded by backslash.
    line = re.sub(
        r'(?<!\\)\b([a-zA-Z0-9][a-zA-Z]*(?:[a-z0-9][a-zA-Z]*[A-Z]|[A-Z][a-zA-Z]*[a-z])[a-zA-Z0-9]*)\b',
        _split_camel,
        line
    )

    if is_latex:
        # ── Subscripts starting with "_{" ────────────────────────────────────
        new_line = ''
        rest = line
        while True:
            m = re.search(r'_\{', rest)
            if not m:
                break
            new_line += rest[:m.start() + 1]   # keep up to and including '_'
            braced, end = _extract_balanced(rest, m.start() + 1)
            if braced is None:
                new_line += rest[m.start() + 1:]
                rest = ''
                break
            inner = braced[1:-1]
            # Keep only \mathrm{…} sub-contents that are longer than 2 chars.
            kept = ''
            sub = inner
            while True:
                m2 = re.search(r'\\mathrm(\{)', sub)
                if not m2:
                    break
                mc, e2 = _extract_balanced(sub, m2.start(1))
                if mc is None:
                    break
                content = mc[1:-1]
                kept += '\\mathrm{' + (content if len(content) > 2 else '') + '}'
                sub = sub[m2.start(1) + len(mc):]
            new_line += '{' + kept + '}'
            rest = rest[end:]
        line = new_line + rest

        # ── Subscripts of the form "_\mathrm{…}" ────────────────────────────
        new_line = ''
        rest = line
        while True:
            m = re.search(r'_\\mathrm\{', rest)
            if not m:
                break
            tag_end = m.start() + len(r'_\mathrm')
            new_line += rest[:tag_end]
            braced, end = _extract_balanced(rest, tag_end)
            if braced is None:
                new_line += rest[tag_end:]
                rest = ''
                break
            inner = braced[1:-1]
            new_line += '{' + (inner if len(inner) > 2 else '') + '}'
            rest = rest[end:]
        line = new_line + rest

        # ── Glossary macros ──────────────────────────────────────────────────
        # Replace the argument of \gls, \glslink, etc. with an empty group so
        # hunspell does not try to spell-check the label.
        for macro in (r'\\gls', r'\\glslink'):
            new_line = ''
            rest = line
            while True:
                m = re.search(macro + r'\{', rest)
                if not m:
                    break
                new_line += rest[:m.start() + len(macro.replace('\\\\', '\\'))]
                braced, end = _extract_balanced(rest, m.end() - 1)
                if braced is None:
                    new_line += rest[m.end() - 1:]
                    rest = ''
                    break
                new_line += '{}'
                rest = rest[end:]
            line = new_line + rest

        # \newacronym{}{} — drop both brace groups
        new_line = ''
        rest = line
        while True:
            m = re.search(r'\\newacronym\{', rest)
            if not m:
                break
            new_line += rest[:m.start() + len(r'\newacronym')]
            pos = m.end() - 1
            b1, pos = _extract_balanced(rest, pos)
            b2, pos = _extract_balanced(rest, pos) if b1 else (None, pos)
            if b1 is None:
                new_line += rest[m.end() - 1:]
                rest = ''
                break
            new_line += '{}{}'
            rest = rest[pos:]
        line = new_line + rest

        # \newglossaryentry{}{name=…} — drop first group, keep rest
        line = re.sub(r'\\newglossaryentry\{[^}]*\}\{name=\{[^}]*\}', r'\\newglossaryentry{}{}', line)
        line = line.replace('firstplural=', 'first plural=')

        # Accent translations (we use the unaccented character so hunspell
        # personal dictionaries work reliably).
        line = re.sub(r"\\'e", 'e', line)
        line = re.sub(r'\\"o',  'o', line)

        # Remove \href{URL}{…} — drop the URL argument entirely.
        line = re.sub(r'\\href\{[^}]+\}', '', line)
    else:
        # Plain text: remove bare URLs.
        line = re.sub(r'https?://\S+(\s|\))', r'\1', line)

    return line


def spell_check(text, text_type, file_name_original):
    """Spell-check a text fragment (string).  text_type is 'latex' or 'text'."""
    is_latex = text_type == 'latex'
    suffix = '.tex' if is_latex else '.txt'
    fd, tmp_name = tempfile.mkstemp(suffix=suffix)
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write(text)
        return spell_check_file(tmp_name, file_name_original)
    finally:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass


def spell_check_file(file_name, file_name_original):
    """Spell-check a file on disk.  Returns a warnings string (may be empty)."""
    is_latex  = file_name.endswith('.tex')
    words     = _load_spell_words()
    dic_file  = file_name + '.dic'
    aff_file  = file_name + '.aff'
    spell_file = file_name + '.spell'

    # Write the personal dictionary.
    unique_words = sorted(set(w for w in words if w))
    with open(dic_file, 'w') as fh:
        fh.write(f'{len(unique_words)}\n')
        fh.write('\n'.join(unique_words) + '\n')

    # Copy the affix file (controls morphological rules).
    if os.path.exists('aux/words.aff'):
        shutil.copy('aux/words.aff', aff_file)
    else:
        open(aff_file, 'w').close()

    # Pre-process the file for spell checking.
    with open(file_name, 'r', errors='replace') as fi, \
         open(spell_file, 'w') as fo:
        for line in fi:
            fo.write(_preprocess_line(line, is_latex, unique_words))

    # Run hunspell in list mode (-l).
    cmd = ['hunspell', '-l']
    if is_latex:
        cmd.append('-t')
    cmd += ['-i', 'utf-8', '-d', f'en_US,{file_name}', spell_file]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        raw_output = result.stdout
    except (OSError, subprocess.SubprocessError):
        raw_output = ''

    # Collect unique misspelled words (case-folded).
    word_counts = {}
    for w in raw_output.splitlines():
        w = w.strip().lower()
        if w:
            word_counts[w] = word_counts.get(w, 0) + 1

    for path in (spell_file, dic_file, aff_file):
        try:
            os.unlink(path)
        except OSError:
            pass

    warnings = ''
    for word in sorted(word_counts):
        count = word_counts[word]
        count_str = f' ({count} instances)' if count > 1 else ''
        warnings += (
            f":warning: Possible misspelled word '{word}'{count_str}"
            f" in file '{file_name_original}'\n"
        )
    return warnings
