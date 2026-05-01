# Shared LaTeX encoding utilities.
# Andrew Benson (ported to Python 2026)
from __future__ import annotations

import re

__all__ = ['latex_encode']


_LATEX_SPECIAL = [
    ('\\', r'\textbackslash{}'),
    ('{',  r'\{'),
    ('}',  r'\}'),
    ('$',  r'\$'),
    ('&',  r'\&'),
    ('%',  r'\%'),
    ('#',  r'\#'),
    ('^',  r'\^{}'),
    ('_',  r'\_'),
    ('~',  r'\~{}'),
]

_LATEX_SPECIAL_MAP = dict(_LATEX_SPECIAL)
_LATEX_SPECIAL_RE  = re.compile(r'[\\{}$&%#^_~]')


def latex_encode(text: str) -> str:
    """Escape special LaTeX characters in text (equivalent to LaTeX::Encode)."""
    return _LATEX_SPECIAL_RE.sub(
        lambda m: _LATEX_SPECIAL_MAP[m.group(0)],
        text,
    )
