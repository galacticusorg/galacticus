"""Provides extra list/dict utility functions.

Andrew Benson (ported to Python 2026)
"""
from __future__ import annotations

from typing import Any

__all__ = ['smart_push', 'as_array', 'hash_list', 'sorted_keys']


def smart_push(array: list, item: Any) -> None:
    """Intelligently append item(s) onto array.

    - If item is None, does nothing.
    - If item is a list, extends array with it.
    - Otherwise appends the single item.
    """
    if item is None:
        return
    if isinstance(item, list):
        array.extend(item)
    else:
        array.append(item)


def as_array(item: Any) -> list:
    """Return item coerced to a list."""
    result: list = []
    smart_push(result, item)
    return result


def hash_list(d: dict | None, key_as: str | None = None) -> list:
    """Return a list of dict values, sorted by key.

    If key_as is given, stamps the key into each value dict under that field name
    (modifying the values in-place).
    """
    if not d:
        return []
    result: list = []
    for key in sorted(d.keys()):
        val = d[key]
        if key_as is not None:
            if not isinstance(val, dict):
                val = {}
            val[key_as] = key
        result.append(val)
    return result


def sorted_keys(d: dict | None) -> list:
    """Return sorted list of dict keys.

    Returns [] if d is None or empty.
    """
    if not d:
        return []
    return sorted(d.keys())
