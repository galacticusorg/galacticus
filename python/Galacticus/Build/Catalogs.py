"""Shared, cached loading of the XML build catalogs in `$BUILDPATH`.

The source-tree hooks consult the same few catalogs (`directiveLocations.xml`,
`stateStorables.xml`, `deepCopyActions.xml`, ...).  Each hook used to parse
them itself, so preprocessing one file parsed `directiveLocations.xml` (1.5 MB)
four times and `stateStorables.xml` seven, and a functionClass parent, which
processes every implementation, around forty times.  `load()` parses each
catalog once per process.

The cache is validated on every call by comparing the file's content with the
content it was parsed from (reading and comparing even the 1.5 MB
`directiveLocations.xml` costs ~0.1 ms, against ~13 ms to parse it), so a
catalog rewritten during the life of the process is always re-read.  A `stat()`
signature is not enough: a same-size rewrite in place within one file-system
timestamp tick leaves it unchanged.  The returned structure is *shared* between
every caller and must be treated as read-only.

Andrew Benson (2026)
"""
from __future__ import annotations

import os
import xml.etree.ElementTree as ET

from XML.Utils import xml_to_dict

__all__ = ['load']

# `(path, keyed_tags) -> (raw file content, converted catalog)`.
_cache: dict = {}


def load(path: str, keyed_tags: set[str] | None = None):
    """Return `xml_to_dict()` of the XML file at `path`, parsing it at most once
    per process while the file is unchanged.

    Raises `FileNotFoundError` if `path` does not exist, and `ET.ParseError` if
    it is not well-formed, exactly as parsing it directly would.
    """
    with open(path, 'rb') as fh:
        content = fh.read()
    key    = (os.path.abspath(path), frozenset(keyed_tags or ()))
    cached = _cache.get(key)
    if cached is not None and cached[0] == content:
        return cached[1]
    catalog     = xml_to_dict(ET.fromstring(content), keyed_tags=keyed_tags)
    _cache[key] = (content, catalog)
    return catalog
