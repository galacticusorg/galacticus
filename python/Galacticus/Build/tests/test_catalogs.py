"""Tests for `Galacticus.Build.Catalogs`, the shared cached loader of the XML
build catalogs.

The cache must return the parsed catalog without re-parsing while the file is
unchanged, but must never serve stale content once the file is rewritten — the
test suite, and any long-lived process, can see a catalog change underneath it.
"""

import os
import xml.etree.ElementTree as ET

import pytest

from Galacticus.Build import Catalogs


def _write(path, text):
    path.write_text(text)


def test_load_matches_xml_to_dict_and_is_cached(tmp_path):
    path = tmp_path / 'stateStorables.xml'
    _write(path, '<storables><functionClasses><name>aClass</name></functionClasses></storables>')
    first = Catalogs.load(str(path))
    assert first == {'functionClasses': {'name': 'aClass'}}
    # A second call returns the very same object, without re-parsing.
    assert Catalogs.load(str(path)) is first


def test_load_rereads_a_rewritten_file(tmp_path):
    # A same-size rewrite in place, with the original timestamp restored: the
    # file's `stat()` signature is unchanged, so only its content tells.
    path = tmp_path / 'directiveLocations.xml'
    _write(path, '<directives><a><file>one.F90</file></a></directives>')
    status = os.stat(path)
    assert Catalogs.load(str(path)) == {'a': {'file': 'one.F90'}}
    _write(path, '<directives><a><file>two.F90</file></a></directives>')
    os.utime(path, ns=(status.st_atime_ns, status.st_mtime_ns))
    assert Catalogs.load(str(path)) == {'a': {'file': 'two.F90'}}


def test_load_caches_conversion_options_separately(tmp_path):
    path = tmp_path / 'openMPCriticalSections.xml'
    _write(path, '<criticalSections><critical name="x"><id>1</id></critical></criticalSections>')
    keyed = Catalogs.load(str(path), keyed_tags={'critical'})
    plain = Catalogs.load(str(path))
    assert keyed == {'critical': {'x': {'name': 'x', 'id': '1'}}}
    assert plain == {'critical': {'name': 'x', 'id': '1'}}


def test_load_errors_match_direct_parsing(tmp_path):
    with pytest.raises(FileNotFoundError):
        Catalogs.load(str(tmp_path / 'missing.xml'))
    path = tmp_path / 'broken.xml'
    _write(path, '<unclosed>')
    with pytest.raises(ET.ParseError):
        Catalogs.load(str(path))
