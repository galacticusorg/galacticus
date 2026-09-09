"""Tests for the source scanning of scripts/aux/linkChecker.py.

`source/` is a directory hierarchy, and essentially every URL embedded in the
source lives in a file below its top level. The sweep must therefore walk the
tree: an earlier non-recursive listing saw only the three files at the top of
`source/`, so the ~200 URLs in the rest of the tree went unchecked entirely.
The cases below pin that behavior, along with the file selection and the source
locations recorded for each URL (which the checker prints when reporting a
broken link, and which are useless if they do not name the real file).

Walking the tree also exposed how a URL is extracted from a line, which had
never been exercised against the bulk of the source. Three forms there produce
a "URL" that is not one, and each was reported as a broken link:

  * an RST link target ``<url>`` inside an XML directive writes its delimiters
    as ``&lt;``/``&gt;``, which the URL pattern cannot see, so unescaping left
    a ``>`` stuck to the end of the URL;
  * trailing punctuation and an enclosing ``)`` interleave (``…/meta.)``), and
    a single pass of each left one of them behind; and
  * a URL assembled by Fortran string concatenation is only a fragment, and no
    complete URL exists in the source to check.
"""

import importlib.util
import os
import sys

import pytest

_CHECKER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'linkChecker.py')


@pytest.fixture(scope='module')
def linkChecker():
    spec = importlib.util.spec_from_file_location('linkChecker', _CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def sourceTree(tmp_path):
    """A miniature `source/` tree, with URLs at several depths."""
    def write(relativePath, url):
        path = tmp_path / relativePath
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('  !![\n'
                        f'  <workaround type="gfortran" PR="1" url="{url}"/>\n'
                        '  !!]\n')
    write('top.F90',                              'https://example.com/top')
    write('utility/nested.F90',                   'https://example.com/nested')
    write('utility/IO/HDF5/deeplyNested.F90',     'https://example.com/deep')
    write('numerical/included.Inc',               'https://example.com/included')
    write('notSource.txt',                        'https://example.com/ignored')
    write('utility/notSource.md',                 'https://example.com/alsoIgnored')
    return tmp_path


def test_urls_are_found_at_every_depth(linkChecker, sourceTree):
    urls = {}
    linkChecker.scan_sources(str(sourceTree), urls)
    assert set(urls) == {
        'https://example.com/top',
        'https://example.com/nested',
        'https://example.com/deep',
        'https://example.com/included',
    }


def test_non_source_files_are_ignored(linkChecker, sourceTree):
    urls = {}
    linkChecker.scan_sources(str(sourceTree), urls)
    assert 'https://example.com/ignored'     not in urls
    assert 'https://example.com/alsoIgnored' not in urls


def test_the_recorded_location_names_the_real_file(linkChecker, sourceTree):
    urls = {}
    linkChecker.scan_sources(str(sourceTree), urls)
    source, = urls['https://example.com/deep']
    assert source['file']       == 'deeplyNested.F90'
    assert source['lineNumber'] == 2
    assert os.path.join(source['path'], source['file']) == \
        str(sourceTree / 'utility' / 'IO' / 'HDF5' / 'deeplyNested.F90')


def test_a_missing_source_directory_is_not_an_error(linkChecker, tmp_path):
    urls = {}
    linkChecker.scan_sources(str(tmp_path / 'doesNotExist'), urls)
    assert urls == {}


@pytest.mark.parametrize('line,expected', [
    # An RST link target inside an XML directive: delimiters entity-encoded.
    ('`web page &lt;https://example.net/a&gt;`_',        'https://example.net/a'),
    # Trailing sentence punctuation, and an enclosing parenthesis, interleaved.
    ('see https://example.net/a/meta.)',                 'https://example.net/a/meta'),
    ('see https://example.net/a.',                       'https://example.net/a'),
    ('(see https://example.net/a)',                      'https://example.net/a'),
    # A balanced parenthesis belongs to the URL and is kept.
    ('https://example.net/Foo_(bar)',                    'https://example.net/Foo_(bar)'),
    # An apostrophe inside a URL is kept; one wrapping it is not.
    ("url='https://example.net/a'",                      'https://example.net/a'),
    ("https://example.net/Claudia's_Model.html",         "https://example.net/Claudia's_Model.html"),
])
def test_url_extraction(linkChecker, tmp_path, line, expected):
    (tmp_path / 'test.F90').write_text(line + '\n')
    urls = {}
    linkChecker.scan_sources(str(tmp_path), urls)
    assert list(urls) == [expected]


@pytest.mark.parametrize('line', [
    # Fortran concatenation, the URL fragment ending at a double quote ...
    'url="https://example.net/archive/v"//version//".tar.gz"',
    # ... at a single quote, with the concatenation inside the matched text ...
    "call download('https://example.net/releases/c'//char(version)//'.tar.gz')",
    # ... and with the operator separated from the literal.
    'url="https://example.net/archive/v" // version',
])
def test_concatenated_url_fragments_are_skipped(linkChecker, tmp_path, line):
    (tmp_path / 'test.F90').write_text(line + '\n')
    urls = {}
    linkChecker.scan_sources(str(tmp_path), urls)
    assert urls == {}


@pytest.mark.parametrize('url,excluded', [
    # RFC 2606 reserved names, used in the source for URLs that are
    # deliberately not real (the download test's injection payloads).
    ('https://invalid.invalid/x',      True),
    ('http://foo.invalid',             True),
    ('https://example.org/data.txt',   True),
    ('https://www.example.com/a',      True),
    # A real host that merely contains the reserved name is not excluded.
    ('https://not-invalid.example.io/a', False),
    ('https://example.org.uk/a',        False),
    ('https://gcc.gnu.org/bugzilla/',   False),
])
def test_reserved_domains_are_excluded(linkChecker, url, excluded):
    assert linkChecker.is_excluded(url) is excluded


@pytest.mark.parametrize('url,tolerated', [
    ('https://www.openmp.org/specifications/',                     True),
    ('http://math.stackexchange.com/questions/40713/x',            True),
    ('https://goldbook.iupac.org/terms/view/B00598',               True),
    ('https://journals.aps.org/prl/abstract/10.1103/PhysRevLett',  True),
    ('https://gcc.gnu.org/bugzilla/',                              False),
])
def test_hosts_that_block_the_checker_tolerate_a_forbidden(linkChecker, url,
                                                           tolerated):
    """These hosts answer 403 to the checker (bot/datacenter-IP blocking by a
    CDN) while serving the page to ordinary clients, so a 403 from them is the
    expected result rather than a broken link."""
    assert linkChecker.tolerates_forbidden(url) is tolerated


def test_urls_are_collected_across_calls(linkChecker, sourceTree, tmp_path):
    """`main` passes one dictionary through the source, docs and wiki sweeps."""
    urls = {'https://example.com/preexisting': [{'file': 'a', 'path': 'b',
                                                 'lineNumber': 1}]}
    linkChecker.scan_sources(str(sourceTree), urls)
    assert 'https://example.com/preexisting' in urls
    assert 'https://example.com/deep'        in urls


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
