#!/usr/bin/env python3
"""Audit the documentation embedded in Galacticus source against the code.

Galacticus carries its reference documentation inside its sources: ``!!{ ... !!}``
blocks hold reStructuredText, and the ``!![ ... !!]`` directive blocks carry
``<description>`` elements which are rendered into the manuals. Both refer to
parameters, classes, implementations, and bibliography keys by name, and nothing
in the build checks that those names still exist---so a parameter rename leaves
the prose which mentions it silently stale, with no warning from the
documentation build.

The following checks are made, over both the embedded documentation and the
hand-written manuals under `docs/manuals/`:

  * every ``[parameter]`` reference names a parameter, class, or implementation
    which exists;
  * every ``:cite:`` key has an entry in `docs/Galacticus.bib`;
  * every ``:galacticus-class:`` reference names a type which exists;
  * every implementation named in an XML example exists;
  * no interpreted-text role runs on into the word which follows it.

The vocabulary of valid parameter, class, and implementation names is taken from
the parameter catalog (see `scripts/build/parameterCatalog.py`), supplemented by
the component classes---which the catalog does not cover---and by the derived
types declared in the source, which `:galacticus-class:` is also applied to.

Usage:
    docBlockAudit.py [--check]

With `--check` the script exits with a non-zero status if any problem is found,
so that it may be used as a lint.

Andrew Benson (2026).
"""

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

# Make the in-tree Python packages importable when run directly.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, 'python'))

from Galacticus.Parameters.catalog import build_catalog

# Directories excluded from the scan: vendored sources.
EXCLUDED = ('source/external/',)

# A documentation block is `!!{ ... !!}`; a directive block is `!![ ... !!]`.
DOC_BLOCK = re.compile(r'!!\{(.*?)!!\}', re.DOTALL)
DIRECTIVE_BLOCK = re.compile(r'^\s*!!\[\s*$(.*?)^\s*!!\]\s*$', re.MULTILINE | re.DOTALL)
DOC_FORMAT = re.compile(r'^[ \t]*(?:RST|latex)[ \t]*\n')

# Markup scanned for in documentation text.
PARAMETER_REFERENCE = re.compile(r'``\[([^\]`]+)\]``')
CITE_REFERENCE = re.compile(r':cite:[a-z]*:`([^`]+)`')
CLASS_REFERENCE = re.compile(r':galacticus-class:`([^`]+)`')
XML_EXAMPLE = re.compile(r'<(\w+)\s+value="([^"]*)"')
# A role must be followed by whitespace or punctuation: a letter immediately
# after the closing backtick (":math:`5500\,`Å") leaves the role unrecognized,
# and reStructuredText reports it as an error.
ROLE_RUN_ON = re.compile(r':[a-z][a-z:-]*:`[^`]*`(?=[^\W\d_])')
INLINE_LITERAL = re.compile(r'``.*?``', re.DOTALL)

# Derived types, which `:galacticus-class:` is also applied to.
DERIVED_TYPE = re.compile(r'^\s*type\s*(?:,[^:]*)?::\s*(\w+)', re.MULTILINE | re.IGNORECASE)

# Not every parameter is declared by an `<inputParameter>` directive: those read
# directly from the parameter object are found here instead.
DIRECT_PARAMETER = re.compile(r"%(?:value|copiesCount|subParameters|copyInstance)\s*\(\s*'([^']+)'")

# A `[reference]` containing any of these is prose or a placeholder---
# ``[L_1,...,L_N]``, or ``[%5.2f|darkMatterParticle/mass]``, for example---and
# not the name of a parameter.
NOT_A_NAME = re.compile(r'[\s%|,:=<>()\[\]{}…]|^[-+.\d]')

# A single character is never a parameter name: `` [N] `` in the documentation
# of XML handling is a positional XPath predicate, not a reference.
MINIMUM_NAME_LENGTH = 2

# Separators used to write a sub-parameter of a class, as in
# ``[componentBlackHole--mass]``.
SUB_PARAMETER = re.compile(r'/|--')

BIBLIOGRAPHY_ENTRY = re.compile(r'@\w+\s*\{\s*([^,\s]+)\s*,')


def repositoryRoot():
    """Return the root of the Galacticus source tree."""
    return os.environ.get('GALACTICUS_EXEC_PATH', os.getcwd())


def sourceFiles(root):
    """Iterate over the paths, relative to `root`, of all non-vendored Fortran sources."""
    for directoryPath, directoryNames, fileNames in os.walk(os.path.join(root, 'source')):
        for fileName in sorted(fileNames):
            if not fileName.endswith(('.F90', '.Inc')):
                continue
            path = os.path.relpath(os.path.join(directoryPath, fileName), root)
            if not any(path.startswith(excluded) for excluded in EXCLUDED):
                yield path


def manualFiles(root):
    """Iterate over the paths, relative to `root`, of all hand-written manual sources."""
    for directoryPath, directoryNames, fileNames in os.walk(os.path.join(root, 'docs', 'manuals')):
        for fileName in sorted(fileNames):
            if fileName.endswith('.rst'):
                yield os.path.relpath(os.path.join(directoryPath, fileName), root)


def readFile(root, path):
    """Return the content of a file."""
    with open(os.path.join(root, path), encoding='utf-8', errors='replace') as file:
        return file.read()


def parseDirectives(block):
    """Parse a directive block into a list of elements, or return an empty list if it is not well formed.

    Directive blocks contain a sequence of XML fragments rather than a single
    document, so a root element is supplied. A block which fails to parse is
    skipped: validating the directives themselves is the job of the build, not
    of this audit.
    """
    try:
        return list(ET.fromstring('<directives>' + block + '</directives>'))
    except ET.ParseError:
        return []


def docText(content):
    """Return the documentation text of a source file: its doc blocks, plus the descriptions in its directives."""
    texts = [DOC_FORMAT.sub('', block, count=1) for block in DOC_BLOCK.findall(content)]
    for block in DIRECTIVE_BLOCK.findall(content):
        for element in parseDirectives(block):
            for description in element.iter('description'):
                texts.append(''.join(description.itertext()))
    return texts


class Vocabulary:
    """The names which documentation may legitimately refer to."""

    def __init__(self, root, paths):
        catalog = build_catalog(os.path.join(root, 'source'))
        self.classes = set(catalog['functionClasses'].keys())
        self.implementations = set(catalog['implementations'].keys())
        # Implementations are named in parameter files by their label---the
        # implementation name with the class name prefix removed. Labels are
        # stored lower-cased, since the case of the leading character is not
        # significant when a parameter file is read.
        self.labelsByClass = {
            name: {label.lower() for label in functionClass['implementations']}
            for name, functionClass in catalog['functionClasses'].items()
        }
        self.parameters = set()
        for implementation in catalog['implementations'].values():
            self.parameters.update(parameter['name'] for parameter in implementation['parameters'])
            self.parameters.update(object_['parameterName'] for object_ in implementation['objects'])
        self.derivedTypes = set()
        self.parametersByFile = {}
        for path in paths:
            self._addSource(root, path)

    def _addSource(self, root, path):
        """Add the names declared by a single source file."""
        content = readFile(root, path)
        self.derivedTypes.update(DERIVED_TYPE.findall(content))
        fileParameters = set(DIRECT_PARAMETER.findall(content))
        for block in DIRECTIVE_BLOCK.findall(content):
            for element in parseDirectives(block):
                for descendant in element.iter():
                    if descendant.tag == 'inputParameter':
                        name = descendant.find('name')
                        if name is not None and name.text is not None:
                            fileParameters.add(name.text.strip())
                    elif descendant.tag == 'objectBuilder':
                        name = descendant.get('parameterName', descendant.get('class'))
                        if name is not None:
                            fileParameters.add(name.strip())
                    elif descendant.tag == 'inputParametersValidate':
                        # Parameters which are read directly, and so are
                        # declared here as permissible rather than by an
                        # `<inputParameter>` directive.
                        for attribute in ('multiParameters', 'extraAllowedNames'):
                            for name in re.split(r'[,\s]+', descendant.get(attribute, '')):
                                if name:
                                    fileParameters.add(name)
                # Components are not covered by the parameter catalog, so
                # collect their classes, implementations, and properties here.
                if element.tag == 'component':
                    self._addComponent(element, fileParameters)
        self.parametersByFile[path] = fileParameters
        self.parameters |= fileParameters

    def _addComponent(self, element, fileParameters):
        """Add the names declared by a `<component>` directive."""
        className, implementationName = element.find('class'), element.find('name')
        if className is None or className.text is None:
            return
        className = className.text.strip()
        className = 'component' + className[0].upper() + className[1:]
        self.classes.add(className)
        self.parameters.add(className)
        # The generated types for the component class and this implementation.
        self.derivedTypes.add('node' + className[0].upper() + className[1:])
        if implementationName is not None and implementationName.text is not None:
            implementationName = implementationName.text.strip()
            self.labelsByClass.setdefault(className, set()).add(implementationName.lower())
            self.derivedTypes.add('node' + className[0].upper() + className[1:]
                                  + implementationName[0].upper() + implementationName[1:])
        for property_ in element.iter('property'):
            name = property_.find('name')
            if name is not None and name.text is not None:
                fileParameters.add(name.text.strip())
        # Every component class also accepts the `null` implementation, which is
        # generated rather than declared in the source.
        self.labelsByClass.setdefault(className, set()).add('null')

    def isParameter(self, name, fileParameters):
        """Return true if `name` is a parameter, class, or implementation which exists."""
        parts = [part for part in SUB_PARAMETER.split(name) if part]
        return all(
            part in fileParameters or part in self.parameters
            or part in self.classes or part in self.implementations
            for part in parts
        )

    def isType(self, name):
        """Return true if `name` is a type which exists.

        The `:galacticus-class:` role resolves the abstract base type `xClass`
        to the page for class `x` (see `docs/conf.py`), and renders anything for
        which no page exists---a utility type, for example---as inline code,
        which is harmless. A reference is therefore stale only if it names no
        type at all.
        """
        candidates = [name] + ([name[:-len('Class')]] if name.endswith('Class') else [])
        return any(candidate in self.classes or candidate in self.implementations
                   or candidate in self.derivedTypes for candidate in candidates)

    def isImplementation(self, className, label):
        """Return true if `label` names an implementation of class `className`."""
        return label.lower() in self.labelsByClass.get(className, set())


def checkText(path, text, fileParameters, vocabulary, keys, problems):
    """Check a single piece of documentation text, appending any problems found."""
    for match in PARAMETER_REFERENCE.finditer(text):
        name = match.group(1)
        if len(name) >= MINIMUM_NAME_LENGTH and not NOT_A_NAME.search(name) \
           and not vocabulary.isParameter(name, fileParameters):
            problems.append((path, 'unknown parameter reference: ' + name))
    for match in XML_EXAMPLE.finditer(text):
        className, label = match.group(1), match.group(2)
        if className in vocabulary.classes and not vocabulary.isImplementation(className, label):
            problems.append((path, 'XML example names unknown implementation: ' + className + '=' + label))
    for match in ROLE_RUN_ON.finditer(text):
        problems.append((path, 'role runs on into the following word (insert an escaped space): '
                         + match.group(0)))
    # Roles are not interpreted inside an inline literal, so text which merely
    # shows the syntax of a role---as the coding guide does---is not a reference.
    text = INLINE_LITERAL.sub('', text)
    for match in CITE_REFERENCE.finditer(text):
        for key in match.group(1).split(','):
            key = key.strip()
            if key and key not in keys:
                problems.append((path, 'unknown citation key: ' + key))
    for match in CLASS_REFERENCE.finditer(text):
        name = match.group(1)
        if not vocabulary.isType(name):
            problems.append((path, 'unknown class reference: ' + name))


def citationKeys(root):
    """Return the set of keys defined in the bibliography."""
    return set(BIBLIOGRAPHY_ENTRY.findall(readFile(root, os.path.join('docs', 'Galacticus.bib'))))


def main():
    parser = argparse.ArgumentParser(description='Audit embedded documentation blocks against the code.')
    parser.add_argument('--check', action='store_true',
                        help='exit with non-zero status if any problem is found')
    arguments = parser.parse_args()

    root = repositoryRoot()
    paths = list(sourceFiles(root))
    vocabulary = Vocabulary(root, paths)
    keys = citationKeys(root)

    problems = []
    for path in paths:
        for text in docText(readFile(root, path)):
            checkText(path, text, vocabulary.parametersByFile[path], vocabulary, keys, problems)
    for path in manualFiles(root):
        checkText(path, readFile(root, path), set(), vocabulary, keys, problems)

    for path, message in sorted(problems):
        print(path + ': ' + message)
    print('Found ' + str(len(problems)) + ' problems in ' + str(len(paths))
          + ' source files and the manuals.')

    sys.exit(1 if arguments.check and problems else 0)


if __name__ == '__main__':
    main()
