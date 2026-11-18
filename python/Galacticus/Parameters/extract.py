"""Reconstruct a Galacticus parameter file from the ``Parameters`` group of a Galacticus output file.

Andrew Benson (17-November-2026)

Galacticus records every parameter used by a model (including those which took default values) as attributes of a hierarchy
of groups beneath ``Parameters`` in its output file. This module reverses that, giving an XML parameter file which describes
the model which generated the output. It is used by ``scripts/parameters/parametersExtract.py``, and by the launcher's
``postprocess`` command, which starts from the parameters of the model whose output is being postprocessed.
"""

import copy
import re

import h5py
import lxml.etree as ET
import numpy as np

# Suffix used to mark, in the output file's `Parameters` group, a parameter which took a default value. This must match
# `defaultedMarkerSuffix` in `source/utility/input_parameters.F90`.
defaultedMarkerSuffix = "{defaulted}"

__all__ = ['parameters_from_hdf5', 'dereference_parameters']


def _value_string(value):
    """Convert an HDF5 attribute value to the string form used in a parameter file."""
    if isinstance(value, bytes):
        value = value.decode()
    elif isinstance(value, list):
        value = " ".join(map(lambda x: x.decode(), value))
    elif isinstance(value, np.ndarray):
        if isinstance(value[0], bytes):
            value = " ".join(map(lambda x: x.decode(), value))
        else:
            value = " ".join(map(lambda x: str(x), value))
    return str(value)


def parameters_from_hdf5(source, dereference=False, include_defaulted=True):
    """Return the parameters recorded in a Galacticus output file, as an ``lxml`` ``parameters`` element.

    ``source`` is the name of the file, or an open ``h5py.File``. If ``dereference`` is true, any ``idRef`` pointer is
    replaced by a copy of the parameter it refers to, and ``id`` attributes are then removed. A ``ValueError`` is raised if a
    reference can not be resolved uniquely.

    If ``include_defaulted`` is false, parameters which took their default values in the model (along with any
    sub-parameters) are omitted, so that a model run from the result takes the defaults afresh. This is needed where defaults
    depend on one another: a parameter whose default is computed from another parameter, and which may not be given together
    with it (``massEnvironment`` and ``radiusEnvironment`` of ``haloEnvironmentNormal``, for example), would otherwise be set
    explicitly, and the model would refuse to run.
    """
    if isinstance(source, h5py.File):
        return _extract(source, dereference, include_defaulted)
    with h5py.File(source, "r") as fileIn:
        return _extract(fileIn, dereference, include_defaulted)


def _extract(fileIn, dereference, include_defaulted):
    parameters = ET.Element("parameters")
    # Elements for parameters which took default values, to be removed if defaulted parameters are not to be included.
    defaulted = []

    def createStructure(name, node):
        # Create the element for each group, so that the structure of the parameter file exists before values are assigned.
        if not isinstance(node, h5py.Group):
            return
        # Split the parameter name into individual elements.
        elementNames = name.split("/")
        # Find the parent element for this new parameter.
        if len(elementNames) > 1:
            parent = parameters.find("/".join(elementNames[0:-1]))
        else:
            parent = parameters
        # Array parameters (i.e. multiple copies of the same named parameter) are defined by a `[N]` suffix. Detect these.
        isArray = re.match(r'(.*)\[(\d+)\]', elementNames[-1])
        if isArray:
            # For array parameters, construct a list of elements of length equal to the current array index.
            elementName = isArray.group(1)
            elementIndex = int(isArray.group(2))
            members = parent.findall(elementName)
            if len(members) <= elementIndex:
                for i in range(len(members), elementIndex):
                    ET.SubElement(parent, elementName)
        else:
            # For non-array parameters simply add an appropriately-named child element.
            ET.SubElement(parent, elementNames[-1])

    def assignValues(name, node):
        # Assign the values of parameters held as attributes of this group.
        if not isinstance(node, h5py.Group):
            return
        # Find the parent element for this group.
        if name == "parameters":
            parent = parameters
        else:
            parent = parameters.find(name)
        # Iterate over all attributes in this group.
        for attributeName in node.attrs.keys():
            # Skip markers recording that a parameter took a default value. These are not parameters themselves - they
            # annotate the parameter of the same name, minus the suffix - and the braces are not valid in an XML element
            # name, so they must not be emitted here.
            if attributeName.endswith(defaultedMarkerSuffix):
                continue
            isDefaulted = (attributeName + defaultedMarkerSuffix) in node.attrs
            value = _value_string(node.attrs[attributeName])
            isReference = re.match(r'\{idRef:([a-zA-Z0-9_]+)\}', value)
            idReference = isReference.group(1) if isReference else None
            isTarget = re.match(r'(.*)\{id:([a-zA-Z0-9_]+)\}', attributeName)
            if isTarget:
                attributeName = isTarget.group(1)
                idTarget = isTarget.group(2)
            # A prefix "meta:" indicates a non-parameter that was added as metadata. Include it but add an `ignoreWarnings`
            # attribute to avoid triggering warnings.
            isMeta = re.match(r'^meta:(.+)', attributeName)
            if isMeta:
                attributeName = isMeta.group(1)
            # Set the value of the parameter (creating the element first if needed).
            if attributeName in node:
                child = parent.find(attributeName)
            else:
                child = ET.SubElement(parent, attributeName)
            if isReference:
                child.set('idRef', idReference)
            else:
                child.set('value', value)
                if isTarget:
                    child.set('id', idTarget)
            if isMeta:
                child.set('ignoreWarnings', 'true')
            if isDefaulted:
                defaulted.append(child)

    parametersGroup = fileIn['Parameters']
    # Add a `formatVersion` element.
    formatVersion = ET.SubElement(parameters, 'formatVersion')
    formatVersion.text = "2"
    # Add a `lastModified` element.
    version = fileIn['Version']
    if 'gitHash' in version.attrs.keys():
        lastModified = ET.SubElement(parameters, 'lastModified')
        lastModified.set('revision', version.attrs['gitHash'].decode())
    # Construct the structure of the parameter file by visiting all groups in `Parameters`.
    parametersGroup.visititems(createStructure)
    # Assign top-level parameter values.
    assignValues('parameters', parametersGroup)
    # Assign parameter values in all sub-parameters.
    parametersGroup.visititems(assignValues)
    # Remove defaulted parameters if requested.
    if not include_defaulted:
        for element in defaulted:
            parent = element.getparent()
            if parent is not None:
                parent.remove(element)
    # Dereference `idRef`s if requested.
    if dereference:
        dereference_parameters(parameters)
    return parameters


def dereference_parameters(parameters, remove_ids=True):
    """Replace, in place, each ``idRef`` pointer in ``parameters`` (an lxml element) by a copy of the parameter it refers to.

    A ``ValueError`` is raised if a reference can not be resolved uniquely. Copies do not carry the ``id`` attributes of the
    parameters copied, so that each ``id`` remains unique. References within a copied parameter are themselves resolved. If
    ``remove_ids`` is true, all ``id`` attributes are then removed; otherwise they are kept, so that references added later (in
    parameters appended to ``parameters``, for example) can still be resolved by a further call.
    """
    # Each pass resolves the references present at its start. A copied parameter may contain references of its own, so passes
    # continue until none remain. A chain of references can be no longer than the number of parameters which are targets, so
    # more passes than that indicate a cycle.
    passesMaximum = len(parameters.findall('.//*[@id]')) + 1
    for _ in range(passesMaximum):
        references = parameters.findall('.//*[@idRef]')
        if not references:
            break
        for reference in references:
            name = reference.tag
            id = reference.attrib['idRef']
            target = parameters.findall(".//" + name + "[@id='" + id + "']")
            if len(target) == 0:
                raise ValueError('`<' + name + ' id="' + id + '"/>` not found')
            elif len(target) > 1:
                raise ValueError('multiple `<' + name + ' id="' + id + '"/>` found')
            replacement = copy.deepcopy(target[0])
            for element in replacement.iter():
                element.attrib.pop('id', None)
            reference.getparent().replace(reference, replacement)
    else:
        raise ValueError('references between parameters form a cycle')
    if remove_ids:
        for target in parameters.findall('.//*[@id]'):
            del target.attrib['id']
    return parameters
