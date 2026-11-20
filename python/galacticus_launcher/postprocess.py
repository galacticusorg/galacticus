"""Postprocess the galaxies of a Galacticus output file with new property extractors.

Andrew Benson (17-November-2026)

A Galacticus output file records the parameters of the model which produced it. This module starts from those parameters,
replaces the task with ``postprocessOutputs`` and the property extractors with those requested, and runs the executable. The
galaxies are rebuilt from the output file (star formation histories, component masses, abundances, sizes, and inclinations), so
that new SEDs, emission lines, or dust-attenuated luminosities can be computed for them without evolving the model again. Since
the star formation history class, stellar populations, and node components come from the original model, they are configured
exactly as they were when the file was written.

The pieces are usable separately:

* :func:`build_parameters` returns the parameter file (as an ``lxml`` element) for postprocessing;
* :func:`write_selection` writes the HDF5 file listing which galaxies to postprocess;
* :func:`postprocess` does both, runs Galacticus, and returns the name of the output file.

The same is available from the command line as ``galacticus postprocess``.

The output file has the same structure as a Galacticus output: for each output group postprocessed,
``[outputsGroupName]/Output[i]/nodeData`` holds one row per galaxy postprocessed. The dataset
``nodeIndexPostprocessSource`` gives, for each row, the index (counting from zero) of the galaxy in the ``nodeData`` datasets
of the original file. Galaxies are *not* written in the order in which they were selected, so use this dataset to match rows
back to the original file.
"""

import copy
import os
import subprocess
import tempfile
from pathlib import Path

import h5py
import lxml.etree as ET
import numpy as np

__all__ = ["build_parameters", "write_selection", "postprocess"]


def _element(value, tag):
    """Return `value` as an lxml element, checking its tag.

    `value` may be an lxml element, a string of XML, or the name of a file containing XML.
    """
    if isinstance(value, ET._Element):
        # Copy the element, so that neither placing it in the parameters nor resolving its references alters the caller's tree.
        element = copy.deepcopy(value)
    elif isinstance(value, (str, os.PathLike)) and not str(value).lstrip().startswith("<"):
        element = ET.parse(str(value)).getroot()
    else:
        element = ET.fromstring(value)
    if element.tag != tag:
        raise ValueError(f"expected a <{tag}> element, but found <{element.tag}>")
    return element


def _replace(parameters, element):
    """Replace all top-level parameters with the tag of `element` by `element`."""
    for existing in parameters.findall(element.tag):
        parameters.remove(existing)
    parameters.append(element)


def _parameter(parent, name, value):
    """Append a `<name value="value"/>` parameter to `parent`."""
    child = ET.SubElement(parent, name)
    child.set("value", str(value))
    return child


def build_parameters(catalog, extractors, output, *, outputs_group="Outputs", selection_file=None, count_batch=1024,
                     count_nodes_buffered=1000, galactic_filter="always"):
    """Return the parameters (an lxml ``parameters`` element) to postprocess the galaxies of a Galacticus output file.

    Parameters
    ----------
    catalog
        The Galacticus output file whose galaxies are to be postprocessed. Its ``Parameters`` group is the starting point.
    extractors
        The property extractors to apply: a ``<nodePropertyExtractor>`` element (an lxml element, a string of XML, or the
        name of a file containing it), or a list of such elements. They are combined, together with the extractor of
        ``nodeIndexPostprocessSource``, in a ``multi`` extractor. They may refer (with ``idRef``) to parameters of the model
        (those given an ``id`` in its parameter file) or to each other; each reference is replaced by a copy of the parameter
        to which it refers.
    output
        The name of the file to which results are to be written.
    outputs_group
        The group in the output file containing the ``Output[i]`` groups to postprocess (``Lightcone`` for a lightcone, for
        example). The same name is used in the file written.
    selection_file
        An HDF5 file listing the galaxies to postprocess (see :func:`write_selection`), or ``None`` to postprocess all galaxies.
    count_batch
        The number of galaxies read from file at a time. Batches are processed in parallel.
    count_nodes_buffered
        The number of galaxies whose properties are accumulated before being written to file.
    galactic_filter
        The filter applied to postprocessed galaxies: ``"always"`` (the default) accepts all, ``"keep"`` keeps the filter
        of the original model, and a ``<galacticFilter>`` element (in any form accepted for ``extractors``) replaces it. The
        filter of the original model was applied when the file was written, and may need properties (positions in a
        lightcone, for example) which are not restored by postprocessing, so it is not kept by default.
    """
    from Galacticus.Parameters.extract import dereference_parameters, parameters_from_hdf5
    # Parameters which took default values are omitted, so that they are defaulted afresh: some defaults depend on other
    # parameters, and can not be given explicitly alongside them. References between the model's parameters are replaced by
    # copies of the parameters to which they refer, as the target of a reference may lie within a parameter which is replaced
    # below. The `id`s of the model's parameters are kept, so that the extractors and filter given can refer to them.
    parameters = parameters_from_hdf5(str(catalog), include_defaulted=False)
    dereference_parameters(parameters, remove_ids=False)
    # Task.
    task = ET.Element("task")
    task.set("value", "postprocessOutputs")
    _parameter(task, "fileName", Path(catalog).resolve())
    _parameter(task, "outputsGroupName", outputs_group)
    if selection_file is not None:
        _parameter(task, "selectionFileName", Path(selection_file).resolve())
    _parameter(task, "countBatch", int(count_batch))
    _replace(parameters, task)
    # Outputter. Properties of the rebuilt node components would only duplicate (or misrepresent) those in the original file,
    # so are not output.
    outputter = ET.Element("mergerTreeOutputter")
    outputter.set("value", "standard")
    _parameter(outputter, "outputsGroupName", outputs_group)
    _parameter(outputter, "outputComponentProperties", "false")
    _parameter(outputter, "countNodesBuffered", int(count_nodes_buffered))
    _replace(parameters, outputter)
    # Property extractors.
    if isinstance(extractors, (list, tuple)):
        extractors = [_element(extractor, "nodePropertyExtractor") for extractor in extractors]
    else:
        extractors = [_element(extractors, "nodePropertyExtractor")]
    if not extractors:
        raise ValueError("at least one property extractor is required")
    multi = ET.Element("nodePropertyExtractor")
    multi.set("value", "multi")
    ET.SubElement(multi, "nodePropertyExtractor").set("value", "indexPostprocessSource")
    for extractor in extractors:
        multi.append(extractor)
    _replace(parameters, multi)
    # Galactic filter.
    if isinstance(galactic_filter, str) and galactic_filter == "always":
        filter_ = ET.Element("galacticFilter")
        filter_.set("value", "always")
        _replace(parameters, filter_)
    elif not (isinstance(galactic_filter, str) and galactic_filter == "keep"):
        _replace(parameters, _element(galactic_filter, "galacticFilter"))
    # Output file.
    outputFileName = ET.Element("outputFileName")
    outputFileName.set("value", str(Path(output).resolve()))
    _replace(parameters, outputFileName)
    # Resolve references in the extractors and filter given, which may refer to the model's parameters, or to each other.
    dereference_parameters(parameters)
    return parameters


def _output_counts(catalog, outputs_group):
    """Return a dictionary mapping the index of each output in `outputs_group` to its number of galaxies."""
    counts = {}
    with h5py.File(str(catalog), "r") as file:
        if outputs_group not in file:
            raise ValueError(f"'{catalog}' has no group '{outputs_group}'")
        for name, group in file[outputs_group].items():
            if not name.startswith("Output") or not name[6:].isdigit():
                continue
            count = 0
            if "nodeData" in group:
                for dataset in group["nodeData"].values():
                    if isinstance(dataset, h5py.Dataset) and dataset.ndim >= 1:
                        count = dataset.shape[0]
                        break
            counts[int(name[6:])] = count
    return counts


def write_selection(path, selection, *, outputs_group="Outputs", catalog=None):
    """Write the file listing the galaxies to postprocess.

    ``selection`` maps the index of each output to postprocess (the ``i`` in ``Output[i]``) to an array of the indices
    (counting from zero) in the ``nodeData`` datasets of the galaxies to postprocess. These may be in any order, and may repeat
    (a repeated galaxy is postprocessed, and output, each time). Outputs not listed are skipped. Where the file has a single
    output, ``selection`` may instead be just the array of indices.

    If ``catalog`` (the output file to be postprocessed) is given, the selection is checked against it.
    """
    counts = _output_counts(catalog, outputs_group) if catalog is not None else None
    if not isinstance(selection, dict):
        if counts is not None and len(counts) != 1:
            raise ValueError(f"'{catalog}' has {len(counts)} outputs in '{outputs_group}', so the selection must be a "
                             "dictionary mapping output index to galaxy indices")
        selection = {next(iter(counts)) if counts else 1: selection}
    with h5py.File(str(path), "w") as file:
        for indexOutput, indices in selection.items():
            indices = np.asarray(indices)
            if indices.ndim != 1 or (indices.size > 0 and not np.issubdtype(indices.dtype, np.integer)):
                raise ValueError(f"the selection for output {indexOutput} must be a one-dimensional array of integers")
            indices = indices.astype(np.int64)
            if counts is not None:
                if int(indexOutput) not in counts:
                    raise ValueError(f"'{catalog}' has no output '{outputs_group}/Output{indexOutput}'")
                if indices.size > 0 and (indices.min() < 0 or indices.max() >= counts[int(indexOutput)]):
                    raise ValueError(f"selected indices for output {indexOutput} must lie in the range "
                                     f"0 to {counts[int(indexOutput)] - 1}")
            file.create_dataset(f"{outputs_group}/Output{int(indexOutput)}/selection", data=indices)
    return Path(path)


def postprocess(catalog, extractors, output, *, selection=None, outputs_group="Outputs", change_files=(), threads=None,
                count_batch=1024, count_nodes_buffered=1000, galactic_filter="always", parameter_file=None, log_file=None,
                install=None):
    """Postprocess the galaxies of a Galacticus output file, returning the name of the file written.

    ``catalog``, ``extractors``, ``output``, ``outputs_group``, ``count_batch``, ``count_nodes_buffered``, and
    ``galactic_filter`` are as for :func:`build_parameters`. ``selection`` is as for :func:`write_selection`, or the name of
    an existing selection file, or ``None`` to postprocess all galaxies.

    ``change_files`` are Galacticus change files applied (in order) to the parameters, for example to alter a parameter of the
    original model. ``threads`` sets the number of OpenMP threads (by default, as set by ``OMP_NUM_THREADS``). The parameter
    file used is written to ``parameter_file`` if given (otherwise to a temporary file which is removed), and the output of
    Galacticus to ``log_file`` (by default, the output file name with the suffix ``.log``).

    The Galacticus executable is that found by ``galacticus run`` (see ``galacticus info``), unless ``install`` (a
    :class:`galacticus_launcher.paths.Install`) is given. A ``RuntimeError`` is raised if Galacticus fails.
    """
    from . import cli, paths
    if install is None:
        install = paths.resolve()
    cli._ensure(install)
    output = Path(output).resolve()
    log_file = Path(log_file) if log_file is not None else output.with_suffix(".log")
    with tempfile.TemporaryDirectory(prefix="galacticusPostprocess") as work:
        # Selection.
        if selection is None or isinstance(selection, (str, os.PathLike)):
            selection_file = selection
        else:
            selection_file = write_selection(Path(work) / "selection.hdf5", selection, outputs_group=outputs_group,
                                             catalog=catalog)
        # Parameters.
        parameters = build_parameters(catalog, extractors, output, outputs_group=outputs_group, selection_file=selection_file,
                                      count_batch=count_batch, count_nodes_buffered=count_nodes_buffered,
                                      galactic_filter=galactic_filter)
        parameter_file = Path(parameter_file) if parameter_file is not None else Path(work) / "postprocess.xml"
        tree = ET.ElementTree(parameters)
        ET.indent(tree, space="  ", level=0)
        tree.write(str(parameter_file), xml_declaration=True, encoding="UTF-8")
        # Run.
        env = install.environ()
        if threads is not None:
            env["OMP_NUM_THREADS"] = str(int(threads))
        command = [str(install.binary), str(parameter_file), *[str(Path(c).resolve()) for c in change_files]]
        with open(log_file, "w") as log:
            status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
    if status.returncode != 0:
        raise RuntimeError(f"Galacticus failed (exit status {status.returncode}) postprocessing '{catalog}'; see '{log_file}'")
    return output
