"""Read and write design files, in the layout written by the Galacticus ``emulatorDesign`` task.

A design file holds a ``design`` group, in the same layout as the ``design`` group of an emulator file (see
:mod:`Galacticus.Emulation.emulatorFile`), and a ``runs`` group listing the runs of the model: one per realization per
design point, each with its change file, output file and (optionally) random number seed.

Designs whose runs are complete can be split into parts (e.g. the even and odd points, to measure the systematic error of
the emulators from the spread between emulators trained on each), or pooled (e.g. two stages run with the same base
parameter files, to train one emulator on all their points), reusing the outputs of their runs.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field, replace

import h5py
import numpy as np

from Galacticus.Emulation.emulatorFile import Design, Prior, _decode, _encode

__all__ = [
    'FORMAT_NAME',
    'FORMAT_VERSION',
    'DesignFile',
    'parse_descriptor',
    'pool_designs',
    'prior_cumulative',
    'read_design',
    'resolve_paths',
    'split_design',
    'serialize_descriptor',
    'write_design',
]

FORMAT_NAME = 'galacticusDesign'
FORMAT_VERSION = 1


def _convert(text):
    """Convert a descriptor value to a number if it looks like one, else leave it as a string."""
    try:
        return float(text)
    except ValueError:
        return text


def parse_descriptor(descriptor):
    """Parse a serialized parameter set, as written by Galacticus' ``inputParameters%serializeToString()``.

    The format is ``name:value`` pairs separated by ``_``, with any sub-parameters of a parameter following its value in
    braces, e.g. ``x0:250.0_sigma:0.5_randomNumberGenerator:GSL{seed:42}``. Returns a dictionary mapping each name to its
    value (converted to ``float`` where possible). A parameter with sub-parameters maps to a pair ``(value, dict)``.
    """
    result = {}
    position = 0
    length = len(descriptor)
    while position < length:
        colon = descriptor.index(':', position)
        name = descriptor[position:colon]
        position = colon + 1
        # The value runs to the next separator or opening brace at this level.
        end = position
        while end < length and descriptor[end] not in '_{':
            end += 1
        value = _convert(descriptor[position:end])
        position = end
        if position < length and descriptor[position] == '{':
            depth = 1
            close = position + 1
            while depth > 0:
                if descriptor[close] == '{':
                    depth += 1
                elif descriptor[close] == '}':
                    depth -= 1
                close += 1
            value = (value, parse_descriptor(descriptor[position + 1:close - 1]))
            position = close
        result[name] = value
        if position < length:
            if descriptor[position] != '_':
                raise ValueError(f"malformed descriptor '{descriptor}' at position {position}")
            position += 1
    return result


def serialize_descriptor(parameters):
    """Serialize a (flat) set of parameters as ``name:value`` pairs separated by ``_``, the inverse of
    :func:`parse_descriptor` for parameters without sub-parameters. Numbers are written in full precision."""
    return '_'.join(f'{name}:{value!r}' if isinstance(value, float) else f'{name}:{value}' for name, value in parameters.items())


@dataclass
class DesignFile:
    """The content of a design file.

    ``design`` is a :class:`Galacticus.Emulation.emulatorFile.Design`; the ``descriptor`` of each prior is kept in its
    parameters, alongside the parsed values. ``runs`` maps each dataset of the ``runs`` group (``pointIndex``,
    ``realizationIndex``, ``changeFileName``, ``outputFileName``, and ``seed`` if set) to its values.
    """

    design: Design
    runs: dict
    attributes: dict = field(default_factory=dict)

    @property
    def count_runs(self):
        return len(self.runs['pointIndex'])


def read_design(path):
    """Read a design file written by the ``emulatorDesign`` task, returning a :class:`DesignFile`."""
    with h5py.File(path, 'r') as file:
        attributes = {name: _decode(value) for name, value in file.attrs.items()}
        if attributes.get('format') != FORMAT_NAME:
            raise ValueError(f"'{path}' is not a Galacticus design file")
        if attributes.get('formatVersion') != FORMAT_VERSION:
            raise ValueError(f"'{path}' has format version {attributes.get('formatVersion')}; this reader supports version {FORMAT_VERSION}")
        group = file['design']
        names = _decode(group['parameterNames'][()])
        priors = []
        for i in range(len(names)):
            prior_attributes = {name: _decode(value) for name, value in group['priors'][f'prior{i + 1}'].attrs.items()}
            parameters = parse_descriptor(prior_attributes['descriptor'])
            parameters['descriptor'] = prior_attributes['descriptor']
            priors.append(Prior(prior_attributes['class'], parameters))
        design_attributes = {name: _decode(value) for name, value in group.attrs.items() if name != 'comment'}
        design = Design(
            names=names,
            priors=priors,
            mappers=_decode(group['mappers'][()]),
            quantiles=group['quantiles'][()],
            values=group['values'][()],
            attributes=design_attributes,
        )
        runs = {}
        for name, dataset in file['runs'].items():
            data = dataset[()]
            runs[name] = _decode(data) if data.dtype.kind == 'S' else np.asarray(data)
    design.validate()
    return DesignFile(design=design, runs=runs, attributes=attributes)


def write_design(path, design_file, overwrite=False):
    """Write a :class:`DesignFile` in the layout written by the ``emulatorDesign`` task (for example, to import a design
    generated elsewhere). The ``descriptor`` of each prior is taken from its parameters if present there, and is otherwise
    serialized from them."""
    design = design_file.design
    design.validate()
    count_points = design.quantiles.shape[0]
    with h5py.File(path, 'w' if overwrite else 'w-') as file:
        attributes = dict(design_file.attributes)
        attributes.update({'format': FORMAT_NAME, 'formatVersion': np.int32(FORMAT_VERSION)})
        for name, value in attributes.items():
            file.attrs[name] = _encode(value)
        group = file.create_group('design')
        design_attributes = {'countPoints': np.int32(count_points), 'realizationsPerPoint': np.int32(1), 'seedPerPoint': np.int32(0)}
        design_attributes.update(design.attributes)
        for name, value in design_attributes.items():
            group.attrs[name] = _encode(value)
        group.create_dataset('parameterNames', data=_encode(list(design.names)))
        group.create_dataset('mappers', data=_encode(list(design.mappers)))
        group.create_dataset('quantiles', data=design.quantiles)
        group.create_dataset('values', data=design.values)
        priors = group.create_group('priors')
        for i, prior in enumerate(design.priors):
            parameters = {name: value for name, value in prior.parameters.items() if name != 'descriptor'}
            prior_group = priors.create_group(f'prior{i + 1}')
            prior_group.attrs['class'] = _encode(prior.class_name)
            prior_group.attrs['descriptor'] = _encode(prior.parameters.get('descriptor', serialize_descriptor(parameters)))
        runs = file.create_group('runs')
        for name, values in design_file.runs.items():
            values = list(values) if isinstance(values, (list, tuple)) else values
            if isinstance(values, list) and values and isinstance(values[0], str):
                runs.create_dataset(name, data=_encode(values))
            else:
                runs.create_dataset(name, data=np.asarray(values))


def prior_cumulative(prior, values):
    """The cumulative probability of ``values`` under a prior (a :class:`Galacticus.Emulation.emulatorFile.Prior`), as
    computed by the Galacticus ``distributionFunction1D`` class of that name: ``uniform``, ``logUniform``, ``normal`` or
    ``logNormal`` (with ``mean`` and ``variance``, or ``x0`` and ``sigma``, and optional ``limitLower`` and ``limitUpper``).
    Values outside the limits of the prior are clipped to 0 or 1."""
    from scipy.special import ndtr
    values = np.asarray(values, dtype=float)
    parameters = {name: float(value) for name, value in prior.parameters.items() if name != 'descriptor'}
    if prior.class_name == 'uniform':
        lower, upper = parameters['limitLower'], parameters['limitUpper']
        return np.clip((values - lower) / (upper - lower), 0.0, 1.0)
    if prior.class_name == 'logUniform':
        lower, upper = parameters['limitLower'], parameters['limitUpper']
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.clip(np.log(values / lower) / np.log(upper / lower), 0.0, 1.0)
    if prior.class_name in ('normal', 'logNormal'):
        if prior.class_name == 'normal':
            mean, variance = parameters['mean'], parameters['variance']
            x, transform = values, (lambda limit: limit)
        else:
            if 'x0' in parameters:
                mean, variance = np.log(parameters['x0']), parameters['sigma'] ** 2
            else:
                variance = np.log(parameters['variance'] / parameters['mean'] ** 2 + 1.0)
                mean = np.log(parameters['mean']) - 0.5 * variance
            with np.errstate(divide='ignore', invalid='ignore'):
                x = np.log(values)
            transform = np.log
        cdf = lambda y: ndtr((y - mean) / np.sqrt(variance))  # noqa: E731
        cdfLower = cdf(transform(parameters['limitLower'])) if 'limitLower' in parameters else 0.0
        cdfUpper = cdf(transform(parameters['limitUpper'])) if 'limitUpper' in parameters else 1.0
        return np.clip((cdf(x) - cdfLower) / (cdfUpper - cdfLower), 0.0, 1.0)
    raise NotImplementedError(f"prior class '{prior.class_name}' is not supported")


def resolve_paths(design_file, directory):
    """Return a copy of a :class:`DesignFile` with the change and output file names of its runs made absolute, resolving
    relative names against ``directory`` (the directory from which the runs were made)."""
    runs = dict(design_file.runs)
    for name in ('changeFileName', 'outputFileName'):
        if name in runs:
            runs[name] = [path if os.path.isabs(path) else os.path.abspath(os.path.join(directory, path)) for path in runs[name]]
    return replace(design_file, runs=runs)


def _select_runs(runs, selected):
    return {name: ([value for value, keep in zip(values, selected) if keep] if isinstance(values, list)
                   else np.asarray(values)[selected]) for name, values in runs.items()}


def split_design(design_file, parts=2):
    """Split a :class:`DesignFile` into ``parts`` designs, by design point: part ``k`` holds the points whose index modulo
    ``parts`` is ``k`` (for two parts, the even and the odd points), with all their runs. A low-discrepancy design split
    this way gives parts which each still fill the space evenly. Returns a list of :class:`DesignFile`."""
    if parts < 2:
        raise ValueError('a design must be split into at least two parts')
    design = design_file.design
    points = np.arange(design.quantiles.shape[0])
    if points.size < parts:
        raise ValueError(f'a design of {points.size} points can not be split into {parts} parts')
    pointIndex = np.asarray(design_file.runs['pointIndex'])
    result = []
    for part in range(parts):
        keep = points[points % parts == part]
        remap = np.full(points.size, -1, dtype=np.int64)
        remap[keep] = np.arange(keep.size)
        runs = _select_runs(design_file.runs, np.isin(pointIndex, keep))
        runs['pointIndex'] = remap[np.asarray(runs['pointIndex'])].astype(np.int32)
        attributes = dict(design.attributes)
        attributes.update({'countPoints': np.int32(keep.size), 'splitPart': np.int32(part), 'splitParts': np.int32(parts)})
        result.append(replace(design_file, runs=runs,
                              design=replace(design, quantiles=design.quantiles[keep], values=design.values[keep], attributes=attributes)))
    return result


def pool_designs(design_files, priors=None):
    """Pool several :class:`DesignFile` into one, holding all of their points and runs.

    The designs must vary the same parameters, with the same mappers, and must have been run with the same base parameter
    files and executable (which can not be checked here). The quantiles of every point are recomputed from its parameter
    values under ``priors`` (a list of :class:`Galacticus.Emulation.emulatorFile.Prior`, one per parameter; by default the
    priors of the first design), which must contain every point: the union of the boxes of the designs, for example. The
    pooled design is then a valid training design for an emulated likelihood with those priors. The run file names should
    first be made absolute (see :func:`resolve_paths`) if the designs were run from different directories.
    """
    if len(design_files) < 2:
        raise ValueError('at least two designs must be given to pool')
    first = design_files[0].design
    for design_file in design_files[1:]:
        if list(design_file.design.names) != list(first.names):
            raise ValueError('designs to pool must vary the same parameters, in the same order')
        if list(design_file.design.mappers) != list(first.mappers):
            raise ValueError('designs to pool must have the same mappers')
        if set(design_file.runs) != set(design_files[0].runs):
            raise ValueError('designs to pool must record the same information about their runs (e.g. seeds)')
    priors = list(first.priors) if priors is None else list(priors)
    if len(priors) != len(first.names):
        raise ValueError(f'{len(priors)} priors given for {len(first.names)} parameters')
    values = np.concatenate([design_file.design.values for design_file in design_files])
    quantiles = np.column_stack([prior_cumulative(prior, values[:, k]) for k, prior in enumerate(priors)])
    for k, prior in enumerate(priors):
        outside = (quantiles[:, k] <= 0.0) | (quantiles[:, k] >= 1.0)
        if outside.any():
            raise ValueError(f"parameter '{first.names[k]}': {int(outside.sum())} points lie outside the support of the given prior")
    runs, offset = {name: [] for name in design_files[0].runs}, 0
    for design_file in design_files:
        for name, value in design_file.runs.items():
            value = list(value) if isinstance(value, list) else list(np.asarray(value))
            runs[name] += [entry + offset for entry in value] if name == 'pointIndex' else value
        offset += design_file.design.quantiles.shape[0]
    for name in runs:
        if name in ('pointIndex', 'realizationIndex'):
            runs[name] = np.array(runs[name], dtype=np.int32)
        elif name == 'seed':
            runs[name] = np.array(runs[name], dtype=np.int64)
    if 'seed' in runs:
        seeds, counts = np.unique(runs['seed'], return_counts=True)
        if np.any((counts > 1) & (seeds >= 0)):
            warnings.warn(f'{int(np.sum(counts[seeds >= 0] - 1))} runs of the pooled designs share a random number seed with another '
                          'run, so their noise is correlated (give each design its own `seedBase`)', stacklevel=2)
    seedPerPoint = {int(design_file.design.attributes.get('seedPerPoint', 0)) for design_file in design_files}
    attributes = {'type': 'pooled', 'countPoints': np.int32(values.shape[0]), 'pooledDesigns': np.int32(len(design_files)),
                  'realizationsPerPoint': np.int32(max(int(design_file.design.attributes.get('realizationsPerPoint', 1))
                                                       for design_file in design_files)),
                  'seedPerPoint': np.int32(max(seedPerPoint))}
    design = Design(names=list(first.names), priors=priors, mappers=list(first.mappers), quantiles=quantiles, values=values,
                    attributes=attributes)
    return DesignFile(design=design, runs=runs, attributes=dict(design_files[0].attributes))
