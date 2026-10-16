"""Read and write design files, in the layout written by the Galacticus ``emulatorDesign`` task.

A design file holds a ``design`` group, in the same layout as the ``design`` group of an emulator file (see
:mod:`Galacticus.Emulation.emulatorFile`), and a ``runs`` group listing the runs of the model: one per realization per
design point, each with its change file, output file and (optionally) random number seed.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import h5py
import numpy as np

from Galacticus.Emulation.emulatorFile import Design, Prior, _decode, _encode

__all__ = [
    'FORMAT_NAME',
    'FORMAT_VERSION',
    'DesignFile',
    'parse_descriptor',
    'read_design',
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
