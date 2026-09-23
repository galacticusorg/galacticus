"""Read design files written by the Galacticus ``emulatorDesign`` task.

A design file holds a ``design`` group, in the same layout as the ``design`` group of an emulator file (see
:mod:`Galacticus.Emulation.emulatorFile`), and a ``runs`` group listing the runs of the model: one per realization per
design point, each with its change file, output file and (optionally) random number seed.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import h5py
import numpy as np

from Galacticus.Emulation.emulatorFile import Design, Prior, _decode

__all__ = [
    'FORMAT_NAME',
    'FORMAT_VERSION',
    'DesignFile',
    'parse_descriptor',
    'read_design',
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
