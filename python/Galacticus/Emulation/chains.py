"""Read the chains of a Galacticus posterior sampling simulation, and compare posteriors.

The comparison measures the systematic error of emulators: emulators trained on parts of a design (e.g. its even and odd
points; see :func:`Galacticus.Emulation.design.split_design`) give independent estimates of the posterior, and the shifts of
their medians from that of the emulators trained on the whole design measure how far the posterior is limited by the
emulators rather than by the data.

Andrew Benson, Claude (28-October-2026).
"""

from __future__ import annotations

import glob
import re

import numpy as np

__all__ = [
    'read_chains',
    'summarize',
    'compare_posteriors',
]

# The 0-based columns of the convergence flag and of the first parameter in a chain log.
_COLUMN_CONVERGED = 3
_COLUMN_PARAMETER_FIRST = 6


def read_chains(root, converged_only=True, burn_fraction=0.5):
    """Read the chain logs ``<root>_NNNN.log`` of a posterior sampling simulation.

    Returns ``(names, states)``: the parameter names (from the header of the logs) and an array of shape ``(N, d)`` of the
    parameter values of the states used. If ``converged_only`` is true and any state is flagged as converged, only the
    converged states are used; otherwise the first ``burn_fraction`` of the states of each chain are discarded.
    """
    files = sorted(path for path in glob.glob(f'{glob.escape(root)}_*.log') if re.fullmatch(re.escape(root) + r'_\d{4}\.log', path))
    if not files:
        raise FileNotFoundError(f"no chain logs found matching '{root}_NNNN.log'")
    names = {}
    with open(files[0]) as handle:
        for line in handle:
            if not line.startswith('#'):
                break
            match = re.match(r'#\s*(\d+)\s*=\s*Parameter\s+`([^`]+)`', line)
            if match:
                names[int(match.group(1))] = match.group(2)
    if not names:
        raise ValueError(f"chain log '{files[0]}' does not name its parameters in its header")
    columns = sorted(names)
    chains = []
    for path in files:
        rows = [line.split() for line in open(path) if not line.startswith('#') and line.strip()]
        if not rows:
            continue
        converged = np.array([row[_COLUMN_CONVERGED] == 'T' for row in rows])
        values = np.array([[float(value) for value in row[_COLUMN_PARAMETER_FIRST:_COLUMN_PARAMETER_FIRST + len(columns)]] for row in rows])
        chains.append((converged, values))
    if not chains:
        raise ValueError(f"the chain logs of '{root}' hold no states")
    if converged_only and any(converged.any() for converged, _ in chains):
        states = np.concatenate([values[converged] for converged, values in chains])
    else:
        states = np.concatenate([values[int(burn_fraction * values.shape[0]):] for _, values in chains])
    return [names[column] for column in columns], states


def summarize(states):
    """The median, and 16th and 84th percentiles, of each parameter: three arrays."""
    lower, median, upper = np.percentile(states, [15.865525393145708, 50.0, 84.13447460685429], axis=0)
    return median, lower, upper


def compare_posteriors(reference, others):
    """Compare posteriors with a reference posterior.

    ``reference`` and each of ``others`` are ``(names, states)`` (as returned by :func:`read_chains`); the others must
    hold the parameters of the reference (in any order). Returns a list with, for each parameter of the reference, a
    dictionary of its ``name``, the reference ``median`` and its half-width ``sigma`` (half the 16-84 percentile range), the
    ``shifts`` of the median of each other posterior from the reference median, their ``sigmas``, and the ``systematic``
    (the largest absolute shift) and its ratio to ``sigma``.
    """
    names, states = reference
    median, lower, upper = summarize(states)
    summaries = []
    for otherNames, otherStates in others:
        missing = [name for name in names if name not in otherNames]
        if missing:
            raise ValueError(f"posterior lacks parameters: {', '.join(missing)}")
        order = [otherNames.index(name) for name in names]
        summaries.append(summarize(otherStates[:, order]))
    result = []
    for k, name in enumerate(names):
        sigma = 0.5 * (upper[k] - lower[k])
        shifts = [float(summary[0][k] - median[k]) for summary in summaries]
        systematic = max(abs(shift) for shift in shifts) if shifts else 0.0
        result.append({
            'name': name,
            'median': float(median[k]),
            'sigma': float(sigma),
            'shifts': shifts,
            'sigmas': [float(0.5 * (summary[2][k] - summary[1][k])) for summary in summaries],
            'systematic': float(systematic),
            'ratio': float(systematic / sigma) if sigma > 0.0 else float('inf'),
        })
    return result
