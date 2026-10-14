"""Provides utilities for working with Galacticus MCMC parameter chains.

Andrew Benson (ported to Python 2026)
"""

import os
import re

import numpy as np

from List.ExtraUtils import as_array

__all__ = [
    # Canonical snake_case names.
    'log_file_root',
    'maximum_posterior_parameter_vector',
    'maximum_likelihood_parameter_vector',
    'maximum_parameter_dict',
    'chain_parameter_names',
    'parameter_names',
    # camelCase compatibility aliases (defined further down in this module).
    'logFileRoot',
    'maximumPosteriorParameterVector',
    'maximumLikelihoodParameterVector',
    'parameterNames',
]

# The 1-based column of the first parameter in a chain log (columns 1-6 are bookkeeping).
_PARAMETER_COLUMN_FIRST = 7


def log_file_root(config, options=None):
    """Return the chain log file root path.

    Respects an options['chainRoot'] override; otherwise reads from the config
    dict's posteriorSampleSimulation/logFileRoot/value.
    """
    options = options or {}
    return (options.get('chainRoot')
            or config['posteriorSampleSimulation']['logFileRoot']['value'])


def _maximum_parameter_vector(config, select_column, options=None):
    """Scan all chain log files and return the parameter vector for the row
    with the highest value in select_column (4 = log-posterior, 5 = log-likelihood).

    Chain log format (space-separated columns):
      col 0   step number
      col 1-2 (unused)
      col 3   convergence flag 'T'/'F'
      col 4   log-posterior
      col 5   log-likelihood
      col 6+  parameter values

    Lines starting with '#' are skipped. options['burnCount'] rows are skipped
    from the start of each file. options['chain'] restricts to a single chain index.

    Returns (np.ndarray of parameter values, float best_value).
    """
    options    = options or {}
    root       = log_file_root(config, options)
    burn_count = int(options.get('burnCount', 0))
    chain_only = options.get('chain')

    best_val    = None
    best_params = None

    i = 0
    while True:
        chain_file = f'{root}_{i:04d}.log'
        if not os.path.exists(chain_file):
            break
        if chain_only is not None and chain_only != 'all' and int(chain_only) != i:
            i += 1
            continue
        step = 0
        with open(chain_file) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                step += 1
                cols = line.split()
                if len(cols) < 7:
                    continue
                if step <= burn_count:
                    continue
                val = float(cols[select_column])
                if best_val is None or val > best_val:
                    best_val    = val
                    best_params = [float(c) for c in cols[6:]]
        i += 1

    if best_params is None:
        raise RuntimeError('No acceptable states found in chain log files')

    # particleSwarm simulation: only the first half of columns are position parameters.
    sim_type = config['posteriorSampleSimulation'].get('value', '')
    if sim_type == 'particleSwarm':
        best_params = best_params[:len(best_params) // 2]

    return np.array(best_params), best_val


def maximum_posterior_parameter_vector(config, options=None):
    """Return (parameter vector, log-posterior value) for the best-posterior state."""
    return _maximum_parameter_vector(config, 4, options)


def maximum_likelihood_parameter_vector(config, options=None):
    """Return (parameter vector, log-likelihood value) for the best-likelihood state."""
    return _maximum_parameter_vector(config, 5, options)


def parameter_names(config):
    """Return the list of active parameter names from posteriorSampleSimulation.

    Each name comes from a <modelParameter value="active"><name value="..."/></modelParameter>
    block in the config XML.
    """
    params = as_array(config['posteriorSampleSimulation']['modelParameter'])
    return [p['name']['value'] for p in params if p.get('value') == 'active']


def chain_parameter_names(config, options=None):
    """Return the parameter names recorded in the chain logs' own header, in column order.

    Galacticus writes a header to each chain log naming every column, e.g.::

        #    7 = Parameter `haloMassFunctionParameters/a`

    These names are authoritative: they describe the chain that was actually run, whereas the
    config may since have been regenerated with a different set of active parameters. The
    returned list is indexed as the parameter vector is, so element *i* names the vector's
    element *i*.

    Returns an empty list for a chain log carrying no such header (an older log), leaving the
    caller to fall back to matching by position.
    """
    options    = options or {}
    chain_file = f'{log_file_root(config, options)}_0000.log'
    if not os.path.exists(chain_file):
        return []
    columns = {}
    with open(chain_file) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            match = re.match(r'#\s*(\d+)\s*=\s*Parameter\s+`([^`]+)`', line)
            if match:
                columns[int(match.group(1))] = match.group(2)
    if not columns:
        return []
    expected = list(range(_PARAMETER_COLUMN_FIRST, _PARAMETER_COLUMN_FIRST + len(columns)))
    if sorted(columns) != expected:
        raise RuntimeError(f'chain log {chain_file} names a non-contiguous set of parameter columns')
    return [columns[column] for column in expected]


def maximum_parameter_dict(config, maximum='posterior', options=None):
    """Return the {name: value} best-fit parameters of the config's active parameters.

    *maximum* selects what is maximized over the chains: ``'posterior'`` (the default) or
    ``'likelihood'``.

    Values are matched to names through the chain logs' own column header rather than by
    position, so a config whose active parameter list differs from the one the chains were run
    with - for instance one regenerated for a narrower simulation selection - cannot silently
    mislabel them. A config-active parameter the chains do not carry is an error rather than a
    value taken from the wrong column.
    """
    if maximum == 'posterior':
        values, _ = maximum_posterior_parameter_vector(config, options)
    elif maximum == 'likelihood':
        values, _ = maximum_likelihood_parameter_vector(config, options)
    else:
        raise ValueError(f"unknown maximum '{maximum}': expected 'posterior' or 'likelihood'")
    names       = parameter_names(config)
    names_chain = chain_parameter_names(config, options)
    if not names_chain:
        # An older chain log with no column header: fall back to matching by position, which is
        # only meaningful if the chains and the config describe the same number of parameters.
        if len(values) != len(names):
            raise RuntimeError(
                f'the chain logs hold {len(values)} parameters but the config declares '
                f'{len(names)} active, and the logs carry no column header to match them by name'
            )
        return {name: float(values[i]) for i, name in enumerate(names)}
    # Keep the first column of each name: a particleSwarm chain repeats the parameter names for
    # its velocity columns, which `_maximum_parameter_vector` has already trimmed away.
    indices = {}
    for index, name in enumerate(names_chain):
        if index < len(values):
            indices.setdefault(name, index)
    missing = [name for name in names if name not in indices]
    if missing:
        raise RuntimeError(
            'the chain logs do not carry these parameters, which the config declares active: '
            + ', '.join(missing)
        )
    return {name: float(values[indices[name]]) for name in names}


# ---------------------------------------------------------------------------
# Public legacy camelCase API aliases
# ---------------------------------------------------------------------------

logFileRoot                     = log_file_root
maximumPosteriorParameterVector = maximum_posterior_parameter_vector
maximumLikelihoodParameterVector = maximum_likelihood_parameter_vector
parameterNames                  = parameter_names
