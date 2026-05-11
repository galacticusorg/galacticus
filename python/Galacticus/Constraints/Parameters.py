"""Provides utilities for working with Galacticus MCMC parameter chains.

Python port of the subset of perl/Galacticus/Constraints/Parameters.pm
used by the dark-matter constraint pipeline.
Andrew Benson (ported to Python 2026)
"""

import os

import numpy as np

from List.ExtraUtils import as_array

__all__ = [
    # Canonical snake_case names.
    'log_file_root',
    'maximum_posterior_parameter_vector',
    'maximum_likelihood_parameter_vector',
    'parameter_names',
    # camelCase compatibility aliases (defined further down in this module).
    'logFileRoot',
    'maximumPosteriorParameterVector',
    'maximumLikelihoodParameterVector',
    'parameterNames',
]


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


# ---------------------------------------------------------------------------
# Public camelCase aliases (matches Perl module naming convention)
# ---------------------------------------------------------------------------

logFileRoot                     = log_file_root
maximumPosteriorParameterVector = maximum_posterior_parameter_vector
maximumLikelihoodParameterVector = maximum_likelihood_parameter_vector
parameterNames                  = parameter_names
