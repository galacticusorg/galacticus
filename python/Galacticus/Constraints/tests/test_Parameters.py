"""Tests for extracting best-fit parameters from Galacticus MCMC chain logs.

The regression these guard is a real one: a pipeline stage whose config had been regenerated for
a narrower simulation selection (21 active parameters) was used to label a parameter vector read
from chains run with the full selection (47 parameters). The positional zip silently took the
first 21 columns, so five calibrated values named a different parameter than they held.

Andrew Benson, Claude (13-October-2026)
"""

import pytest

from Galacticus.Constraints.Parameters import (
    chain_parameter_names,
    maximum_parameter_dict,
)

# The parameters the chains in these tests were run with, in column order.
NAMES_CHAIN = [
    'haloMassFunctionParameters/a',
    'haloMassFunctionParameters/normalization',
    'varianceFractionalModelDiscrepancy',
    'haloMassFunctionParameters/perturbationCubeBigMDPL',
    'haloMassFunctionParameters/perturbationCubeHugeMDPL',
    'haloMassFunctionParameters/perturbationCubeMDPL2',
    'haloMassFunctionParameters/massMinimumParticleCountMDPL3',
]
# One row per step; the best log posterior is in the second row.
ROWS = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1],
]
LOG_POSTERIORS  = [-100.0, -10.0]
LOG_LIKELIHOODS = [-11.0, -12.0]


def _write_chain(directory, names=NAMES_CHAIN, rows=ROWS, header=True, first_column=7):
    """Write a single chain log, returning its log file root."""
    root  = str(directory / 'chains')
    lines = []
    if header:
        lines.append('# Simulation state chain file')
        lines.append('# Columns:')
        for i, name in enumerate(names):
            lines.append(f'#   {first_column + i} = Parameter `{name}`')
    for step, row in enumerate(rows):
        columns = [str(step), '0', '0.0', 'F', str(LOG_POSTERIORS[step]), str(LOG_LIKELIHOODS[step])]
        columns += [repr(value) for value in row]
        lines.append(' '.join(columns))
    with open(f'{root}_0000.log', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return root


def _config(root, names):
    """A parsed-config dict declaring *names* as the active parameters."""
    return {
        'posteriorSampleSimulation': {
            'logFileRoot':    {'value': root},
            'modelParameter': [{'value': 'active', 'name': {'value': name}} for name in names],
        }
    }


def test_chain_parameter_names_reads_the_header_in_column_order(tmp_path):
    root = _write_chain(tmp_path)
    assert chain_parameter_names(_config(root, NAMES_CHAIN)) == NAMES_CHAIN


def test_chain_parameter_names_is_empty_without_a_header(tmp_path):
    root = _write_chain(tmp_path, header=False)
    assert chain_parameter_names(_config(root, NAMES_CHAIN)) == []


def test_chain_parameter_names_is_empty_without_a_chain_log(tmp_path):
    assert chain_parameter_names(_config(str(tmp_path / 'absent'), NAMES_CHAIN)) == []


def test_chain_parameter_names_rejects_columns_not_starting_where_parameters_do(tmp_path):
    root = _write_chain(tmp_path, first_column=9)
    with pytest.raises(RuntimeError, match='non-contiguous'):
        chain_parameter_names(_config(root, NAMES_CHAIN))


def test_maximum_parameter_dict_takes_the_best_posterior_row(tmp_path):
    root   = _write_chain(tmp_path)
    values = maximum_parameter_dict(_config(root, NAMES_CHAIN))
    assert values == dict(zip(NAMES_CHAIN, ROWS[1]))


def test_maximum_parameter_dict_can_maximize_the_likelihood_instead(tmp_path):
    root   = _write_chain(tmp_path)
    values = maximum_parameter_dict(_config(root, NAMES_CHAIN), maximum='likelihood')
    assert values == dict(zip(NAMES_CHAIN, ROWS[0]))


def test_maximum_parameter_dict_rejects_an_unknown_maximum(tmp_path):
    root = _write_chain(tmp_path)
    with pytest.raises(ValueError, match='unknown maximum'):
        maximum_parameter_dict(_config(root, NAMES_CHAIN), maximum='mode')


def test_maximum_parameter_dict_matches_a_narrower_config_by_name(tmp_path):
    """The regression: a config declaring a subset of the chains' parameters, at other positions.

    Matched by position, `perturbationCubeMDPL2` would take the value of
    `perturbationCubeBigMDPL` (column 7, the first) rather than its own (column 12).
    """
    root   = _write_chain(tmp_path)
    names  = [
        'haloMassFunctionParameters/a',
        'haloMassFunctionParameters/perturbationCubeMDPL2',
        'haloMassFunctionParameters/massMinimumParticleCountMDPL3',
    ]
    values = maximum_parameter_dict(_config(root, names))
    assert values == {
        'haloMassFunctionParameters/a':                            ROWS[1][0],
        'haloMassFunctionParameters/perturbationCubeMDPL2':         ROWS[1][5],
        'haloMassFunctionParameters/massMinimumParticleCountMDPL3': ROWS[1][6],
    }


def test_maximum_parameter_dict_rejects_a_parameter_absent_from_the_chains(tmp_path):
    root = _write_chain(tmp_path)
    with pytest.raises(RuntimeError, match='exponentMassDetectionMDPL3'):
        maximum_parameter_dict(_config(root, ['haloMassFunctionParameters/exponentMassDetectionMDPL3']))


def test_maximum_parameter_dict_falls_back_to_position_without_a_header(tmp_path):
    root   = _write_chain(tmp_path, header=False)
    values = maximum_parameter_dict(_config(root, NAMES_CHAIN))
    assert values == dict(zip(NAMES_CHAIN, ROWS[1]))


def test_maximum_parameter_dict_rejects_a_count_mismatch_without_a_header(tmp_path):
    root = _write_chain(tmp_path, header=False)
    with pytest.raises(RuntimeError, match='no column header'):
        maximum_parameter_dict(_config(root, NAMES_CHAIN[:3]))
