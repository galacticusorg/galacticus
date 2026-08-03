#!/usr/bin/env python3
"""Resolution-convergence harness for the dark-matter structure models.

Arm C of the harness described in `darkMatterPipelineReview.md` §10: an MDPL2
*self-convergence* test. It runs the current spin/concentration model at a ladder
of fixed mass resolutions and measures how far the predicted median c(M), the
spin distribution P(λ), and the Johnson et al. (2021) energy-model failure rate
drift as resolution changes. Flat curves versus resolution are the pass
criterion. This is the acceptance test for the physics improvements planned in
`darkMatterPipelineNextStages.md` §4–5, run against the *current* (unimproved)
code to baseline the drift.

Self-convergence only: the model is compared against itself at different
resolutions, not against N-body targets (that is Arm A, Symphony Halo004). No
target files are needed — the concentration and spin analyses are configured
entirely from inline parameters, so this sidesteps the 2-D/1-D target-format
reconciliation (nextStages §2.2 step 2).

Andrew Benson (2026).
"""

import argparse
import glob
import os
import re

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from lxml import etree

import queueManager
from queueManager import submit_jobs

# ---------------------------------------------------------------------------
# Harness configuration
# ---------------------------------------------------------------------------

# The resolution ladder: absolute mass resolution (M_res, in M_sun), swept via
# `mergerTreeMassResolutionFixed`. Fixed (not fractional) resolution holds M_res
# constant across tree mass, removing the tree-mass confound that `scaled`
# introduces. The finest rung is a few times the MDPL2 particle mass; the
# coarsest is ~1 decade below the lowest measured halo mass. At the default
# `factorMassResolution = 100` the Johnson2021 energy model applies only above
# 100*M_res, so this ladder sweeps the model from "fall-back seeded" (coarse) to
# "energy-model applied throughout the measured range" (fine).
MASS_RESOLUTION_LADDER = [1.0e12, 3.16e11, 1.0e11, 3.16e10, 1.0e10]

# --- Second axis: the Johnson2021 seeding threshold ------------------------
#
# `factorMassResolution` (f) is the Johnson2021 parameter that decides *where*
# the energy model hands over to the fall-back: the energy model is applied only
# to nodes with mass > f*M_res; below that the whole sub-branch is seeded from
# the fall-back c(M) relation plus correlated scatter (Johnson2021.F90:487). So
# the seeding mass is
#
#     M_seed = f * M_res
#
# and the M_res ladder above (run at the default f = 100) necessarily moved
# M_seed and tree depth *together*. This second axis separates them, because f
# changes M_seed while leaving tree construction — and therefore tree depth and
# the resolved merger spectrum — completely untouched.
#
# The diagnostic is a collapse test: plot every run against M_seed. If the runs
# collapse onto a single curve, the resolution drift seen on axis 1 is entirely a
# seeding-threshold effect and is cured by pushing f*M_res down (or by §5.2's
# self-consistent seeding). If they do not collapse, tree depth carries drift of
# its own, and the §5.3 per-merger scatter is implicated.
#
# f = 1 is the "no fall-back" limit: every node with progenitors gets the energy
# model, and only leaf nodes are seeded. It is also where the §5.1 positive-energy
# failure path is most likely to be exercised, since the least-resolved nodes are
# the ones newly handed to the energy model.
FACTOR_MASS_RESOLUTION_DEFAULT = 100.0

# The f ladder is swept at this (middle) rung, whose f = 100 run already exists
# from axis 1 and serves as the free anchor point.
MASS_RESOLUTION_FACTOR_AXIS = 1.0e11
FACTOR_MASS_RESOLUTION_LADDER = [1.0, 3.0, 10.0, 30.0]

# Degeneracy checks: (M_res, f) pairs chosen so that f*M_res reproduces a seeding
# mass already sampled at a *different* M_res. Together with the runs above these
# give a triple at M_seed = 1e12 -- (1e10, 100), (1e11, 10), (1e12, 1) -- and a
# pair at M_seed = 1e13 -- (1e11, 100), (1e12, 10). Within each set M_seed is
# identical and only tree depth differs, so any spread is pure depth dependence.
DEGENERACY_POINTS = [(1.0e12, 1.0), (1.0e12, 10.0)]

# The M_res ladder is run at each of these f values. f = 100 is the code default
# and gives the original (§10.7) baseline; f = 1 repeats the ladder with fall-back
# seeding minimized, so that the drift which remains is attributable to tree depth
# acting on the energy model itself rather than to the moving seeding threshold.
# The §10.9 result makes the f = 100 ladder hard to interpret on its own: at
# f = 100 the lowest mass bin lies entirely below M_seed at every rung, so its
# apparent convergence measured the fall-back relation, not the energy model.
LADDER_FACTORS = [FACTOR_MASS_RESOLUTION_DEFAULT, 1.0]

# MDPL2 particle mass (M_sun), from the box metadata. Used for the N-body error
# model and reported for reference against the ladder floor.
MASS_PARTICLE = 2.228124538881511e9

# Tree mass range and sampling. Trees are built (and their z=0 roots measured)
# across this range; the lowest bin edge equals `MASS_TREE_MINIMUM`.
MASS_TREE_MINIMUM = 3.16e12   # 10^12.5
MASS_TREE_MAXIMUM = 1.00e15   # 10^15.0
TREES_PER_DECADE  = 300       # diagnostic sampling (base config uses 2000)

# Concentration distribution: one analysis per 0.5-dex mass bin. Concentration
# binning matches the target-construction config
# (`concentrationDistributionFunctionCompute.xml`): 1 <= c <= 300 at 5/decade.
MASS_BINS_LOG10       = [(12.5, 13.0), (13.0, 13.5), (13.5, 14.0),
                         (14.0, 14.5), (14.5, 15.0)]
CONCENTRATION_MINIMUM = 1.0
CONCENTRATION_MAXIMUM = 300.0
CONCENTRATION_PER_DECADE = 5

# Spin distribution: a single analysis over the full measured mass range. Spin
# binning matches the MDPL2 spin target (~5/decade over 0.0013 <= λ <= 0.79).
SPIN_MINIMUM       = 1.26e-3
SPIN_MAXIMUM       = 0.794
SPIN_PER_DECADE    = 5
LOGNORMAL_RANGE    = 1.420    # spin measurement-error convolution width (base config)

# Held-fixed analysis constants. These affect the absolute level of the
# distributions but not the *drift* with resolution (the quantity of interest),
# since they are identical across every rung.
TIME_RECENT                          = 0.0   # no recent-major-merger exclusion (revisit for Arm A)
PARTICLE_COUNT_MINIMUM               = 300    # integer parameter
ENERGY_ESTIMATE_PARTICLE_COUNT_MAX   = 1000   # integer parameter

REDSHIFT = 0.0


def _run_grid(ladder_only=False):
    """The full set of (M_res, f) run points, axis 1 first.

    Axis 1 is the M_res ladder, run at each f in `LADDER_FACTORS`; axis 2 is the f ladder at
    `MASS_RESOLUTION_FACTOR_AXIS`, plus the degeneracy points. Duplicates are
    dropped so the shared anchor point is run only once.

    `ladder_only` restricts the grid to the plain M_res ladder at the default seeding
    threshold (the f = 100 axis-1 points). This is the right subset when the swept
    quantity affects only spin (e.g. the Vitvitska sub-resolution method), since the
    factorMassResolution axis touches only the Johnson2021 concentration seeding and
    would just reproduce the same spins redundantly."""
    if ladder_only:
        return [(m, FACTOR_MASS_RESOLUTION_DEFAULT) for m in MASS_RESOLUTION_LADDER]
    points = [(m, f) for f in LADDER_FACTORS for m in MASS_RESOLUTION_LADDER]
    points += [(MASS_RESOLUTION_FACTOR_AXIS, f) for f in FACTOR_MASS_RESOLUTION_LADDER]
    points += DEGENERACY_POINTS
    unique = []
    for point in points:
        if point not in unique:
            unique.append(point)
    return unique

# Regexes for the Johnson2021 energy-model statistics reported at output-file
# close (see `Johnson2021_statistics.F90`).
_FAILURE_RE      = re.compile(r'energy model failed.*?for\s+(\d+)\s+of the\s+(\d+)\s+nodes')
_NEVER_APPLIED_RE = re.compile(r'energy model was constructed but never applied')


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Resolution-convergence harness (Arm C: MDPL2 self-convergence).'
    )
    parser.add_argument('--pipelinePath',    required=True,
                        help='Path to pipeline configuration files (holds spinConcentrationBaseMDPL2.xml)')
    parser.add_argument('--outputDirectory', required=True,
                        help='Directory for generated parameter files, model outputs, and plots')
    parser.add_argument('--force',           default='no',
                        help='Re-run models even if outputs exist (default: no)')
    parser.add_argument('--plotOnly',        action='store_true',
                        help='Skip generation and submission; only (re)build plots from existing outputs')
    parser.add_argument('--treesPerDecade',  default=TREES_PER_DECADE, type=int,
                        help=f'Trees per decade of tree mass (default: {TREES_PER_DECADE})')
    parser.add_argument('--subresolutionMethod', default=None,
                        help='Override the Vitvitska sub-resolution angular-momentum method '
                             '(resolutionScaled | massScaled | original); default: inherit the base config.')
    parser.add_argument('--ladderOnly',      action='store_true',
                        help='Run only the M_res ladder at the default seeding threshold '
                             '(skip the factorMassResolution axis); use when sweeping a spin-only quantity.')
    parser.add_argument('--ppn',             default=16, type=int,
                        help='Processors (OpenMP threads) per model run (default: 16)')
    parser.add_argument('--nodes',           default=1, type=int,
                        help='Nodes per model run (default: 1)')
    parser.add_argument('--walltime',        default='24:00:00',
                        help='Walltime per model run (default: 24:00:00)')
    parser.add_argument('--memory',          default=None, type=int,
                        help='Memory per model run in MB (SLURM --mem). Default: cluster default '
                             '(~1 GB/CPU). The finest rung builds the deepest trees and may need '
                             'more; e.g. 120000 for a 128 GB node.')
    # queueManager pass-through options.
    parser.add_argument('--partition',       default=None, help='SLURM partition override')
    parser.add_argument('--jobMaximum',      default=None, type=int, help='Maximum concurrent jobs')
    parser.add_argument('--waitOnSubmit',    default=None, type=int, help='Seconds to wait after submitting')
    parser.add_argument('--waitOnActive',    default=None, type=int, help='Seconds between polling active jobs')
    args = parser.parse_args()
    for key in ('pipelinePath', 'outputDirectory'):
        val = getattr(args, key)
        if val and not val.endswith('/'):
            setattr(args, key, val + '/')
    return args


# ---------------------------------------------------------------------------
# Labels and paths
# ---------------------------------------------------------------------------

def _resolution_label(mass_resolution):
    """Filename-safe label for a mass-resolution rung, e.g. 3.16e11 -> '3p16e11'."""
    return f'{mass_resolution:.2e}'.replace('.', 'p').replace('+', '')


def _concentration_label(log_min, log_max):
    return f'c_{log_min:.1f}_{log_max:.1f}'.replace('.', 'p')


def _factor_label(factor):
    """Filename-safe label for a seeding-threshold value, e.g. 3.0 -> '3'."""
    return f'{factor:g}'.replace('.', 'p')


def _rung_prefix(output_dir, mass_resolution, factor=FACTOR_MASS_RESOLUTION_DEFAULT):
    """Prefix for one (M_res, f) run point.

    Runs at the default f carry no f suffix, so the axis-1 outputs generated
    before this second axis existed are found unchanged and are reused as
    anchor points rather than being re-run."""
    prefix = output_dir + f'convergenceMDPL2_mRes{_resolution_label(mass_resolution)}'
    if factor != FACTOR_MASS_RESOLUTION_DEFAULT:
        prefix += f'_f{_factor_label(factor)}'
    return prefix


def _model_output(prefix):
    """Return the model output file for a rung, tolerating an MPI0000 suffix."""
    for candidate in (prefix + '.hdf5', prefix + '.hdf5:MPI0000', prefix + ':MPI0000.hdf5'):
        if os.path.exists(candidate):
            return candidate
    return prefix + '.hdf5'


# ---------------------------------------------------------------------------
# Parameter-file generation
# ---------------------------------------------------------------------------

def _set_scalar(parent, tag, value):
    """Set (creating if absent) a single-valued child element `<tag value=.../>`."""
    node = parent.find(tag)
    if node is None:
        node = etree.SubElement(parent, tag)
    node.set('value', value)
    return node


def _concentration_block(log_min, log_max):
    """Build one inline `concentrationDistribution` analysis for a mass bin."""
    node = etree.Element('outputAnalysis')
    node.set('value', 'concentrationDistribution')
    _set_scalar(node, 'label',                        _concentration_label(log_min, log_max))
    _set_scalar(node, 'comment',
                f'MDPL2 self-convergence; {log_min:.1f} < log10(M/Msun) <= {log_max:.1f}; z={REDSHIFT:.1f}')
    _set_scalar(node, 'redshift',                     f'{REDSHIFT}')
    _set_scalar(node, 'massMinimum',                  f'{10.0**log_min:.5e}')
    _set_scalar(node, 'massMaximum',                  f'{10.0**log_max:.5e}')
    _set_scalar(node, 'concentrationMinimum',         f'{CONCENTRATION_MINIMUM}')
    _set_scalar(node, 'concentrationMaximum',         f'{CONCENTRATION_MAXIMUM}')
    _set_scalar(node, 'countConcentrationsPerDecade', f'{CONCENTRATION_PER_DECADE}')
    _set_scalar(node, 'timeRecent',                   f'{TIME_RECENT}')
    _set_scalar(node, 'massParticle',                 f'{MASS_PARTICLE:.5e}')
    return node


def _spin_block():
    """Build the inline `spinDistribution` analysis over the full mass range."""
    node = etree.Element('outputAnalysis')
    node.set('value', 'spinDistribution')
    _set_scalar(node, 'label',                              'spin')
    _set_scalar(node, 'comment',
                f'MDPL2 self-convergence; {MASS_TREE_MINIMUM:.2e} < M/Msun <= {MASS_TREE_MAXIMUM:.2e}; z={REDSHIFT:.1f}')
    _set_scalar(node, 'redshift',                           f'{REDSHIFT}')
    _set_scalar(node, 'massMinimum',                        f'{MASS_TREE_MINIMUM:.5e}')
    _set_scalar(node, 'massMaximum',                        f'{MASS_TREE_MAXIMUM:.5e}')
    _set_scalar(node, 'spinMinimum',                        f'{SPIN_MINIMUM:.5e}')
    _set_scalar(node, 'spinMaximum',                        f'{SPIN_MAXIMUM:.5e}')
    _set_scalar(node, 'countSpinsPerDecade',               f'{SPIN_PER_DECADE}')
    _set_scalar(node, 'timeRecent',                         f'{TIME_RECENT}')
    _set_scalar(node, 'massParticle',                       f'{MASS_PARTICLE:.5e}')
    _set_scalar(node, 'particleCountMinimum',               f'{PARTICLE_COUNT_MINIMUM:d}')
    _set_scalar(node, 'energyEstimateParticleCountMaximum', f'{ENERGY_ESTIMATE_PARTICLE_COUNT_MAX:d}')
    _set_scalar(node, 'logNormalRange',                     f'{LOGNORMAL_RANGE}')
    _set_scalar(node, 'errorTolerant',                      'true')
    return node


def _generate_parameter_file(base_tree, mass_resolution, factor, options):
    """Emit a run-point parameter file from the MDPL2 base tree; return its path."""
    root   = base_tree.getroot()
    tree   = etree.ElementTree(etree.fromstring(etree.tostring(root)))  # deep copy
    params = tree.getroot()
    prefix = _rung_prefix(options['outputDirectory'], mass_resolution, factor)

    # --- Fixed mass resolution (replace the `scaled` block). ---
    res_node = params.find('mergerTreeMassResolution')
    for child in list(res_node):
        res_node.remove(child)
    res_node.set('value', 'fixed')
    _set_scalar(res_node, 'massResolution', f'{mass_resolution:.5e}')

    # --- Seeding threshold. The base config leaves `factorMassResolution` at its
    #     default (100), so set it explicitly on the johnson2021 element, which
    #     is nested inside the `concentrationLimiter` wrapper. ---
    johnson = params.find(".//darkMatterProfileScaleRadius[@value='johnson2021']")
    if johnson is None:
        raise RuntimeError('no johnson2021 scale-radius element found in the base configuration')
    _set_scalar(johnson, 'factorMassResolution', f'{factor:.5e}')

    # --- Tree mass range and sampling. ---
    build_node = params.find('mergerTreeBuildMasses')
    _set_scalar(build_node, 'massTreeMinimum', f'{MASS_TREE_MINIMUM:.5e}')
    _set_scalar(build_node, 'massTreeMaximum', f'{MASS_TREE_MAXIMUM:.5e}')
    _set_scalar(build_node, 'treesPerDecade',  f'{options["treesPerDecade"]}')

    # --- Satellite orbits. The base config (2021) predates orbit-setting being
    #     moved into `nodeOperatorSatelliteOrbit`; the Vitvitska operator now
    #     *reads* a stored virial orbit at tree initialization and aborts
    #     ("radius has not been set") if none was set. Use the `orbiting`
    #     satellite component and prepend a `satelliteOrbit` operator that only
    #     initializes orbits (no orbital evolution, which the spin/concentration
    #     work does not need). It must precede the Vitvitska operator so the
    #     orbit exists before it is read. ---
    _set_scalar(params, 'componentSatellite', 'orbiting')
    operators = params.find("nodeOperator[@value='multi']")
    orbit_op = etree.Element('nodeOperator')
    orbit_op.set('value', 'satelliteOrbit')
    _set_scalar(orbit_op, 'acceptUnboundOrbits', 'false')
    _set_scalar(orbit_op, 'initializeOnly',      'true')
    operators.insert(0, orbit_op)

    # --- Sub-resolution angular-momentum method (Vitvitska stochastic term). If set,
    #     override the base config's choice; used to compare the resolution-convergent
    #     `resolutionScaled` method (nextStages §4.1) against the legacy `massScaled`. ---
    if options.get('subresolutionMethod'):
        vitvitska = params.find("nodeOperator[@value='multi']/nodeOperator[@value='haloAngularMomentumVitvitska2002']")
        _set_scalar(vitvitska, 'subresolutionAngularMomentumMethod', options['subresolutionMethod'])

    # --- Output file. ---
    _set_scalar(params, 'outputFileName', prefix + '.hdf5')

    # --- Raise the task's internal walltime cap so long fine-resolution runs
    #     are not aborted mid-forest. ---
    task_node = params.find('task')
    if task_node is not None and task_node.get('value') == 'evolveForests':
        _set_scalar(task_node, 'walltimeMaximum', '360000')

    # --- Replace the analyses: keep virialDensityContrastDefinition, drop the
    #     file-based analyses, add inline ones. ---
    multi = params.find("outputAnalysis[@value='multi']")
    for child in list(multi):
        if child.get('value') in ('spinDistribution', 'concentrationDistribution'):
            multi.remove(child)
    for (log_min, log_max) in MASS_BINS_LOG10:
        multi.append(_concentration_block(log_min, log_max))
    multi.append(_spin_block())

    path = prefix + '.xml'
    tree.write(path, pretty_print=True, xml_declaration=True, encoding='UTF-8')
    return path


def _generate(options):
    """Generate all run-point parameter files; return (M_res, f, path) triples."""
    os.makedirs(options['outputDirectory'], exist_ok=True)
    base_path = options['pipelinePath'] + 'spinConcentrationBaseMDPL2.xml'
    parser    = etree.XMLParser(remove_blank_text=True)
    base_tree = etree.parse(base_path, parser)
    generated = []
    for mass_resolution, factor in _run_grid(options.get('ladderOnly', False)):
        path = _generate_parameter_file(base_tree, mass_resolution, factor, options)
        generated.append((mass_resolution, factor, path))
        print(f'  generated {os.path.basename(path)}  '
              f'(M_res = {mass_resolution:.2e}, f = {factor:g}, M_seed = {factor * mass_resolution:.2e} Msun)')
    return generated


# ---------------------------------------------------------------------------
# Job construction
# ---------------------------------------------------------------------------

def _build_jobs(generated, options):
    galacticus_exe = os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe'
    jobs = []
    for mass_resolution, factor, param_path in generated:
        prefix = _rung_prefix(options['outputDirectory'], mass_resolution, factor)
        if os.path.exists(_model_output(prefix)) and options['force'] == 'no':
            continue
        job = {
            'command':    f'{galacticus_exe} {param_path}',
            'launchFile': prefix + '.sh',
            'logFile':    prefix + '.log',
            'label':      os.path.basename(prefix),
            'ppn':        options['ppn'],
            'nodes':      options['nodes'],
            'walltime':   options['walltime'],
            'mpi':        'no',
        }
        if options.get('memory') is not None:
            job['memory'] = options['memory']
        jobs.append(job)
    return jobs


# ---------------------------------------------------------------------------
# Reading model outputs
# ---------------------------------------------------------------------------

def _find_analysis_group(analyses, kind, label):
    """Find the analyses subgroup for a given analysis `kind` ('concentration' or
    'spin') and user `label`. The analyzer prefixes the stored group name with the
    analysis type (e.g. label 'c_12p5_13p0' -> 'concentrationDistributionc_12p5_13p0'),
    so match by suffix and kind rather than assuming the exact prefix."""
    for name in analyses:
        if kind in name.lower() and name.endswith(label):
            return analyses[name]
    return None


def _read_distribution(analysis_group):
    """Return (x, y) for a volumeFunction1D analysis, using the group attributes
    that name the x and y datasets."""
    x_name = analysis_group.attrs['xDataset']
    y_name = analysis_group.attrs['yDataset']
    if isinstance(x_name, bytes):
        x_name = x_name.decode()
    if isinstance(y_name, bytes):
        y_name = y_name.decode()
    return analysis_group[x_name][:], analysis_group[y_name][:]


def _weighted_median_log10(x, weights):
    """Median of a distribution sampled at points `x` with (non-negative)
    `weights`, interpolated linearly in log10(x). Returns NaN if no weight."""
    x       = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good    = np.isfinite(weights) & (weights > 0.0) & (x > 0.0)
    if good.sum() < 2:
        return np.nan
    lx  = np.log10(x[good])
    w   = weights[good]
    order = np.argsort(lx)
    lx, w = lx[order], w[order]
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return 10.0 ** np.interp(0.5, cdf, lx)


def _percentile_log10(x, weights, percentile):
    x       = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good    = np.isfinite(weights) & (weights > 0.0) & (x > 0.0)
    if good.sum() < 2:
        return np.nan
    lx  = np.log10(x[good])
    w   = weights[good]
    order = np.argsort(lx)
    lx, w = lx[order], w[order]
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return 10.0 ** np.interp(percentile, cdf, lx)


def _read_failure_rate(log_path):
    """Parse the Johnson2021 energy-model statistics from a run log.

    Returns (applied, failed) as ints, (0, 0) if the model was never applied, or
    None if the statistics line was not found."""
    if not os.path.exists(log_path):
        return None
    applied = failed = None
    never   = False
    with open(log_path, 'r', errors='replace') as handle:
        for line in handle:
            match = _FAILURE_RE.search(line)
            if match:
                failed, applied = int(match.group(1)), int(match.group(2))
            elif _NEVER_APPLIED_RE.search(line):
                never = True
    if applied is not None:
        return applied, failed
    if never:
        return 0, 0
    return None


def _collect(options):
    """Read median c per mass bin, spin median/width, and failure rate per rung."""
    results = []
    for mass_resolution, factor in _run_grid(options.get('ladderOnly', False)):
        prefix      = _rung_prefix(options['outputDirectory'], mass_resolution, factor)
        output_file = _model_output(prefix)
        record      = {'massResolution': mass_resolution,
                       'factor':         factor,
                       'massSeed':       factor * mass_resolution,
                       'concentration':  {},
                       'spin':           None}
        # A rung that aborted mid-run (e.g. OOM-killed) can leave a missing,
        # partial, or truncated/corrupt HDF5 file; treat any read failure as
        # missing data for that rung rather than aborting the whole harness.
        if not os.path.exists(output_file):
            # A run point that has not been run yet (e.g. a newly added grid
            # point, or one still queued): reported as missing, without noise.
            print(f'  (not yet run: M_res = {mass_resolution:.2e}, f = {factor:g})')
            results.append(record | {'failure': _read_failure_rate(prefix + '.log')})
            continue
        try:
            with h5py.File(output_file, 'r') as model:
                analyses = model['analyses'] if 'analyses' in model else {}
                for (log_min, log_max) in MASS_BINS_LOG10:
                    group = _find_analysis_group(analyses, 'concentration',
                                                 _concentration_label(log_min, log_max))
                    if group is not None:
                        x, y = _read_distribution(group)
                        record['concentration'][(log_min, log_max)] = _weighted_median_log10(x, y)
                spin_group = _find_analysis_group(analyses, 'spin', 'spin')
                if spin_group is not None:
                    x, y = _read_distribution(spin_group)
                    record['spin'] = {
                        'median': _weighted_median_log10(x, y),
                        'lo':     _percentile_log10(x, y, 0.16),
                        'hi':     _percentile_log10(x, y, 0.84),
                    }
        except (OSError, RuntimeError, KeyError) as error:
            print(f'  WARNING: could not read {output_file} ({error}); rung treated as missing.')
        record['failure'] = _read_failure_rate(prefix + '.log')
        results.append(record)
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _log_scale_y(ax, values):
    """Set a log y-scale only if there is positive data to scale; matplotlib
    raises otherwise, which would abort a plot of a partially-complete grid."""
    values = np.asarray([v for v in values], dtype=float)
    if np.any(np.isfinite(values) & (values > 0.0)):
        ax.set_yscale('log')


def _plot(all_results, options):
    """Axis 1: drift with M_res, one ladder per seeding threshold in `LADDER_FACTORS`.

    The f = 100 ladder is the code-default baseline; the f = 1 ladder repeats it with
    fall-back seeding minimized, so the two together separate drift caused by the moving
    seeding threshold from drift caused by tree depth acting on the energy model itself."""
    ladders = []
    for factor in LADDER_FACTORS:
        runs = sorted((r for r in all_results
                       if r['factor'] == factor and r['massResolution'] in MASS_RESOLUTION_LADDER),
                      key=lambda r: r['massResolution'])
        if runs:
            ladders.append((factor, runs, '--' if factor == FACTOR_MASS_RESOLUTION_DEFAULT else '-'))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    colors    = plt.cm.viridis(np.linspace(0.0, 0.9, len(MASS_BINS_LOG10)))

    # Panel 1: median c versus resolution, one line per mass bin, per ladder.
    ax = axes[0]
    for index, (log_min, log_max) in enumerate(MASS_BINS_LOG10):
        for (factor, runs, style) in ladders:
            resolutions = [r['massResolution'] for r in runs]
            c           = [r['concentration'].get((log_min, log_max), np.nan) for r in runs]
            ax.plot(resolutions, c, linestyle=style, marker='o', markersize=3,
                    color=colors[index],
                    label=f'{log_min:.1f}–{log_max:.1f}' if style == '-' else None)
    ax.axvline(MASS_PARTICLE, color='dimgray', linestyle=':', linewidth=1)
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'median concentration $c$')
    ax.set_title('Concentration drift', fontsize=9)
    ax.legend(title=r'$\log_{10} M$ bin', fontsize=6, title_fontsize=6)

    # Panel 2: spin median with 16-84% band versus resolution. The band is drawn only
    # for the f = 1 ladder, to keep the comparison of the two medians readable.
    ax = axes[1]
    spin_values = []
    for (factor, runs, style) in ladders:
        resolutions = [r['massResolution'] for r in runs]
        median      = np.array([r['spin']['median'] if r['spin'] else np.nan for r in runs])
        lo          = np.array([r['spin']['lo']     if r['spin'] else np.nan for r in runs])
        hi          = np.array([r['spin']['hi']     if r['spin'] else np.nan for r in runs])
        ax.plot(resolutions, median, linestyle=style, marker='o', markersize=3,
                color='orangered', label=f'$f={factor:g}$')
        if style == '-':
            ax.fill_between(resolutions, lo, hi, color='orangered', alpha=0.2)
        spin_values.append(np.concatenate([median, lo, hi]))
    ax.axvline(MASS_PARTICLE, color='dimgray', linestyle=':', linewidth=1)
    ax.set_xscale('log')
    _log_scale_y(ax, np.concatenate(spin_values) if spin_values else [])
    ax.set_xlabel(r'$M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'spin $\lambda$ (median, 16–84%)')
    ax.set_title('Spin drift', fontsize=9)
    ax.legend(fontsize=6)

    # Panel 3: energy-model failure rate versus resolution.
    ax = axes[2]
    for (factor, runs, style) in ladders:
        resolutions = [r['massResolution'] for r in runs]
        rate        = [r['failure'][1] / r['failure'][0]
                       if r['failure'] is not None and r['failure'][0] > 0 else np.nan
                       for r in runs]
        ax.plot(resolutions, rate, linestyle=style, marker='o', markersize=3,
                color='darkgreen', label=f'$f={factor:g}$')
    ax.axvline(MASS_PARTICLE, color='dimgray', linestyle=':', linewidth=1)
    ax.set_xscale('log')
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r'$M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'energy-model failure fraction')
    ax.set_title('Failure rate', fontsize=9)
    ax.legend(fontsize=6)

    fig.suptitle('MDPL2 self-convergence (Arm C), z=0', fontsize=10)
    fig.tight_layout()
    output_pdf = options['outputDirectory'] + 'convergenceMDPL2.pdf'
    fig.savefig(output_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {output_pdf}')


# Families of runs shown on the seeding-mass plot. Each is a set of runs sharing
# one held-fixed quantity, drawn as a line against M_seed:
#   - a constant-M_res family varies f alone, so tree depth is held fixed and
#     only the seeding threshold moves;
#   - the constant-f family is axis 1, where depth and threshold moved together.
# If the drift is purely a seeding-threshold effect, all families of a given
# colour lie on top of one another.
_FAMILIES = [
    ('massResolution', 1.0e11, '-',  'o', r'$M_\mathrm{res}=10^{11}$, vary $f$'),
    ('massResolution', 1.0e12, '--', 's', r'$M_\mathrm{res}=10^{12}$, vary $f$'),
    ('factor',         FACTOR_MASS_RESOLUTION_DEFAULT, ':', '^', r'$f=100$, vary $M_\mathrm{res}$'),
]


def _family_members(results, key, value):
    """Runs belonging to one family, ordered by seeding mass."""
    return sorted((r for r in results if r[key] == value), key=lambda r: r['massSeed'])


def _plot_factor(results, options):
    """Axis 2: the collapse test against seeding mass M_seed = f*M_res."""
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4))
    colors    = plt.cm.viridis(np.linspace(0.0, 0.9, len(MASS_BINS_LOG10)))

    # Panel 1: median c versus seeding mass, per mass bin, per family.
    ax = axes[0]
    for index, (log_min, log_max) in enumerate(MASS_BINS_LOG10):
        for (key, value, style, marker, _label) in _FAMILIES:
            members = _family_members(results, key, value)
            if len(members) < 2:
                continue
            seed = [r['massSeed'] for r in members]
            c    = [r['concentration'].get((log_min, log_max), np.nan) for r in members]
            ax.plot(seed, c, linestyle=style, marker=marker, markersize=3,
                    color=colors[index],
                    label=f'{log_min:.1f}–{log_max:.1f}' if style == '-' else None)
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_\mathrm{seed} = f\,M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'median concentration $c$')
    ax.set_title('Concentration vs. seeding mass', fontsize=9)
    ax.legend(title=r'$\log_{10} M$ bin', fontsize=6, title_fontsize=6)

    # Panel 2: spin versus seeding mass. `factorMassResolution` is a Johnson2021
    # parameter only -- the Vitvitska spin model has no seeding threshold -- so
    # the constant-M_res families are expected to be flat here. Any slope along
    # them would mean spin is picking up the scale-radius model, and the contrast
    # with the (steep) constant-f family isolates spin's drift as tree depth.
    ax = axes[1]
    for (key, value, style, marker, label) in _FAMILIES:
        members = _family_members(results, key, value)
        if len(members) < 2:
            continue
        seed   = [r['massSeed'] for r in members]
        median = [r['spin']['median'] if r['spin'] else np.nan for r in members]
        ax.plot(seed, median, linestyle=style, marker=marker, markersize=3,
                color='orangered', label=label)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$M_\mathrm{seed} = f\,M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'median spin $\lambda$')
    ax.set_title('Spin vs. seeding mass (null test)', fontsize=9)
    ax.legend(fontsize=6)

    # Panel 3: energy-model failure fraction. Lowering f hands the least-resolved
    # nodes to the energy model, so this is where the §5.1 positive-energy path
    # is most likely to fire.
    ax = axes[2]
    for (key, value, style, marker, label) in _FAMILIES:
        members = _family_members(results, key, value)
        if len(members) < 2:
            continue
        seed = [r['massSeed'] for r in members]
        rate = [r['failure'][1] / r['failure'][0]
                if r['failure'] is not None and r['failure'][0] > 0 else np.nan
                for r in members]
        ax.plot(seed, rate, linestyle=style, marker=marker, markersize=3,
                color='darkgreen', label=label)
    ax.set_xscale('log')
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r'$M_\mathrm{seed} = f\,M_\mathrm{res}$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'energy-model failure fraction')
    ax.set_title('Failure rate vs. seeding mass', fontsize=9)
    ax.legend(fontsize=6)

    fig.suptitle(r'MDPL2 seeding-threshold collapse test (Arm C, axis 2), $z=0$', fontsize=10)
    fig.tight_layout()
    output_pdf = options['outputDirectory'] + 'convergenceMDPL2_factor.pdf'
    fig.savefig(output_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {output_pdf}')


def _report(results):
    print('\nConvergence summary (median c per mass bin; spin median; failure rate):')
    header = (['M_res', 'f', 'M_seed']
              + [f'c[{a:.1f}-{b:.1f}]' for (a, b) in MASS_BINS_LOG10] + ['lambda', 'fail'])
    print('  ' + '  '.join(f'{h:>11}' for h in header))
    for r in sorted(results, key=lambda x: (-x['massResolution'], -x['factor'])):
        row = [f'{r["massResolution"]:.2e}', f'{r["factor"]:g}', f'{r["massSeed"]:.2e}']
        for key in MASS_BINS_LOG10:
            c = r['concentration'].get(key, np.nan)
            row.append(f'{c:.3f}' if np.isfinite(c) else '   -')
        lam = r['spin']['median'] if r['spin'] else np.nan
        row.append(f'{lam:.4f}' if np.isfinite(lam) else '   -')
        if r['failure'] is not None and r['failure'][0] > 0:
            row.append(f'{r["failure"][1] / r["failure"][0]:.2e}')
        elif r['failure'] is not None:
            row.append('n/a')
        else:
            row.append('   -')
        print('  ' + '  '.join(f'{v:>11}' for v in row))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args    = _parse_args()
    options = vars(args)

    if not options['plotOnly']:
        print('Generating rung parameter files:')
        generated = _generate(options)
        jobs      = _build_jobs(generated, options)
        if jobs:
            print(f'Submitting {len(jobs)} model run(s):')
            manager = queueManager.factory(args)
            submit_jobs(manager, jobs)
        else:
            print('All model outputs already present; nothing to submit (use --force to re-run).')

    results = _collect(options)
    _report(results)
    _plot(results, options)
    _plot_factor(results, options)


if __name__ == '__main__':
    main()
