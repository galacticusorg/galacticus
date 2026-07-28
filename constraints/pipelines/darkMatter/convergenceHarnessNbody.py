#!/usr/bin/env python3
"""Converge-to-truth harness for the dark-matter structure models (Arms A & B).

Arms A and B of the harness described in `darkMatterPipelineReview.md` §10. Unlike
Arm C (`convergenceHarness.py`, MDPL2 self-convergence), these arms compare the
current spin/concentration model against *real N-body targets* and ask whether the
model converges to that truth as the N-body resolution improves:

  Arm A  Symphony Milky-Way CDM, resolutions X1 / X8 / X64 (realization Halo004).
  Arm B  COZMIC Milky-Way non-CDM, a representative spanning set of power spectra,
         resolutions X1 / X8.

Both suites are zoom-in simulations of a single Milky-Way host; the concentration
and spin targets are measured for the *field* halos in the surrounding zoom region
(not bound subhalos), so no orbital/tidal evolution is needed and the model setup is
Arm-C style (isolated halos). What differs from Arm C:

  * the Symphony/COZMIC cosmology and power spectrum replace MDPL's;
  * the resolution is set by the actual simulation particle mass at each X-level,
    so the model tracks the N-body's own resolution ladder;
  * for COZMIC each model swaps in a tabulated non-CDM transfer function;
  * the analyses read *real* N-body target distributions (file-based path) rather
    than being configured inline.

The §2.2-step-2 target-format reconciliation is implemented here: the Symphony and
COZMIC targets are single 2-D files (distribution x mass), whereas the file-based
`concentrationDistribution` / `spinDistribution` analyses expect the MDPL per-mass-
bin 1-D format (`_reformat_target`, below). Each 2-D target is sliced into one 1-D
per-mass-bin file carrying the attributes the analysis constructor reads.

Andrew Benson (2026).
"""

import argparse
import os

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from lxml import etree

import queueManager
from queueManager import submit_jobs

# Reuse the distribution-reading and weighted-median helpers from the Arm C harness.
from convergenceHarness import (_set_scalar, _read_distribution, _find_analysis_group,
                                _weighted_median_log10, _percentile_log10, _read_failure_rate,
                                _model_output)

# ---------------------------------------------------------------------------
# Harness configuration
# ---------------------------------------------------------------------------

DATASETS = '/scratch/abenson/datasets/static/darkMatter'
SYMPHONY_DIR = DATASETS + '/Symphony/MilkyWay'
COZMIC_DIR   = DATASETS + '/COZMIC/MilkyWay'

# Milky-Way zoom particle mass (M_sun) at each resolution level. Symphony and COZMIC
# share the same MW zoom initial conditions, so these are identical (see
# `simulation_Symphony_MilkyWay.xml` / `simulation_COZMIC_MilkyWay.xml`).
PARTICLE_MASS = {
    'resolutionX1' : 4.02830e5,
    'resolutionX8' : 5.03538e4,
    'resolutionX64': 6.29448e3,
}

# The single realization for which the higher-resolution runs exist.
REALIZATION = 'Halo004'

# The model's fixed mass resolution is tied to the simulation particle mass:
#     M_res = MODEL_RES_PARTICLES * massParticle(resolution).
# This makes the model's tree resolution track the N-body's, so X1/X8/X64 give three
# model resolutions matching the three N-body resolutions -- the convergence axis.
# 10 particles is a coarse floor mirroring the N-body's own halo-detection threshold;
# it is a documented knob worth revisiting (a smaller value resolves more history at
# fixed cost only for the low-mass trees).
MODEL_RES_PARTICLES = 10.0

# Symphony/COZMIC (Mao et al. 2015 MW zoom ICs; Nadler et al. 2023) cosmology and
# power spectrum. NB: the *target HDF5* files carry a mislabeled Planck cosmology in
# their metadata; the values below are the true simulation cosmology and match the
# pipeline's `cosmology_Symphony.xml` / `powerSpectrum_Symphony.xml`.
COSMOLOGY = {
    'HubbleConstant' : '70.000000',
    'OmegaMatter'    : ' 0.286000',
    'OmegaDarkEnergy': ' 0.714000',
    'OmegaBaryon'    : ' 0.047000',
    'temperatureCMB' : ' 2.725480',
}
SIGMA_8            = '0.820'
POWER_SPECTRUM_INDEX = '0.96'

# Arm B: a representative spanning set of COZMIC non-CDM models, one per physics class
# and cutoff strength, each with both X1 and X8 for Halo004. Each entry maps the model
# name (its target-directory name) to its tabulated transfer-function parameter file.
COZMIC_MODELS = {
    'WDM:3keV'         : 'transferFunction_COZMIC_WDM:3keV.xml',           # strong thermal cutoff
    'WDM:6keV'         : 'transferFunction_COZMIC_WDM:6keV.xml',           # mild thermal cutoff
    'FDM:25.9e-22eV'   : 'transferFunction_COZMIC_FDM:25.9e-22eV.xml',     # fuzzy dark matter
    'IDM:1GeV:halfmode': 'transferFunction_COZMIC_IDM:1GeV:halfmode.xml',  # interacting dark matter
    'WDM:3keV:bump'    : 'transferFunction_COZMIC_WDM:3keV:bump.xml',      # cutoff + power bump
}

# Only mass bins with at least this many N-body halos are turned into targets; below
# this the measured distribution is too noisy to be a meaningful convergence target.
COUNT_MINIMUM = 100

# Spin N-body error-model settings (as in Arm C); the per-bin particle-count minimum
# is set from the bin mass and the simulation resolution, so it is resolution-aware.
ENERGY_ESTIMATE_PARTICLE_COUNT_MAX = 1000
LOGNORMAL_RANGE = 1.420

REDSHIFT = 0.0
TIME_RECENT = 0.0
HALF_DEX = 10.0 ** 0.25   # half of a 0.5-dex mass bin, for bin edges from centers


# ---------------------------------------------------------------------------
# Run grid
# ---------------------------------------------------------------------------

def _run_grid(options):
    """The list of runs as dicts. Each run is one (suite, model, resolution) point."""
    runs = []
    # Arm A: Symphony CDM at three resolutions.
    for resolution in ('resolutionX1', 'resolutionX8', 'resolutionX64'):
        runs.append({
            'arm': 'A', 'suite': 'Symphony', 'model': 'CDM', 'resolution': resolution,
            'targetDir': f'{SYMPHONY_DIR}/{resolution}/CDM/{REALIZATION}',
            'transferFile': None,   # CDM: CAMB, set in the base config
            'label': f'Symphony_CDM_{resolution[len("resolution"):]}',
        })
    # Arm B: COZMIC non-CDM at two resolutions.
    if not options['armAOnly']:
        for model, transfer in COZMIC_MODELS.items():
            for resolution in ('resolutionX1', 'resolutionX8'):
                safe = model.replace(':', '-')
                runs.append({
                    'arm': 'B', 'suite': 'COZMIC', 'model': model, 'resolution': resolution,
                    'targetDir': f'{COZMIC_DIR}/{resolution}/{model}/{REALIZATION}',
                    'transferFile': transfer,
                    'label': f'COZMIC_{safe}_{resolution[len("resolution"):]}',
                })
    if options['armAOnly']:
        return [r for r in runs if r['arm'] == 'A']
    return runs


# ---------------------------------------------------------------------------
# Target reformatting: 2-D N-body target -> per-mass-bin 1-D MDPL-format files
# ---------------------------------------------------------------------------

def _reformat_target(source_file, kind, mass_particle, out_dir, run_label):
    """Slice a 2-D Symphony/COZMIC target into per-mass-bin 1-D files.

    `kind` is 'concentration' or 'spin'. Returns a list of
    (massMinimum, massMaximum, path) for the bins with sufficient counts, in
    increasing mass order. Each output file mirrors the MDPL per-bin format the
    file-based analysis constructor reads (`concentration_distribution/_class.F90`,
    `spin_distribution/_class.F90`)."""
    axis_dataset = 'concentration' if kind == 'concentration' else 'spin'
    dist_dataset = f'{kind}DistributionFunction'
    bins = []
    with h5py.File(source_file, 'r') as source:
        group = source['simulation0001']
        axis  = group[axis_dataset][:]                 # (nAxis,)  concentration or spin centers
        dist  = group[dist_dataset][:]                 # (nAxis, nMass)
        count = group['count'][:]                      # (nAxis, nMass)
        mass  = group['mass'][:]                        # (nMass,)  bin centers, 0.5 dex
    os.makedirs(out_dir, exist_ok=True)
    for j in range(mass.size):
        total = int(count[:, j].sum())
        if total < COUNT_MINIMUM:
            continue
        mass_minimum = mass[j] / HALF_DEX
        mass_maximum = mass[j] * HALF_DEX
        path = f'{out_dir}/{run_label}_{kind}_bin{j:02d}.hdf5'
        with h5py.File(path, 'w') as target:
            sim = target.create_group('simulation0001')
            sim.create_dataset(axis_dataset, data=axis)
            sim.create_dataset(dist_dataset, data=dist[:, j])
            sim.create_dataset('count',      data=count[:, j])
            attributes = sim.create_group('simulation')
            attributes.attrs['labelTarget']  = f'{run_label} {mass_minimum:.2e}-{mass_maximum:.2e}'
            attributes.attrs['massMinimum']  = float(mass_minimum)
            attributes.attrs['massMaximum']  = float(mass_maximum)
            attributes.attrs['massParticle'] = float(mass_particle)
            attributes.attrs['redshift']     = float(REDSHIFT)
            attributes.attrs['timeRecent']   = float(TIME_RECENT)
            if kind == 'spin':
                # Resolution-aware minimum particle count: the number of particles in a
                # halo at the low edge of this mass bin, floored at 1.
                particle_count_minimum = max(1, int(round(mass_minimum / mass_particle)))
                # Types must match the Fortran constructor: `particleCountMinimum` is a
                # default (32-bit) integer; `energyEstimateParticleCountMaximum` is read
                # as a double (`spin_distribution/_class.F90:371-373`).
                attributes.attrs['particleCountMinimum']               = np.int32(particle_count_minimum)
                attributes.attrs['energyEstimateParticleCountMaximum'] = float(ENERGY_ESTIMATE_PARTICLE_COUNT_MAX)
        bins.append((float(mass_minimum), float(mass_maximum), path))
    return bins


# ---------------------------------------------------------------------------
# Parameter-file generation
# ---------------------------------------------------------------------------

def _apply_cosmology(params):
    """Swap the MDPL cosmology/power-spectrum blocks for the Symphony/COZMIC ones."""
    cosmology = params.find('cosmologyParameters')
    for key, value in COSMOLOGY.items():
        _set_scalar(cosmology, key, value)
    _set_scalar(params.find('powerSpectrumPrimordial'), 'index', POWER_SPECTRUM_INDEX)
    _set_scalar(params.find('cosmologicalMassVariance'), 'sigma_8', SIGMA_8)


def _apply_transfer_function(params, transfer_file, pipeline_path):
    """CDM keeps the base config's CAMB transfer; non-CDM reads a tabulated file."""
    if transfer_file is None:
        return
    transfer_tree = etree.parse(pipeline_path + transfer_file,
                                etree.XMLParser(remove_blank_text=True))
    new_transfer = transfer_tree.find('transferFunction')
    old_transfer = params.find('transferFunction')
    # The tabulated file paths in the transfer-function XML are pipeline-relative; make
    # them absolute against the execution path so the run finds them from any cwd.
    file_name = new_transfer.find('fileName')
    if file_name is not None and not file_name.get('value').startswith('/'):
        exec_path = os.environ['GALACTICUS_EXEC_PATH'].rstrip('/')
        file_name.set('value', f"{exec_path}/{file_name.get('value')}")
    old_transfer.getparent().replace(old_transfer, new_transfer)


def _apply_orbit_fix(params):
    """Arm-C-style isolated-halo orbit setup (initialize only, no orbital evolution)."""
    _set_scalar(params, 'componentSatellite', 'orbiting')
    operators = params.find("nodeOperator[@value='multi']")
    orbit_op = etree.Element('nodeOperator')
    orbit_op.set('value', 'satelliteOrbit')
    _set_scalar(orbit_op, 'acceptUnboundOrbits', 'false')
    _set_scalar(orbit_op, 'initializeOnly',      'true')
    operators.insert(0, orbit_op)


def _concentration_analysis(mass_min, mass_max, path):
    node = etree.Element('outputAnalysis')
    node.set('value', 'concentrationDistribution')
    _set_scalar(node, 'fileName', path)
    _set_scalar(node, 'label',    f'c_{np.log10(mass_min):.2f}'.replace('.', 'p'))
    _set_scalar(node, 'comment',  f'{mass_min:.3e} < M/Msun <= {mass_max:.3e}; z={REDSHIFT}')
    return node


def _spin_analysis(mass_min, mass_max, path):
    node = etree.Element('outputAnalysis')
    node.set('value', 'spinDistribution')
    _set_scalar(node, 'fileName',       path)
    _set_scalar(node, 'label',          f's_{np.log10(mass_min):.2f}'.replace('.', 'p'))
    _set_scalar(node, 'comment',        f'{mass_min:.3e} < M/Msun <= {mass_max:.3e}; z={REDSHIFT}')
    _set_scalar(node, 'logNormalRange', f'{LOGNORMAL_RANGE}')
    _set_scalar(node, 'errorTolerant',  'true')
    return node


def _generate_parameter_file(base_tree, run, options):
    """Emit a parameter file for one (suite, model, resolution) run; return its path."""
    tree   = etree.ElementTree(etree.fromstring(etree.tostring(base_tree.getroot())))
    params = tree.getroot()
    prefix = options['outputDirectory'] + run['label']
    mass_particle = PARTICLE_MASS[run['resolution']]

    _apply_cosmology(params)
    _apply_transfer_function(params, run['transferFile'], options['pipelinePath'])
    _apply_orbit_fix(params)

    # Fixed mass resolution tied to the simulation particle mass.
    res_node = params.find('mergerTreeMassResolution')
    for child in list(res_node):
        res_node.remove(child)
    res_node.set('value', 'fixed')
    _set_scalar(res_node, 'massResolution', f'{MODEL_RES_PARTICLES * mass_particle:.5e}')

    # N-body error model uses this simulation's particle mass.
    _set_scalar(params.find('nbodyHaloMassError'), 'massParticle', f'{mass_particle:.5e}')

    # Reformat this run's 2-D targets into per-mass-bin 1-D files.
    target_dir = options['outputDirectory'] + 'targets/' + run['label']
    conc_bins = _reformat_target(f"{run['targetDir']}/concentrationDistributionFunction_z0.000.hdf5",
                                 'concentration', mass_particle, target_dir, run['label'])
    spin_bins = _reformat_target(f"{run['targetDir']}/spinDistributionFunction_z0.000.hdf5",
                                 'spin', mass_particle, target_dir, run['label'])

    # Tree mass range spans the union of the populated target bins.
    all_edges = [e for (lo, hi, _p) in conc_bins + spin_bins for e in (lo, hi)]
    _set_scalar(params.find('mergerTreeBuildMasses'), 'massTreeMinimum', f'{min(all_edges):.5e}')
    _set_scalar(params.find('mergerTreeBuildMasses'), 'massTreeMaximum', f'{max(all_edges):.5e}')
    _set_scalar(params.find('mergerTreeBuildMasses'), 'treesPerDecade',  f'{options["treesPerDecade"]}')

    _set_scalar(params, 'outputFileName', prefix + '.hdf5')
    task_node = params.find('task')
    if task_node is not None and task_node.get('value') == 'evolveForests':
        _set_scalar(task_node, 'walltimeMaximum', '360000')

    # Replace the file-based MDPL analyses with the reformatted per-bin ones.
    multi = params.find("outputAnalysis[@value='multi']")
    for child in list(multi):
        if child.get('value') in ('spinDistribution', 'concentrationDistribution'):
            multi.remove(child)
    for (lo, hi, path) in conc_bins:
        multi.append(_concentration_analysis(lo, hi, path))
    for (lo, hi, path) in spin_bins:
        multi.append(_spin_analysis(lo, hi, path))

    path = prefix + '.xml'
    tree.write(path, pretty_print=True, xml_declaration=True, encoding='UTF-8')
    run['concentrationBins'] = conc_bins
    run['spinBins']          = spin_bins
    return path


def _generate(options):
    os.makedirs(options['outputDirectory'], exist_ok=True)
    base_tree = etree.parse(options['pipelinePath'] + 'spinConcentrationBaseMDPL2.xml',
                            etree.XMLParser(remove_blank_text=True))
    runs = _run_grid(options)
    for run in runs:
        run['paramPath'] = _generate_parameter_file(base_tree, run, options)
        print(f"  generated {os.path.basename(run['paramPath'])}  "
              f"({len(run['concentrationBins'])} c-bins, {len(run['spinBins'])} spin-bins, "
              f"M_res = {MODEL_RES_PARTICLES * PARTICLE_MASS[run['resolution']]:.2e})")
    return runs


# ---------------------------------------------------------------------------
# Job construction
# ---------------------------------------------------------------------------

def _build_jobs(runs, options):
    galacticus_exe = os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe'
    jobs = []
    for run in runs:
        prefix = options['outputDirectory'] + run['label']
        if os.path.exists(_model_output(prefix)) and options['force'] == 'no':
            continue
        job = {
            'command':    f"{galacticus_exe} {run['paramPath']}",
            'launchFile': prefix + '.sh',
            'logFile':    prefix + '.log',
            'label':      run['label'],
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
# Collection
# ---------------------------------------------------------------------------

def _target_median(path, kind):
    """Median of a reformatted 1-D target file (model-independent truth)."""
    axis_dataset = 'concentration' if kind == 'concentration' else 'spin'
    dist_dataset = f'{kind}DistributionFunction'
    with h5py.File(path, 'r') as target:
        group = target['simulation0001']
        return _weighted_median_log10(group[axis_dataset][:], group[dist_dataset][:])


def _collect(runs, options):
    """For each run and mass bin, read model and target median c and lambda."""
    for run in runs:
        prefix      = options['outputDirectory'] + run['label']
        output_file = _model_output(prefix)
        run['concentration'] = {}   # massMinimum -> {'model':.., 'truth':..}
        run['spin']          = {}
        run['failure']       = _read_failure_rate(prefix + '.log')
        if not os.path.exists(output_file):
            print(f"  (not yet run: {run['label']})")
            continue
        try:
            with h5py.File(output_file, 'r') as model:
                analyses = model['analyses'] if 'analyses' in model else {}
                for (lo, hi, path) in run['concentrationBins']:
                    label = f'c_{np.log10(lo):.2f}'.replace('.', 'p')
                    group = _find_analysis_group(analyses, 'concentration', label)
                    model_median = (_weighted_median_log10(*_read_distribution(group))
                                    if group is not None else np.nan)
                    run['concentration'][lo] = {
                        'model': model_median,
                        'truth': _target_median(path, 'concentration'),
                        'mass':  np.sqrt(lo * hi),
                    }
                for (lo, hi, path) in run['spinBins']:
                    label = f's_{np.log10(lo):.2f}'.replace('.', 'p')
                    group = _find_analysis_group(analyses, 'spin', label)
                    model_median = (_weighted_median_log10(*_read_distribution(group))
                                    if group is not None else np.nan)
                    run['spin'][lo] = {
                        'model': model_median,
                        'truth': _target_median(path, 'spin'),
                        'mass':  np.sqrt(lo * hi),
                    }
        except (OSError, RuntimeError, KeyError) as error:
            print(f'  WARNING: could not read {output_file} ({error}); run treated as missing.')
    return runs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

_RES_STYLE = {'resolutionX1': (':', 'o'), 'resolutionX8': ('--', 's'), 'resolutionX64': ('-', '^')}


def _plot_suite(runs, suite_runs, title, output_pdf):
    """Model vs. truth median c(M) and lambda(M), one line per resolution.

    Solid/dashed/dotted encodes resolution; model is coloured, truth is gray. The
    convergence question is whether the coloured (model) lines approach the gray
    (truth) lines as resolution goes X1 -> X64, and whether both flatten."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for (kind, key, ax) in (('concentration', 'concentration', axes[0]),
                            ('spin',          'spin',          axes[1])):
        for run in suite_runs:
            style, marker = _RES_STYLE[run['resolution']]
            data = sorted(run[key].values(), key=lambda d: d['mass'])
            if not data:
                continue
            mass  = [d['mass']  for d in data]
            model = [d['model'] for d in data]
            truth = [d['truth'] for d in data]
            res = run['resolution'][len('resolution'):]
            ax.plot(mass, model, style, marker=marker, markersize=3, color='C0',
                    label=f'model {res}')
            ax.plot(mass, truth, style, marker=marker, markersize=3, color='gray',
                    label=f'truth {res}')
        ax.set_xscale('log')
        ax.set_xlabel(r'$M$ [$\mathrm{M}_\odot$]')
        if kind == 'spin':
            ax.set_yscale('log')
            ax.set_ylabel(r'median spin $\lambda$')
        else:
            ax.set_ylabel(r'median concentration $c$')
        ax.set_title(kind, fontsize=9)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(output_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {output_pdf}')


def _plot(runs, options):
    arm_a = [r for r in runs if r['arm'] == 'A']
    if arm_a:
        _plot_suite(runs, arm_a, r'Arm A: Symphony MW CDM, model vs. N-body truth, $z=0$',
                    options['outputDirectory'] + 'convergenceSymphony.pdf')
    for model in COZMIC_MODELS:
        arm_b = [r for r in runs if r['arm'] == 'B' and r['model'] == model]
        if arm_b:
            safe = model.replace(':', '-')
            _plot_suite(runs, arm_b, f'Arm B: COZMIC {model}, model vs. N-body truth, $z=0$',
                        options['outputDirectory'] + f'convergenceCOZMIC_{safe}.pdf')


def _report(runs):
    print('\nConverge-to-truth summary (median c and lambda; model | truth per mass bin):')
    for run in runs:
        if not run['concentration'] and not run['spin']:
            continue
        print(f"\n  {run['label']}  (failure={run['failure']})")
        for lo in sorted(run['concentration']):
            d = run['concentration'][lo]
            print(f"    c  M={d['mass']:.2e}  model={d['model']:.3f}  truth={d['truth']:.3f}")
        for lo in sorted(run['spin']):
            d = run['spin'][lo]
            print(f"    la M={d['mass']:.2e}  model={d['model']:.4f}  truth={d['truth']:.4f}")


# ---------------------------------------------------------------------------
# Argument parsing and main
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Converge-to-truth harness (Arms A & B: Symphony/COZMIC vs N-body).')
    parser.add_argument('--pipelinePath',    required=True)
    parser.add_argument('--outputDirectory', required=True)
    parser.add_argument('--force',           default='no')
    parser.add_argument('--plotOnly',        action='store_true')
    parser.add_argument('--armAOnly',        action='store_true',
                        help='Run only Arm A (Symphony); skip the COZMIC non-CDM arm.')
    parser.add_argument('--treesPerDecade',  default=300, type=int)
    parser.add_argument('--ppn',             default=16, type=int)
    parser.add_argument('--nodes',           default=1, type=int)
    parser.add_argument('--walltime',        default='24:00:00')
    parser.add_argument('--memory',          default=None, type=int)
    parser.add_argument('--partition',       default=None)
    parser.add_argument('--jobMaximum',      default=None, type=int)
    parser.add_argument('--waitOnSubmit',    default=None, type=int)
    parser.add_argument('--waitOnActive',    default=None, type=int)
    args = parser.parse_args()
    for key in ('pipelinePath', 'outputDirectory'):
        value = getattr(args, key)
        if value and not value.endswith('/'):
            setattr(args, key, value + '/')
    return args


def main():
    args    = _parse_args()
    options = vars(args)
    if not options['plotOnly']:
        print('Generating run parameter files:')
        runs = _generate(options)
        jobs = _build_jobs(runs, options)
        if jobs:
            print(f'Submitting {len(jobs)} model run(s):')
            submit_jobs(queueManager.factory(args), jobs)
        else:
            print('All model outputs already present; nothing to submit (use --force to re-run).')
    else:
        runs = _generate(options)   # regenerate to populate bin metadata for collection
    _collect(runs, options)
    _report(runs)
    _plot(runs, options)


if __name__ == '__main__':
    main()
