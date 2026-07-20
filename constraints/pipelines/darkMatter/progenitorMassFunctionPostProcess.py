#!/usr/bin/env python3
# Script to generate post-processed progenitor mass function models and plots.
# Andrew Benson (05-May=2026)

import argparse
import os
import sys
import re

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import queueManager
from queueManager import translate_job, submit_jobs
from Galacticus.Constraints.Simulations import iterate, parse_simulations_xml
from Galacticus._logging                 import configure_default as _configure_default

# Show INFO-level diagnostic output from library modules.
_configure_default()


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Generate post-processed progenitor mass function models and plots.'
    )
    parser.add_argument('--pipelinePath',           required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--outputDirectory',        required=True,
                        help='Directory for generated output files')
    parser.add_argument('--force',                  default='no',
                        help='Re-run models even if outputs exist (default: no)')
    parser.add_argument('--countParticlesMinimum',  default=300, type=int,
                        help='Minimum particles per halo for constraints (default: 300)')
    parser.add_argument('--progenitorMassFunctionPPN',  default=1, type=int,
                        help='Number of processors per node to use (default: 1)')
    parser.add_argument('--select',                 default=None, action='append',
                        help='Simulation selection filter (suite::group::resolution::...); may be repeated')
    parser.add_argument('--partition',              default=None,
                        help='SLURM partition override')
    parser.add_argument('--jobMaximum',             default=None, type=int,
                        help='Maximum number of concurrent jobs')
    parser.add_argument('--waitOnSubmit',       default=None, type=int,
                        help='Time (in seconds) to wait after submitting')
    parser.add_argument('--waitOnActive',       default=None, type=int,
                        help='Time (in seconds) to wait between polling active jobs')
    args = parser.parse_args()
    for key in ('pipelinePath', 'outputDirectory'):
        val = getattr(args, key)
        if val and not val.endswith('/'):
            setattr(args, key, val + '/')
    return args


# ---------------------------------------------------------------------------
# Entry grouping over redshifts
# ---------------------------------------------------------------------------

def _group_entries(simulations, options):
    """Group entries that differ only by redshift into ordered lists.

    Returns a dict keyed by suite::group::res::sim::realization, each value
    being a list of entry dicts (one per redshift, in iterate() order).
    Each entry gets a 'countRealizations' key attached.
    """
    entry_groups = {}
    have_models  = False

    for entry in iterate(simulations, options):
        have_models = True
        group_key = '::'.join([
            entry['suite'     ]['name'],
            entry['group'     ]['name'],
            entry['resolution']['name'],
            entry['simulation']['name'],
            entry['realization'],
        ])
        entry_groups.setdefault(group_key, []).append(entry)

    if not have_models:
        raise RuntimeError('no models match this selection')

    return entry_groups

# ---------------------------------------------------------------------------
# Job construction
# ---------------------------------------------------------------------------

def _build_jobs(entry_groups, options):
    """Return list of job dicts for the two model variants."""
    galacticus_exe = os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe'
    pipeline_path  = os.environ['GALACTICUS_EXEC_PATH'] + '/constraints/pipelines/darkMatter/'
    force          = options['force']
    jobs           = []

    variants = [
        # (suffix, extra_xml, output_check_suffix)
        ('',      None                                   , ':MPI0000.hdf5'      ),
        ('_true', 'changesProgenitorMassFunctionTrue.xml', '.hdf5_true:MPI0000' ),
    ]

    for identifier in sorted(entry_groups):
        entry_group = entry_groups[identifier]
        entry       = entry_group[0]
        redshift_parent = sorted([e['redshift'] for e in entry_group],key=lambda x: float(x)).pop(0)
        identifier_ = (
            f"{entry['suite']['name']}_{entry['group']['name']}_"
            f"{entry['resolution']['name']}_{entry['simulation']['name']}_"
            f"{entry['realization']}_z{redshift_parent}"
        )
        file_prefix  = options['outputDirectory'] + f'progenitorMassFunction_{identifier_}'
        file_params  = options['outputDirectory'] + f'progenitorMassFunctionBase_{identifier_}.xml'
        file_changes = options['outputDirectory'] + f'progenitorMassFunctionChanges_{identifier_}.xml'
        label_base   = f"progenitorMassFunction_{identifier_}"
        
        for suffix, changes_xml, check_suffix in variants:
            output_file = file_prefix + check_suffix
            if os.path.exists(output_file) and force == 'no':
                continue
            command = f'{galacticus_exe} {file_params}'
            if os.path.exists(file_changes) :
                command += f' {file_changes}'
            if changes_xml is not None:
                command += f' {pipeline_path}{changes_xml}'
            ompThreads = options['progenitorMassFunctionPPN'] if "progenitorMassFunctionPPN" in options else "1"
            jobs.append({
                'command':      command,
                'launchFile':   f'{file_prefix}{suffix}.sh',
                'logFile':      f'{file_prefix}{suffix}.log',
                'label':        f'{label_base}{suffix}',
                'ompThreads':   ompThreads,
                'tasksPerNode': ompThreads,
                'nodes':        1,
                'walltime':     '2:00:00',
                'mpi':          'no',
            })

    return jobs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _read_pmf(hdf5_path, analysis_name):
    """Read massRatio, progenitorMassFunction, and  from analyses."""
    with h5py.File(hdf5_path, 'r') as f:
        grp            = f['analyses/' + analysis_name]
        mass_ratio     = grp['massRatio'][:]
        pmf            = grp['progenitorMassFunction'][:]
        pmf_target     = grp['progenitorMassFunctionTarget'][:]
        pmf_cov        = grp['progenitorMassFunctionCovariance'][:]
        pmf_cov_target = grp['progenitorMassFunctionCovarianceTarget'][:]
    return mass_ratio, pmf, pmf_cov, pmf_target, pmf_cov_target


def _plot_group(group, options):
    entry           = group[0]
    redshift_progenitors = sorted([e['redshift'] for e in group],key=lambda x: float(x))
    redshift_parent = redshift_progenitors.pop(0)

    identifier_ = (
        f"{entry['suite']['name']}_{entry['group']['name']}_"
        f"{entry['resolution']['name']}_{entry['simulation']['name']}_"
        f"{entry['realization']}"
    )
    file_prefix = options['outputDirectory'] + f'progenitorMassFunction_{identifier_}_z{redshift_parent}'
    model = h5py.File(file_prefix + ':MPI0000.hdf5','r')
    analyses = model['analyses']

    for name, analysis in analyses.items():
        if not isinstance(analysis, h5py.Group) or not re.match(r'progenitorMassFunction',name):
            continue
        redshift_progenitor =          analysis.attrs['redshiftProgenitor'        ]
        mass_ratio_min      =          analysis.attrs['massRatioLikelihoodMinimum']
        mass_ratio_max      =          analysis.attrs['massRatioLikelihoodMaximum']
        mass_parent_min     = np.log10(analysis.attrs['massParentMinimum'         ])
        mass_parent_max     = np.log10(analysis.attrs['massParentMaximum'         ])
        output_pdf          = file_prefix + f'_z{redshift_progenitor:.3f}_logM{mass_parent_min:.1f}-{mass_parent_max:.1f}.pdf'

        if os.path.exists(output_pdf) and options['force'] == 'no':
            return

        # Read model outputs.
        mass_ratio_model, pmf_model, pmf_cov_model, pmf_target, pmf_cov_target = _read_pmf(file_prefix + ':MPI0000.hdf5'     , name)
        mass_ratio_true,  pmf_true , pmf_cov_true, *_                          = _read_pmf(file_prefix + '.hdf5_true:MPI0000', name)
      
        # Determine axis ranges (decade-rounded).
        non_zero = np.where((pmf_target > 0) & (mass_ratio_model >= mass_ratio_min) & (mass_ratio_model <= mass_ratio_max))
        x_min = 10 ** np.floor(np.log10(mass_ratio_model[non_zero].min()))
        x_max = 10 ** np.ceil( np.log10(mass_ratio_model[non_zero].max()))
        y_min = 10 ** np.floor(np.log10(pmf_target[non_zero].min()))
        y_max = 10 ** np.ceil( np.log10(pmf_target[non_zero].max()))

        # Build the figure.
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(r'$M_\mathrm{prog}/M_\mathrm{par}$')
        ax.set_ylabel(r'$\mathrm{d}f / \mathrm{d}\log M_\mathrm{prog}$')
        ax.set_title(
            f"{entry['suite']['name']} {entry['group']['name']} "
            f"{entry['resolution']['name']} {entry['simulation']['name']} "
            f"{entry['realization']} $\\log_{{10}}M_\\mathrm{{par}}={mass_parent_min:.1f}$–${mass_parent_max:.1f}$ $z={entry['redshift']}$ → ${redshift_progenitor:.3f}$",
            fontsize=8
        )

        # 1. Vertical lines at min/max mass ratios.
        ax.axvline(mass_ratio_min, color='dimgray', linestyle='--', linewidth=1)
        ax.axvline(mass_ratio_max, color='dimgray', linestyle='--', linewidth=1)

        # 2. Target data outside fit range (transparent).
        non_fit = (mass_ratio_model < mass_ratio_min) | (mass_ratio_model > mass_ratio_max)
        if non_fit.size > 0:
            ax.errorbar(
                mass_ratio_model[non_fit], pmf_target[non_fit],
                yerr=np.sqrt(np.diag(pmf_cov_target))[non_fit],
                fmt='o', markersize=2, color='mediumseagreen', alpha=0.25,
                elinewidth=0.5, capsize=0
            )

        # 3. Target data in fit region.
        fit = (mass_ratio_model >= mass_ratio_min) & (mass_ratio_model <= mass_ratio_max)
        if fit.size > 0:
            ax.errorbar(
                mass_ratio_model[fit], pmf_target[fit],
                yerr=np.sqrt(np.diag(pmf_cov_target))[fit],
                fmt='o', markersize=2, color='mediumseagreen',
                elinewidth=0.5, capsize=0
            )

        # 4. True model (no numerical artifacts).
        nz = pmf_true > 0
        if nz.any():
            ax.plot(mass_ratio_true[nz], pmf_true[nz],
                    color='orangered', linestyle='--', linewidth=1)

        # 5. Base (fitted) model.
        nz = pmf_model > 0
        if nz.any():
            pmf_err = np.sqrt(np.diag(pmf_cov_model))
            ax.fill_between(mass_ratio_model[nz], pmf_model[nz] - pmf_err[nz], pmf_model[nz] + pmf_err[nz],
                            color='orangered', alpha=0.25)
            ax.plot(mass_ratio_model[nz], pmf_model[nz],
                    color='orangered', linestyle='-', linewidth=1)

        fig.savefig(output_pdf, bbox_inches='tight')
        plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args    = _parse_args()
    options = vars(args)

    simulations = parse_simulations_xml(options['pipelinePath'] + 'simulations.xml')

    # Force re-run when --select is given.
    if options.get('select') is not None:
        options['force'] = 'yes'

    # Group entries over redshift.
    entry_groups = _group_entries(simulations, options)

    # Submit Galacticus jobs for all model variants.
    jobs = _build_jobs(entry_groups, options)
    if jobs:
        manager = queueManager.factory(args)
        submit_jobs(manager, jobs)

    # Generate comparison plots.
    for identifier in sorted(entry_groups):
        group = entry_groups[identifier]
        _plot_group(group, options)


if __name__ == '__main__':
    main()
