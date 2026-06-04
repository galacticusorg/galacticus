#!/usr/bin/env python3
"""Script to generate post-processed halo mass function models and plots.

Python port of constraints/pipelines/darkMatter/haloMassFunctionPostProcess.pl
Andrew Benson (ported to Python 2026)
"""

import argparse
import os

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
        description='Generate post-processed halo mass function models and plots.'
    )
    parser.add_argument('--pipelinePath',           required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--outputDirectory',        required=True,
                        help='Directory for generated output files')
    parser.add_argument('--force',                  default='no',
                        help='Re-run models even if outputs exist (default: no)')
    parser.add_argument('--countParticlesMinimum',  default=300, type=int,
                        help='Minimum particles per halo for constraints (default: 300)')
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
# Job construction
# ---------------------------------------------------------------------------

def _build_jobs(entries, options):
    """Return list of Perl-style job dicts for the four model variants."""
    galacticus_exe = os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe'
    pipeline_path  = os.environ['GALACTICUS_EXEC_PATH'] + '/constraints/pipelines/darkMatter/'
    force          = options['force']
    jobs           = []

    variants = [
        # (suffix, extra_xml, output_check_suffix)
        ('',           None,                    ':MPI0000.hdf5'           ),
        ('_true',      'changesTrue.xml',        '.hdf5_true:MPI0000'      ),
        ('_despali2015', 'changesDespali2015.xml', '.hdf5_Despali2015:MPI0000'),
        ('_bohr2021',  'changesBohr2021.xml',    '.hdf5_Bohr2021:MPI0000'  ),
    ]

    for entry in entries:
        file_prefix = entry['filePrefix']
        file_params = entry['fileParameters']
        label_base  = (
            f"haloMassFunction_{entry['suite']['name']}_{entry['group']['name']}_"
            f"{entry['resolution']['name']}_{entry['simulation']['name']}_"
            f"{entry['realization']}_z{entry['redshift']}"
        )

        for suffix, changes_xml, check_suffix in variants:
            output_file = file_prefix + check_suffix
            if os.path.exists(output_file) and force == 'no':
                continue
            command = f'{galacticus_exe} {file_params}'
            if changes_xml is not None:
                command += f' {pipeline_path}{changes_xml}'
            jobs.append({
                'command':    command,
                'launchFile': f'{file_prefix}{suffix}.sh',
                'logFile':    f'{file_prefix}{suffix}.log',
                'label':      f'{label_base}{suffix}',
                'ppn':        1,
                'nodes':      1,
                'walltime':   '2:00:00',
                'mpi':        'no',
            })

    return jobs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _read_hmf(hdf5_path):
    """Read haloMass and haloMassFunctionLnMBinAveraged from Outputs/Output1."""
    with h5py.File(hdf5_path, 'r') as f:
        grp = f['Outputs/Output1']
        mass = grp['haloMass'][:]
        hmf  = grp['haloMassFunctionLnMBinAveraged'][:]
    return mass, hmf


def _plot_entry(entry, options):
    file_prefix  = entry['filePrefix']
    output_pdf   = file_prefix + '.pdf'
    if os.path.exists(output_pdf) and options['force'] == 'no':
        return

    # Read model outputs.
    mass_model,       hmf_model       = _read_hmf(file_prefix + ':MPI0000.hdf5')
    mass_true,        hmf_true        = _read_hmf(file_prefix + '.hdf5_true:MPI0000')
    mass_despali2015, hmf_despali2015 = _read_hmf(file_prefix + '.hdf5_Despali2015:MPI0000')
    mass_bohr2021,    hmf_bohr2021    = _read_hmf(file_prefix + '.hdf5_Bohr2021:MPI0000')

    # Read target data.
    with h5py.File(entry['fileTargetData'], 'r') as f:
        grp          = f['simulation0001']
        mass_target  = grp['mass'][:]
        hmf_target   = grp['massFunction'][:]
        count_target = grp['count'][:]

    # Derive minimum halo mass threshold and filter masks.
    mass_halo_min   = options['countParticlesMinimum'] * entry['resolution']['massParticle']
    hmf_err         = hmf_target.copy()
    non_zero        = np.where((count_target > 0) & (mass_target >= mass_halo_min))[0]
    non_fit         = np.where((count_target > 0) & (mass_target <  mass_halo_min))[0]
    hmf_err[non_zero] /= np.sqrt(count_target[non_zero].astype(float))

    # Determine axis ranges (decade-rounded, matching Perl floor/ceil logic).
    x_min = 10 ** np.floor(np.log10(mass_target[non_zero].min()))
    x_max = 10 ** np.ceil( np.log10(mass_target[non_zero].max()))
    y_min = 10 ** np.floor(np.log10(hmf_target[non_zero].min()))
    y_max = 10 ** np.ceil( np.log10(hmf_target[non_zero].max()))

    # Build the figure.
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(r'$M$ [$\mathrm{M}_\odot$]')
    ax.set_ylabel(r'$\mathrm{d}n / \mathrm{d}\log M$ [Mpc$^{-3}$]')
    ax.set_title(
        f"{entry['suite']['name']} {entry['group']['name']} "
        f"{entry['resolution']['name']} {entry['simulation']['name']} "
        f"{entry['realization']} $z={entry['redshift']}$",
        fontsize=8
    )

    # 1. Vertical line at minimum halo mass.
    ax.axvline(mass_halo_min, color='dimgray', linestyle='--', linewidth=1)

    # 2. Target data below fit threshold (transparent).
    if non_fit.size > 0:
        ax.errorbar(
            mass_target[non_fit], hmf_target[non_fit],
            yerr=hmf_err[non_fit],
            fmt='o', markersize=2, color='mediumseagreen', alpha=0.25,
            elinewidth=0.5, capsize=0
        )

    # 3. Target data in fit region.
    if non_zero.size > 0:
        ax.errorbar(
            mass_target[non_zero], hmf_target[non_zero],
            yerr=hmf_err[non_zero],
            fmt='o', markersize=2, color='mediumseagreen',
            elinewidth=0.5, capsize=0
        )

    # 4. Despali et al. (2015) model.
    nz = hmf_despali2015 > 0
    if nz.any():
        ax.plot(mass_despali2015[nz], hmf_despali2015[nz],
                color='cornflowerblue', linestyle='--', linewidth=1)

    # 5. Bohr et al. (2021) model.
    nz = hmf_bohr2021 > 0
    if nz.any():
        ax.plot(mass_bohr2021[nz], hmf_bohr2021[nz],
                color='lightgray', linestyle='--', linewidth=1)

    # 6. True model (no numerical artifacts).
    nz = hmf_true > 0
    if nz.any():
        ax.plot(mass_true[nz], hmf_true[nz],
                color='orangered', linestyle='--', linewidth=1)

    # 7. Base (fitted) model.
    nz = hmf_model > 0
    if nz.any():
        ax.plot(mass_model[nz], hmf_model[nz],
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

    # Perl: force re-run when --select is given.
    if options.get('select') is not None:
        options['force'] = 'yes'

    # Build full entry list and attach derived file paths.
    entries = iterate(simulations, options)
    for entry in entries:
        identifier = (
            f"{entry['suite']['name']}_{entry['group']['name']}_"
            f"{entry['resolution']['name']}_{entry['simulation']['name']}_"
            f"{entry['realization']}_z{entry['redshift']}"
        )
        entry['fileParameters'] = options['outputDirectory'] + f'haloMassFunctionBase_{identifier}.xml'
        entry['filePrefix']     = options['outputDirectory'] + f'haloMassFunction_{identifier}'

    # Group by (suite, group, resolution, simulation, redshift) — mirrors Perl %setsRealization.
    sets_realization = {}
    for entry in entries:
        key = (
            f"{entry['suite']['name']}_{entry['group']['name']}_"
            f"{entry['resolution']['name']}_{entry['simulation']['name']}_"
            f"z{entry['redshift']}"
        )
        sets_realization.setdefault(key, []).append(entry)

    # Submit Galacticus jobs for all four model variants.
    jobs = _build_jobs(entries, options)
    if jobs:
        manager = queueManager.factory(args)
        submit_jobs(manager, jobs)

    # Generate comparison plots.
    for entry in entries:
        _plot_entry(entry, options)


if __name__ == '__main__':
    main()
