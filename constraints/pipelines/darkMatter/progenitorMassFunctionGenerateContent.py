#!/usr/bin/env python3
# Script to generate content for progenitor mass function constraint pipeline.
# Andrew Benson (01-May-2026)

import argparse
import math
import os
import re
import sys

import h5py
import lxml.etree as ET
import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))
from Galacticus.Constraints.Simulations import iterate, parse_simulations_xml
from XML.Utils import xml_to_dict

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Generate content for the progenitor mass function constraint pipeline.'
    )
    parser.add_argument('--pipelinePath',      required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--outputDirectory',   required=True,
                        help='Directory for generated output files')
    parser.add_argument('--countParticlesMinimum',       default=300, type=int,
                        help='Minimum particles per halo for constraints (default: 300)')
    parser.add_argument('--convergeAfterCount',   default=-1, type=int,
                        help='Number of steps after which to declare convergence (-1 to never stop; default: 1000)')
    parser.add_argument('--correlationStopAfterCount',   default=1000, type=int,
                        help='Number of correlation length after which to stop (default: 1000)')
    parser.add_argument('--select',            default=None, action='append',
                        help='Simulation selection filter (suite::group::resolution::...); may be repeated')
    parser.add_argument('--initializeToPosteriorMaximum', default=None,
                        help='Log file root for initializing from a prior posterior maximum')
    parser.add_argument('--randomize', default='no',
                        help='Randomize merger tree construction on each step (default: no)')
    parser.add_argument('--massRatioDistribution', default='hinkley',
                        choices=['normal', 'hinkley'],
                        help='Distribution model for the progenitor-to-parent mass ratio under N-body '
                             'mass errors: "hinkley" (exact ratio-of-normals) or "normal" (linearized '
                             'Gaussian approximation) (default: hinkley)')
    args = parser.parse_args()
    # Normalise paths to end with '/'.
    for key in ('pipelinePath', 'outputDirectory'):
        val = getattr(args, key)
        if val and not val.endswith('/'):
            setattr(args, key, val + '/')
    return vars(args)


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

def _base_files(entry_groups,options):
    """Generate progenitorMassFunctionBase_*.xml files and accumulate config_likelihood.

    Returns config_likelihood.
    """
    pipeline_path          = options['pipelinePath']
    output_dir             = options['outputDirectory']
    count_particles_min    = options['countParticlesMinimum']
    fraction_mass_primary  = 5.0
    data_path              = os.environ.get('GALACTICUS_DATA_PATH', '')

    config_likelihood      = ''
    for identifier in sorted(entry_groups):
        group   = entry_groups[identifier]
        entry   = group[0]
        suite   = entry['suite']
        grp     = entry['group']
        res     = entry['resolution']
        sim     = entry['simulation']
        real    = entry['realization']
        redshifts            = sorted([e['redshift'] for e in group],key=lambda x: float(x))
        redshift_progenitors = redshifts
        redshift_parent      = redshift_progenitors.pop(0)
        redshift_high        = np.ceil(float(redshift_progenitors[-1]))+1

        print(f"  {suite['name']}\t{grp['name']}\t{res['name']}\t{sim['name']}\t{real}\tz={redshift_parent} ⟶  {', '.join(redshift_progenitors)}")

        # --- Mass resolution ---
        mass_resolution = f"{count_particles_min / 5.0 * res['massParticle']:11.5e}"
        mass_resolution_comment = "1/5th of the minimum progenitor halo mass - such that our resolution is well below any mass scale that we care about"
        
        # --- Mass limits ---
        mass_halo_minimum = f"{count_particles_min * res['massParticle']:11.5e}"
        limits_str = suite.get('limitMassMaximum', {}).get('value', '')
        limits     = limits_str.split(':') if limits_str else []

        mass_maxima = []
        if 'primaryFraction' in limits:
            mass_maxima.append({
                'mass':   fraction_mass_primary * entry['massPrimary'],
                'reason': f'{fraction_mass_primary} of the target halo mass',
            })

        if sim['name'] != 'CDM' and 'matchedPhaseICs' in limits:
            redshift = sorted(redshifts, key=float)[0]
            matched  = suite['matchedICs']
            data_path_static = data_path + '/static/darkMatter/'
            file_self = (data_path_static
                         + f"{suite['name']}/{grp['name']}/{res['name']}"
                         + f"/{sim['name']}/{real}/progenitorMassFunction_z{redshift}.hdf5")
            file_ref  = (data_path_static
                         + f"{matched['suite']}/{grp['name']}/{res['name']}"
                         + f"/{matched['simulation']}/{real}/progenitorMassFunction_z{redshift}.hdf5")
            with h5py.File(file_self, 'r') as hf:
                count_self = hf['simulation0001/count'][:]
            with h5py.File(file_ref, 'r') as hf:
                count_ref  = hf['simulation0001/count'][:]
                mass       = hf['simulation0001/mass'][:]

            count_combined = count_ref.copy()
            count_combined[count_ref == 0] = count_self[count_ref == 0]

            non_zero = count_combined > 0
            delta    = (np.abs(count_self[non_zero].astype(float)
                               - count_ref[non_zero].astype(float))
                        / np.sqrt(count_combined[non_zero].astype(float)))

            i_maximum = len(delta) - 1
            for i, d in enumerate(delta):
                if d <= 1.0:
                    i_maximum = i
                    break

            mass_non_zero = mass[non_zero]
            mass_maximum  = mass_non_zero[i_maximum] / np.sqrt(mass[1] / mass[0])
            mass_maxima.append({
                'mass':   float(mass_maximum),
                'reason': f"Matched phases with {matched['suite']} {matched['simulation']}",
            })

        mass_halo_maximum        = '1.0e16'
        mass_halo_maximum_reason = ''
        if mass_maxima:
            best = min(mass_maxima, key=lambda m: m['mass'])
            mass_halo_maximum        = best['mass']
            mass_halo_maximum_reason = best['reason']

        # --- File names ---
        file_name_base = f"{output_dir}progenitorMassFunctionBase_{suite['name']}_{grp['name']}_{res['name']}_{sim['name']}_{real}_z{redshift_parent}.xml"
        file_name_target = f"%DATASTATICPATH%/darkMatter/{suite['name']}/{grp['name']}/{res['name']}/{sim['name']}/{real}/progenitorMassFunction_z{redshift_parent}.hdf5"
        # --- Non-CDM modifications ---
        changes = None
        changes_file_name = f"{output_dir}progenitorMassFunctionChanges_{suite['name']}_{grp['name']}_{res['name']}_{sim['name']}_{real}_z{redshift_parent}.xml"
        if sim['transferFunction']['value'] != "CDM":
            changes = (
                f'<changes>\n'
                f'  <change type="update" path="mergerTreeBranchingProbability/cdmAssumptions" value="false"/>\n'
                f'</changes>\n'
                )
        
        # --- Accumulate config_likelihood ---
        mass_minimum_comment = (f'{count_particles_min} times'
                                f" {suite['name']} {grp['name']} particle mass")
        config_likelihood += (
            f"    <!-- Suite: {suite['name']}; Group: {grp['name']};"
            f" Simulation: {sim['name']} -->\n"
            f'    <parameterMap value="mergerTreeBranchingProbability/G0'
            f'                                mergerTreeBranchingProbability/gamma1\n'
            f'                         mergerTreeBranchingProbability/gamma2'
            f'                                mergerTreeBranchingProbability/gamma3\n'
            f'\n'
            f'                         progenitorMassFunctionParameters/rootVarianceTargetFractional\n'
            f'\n'
            f'                         progenitorMassFunctionParameters/errorAmplitude\n'
            f'                         progenitorMassFunctionParameters/errorExponent\n'
            f'                         progenitorMassFunctionParameters/errorFractionalHighMass\n'
            f'                        "/>\n'
            )
        if options['randomize'] == 'yes':
            config_likelihood += f'    <parameterInactiveMap value="randomNumberGenerator/seed"/>\n'
        config_likelihood += (
            f'    <posteriorSampleLikelihood value="galaxyPopulation">\n'
            f'      <baseParametersFileName    value="{file_name_base}" />\n'
            )
        if changes is not None:
            config_likelihood += (
                f'      <changeParametersFileNames    value="{changes_file_name}" />\n'
                )
        config_likelihood += (
            f'      <evolveForestsVerbosity    value="silent"           />\n'
            f'      <reportFileName            value="false"            />\n'
            f'      <reportState               value="false"            />\n'
            f'      <countCollaborativeGroups  value="-1"               />\n'
            f'    </posteriorSampleLikelihood>\n'
        )

        # --- Write one base parameter file ---
        xi = 'xmlns:xi="http://www.w3.org/2001/XInclude"'
        xp = 'xpointer="xpointer(parameters/*)"'
        suite_n = suite['name']
        grp_n   = grp['name']
        res_n   = res['name']
        sim_n   = sim['name']
        base = (
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<parameters>\n'
            f'  <formatVersion>2</formatVersion>\n'
            f'\n'
            f'  <!-- Include cosmology and mass function parameters -->\n'
            f'  <xi:include href="progenitorMassFunctionParameters.xml"'
            f'                                                                       {xp} {xi}/>\n'
            f'  <xi:include href="{pipeline_path}simulation_{suite_n}_{grp_n}.xml"'
            f'            {xp} {xi}/>\n'
            f'  <xi:include href="{pipeline_path}cosmology_{suite_n}.xml"'
            f'                                           {xp} {xi}/>\n'
            f'  <xi:include href="{pipeline_path}powerSpectrum_{suite_n}.xml"'
            f'                              {xp} {xi}/>\n'
            f'  <xi:include href="{pipeline_path}transferFunction_{suite_n}_{sim_n}.xml"'
            f'                               {xp} {xi}/>\n'
            f'  <xi:include href="{pipeline_path}mergerTreeEvolution.xml"'
            f'                               {xp} {xi}/>\n'
            f'  <xi:include href="{output_dir}haloMassFunctionParameters.xml"'
            f'                                 {xp} {xi}/>\n'
            f'  <xi:include href="{output_dir}progenitorMassFunctionParameters.xml"'
            f'                                 {xp} {xi}/>\n'
            f'\n'
        )
        if options['randomize'] == 'yes':
            base += (
                f'<!-- Random number generator -->\n'
                f'<randomNumberGenerator value="GSL">\n'
                f'  <seed value="9372"/>\n'
                f'</randomNumberGenerator>\n'
                f'\n'
                )
        base += (
            f'<!-- Task control -->\n'
            f'<evolveForestsWorkShare value="cyclic"/>\n'
            f'\n'
            f'<!-- Mass resolution -->\n'
            f'<mergerTreeMassResolution value="fixed">\n'
            f'  <!-- {mass_resolution_comment} -->\n'
            f'  <massResolution value="{mass_resolution}"/>\n'
            f'</mergerTreeMassResolution>\n'
            f'\n'
            f'  <!-- Particle mass at the current resolution -->\n'
            f'  <massParticleAtResolution value="=[simulation/massParticle/{res_n}]"'
            f' ignoreWarnings="true"/>\n'
            f'\n'
            f'  <!-- N-body halo mass uncertainties: power-law fractional-error model, calibrated\n'
            f'       in the MCMC. The fractional error is written as a function of particle number\n'
            f'       N = M/m_p so that its amplitude is fixed at fixed N across resolutions:\n'
            f'         sigma_pn(M) = A * (N/1000)^gamma,   with N = M/m_p,\n'
            f'       reproduced by the powerLaw class (sigma = [sigma_12^2 (M/1e12)^2gamma + sigma_inf^2]^1/2)\n'
            f'       by setting the internal normalization sigma_12 = A * (1e9/m_p)^gamma. Here A\n'
            f'       (errorAmplitude) is the fractional error at N=1000 particles (Trenti et al. 2010: 0.135),\n'
            f'       gamma (errorExponent) the mass slope, and sigma_inf (errorFractionalHighMass) a\n'
            f'       resolution-independent floor. Correlation coefficients are fixed at Trenti values. -->\n'
            f'  <nbodyHaloMassError value="powerLaw">\n'
            f'    <normalization               value="=[progenitorMassFunctionParameters/errorAmplitude]*(1.0e9/[massParticleAtResolution])^([progenitorMassFunctionParameters/errorExponent])"/>\n'
            f'    <exponent                    value="=[progenitorMassFunctionParameters/errorExponent]"          />\n'
            f'    <fractionalErrorHighMass     value="=[progenitorMassFunctionParameters/errorFractionalHighMass]" />\n'
            f'    <correlationModelTrivial     value="false"/>\n'
            f'    <correlationNormalization    value="1.0"  />\n'
            f'    <correlationMassExponent     value="1.0"  />\n'
            f'    <correlationRedshiftExponent value="0.0"  />\n'
            f'  </nbodyHaloMassError>\n'
            f'\n'
            f'  <!-- Output control -->\n'
            f'  <outputFileName value="{output_dir}progenitorMassFunction'
            f'_{suite_n}_{grp_n}_{res_n}_{sim_n}_{real}_z{redshift_parent}.hdf5"/>\n'
            f'  <outputTimes value="list">\n'
            f'    <redshifts value="{" ".join(redshift_progenitors)}"/>\n'
            f'  </outputTimes>  \n'
            f'  <mergerTreeOutputter value="analyzer">\n'
            f'  </mergerTreeOutputter>\n'
            f'  <mergerTreeOperator value="regridTimes">\n'
            f'    <!-- Regrid the trees to massively speed up searching for parents etc. -->\n'
            f'    <outputTimes value="list">\n'
            f'      <redshifts value="{redshift_parent} {" ".join(redshift_progenitors)} {redshift_high}"/>\n'
            f'    </outputTimes>\n'
            f'  </mergerTreeOperator>\n'
        )
        # Add tree mass distribution settings.
        tree_masses_method = entry['resolution']['progenitorMassFunction']['treeMasses']['value']
        if tree_masses_method == "sampled":
            trees_per_decade = entry['resolution']['progenitorMassFunction']['treesPerDecade']['value']
            base += (
                f'<!-- Distribution of merger tree masses to build -->\n'
                f'<mergerTreeBuildMasses value="sampledDistributionUniform">\n'
                f'  <massTreeMinimum value="{mass_halo_minimum}"/>  <!-- {mass_minimum_comment} -->\n'
                f'  <massTreeMaximum value="{mass_halo_maximum}"/>  <!-- {mass_halo_maximum_reason} -->\n'
                f'  <treesPerDecade  value="{trees_per_decade}" />\n'
                f'</mergerTreeBuildMasses>\n'
                f'<mergerTreeBuildMassDistribution value="powerLaw">\n'
                f'  <exponent value="1.0"/>\n'
                f'</mergerTreeBuildMassDistribution>\n'
                )
        elif tree_masses_method == "file":
            file_name_mass = f'{data_path}/static/darkMatter/{suite_n}/{grp_n}/{res_n}/{sim_n}/hostHaloMasses_z{redshift_parent}.hdf5'
            count_replications = entry['resolution']['progenitorMassFunction']['countReplications']['value']
            base += (
                f'<!-- Replicated list of merger tree masses to build -->\n'
                f'<mergerTreeBuildMasses value="replicate">\n'
                f'  <replicationCount value="{count_replications}"/>\n'
                f'  <mergerTreeBuildMasses value="readHDF5">\n'
                f'    <fileName value="{file_name_mass}"/>\n'
                f'  </mergerTreeBuildMasses>\n'
                f'</mergerTreeBuildMasses>\n'
                )
        else:
            raise RuntimeError('tree masses method should be `sampled` or `file`; found `'+tree_masses_method+'`')
        # Access the target data file.
        file_name_target = file_name_target.replace('%DATASTATICPATH%', data_path + '/static/')
        with h5py.File(file_name_target, 'r') as data_target:
            mass_bins   = data_target['simulation0001/massParent'         ][:]
            mass_ratios = data_target['simulation0001/massRatioProgenitor'][:]
            count_bins  = data_target['simulation0001/count'              ][:]
        mass_bin_half_width_geometric = np.sqrt(mass_bins[1]/mass_bins[0])
        # Iterate over mass and redshift bins. Include only cases that have non-zero constraint.
        comment_range = '&lt; \log_{{10}}(M_\mathrm{{parent}}/\mathrm{{M}}_\odot) \le'
        base += (
            f'  <outputAnalysis value="multi">\n'
            f'    <virialDensityContrastDefinition value="sphericalCollapseClsnlssMttrCsmlgclCnstnt"/>\n'
            )
        for i in range(len(redshift_progenitors)):
            for j in range(len(mass_bins)):
                mass_bin_min     = mass_bins[j]/mass_bin_half_width_geometric
                mass_bin_max     = mass_bins[j]*mass_bin_half_width_geometric
                mass_ratio_min   = count_particles_min * float(res['massParticle']) / mass_bin_max
                mass_bin_min_log = np.log10(mass_bin_min)
                mass_bin_max_log = np.log10(mass_bin_max)
                in_range         = (mass_ratios >= mass_ratio_min) & (mass_ratios <= 1.0)
                if np.count_nonzero(count_bins[i,in_range,j]) == 0:
                    # Skip bins (of progenitor mass, redshift) that have no useable target data.
                    continue
                base += (
                    f'    <outputAnalysis value="progenitorMassFunction">\n'
                    f'      <massRatioDistribution        value="{options["massRatioDistribution"]}"                                                             />\n'
                    f'      <fileName                     value="{file_name_target}"                                                                            />\n'
                    f'      <indexRedshift                value="{i}"                                                                                           />\n'
                    f'      <indexParent                  value="{j}"                                                                                           />\n'
                    f'      <massRatioLikelihoodMinimum   value="{mass_ratio_min}"                                                                              />\n'
                    f'      <massRatioLikelihoodMaximum   value="1.00"                                                                                          />\n'
                    f'      <covarianceDiagonalize        value="true"                                                                                          />\n'
                    f'      <covarianceTargetOnly         value="true"                                                                                          />\n'
                    f'      <likelihoodInCounts           value="true"                                                                                          />\n'
                    f'      <likelihoodInLog              value="false"                                                                                         />\n'
                    f'      <fillInZeroBins               value="=[%s|progenitorMassFunctionParameters/fillInZeroBins]"                                         />\n'
                    f'      <rootVarianceTargetFractional value="=[progenitorMassFunctionParameters/rootVarianceTargetFractional]"                              />\n'
                    f'      <alwaysIsolatedOnly           value="true"                                                                                          />\n'
                    f'      <redshiftParent               value="{redshift_parent}"                                                                             />\n'
                    f'      <label                        value="M{j}Z{i}"                                                                                      />\n'
                    f'      <targetLabel                  value="{suite_n}; {grp_n}; {res_n}; {sim_n}"                                                          />\n'
                    f'      <comment                      value="${mass_bin_min_log:.2f} {comment_range} {mass_bin_max_log:.2f}$; $z={redshift_progenitors[i]}$"/>\n'
                    f'    </outputAnalysis>\n'
                    f'\n'
                )
        base += (
            f'  </outputAnalysis>\n'
            f'\n'
        )
        if suite.get('includeEnvironment', {}).get('value') == 'true':
            env = entry['environment']
            base += (
                f'\n'
                f'  <!-- Halo environments -->\n'
                f'  <haloEnvironment value="fixed">\n'
                f'    <massEnvironment value="{env["massEnvironment"]}"/>\n'
                f'    <overdensity     value="{env["overdensityEnvironment"]}"/>\n'
                f'  </haloEnvironment>\n'
            )
        base += '\n</parameters>\n'

        with open(file_name_base, 'w') as fh:
            fh.write(base)

    print('...done')
    return config_likelihood


# ---------------------------------------------------------------------------
# Assemble MCMC config XML strings
# ---------------------------------------------------------------------------

def _config_strings(options):
    """Return (config_opener, config_initializer, config_resumer,
               config_closer, parameters_opener, parameters_closer)."""
    import math

    output_dir   = options['outputDirectory']
    pipeline_path = options['pipelinePath']

    # Active model parameters. Assembled as a single block so that the parameter
    # count used to set the differential-evolution proposal size is *derived* from
    # the actual number of active parameters, never hardcoded. Add or remove a
    # <modelParameter value="active"> block here and gammaInitial tracks it.
    active_parameters = (
        '    <modelParameter value="active">\n'
        '      <name value="mergerTreeBranchingProbability/G0"/>\n'
        '      <label value="G_0" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="0.10"/>\n'
        '        <limitUpper value="3.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="mergerTreeBranchingProbability/gamma1"/>\n'
        '      <label value="\gamma_1" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="-1.00"/>\n'
        '        <limitUpper value="+0.99"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '     <modelParameter value="active">\n'
        '      <name value="mergerTreeBranchingProbability/gamma2" />\n'
        '      <label value="\gamma_2" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="-1.00" />\n'
        '        <limitUpper value="+1.00" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity" />\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="mergerTreeBranchingProbability/gamma3"/>\n'
        '      <label value="\gamma_3" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="+0.00"/>\n'
        '        <limitUpper value="+1.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="progenitorMassFunctionParameters/rootVarianceTargetFractional"/>\n'
        r'      <label value="\mathcal{C}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="+1.00e-4"/>\n'
        '        <limitUpper value="+1.00e+0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        # --- N-body halo mass-error (power-law) parameters ---
        '    <modelParameter value="active">\n'
        '      <name value="progenitorMassFunctionParameters/errorAmplitude"/>\n'
        r'      <label value="\sigma_{1000}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="+1.00e-2"/>\n'
        '        <limitUpper value="+1.00e+0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="progenitorMassFunctionParameters/errorExponent"/>\n'
        r'      <label value="\gamma_\sigma" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="-1.00"/>\n'
        '        <limitUpper value="+0.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="progenitorMassFunctionParameters/errorFractionalHighMass"/>\n'
        r'      <label value="\sigma_\infty" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="+1.00e-4"/>\n'
        '        <limitUpper value="+1.00e+0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0"/>\n'
        '        <scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
    )
    count_parameters = active_parameters.count('<modelParameter value="active">')
    gamma_initial    = 2.35 / math.sqrt(count_parameters)

    # ------------------------------------------------------------------ opener
    config_opener = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<parameters>\n'
        '  <!-- Posterior sampling simulation parameter file for constraining to dark matter progenitor mass functions. -->\n'
        '  <!-- Andrew Benson (01-May-2026)                                                                             -->\n'
        '  <formatVersion>2</formatVersion>\n'
        '\n'
        '  <verbosityLevel value="standard"/>\n'
        '\n'
        '  <task value="posteriorSample">\n'
        '    <initializeNodeClassHierarchy value="true"/>\n'
        '  </task>\n'
        '\n'
        f'  <outputFileName value="{output_dir}progenitorMassFunction.hdf5"/>\n'
        '\n'
        '  <!-- Likelihood -->\n'
        '  <posteriorSampleLikelihood value="independentLikelihoods">\n'
        '    <orderRotation value="byRankOnNode"/> <!-- Rotate likelihoods by the on-node rank to ensure that each process begins with a different power spectrum class. -->\n'
    )

    # --------------------------------------------------------------- initializer
    config_initializer = '  </posteriorSampleLikelihood>\n\n'

    prior_log = options.get('initializeToPosteriorMaximum')
    if prior_log is not None:
        config_initializer += (
            f'  <posteriorSampleStateInitialize value="posteriorMaximumGaussianSphere">\n'
            f'     <logFileRoot      value="{prior_log}"/>\n'
            f'     <radiusSphere     value="1.0e-6"/>\n'
            f'     <radiusIsRelative value="true"  />\n'
            f'     <posteriorSampleStateInitialize value="gaussianSphere">\n'
            f'        <radiusSphere     value="1.0e-6"/>\n'
            f'        <radiusIsRelative value="true"  />\n'
            f'     </posteriorSampleStateInitialize>\n'
            f'  </posteriorSampleStateInitialize>   \n\n'
        )
    else:
        config_initializer += (
            '  <posteriorSampleStateInitialize value="gaussianSphere">\n'
            '     <radiusSphere     value="1.0e-6"/>\n'
            '     <radiusIsRelative value="true"  />\n'
            '     <!-- Initialize to the Parkinson, Cole, and Helly (2008) fit for the branching\n'
            '          parameters, and to the Trenti et al. (2010)-equivalent error model. Order must\n'
            '          match the active <modelParameter> declarations: G0 gamma1 gamma2 gamma3\n'
            '          rootVarianceTargetFractional errorAmplitude errorExponent errorFractionalHighMass -->\n'
            '     <position         value="0.57 0.38 -0.01 0.00 1.0e-2 0.135 -0.333 0.01"/>\n'
            '  </posteriorSampleStateInitialize>   \n\n'
        )

    config_initializer += (
        f'  <posteriorSampleDffrntlEvltnProposalSize value="adaptive" >\n'
        f'    <gammaInitial          value="{gamma_initial:.5e}"/>\n'
        f'    <gammaAdjustFactor     value="1.10e+0"/>\n'
        f'    <gammaMinimum          value="1.00e-4"/>\n'
        f'    <gammaMaximum          value="3.00e+0"/>\n'
        f'    <acceptanceRateMinimum value="0.10e+0"/>\n'
        f'    <acceptanceRateMaximum value="0.90e+0"/>\n'
        f'    <updateCount           value="10"     />\n'
        f'    <appendLog             value="false"  />\n'
        f'    <flushLog              value="true"   />\n'
        f'    <restoreFromLog        value="false"  />\n'
        f'    <logFileName           value="{output_dir}progenitorMassFunctionGamma.log"/>\n'
        f'  </posteriorSampleDffrntlEvltnProposalSize>\n'
        f'\n'
        f'  <!-- MCMC -->\n'
        f'  <posteriorSampleSimulation value="differentialEvolution">\n'
        f'    <stepsMaximum           value="100000"                                  />\n'
        f'    <acceptanceAverageCount value="    10"                                  />\n'
        f'    <stateSwapCount         value="    11"                                  /> <!-- Offset swaps from reporting, otherwise we only get reports for swap steps, which gives a biased view of progress. -->\n'
        f'    <slowStepCount          value="    12"                                  />\n'
        f'    <logFileRoot            value="{output_dir}progenitorMassFunctionChains"/>\n'
        f'    <reportCount            value="    10"                                  />\n'
        f'    <sampleOutliers         value="false"                                   />\n'
        f'    <logFlushCount          value="     1"                                  />\n'
        f'    <appendLogs             value="false"                                   />\n'
        f'    <loadBalance            value="false"                                   />\n'
    )

    # ----------------------------------------------------------------- resumer
    config_resumer = (
        '  </posteriorSampleLikelihood>\n'
        '\n'
        f'  <posteriorSampleStateInitialize value="resume">\n'
        f'    <logFileRoot  value="{output_dir}progenitorMassFunctionChains"/>\n'
        f'    <restoreState value="true"                                    />\n'
        f'  </posteriorSampleStateInitialize>   \n'
        f'\n'
        f'  <posteriorSampleDffrntlEvltnProposalSize value="adaptive" >\n'
        f'    <gammaInitial          value="{gamma_initial:.5e}"/>\n'
        f'    <gammaAdjustFactor     value="1.10e+0"/>\n'
        f'    <gammaMinimum          value="1.00e-4"/>\n'
        f'    <gammaMaximum          value="3.00e+0"/>\n'
        f'    <acceptanceRateMinimum value="0.10e+0"/>\n'
        f'    <acceptanceRateMaximum value="0.90e+0"/>\n'
        f'    <updateCount           value="10"     />\n'
        f'    <appendLog             value="true"   />\n'
        f'    <flushLog              value="true"   />\n'
        f'    <restoreFromLog        value="true"   />\n'
        f'    <logFileName           value="{output_dir}progenitorMassFunctionGamma.log"/>\n'
        f'  </posteriorSampleDffrntlEvltnProposalSize>\n'
        f'\n'
        f'  <!-- MCMC -->\n'
        f'  <posteriorSampleSimulation value="differentialEvolution">\n'
        f'    <stepsMaximum           value="100000"                                  />\n'
        f'    <acceptanceAverageCount value="    10"                                  />\n'
        f'    <stateSwapCount         value="    11"                                  /> <!-- Offset swaps from reporting, otherwise we only get reports for swap steps, which gives a biased view of progress. -->\n'
        f'    <slowStepCount          value="    12"                                  />\n'
        f'    <logFileRoot            value="{output_dir}progenitorMassFunctionChains"/>\n'
        f'    <reportCount            value="    10"                                  />\n'
        f'    <sampleOutliers         value="false"                                   />\n'
        f'    <logFlushCount          value="     1"                                  />\n'
        f'    <appendLogs             value="true"                                    />\n'
        f'    <loadBalance            value="false"                                   />\n'
        f'\n'
    )

    # ----------------------------------------------------------------- closer
    config_closer = (
        '    <posteriorSampleState value="correlation">\n'
        '      <acceptedStateCount value="100"/>\n'
        '    </posteriorSampleState>\n'
        '\n'
        )
    if options['convergeAfterCount'] == -1:
        config_closer += (
            '    <posteriorSampleConvergence value="never"/>\n'
            '    \n'
        )
    elif options['convergeAfterCount'] > 0:
        config_closer += (
            '    <posteriorSampleConvergence value="stepCount">\n'
            f'       <countSteps value="{options["convergeAfterCount"]}"/>\n'
            '    </posteriorSampleConvergence>\n'
            '    \n'
        )
    else:
       raise RuntimeError('`convergeAfterCount` should be -1 or > 0')
    config_closer += (
        '    <posteriorSampleStoppingCriterion value="correlationLength">\n'
        f'      <stopAfterCount value="{options["correlationStopAfterCount"]}"/>\n'
        '    </posteriorSampleStoppingCriterion>\n'
        '\n'
        '    <posteriorSampleDffrntlEvltnRandomJump   value="adaptive"/>\n'
        '\n'
        '    <!-- Parameters of the dark matter progenitor mass function. -->\n'
        )
    config_closer += active_parameters
    config_closer += (
        '\n'
        )
    if options['randomize'] == 'yes':
        config_closer += (
            '    <modelParameter value="derived">\n'
            '      <name value="randomNumberGenerator/seed"/>\n'
            '      <definition value="1234+%[posteriorSimulationStep]"/>\n'
            '      <isInteger value="true"/>\n'
            '    </modelParameter>\n'
            '\n'
        )

    # Finish the closer
    config_closer += (
        '\n'
        '  </posteriorSampleSimulation>\n'
        '\n'
        '  <!-- Random seed -->\n'
        '  <randomNumberGenerator value="GSL">\n'
        '    <seed          value="9372"/>\n'
        '    <mpiRankOffset value="true"/>\n'
        '  </randomNumberGenerator>\n'
        '\n'
        '  <!-- Component selection -->\n'
        '  <componentBasic              value="standard"  />\n'
        '  <componentBlackHole          value="null"      />\n'
        '  <componentDarkMatterProfile  value="scaleFree" />\n'
        '  <componentDisk               value="null"      />\n'
        '  <componentHotHalo            value="null"      />\n'
        '  <componentSatellite          value="null"      />\n'
        '  <componentSpheroid           value="null"      />\n'
        '  <componentSpin               value="null"      />\n'
        '  <darkMatterProfileDMO        value="isothermal"/>\n'
        '\n'
        '</parameters>\n'
    )

    # -------------------------------------------------------- parameters_opener
    parameters_opener = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<parameters>\n'
        '  <formatVersion>2</formatVersion>\n'
        '\n'
        '  <!-- Controllable parameters of the progenitor mass function -->\n'
        '  <progenitorMassFunctionParameters value="" ignoreWarnings="true">\n'
        '     <!-- Parameter controlling the fractional variance due to model discrepancy -->\n'
        '     <rootVarianceTargetFractional value="0.1" />\n'
        '     <!-- Parameter controlling if empty bins should be filled in using extrapolation -->\n'
        '     <fillInZeroBins               value="true"/>\n'
        '     <!-- Power-law N-body halo mass-error model parameters (calibrated in the MCMC).\n'
        '          Parametrized by particle number so the amplitude is fixed at fixed N=M/m_p across\n'
        '          resolutions; the internal sigma_12 is derived per-resolution in the base file.\n'
        '          Defaults reproduce the Trenti et al. (2010) error. -->\n'
        '     <errorAmplitude              value="0.135"     /> <!-- A        : fractional error at N=1000 particles -->\n'
        '     <errorExponent               value="-0.333333" /> <!-- gamma    : power-law slope with particle number -->\n'
        '     <errorFractionalHighMass     value="0.01"      /> <!-- sigma_inf: resolution-independent floor         -->\n'
    )

    # -------------------------------------------------------- parameters_closer
    parameters_closer = (
        '  </progenitorMassFunctionParameters>\n'
        '\n'
        '</parameters>\n'
    )

    return (config_opener, config_initializer, config_resumer,
            config_closer, parameters_opener, parameters_closer)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    options = _parse_args()
    simulations = parse_simulations_xml(
        options['pipelinePath'] + 'simulations.xml'
    )

    print('Building entry groups...')
    entry_groups = _group_entries(simulations, options)

    print('Generating base parameter files...')
    config_likelihood = _base_files(entry_groups,options)

    print('Assembling MCMC config strings...')
    (config_opener, config_initializer, config_resumer,
     config_closer, parameters_opener, parameters_closer) = _config_strings(options)

    output_dir = options['outputDirectory']

    with open(output_dir + 'progenitorMassFunctionConfig.xml', 'w') as fh:
        fh.write(config_opener)
        fh.write(config_likelihood)
        fh.write(config_initializer)
        fh.write(config_closer)

    with open(output_dir + 'progenitorMassFunctionConfigResume.xml', 'w') as fh:
        config_likelihood_resume = re.sub(r'(<appendSamples\s+value=")false("\/>)',r'\1true\2',config_likelihood)
        fh.write(config_opener)
        fh.write(config_likelihood_resume)
        fh.write(config_resumer)
        fh.write(config_closer)

    with open(output_dir + 'progenitorMassFunctionParameters.xml', 'w') as fh:
        fh.write(parameters_opener)
        fh.write(parameters_closer)

    print('Done.')


if __name__ == '__main__':
    main()
