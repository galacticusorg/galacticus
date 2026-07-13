#!/usr/bin/env python3
"""Script to generate content for halo mass function constraint pipeline.

Andrew Benson (2026)
"""

import argparse
import math
import os
import re

import h5py
import lxml.etree as ET
import numpy as np

from Galacticus.Constraints.Simulations import iterate, parse_simulations_xml
from Galacticus._logging                 import configure_default as _configure_default
from XML.Utils import xml_to_dict

# Show INFO-level diagnostic output from library modules.
_configure_default()


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def _remove_hmf_nodes(root, removal_value):
    """Splice out <haloMassFunction value=removal_value> nodes from the chain.

    The removed node is discarded along with any non-haloMassFunction
    children it carries; only its inner <haloMassFunction> child is promoted
    into its place.
    """
    parent = root
    child  = parent.find('haloMassFunction')
    while child is not None:
        next_child = child.find('haloMassFunction')
        if child.get('value') == removal_value:
            idx = list(parent).index(child)
            parent.remove(child)
            if next_child is not None:
                parent.insert(idx, next_child)
            child = next_child          # parent stays; advance along chain
        else:
            parent = child
            child  = next_child


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Generate content for the halo mass function constraint pipeline.'
    )
    parser.add_argument('--pipelinePath',      required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--outputDirectory',   required=True,
                        help='Directory for generated output files')
    parser.add_argument('--binAverage',                  default='true',
                        help='Average halo mass functions over each bin (default: true)')
    parser.add_argument('--includeCorrelations',         default='true',
                        help='Include correlations between HMFs at different redshifts (default: true)')
    parser.add_argument('--removeAccelerator',           default='false',
                        help='Remove accelerator haloMassFunction components (default: false)')
    parser.add_argument('--removeDetectionEfficiency',   default='false',
                        help='Remove detectionEfficiency haloMassFunction components (default: false)')
    parser.add_argument('--removeErrorConvolved',        default='false',
                        help='Remove errorConvolved haloMassFunction components (default: false)')
    parser.add_argument('--removeSimulationVariance',    default='false',
                        help='Remove simulationVariance haloMassFunction components (default: false)')
    parser.add_argument('--removeMultiplier',            default='false',
                        help='Remove multiplier haloMassFunction components (default: false)')
    parser.add_argument('--heatRepeats',                 default='true',
                        help='Heat repeated realizations to account for correlations (default: true)')
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
    args = parser.parse_args()
    # Normalise paths to end with '/'.
    for key in ('pipelinePath', 'outputDirectory'):
        val = getattr(args, key)
        if val and not val.endswith('/'):
            setattr(args, key, val + '/')
    return vars(args)


# ---------------------------------------------------------------------------
# Steps A & B
# ---------------------------------------------------------------------------

def _step_a_group_labels(simulations, options):
    """Step A: build isolation-bias and perturbation label dicts (stop_after='group').

    Returns (isolation_biases, perturbations) and mutates each entry['group']
    to carry 'labelIsolationBias', 'isolationBias', 'labelPerturbation', 'parameterPerturbation;
    and 'perturbation' keys used downstream.
    """
    isolation_biases = {}
    perturbations    = {}
    have_groups      = False

    for entry in iterate(simulations, options, stop_after='group'):
        have_groups = True
        suite = entry['suite']
        group = entry['group']
        print(f"  {suite['name']}\t{group['name']}")

        # Isolation bias label.
        if (suite.get('includeIsolationBias', {}).get('value') == 'true'
                and options.get('removeMultiplier') != 'true'):
            group_short = group.get('shortName', group['name'])
            if 'matchedIsolation' in suite:
                label = suite['matchedIsolation']['suite'] + group_short
            else:
                label = suite['name'] + group_short
            label = label.replace(':', '_')
            group['labelIsolationBias'] = label
            group['isolationBias'] = (
                f'haloMassFunctionParameters/isolationBias{label}'
                f'  haloMassFunctionParameters/isolationBiasExponent{label}'
            )
            isolation_biases[label] = {'groupName': group_short}
        else:
            group['labelIsolationBias'] = ''
            group['isolationBias']      = ''

        # Perturbation label.
        if (suite.get('includePerturbation', {}).get('value') == 'true'
                and options.get('removeSimulationVariance') != 'true'):
            label = 'cube,' + group['name']
            labelParameter = 'Cube' + group['name']
            group['labelPerturbation']     = label
            group['parameterPerturbation'] = labelParameter
            group['perturbation']          = f'haloMassFunctionParameters/perturbation{labelParameter}'
            perturbations[labelParameter]  = perturbations.get(label, 0) + 1
        else:
            group['labelPerturbation']     = ''
            group['parameterPerturbation'] = ''
            group['perturbation']          = ''

    if not have_groups:
        raise RuntimeError('no groups match this selection')

    return isolation_biases, perturbations


def _step_b_suite_hmf_xml(simulations, options):
    """Step B: copy and optionally prune haloMassFunction_{suite}.xml files (stop_after='suite')."""
    pipeline_path  = options['pipelinePath']
    output_dir     = options['outputDirectory']
    removals = {
        'accelerator':         options.get('removeAccelerator')         == 'true',
        'detectionEfficiency': options.get('removeDetectionEfficiency') == 'true',
        'errorConvolved':      options.get('removeErrorConvolved')      == 'true',
        'simulationVariance':  options.get('removeSimulationVariance')  == 'true',
        'multiplier':          options.get('removeMultiplier')          == 'true',
    }

    for entry in iterate(simulations, options, stop_after='suite'):
        suite_name = entry['suite']['name']
        src  = pipeline_path  + f'haloMassFunction_{suite_name}.xml'
        dst  = output_dir     + f'haloMassFunction_{suite_name}.xml'
        root = ET.parse(src).getroot()
        for removal_value, enabled in removals.items():
            if enabled:
                _remove_hmf_nodes(root, removal_value)
        ET.indent(root)
        ET.ElementTree(root).write(dst, xml_declaration=True, encoding='utf-8')

# ---------------------------------------------------------------------------
# Steps C & D
# ---------------------------------------------------------------------------

def _step_c_count_realizations(simulations, options):
    """Step C: count how many simulations share each (suite,group,res,realization,redshift) key."""
    count_realizations = {}
    for entry in iterate(simulations, options):
        key = '::'.join([
            entry['suite'     ]['name'],
            entry['group'     ]['name'],
            entry['resolution']['name'],
            entry['realization'],
            entry['redshift'],
        ])
        count_realizations[key] = count_realizations.get(key, 0) + 1
    return count_realizations


def _step_d_group_entries(simulations, options, count_realizations):
    """Step D: group entries that differ only by redshift into ordered lists.

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
        real_key = '::'.join([
            entry['suite'     ]['name'],
            entry['group'     ]['name'],
            entry['resolution']['name'],
            entry['realization'],
            entry['redshift'],
        ])
        entry['countRealizations'] = count_realizations.get(real_key, 1)
        entry_groups.setdefault(group_key, []).append(entry)

    if not have_models:
        raise RuntimeError('no models match this selection')

    return entry_groups


# ---------------------------------------------------------------------------
# Step E
# ---------------------------------------------------------------------------

def _step_e_base_files(entry_groups, options):
    """Step E: generate haloMassFunctionBase_*.xml files and accumulate config_likelihood.

    Returns (config_likelihood, detection_efficiency_classes).
    """
    pipeline_path          = options['pipelinePath']
    output_dir             = options['outputDirectory']
    count_particles_min    = options['countParticlesMinimum']
    fraction_mass_primary  = 0.1
    bin_average            = options['binAverage']
    include_correlations   = options['includeCorrelations']
    data_path              = os.environ.get('GALACTICUS_DATA_PATH', '')

    config_likelihood          = ''
    detection_efficiency_classes = {}

    for identifier in sorted(entry_groups):
        group   = entry_groups[identifier]
        entry   = group[0]
        suite   = entry['suite']
        grp     = entry['group']
        res     = entry['resolution']
        sim     = entry['simulation']
        real    = entry['realization']
        redshifts = [e['redshift'] for e in group]

        print(f"  {suite['name']}\t{grp['name']}\t{res['name']}\t{sim['name']}\t{real}\tz={', '.join(redshifts)}")

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
                         + f"/{sim['name']}/{real}/haloMassFunction_z{redshift}.hdf5")
            file_ref  = (data_path_static
                         + f"{matched['suite']}/{grp['name']}/{res['name']}"
                         + f"/{matched['simulation']}/{real}/haloMassFunction_z{redshift}.hdf5")
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

        mass_halo_maximum        = None
        mass_halo_maximum_reason = ''
        if mass_maxima:
            best = min(mass_maxima, key=lambda m: m['mass'])
            mass_halo_maximum        = best['mass']
            mass_halo_maximum_reason = best['reason']

        # --- File names ---
        file_names_base = [
            output_dir + f"haloMassFunctionBase_{suite['name']}_{grp['name']}"
            + f"_{res['name']}_{sim['name']}_{real}_z{z}.xml"
            for z in redshifts
        ]
        file_names_target = [
            f"%DATASTATICPATH%/darkMatter/{suite['name']}/{grp['name']}"
            + f"/{res['name']}/{sim['name']}/{real}/haloMassFunction_z{z}.hdf5"
            for z in redshifts
        ]

        # --- Detection efficiency class ---
        suite_name_det = suite.get('matchedDetection', {}).get('suite', suite['name'])
        suite_name_det = suite_name_det.replace(':', '')
        class_name     = suite_name_det + grp.get('detectionEfficiencyClass', '')
        if options.get('removeDetectionEfficiency') != 'true':
            detection_efficiency_classes[class_name] = \
                detection_efficiency_classes.get(class_name, 0) + 1
            detection_efficiency1 = (
                f'haloMassFunctionParameters/massMinimumParticleCount{class_name}'
                f' haloMassFunctionParameters/efficiencyAtMassMinimum{class_name}'
            )
            detection_efficiency2 = (
                f'haloMassFunctionParameters/exponentMassDetection{class_name}'
                f'    haloMassFunctionParameters/exponentRedshiftDetection{class_name}'
            )
        else:
            detection_efficiency1 = ''
            detection_efficiency2 = ''

        # --- Likelihood heating ---
        likelihood_heating_opener = ''
        likelihood_heating_closer = ''
        likelihood_temperature    = entry.get('options', {}).get('likelihoodTemperature')
        if likelihood_temperature is not None or options.get('heatRepeats') == 'true':
            temperature = 1.0
            if likelihood_temperature is not None:
                temperature *= float(likelihood_temperature)
            if options.get('heatRepeats') == 'true':
                temperature *= entry['countRealizations']
            likelihood_heating_opener = (
                f'    <posteriorSampleLikelihood value="heated">\n'
                f'       <temperature value="{temperature}"/>\n'
            )
            likelihood_heating_closer = '</posteriorSampleLikelihood>\n'

        # --- Accumulate config_likelihood ---
        mass_minimum_comment = (f'{count_particles_min} times zoom-in'
                                f" {suite['name']} {grp['name']} particle mass")
        config_likelihood += (
            f"    <!-- Suite: {suite['name']}; Group: {grp['name']};"
            f" Simulation: {sim['name']} -->\n"
            f'    <parameterMap value="haloMassFunctionParameters/a'
            f'                                haloMassFunctionParameters/b\n'
            f'                         haloMassFunctionParameters/c'
            f'                                haloMassFunctionParameters/p\n'
            f'                         haloMassFunctionParameters/q'
            f'                                haloMassFunctionParameters/normalization \n'
            f'\n'
            f'                         haloMassFunctionParameters/cW0'
            f'                              haloMassFunctionParameters/beta0\n'
            f'                         haloMassFunctionParameters/cW1'
            f'                              haloMassFunctionParameters/beta1\n'
            f'                         haloMassFunctionParameters/wavenumberScaledMinimum'
            f'          haloMassFunctionParameters/powerSpectrumSmoothingWidth\n'
            f'\n'
            f'                         haloMassFunctionParameters/artificialExponentMass'
            f'           haloMassFunctionParameters/artificialExponentGrowthFactor\n'
            f'                         haloMassFunctionParameters/artificialCountParticles\n'
            f'\n'
            f'                         {detection_efficiency1}\n'
            f'                         {detection_efficiency2}\n'
            f' \n'
            f'                         varianceFractionalModelDiscrepancy\n'
            f'\n'
            f"                         {grp['perturbation']}\n"
            f"                         {grp['isolationBias']}\n"
            f'                        "/>\n'
            f'    <parameterInactiveMap value=""     ignoreWarnings="true"/>\n'
            f'{likelihood_heating_opener}'
            f'    <posteriorSampleLikelihood value="haloMassFunction">\n'
            f'      <!-- Options matched to those of Benson (2017;'
            f' https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.3454B) -->\n'
            f'      <baseParametersFileName value="{file_names_base[0]}"/>\n'
            f'      <fileNames              value="{" ".join(file_names_target)}"/>\n'
            f'      <pathSamples            value="{output_dir}samples"/>\n'
            f'      <appendSamples          value="false"/>\n'
            f'      <redshifts              value="{" ".join(redshifts)}"/>\n'
            f'      <allowEmptyMassFunction value="true"                           />\n'
            f'      <massRangeMinimum       value="{mass_halo_minimum}"'
            f'             /> <!-- {mass_minimum_comment} -->\n'
        )
        if limits_str:
            config_likelihood += (
                f'      <massRangeMaximum       value="{mass_halo_maximum:.5e}"/>'
                f' <!-- {mass_halo_maximum_reason} -->\n'
            )
        config_likelihood += (
            f'      <binCountMinimum        value="0"                     />\n'
            f'      <likelihoodPoisson      value="true"                  />\n'
            f'      <binAverage             value="{bin_average}"         />\n'
            f'      <includeCorrelations    value="{include_correlations}"/>\n'
            f'    </posteriorSampleLikelihood>\n'
            f'{likelihood_heating_closer}'
        )

        # --- Write one base parameter file per redshift ---
        xi = 'xmlns:xi="http://www.w3.org/2001/XInclude"'
        xp = 'xpointer="xpointer(parameters/*)"'
        for i_z, redshift in enumerate(redshifts):
            suite_n = suite['name']
            grp_n   = grp['name']
            res_n   = res['name']
            sim_n   = sim['name']
            base = (
                f'<?xml version="1.0" encoding="UTF-8"?>\n'
                f'<parameters>\n'
                f'  <formatVersion>2</formatVersion>\n'
                f'\n'
                f'  <!-- Output control -->\n'
                f'  <outputFileName value="{output_dir}haloMassFunction'
                f'_{suite_n}_{grp_n}_{res_n}_{sim_n}_{real}_z{redshift}.hdf5"/>\n'
                f'  <outputTimes value="list">\n'
                f'    <redshifts value="{redshift}"/>\n'
                f'  </outputTimes>  \n'
                f'\n'
                f'  <!-- Include cosmology and mass function parameters -->\n'
                f'  <xi:include href="haloMassFunctionParameters.xml"'
                f'                                                                       {xp} {xi}/>\n'
                f'  <xi:include href="{pipeline_path}simulation_{suite_n}_{grp_n}.xml"'
                f'            {xp} {xi}/>\n'
                f'  <xi:include href="{pipeline_path}cosmology_{suite_n}.xml"'
                f'                                           {xp} {xi}/>\n'
                f'  <xi:include href="{pipeline_path}powerSpectrumEffective_{suite_n}.xml"'
                f'                              {xp} {xi}/>\n'
                f'  <xi:include href="{output_dir}haloMassFunction_{suite_n}.xml"'
                f'                                 {xp} {xi}/>\n'
                f'  <xi:include href="{pipeline_path}transferFunction_{suite_n}_{sim_n}.xml"'
                f' {xp} {xi}/>\n'
                f'\n'
                f'  <!-- Particle mass at the current resolution -->\n'
                f'  <massParticleAtResolution value="=[simulation/massParticle/{res_n}]"'
                f' ignoreWarnings="true"/>\n'
                f'\n'
                f'  <!-- Detection efficiency -->\n'
                f'  <detectionMassMinimumParticleCount'
                f' value="=[haloMassFunctionParameters/massMinimumParticleCount{class_name}]"'
                f'  ignoreWarnings="true"/>\n'
                f'  <detectionEfficiencyAtMassMinimum '
                f' value="=[haloMassFunctionParameters/efficiencyAtMassMinimum{class_name}]"'
                f'   ignoreWarnings="true"/>\n'
                f'  <detectionExponentMass            '
                f' value="=[haloMassFunctionParameters/exponentMassDetection{class_name}]"'
                f'     ignoreWarnings="true"/>\n'
                f'  <detectionExponentRedshift        '
                f' value="=[haloMassFunctionParameters/exponentRedshiftDetection{class_name}]"'
                f' ignoreWarnings="true"/>\n'
            )
            if (suite.get('includePerturbation', {}).get('value') == 'true'
                    and options.get('removeSimulationVariance') != 'true'):
                lbl = grp['parameterPerturbation']
                base += (
                    f'\n'
                    f'  <!-- Use the simulation variance perturbation relevant to this simulation -->\n'
                    f'  <perturbationFractional'
                    f' value="=[haloMassFunctionParameters/perturbation{lbl}]"/>\n'
                )
            if (suite.get('includeIsolationBias', {}).get('value') == 'true'
                    and options.get('removeMultiplier') != 'true'):
                lbl = grp['labelIsolationBias']
                base += (
                    f'\n'
                    f'  <!-- Isolation bias -->\n'
                    f'  <isolationBias        '
                    f' value="=[haloMassFunctionParameters/isolationBias{lbl}]"'
                    f'         ignoreWarnings="true"/>\n'
                    f'  <isolationBiasExponent'
                    f' value="=[haloMassFunctionParameters/isolationBiasExponent{lbl}]"'
                    f' ignoreWarnings="true"/>\n'
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

            with open(file_names_base[i_z], 'w') as fh:
                fh.write(base)

    print('...done')
    return config_likelihood, detection_efficiency_classes


# ---------------------------------------------------------------------------
# Step F — assemble MCMC config XML strings
# ---------------------------------------------------------------------------

def _step_f_config_strings(perturbations, isolation_biases,
                           detection_efficiency_classes, options):
    """Return (config_opener, config_initializer, config_resumer,
               config_closer, parameters_opener, parameters_closer)."""
    import math

    output_dir   = options['outputDirectory']
    pipeline_path = options['pipelinePath']

    count_parameters = (16
                        + len(perturbations)
                        + 2 * len(isolation_biases)
                        + 4 * len(detection_efficiency_classes))
    gamma_initial   = 2.35 / math.sqrt(count_parameters)
    recompute_count = ('<recomputeCount value="13"/>'
                       if options.get('includeCorrelations') == 'true'
                       else '')

    # ------------------------------------------------------------------ opener
    config_opener = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<parameters>\n'
        '  <!-- Posterior sampling simulation parameter file for constraining to dark matter halo mass functions. -->\n'
        '  <!-- Andrew Benson (17-September-2020)                                                                 -->\n'
        '  <formatVersion>2</formatVersion>\n'
        '\n'
        '  <verbosityLevel value="standard"/>\n'
        '\n'
        '  <task value="posteriorSample">\n'
        '    <initializeNodeClassHierarchy value="true"/>\n'
        '  </task>\n'
        '\n'
        f'  <outputFileName value="{output_dir}haloMassFunction.hdf5"/>\n'
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
        f'    <logFileName           value="{output_dir}haloMassFunctionGamma.log"/>\n'
        f'  </posteriorSampleDffrntlEvltnProposalSize>\n'
        f'\n'
        f'  <!-- MCMC -->\n'
        f'  <posteriorSampleSimulation value="differentialEvolution">\n'
        f'    <stepsMaximum           value="100000"                                  />\n'
        f'    <acceptanceAverageCount value="    10"                                  />\n'
        f'    <stateSwapCount         value="    11"                                  /> <!-- Offset swaps from reporting, otherwise we only get reports for swap steps, which gives a biased view of progress. -->\n'
        f'    <slowStepCount          value="    12"                                  />\n'
        f'    {recompute_count}\n'
        f'    <logFileRoot            value="{output_dir}haloMassFunctionChains"/>\n'
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
        f'    <logFileRoot  value="{output_dir}haloMassFunctionChains"/>\n'
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
        f'    <logFileName           value="{output_dir}haloMassFunctionGamma.log"/>\n'
        f'  </posteriorSampleDffrntlEvltnProposalSize>\n'
        f'\n'
        f'  <!-- MCMC -->\n'
        f'  <posteriorSampleSimulation value="differentialEvolution">\n'
        f'    <stepsMaximum           value="100000"                                  />\n'
        f'    <acceptanceAverageCount value="    10"                                  />\n'
        f'    <stateSwapCount         value="    11"                                  /> <!-- Offset swaps from reporting, otherwise we only get reports for swap steps, which gives a biased view of progress. -->\n'
        f'    <slowStepCount          value="    12"                                  />\n'
        f'    <logFileRoot            value="{output_dir}haloMassFunctionChains"/>\n'
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
        '    <!-- Parameters of the dark matter halo mass function. -->\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/a"/>\n'
        '      <label value="a" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value=" 0.03"/>\n'
        '\t<limitUpper value="10.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/b"/>\n'
        '      <label value="b" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value="-3.00"/>\n'
        '\t<limitUpper value="+3.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '     <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/c" />\n'
        '      <label value="c" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="+0.10" />\n'
        '        <limitUpper value="+5.00" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/p"/>\n'
        '      <label value="p" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value="-3.0"/>\n'
        '\t<limitUpper value="+3.0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/q"/>\n'
        '      <label value="q" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value="-3.00"/>\n'
        '\t<limitUpper value="+3.00"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/normalization"/>\n'
        '      <label value="A" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '\t<limitLower value="1.0e-3"/>\n'
        '\t<limitUpper value="1.0e+3"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '\n'
        '    <!-- Artificial halo model parameters -->\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/artificialExponentMass"/>\n'
        r'      <label value="\alpha" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value="-3.0"/>\n'
        '\t<limitUpper value="+0.0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/artificialExponentGrowthFactor"/>\n'
        r'      <label value="\beta" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '\t<limitLower value="-3.0"/>\n'
        '\t<limitUpper value="+3.0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/artificialCountParticles"/>\n'
        '      <label value="N" ignoreWarnings="true"/>\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '\t<limitLower value="1.0e+0"/>\n'
        '\t<limitUpper value="1.0e+3"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '\t<median value="0.0"/>\n'
        '\t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
        '\n'
        '    <!-- Window function parameters -->\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/cW0" />\n'
        r'      <label value="c_\mathrm{W,0}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="0.50" />\n'
        '        <limitUpper value="6.00" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/beta0" />\n'
        r'      <label value="\beta_0" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value=" 0.50" />\n'
        '        <limitUpper value="10.00" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/cW1" />\n'
        r'      <label value="c_\mathrm{W,1}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="0.1" />\n'
        '        <limitUpper value="5.0" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/beta1" />\n'
        r'      <label value="\beta_1" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="0.1" />\n'
        '        <limitUpper value="5.0" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/wavenumberScaledMinimum" />\n'
        r'      <label value="x_\mathrm{min}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="uniform">\n'
        '        <limitLower value="0.0" />\n'
        '        <limitUpper value="5.0" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="identity" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '    <modelParameter value="active">\n'
        '      <name value="haloMassFunctionParameters/powerSpectrumSmoothingWidth" />\n'
        r'      <label value="\log k_\mathrm{width}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '        <median value="0.0" />\n'
        '        <scale value="1.0e-4" />\n'
        '      </distributionFunction1DPerturber>\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '        <limitLower value="0.01" />\n'
        '        <limitUpper value="10.0" />\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm" />\n'
        '      <slow value="false" />\n'
        '    </modelParameter>\n'
        '\n'
        '    <modelParameter value="active">\n'
        '      <name value="varianceFractionalModelDiscrepancy"/>\n'
        r'      <label value="\mathcal{C}_\mathrm{disc}" ignoreWarnings="true"/>' + '\n'
        '      <distributionFunction1DPrior value="logUniform">\n'
        '    \t<limitLower value="1.0e-6"/>\n'
        '    \t<limitUpper value="1.0e+0"/>\n'
        '      </distributionFunction1DPrior>\n'
        '      <operatorUnaryMapper value="logarithm"/>\n'
        '      <distributionFunction1DPerturber value="cauchy">\n'
        '    \t<median value="0.0"/>\n'
        '    \t<scale value="1.0e-4"/>\n'
        '      </distributionFunction1DPerturber>\n'
        '    </modelParameter>\n'
    )

    # Perturbation parameters
    if perturbations:
        config_closer += '\n    <!-- Perturbation model parameters -->\n'
        for lbl in sorted(perturbations):
            config_closer += (
                f'    <modelParameter value="active">\n'
                f'      <name value="haloMassFunctionParameters/perturbation{lbl}" />\n'
                r'      <label value="\epsilon_\mathrm{' + lbl + r'}" ignoreWarnings="true"/>' + '\n'
                f'      <distributionFunction1DPerturber value="cauchy">\n'
                f'        <median value="0.0" />\n'
                f'        <scale value="1.0e-4" />\n'
                f'      </distributionFunction1DPerturber>\n'
                f'      <distributionFunction1DPrior value="normal">\n'
                f'        <limitLower value="-10.0" />\n'
                f'        <limitUpper value="+10.0" />\n'
                f'        <mean value="0.0" />\n'
                f'        <variance value="1.0" />\n'
                f'      </distributionFunction1DPrior>\n'
                f'      <operatorUnaryMapper value="identity" />\n'
                f'    </modelParameter>\n'
            )

    # Isolation bias parameters
    if isolation_biases:
        config_closer += '\n    <!-- Isolation bias model parameters -->\n'
        for lbl in sorted(isolation_biases):
            group_name = isolation_biases[lbl]['groupName']
            config_closer += (
                f'    <modelParameter value="active">\n'
                f'      <name value="haloMassFunctionParameters/isolationBias{lbl}" />\n'
                r'      <label value="\mathcal{I}_\mathrm{' + group_name + r'}" ignoreWarnings="true"/>' + '\n'
                f'      <distributionFunction1DPerturber value="cauchy">\n'
                f'        <median value="0.0" />\n'
                f'        <scale value="1.0e-4" />\n'
                f'      </distributionFunction1DPerturber>\n'
                f'      <distributionFunction1DPrior value="logNormal">\n'
                f'        <limitLower value="1.0e-3" />\n'
                f'        <limitUpper value="1.0e+3" />\n'
                f'        <mean value="1.0" />\n'
                f'        <variance value="0.25" />\n'
                f'      </distributionFunction1DPrior>\n'
                f'      <operatorUnaryMapper value="identity" />\n'
                f'    </modelParameter>    \n'
                f'    <modelParameter value="active">\n'
                f'      <name value="haloMassFunctionParameters/isolationBiasExponent{lbl}" />\n'
                r'      <label value="\alpha_\mathrm{iso,' + group_name + r'}" ignoreWarnings="true"/>' + '\n'
                f'      <distributionFunction1DPerturber value="cauchy">\n'
                f'        <median value="0.0" />\n'
                f'        <scale value="1.0e-4" />\n'
                f'      </distributionFunction1DPerturber>\n'
                f'      <distributionFunction1DPrior value="uniform">\n'
                f'        <limitLower value="-3.0" />\n'
                f'        <limitUpper value="+0.0" />\n'
                f'      </distributionFunction1DPrior>\n'
                f'      <operatorUnaryMapper value="identity" />\n'
                f'    </modelParameter>\n'
            )

    # Detection efficiency parameters
    if detection_efficiency_classes:
        config_closer += '\n    <!-- Detection efficiency model parameters -->\n'
    for cls in sorted(detection_efficiency_classes):
        config_closer += (
            f'    <modelParameter value="active">\n'
            f'      <name value="haloMassFunctionParameters/massMinimumParticleCount{cls}" />\n'
            r'      <label value="N_\mathrm{min,' + cls + r'}" ignoreWarnings="true"/>' + '\n'
            f'      <distributionFunction1DPerturber value="cauchy">\n'
            f'        <median value="0.0" />\n'
            f'        <scale value="1.0e-4" />\n'
            f'      </distributionFunction1DPerturber>\n'
            f'      <distributionFunction1DPrior value="logUniform">\n'
            f'        <limitLower value="  1.00" />\n'
            f'        <limitUpper value="100.00" />\n'
            f'      </distributionFunction1DPrior>\n'
            f'      <operatorUnaryMapper value="logarithm" />\n'
            f'    </modelParameter>\n'
            f'    <modelParameter value="active">\n'
            f'      <name value="haloMassFunctionParameters/efficiencyAtMassMinimum{cls}" />\n'
            r'      <label value="\epsilon_\mathrm{min,' + cls + r'}" ignoreWarnings="true"/>' + '\n'
            f'      <distributionFunction1DPerturber value="cauchy">\n'
            f'        <median value="0.0" />\n'
            f'        <scale value="1.0e-4" />\n'
            f'      </distributionFunction1DPerturber>\n'
            f'      <distributionFunction1DPrior value="uniform">\n'
            f'        <limitLower value="+0.00" />\n'
            f'        <limitUpper value="+1.00" />\n'
            f'      </distributionFunction1DPrior>\n'
            f'      <operatorUnaryMapper value="identity" />\n'
            f'    </modelParameter>\n'
            f'    <modelParameter value="active">\n'
            f'      <name value="haloMassFunctionParameters/exponentMassDetection{cls}" />\n'
            r'      <label value="\alpha_\mathrm{det,' + cls + r'}" ignoreWarnings="true"/>' + '\n'
            f'      <distributionFunction1DPerturber value="cauchy">\n'
            f'        <median value="0.0" />\n'
            f'        <scale value="1.0e-4" />\n'
            f'      </distributionFunction1DPerturber>\n'
            f'      <distributionFunction1DPrior value="uniform">\n'
            f'        <limitLower value="-3.00" />\n'
            f'        <limitUpper value="+0.00" />\n'
            f'      </distributionFunction1DPrior>\n'
            f'      <operatorUnaryMapper value="identity" />\n'
            f'    </modelParameter>\n'
            f'    <modelParameter value="active">\n'
            f'      <name value="haloMassFunctionParameters/exponentRedshiftDetection{cls}" />\n'
            r'      <label value="\beta_\mathrm{det,' + cls + r'}" ignoreWarnings="true"/>' + '\n'
            f'      <distributionFunction1DPerturber value="cauchy">\n'
            f'        <median value="0.0" />\n'
            f'        <scale value="1.0e-4" />\n'
            f'      </distributionFunction1DPerturber>\n'
            f'      <distributionFunction1DPrior value="uniform">\n'
            f'        <limitLower value="-3.00" />\n'
            f'        <limitUpper value="+3.00" />\n'
            f'      </distributionFunction1DPrior>\n'
            f'      <operatorUnaryMapper value="identity" />\n'
            f'    </modelParameter>\n'
        )

    # Finish the closer
    config_closer += (
        '\n'
        '  </posteriorSampleSimulation>\n'
        '\n'
        '  <!-- Random seed -->\n'
        '  <randomNumberGenerator value="GSL">\n'
        '    <seed          value="219" />\n'
        '    <mpiRankOffset value="true"/>\n'
        '  </randomNumberGenerator>\n'
        '\n'
        '</parameters>\n'
    )

    # -------------------------------------------------------- parameters_opener
    parameters_opener = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<parameters>\n'
        '  <formatVersion>2</formatVersion>\n'
        '\n'
        '  <!-- Controllable parameters of the halo mass function -->\n'
        '  <haloMassFunctionParameters value="" ignoreWarnings="true">\n'
        '    <!-- Halo mass function -->\n'
        '    <a                                value="+0.75"/>\n'
        '    <b                                value="+0.84"/>\n'
        '    <c                                value="+2.26"/>\n'
        '    <p                                value="-1.09"/>\n'
        '    <q                                value="+0.51"/>\n'
        '    <normalization                    value="+0.57"/>\n'
        '    <!-- Artificial halo variance -->\n'
        '    <artificialExponentMass           value="-0.50e+0"/>\n'
        '    <artificialExponentGrowthFactor   value="+1.00e+0"/>\n'
        '    <artificialCountParticles         value="+1.00e+2"/>\n'
        '    <!-- ETHOS window function -->\n'
        '    <cW0                              value="+2.59"/>\n'
        '    <beta0                            value="+3.51"/>\n'
        '    <cW1                              value="+0.74"/>\n'
        '    <beta1                            value="+4.83"/>\n'
        '    <wavenumberScaledMinimum          value="+0.00"/>\n'
        '    <powerSpectrumSmoothingWidth      value="+1.00"/>\n'
    )

    # -------------------------------------------------------- parameters_closer
    parameters_closer = ''
    if perturbations:
        parameters_closer += '    <!-- Perturbation model parameters -->\n'
        for lbl in sorted(perturbations):
            parameters_closer += f'    <perturbation{lbl} value="0.0"/>\n'

    if isolation_biases:
        parameters_closer += '    <!-- Isolation bias model parameters -->\n'
        for lbl in sorted(isolation_biases):
            parameters_closer += (
                f'    <isolationBias{lbl}         value="1.0"/>\n'
                f'    <isolationBiasExponent{lbl} value="0.0"/>\n'
            )

    if detection_efficiency_classes:
        parameters_closer += '    <!-- Detection efficiency model parameters -->\n'
    for cls in sorted(detection_efficiency_classes):
        parameters_closer += (
            f'    <massMinimumParticleCount{cls}  value="10.0"/>\n'
            f'    <efficiencyAtMassMinimum{cls}   value=" 1.0"/>\n'
            f'    <exponentMassDetection{cls}     value=" 0.0"/>\n'
            f'    <exponentRedshiftDetection{cls} value=" 0.0"/>\n'
        )

    parameters_closer += (
        '  </haloMassFunctionParameters>\n'
        '\n'
        '  <!-- Parameter controlling the fractional variance due to model discrepancy -->\n'
        '  <varianceFractionalModelDiscrepancy value="0.0"/>\n'
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

    print('Creating parameters...')
    isolation_biases, perturbations = _step_a_group_labels(simulations, options)

    print('Modifying halo mass function XML files...')
    _step_b_suite_hmf_xml(simulations, options)

    print('Counting realizations...')
    count_realizations = _step_c_count_realizations(simulations, options)

    print('Building entry groups...')
    entry_groups = _step_d_group_entries(simulations, options, count_realizations)

    print('Generating base parameter files...')
    config_likelihood, detection_efficiency_classes = _step_e_base_files(
        entry_groups, options
    )

    print('Assembling MCMC config strings...')
    (config_opener, config_initializer, config_resumer,
     config_closer, parameters_opener, parameters_closer) = _step_f_config_strings(
        perturbations, isolation_biases, detection_efficiency_classes, options
    )

    output_dir = options['outputDirectory']

    with open(output_dir + 'haloMassFunctionConfig.xml', 'w') as fh:
        fh.write(config_opener)
        fh.write(config_likelihood)
        fh.write(config_initializer)
        fh.write(config_closer)

    with open(output_dir + 'haloMassFunctionConfigResume.xml', 'w') as fh:
        config_likelihood_resume = re.sub(r'(<appendSamples\s+value=")false("\/>)',r'\1true\2',config_likelihood)
        fh.write(config_opener)
        fh.write(config_likelihood_resume)
        fh.write(config_resumer)
        fh.write(config_closer)

    with open(output_dir + 'haloMassFunctionParameters.xml', 'w') as fh:
        fh.write(parameters_opener)
        fh.write(parameters_closer)

    print('Done.')


if __name__ == '__main__':
    main()
