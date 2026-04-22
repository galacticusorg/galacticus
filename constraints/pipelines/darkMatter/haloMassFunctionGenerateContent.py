#!/usr/bin/env python3
# Script to generate content for halo mass function constraint pipeline.
# Python port of constraints/pipelines/darkMatter/haloMassFunctionGenerateContent.pl
# Andrew Benson (ported to Python 2026)

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
# XML helpers
# ---------------------------------------------------------------------------

def _remove_hmf_nodes(root, removal_value):
    """Splice out <haloMassFunction value=removal_value> nodes from the chain.

    Mirrors the Perl linked-list traversal: the removed node is discarded along
    with any non-haloMassFunction children it carries; only its inner
    <haloMassFunction> child is promoted into its place.
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
    parser.add_argument('--select',            default=None,
                        help='Simulation selection filter (suite::group::resolution::...)')
    parser.add_argument('--initializeToPosteriorMaximum', default=None,
                        help='Log file root for initializing from a prior posterior maximum')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Steps A & B
# ---------------------------------------------------------------------------

def _step_a_group_labels(simulations, options):
    """Step A: build isolation-bias and perturbation label dicts (stop_after='group').

    Returns (isolation_biases, perturbations) and mutates each entry['group']
    to carry 'labelIsolationBias', 'isolationBias', 'labelPerturbation',
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
            group['labelPerturbation'] = label
            group['perturbation']      = f'haloMassFunctionParameters/perturbation{label}'
            perturbations[label]       = perturbations.get(label, 0) + 1
        else:
            group['labelPerturbation'] = ''
            group['perturbation']      = ''

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

    print('Generating base parameter files...')

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
                         + f"haloMassFunction_{suite['name']}_{grp['name']}_{res['name']}"
                         + f"_{sim['name']}_{real}_z{redshift}.hdf5")
            file_ref  = (data_path_static
                         + f"haloMassFunction_{matched['suite']}_{grp['name']}_{res['name']}"
                         + f"_{matched['simulation']}_{real}_z{redshift}.hdf5")
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
            f"%DATASTATICPATH%/darkMatter/haloMassFunction_{suite['name']}_{grp['name']}"
            + f"_{res['name']}_{sim['name']}_{real}_z{z}.hdf5"
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
                f'  <version>0.9.4</version>\n'
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
                lbl = grp['labelPerturbation']
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
