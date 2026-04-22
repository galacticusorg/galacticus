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
