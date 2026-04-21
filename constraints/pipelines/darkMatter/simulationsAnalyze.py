#!/usr/bin/env python3
# Analyze a variety of cosmological N-body simulations to extract statistics of interest.
# Python port of constraints/pipelines/darkMatter/simulationsAnalyze.pl
# Andrew Benson (ported to Python 2026)

import argparse
import math
import os
import re
import shutil
import sys
from datetime import datetime

import h5py
import lxml.etree as ET
import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))
import queueManager
from Galacticus.Constraints.Simulations import iterate


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------

def _element_to_dict_keyed(el, keyed_tags):
    """Recursively convert an lxml element to a nested dict.

    Elements whose tag is in keyed_tags are converted to dicts keyed by their
    'name' attribute, matching XML::Simple's KeyAttr behaviour. All other
    children with a unique tag become plain dict values; repeating tags become
    lists. An element's own XML attributes are included in the result dict.
    """
    result = dict(el.attrib)
    children_by_tag = {}
    for child in el:
        children_by_tag.setdefault(child.tag, []).append(child)
    for tag, children in children_by_tag.items():
        if tag in keyed_tags:
            result[tag] = {
                child.get('name'): _element_to_dict_keyed(child, keyed_tags)
                for child in children
            }
        elif len(children) == 1:
            result[tag] = _element_to_dict_keyed(children[0], keyed_tags)
        else:
            result[tag] = [_element_to_dict_keyed(c, keyed_tags) for c in children]
    return result


def parse_simulations_xml(path):
    """Parse simulations.xml into a nested dict matching XML::Simple KeyAttr output.

    suite, group, resolution, and simulation elements are keyed by their
    'name' attribute in the result dict.
    """
    keyed_tags = {'suite', 'group', 'resolution', 'simulation'}
    tree = ET.parse(path)
    root = tree.getroot()
    result = {}
    suites = root.findall('suite')
    if suites:
        result['suite'] = {
            el.get('name'): _element_to_dict_keyed(el, keyed_tags)
            for el in suites
        }
    return result


def _parse_param_xml(path):
    """Return the root element of a Galacticus parameter XML file."""
    return ET.parse(path).getroot()


def _write_param_xml(root, path):
    """Serialise a Galacticus parameter XML element tree to file."""
    ET.indent(root)
    tree = ET.ElementTree(root)
    tree.write(path, xml_declaration=True, encoding='unicode')


# ---------------------------------------------------------------------------
# Job submission helpers
# ---------------------------------------------------------------------------

def _translate_job(job):
    """Translate Perl-style job dict keys to Python queueManager keys."""
    j = dict(job)
    # OMP thread count
    if 'ompThreads' in j:
        j['countOpenMPThreads'] = j.pop('ompThreads')
    # Memory: Perl uses 'mem', Python manager uses 'memory'
    if 'mem' in j:
        j['memory'] = j.pop('mem')
    # Log file: map to both output and error
    if 'logFile' in j:
        log = j.pop('logFile')
        j.setdefault('logOutput', log)
        j.setdefault('logError',  log)
    # ppn maps to tasksPerNode
    if 'ppn' in j:
        j['tasksPerNode'] = j.pop('ppn')
    return j


def _submit_jobs(manager, jobs):
    """Submit a list of Perl-style job dicts via the Python queue manager."""
    if jobs:
        manager.submitJobs([_translate_job(j) for j in jobs])


# ---------------------------------------------------------------------------
# Active-step resolution
# ---------------------------------------------------------------------------

def _build_active_steps(active_analyses):
    """Return a dict mapping step IDs to the analyses that require them."""
    active_steps = {}
    if 'haloMassFunction' in active_analyses:
        for step in ('identifyAlwaysIsolated', 'extractHalos', 'massFunctions'):
            active_steps.setdefault(step, []).append('haloMassFunction')
    if 'subhaloStatistics' in active_analyses:
        for step in ('identifyAlwaysIsolated', 'extractSubhalos', 'subhaloFunctions'):
            active_steps.setdefault(step, []).append('subhaloStatistics')
    return active_steps


# ---------------------------------------------------------------------------
# Pre/post-process dispatch
# ---------------------------------------------------------------------------

def _run_hooks(hook_key, step_id, entries, suites_cfg, manager, options):
    """Run all preprocess or postprocess hooks for a step in iteration order.

    Mirrors the Perl while(workDone) loop: iterates through entries calling
    the i-th hook function until no entry has an i-th hook, then submits jobs.
    """
    iteration = 0
    while True:
        jobs = []
        work_done = False
        for entry in entries:
            suite_name = entry['suite']['name']
            hooks = (suites_cfg
                     .get(suite_name, {})
                     .get('steps', {})
                     .get(step_id, {})
                     .get(hook_key, []))
            if iteration < len(hooks):
                work_done = True
                hooks[iteration](entry, jobs, options)
        _submit_jobs(manager, jobs)
        if not work_done:
            break
        iteration += 1


# ---------------------------------------------------------------------------
# Symphony / COZMIC processor stubs  (implemented in later steps)
# ---------------------------------------------------------------------------

def symphony_process_identify_always_isolated(entry, expansion_factor, parameters, options):
    raise NotImplementedError


def symphony_preprocess_extract_halos_locate(entry, jobs, options):
    raise NotImplementedError


def symphony_preprocess_extract_halos_uncontaminated(entry, jobs, options):
    raise NotImplementedError


def symphony_process_extract_halos(entry, expansion_factor, parameters, options):
    raise NotImplementedError


def symphony_process_extract_subhalos(entry, expansion_factor, parameters, options):
    raise NotImplementedError


def symphony_process_subhalo_functions(entry, expansion_factor, parameters, options):
    raise NotImplementedError


def symphony_postprocess_select_in_sphere(entry, jobs, options):
    raise NotImplementedError


def symphony_postprocess_select_in_ics(entry, jobs, options):
    raise NotImplementedError


def symphony_postprocess_analyze(entry, jobs, options):
    raise NotImplementedError


def symphony_postprocess_set_volume(entry, jobs, options):
    raise NotImplementedError


def symphony_postprocess_mass_function(entry, jobs, options):
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Suite configuration
# ---------------------------------------------------------------------------

SUITES = {
    'MDPL': {
        'analyses': ['haloMassFunction'],
        'steps':    {},
    },
    'Symphony': {
        'analyses': ['haloMassFunction', 'subhaloStatistics'],
        'steps': {
            'identifyAlwaysIsolated': {
                'processParameters': symphony_process_identify_always_isolated,
            },
            'extractHalos': {
                'preprocess': [
                    symphony_preprocess_extract_halos_locate,
                    symphony_preprocess_extract_halos_uncontaminated,
                ],
                'processParameters': symphony_process_extract_halos,
                'postprocess': [
                    symphony_postprocess_select_in_sphere,
                    symphony_postprocess_select_in_ics,
                    symphony_postprocess_analyze,
                    symphony_postprocess_set_volume,
                ],
            },
            'extractSubhalos': {
                'preprocess': [
                    symphony_preprocess_extract_halos_locate,
                    symphony_preprocess_extract_halos_uncontaminated,
                ],
                'processParameters': symphony_process_extract_subhalos,
            },
            'massFunctions': {
                'postprocess': [symphony_postprocess_mass_function],
            },
            'subhaloFunctions': {
                'processParameters': symphony_process_subhalo_functions,
            },
        },
    },
}
# COZMIC uses exactly the same processing pipeline as Symphony.
SUITES['COZMIC'] = SUITES['Symphony']


# ---------------------------------------------------------------------------
# Step function stubs  (implemented in later steps)
# ---------------------------------------------------------------------------

def step_identify_always_isolated(entries, suites_cfg, active_steps, manager, options, omp_threads):
    raise NotImplementedError


def step_extract_halos(entries, suites_cfg, active_steps, manager, options, omp_threads):
    raise NotImplementedError


def step_extract_subhalos(entries, suites_cfg, active_steps, manager, options, omp_threads):
    raise NotImplementedError


def step_mass_functions(entries, suites_cfg, active_steps, manager, options, omp_threads):
    raise NotImplementedError


def step_subhalo_functions(entries, suites_cfg, active_steps, manager, options, omp_threads):
    raise NotImplementedError


STEP_FUNCTIONS = {
    'identifyAlwaysIsolated': step_identify_always_isolated,
    'extractHalos':           step_extract_halos,
    'extractSubhalos':        step_extract_subhalos,
    'massFunctions':          step_mass_functions,
    'subhaloFunctions':       step_subhalo_functions,
}

STEP_IDS = [
    'identifyAlwaysIsolated',
    'extractHalos',
    'extractSubhalos',
    'massFunctions',
    'subhaloFunctions',
]


# ---------------------------------------------------------------------------
# Results storage stub  (implemented in later steps)
# ---------------------------------------------------------------------------

def store_results(active_analyses, entries, options):
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze cosmological N-body simulations to extract statistics of interest.'
    )
    parser.add_argument('--simulationDataPath', required=True,
                        help='Path to the simulation data directory')
    parser.add_argument('--pipelinePath', required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--analyses', default='haloMassFunction,subhaloStatistics',
                        help='Comma-separated list of analyses to perform')
    parser.add_argument('--ompThreads', default='max',
                        help='Number of OpenMP threads, or "max" to use the queue config value')
    parser.add_argument('--jobMaximum', type=int, default=64,
                        help='Maximum number of concurrent jobs')
    parser.add_argument('--waitOnSubmit', type=int, default=5,
                        help='Seconds to wait between job submissions')
    parser.add_argument('--waitOnActive', type=int, default=30,
                        help='Seconds to wait when polling active jobs')
    parser.add_argument('--partition', default=None,
                        help='SLURM partition override')
    parser.add_argument('--select', default=None,
                        help='Simulation selection filter (suite::group::resolution::...)')
    parser.add_argument('--reOrder', default='no',
                        help='Re-order simulations by power spectrum class ("yes"/"no")')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = _parse_args()

    # Normalise paths to end with '/'.
    for attr in ('simulationDataPath', 'pipelinePath'):
        val = getattr(args, attr)
        if not val.endswith('/'):
            setattr(args, attr, val + '/')

    # Build the options dict that Simulations.iterate() expects.
    options = vars(args)

    # Determine which analyses and steps are active.
    active_analyses = {a: True for a in args.analyses.split(',')}
    active_steps    = _build_active_steps(active_analyses)

    # Initialise the queue manager.
    manager = queueManager.factory(args)

    # Resolve the number of OMP threads.
    if args.ompThreads == 'max':
        omp_threads = manager.options.get('tasksPerNode', 1)
    else:
        omp_threads = int(args.ompThreads)

    # Parse the simulation definition file.
    simulations = parse_simulations_xml(args.pipelinePath + 'simulations.xml')

    # Build the list of simulation entries to process (up to realization level).
    entries = iterate(simulations, options, stop_after='realization')

    # Execute each active step.
    for step_id in STEP_IDS:
        if step_id not in active_steps:
            continue

        # Preprocess hooks.
        _run_hooks('preprocess', step_id, entries, SUITES, manager, options)

        # Main step function.
        STEP_FUNCTIONS[step_id](entries, SUITES, active_steps, manager, options, omp_threads)

        # Postprocess hooks.
        _run_hooks('postprocess', step_id, entries, SUITES, manager, options)

    # Store results to the datasets repository.
    store_results(active_analyses, entries, options)


if __name__ == '__main__':
    main()
