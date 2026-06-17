#!/usr/bin/env python3
"""Top-level orchestration script for the dark matter constraint pipeline.

Python port of constraints/pipelines/darkMatter/pipeline.pl
Andrew Benson (ported to Python 2026)
"""

import argparse
import datetime
import os
import re
import subprocess
import sys

import lxml.etree as ET
import numpy as np

import queueManager
from queueManager import submit_jobs
from List.ExtraUtils import as_array
from XML.Utils import xml_to_dict
from Galacticus.Constraints.Parameters import (
    log_file_root,
    maximum_likelihood_parameter_vector,
    maximum_posterior_parameter_vector,
    parameter_names,
)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    _default_pipeline = os.path.join(
        os.environ.get('GALACTICUS_EXEC_PATH', ''),
        'constraints', 'pipelines', 'darkMatter', ''
    )
    parser = argparse.ArgumentParser(
        description='Dark matter constraint pipeline orchestration script.'
    )
    parser.add_argument('--outputDirectory',    default=os.path.join(os.getcwd(), 'pipeline'),
                        help='Output directory (default: ./pipeline)')
    parser.add_argument('--updateResults',      default='yes',
                        help='Extract best-fit parameters after MCMC (default: yes)')
    parser.add_argument('--maximum',            default='posterior',
                        choices=['posterior', 'likelihood'],
                        help='Whether to extract posterior or likelihood maximum (default: posterior)')
    parser.add_argument('--writeParameters',    default='yes',
                        help='Write updated parameter files to disk (default: yes)')
    parser.add_argument('--callgrind',          default='no',
                        help='Run Galacticus under callgrind (default: no)')
    parser.add_argument('--generateContent',    default='yes',
                        help='Run GenerateContent script before MCMC (default: yes)')
    parser.add_argument('--haloMassFunctionNodes', default=16, type=int,
                        help='Number of nodes for the haloMassFunction MCMC job (default: 16)')
    parser.add_argument('--haloMassFunctionPPN',   default=32, type=int,
                        help='Processors per node for the haloMassFunction MCMC job (default: 32)')
    parser.add_argument('--progenitorMassFunctionNodes', default=16, type=int,
                        help='Number of nodes for the progenitorMassFunction MCMC job (default: 16)')
    parser.add_argument('--progenitorMassFunctionPPN',   default=32, type=int,
                        help='Processors per node for the progenitorMassFunction MCMC job (default: 32)')
    parser.add_argument('--countParticlesMinimum', default=300, type=int,
                        help='Minimum particles per halo (default: 300); forwarded to sub-scripts')
    parser.add_argument('--force',              default='yes',
                        help='Force generation of plots, even if they already exist (default: yes)')
    parser.add_argument('--select',             default=None, action='append',
                        help='Simulation selection filter; may be repeated; forwarded to sub-scripts')
    parser.add_argument('--partition',          default=None,
                        help='SLURM partition override; forwarded to PostProcess')
    parser.add_argument('--jobMaximum',         default=None, type=int,
                        help='Maximum concurrent jobs; forwarded to PostProcess')
    parser.add_argument('--waitOnSubmit',       default=5, type=int,
                        help='Time (in seconds) to wait after submitting; forwarded to PostProcess (default: 5)')
    parser.add_argument('--waitOnActive',       default=60, type=int,
                        help='Time (in seconds) to wait between polling active jobs; forwarded to PostProcess (default: 60)')
    parser.add_argument('--pipelinePath',       default=_default_pipeline,
                        help='Pipeline directory (default: $GALACTICUS_EXEC_PATH/constraints/pipelines/darkMatter/)')
    parser.add_argument('--initializeToPosteriorMaximum', default=None,
                        help='Log file root for initializing MCMC from prior posterior maximum; '
                             'forwarded to GenerateContent')
    args, args_extra = parser.parse_known_args()
    for attr in ('outputDirectory', 'pipelinePath'):
        val = getattr(args, attr)
        if val and not val.endswith('/'):
            setattr(args, attr, val + '/')
    return args, args_extra


# ---------------------------------------------------------------------------
# Config XML helpers
# ---------------------------------------------------------------------------

def _parse_config(config_file_name):
    """Parse a Galacticus config XML file into a nested dict (no XInclude)."""
    return xml_to_dict(ET.parse(config_file_name).getroot())


def _innermost_likelihood_model(model):
    """Walk the nested posteriorSampleLikelihood chain to the innermost model."""
    while 'posteriorSampleLikelihood' in model:
        model = model['posteriorSampleLikelihood']
    return model


def _load_base_parameter_trees(config):
    """Load per-redshift base parameter XML files as lxml element trees.

    Walks each outer posteriorSampleLikelihood down to the innermost
    model, derives filenames (substituting the z-suffix in
    baseParametersFileName if needed), expands XInclude, and stores a
    list of (path, root_element) tuples in inner['_trees']. Mutates the
    config dict in-place.
    """
    outer  = config['posteriorSampleLikelihood']
    models = as_array(outer.get('posteriorSampleLikelihood', []))
    for model in models:
        inner     = _innermost_likelihood_model(model)
        trees = []
        if 'redshifts' in inner:
            redshifts = inner['redshifts']['value'].split()
            pattern   = inner['baseParametersFileName']['value']
            for z in redshifts:
                path = re.sub(r'_z\d+\.\d+', f'_z{z}', pattern)
                tree = ET.parse(path)
                tree.xinclude()
                trees.append((path, tree.getroot()))
        else:
            path = inner['baseParametersFileName']['value']
            tree = ET.parse(path)
            tree.xinclude()
            trees.append((path, tree.getroot()))
        inner['_trees'] = trees


def _find_param_element(root_elem, param_name):
    """Return the lxml element reached by following the '/'-separated param_name path.

    Each path component is the tag of the next child element. Returns None if
    any component is absent.
    """
    elem = root_elem
    for part in param_name.split('/'):
        elem = elem.find(part)
        if elem is None:
            return None
    return elem


def _apply_parameters(config, params_determined):
    """Inject best-fit parameter values into the in-memory base parameter trees.

    Respects parameterMap entries on the outer posteriorSampleLikelihood
    to restrict which active parameters are applied to each inner
    likelihood model.
    """
    if not params_determined:
        return
    simulation = config['posteriorSampleSimulation']
    outer      = config['posteriorSampleLikelihood']
    models     = as_array(outer.get('posteriorSampleLikelihood', []))
    param_maps = as_array(outer.get('parameterMap', []))
    parameters = as_array(simulation.get('modelParameter', []))
    parameter_names = [parameter['name']['value'] for parameter in parameters]
    for i, model in enumerate(models):
        allowed = (param_maps[i]['value'].split() if i < len(param_maps) else [])
        inner   = _innermost_likelihood_model(model)
        for path, elem_root in inner.get('_trees', []):
            for param_name, value in params_determined.items():
                if allowed and param_name not in allowed and param_name in parameter_names:
                    continue
                elem = _find_param_element(elem_root, param_name)
                if elem is not None:
                    elem.set('value', str(value))


def _write_parameters(config):
    """Write each in-memory base parameter element tree back to its source file."""
    outer  = config['posteriorSampleLikelihood']
    models = as_array(outer.get('posteriorSampleLikelihood', []))
    for model in models:
        inner = _innermost_likelihood_model(model)
        for path, elem_root in inner.get('_trees', []):
            ET.indent(elem_root)
            ET.ElementTree(elem_root).write(
                path, xml_declaration=True, encoding='utf-8'
            )


# ---------------------------------------------------------------------------
# Subprocess helper
# ---------------------------------------------------------------------------

def _call_subscript(script_path, flag_pairs):
    """Invoke a Python sub-script as a subprocess with the given flag pairs.

    flag_pairs is a list of (flag_name, value) tuples. None values are omitted;
    list values produce one --flag value per element.
    """
    cmd = [sys.executable, script_path]
    for key, val in flag_pairs:
        if val is None:
            continue
        if isinstance(val, list):
            for v in val:
                cmd += [f'--{key}', str(v)]
        else:
            cmd += [f'--{key}', str(val)]
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args, args_extra = _parse_args()
    options          = vars(args)
    pipeline_path    = options['pipelinePath']
    output_dir       = options['outputDirectory']
    galacticus       = os.path.join(
        os.environ.get('GALACTICUS_EXEC_PATH', ''), 'Galacticus.exe'
    )

    os.makedirs(output_dir, exist_ok=True)

    with open(output_dir + 'command.txt', 'w') as fh:
        fh.write(' '.join(sys.argv) + '\n')

    tasks = [
        {
            'label':       'haloMassFunction',
            'ppn':         options['haloMassFunctionPPN'],
            'nodes':       options['haloMassFunctionNodes'],
            'postprocess': True,
        },
        {
            'label':       'progenitorMassFunction',
            'ppn':         options['progenitorMassFunctionPPN'],
            'nodes':       options['progenitorMassFunctionNodes'],
            'postprocess': True,
        }, 
   ]

    params_determined = {}

    for task in tasks:
        label       = task['label']
        config_file = output_dir + f'{label}Config.xml'
        print(f'Begin task: {label}')

        # Step 1: generate config and base parameter XML files.
        if options['generateContent'] == 'yes':
            print('  Generating content...')
            subscript_options = [
                    ('pipelinePath',                  pipeline_path),
                    ('outputDirectory',               output_dir),
                    ('countParticlesMinimum',         options['countParticlesMinimum']),
                    ('select',                        options['select']),
                    ('initializeToPosteriorMaximum',  options.get('initializeToPosteriorMaximum')),
                ]
            for i in range(0, len(args_extra), 2):
                m = re.match(r'\-\-'+label+r'_([a-zA-Z]+)',args_extra[i])
                if m:
                    subscript_options.append((m.group(1), args_extra[i+1]))
            _call_subscript(
                pipeline_path + f'{label}GenerateContent.py',
                subscript_options
            )
            print('  ...done')

        # Step 2: record postprocess command to file.
        postprocess_cmd = [
            sys.executable,
            pipeline_path + f'{label}PostProcess.py',
            '--pipelinePath',    pipeline_path,
            '--outputDirectory', output_dir,
        ]
        for key in ('countParticlesMinimum', 'partition', 'jobMaximum', 'waitOnActive', 'waitOnSubmit', 'force'):
            if options.get(key) is not None:
                postprocess_cmd += [f'--{key}', str(options[key])]
        if options.get('select'):
            for s in options['select']:
                postprocess_cmd += ['--select', s]
        with open(output_dir + 'postProcessCommand.txt', 'w') as fh:
            fh.write(' '.join(postprocess_cmd) + '\n')

        # Step 3: parse config; load base parameter element trees.
        print('  Processing base parameter files...')
        config = _parse_config(config_file)
        _load_base_parameter_trees(config)

        # Step 4: inject parameters determined by earlier tasks.
        _apply_parameters(config, params_determined)
        print('  ...done')

        # Step 5: write base parameter files to disk.
        if options['writeParameters'] == 'yes':
            print('  Writing updated parameter files...')
            _write_parameters(config)
            print('  ...done')

        # Step 6: submit the MCMC job (skip if chain log already exists).
        chain_log_0 = log_file_root(config) + '_0000.log'
        if not os.path.exists(chain_log_0):
            n_proc = task["ppn"]*task['nodes']
            omp_prefix = (
                f'export OMP_NUM_THREADS=1; '
                )
            mpi_prefix = (
                f'mpirun --oversubscribe --n {n_proc} '
                )
            callgrind_prefix = (
                f'--output-filename {output_dir}{label}.vlog '
                f'valgrind --tool=callgrind '
                if options['callgrind'] == 'yes' else ''
            )
            job = {
                'command':    f'{omp_prefix} {mpi_prefix}{callgrind_prefix}{galacticus} {config_file}',
                'launchFile': f'{output_dir}{label}.sh',
                'logFile':    f'{output_dir}{label}.log',
                'label':      f'darkMatterPipeline{label[0].upper()}{label[1:]}',
                'ppn':        task['ppn'],
                'nodes':      task['nodes'],
                'walltime':   '7-00:00:00',
            }
            print(f"  Running MCMC for '{label}'  [{datetime.datetime.now()}]")
            manager = queueManager.factory(args)
            submit_jobs(manager, [job])
            print(f'  ...done [{datetime.datetime.now()}]')

        # Step 7: extract best-fit parameters from chain logs.
        results_file = output_dir + 'results.txt'
        if options['updateResults'] == 'yes':
            print(f"  Extracting maximum {options['maximum']} parameter vector...")
            if options['maximum'] == 'likelihood':
                params_vec, _ = maximum_likelihood_parameter_vector(config)
            else:
                params_vec, _ = maximum_posterior_parameter_vector(config)
            names = parameter_names(config)
            for i, name in enumerate(names):
                params_determined[name] = float(params_vec[i])
            print('  ...done')

            print('  Outputting maximum likelihood parameters...')
            with open(results_file, 'w') as fh:
                for name in sorted(params_determined):
                    fh.write(f'{name}\t{params_determined[name]}\n')
            print('  ...done')
        else:
            if not os.path.exists(results_file):
                raise FileNotFoundError(f'No results.txt found at {results_file}')
            print('  Reading prior maximum likelihood parameters...')
            with open(results_file) as fh:
                for line in fh:
                    m = re.match(r'([a-zA-Z0-9:/]+)\s+([\+\-\d\.e]+)', line)
                    if not m:
                        raise ValueError(f'Cannot parse results.txt line: {line!r}')
                    params_determined[m.group(1)] = float(m.group(2))
            print('  ...done')

        # Step 8: re-apply and re-write with the newly extracted parameters.
        if options['writeParameters'] == 'yes':
            print('  Updating parameter files...')
            _apply_parameters(config, params_determined)
            _write_parameters(config)
            print('  ...done')

        # Step 9: postprocess.
        if task['postprocess']:
            print('  Postprocessing...')
            subprocess.run(postprocess_cmd, check=True)
            print('  ...done')


if __name__ == '__main__':
    main()
