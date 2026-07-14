#!/usr/bin/env python3
"""Top-level orchestration script for the dark matter constraint pipeline.

Python port of constraints/pipelines/darkMatter/pipeline.pl
Andrew Benson (ported to Python 2026)
"""

import argparse
import datetime
import json
import os
import re
import subprocess
import sys
import time

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
    parser.add_argument('--markConverged', default=None, metavar='STAGE',
                        help='Human judgement gate: mark STAGE (e.g. haloMassFunction) as '
                             'converged — scancel its running job and record a durable marker — '
                             'then exit. The next pipeline run advances past it.')
    parser.add_argument('--diagnose', default=None, metavar='STAGE',
                        help='Print DE-appropriate convergence diagnostics for STAGE (via dendros) '
                             'and exit. Aids the human judgement; does not make it.')
    parser.add_argument('--diagnoseSteps', default=1000, type=int,
                        help='Recent-window size (steps per chain) for --diagnose (default 1000; '
                             '0 = full history, slower).')
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


def _write_calibrated_containers(output_dir, params_determined):
    """Write calibrated parameter values into the on-disk container files.

    The cross-stage coupling is file-based: a stage's base files ``xi:include`` a
    container such as ``haloMassFunctionParameters.xml``, and a *downstream* stage
    ``xi:include``s the *same* file to reuse the calibrated model (e.g. the
    progenitor stage reusing the calibrated HMF window function). The generators
    write those containers with prior/default values; this updates them in place
    with the extracted best-fit values so the coupling survives standalone or
    partial runs — a fresh ``xinclude`` of the container then yields calibrated
    values without relying on the in-memory ``params_determined`` handoff.

    Each parameter ``root/leaf/...`` whose container ``{output_dir}{root}.xml``
    exists has its element at that path set to the calibrated value; parameters
    with no matching container file are left to the base-file injection.
    """
    roots = {}
    for name in params_determined:
        roots.setdefault(name.split('/', 1)[0], []).append(name)
    for root, names in sorted(roots.items()):
        container_path = f'{output_dir}{root}.xml'
        if not os.path.exists(container_path):
            continue
        tree = ET.parse(container_path)
        container_root = tree.getroot()
        updated = 0
        for name in names:
            elem = _find_param_element(container_root, name)
            if elem is not None:
                elem.set('value', str(params_determined[name]))
                updated += 1
        ET.indent(container_root)
        tree.write(container_path, xml_declaration=True, encoding='utf-8')
        print(f'    Wrote {updated}/{len(names)} calibrated values -> {os.path.basename(container_path)}')


# ---------------------------------------------------------------------------
# Stage state machine: manifest + SLURM introspection + human convergence gate
# ---------------------------------------------------------------------------
#
# Convergence is a deliberate human judgement (the stopping criterion is set
# unreachably high on purpose), so the driver separates the one human bit from
# the machine-decidable bits. Per stage, on each entry:
#
#   converged marker | chain logs | job active |  action
#   ---------------- | ---------- | ---------- |  -------------------------------
#   yes (human-set)  |     -      |     -      |  extract -> handoff -> postprocess
#   no               |    no      |    no      |  submit FRESH  ({stage}Config.xml)
#   no               |   yes      |    no      |  submit RESUME ({stage}ConfigResume.xml)
#   no               |     -      |   yes      |  in progress -- wait for it
#
# The human flips the marker with `--markConverged {stage}` (which scancels the
# still-running chain job); that unblocks the waiting driver, which then advances.

def _manifest_path(output_dir, stage):
    return f'{output_dir}{stage}.state.json'


def _manifest_read(output_dir, stage):
    """Return the stage manifest dict ({} if none written yet)."""
    path = _manifest_path(output_dir, stage)
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def _manifest_write(output_dir, stage, **updates):
    """Merge *updates* into the stage manifest and persist it."""
    manifest = _manifest_read(output_dir, stage)
    manifest.update(updates)
    with open(_manifest_path(output_dir, stage), 'w') as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write('\n')
    return manifest


def _job_label(stage):
    """The deterministic SLURM job name for a stage (matches submission)."""
    return f'darkMatterPipeline{stage[0].upper()}{stage[1:]}'


def _job_active(stage):
    """True if a queued/running SLURM job for this stage exists (by name)."""
    try:
        result = subprocess.run(
            ['squeue', '--me', '--name', _job_label(stage), '-h', '-o', '%i'],
            capture_output=True, text=True,
        )
    except FileNotFoundError:
        return False
    return result.returncode == 0 and bool(result.stdout.strip())


def _wait_for_job(stage, interval=60):
    """Block until no SLURM job for this stage remains active."""
    while _job_active(stage):
        time.sleep(interval)


def _git_revision():
    try:
        return subprocess.run(
            ['git', '-C', os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'rev-parse', 'HEAD'],
            capture_output=True, text=True,
        ).stdout.strip() or None
    except Exception:
        return None


def _do_mark_converged(output_dir, stage):
    """Record the human 'converged' decision for a stage and stop its job."""
    if _job_active(stage):
        print(f"Cancelling active '{stage}' job ({_job_label(stage)})...")
        subprocess.run(['scancel', '--name', _job_label(stage)], check=False)
    _manifest_write(
        output_dir, stage,
        converged=True,
        convergedAt=datetime.datetime.now().isoformat(timespec='seconds'),
        gitRevision=_git_revision(),
    )
    print(f"Marked '{stage}' converged (manifest: {_manifest_path(output_dir, stage)}).")
    print(f"The next `pipeline.py` run will extract the maximum, write the calibrated "
          f"handoff, postprocess, and advance to the next stage.")


def _do_diagnose(output_dir, stage, steps):
    """Emit DE-appropriate convergence aids for a stage via dendros, then exit.

    Deliberately uses the effect-size ensemble drift + ESS/autocorrelation, and
    *not* per-walker Gelman-Rubin / Geweke, which are inflated for an interacting
    differentialEvolution ensemble (see plan review Section 3).
    """
    from dendros import open_mcmc

    config = f'{output_dir}{stage}Config.xml'
    if not os.path.exists(config):
        raise FileNotFoundError(f"No config for stage '{stage}': {config}")
    window = steps if steps and steps > 0 else None
    run = open_mcmc(config, max_steps=window)
    n_chains = len(run.chains)
    n_steps = run.chains[0].n_steps if n_chains else 0
    outliers = run.outlier_chains() if n_chains else ()
    manifest = _manifest_read(output_dir, stage)

    print(f"[diagnose {stage}]  chains={n_chains}  steps/chain kept={n_steps}"
          f"  ({'last '+str(window) if window else 'full history'})  outliers={len(outliers)}")
    print(f"  converged marker: {manifest.get('converged', False)}"
          f"{'  (set '+manifest['convergedAt']+')' if manifest.get('convergedAt') else ''}")
    if n_chains == 0 or n_steps < 4:
        print("  (too few samples yet for diagnostics)")
        return
    try:
        drift = run.ensemble_drift(drop_chains=outliers)
        print(f"  ensemble-mean drift (effect size, recent window halves):"
              f"  max={drift.max_drift():.3f} sigma   worst={drift.worst_parameter()}")
        print(f"    -> a well-mixed tail sits well below ~0.1 sigma once the window spans")
        print(f"       several autocorrelation times; widen --diagnoseSteps if it's short.")
    except Exception as e:
        print(f"  ensemble_drift unavailable: {e}")
    print("  (effect-size drift is the DE-appropriate signal; per-walker R-hat/Geweke are")
    print("   omitted — inflated for a DE ensemble — see plan Section 3. Slower autocorrelation")
    print("   diagnostics: dendros effective_sample_size / acceptance_rate on this config.)")
    print(f"  When you judge it done:  {sys.argv[0]} --markConverged {stage}")


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


def _build_job(label, task, output_dir, galacticus, config_file, callgrind):
    """Build the queueManager job dict for a stage's MCMC, running *config_file*."""
    n_proc     = task['ppn'] * task['nodes']
    omp_prefix = 'export OMP_NUM_THREADS=1; '
    mpi_prefix = f'mpirun --oversubscribe --n {n_proc} '
    callgrind_prefix = (
        f'--output-filename {output_dir}{label}.vlog valgrind --tool=callgrind '
        if callgrind == 'yes' else ''
    )
    return {
        'command':    f'{omp_prefix} {mpi_prefix}{callgrind_prefix}{galacticus} {config_file}',
        'launchFile': f'{output_dir}{label}.sh',
        'logFile':    f'{output_dir}{label}.log',
        'label':      _job_label(label),
        'ppn':        task['ppn'],
        'nodes':      task['nodes'],
        'walltime':   '7-00:00:00',
    }


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

    # Subcommands act on an existing run directory and then exit.
    if options.get('diagnose'):
        _do_diagnose(output_dir, options['diagnose'], options['diagnoseSteps'])
        return
    if options.get('markConverged'):
        _do_mark_converged(output_dir, options['markConverged'])
        return

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

    for task_index, task in enumerate(tasks):
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

        # Step 6 (stage state machine): the stage is "done" only when a human has
        # marked it converged (not merely because a chain log exists). Interlock on
        # upstream stages, then submit fresh/resume and wait, resubmitting the resume
        # config on each walltime kill, until the marker is set (which scancels the
        # running job and unblocks us).
        for prior in tasks[:task_index]:
            if not _manifest_read(output_dir, prior['label']).get('converged'):
                print(f"  Refusing to start '{label}': upstream stage "
                      f"'{prior['label']}' is not marked converged.")
                print(f"    Inspect:  {sys.argv[0]} --diagnose {prior['label']}")
                print(f"    Accept:   {sys.argv[0]} --markConverged {prior['label']}")
                sys.exit(1)

        chain_log_0   = log_file_root(config) + '_0000.log'
        resume_config = output_dir + f'{label}ConfigResume.xml'
        resubmits     = 0
        while not _manifest_read(output_dir, label).get('converged'):
            if _job_active(label):
                print(f"  A '{label}' job is already active; waiting for it to end...")
                _wait_for_job(label)
                continue
            logs_exist    = os.path.exists(chain_log_0)
            submit_config = resume_config if logs_exist else config_file
            kind          = 'resume' if logs_exist else 'fresh'
            print(f"  Running {kind} MCMC for '{label}' "
                  f"({os.path.basename(submit_config)})  [{datetime.datetime.now()}]")
            job     = _build_job(label, task, output_dir, galacticus, submit_config,
                                 options['callgrind'])
            manager = queueManager.factory(args)
            t_start = time.time()
            submit_jobs(manager, [job])
            elapsed = time.time() - t_start
            if _manifest_read(output_dir, label).get('converged'):
                break                       # human marked it done (scancel unblocked us)
            resubmits += 1
            if elapsed < 120:
                print(f"  '{label}' job ended after {elapsed:.0f}s but is not marked "
                      f"converged — this looks like a crash, not a walltime kill. Stopping.")
                print(f"    Check {output_dir}{label}.log, fix, and re-run.")
                sys.exit(1)
            print(f"  '{label}' job ended (walltime?) and is not marked converged; "
                  f"resubmitting resume [#{resubmits}].")
            print(f"    When you judge it done:  {sys.argv[0]} --markConverged {label}")
        print(f"  '{label}' is marked converged  [{datetime.datetime.now()}]")

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

        # Step 8: re-apply and re-write with the newly extracted parameters, and
        # persist the calibrated container files for the file-based cross-stage
        # handoff (so a downstream stage's xi:include picks up calibrated values).
        if options['writeParameters'] == 'yes':
            print('  Updating parameter files...')
            _apply_parameters(config, params_determined)
            _write_parameters(config)
            _write_calibrated_containers(output_dir, params_determined)
            print('  ...done')

        # Step 9: postprocess.
        if task['postprocess']:
            print('  Postprocessing...')
            subprocess.run(postprocess_cmd, check=True)
            print('  ...done')


if __name__ == '__main__':
    main()
