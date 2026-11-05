"""The emulated variant of the progenitor mass function stage of the dark-matter constraint pipeline.

With ``--progenitorMassFunctionEmulated yes`` (in the detached driver), the progenitor mass function stage samples the
posterior of the merger-tree branching parameters with emulated likelihoods (see the "Emulator-Assisted Calibration"
chapter of the user guide) in place of running the model at every step of every chain. It starts from the content the
direct stage generates (its MCMC configuration and one base parameter file per simulation component), and passes through
these sub-states, taking one action per invocation of the driver:

``generate``
    Resolve each component's base parameter file (with its changes, and ``mergerTreeSeeds`` set to ``random``, so that
    runs with nearby seeds build independent forests), record a hash of them all (the hand-off from upstream stages), and
    write the profile configuration.
``profileStar``, ``profileLine``
    Fit the parameters of the model of N-body halo mass errors, which change only the analyses and which the emulators
    can not resolve, by direct runs (``scripts/emulation/emulatorProfile.py``): a star design at one seed, then a line
    search at several seeds. Their runs are made by stage workers.
``profileAccept``
    Wait for the user to accept the profile (``--markProfiled progenitorMassFunction``), or to recenter it on its optimum
    and repeat it (``--reprofile progenitorMassFunction``).
``fix``
    Write each component's base parameter file with the profiled parameters fixed.
``design``
    Generate a design over the branching parameters (the ``emulatorDesign`` task, with a seed per point), and, for each
    component, a campaign over it, with training configurations for emulators on the whole design and on each half.
``campaign``
    Run the campaigns and train the emulators with stage workers (``scripts/emulation/emulatorStageWorker.py``).
``mcmc``
    Sample the posterior three times: with the emulators trained on the whole design, and on each half.
``systematic``
    Compare the posteriors from the halves with that from the whole design: the largest shift of each parameter's
    median is the systematic error of the emulators (``scripts/emulation/emulatorChainCompare.py``).
``ready``
    Wait for the user to judge the posterior converged (``--markConverged progenitorMassFunction``), after which the
    pipeline extracts the best fit (with the profiled parameters) and hands it on, as for the direct stage.

Whenever the hand-off changes (the hash of the resolved base parameter files differs from that recorded), the emulated
stage is invalidated: its files are moved aside, and it restarts from ``generate``.

Andrew Benson, Claude (28-October-2026).
"""

import datetime
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

import h5py
import lxml.etree as ET
import numpy as np

STAGE = 'progenitorMassFunction'
STATES = ('generate', 'profileStar', 'profileLine', 'profileAccept', 'fix', 'design', 'campaign', 'mcmc', 'systematic', 'ready')
VARIANTS = ('full', 'half0', 'half1')
# Parameters emulated (the branching parameters), and the parameter of the likelihood (the fractional root variance of
# counts). Every other active parameter of the direct stage is profiled.
PREFIX_EMULATED = 'mergerTreeBranchingProbability/'
NAME_ROOT_VARIANCE = 'progenitorMassFunctionParameters/rootVarianceTargetFractional'
# The observables: every progenitor mass function analysis, emulated in its logarithm (as in the calibration of
# 2026-10, which validated these settings), and compared to counts.
OBSERVABLE = ('label="progenitorMassFunction.*" pattern="true" transform="log10" floor="-13.6" rootVarianceFloored="1.0" '
              'likelihoodRangeOnly="true"')
RESUBMISSIONS_MAXIMUM = 3


def add_arguments(parser):
    """Add the options of the emulated stage to the pipeline's argument parser."""
    group = parser.add_argument_group('emulated progenitor mass function stage')
    group.add_argument('--progenitorMassFunctionEmulated', default='no', choices=['no', 'yes'],
                       help='sample the progenitor mass function stage with emulated likelihoods (detached driver only)')
    group.add_argument('--progenitorMassFunctionProfile', default='yes', choices=['no', 'yes'],
                       help='fit the parameters of the N-body halo mass error model by direct runs before emulating (default: yes)')
    group.add_argument('--progenitorMassFunctionProfileModelVariance', default='no', choices=['no', 'yes'],
                       help='include the variance of the model (from its finite number of trees) in the likelihoods of the profile, '
                            'reducing their bias at few seeds (default: no)')
    group.add_argument('--progenitorMassFunctionProfileStepFraction', default=0.03, type=float,
                       help='the step of each profiled parameter, as a fraction of the width of its prior (in its mapped space)')
    group.add_argument('--progenitorMassFunctionProfileSteps', default=None,
                       help='steps of individual profiled parameters, overriding the fraction: `name=step,...` (names without their container)')
    group.add_argument('--progenitorMassFunctionEmulatedChanges', default=None,
                       help='a change file applied to the base parameter file of every component (e.g. to reduce the number of trees for a test)')
    group.add_argument('--progenitorMassFunctionDesignPoints', default=256, type=int, help='the number of points of the design')
    group.add_argument('--progenitorMassFunctionSeedBase', default=1000, type=int, help='the seed of the first run of the design')
    group.add_argument('--progenitorMassFunctionWorkers', default=4, type=int, help='the number of stage worker jobs')
    group.add_argument('--progenitorMassFunctionWorkerCPUs', default=128, type=int, help='the CPUs (runs at once) of each worker')
    group.add_argument('--progenitorMassFunctionWorkerWalltime', default='24:00:00', help='the walltime of each worker job')
    group.add_argument('--progenitorMassFunctionWorkerMemory', default=None, type=int, help='the memory of each worker job, in MB')
    group.add_argument('--progenitorMassFunctionTrainingProcesses', default=8, type=int, help='the processes of each training')
    group.add_argument('--progenitorMassFunctionChains', default=32, type=int, help='the chains (MPI processes) of each emulated MCMC')
    group.add_argument('--progenitorMassFunctionSteps', default=20000, type=int, help='the steps of each chain of each emulated MCMC')
    group.add_argument('--progenitorMassFunctionLocal', default='no', choices=['no', 'yes'],
                       help='run workers and MCMCs locally (blocking), rather than submitting them to Slurm (for testing)')
    group.add_argument('--markProfiled', default=None, metavar='STAGE', help='accept the profile of an emulated stage, then exit')
    group.add_argument('--reprofile', default=None, metavar='STAGE',
                       help='recenter the profile of an emulated stage on its optimum and repeat it, then exit')


class EmulatedStage:
    """The emulated progenitor mass function stage of a pipeline run, with its state."""

    def __init__(self, output_dir, options, galacticus):
        self.output_dir = output_dir
        self.options = options
        self.galacticus = galacticus
        self.directory = os.path.join(output_dir, f'{STAGE}Emulated')
        self.state_file = os.path.join(self.directory, 'state.json')
        self.exec_path = os.environ.get('GALACTICUS_EXEC_PATH', '.')
        sys.path.insert(0, os.path.join(self.exec_path, 'scripts', 'emulation'))

    # State.

    def state(self):
        if not os.path.exists(self.state_file):
            return {'state': 'generate'}
        with open(self.state_file) as handle:
            return json.load(handle)

    def update(self, **updates):
        os.makedirs(self.directory, exist_ok=True)
        state = self.state()
        state.update(updates)
        temporary = self.state_file + '.tmp'
        with open(temporary, 'w') as handle:
            json.dump(state, handle, indent=1)
        os.replace(temporary, self.state_file)
        return state

    def path(self, *names):
        return os.path.join(self.directory, *names)

    def local(self):
        return self.options['progenitorMassFunctionLocal'] == 'yes'

    # The direct stage's content.

    def direct_config(self):
        return os.path.join(self.output_dir, f'{STAGE}Config.xml')

    def components(self):
        """The components of the direct stage: (name, base parameter file, change files), one per likelihood."""
        root = ET.parse(self.direct_config()).getroot()
        components = []
        for element in root.find('posteriorSampleLikelihood').findall('posteriorSampleLikelihood'):
            base = element.find('baseParametersFileName').get('value')
            changes = element.find('changeParametersFileNames')
            changes = [] if changes is None else changes.get('value').split()
            name = os.path.splitext(os.path.basename(base))[0].replace(f'{STAGE}Base_', '')
            components.append((name, base, changes))
        return components

    def active_parameters(self):
        """The active `modelParameter` elements of the direct stage, by name."""
        root = ET.parse(self.direct_config()).getroot()
        return {element.find('name').get('value'): element
                for element in root.find('posteriorSampleSimulation').findall('modelParameter') if element.get('value') == 'active'}

    def profiled_names(self):
        return [name for name in self.active_parameters() if not name.startswith(PREFIX_EMULATED) and name != NAME_ROOT_VARIANCE]

    def emulated_names(self):
        return [name for name in self.active_parameters() if name.startswith(PREFIX_EMULATED)]

    def resolve_bases(self):
        """Write each component's resolved base parameter file (with its changes, and random merger tree seeds), and return
        the hash of them all."""
        from Galacticus.Parameters.resolve import resolve_file
        os.makedirs(self.path('bases'), exist_ok=True)
        seeds = self.path('mergerTreeSeedsRandom.xml')
        with open(seeds, 'w') as handle:
            handle.write('<changes>\n  <change type="replaceOrAppend" path="mergerTreeSeeds">\n    <mergerTreeSeeds value="random"/>\n'
                         '  </change>\n</changes>\n')
        digest = hashlib.sha256()
        for name, base, changes in self.components():
            output = self.path('bases', f'{name}.xml')
            self.resolve_base(base, changes, output)
            with open(output, 'rb') as handle:
                digest.update(handle.read())
        return digest.hexdigest()

    def resolve_base(self, base, changes, output):
        """Resolve one component's base parameter file, with its changes, random merger tree seeds, and any extra changes,
        and give it a random number generator if it has none (each run of a design sets its seed)."""
        from Galacticus.Parameters.resolve import resolve_file
        tree = resolve_file(base, change_files=list(changes) + [self.path('mergerTreeSeedsRandom.xml')] + self.extra_changes())
        root = tree.getroot()
        if root.find('randomNumberGenerator') is None:
            generator = ET.SubElement(root, 'randomNumberGenerator', value='GSL')
            ET.SubElement(generator, 'seed', value='1')
        ET.indent(root)
        tree.write(output, xml_declaration=True, encoding='utf-8')

    def extra_changes(self):
        changes = self.options.get('progenitorMassFunctionEmulatedChanges')
        return [os.path.abspath(changes)] if changes else []

    def handoff_changed(self):
        """Whether the hand-off has changed since the stage was generated (recomputing the hash into a scratch copy)."""
        recorded = self.state().get('handoff')
        if recorded is None or not os.path.exists(self.direct_config()):
            return False
        digest = hashlib.sha256()
        scratch = self.path('handoffCheck.xml')
        for name, base, changes in self.components():
            self.resolve_base(base, changes, scratch)
            with open(scratch, 'rb') as handle:
                digest.update(handle.read())
        os.remove(scratch)
        return digest.hexdigest() != recorded

    def invalidate(self):
        """Move the emulated stage aside, so that it restarts from `generate`."""
        stamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        aside = f'{self.directory}.invalidated.{stamp}'
        os.rename(self.directory, aside)
        print(f'[{STAGE}] the hand-off from upstream stages has changed: the emulated stage was moved to {aside}, and restarts.')

    # Profile.

    def write_profile(self, centers=None):
        """Write the profile configuration over the profiled parameters, centered on their values in the resolved base files
        (or on ``centers``)."""
        from Galacticus.Parameters.resolve import load, resolve_path
        active = self.active_parameters()
        steps = {}
        if self.options.get('progenitorMassFunctionProfileSteps'):
            for entry in self.options['progenitorMassFunctionProfileSteps'].split(','):
                name, value = entry.split('=')
                steps[name.strip()] = float(value)
        components = self.components()
        tree = load(self.path('bases', f'{components[0][0]}.xml'))
        lines = ['<emulatorProfile>', f'  <directory value="{self.path("profile")}"/>',
                 f'  <executable value="{self.executable()}"/>',
                 f'  <rootVarianceTargetFractional value="profile" minimum="{_limits(active[NAME_ROOT_VARIANCE])[0]}" '
                 f'maximum="{_limits(active[NAME_ROOT_VARIANCE])[1]}"/>',
                 f'  <modelVariance value="{"true" if self.options["progenitorMassFunctionProfileModelVariance"] == "yes" else "false"}"/>']
        for name in self.profiled_names():
            element = active[name]
            lower, upper = _limits(element)
            # Step in the logarithm for a prior uniform (or normal) in the logarithm, or a parameter sampled in its logarithm.
            logarithmic = (element.find('distributionFunction1DPrior').get('value') in ('logUniform', 'logNormal')
                           or element.find('operatorUnaryMapper').get('value') == 'logarithm')
            width = math.log(upper / lower) if logarithmic else upper - lower
            if centers is not None:
                center = centers[name]
            else:
                nodes = resolve_path(tree.getroot(), name, tree.getroot())
                if len(nodes) != 1:
                    raise RuntimeError(f"parameter '{name}' is not (uniquely) in the base parameter files")
                center = float(nodes[0].get('value'))
            step = steps.get(name.split('/')[-1], self.options['progenitorMassFunctionProfileStepFraction'] * width)
            lines.append(f'  <parameter name="{name}" center="{center!r}" step="{step!r}" scale="{"log" if logarithmic else "linear"}"/>')
        for name, _, _ in components:
            lines += [f'  <component name="{name}" baseParameters="{self.path("bases", name + ".xml")}">',
                      f'    <observable {OBSERVABLE} likelihoodForm="poissonCounts"/>', '  </component>']
        lines.append('</emulatorProfile>')
        with open(self.path('profile.xml'), 'w') as handle:
            handle.write('\n'.join(lines) + '\n')

    def profile(self):
        import emulatorProfile
        return emulatorProfile.read_configuration(self.path('profile.xml'))

    # Running stages of runs, and MCMCs.

    def executable(self):
        """The copy of the executable used by every run of the stage (so that a relink of the original does not disturb
        them)."""
        return self.path('Galacticus.exe')

    def run_stage(self, stage_file, label):
        """Run (or check on) the workers of a stage. Returns True once all of its work has ended."""
        from Galacticus.Emulation.stage import Stage
        script = os.path.join(self.exec_path, 'scripts', 'emulation', 'emulatorStageWorker.py')
        state = self.state()
        if not self.local() and _job_active(state.get('jobs', {}).get(label)):
            print(f'[{STAGE}] workers ({label}) are running.')
            return False
        stage = Stage.load(stage_file)
        summary = stage.plan()
        for component, index, message in summary['failed']:
            print(f'[{STAGE}]   {component} run {index} failed: {message}')
        if summary['runs'] == 0 and not stage.check():
            return True
        submissions = state.get('submissions', {}).get(label, 0)
        if submissions >= RESUBMISSIONS_MAXIMUM:
            raise RuntimeError(f'the workers of {label} have been submitted {submissions} times without finishing the work '
                               f'(see {stage.directory}/state and the worker logs); fix the problem, then clear `submissions` '
                               f'in {self.state_file}')
        self.update(submissions={**state.get('submissions', {}), label: submissions + 1})
        if self.local():
            print(f"[{STAGE}] running a local worker on {label} ({summary['runs']} runs, {summary['trainings']} trainings)...")
            subprocess.run([sys.executable, script, 'work', stage_file, '--slots', str(self.options['progenitorMassFunctionWorkerCPUs'])],
                           check=True)
            return False
        command = [sys.executable, script, 'submit', stage_file, '--label', label,
                   '--workers', str(self.options['progenitorMassFunctionWorkers']),
                   '--cpus', str(self.options['progenitorMassFunctionWorkerCPUs']),
                   '--walltime', self.options['progenitorMassFunctionWorkerWalltime']]
        if self.options.get('progenitorMassFunctionWorkerMemory'):
            command += ['--memory', str(self.options['progenitorMassFunctionWorkerMemory'])]
        if self.options.get('partition'):
            command += ['--partition', self.options['partition']]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout, end='')
        match = re.search(r'as array job (\S+)', result.stdout)
        if match is None and 'nothing to do' in result.stdout:
            return False
        if match is None:
            raise RuntimeError(f'the job ID of the submitted workers ({label}) was not reported')
        self.update(jobs={**self.state().get('jobs', {}), label: match.group(1)})
        return False

    def profile_stage_file(self, kind):
        """The stage file of the runs of a stage of the profile (one component per profile component, with no training)."""
        from Galacticus.Emulation.stage import Stage, StageComponent
        path = self.path('profile', f'{kind}Stage.json')
        if not os.path.exists(path):
            profile = self.profile()
            Stage.create(path, [StageComponent(name=component.name, manifest=profile.manifest(kind, component))
                                for component in profile.components])
        return path

    # Design and campaign.

    def write_design(self):
        """Generate the design over the emulated parameters, and a campaign, training configurations and a stage file for each
        component."""
        from Galacticus.Emulation.campaign import Campaign
        from Galacticus.Emulation.design import read_design, split_design, write_design
        from Galacticus.Emulation.stage import Stage, StageComponent, Training
        os.makedirs(self.path('design'), exist_ok=True)
        active = self.active_parameters()
        blocks = ''.join('    ' + ET.tostring(active[name], encoding='unicode').strip() + '\n' for name in self.emulated_names())
        task = self.path('design', 'designTask.xml')
        with open(task, 'w') as handle:
            handle.write(f'''<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <formatVersion>2</formatVersion>
  <task value="emulatorDesign">
    <designFileName       value="{self.path('design', 'design.hdf5')}"/>
    <changeFilesRoot      value="{self.path('design', 'changes', 'point')}"/>
    <outputFileNameRoot   value="{self.path('design', 'models', 'point')}"/>
    <seedPerPoint         value="true"/>
    <seedBase             value="{self.options['progenitorMassFunctionSeedBase']}"/>
{blocks}    <posteriorSamples value="sobol">
      <countSamples value="{self.options['progenitorMassFunctionDesignPoints']}"/>
    </posteriorSamples>
  </task>
</parameters>
''')
        with open(self.path('design', 'designTask.log'), 'w') as log:
            status = subprocess.run([self.executable(), task], stdout=log, stderr=subprocess.STDOUT).returncode
        if status != 0:
            raise RuntimeError(f"the emulatorDesign task failed (see {self.path('design', 'designTask.log')})")
        design = read_design(self.path('design', 'design.hdf5'))
        names = design.design.names
        components = []
        for name, _, _ in self.components():
            directory = self.path('campaign', name)
            for subdirectory in ('changes', 'models', 'logs'):
                os.makedirs(os.path.join(directory, subdirectory), exist_ok=True)
            runs = dict(design.runs)
            runs['changeFileName'], runs['outputFileName'] = [], []
            for index, (point, seed) in enumerate(zip(design.runs['pointIndex'], design.runs['seed'])):
                change = os.path.join(directory, 'changes', f'point{index:04d}.xml')
                output = os.path.join(directory, 'models', f'point{index:04d}.hdf5')
                with open(change, 'w') as handle:
                    handle.write('<changes>\n')
                    for k, parameter in enumerate(names):
                        handle.write(f'  <change type="update" path="{parameter}" value="{float(design.design.values[point, k])!r}"/>\n')
                    handle.write(f'  <change type="replaceOrAppend" path="outputFileName">\n    <outputFileName value="{output}"/>\n  </change>\n')
                    handle.write(f'  <change type="replaceOrAppend" path="randomNumberGenerator/seed">\n    <seed value="{int(seed)}"/>\n  </change>\n')
                    handle.write('</changes>\n')
                runs['changeFileName'].append(change)
                runs['outputFileName'].append(output)
            content = design.__class__(design=design.design, runs=runs, attributes=design.attributes)
            write_design(os.path.join(directory, 'design.hdf5'), content, overwrite=True)
            for part, half in enumerate(split_design(content)):
                write_design(os.path.join(directory, f'half{part}.hdf5'), half, overwrite=True)
            manifest = os.path.join(directory, 'campaign.json')
            Campaign.create(manifest, os.path.join(directory, 'design.hdf5'), self.path('fixed', f'{name}.xml'),
                            executable=self.executable(), working_directory=directory, log_directory=os.path.join(directory, 'logs'),
                            overwrite=True)
            trainings = []
            for variant, design_file in (('full', 'design.hdf5'), ('half0', 'half0.hdf5'), ('half1', 'half1.hdf5')):
                configuration = os.path.join(directory, f'train_{variant}.xml')
                with open(configuration, 'w') as handle:
                    handle.write(f'''<emulatorTrain>
  <designFileName       value="{os.path.join(directory, design_file)}"/>
  <emulatorFileName     value="{os.path.join(directory, f'emulator_{variant}.hdf5')}"/>
  <foldsCrossValidation value="{5 if variant == 'full' else 0}"/>
  <storeCholesky        value="false"/>
  <observable {OBSERVABLE}/>
</emulatorTrain>
''')
                trainings.append(Training(variant, configuration, os.path.join(directory, f'emulator_{variant}.hdf5')))
            components.append(StageComponent(name=name, manifest=manifest, trainings=trainings))
        Stage.create(self.path('campaign', 'stage.json'), components, overwrite=True,
                     training_processes=self.options['progenitorMassFunctionTrainingProcesses'])

    # MCMC.

    def write_mcmc(self, variant):
        """Write the emulated MCMC configuration for the emulators of one variant (full, half0, half1)."""
        root = ET.parse(self.direct_config()).getroot()
        emulators = []
        for name, _, _ in self.components():
            path = self.path('campaign', name, f'emulator_{variant}.hdf5')
            with h5py.File(path, 'r') as file:
                emulators += [(path, label) for label in sorted(file['emulators'].keys())]
        likelihood = root.find('posteriorSampleLikelihood')
        replacement = ET.Element('posteriorSampleLikelihood', value='emulated')
        for path, label in emulators:
            emulator = ET.SubElement(replacement, 'emulator', value='gaussianProcess')
            ET.SubElement(emulator, 'fileName', value=path)
            ET.SubElement(emulator, 'label', value=label)
        ET.SubElement(replacement, 'likelihoodForms', value='poissonCounts')
        ET.SubElement(replacement, 'includeEmulatorVariance', value='false')
        ET.SubElement(replacement, 'rootVarianceTargetFractionalParameter', value=NAME_ROOT_VARIANCE)
        root.replace(likelihood, replacement)
        simulation = root.find('posteriorSampleSimulation')
        profiled = set(self.profiled_names())
        for element in list(simulation.findall('modelParameter')):
            if element.get('value') == 'derived' or element.find('name').get('value') in profiled:
                simulation.remove(element)
        countActive = sum(1 for element in simulation.findall('modelParameter') if element.get('value') == 'active')
        directory = self.path('mcmc')
        os.makedirs(directory, exist_ok=True)
        _set(simulation, 'logFileRoot', os.path.join(directory, f'{variant}Chains'))
        _set(simulation, 'stepsMaximum', str(self.options['progenitorMassFunctionSteps']))
        _set(simulation, 'appendLogs', 'false')
        for element in simulation.findall('posteriorSampleStoppingCriterion'):
            simulation.remove(element)
        _set(root, 'outputFileName', os.path.join(directory, f'{variant}.hdf5'))
        initializer = root.find('posteriorSampleStateInitialize')
        root.replace(initializer, ET.fromstring('<posteriorSampleStateInitialize value="latinHypercube">'
                                                '<maximinTrialCount value="1000"/></posteriorSampleStateInitialize>'))
        proposal = root.find('posteriorSampleDffrntlEvltnProposalSize')
        _set(proposal, 'gammaInitial', f'{2.38 / math.sqrt(countActive):.5e}')
        _set(proposal, 'logFileName', os.path.join(directory, f'{variant}Gamma.log'))
        path = os.path.join(directory, f'{variant}Config.xml')
        ET.indent(root)
        ET.ElementTree(root).write(path, xml_declaration=True, encoding='utf-8')
        return path

    def run_mcmc(self, submit):
        """Run (or check on) the three emulated MCMCs. Returns True once all have ended."""
        chains = self.options['progenitorMassFunctionChains']
        pending = []
        for variant in VARIANTS:
            label = _job_label(variant)
            root = self.path('mcmc', f'{variant}Chains')
            if not self.local() and _job_active(self.state().get('jobs', {}).get(label)):
                pending.append(variant)
                continue
            if self.state().get('mcmcSubmitted', {}).get(variant):
                if not os.path.exists(f'{root}_0000.log'):
                    raise RuntimeError(f"the {variant} emulated MCMC ended without writing chains (see {self.path('mcmc', f'{variant}.log')})")
                continue
            config = self.write_mcmc(variant)
            command = f'export OMP_NUM_THREADS=1; mpirun --oversubscribe --n {chains} {self.executable()} {config}'
            self.update(mcmcSubmitted={**self.state().get('mcmcSubmitted', {}), variant: True})
            if self.local():
                print(f'[{STAGE}] running the {variant} emulated MCMC locally ({chains} chains)...')
                with open(self.path('mcmc', f'{variant}.log'), 'w') as log:
                    status = subprocess.run(['bash', '-c', command], stdout=log, stderr=subprocess.STDOUT).returncode
                if status != 0:
                    raise RuntimeError(f"the {variant} emulated MCMC failed (see {self.path('mcmc', f'{variant}.log')})")
            else:
                job_id = submit({'command': command, 'launchFile': self.path('mcmc', f'{variant}.sh'),
                                 'logFile': self.path('mcmc', f'{variant}.log'), 'label': label, 'ppn': chains, 'nodes': 1,
                                 'walltime': '2-00:00:00'})
                self.update(jobs={**self.state().get('jobs', {}), label: str(job_id)})
                print(f'[{STAGE}] submitted the {variant} emulated MCMC (job {job_id}).')
                pending.append(variant)
        return not pending

    # The state machine.

    def advance(self, generate, submit):
        """Take one action of the emulated stage. ``generate()`` generates the direct stage's content; ``submit(job)``
        submits a queueManager job without waiting, returning its ID. Returns True once the stage is ready for the user's
        convergence decision."""
        if os.path.isdir(self.directory) and self.handoff_changed():
            self.invalidate()
        state = self.state()['state']
        if state == 'generate':
            if not os.path.exists(self.direct_config()):
                print(f'[{STAGE}] generating the content of the direct stage...')
                generate()
            os.makedirs(self.directory, exist_ok=True)
            shutil.copy2(self.galacticus, self.executable())
            handoff = self.resolve_bases()
            profiling = self.options['progenitorMassFunctionProfile'] == 'yes' and self.profiled_names()
            if profiling:
                self.write_profile()
            self.update(state='profileStar' if profiling else 'fix', handoff=handoff,
                        generatedAt=datetime.datetime.now().isoformat(timespec='seconds'))
            print(f'[{STAGE}] emulated stage generated ({len(self.components())} components); '
                  + ('profiling the error parameters next.' if profiling else 'no parameters to profile.'))
            return False
        if state in ('profileStar', 'profileLine'):
            import emulatorProfile
            kind = 'star' if state == 'profileStar' else 'line'
            profile = self.profile()
            if not os.path.exists(profile.manifest(kind, profile.components[0])):
                profile.prepare(kind)
            if not self.run_stage(self.profile_stage_file(kind), _job_label(f'profile{kind.capitalize()}')):
                return False
            profile = self.profile()
            result = profile.analyze(kind)
            emulatorProfile._report_result(profile, result)
            if kind == 'star':
                self.update(state='profileLine')
                return False
            self.update(state='profileAccept')
            state = 'profileAccept'
        if state == 'profileAccept':
            print(f'[{STAGE}] the profile of the error parameters is done ({self.path("profile", "line", "result.json")}). Either:')
            print(f'  accept it:  --markProfiled {STAGE}')
            print(f'  or recenter on its optimum and repeat it:  --reprofile {STAGE}')
            return False
        if state == 'fix':
            from Galacticus.Parameters.resolve import resolve_file
            os.makedirs(self.path('fixed'), exist_ok=True)
            changes = []
            if os.path.exists(self.path('profile.xml')):
                profile = self.profile()
                result = profile.read_result('line')
                profile.write_fixed(self.path('fixedParameters.xml'), result['optimum'])
                changes = [self.path('fixedParameters.xml')]
                self.update(fixedParameters=result['optimum'])
            for name, _, _ in self.components():
                resolve_file(self.path('bases', f'{name}.xml'), change_files=changes, output=self.path('fixed', f'{name}.xml'))
            self.update(state='design')
            print(f'[{STAGE}] base parameter files written with the profiled parameters fixed.')
            return False
        if state == 'design':
            self.write_design()
            self.update(state='campaign')
            print(f"[{STAGE}] design of {self.options['progenitorMassFunctionDesignPoints']} points written, with a campaign per component.")
            return False
        if state == 'campaign':
            if not self.run_stage(self.path('campaign', 'stage.json'), _job_label('campaign')):
                return False
            self.update(state='mcmc')
            print(f'[{STAGE}] every campaign is complete, and every emulator trained.')
            return False
        if state == 'mcmc':
            if not self.run_mcmc(submit):
                print(f'[{STAGE}] emulated MCMCs are running.')
                return False
            self.update(state='systematic')
            state = 'systematic'
        if state == 'systematic':
            from Galacticus.Emulation.chains import compare_posteriors, read_chains
            roots = [self.path('mcmc', f'{variant}Chains') for variant in VARIANTS]
            rows = compare_posteriors(read_chains(roots[0]), [read_chains(root) for root in roots[1:]])
            with open(self.path('systematic.json'), 'w') as handle:
                json.dump(rows, handle, indent=1)
            self.update(state='ready')
            state = 'ready'
        if state == 'ready':
            self.report()
            print(f'  If the posterior is converged:  --markConverged {STAGE}')
            return True
        raise RuntimeError(f"unknown state '{state}' of the emulated {STAGE} stage")

    def report(self):
        """Print the systematic error of the emulators."""
        with open(self.path('systematic.json')) as handle:
            rows = json.load(handle)
        print(f'[{STAGE}] emulated posterior (whole design), and the systematic error of the emulators (largest shift of the '
              'median between the halves and the whole):')
        width = max(len(row['name']) for row in rows)
        for row in rows:
            print(f"  {row['name']:{width}}  median {row['median']:12.5g}  sigma {row['sigma']:10.3g}  systematic {row['systematic']:10.3g}"
                  f"  ({row['ratio']:.2f} sigma)")

    def config_for_extraction(self):
        """The configuration whose chains give the stage's best fit: that of the emulators trained on the whole design."""
        return self.path('mcmc', 'fullConfig.xml')

    def fixed_parameters(self):
        return self.state().get('fixedParameters', {})

    def mark_profiled(self):
        if self.state()['state'] != 'profileAccept':
            raise RuntimeError(f"the emulated {STAGE} stage is not waiting for its profile to be accepted (state: {self.state()['state']})")
        self.update(state='fix')
        print(f'[{STAGE}] profile accepted; re-invoke the pipeline to continue.')

    def reprofile(self):
        if self.state()['state'] != 'profileAccept':
            raise RuntimeError(f"the emulated {STAGE} stage is not waiting for its profile to be accepted (state: {self.state()['state']})")
        result = self.profile().read_result('line')
        stamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        os.rename(self.path('profile'), self.path(f'profile.{stamp}'))
        os.replace(self.path('profile.xml'), self.path(f'profile.{stamp}.xml'))
        self.write_profile(centers=result['optimum'])
        self.update(state='profileStar', submissions={})
        print(f'[{STAGE}] profile recentered on its optimum (the previous one is in profile.{stamp}); re-invoke the pipeline to repeat it.')


def _limits(element):
    prior = element.find('distributionFunction1DPrior')
    return float(prior.find('limitLower').get('value')), float(prior.find('limitUpper').get('value'))


def _set(parent, tag, value):
    element = parent.find(tag)
    if element is None:
        element = ET.SubElement(parent, tag)
    element.set('value', value)


# Slurm job states (as reported by `sacct`) of jobs which have ended.
JOB_STATES_ENDED = {'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'DEADLINE', 'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'TIMEOUT'}


def _job_label(name):
    return f'darkMatterPipeline{STAGE[0].upper()}{STAGE[1:]}{name[0].upper()}{name[1:]}'


def _job_active(job_id):
    """Whether a submitted Slurm job (None if none was submitted) may still be pending or running. The job is checked by its
    ID, in the accounting records and then in the queue; if neither says that it has ended, it is taken to be active, so that
    nothing (in particular, planning a stage, which clears its work queue) is done while its workers may be running."""
    if job_id is None:
        return False
    try:
        result = subprocess.run(['sacct', '-j', job_id, '-X', '-n', '-P', '-o', 'State'], capture_output=True, text=True, timeout=120)
        states = [line.split()[0] for line in result.stdout.splitlines() if line.strip()] if result.returncode == 0 else []
    except (OSError, subprocess.TimeoutExpired):
        states = []
    if states:
        return any(state not in JOB_STATES_ENDED for state in states)
    # No accounting records: ask the queue, which lists the job (by its ID) until it ends, and soon after forgets it.
    try:
        result = subprocess.run(['squeue', '-j', job_id, '-h', '-o', '%i'], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        result = None
    if result is not None and result.returncode == 0:
        return bool(result.stdout.strip())
    if result is not None and 'Invalid job id' in result.stderr:
        return False
    print(f'[{STAGE}] could not query Slurm about job {job_id}: waiting.')
    return True
