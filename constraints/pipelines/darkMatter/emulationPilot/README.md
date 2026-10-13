# Emulator pilot: progenitor mass function stage, MDPL2 (M4)

Instructions for running the emulator-assisted calibration pilot of the progenitor mass function stage of the
dark-matter pipeline, on the obshpc cluster. They are written for a Claude Code session on the cluster, working with
Andrew Benson. Background and design decisions are in the plan, `~/Dropbox/AI/galacticus-emulator-calibration-plan.md`
(sections "M4 design pass" and "M4 smoke test"), on the workstation; this file is self-contained, so the cluster does
not need it.

**Scope.** MDPL2 only (`--select 'MDPL::MDPL2::resolutionX1::CDM'`), the only simulation with progenitor mass function
data. 36 analyses (`progenitorMassFunctionM{j}Z{i}`: parent-mass bin j, progenitor redshift i), 388 bins in the
likelihood. 11 sampled parameters: 10 emulator inputs (G0, gamma1-3; errorAmplitude, errorExponent,
errorFractionalHighMass; errorTreeNormalization, errorTreeExponentMass, errorTreeExponentTimeLag) and
`rootVarianceTargetFractional` (C), which changes only the likelihood. There are **no direct-calibration chains** to
compare with: validation is by held-out runs here, and later by direct runs at the emulated posterior.

**Stop and report to Andrew** at each point marked **GATE**. Do not skip ahead.

**Code.** This branch, `pilot/m4-mdpl2`, is a local integration branch, not for merging. It combines four branches
under review: PR #1587 (`poissonCounts`), #1601 (`logNormalFloored`, `pathSamples`, `willEvaluate`), #1592 (the
differential evolution fixes and archive), and `feat/emulation-train-likelihood-range` (`likelihoodRangeOnly` and
parallel training). If any of these change before the pilot, rebuild this branch from them.

---

## 1. Setup

Read `~/Dropbox/AI/memory/galacticus.md` and the obshpc section of `machines.md` (synced to the cluster) before
building. In particular: **append** to `GALACTICUS_FCFLAGS`, never set it; compute nodes have 128 cores and ~1 TB, so
the workstation's `-j8` limit does not apply.

```bash
cd /scratch/abenson
git clone https://github.com/galacticusorg/galacticus.git galacticus_m4   # or reuse a checkout
cd galacticus_m4
git fetch origin pilot/m4-mdpl2 && git switch pilot/m4-mdpl2
export GALACTICUS_EXEC_PATH=`pwd`
pip install -e '.[emulation-gp]'            # scikit-learn, for training
```

Build **two** executables, on a compute node (`srun --pty ...`):

```bash
make -j64 Galacticus.exe                    # serial: the campaign runs (one process each, under Slurm)
make -j64 MPI=yes Galacticus.exe            # gives Galacticus.exe_MPI: the MCMC (one chain per process)
```

The campaign must use the serial build: an MPI build started without `mpirun` aborts in `MPI_Init` under Slurm. Note
`make MPI=yes Galacticus.exe_MPI` silently does nothing; the target is always `Galacticus.exe`. Confirm both binaries'
timestamps moved.

## 2. Generate the stage, with the calibrated halo mass function

The emulators are valid only for the halo mass function hand-off they were trained with, so the stage must use the
**calibrated** halo mass function container from the converged halo mass function run. **Ask Andrew for that run's
output directory** (call it `$HMF_RUN`).

```bash
export PILOT=/scratch/abenson/m4pilot
P=$GALACTICUS_EXEC_PATH/constraints/pipelines/darkMatter/
mkdir -p $PILOT/pipeline
for stage in haloMassFunction progenitorMassFunction; do
  python3 $P/${stage}GenerateContent.py --pipelinePath $P --outputDirectory $PILOT/pipeline/ \
      --select 'MDPL::MDPL2::resolutionX1::CDM'
done
cp $HMF_RUN/haloMassFunctionParameters.xml $PILOT/pipeline/
diff $HMF_RUN/haloMassFunction_MDPL.xml $PILOT/pipeline/haloMassFunction_MDPL.xml
diff $HMF_RUN/haloMassFunctionMappings_MDPL_MDPL2.xml $PILOT/pipeline/haloMassFunctionMappings_MDPL_MDPL2.xml
```

`pipelinePath` must be **absolute**. If either `diff` shows differences, use the converged run's copy, and tell Andrew.
Check that the copied `haloMassFunctionParameters.xml` holds calibrated values (not the generator's defaults).

## 3. Step 0: cost and timing (small)

Write the standalone base model (at the stage's own 3000 trees per decade) and a design:

```bash
D=$GALACTICUS_EXEC_PATH/constraints/pipelines/darkMatter/emulationPilot
python3 $D/makeInputs.py --pipelineOutput $PILOT/pipeline --campaignDirectory $PILOT/wave1 --countSamples 512
cd $PILOT/wave1 && $GALACTICUS_EXEC_PATH/Galacticus.exe design.xml
```

This writes `$PILOT/pipeline/emulatorBase.xml`, and `design.hdf5` with one change file per run in `changes/`. Then
time single full-resolution runs on a compute node, each with `/usr/bin/time -v` and `OMP_NUM_THREADS` set:

1. design point 0 with 32 threads;
2. the same with 128 threads (the scaling decides how many CPUs each run gets);
3. the same with 32 threads, but with the `<outputAnalysis>` element removed from a copy of `emulatorBase.xml` (the
   difference from run 1 is the cost of the analyses):
   ```python
   import lxml.etree as ET
   t = ET.parse('emulatorBase.xml'); r = t.getroot(); r.remove(r.find('outputAnalysis'))
   t.write('emulatorBaseNoAnalyses.xml', xml_declaration=True, encoding='UTF-8')
   ```
   (keep it in `$PILOT/pipeline/`, beside `emulatorBase.xml`).

Run as `Galacticus.exe $PILOT/pipeline/emulatorBase.xml $PILOT/wave1/changes/point0000.xml`, writing the outputs
somewhere other than `models/` (edit a copy of the change file's `outputFileName`), so that the campaign's
bookkeeping is not confused. The smoke test (125 trees) took 462 CPU-s and 108 MB, which scales to ~13 CPU-hours per
run at 3000 trees per decade: confirm or correct this. Check that run 1's output has 36 `analyses/progenitorMassFunction*`
groups, each with `countTarget` and `countPerValue`.

**GATE 1.** Report the wall-clock and CPU time and peak memory of each run, the analyses' share of the cost, and the
implied cost of 512 runs. The decision is whether to generate the training data with a *nested* design (each run
carrying several copies of the analyses with different N-body error parameters, so that one tree build serves several
design points). That needs new tooling on the workstation, and pays off only if the analyses are a small share of the
cost. Otherwise proceed with the plain design below.

## 4. Wave 1 campaign (512 runs)

```bash
cd $PILOT/wave1
C=$GALACTICUS_EXEC_PATH/scripts/emulation/emulatorCalibration.py
python3 $C init campaign.json --design design.hdf5 --base $PILOT/pipeline/emulatorBase.xml \
        --executable $GALACTICUS_EXEC_PATH/Galacticus.exe
python3 $C validate campaign.json --parallel 8
python3 $C submit campaign.json --cpus <from step 0> --walltime <~2x the step 0 wall time> --memory <from step 0>
python3 $C status campaign.json          # repeat; `resubmit` re-runs failed runs
```

`submit` submits one Slurm array job (partition from `galacticusConfig.xml`, or `--partition`). Every run has its own
seed (the design's `seedPerPoint`), so its noise is independent. A run that fails is reported by `status` with the
reason. Look at the first failure before resubmitting anything. A handful of failures is acceptable (training excludes
them), but tell Andrew.

## 5. Train

On one compute node, once the campaign is complete:

```bash
python3 $D/makeConfigs.py --pipelineOutput $PILOT/pipeline --campaignDirectory $PILOT/wave1 \
        --modelFile "$(ls $PILOT/wave1/models/point0000*.hdf5 | head -1)"
python3 $GALACTICUS_EXEC_PATH/scripts/emulation/emulatorTrainReference.py $PILOT/wave1/train.xml \
        --processes 128 > train.log 2>&1
```

`train.xml` has one observable per analysis, in log10, with only the likelihood's bins (`likelihoodRangeOnly`), each
floored at an expected count of 0.01, and 5-fold cross-validation. `--processes` fits the 36 emulators and their 180
cross-validation folds in parallel, one thread each, with results identical to a serial run. A single 10-input fit took
16 s at 512 points on the workstation (single-threaded, synthetic data). Real fits with ~5-15 components each will take
longer.

## 6. Validate the emulators

From `train.log` and the `validation` groups of `emulator.hdf5` (read with
`Galacticus.Emulation.emulatorFile.read`), report for each observable:
* the number of principal components;
* the median and 90th percentile of the held-out error in dex, and its fractional equivalent,
  `ln(10) x dex`;
* the RMS standardized residual (well calibrated: ~1) and the 1 and 2 sigma coverage (nominal 0.68 and 0.95);
* the number of floored training values.

The likelihood's width is set by C, the fractional root variance of the counts (the target counts are up to ~10^7 per
bin, so Poisson errors are ~0.04%). So the emulators' fractional errors must be **well below the fitted C**; compare them
with the C posterior from step 7.

## 7. Emulated MCMC

```bash
cd $PILOT/wave1
mpirun -np 32 $GALACTICUS_EXEC_PATH/Galacticus.exe_MPI mcmc.xml > mcmc.log 2>&1
```

This is differential evolution with one chain per process, the DE-MCz archive (default), and chains started from prior
draws. The 11 parameters are the stage's own `modelParameter` blocks, unchanged. In the smoke test each evaluation took
~0.2 ms. With 512 training points expect a few ms, so 20000 steps on 32 chains take minutes to hours. Chain logs go to
`chains/chain_NNNN.log`: column 5 is the log posterior, column 6 the log likelihood, and columns 7-17 the
parameters, in the order of the `modelParameter` blocks in `mcmc.xml`.

Judge convergence by the drift of the ensemble mean of each parameter over the last few thousand steps, as in the
pipeline README; per-chain Gelman-Rubin is inflated for an interacting DE ensemble. Check that all chains
reach the same region (DE chains can stall in a secondary mode).

**GATE 2.** Report:
* the emulator validation of step 6;
* the maximum a posteriori point and the marginal medians and 16-84% ranges of all 11 parameters (C included);
* whether any parameter sits at a prior bound;
* the emulators' fractional errors relative to the fitted C;
* the costs: CPU-hours of the campaign, wall-clock time of training and of the MCMC.

Write these to `$PILOT/wave1/report.md`, and give Andrew a summary to carry back to the workstation session.
The next steps, after review: a second, focused wave around the wave-1 posterior; direct runs of the full model at the
emulated MAP and at ~20 posterior draws, to compare direct and emulated likelihoods; then integration into
`pipeline.py`.

## Pitfalls already met (smoke test, 2026-10-13)

* The generated base file includes `progenitorMassFunctionParameters.xml` by a path relative to its own directory
  (`makeInputs.py` makes it absolute).
* Don't add component selections to the base: `mergerTreeEvolution.xml` already includes them.
* `posteriorSampleDffrntlEvltnProposalSize` requires `logFileName` (`makeConfigs.py` sets it).
* `fillInZeroBins` extrapolates empty model bins downwards by a factor of 2 per bin (values reach 1e-308). That, and
  `countPerValue` of up to 4e11, are why floors are set per observable as an expected count, not as a fixed value.
* On a fatal error Galacticus may sleep (up to a day) "to allow for attachment of a debugger" instead of exiting. Watch
  logs for `Fatal error`, not just for exit.
