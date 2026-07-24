# Dark-matter constraint pipeline

`pipeline.py` orchestrates the calibration of dark-matter merger-tree model
parameters against N-body simulation measurements. It runs a sequence of MCMC
stages, each calibrating one class of quantity, and passes the calibrated result
of each stage forward to the next.

The stages, in order, are:

1. **`haloMassFunction`** — calibrates the halo mass function model (including the
   `ETHOSExtended` window function / `cosmologicalMassVariance`) against the N-body
   halo mass functions.
2. **`progenitorMassFunction`** — calibrates the merger-tree branching model against
   the N-body progenitor mass functions, **reusing the calibrated window function
   from stage 1**.

Each stage is a differential-evolution MCMC run under SLURM, using Galacticus as
the likelihood.

---

## Quick start

```bash
# from $GALACTICUS_EXEC_PATH
constraints/pipelines/darkMatter/pipeline.py \
    --outputDirectory /path/to/run/ \
    --select 'MDPL::MDPL2::resolutionX1::CDM'
```

This generates the parameter files for the first stage, submits its MCMC to the
queue, and **blocks** while it runs. You then drive it to completion with the two
interaction commands below.

The `--select` value filters which simulations/resolutions/dark-matter types are
included; it is a `::`-separated path into `simulations.xml`
(`suite::group::resolution::transferFunction::…`) and may be repeated.

---

## How a run proceeds — the stage state machine

Convergence of these chains is **a deliberate human judgement**, not an automatic
threshold: the stopping criterion is set unreachably high on purpose, so a chain
runs until *you* inspect it and decide it is done. The pipeline therefore separates
the one human decision from everything it *can* decide mechanically.

On each entry to a stage the driver looks at three things — a durable **converged
marker**, whether **chain logs exist**, and whether a **job is active** in SLURM —
and acts:

| converged marker | chain logs | job active | action |
|---|---|---|---|
| **yes** (you set it) | – | – | extract best fit → write calibrated handoff → post-process → advance to next stage |
| no | no | no | submit a **fresh** run (`{stage}Config.xml`) |
| no | yes | no | submit a **resume** run (`{stage}ConfigResume.xml`) — e.g. after a wall-time kill |
| no | – | yes | a job is already running — wait for it |

The driver keeps (re)submitting and waiting until the stage is marked converged.
A downstream stage **refuses to start** until every upstream stage is marked
converged (this is what protects the progenitor stage from calibrating against an
unfinished halo-mass-function fit).

### Interacting with a running stage

Open a second terminal (the first is blocked waiting on the job).

**Inspect convergence:**

```bash
pipeline.py --outputDirectory /path/to/run/ --diagnose haloMassFunction
```

This reads the chains (via [`dendros`](https://github.com/galacticusorg/dendros))
and prints a *differential-evolution-appropriate* mixing diagnostic: the
standardized **ensemble-mean drift** over a recent window (`--diagnoseSteps`,
default 1000 steps/chain) and which parameter is moving most. A well-mixed tail
sits well below ~0.1 σ once the window spans several autocorrelation times.

> Per-walker Gelman-Rubin R̂ and Geweke z-scores are intentionally **not** used
> here: for an interacting DE ensemble they are inflated by the long
> autocorrelation time and over-report non-convergence. See the developer guide.

**Accept it (the human gate):**

```bash
pipeline.py --outputDirectory /path/to/run/ --markConverged haloMassFunction
```

This cancels the still-running chain job, writes a durable `converged` marker (with
a timestamp and git revision) to `haloMassFunction.state.json`, and exits.
Cancelling the job unblocks the waiting driver in the first terminal, which then
extracts the best fit, writes the calibrated hand-off, post-processes, and advances
to the progenitor stage.

If instead the wall-timer kills a job before you mark it converged, the driver
resubmits the **resume** config and keeps going. (A job that dies in under two
minutes is treated as a crash, not a wall-time kill, and the driver stops so you
can investigate.)

### Blocking vs. detached

By default the driver **blocks**: it waits on the queue for the lifetime of each
MCMC (up to the 7-day wall-time). That single long-lived process must survive for
the whole calibration.

With `--detached` the driver instead performs **one action per invocation and
exits**: it submits the next job (fresh, or resume if chain logs exist), records the
job id in the manifest, and returns; a converged-but-unprocessed stage is instead
advanced (extract → hand-off → post-process). You then re-invoke the *same* command
— by hand, from `cron`, or via the `/loop` skill — to take the next step. Nothing
long-lived has to stay alive, which suits HPC. `--diagnose` and `--markConverged`
work identically in both modes. (This relies on the file-based hand-off: each
detached invocation is a fresh process, so the calibrated container on disk — not
in-memory state — is what carries the coupling forward.)

```bash
# submit / advance one step, then exit
pipeline.py --outputDirectory /path/to/run/ --select '…' --detached yes
```

---

## The cross-stage hand-off

Both stages `xi:include` the same container file, `haloMassFunctionParameters.xml`
(written with prior/default values by the stage-1 generator). When stage 1 is
accepted, the driver writes the **calibrated** values back into that container in
place (`_write_calibrated_containers`). The progenitor stage's parameter files then
pick up the calibrated window function directly from disk, so the coupling is
file-based and survives partial or standalone runs — it does not rely on values
held only in memory.

---

## Key files in the output directory

| file | what it is |
|---|---|
| `{stage}Config.xml` / `{stage}ConfigResume.xml` | the MCMC configuration (fresh / resume) |
| `{stage}Base_*.xml` | per-selection Galacticus model files (the likelihood) |
| `{stage}Parameters.xml` | the stage's controllable-parameter container |
| `haloMassFunctionParameters.xml` | the shared HMF container (calibrated in place after stage 1) |
| `{stage}Chains_NNNN.log` | per-rank MCMC chain logs |
| `{stage}Gamma.log` | differential-evolution step-size history |
| `{stage}.state.json` | stage manifest: the `converged` marker + provenance |
| `results.txt` | accumulated best-fit parameter values |
| `command.txt` | the exact command used to launch the run |

---

## Commonly used options

| option | meaning |
|---|---|
| `--outputDirectory DIR` | where everything is written (default `./pipeline`) |
| `--select SEL` | simulation selection filter (repeatable) |
| `--diagnose STAGE` | print convergence diagnostics for a stage and exit |
| `--diagnoseSteps N` | recent-window size for `--diagnose` (default 1000; 0 = full history) |
| `--markConverged STAGE` | mark a stage converged (cancels its job) and exit |
| `--detached yes` | perform one state-machine action and exit (re-invoke to advance) instead of blocking |
| `--{stage}Nodes N`, `--{stage}PPN N` | SLURM node / processors-per-node for a stage's job |
| `--generateContent no` | skip regenerating parameter files (reuse what is on disk) |
| `--maximum likelihood` | extract the maximum-likelihood (not maximum-posterior) point |
| `--{stage}_{param} VALUE` | forward `--{param} VALUE` to that stage's `GenerateContent` script |

The last form is how per-stage generator options are set, e.g.
`--progenitorMassFunction_massRatioDistribution normal` selects the Gaussian (rather
than the default Hinkley) mass-ratio distribution for the progenitor stage.

---

## Notes

- **MPI launch.** The MCMC jobs run Galacticus under `mpirun`. If you evaluate a
  single generated model file by hand with an MPI-built `Galacticus.exe` on a SLURM
  node, launch it as `mpirun -np 1 ./Galacticus.exe model.xml` — a bare
  `./Galacticus.exe` aborts in `MPI_Init` (Open MPI tries to use SLURM's PMI).
  MPI output files gain a `:MPI####.hdf5` suffix.
- **Merger-tree tolerance.** The progenitor stage builds a wide range of tree
  masses; `mergerTreeEvolution.xml` uses a relaxed `cole2000`
  `toleranceTimeEarliest` (`1.0e-4`) to avoid first-run "branch is making no
  progress" aborts during adaptive tabulation.
