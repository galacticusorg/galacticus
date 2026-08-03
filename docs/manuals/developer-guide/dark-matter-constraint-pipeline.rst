.. _manual-sec-dark-matter-constraint-pipeline:

The Dark-Matter Constraint Pipeline
===================================

The dark-matter constraint pipeline (``constraints/pipelines/darkMatter/``)
calibrates merger-tree model parameters against N-body measurements as a sequence
of MCMC stages. The user-facing workflow is documented in that directory's
``README.md``; this section describes the internals for developers who need to
extend or maintain it.

Architecture
------------

``pipeline.py`` is the orchestrator. It runs a fixed list of *stages* (currently
``haloMassFunction`` then ``progenitorMassFunction``); each stage is a
differential-evolution MCMC using Galacticus as the likelihood. For each stage the
driver:

#. runs ``{stage}GenerateContent.py`` to emit the MCMC configuration and the
   per-selection Galacticus model files;
#. submits the MCMC to SLURM (via ``queueManager.submit_jobs``, which blocks until
   the job ends) and, once the stage is accepted, extracts the best-fit parameters;
#. writes the calibrated result into the container file the next stage includes,
   and post-processes.

The set of simulations, resolutions, redshifts and dark-matter types is described
in ``simulations.xml`` and filtered with the ``--select`` option.

Cross-stage coupling
--------------------

The reason for staging halo-mass-function first is that the progenitor stage must
reuse the **calibrated** window function (``cosmologicalMassVariance`` /
``ETHOSExtended``). The coupling is file-based:

* Both stages' base model files ``xi:include`` the same container,
  ``haloMassFunctionParameters.xml``, which the stage-1 generator writes with
  prior/default values.
* ``_load_base_parameter_trees`` expands XIncludes in place (``lxml``'s
  ``tree.xinclude()``), so the container's parameters are inlined into each base
  tree; ``_apply_parameters`` then overwrites them from the in-memory
  ``params_determined`` dict, and ``_write_parameters`` writes the expanded tree.
* Crucially, after each stage's extraction the driver also calls
  ``_write_calibrated_containers``, which writes the calibrated values back into the
  on-disk container itself. A parameter named ``root/leaf/…`` is written to the
  element at that path in ``{output_dir}{root}.xml`` when that file exists.

Persisting the calibrated container (rather than relying only on the in-memory
hand-off) makes the coupling robust to standalone or partial runs: a fresh
``xinclude`` of ``haloMassFunctionParameters.xml`` yields calibrated values even
when ``params_determined`` is empty.

The stage state machine
-----------------------

Convergence of these chains is a **deliberate human judgement** — the stopping
criterion is set unreachably high (``convergeAfterCount = -1`` →
``posteriorSampleConvergence value="never"``). The driver therefore does *not* try
to detect convergence automatically; it makes the human decision an explicit,
durable artifact and makes everything around it machine-decidable.

Each stage has a JSON manifest, ``{output_dir}{stage}.state.json``, read and written
with ``_manifest_read`` / ``_manifest_write``. The one field that matters for
control flow is ``converged``. On each entry to a stage the driver distinguishes
four situations from the manifest, the presence of chain logs, and a SLURM query
(``_job_active`` matches the deterministic job name ``darkMatterPipeline{Stage}``):

.. list-table::
   :header-rows: 1

   * - ``converged``
     - chain logs
     - job active
     - action
   * - yes
     - –
     - –
     - extract → hand-off → post-process → advance
   * - no
     - no
     - no
     - submit fresh (``{stage}Config.xml``)
   * - no
     - yes
     - no
     - submit resume (``{stage}ConfigResume.xml``)
   * - no
     - –
     - yes
     - wait for the running job

The driver loops — submitting fresh, then resume on each wall-time kill — until the
manifest shows ``converged``. A job that ends in under two minutes without being
marked converged is treated as a crash (not a wall-time kill) and stops the driver.
A stage also refuses to start until every upstream stage's manifest shows
``converged`` (the inter-stage interlock).

The human flips the marker with the ``--markConverged {stage}`` subcommand
(``_do_mark_converged``): it ``scancel``\ s the running job, records
``converged: true`` plus provenance (timestamp, git revision), and exits.
Because ``submit_jobs`` blocks, cancelling the job is also what unblocks the waiting
driver so it can advance.

Blocking vs. detached
~~~~~~~~~~~~~~~~~~~~~~~

The default driver is **blocking**: it waits on ``submit_jobs`` (which polls
``squeue``/``sacct``) for the lifetime of each MCMC. That single process must
survive the whole calibration.

``--detached yes`` selects ``_run_detached``, which performs exactly one action from
the table above per invocation and then exits — submit the next job, or advance a
converged stage — recording the job id in the manifest. It is re-invoked (by a
human, ``cron``, or the ``/loop`` skill) to take each subsequent step, removing the
long-lived process entirely. Two things make this work:

* a **non-blocking submission path** in ``queueManager`` — ``submit_job_detached``
  (and ``SLURMManager.submit_detached``) writes the launch file and ``sbatch``\ es
  the job via the shared ``_submitOne`` helper, then returns the job id without
  polling;
* the **file-based hand-off** — because each detached invocation is a fresh process
  with no in-memory ``params_determined``, the calibrated container on disk is what
  carries the coupling forward. ``_run_detached`` regenerates content only for a
  *fresh* submit (never on a resume, so a resuming job's files are not clobbered).

Convergence diagnostics for a DE ensemble
-----------------------------------------

``--diagnose {stage}`` (``_do_diagnose``) prints a mixing diagnostic via
``dendros``. The choice of statistic matters: Galacticus differential-evolution runs
are an **interacting ensemble** of many walkers with a long integrated
autocorrelation time, and per-walker Gelman-Rubin :math:`\hat{R}` and Geweke
z-scores are badly inflated for such chains — both put an independence-assuming
:math:`s^2/n` variance in the denominator, so they over-report non-convergence even
on a well-mixed run. They are therefore deliberately **not** used as the gate.

Instead the diagnostic reports the standardized **ensemble-mean drift**
(``dendros.ensemble_drift``): the shift of the pooled ensemble mean between an early
and a late window, in units of the posterior width. This is an *effect size*, not a
significance test — appropriate because with many walkers the ensemble mean is
estimated so precisely that a significance test rejects on a negligible drift. The
diagnostic reads only a recent window (``--diagnoseSteps``) for speed. The
autocorrelation-based effective sample size is available from ``dendros`` on demand
but is omitted from the fast path because it is expensive over hundreds of chains.

Extending the pipeline
----------------------

* **A new stage** is a dict appended to the ``tasks`` list in ``main`` plus a
  ``{stage}GenerateContent.py`` and ``{stage}PostProcess.py``. If it consumes an
  upstream calibrated model, have its base files ``xi:include`` the relevant
  container; the interlock and hand-off then apply automatically.
* **A new MCMC parameter** is added in the stage's ``GenerateContent.py`` (an active
  ``modelParameter`` block, a ``parameterMap`` entry, an initializer position, and a
  container entry). The differential-evolution proposal size ``gammaInitial`` should
  be derived from the active-parameter count, not hard-coded.
