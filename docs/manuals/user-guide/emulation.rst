.. _manual-sec-Emulation:

Emulator-Assisted Calibration
=============================

.. note::

   Emulator-assisted calibration is under development. This chapter currently specifies the emulator file format; the
   tasks and tools which create and use emulator files will be described here as they are added.

Calibrating a model by Markov Chain Monte Carlo (see :doc:`tutorials/constraining-parameters`) requires the model to be
run at every step of every chain. When a single model run is expensive this is impractical. Emulator-assisted
calibration instead runs the model at a fixed, space-filling set of points in parameter space (a *design*), trains a
Gaussian-process emulator of the model's *predictions* (for example, a binned stellar mass function) from those runs,
and then samples the posterior using the emulator in place of the model. The likelihood is evaluated exactly, from the
emulated predictions, so it can take any form and can depend on parameters (such as those describing systematic errors
in the data) which the emulator does not need to know about.

Designs are generated in the space of prior *quantiles*: each parameter is described by the cumulative probability of
its prior, :math:`u_i = P_i(<\theta_i)`, which lies in :math:`[0,1]` whatever the prior. The
:galacticus-class:`posteriorSamplesSobol` class generates such a design from a Sobol low-discrepancy sequence. The
emulator's inputs are these quantiles.

.. _manual-sec-EmulatorFileFormat:

Emulator File Format
--------------------

An emulator file is an HDF5 file holding everything needed to evaluate, audit, and reproduce a trained emulator. The
reference implementation of the format, including reading, writing, and evaluation of the emulator, is the Python module
``Galacticus.Emulation.emulatorFile`` (in ``python/Galacticus/Emulation/``).

Conventions
~~~~~~~~~~~

* Array shapes below are given in C (row-major) order, as seen by ``h5py`` and NumPy. Fortran, being column-major, sees
  each array with its dimensions reversed.
* :math:`N` is the number of training points, :math:`d` the number of design parameters, :math:`d_\mathrm{e}` the
  number of inputs to a given emulator, :math:`B` the number of bins in an observable, and :math:`K` the number of
  principal components retained for it.
* All strings, in attributes and in datasets, are fixed-length, null-padded ASCII.
* Stored indices (``pointIndex``, ``folds``) are zero-based.
* Logarithms (``logAmplitude``, ``logLengthScales``) are natural logarithms.
* Group names derived from a count (``prior1``, ``component1``, ...) are numbered from 1.

Root group
~~~~~~~~~~

Attributes:

``format``
   Always ``galacticusEmulator``.
``formatVersion``
   The version of this format, currently ``1``. Readers must reject versions they do not support.
``created``
   The time at which the file was written (ISO 8601).

Further attributes (for example ``creator`` and ``gitHash``) may record provenance.

``design``
~~~~~~~~~~

The parameters varied, their priors, and the points at which the model was run. Attributes record how the design was
generated (for example ``type``, ``seed``, ``randomShift``).

``parameterNames`` (:math:`d`)
   The parameter paths, exactly as given in each ``modelParameter``'s ``name``.
``mappers`` (:math:`d`)
   The ``operatorUnaryMapper`` class of each parameter.
``quantiles`` (:math:`N_\mathrm{design} \times d`)
   The prior quantile, :math:`u_i`, of each parameter at each point of the design.
``values`` (:math:`N_\mathrm{design} \times d`)
   The physical value of each parameter at each point of the design.
``priors/prior{i}``
   The prior of parameter :math:`i`, described by its attributes: ``class`` (the ``distributionFunction1D`` class name,
   e.g. ``logNormal``), plus one attribute for each parameter of that class (e.g. ``x0``, ``sigma``, ``limitLower``,
   ``limitUpper``). These must describe the prior completely: they are compared against the priors of any MCMC which uses
   the emulator, since a change of prior between training and sampling silently changes the meaning of the emulator's
   inputs.

``trainingSets/{label}``
~~~~~~~~~~~~~~~~~~~~~~~~

The training set for one observable, after any transforms, floors, and masking, where ``{label}`` identifies the
observable (normally the label of the ``outputAnalysis`` which produced it). Attributes describe the transform applied
(for example ``transform``, ``floor``, ``rootVarianceFloored``, ``undefined``, ``rootVarianceUndefined``,
``xAxisIsLog``, ``yAxisIsLog``).

``x`` (:math:`B`)
   The bin centers.
``pointIndex`` (:math:`N`)
   The design point from which each row comes. Design points which failed, or were excluded, are simply absent.
``y`` (:math:`N \times B`)
   The transformed model prediction.
``rootVariance`` (:math:`N \times B`)
   The uncertainty in ``y`` due to finite sampling in the model.
``mask`` (:math:`N \times B`)
   Nonzero where a value was floored or replaced.
``yTarget`` (:math:`B`), ``covarianceTarget`` (:math:`B \times B`)
   The data to which the observable is compared, and its covariance.
``yRaw`` (:math:`N \times B`), ``covarianceRaw`` (:math:`N \times B \times B`)
   Optional. The model prediction and its covariance before any transform.

``emulators/{label}``
~~~~~~~~~~~~~~~~~~~~~

The emulator for one observable. Every emulator has a training set of the same label. Attributes:

``kernel``
   The covariance function; currently always ``matern52ARD`` (see below).
``jitter``
   A constant added to the diagonal of each covariance matrix, for numerical stability.
``pcaVarianceRetained``
   The fraction of the variance retained by the principal components.
``countComponents``
   :math:`K`, the number of components (and of ``component{k}`` groups).

Datasets:

``inputNames`` (:math:`d_\mathrm{e}`)
   The design parameters used as inputs, in the emulator's input order. Each must appear in ``design/parameterNames``.
``inputs`` (:math:`N \times d_\mathrm{e}`)
   The training inputs: the design quantiles of those parameters at the training set's points.
``binMean``, ``binScale`` (:math:`B`)
   The standardization of each bin, :math:`z_b = (y_b - \mu_b)/s_b`.
``pcaComponents`` (:math:`K \times B`)
   The retained principal components, :math:`E_{kb}`.
``coefficientMean``, ``coefficientScale`` (:math:`K`)
   The standardization of each principal component coefficient, :math:`\tilde{a}_k = (a_k - \bar{a}_k)/q_k`.

Each ``component{k}`` group holds the Gaussian process for coefficient :math:`k`, with attributes ``logAmplitude``
(:math:`\ln A`) and ``logMarginalLikelihood``, and datasets:

``logLengthScales`` (:math:`d_\mathrm{e}`)
   The length scale, :math:`\ln \ell_i`, for each input.
``targets`` (:math:`N`)
   The standardized coefficients, :math:`\tilde{a}_k`, at the training inputs.
``noiseVariance`` (:math:`N`)
   Their variance, added to the diagonal of the covariance matrix.
``alpha`` (:math:`N`)
   The solution, :math:`\alpha`, of :math:`\mathsf{K} \alpha = \tilde{a}_k`.
``choleskyFactor`` (:math:`N \times N`)
   Optional. The lower-triangular Cholesky factor, :math:`\mathsf{L}`, of :math:`\mathsf{K}`. It is fully determined by the
   other data, and occupies :math:`8 N^2` bytes per component, so it may be omitted and recomputed by the reader.

``validation/{label}``
~~~~~~~~~~~~~~~~~~~~~~

Optional cross-validation results for one observable:

``folds`` (:math:`N`)
   The fold in which each training row was held out.
``heldOutPrediction``, ``heldOutVariance`` (:math:`N \times B`)
   The prediction for each row, and its variance, from the emulator trained without that row's fold.
``rmse``, ``r2``, ``rmsStandardizedResidual``, ``coverage1Sigma``, ``coverage2Sigma`` (:math:`B`)
   Per-bin summary statistics. A well-calibrated emulator has an RMS standardized residual near 1, and 1 and 2 sigma
   coverage near 0.68 and 0.95.

Evaluating an emulator
~~~~~~~~~~~~~~~~~~~~~~

The covariance function is the Matérn :math:`\nu=5/2` kernel with a separate length scale for each input
(automatic relevance determination):

.. math::

   k(\mathbf{u},\mathbf{u}^\prime) = A \left(1 + \sqrt{5} r + \frac{5}{3} r^2\right) \exp\left(-\sqrt{5} r\right),
   \quad r^2 = \sum_i \left(\frac{u_i - u^\prime_i}{\ell_i}\right)^2,

and the covariance matrix of the training data for component :math:`k` is
:math:`\mathsf{K}_{nm} = k(\mathbf{u}_n,\mathbf{u}_m) + \delta_{nm} (\sigma^2_{n} + \epsilon)`, where :math:`\sigma^2_n` is
``noiseVariance`` and :math:`\epsilon` is ``jitter``. Given an input vector :math:`\mathbf{u}_*` (in the order of
``inputNames``) and the vector :math:`\mathbf{k}_*` of its covariances with the training inputs, each component gives a
mean and variance

.. math::

   m_k = \mathbf{k}_*^\mathrm{T} \alpha, \quad v_k = A - \left|\mathsf{L}^{-1} \mathbf{k}_*\right|^2,

which are transformed back to coefficients, :math:`c_k = \bar{a}_k + q_k m_k` and :math:`\sigma^2_k = q_k^2 v_k`, and projected onto the
bins:

.. math::

   y_b = \mu_b + s_b \sum_k c_k E_{kb}, \quad \sigma^2_b = s_b^2 \sum_k E_{kb}^2 \sigma^2_k.

This variance neglects covariance between components and the error due to truncating the principal components; the
cross-validation coverage statistics are the check on whether that matters.
