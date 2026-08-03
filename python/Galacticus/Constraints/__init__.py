"""Helpers for the Galacticus constraint pipelines.

These modules support the dark-matter constraint pipeline scripts under
``constraints/pipelines/darkMatter/`` (notably ``pipeline.py`` and
``haloMassFunctionPostProcess.py``).

Submodules:

* :mod:`~Galacticus.Constraints.Parameters` -- read MCMC chain output,
  resolve maximum-likelihood / maximum-posterior parameter vectors, and
  expose parameter names and log-file roots.
* :mod:`~Galacticus.Constraints.Simulations` -- parse and iterate over the
  simulation-set XML descriptors that drive the pipeline.
"""
