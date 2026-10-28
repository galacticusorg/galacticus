"""Emulator-assisted calibration of Galacticus.

* :mod:`Galacticus.Emulation.emulatorFile` -- the emulator file format: reading, writing, and the reference
  implementation of emulator evaluation.
* :mod:`Galacticus.Emulation.design` -- reading design files written by the ``emulatorDesign`` task.
* :mod:`Galacticus.Emulation.campaign` -- running and tracking a campaign of model runs over a design.
* :mod:`Galacticus.Emulation.train` -- collecting training sets from a campaign, and training emulators of them.
* :mod:`Galacticus.Emulation.likelihood` -- the likelihood forms of the emulated likelihood, for direct evaluation.
* :mod:`Galacticus.Emulation.profile` -- fitting parameters which change only the analyses of a model by direct runs.
"""
