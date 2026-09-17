Dust Absorption and Emission
============================

Dust absorbs the ultraviolet and optical light of a galaxy---from its stars, from the gas they ionize, and from an
active nucleus---and re-radiates that energy in the infrared. This tutorial runs a small model which follows both
sides of that exchange, and analyzes its output in a Jupyter notebook.

The model is described by the parameter file ``parameters/tutorials/dustEmission.xml``. It builds a single merger tree
for a halo of :math:`10^{12}\,M_\odot`, and outputs its galaxy at redshift one, when it is actively forming stars and
hosts a supermassive black hole. To run it,
from the directory where you installed Galacticus:

.. code-block:: bash

   $ ./Galacticus.exe parameters/tutorials/dustEmission.xml

which writes ``dustEmissionTutorial.hdf5``. The first run takes some time---tens of minutes---because Galacticus
tabulates the spectra of stellar populations, and the luminosities of their emission lines, at the output time and
spectral resolution of the model. Those tables are cached (in ``$GALACTICUS_DYNAMIC_DATA_PATH``), so later runs with the
same settings take only seconds. The parameter file outputs:

* spectra of the light of the disk and spheroid stars, and of the AGN, and the luminosities of emission lines from
  HII regions and from the narrow-line region of the AGN;
* for each of those, attenuated by the dust model of :galacticus-class:`dustAttenuationCharlotFall2000`, the light
  transmitted and the light absorbed by each phase of dust: birth clouds around young stars, and the diffuse
  interstellar medium;
* the spectrum of thermal emission from the dust (:galacticus-class:`nodePropertyExtractorSEDDustEmission`), in which
  warm dust in birth clouds emits as a modified blackbody at a fixed temperature, and dust in the diffuse medium emits
  as PAHs (:galacticus-class:`dustEmissionSpectrumRichieHensley2026`), heated by the spectrum of the light they absorb,
  and larger grains, whose temperature is set by energy balance.

The accompanying notebook, ``tutorials/models/01-dust-absorption-and-emission.ipynb``, reads that output with the
`dendros <https://github.com/galacticusorg/dendros>`_ package and plots the spectrum of the galaxy before and after
dust, the light each phase of dust absorbs, and the infrared emission it produces, and checks that the dust emits
exactly the energy it absorbs. It can be viewed, with its output, in the :doc:`notebook tutorials
<../../../tutorials/index>`.
