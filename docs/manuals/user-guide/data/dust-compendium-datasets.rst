Dust Compendium Datasets
========================

Datasets containing dust extinctions for simple galactic geometries, computed using the methods described by `Benson (2018) <https://ui.adsabs.harvard.edu/abs/2018RNAAS...2..188B>`_ and in a format compatible with the `DustCompendium <https://github.com/galacticusorg/analysis-python/blob/master/galacticus/dust/dustCompendium.py>`_ class in the `analysis-python <https://github.com/galacticusorg/analysis-python>`_ tools are available for download. The individual datasets are linked to and described below.

Any questions about these datasets should be directed to `Andrew Benson <mailto:abenson@carnegiescience.edu>`_.

Draine (2003) grains
--------------------

Dust grain properties are taken from `Draine (2003) <http://adsabs.harvard.edu/abs/2003ARA%26A..41..241D>`_ for three different values of *R*\ :sub:`V`\ .

Stellar Geometry
~~~~~~~~~~~~~~~~~

Galactic disks follow exponential profiles in the radial direction, and sech2 distributions in the vertical direction. The vertical scale height is set to a multiple, *h*\ :sub:`d`\ , of the stellar disk scale length, and its value is encoded in each file name. Spheroids follow spherical `Hernquist (1990) <http://adsabs.harvard.edu/abs/2003ARA%26A..41..241D>`_ profiles.

Dust Geometry
~~~~~~~~~~~~~

Dust is distributed in the disk, and follows exponential profiles in both radial and vertical directions. The vertical scale height is set to a multiple, *h*\ :sub:`z`\ , of the stellar disk scale length, and its value is encoded in each file name.

Files
~~~~~

* `R_V=3.1 <https://doi.org/10.5281/zenodo.6335021>`_
* `R_V=4.0 <https://doi.org/10.5281/zenodo.6335545>`_
* `R_V=5.5 <https://doi.org/10.5281/zenodo.6335642>`_

Kim, Martin, and Hendry (1994) grains
-------------------------------------

Dust grain properties are taken from `Kim, Martin, and Hendry (1994) <http://adsabs.harvard.edu/abs/1994ApJ...422..164K>`_, specifically their model with *R*\ :sub:`V`\ =3.1 and either a full scattering calculation or Henyey-Greenstein scattering.

Stellar Geometry
~~~~~~~~~~~~~~~~~

Galactic disks follow exponential profiles in the radial direction, and sech2 distributions in the vertical direction. The vertical scale height is set to a multiple, *h*\ :sub:`d`\ , of the stellar disk scale length, and its value is encoded in each file name. Spheroids follow spherical `Hernquist (1990) <http://adsabs.harvard.edu/abs/2003ARA%26A..41..241D>`_ profiles.

Dust Geometry
~~~~~~~~~~~~~

Dust is distributed in the disk, and follows exponential profiles in both radial and vertical directions. The vertical scale height is set to a multiple, *h*\ :sub:`z`\ , of the stellar disk scale length, and its value is encoded in each file name.

Files
~~~~~

* `R_V=3.1 and full scattering calculation <https://doi.org/10.5281/zenodo.6335668>`_
* `R_V=3.1 and Henyey-Greenstein scattering <https://doi.org/10.5281/zenodo.6335670>`_

Grasil-like models
------------------

These models are intended to mimic the geometries used by `Grasil <https://adlibitum.oats.inaf.it/silva/grasil/grasil.html>`_.

Dust grain properties are taken from `Draine (2003) <http://adsabs.harvard.edu/abs/2003ARA%26A..41..241D>`_ - specifically their model with either *R*\ :sub:`V`\ =3.1 or *R*\ :sub:`V`\ =5.5 as encoded in the file name with prefix ``dustD03``.

Stellar Geometry
~~~~~~~~~~~~~~~~~

Galactic disks follow exponential profiles in both radial and vertical directions, with the vertical scale height equal to 0.1 or 0.5 times the radial scale length as encoded in the file name with prefix ``hzStars``. Note that spheroid radii in this work are listed as the scale radius, *r*\ :sub:`s`\ , while Ferrara et al. (1999) listed the corresponding effective radius, *r*\ :sub:`e`\ =\ *r*\ :sub:`s`\ /1.16.

Dust Geometry
~~~~~~~~~~~~~

Dust is distributed in the disk, and follows exponential profiles in both radial and vertical directions. The vertical scale height is set to a multiple, *h*\ :sub:`z`\ , of the stellar disk scale height. The value of *h*\ :sub:`z`\  is encoded in each file name with prefix ``hzDust``.

Files
~~~~~

* All files can be found `here <https://doi.org/10.5281/zenodo.6335951>`_.

Ferrara et al. (1999)-like models
---------------------------------

These models are intended to closely match the models run by `Ferrara et al. (1999) <http://adsabs.harvard.edu/abs/1999ApJS..123..437F>`_ - they use the same dust grain properties and galactic geometry.

Files
~~~~~

* `"Original" <https://doi.org/10.5281/zenodo.6336095>`_ - models tabulated at the same inclinations, optical depths, wavelengths, and morphologies as in `Ferrara et al. (1999) <http://adsabs.harvard.edu/abs/1999ApJS..123..437F>`_.
* `"HiRes" <https://doi.org/10.5281/zenodo.6336097>`_ - models tabulated with a much higher resolution grid of inclinations, optical depths, wavelengths, and morphologies than in `Ferrara et al. (1999) <http://adsabs.harvard.edu/abs/1999ApJS..123..437F>`_, but with the same dust properties and galactic mass distributions.
