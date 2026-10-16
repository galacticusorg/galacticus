Dust Compendium Datasets
========================

Datasets containing dust extinctions for simple galactic geometries, computed using the methods described by :cite:t:`benson_compendium_2018`, are available for download. The individual datasets are linked to and described below.

Galacticus reads these files directly, through the :galacticus-class:`dustAttenuationAtlasCompendium` dust attenuation class: set its ``fileName`` to the name of the file to use, and its ``url`` to the download link given below for the record containing it. The file is fetched once and cached under the dynamic datasets path, so it need not be downloaded by hand.

Spheroid radii
--------------

Spheroid sizes in these datasets are tabulated as the **scale radius**, *r*\ :sub:`s`\ , of the spheroid density profile, in units of the radial scale length of the stellar disk. This differs from the convention of :cite:t:`ferrara_atlas_1999`, whose atlas is tabulated instead against an **effective radius**, *R*\ :sub:`e`\ , and the two are related by

.. math::

   r_\mathrm{s} = 1.16 \, R_\mathrm{e}.

That factor is not a property of any one profile, but a correspondence between two of them, and it is worth understanding where it comes from before using either dataset.

:cite:t:`ferrara_atlas_1999` describe their spheroids as having an effective radius, but realize them in the radiative transfer calculation as Jaffe profiles, which stand in for the *R*\ :sup:`1/4` profiles that the effective radius properly belongs to. The appendix of :cite:t:`bianchi_monte_carlo_1996`, which describes the Monte Carlo code the atlas was computed with, sets the correspondence between the two, and their reasoning fixes the number. The *R*\ :sup:`1/4` and Jaffe distributions are similar in shape, but the relation between the Jaffe scale radius and the effective radius depends on which property of the galaxy one chooses to match best. Their Monte Carlo associates a random number directly with the fraction of the luminosity enclosed within a given radius, and draws photon emission positions from it, so the quantity that must correspond between the two profiles is that enclosed luminosity --- not any particular radius, and not the total light. Under that criterion the best correspondence is *r*\ :sub:`s`\ =1.16\ *R*\ :sub:`e`\ ; matching the half-light radii instead would have given roughly 1.35. The choice is therefore a real one, and not a matter of convention: the two differ by about the same 16 percent as the factor itself.

Two consequences follow for anyone using these datasets.

First, a spheroid described by a given effective radius in :cite:t:`ferrara_atlas_1999` corresponds to a **larger** scale radius here, by 16 percent. The two tabulations are not interchangeable on the same axis values.

Second, matching a galaxy whose spheroid follows some third profile --- a Hernquist profile, say --- onto either tabulation is necessarily approximate, and by the same reasoning it should be done on an integral property rather than on a scale radius. Galacticus matches on the half-mass radius, which is the one measure that means the same thing whatever profile either side assumes, and then converts to whichever radius the tabulation is labeled with: this is what the ``spheroidProfile`` parameter of :galacticus-class:`dustAttenuationAtlasCompendium` selects. The residual uncertainty in doing so is of the same order as the 1.16 factor itself, so spheroid attenuations for a profile other than the tabulated one should not be relied upon at better than the ten percent level.

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

Galactic disks follow exponential profiles in both radial and vertical directions, with the vertical scale height equal to 0.1 or 0.5 times the radial scale length as encoded in the file name with prefix ``hzStars``.

Dust Geometry
~~~~~~~~~~~~~

Dust is distributed in the disk, and follows exponential profiles in both radial and vertical directions. The vertical scale height is set to a multiple, *h*\ :sub:`z`\ , of the stellar disk scale height. The value of *h*\ :sub:`z`\  is encoded in each file name with prefix ``hzDust``.

Files
~~~~~

* All files can be found `here <https://doi.org/10.5281/zenodo.6335951>`_.

Ferrara et al. (1999)-like models
---------------------------------

These models are intended to closely match the models run by :cite:t:`ferrara_atlas_1999` - they use the same dust grain properties and galactic geometry.

Files
~~~~~

* `"Original" <https://doi.org/10.5281/zenodo.6336095>`_ - models tabulated at the same inclinations, optical depths, wavelengths, and morphologies as in :cite:t:`ferrara_atlas_1999`.
* `"HiRes" <https://doi.org/10.5281/zenodo.6336097>`_ - models tabulated with a much higher resolution grid of inclinations, optical depths, wavelengths, and morphologies than in :cite:t:`ferrara_atlas_1999`, but with the same dust properties and galactic mass distributions.
