Computing Broadband Stellar Luminosities
=========================================

Computing, and outputting, broadband luminosities of the stellar continuum can be achieved by adding parameters to your parameter file as in this example:

.. code-block:: xml

   <luminosityFilter   value="RGO_I SDSS_g   SDSS_r   SDSS_i"  />
   <luminosityRedshift value="0.0   0.1      0.1      0.1"     />
   <luminosityType     value="rest  observed observed observed"/>

For each broadband luminosity to be computed we must specify the filter to use (via the ``luminosityFilter`` parameter), the redshift at which the luminosity should be computed (and output; via the ``luminosityRedshift`` parameter), and the "type" of luminosity, either ``rest`` or ``observed`` frame (via the ``luminosityType``) parameter. Each of these parameters can list as many different luminosities as you want (separated by spaces) - they must all have the same number of entries, otherwise an error will be thrown.

Available filters can be found in the `datasets <https://github.com/galacticusorg/datasets>`_ repository `here <https://github.com/galacticusorg/datasets/tree/master/static/filters>`_ - simply use the file name of the filter (with the ``.xml`` suffix removed) in the ``luminosityFilter`` parameter.

The ``luminosityRedshift`` parameter specifies redshift at which each luminosity will be computed and output. Note that if you do not have a corresponding output redshift defined via the ``outputTimes`` parameter, then the luminosity will not be output. If you want to output the same luminosity at every output redshift you can use the special value "``all``" in place of a numerical redshift in this parameter.

The ``luminosityType`` parameter specifies whether the luminosity should be computed in the rest-frame or observed-frame of the galaxy. If you choose ``rest`` then the luminosity will be computed in the galaxy rest-frame. If you choose ``observed`` then the luminosity is computed in the observer frame. Note that, in this case, the calculation is performed by blueshifting the filter to the redshift of the galaxy (rather than redshifting the galaxy spectrum to :math:`z=0`). This means that you must account for the compression of photon frequencies when converting observer-frame luminosities to apparent magnitudes.

Luminosities are output separately for disk and spheroid components. For example, using the above parameters the first luminosity entry would lead to the following two datasets being output:

.. code-block:: text

   diskLuminositiesStellar:RGO_I:rest:z0.000
   spheroidLuminositiesStellar:RGO_I:rest:z0.000

Luminosities are output in units corresponding to the zero-point of the `AB-magnitude system <https://en.wikipedia.org/wiki/AB_magnitude>`_ This allows for convenient conversion into AB magnitudes simply using:

.. math::

   M = -2.5 \log_{10}L

where :math:`L` is the luminosity as output, and :math:`M` is the AB magnitude. Apparent magnitudes can then be computed by simply adding the appropriate distance modulus. Note that, for observed apparent magnitudes you must also add :math:`-2.5 \log_{10} (1+z)` to the absolute magnitude to account for the compression of photon frequencies due to redshifting (which is not included automatically due to the way in which Galacticus computes observed-frame luminosities as described above).
