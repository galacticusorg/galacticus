Input Data
==========

In some configurations, Galacticus requires additional input data to run. For example, if asked to process galaxy formation through a set of externally derived merger trees, then a file describing those trees must be given. In the remainder of this section we describe the structure of external datasets which can be inputs to Galacticus.

Merger Tree Files
-----------------

Galacticus can read merger tree structures from file, and then form and evolve galaxies within those files. The most common use case involving such files is when merger trees are extracted from an N-body simulation. Galacticus can read several formats of merger tree data file, but its preferred format is that described `here <https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format>`_. A tool "``astrosylva``" (`GitHub <https://github.com/galacticusorg/astrosylva>`_; `PyPI <https://pypi.org/project/astrosylva/>`_; `documentation <https://astrosylva.readthedocs.io/en/latest/>`_) is available to convert merger trees from a variety of formats, including the `Rockstar <https://bitbucket.org/gfcstanford/rockstar>`_/`ConsistentTrees <https://bitbucket.org/pbehroozi/consistent-trees>`_ tree builders, into Galacticus's preferred format.

Broadband Filters
-----------------

To compute luminosities through a given filter, Galacticus requires the response function, :math:`R(\lambda)`, of that filter to be defined. Galacticus follows the convention of :cite:t:`hogg_k_2002` in defining the filter response to be the fraction of incident photons received by the detector at a given wavelength, multiplied by the relative photon response (which will be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a bolometer/calorimeter type detector. Filter response files are stored in the ``static/filters/`` subdirectory of the Galacticus datasets. A large number of filters are provided in that location, but it is straightforward to add new filters. The structure of a filter definition file is shown below, with the ``SDSS_g.xml`` filter response file used as an example:

.. code-block:: none

    <filter>
     <description>SDSS g vacuum (filter+CCD +0 air mass)</description>
     <name>SDSS g</name>
     <origin>Michael Blanton</origin>
     <response>
       <datum>   3630.000      0.0000000E+00</datum>
       <datum>   3680.000      2.2690000E-03</datum>
       <datum>   3730.000      5.4120002E-03</datum>
       <datum>   3780.000      9.8719997E-03</datum>
       <datum>   3830.000      2.9449999E-02</datum>
       .
       .
       .
     </response>
     <effectiveWavelength>4727.02994472695</effectiveWavelength>
     <vegaOffset>0.107430167298754</vegaOffset>
   </filter>

The ``description`` element should provide a description of the filter, while the ``name`` element provides a shorter name. The ``origin`` element should describe from where/whom this filter originated. The ``response`` element contains a list of ``datum`` elements each giving a wavelength (in Angstroms) and response pair. The normalization of the response is arbitrary. The ``effectiveWavelength`` element gives the mean, response-weighted wavelength of the filter and is used, for example, in dust attenuation calculations. The ``vegaOffset`` element gives the value (in magnitudes) which must be added to an AB-system magnitude in this system to place it into the Vega system.
