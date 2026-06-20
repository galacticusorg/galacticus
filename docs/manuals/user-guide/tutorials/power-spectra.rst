Power Spectra
=============

In this tutorial we'll use Galacticus to compute and output the cosmological linear theory power spectrum.

If you haven't already installed Galacticus, you can find the installation instructions `here <https://github.com/galacticusorg/galacticus/wiki#how-do-i-install-and-use-galacticus>`_.

Running the calculation
-----------------------

Galacticus works by reading a parameter file that describes what you want it to do, and gives the numerical values for all relevant parameters (e.g. cosmological density parameters). For this tutorial we'll use the ``parameters/tutorials/powerSpectrum.xml`` parameter file. To run Galacticus on this parameter file issue the following command from the source directory where you installed Galacticus (note that the initial ``$`` represents your command line prompt - don't type it in!):

.. code-block:: bash

   $ ./Galacticus.exe parameters/tutorials/powerSpectrum.xml

If everything is working you should see output which looks something like:

.. code-block:: text

                 ##
      ####        #                  #
     #   #        #             #
    #       ###   #  ###   ### ###  ##   ### ## ##   ##
    #       #  #  #  #  # #  #  #    #  #  #  #  #  #
    #   ###  ###  #   ### #     #    #  #     #  #   #
     #   #  #  #  #  #  # #     #    #  #     #  #    #
      ####  #### ### ####  ###   ## ###  ###   #### ##

    © 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
      2017, 2018, 2019, 2020
      - Andrew Benson

   M: Memory:         code +      nodes +       misc =      total
   M:            22.246Mib +   1.000  b +   9.000  b =  22.246Mib
   M: -> Begin task: power spectrum
   M: <- Done task: power spectrum

and a file ``powerSpectrum.hdf5`` will have been created.

Understanding the input
-----------------------

You can look at the entire parameter file for this tutorial `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/powerSpectrum.xml>`_. Below we'll explore each section of the file to understand what it does. Galacticus parameter files use `XML <https://www.w3schools.com/xml/xml_whatis.asp>`_.

.. code-block:: xml

   <!-- Specify tasks to perform -->
   <task value="powerSpectra"/>

This block tells Galacticus what task we want it to perform - here we're just asking it to compute the power spectrum.

.. code-block:: xml

   <!-- Cosmological parameters -->
   <cosmologyFunctions  value="matterLambda"/>
   <cosmologyParameters value="simple">
     <HubbleConstant  value="70.20000"/>
     <OmegaMatter     value=" 0.27250"/>
     <OmegaDarkEnergy value=" 0.72750"/>
     <OmegaBaryon     value=" 0.04550"/>
     <temperatureCMB  value=" 2.72548"/>
   </cosmologyParameters>

In this block we're setting up a cosmological model. The first option specifies that cosmological functions (e.g. the expansion factor as a function of time) are computed assuming a universe containing matter and a cosmological constant ("lambda"). The next option specifies that we want to use the "simple" cosmological parameters type which just allows us to directly specify the values of all cosmological parameters. Within this option we specify the actual values for the Hubble constant, density parameters, and CMB temperature.

.. code-block:: xml

   <!-- Power spectrum options -->
   <transferFunction value="eisensteinHu1999">
     <neutrinoNumberEffective value="3.046"/>
     <neutrinoMassSummed      value="0.000"/>
   </transferFunction>
   <powerSpectrumPrimordial value="powerLaw">
     <index               value="0.961"/>
     <wavenumberReference value="1.000"/>
     <running             value="0.000"/>
   </powerSpectrumPrimordial>
   <powerSpectrumPrimordialTransferred value="simple"/>
   <cosmologicalMassVariance value="filteredPower">
     <sigma_8   value="0.807" />
     <tolerance value="1.0e-3"/>
   </cosmologicalMassVariance>

Next we set up the linear theory power spectrum. We choose to use the `Eisenstein & Hu (1999) <http://adsabs.harvard.edu/abs/1999ApJ...511....5E>`_ CDM transfer function, a power-law primordial power spectrum (with no running of the index). The "``powerSpectrumPrimordialTransferred``" options specifies how to combine these into the linear theory power spectrum. The "``simple``" option just multiplies the primordial power spectrum by the transfer function squared. Finally, we specify a normalization for the power spectrum by choosing an option for computing the variance of the mass density field ("``filteredPower``" just means that we integrate the power spectrum under some window function - a top-hat function by default - to compute the variance), and then giving a value for the sigma_8 parameter. We also specify a tolerance for the integration of the power spectrum here - this determines how accurately the variance of the mass density field will be computed.

.. code-block:: xml

   <!-- Structure formation options -->
   <linearGrowth value="collisionlessMatter"/>

This section deals with structure growth. We choose to model linear growth assuming collisionless matter.

.. code-block:: xml

   <!-- Output options -->
   <galacticusOutputFileName value="haloMassFunction.hdf5"/>
   <outputTimes value="list">
     <redshifts value="0.0 1.0"/>
   </outputTimes>

Finally we specify a file to output the results to, and a set of redshifts at which we want the power spectrum to be computed and output.

Understanding the output
------------------------

The output file ``powerSpectrum.hdf5`` is an `HDF5 <https://docs.hdfgroup.org/documentation/hdf5/latest/_intro_h_d_f5.html>`_ file - it contains datasets organized into a hierarchical tree of "groups" (like files in directories in a filesystem). Most programming languages have tools to interact with HDF5 files. Here we'll explore the output file using the command-line tools `h5ls <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Ls>`_ and `h5dump <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Dump>`_.

To see the content of the file use:

.. code-block:: console

   $ h5ls powerSpectrum.hdf5
   Build                    Group
   Outputs                  Group
   Parameters               Group
   Version                  Group

The ``h5ls`` command shows the content of the top-level group of the output file. We'll ignore most of the groups in this tutorial, and focus just on ``Outputs`` which contains the main output quantities. We can look inside this group using:

.. code-block:: console

   $ h5ls powerSpectrum.hdf5/Outputs
   Output1                  Group
   Output2                  Group

There are two output groups - corresponding to the two output redshifts specified in our parameter file. Output groups are numbered from the earliest to latest times, so in this case ``Output1`` is the earliest time (*z*\ =1) and ``Output2`` is the latest time (*z*\ =0).

Let's look at the output for *z*\ =0:

.. code-block:: console

   $ h5ls powerSpectrum.hdf5/Outputs/Output2
   alpha                    Dataset {61}
   growthFactor             Dataset {61}
   growthFactorLogDerivative Dataset {61}
   mass                     Dataset {61}
   powerSpectrum            Dataset {61}
   sigma                    Dataset {61}
   wavenumber               Dataset {61}

Finally we find the power spectrum information! We'll ignore most of these datasets in this short tutorial and just focus on ``wavenumber`` and ``powerSpectrum``. You can see the content of these datasets using:

.. code-block:: console

   $ h5ls -d powerSpectrum.hdf5/Outputs/Output2/wavenumber
   wavenumber               Dataset {61}
       Data:
           (0) 0.001, 0.00125892541179417, 0.00158489319246111, 0.00199526231496888, 0.00251188643150958, 0.00316227766016838, 0.00398107170553497, 0.00501187233627273, 0.00630957344480193, 0.00794328234724282, 0.01,
           (11) 0.0125892541179417, 0.0158489319246111, 0.0199526231496888, 0.0251188643150958, 0.0316227766016838, 0.0398107170553497, 0.0501187233627272, 0.0630957344480193, 0.0794328234724281, 0.1, 0.125892541179417,
           (22) 0.158489319246111, 0.199526231496888, 0.251188643150958, 0.316227766016838, 0.398107170553497, 0.501187233627272, 0.630957344480193, 0.794328234724282, 1, 1.25892541179417, 1.58489319246111,
           (33) 1.99526231496888, 2.51188643150958, 3.16227766016838, 3.98107170553497, 5.01187233627272, 6.30957344480193, 7.94328234724282, 9.99999999999999, 12.5892541179417, 15.8489319246111, 19.9526231496888,
           (44) 25.1188643150958, 31.6227766016837, 39.8107170553498, 50.1187233627273, 63.0957344480194, 79.4328234724281, 100, 125.892541179417, 158.489319246111, 199.526231496888, 251.188643150958, 316.227766016838,
           (56) 398.107170553497, 501.187233627272, 630.957344480193, 794.328234724282, 1000

and

.. code-block:: console

   $ h5ls -d powerSpectrum.hdf5/Outputs/Output2/powerSpectrum
   powerSpectrum            Dataset {61}
       Data:
           (0) 19752.765289969, 24275.8564961514, 29641.8593656998, 35882.4927362934, 42954.3023662321, 50704.5590320977, 58839.9340448319, 66901.937276083, 74241.9978355227, 79967.821801863, 82844.7851128752,
           (11) 81397.937733453, 75005.5647253165, 65337.3202631851, 54893.6236162408, 44836.2992006317, 35472.7444185093, 27064.7575914213, 19889.6111811954, 14100.9335074837, 9675.51818378156, 6450.72067033403,
           (22) 4195.23476538419, 2670.86480989484, 1669.46464710895, 1026.96411022781, 622.851533931182, 372.984109335854, 220.787916727825, 129.318408458248, 75.008972039714, 43.1187493270437, 24.5824755625454,
           (33) 13.9083861955195, 7.81421081814523, 4.36211699619959, 2.42070001882387, 1.33605767868285, 0.733738111969742, 0.401109051059581, 0.21834693384002, 0.11839620191705, 0.0639681164280232, 0.034446182833649,
           (44) 0.0184916545920113, 0.00989835104433317, 0.00528431342684245, 0.00281404516782447, 0.00149506907802831, 0.000792580963481224, 0.00041931165533301, 0.000221408739155846, 0.000116698731224898,
           (53) 6.14039039448803e-05, 3.22571173377652e-05, 1.69197052567665e-05, 8.8620191155726e-06, 4.63529660381791e-06, 2.42134808751638e-06, 1.26327961359142e-06, 6.58309379694076e-07

The "``-d``" option to ``h5ls`` tells it to actually output the content of the dataset.

As you'd expect, ``wavenumber`` contains a list of wavenumbers (in units of Mpc\ :sup:`-1`) at which the power spectrum was tabulated, while ``powerSpectrum`` contains the power spectrum itself (in units of Mpc\ :sup:`3`).

You can also explore the content of the HDF5 file using the ``h5dump`` tool, which shows some extra, useful information (I'll cut out much of the very long output for brevity):

.. code-block:: console

   $ h5dump -A -g Outputs/Output2 powerSpectrum.hdf5
   HDF5 "haloMassFunction.hdf5" {
      ATTRIBUTE "outputRedshift" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "outputTime" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 13.7915
         }
      }
      .
      .
      .
      DATASET "powerSpectrum" {
      COMMENT "The power spectrum."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( 61 ) / ( 61 ) }
         ATTRIBUTE "unitsInSI" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 2.938e+67
            }
         }
      }
      .
      .
      .
   }
   }

For example, you can see that the output group contains some "attributes" such as the redshift and time of the output, and that each dataset has a description and a "``unitsInSI``" attribute which is the value you would multiply this dataset by to convert to SI units (and from there you can convert it to whatever units you want).

The full list of output datasets is:

* ``alpha``: The logarithmic slope of :math:`\sigma(M)`: :math:`\alpha = \mathrm{d} \ln \sigma / \mathrm{d} \ln M`;
* ``mass``: The mass scale, :math:`M`, corresponding to the given wavenumber, :math:`k`, defined such that :math:`M = 4 \pi \Omega_\mathrm{M} \rho_\mathrm{crit} / 3 k^3` (in units of :math:`M_\odot`);
* ``powerSpectrum``: The linear theory power spectrum at :math:`z=0`: :math:`P(k)` in units of Mpc\ :sup:`3`;
* ``sigma``: The dimensionless linear theory mass fluctuation at :math:`z=0`: :math:`\sigma(M)`;
* ``wavenumber``: The wavenumber in units of Mpc\ :sup:`-1`.
