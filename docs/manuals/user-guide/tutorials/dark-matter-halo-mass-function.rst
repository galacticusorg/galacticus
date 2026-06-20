Dark Matter Halo Mass Function
==============================

In this tutorial we'll use Galacticus to compute and output the mass function of dark matter halos.

If you haven't already installed Galacticus, you can find the installation instructions `here <https://github.com/galacticusorg/galacticus/wiki#how-do-i-install-and-use-galacticus>`_.

Running the calculation
-----------------------

Galacticus works by reading a parameter file that describes what you want it to do, and gives the numerical values for all relevant parameters (e.g. cosmological density parameters). For this tutorial we'll use the ``parameters/tutorials/haloMassFunction.xml`` parameter file. To run Galacticus on this parameter file issue the following command from the source directory where you installed Galacticus (note that the initial ``$`` represents your command line prompt - don't type it in!):

.. code-block:: console

   $ ./Galacticus.exe parameters/tutorials/haloMassFunction.xml

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
   M: -> Begin task: halo mass function
   M: <- Done task: halo mass function

and a file ``haloMassFunction.hdf5`` will have been created.

Understanding the input
-----------------------

You can look at the entire parameter file for this tutorial `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/haloMassFunction.xml>`_. Below we'll explore each section of the file to understand what it does. Galacticus parameter files use `XML <https://www.w3schools.com/xml/xml_whatis.asp>`_.

.. code-block:: xml

   <!-- Specify tasks to perform -->
   <task value="haloMassFunction"/>

This block tells Galacticus what task we want it to perform - here we're just asking it to compute a halo mass function.

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
     <sigma_8 value="0.807"/>
   </cosmologicalMassVariance>

Next we set up the linear theory power spectrum. We choose to use the `Eisenstein & Hu (1999) <http://adsabs.harvard.edu/abs/1999ApJ...511....5E>`_ CDM transfer function, a power-law primordial power spectrum (with no running of the index). The "``powerSpectrumPrimordialTransferred``" options specifies how to combine these into the linear theory power spectrum. The "``simple``" option just multiplies the primordial power spectrum by the transfer function squared. Finally, we specify a normalization for the power spectrum by choosing an option for computing the variance of the mas density field ("``filteredPower``" just means that we integrate the power spectrum under some window function - a top-hat function by default - to compute the variance), and then giving a value for the sigma_8 parameter.

.. code-block:: xml

   <!-- Structure formation options -->
   <linearGrowth          value="collisionlessMatter"                       />
   <haloMassFunction      value="tinker2008"                                />
   <criticalOverdensity   value="sphericalCollapseCllsnlssMttrCsmlgclCnstnt"/>
   <virialDensityContrast value="sphericalCollapseCllsnlssMttrCsmlgclCnstnt"/>

This section deals with structure growth. We choose to model linear growth assuming collisionless matter, specify that critical overdensities for collapse of halos (as used in `Press-Schechter <http://adsabs.harvard.edu/abs/1974ApJ...187..425P>`_-type models) and virial density contrasts of halos are computed using the spherical collapse model, and that the halo mass function should be computed using the fitting function of `Tinker et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...688..709T>`_.

.. code-block:: xml

   <!-- Dark matter halo structure options -->
   <darkMatterProfileDMO           value="NFW"    />
   <darkMatterProfileConcentration value="gao2008"/>

We next specify how dark matter halo density profiles should be computed. We choose to use `NFW <http://adsabs.harvard.edu/abs/1997ApJ...490..493N>`_ halos, and the `Gao et al. (2008) <http://adsabs.harvard.edu/abs/2008MNRAS.387..536G>`_ fitting function to assign them concentrations.

.. code-block:: xml

   <!-- Output options -->
   <galacticusOutputFileName value="haloMassFunction.hdf5"/>
   <outputTimes value="list">
     <redshifts value="0.0 1.0"/>
   </outputTimes>

Finally we specify an file to output the results to, and a set of redshifts at which we want the halo mass function to be computed and output.

Controlling the range of masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The range and resolution of the mass tabulation can be controlled by three optional parameters, passed to the ``haloMassFunction`` task, e.g.

.. code-block:: xml

   <task value="haloMassFunction">
    <haloMassMinimum value="1.0e06"/>
    <haloMassMaximum value="1.0e15"/>
    <pointsPerDecade value="30"    />
   </task>

The parameters have the following meanings:

* ``haloMassMinimum``: The lowest mass halo (in units of :math:`M_\odot`) at which to tabulate;
* ``haloMassMaximum``: The highest mass halo (in units of :math:`M_\odot`) at which to tabulate;
* ``pointsPerDecade``: The number of points per decade of halo mass at which to tabulate.

Understanding the output
------------------------

The output file ``haloMassFunction.hdf5`` is an `HDF5 <https://docs.hdfgroup.org/documentation/hdf5/latest/_intro_h_d_f5.html>`_ file - it contains datasets organized into a hierarchical tree of "groups" (like files in a directories in a filesystem). Most programming languages have tools to interact with HDF5 files. Here we'll explore the output file using the command-line tools `h5ls <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Ls>`_ and `h5dump <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Dump>`_.

To see the content of the file use:

.. code-block:: console

   $ h5ls haloMassFunction.hdf5
   Build                    Group
   Outputs                  Group
   Parameters               Group
   Version                  Group
   cosmology                Group

The ``h5ls`` command shows the content of the top-level group of the output file. We'll ignore most of the groups in this tutorial, and focus just on ``Outputs`` which contains the main output quantities. We can look inside this group using:

.. code-block:: console

   $ h5ls haloMassFunction.hdf5/Outputs
   Output1                  Group
   Output2                  Group

There are two output groups - corresponding to the two output redshifts specified in our parameter file. Output groups are numbered from the earliest to latest times, so in this case ``Output1`` is the earliest time (:math:`z=1`) and ``Output2`` is the latest time (:math:`z=0`).

Let's look at the output for :math:`z=0`:

.. code-block:: console

   $ h5ls haloMassFunction.hdf5/Outputs/Output2
   haloAlpha                Dataset {51}
   haloBias                 Dataset {51}
   haloMass                 Dataset {51}
   haloMassFractionCumulative Dataset {51}
   haloMassFunctionCumulative Dataset {51}
   haloMassFunctionLnM      Dataset {51}
   haloMassFunctionLnMBinAveraged Dataset {51}
   haloMassFunctionM        Dataset {51}
   haloMassFunctionNuFNu    Dataset {51}
   haloPeakHeightNu         Dataset {51}
   haloScaleRadius          Dataset {51}
   haloSigma                Dataset {51}
   haloVelocityMaximum      Dataset {51}
   haloVirialRadius         Dataset {51}
   haloVirialTemperature    Dataset {51}
   haloVirialVelocity       Dataset {51}

Finally we find the halo mass function information! We'll ignore most of these datasets for now (more on them below) and just focus on ``haloMass`` and ``haloMassFunctionLnM``. You can see the content of these datasets using:

.. code-block:: console

   $ h5ls -d haloMassFunction.hdf5/Outputs/Output2/haloMass
   haloMass                 Dataset {51}
       Data:
           (0) 10000000000, 12589254117.9417, 15848931924.6112, 19952623149.6888, 25118864315.0958, 31622776601.6838,
           (6) 39810717055.3498, 50118723362.7273, 63095734448.0193, 79432823472.4281, 100000000000, 125892541179.417,
           (12) 158489319246.111, 199526231496.888, 251188643150.958, 316227766016.838, 398107170553.497, 501187233627.273,
           (18) 630957344480.194, 794328234724.281, 999999999999.999, 1258925411794.16, 1584893192461.11, 1995262314968.88,
           (24) 2511886431509.59, 3162277660168.38, 3981071705534.97, 5011872336272.72, 6309573444801.92, 7943282347242.82,
           (30) 10000000000000, 12589254117941.7, 15848931924611.1, 19952623149688.8, 25118864315095.8, 31622776601683.8,
           (36) 39810717055349.6, 50118723362727.1, 63095734448019.3, 79432823472427.8, 99999999999999.4, 125892541179417,
           (42) 158489319246111, 199526231496888, 251188643150957, 316227766016837, 398107170553498, 501187233627273,
           (48) 630957344480194, 794328234724281, 999999999999999

and

.. code-block:: console

   $ h5ls -d haloMassFunction.hdf5/Outputs/Output2/haloMassFunctionLnM
   haloMassFunctionLnM      Dataset {51}
       Data:
           (0) 0.102089448017817, 0.0828171612887022, 0.0672008803636811, 0.0545437280552853, 0.0442819649125888,
           (5) 0.0359599258177283, 0.0292090625629812, 0.0237311702401966, 0.01928494203462, 0.0156750547925035,
           (10) 0.0127433541026886, 0.0103617588580126, 0.00842650559954858, 0.00685351336731697, 0.00557462703726862,
           (15) 0.00453458700517561, 0.00368857537608279, 0.00300023229695999, 0.00244005052208817, 0.00198407691243919,
           (20) 0.00161286097807425, 0.00131060349717597, 0.00106446708551513, 0.000864017697775054, 0.000700771988079618,
           (25) 0.000567830431014216, 0.000459579860410537, 0.000371452186163999, 0.000299728655123172, 0.000241380985979144,
           (30) 0.000193942397008319, 0.000155402865211002, 0.000124124033905483, 9.87700567388156e-05, 7.82513690942053e-05,
           (35) 6.16789544141317e-05, 4.83271231039452e-05, 3.76032088600067e-05, 2.90228844294305e-05, 2.21900393502974e-05,
           (40) 1.6780375318062e-05, 1.25280265154478e-05, 9.21465310325886e-06, 6.66056119352857e-06, 4.71749332015699e-06,
           (45) 3.26280662152348e-06, 2.19481196482058e-06, 1.42909146602241e-06, 8.95638243190209e-07, 5.3668178845726e-07,
           (50) 3.05048538970939e-07

The "``-d``" option to ``h5ls`` tells it to actually output the content of the dataset.

As you'd expect, ``haloMass`` contains a list of halo masses (in units of Solar masses) at which the halo mass function was tabulated, while ``haloMassFunctionLnM`` contains the halo mass function (specifically the number of halos per unit mass, per natural logarithm of halo mass; in units of :math:`\mathrm{Mpc}^{-3}`).

You can also explore the content of the HDF5 file using the ``h5dump`` tool, which shows some extra, useful information (I'll cut out much of the very long output for brevity):

.. code-block:: console

   $ h5dump -A -g Outputs/Output2 haloMassFunction.hdf5
   HDF5 "haloMassFunction.hdf5" {
   GROUP "Outputs/Output2" {
      COMMENT "Data for output number 2"
      ATTRIBUTE "criticalOverdensity" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 1.67463
         }
      }
      .
      .
      .
      DATASET "haloMass" {
      COMMENT "The mass of the halo."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( 51 ) / ( 51 ) }
         ATTRIBUTE "unitsInSI" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 1.98892e+30
            }
         }
      }
      .
      .
      .
   }
   }

For example, you can see that the output group contains some "attributes" such as the critical linear theory overdensity for collapse at this epoch, and that each dataset has a description and a "``unitsInSI``" attribute which is the value you would multiply this dataset by to convert to SI units (and from there you can convert it to whatever units you want).

More details on the output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Outputs/OutputN`` groups contain atrributes and datasets which give properties at the corresponding output time as follows:

* ``massHaloCharacteristic``: The characteristic mass scale (in units of :math:`M_\odot`), :math:`M_*`, at which :math:`\sigma(M)=\delta_\mathrm{c}(z)`;
* ``criticalOverdensity``: The critical overdensity for collapse of halos, :math:`\delta_\mathrm{c}`;
* ``outputExpansionFactor``: The expansion factor;
* ``growthFactor``: The linear growth factor;
* ``outputRedshift``: The redshift;
* ``outputTime``: The cosmic time (in units of Gyr);
* ``virialDensityContrast``: The virial density contrast of halos.
* ``haloMass``: The mass of the halo, :math:`M_\mathrm{halo}` (in :math:`M_\odot`);
* ``haloSigma``: The root-variance of the mass field smoothed in top-hat spheres, :math:`\sigma(M)`;
* ``haloAlpha``: The loagrithmic gradient of the root-variance of the mass field smoothed in top-hat spheres with mass, :math:`\mathrm{d} \log\sigma(M)/\mathrm{d} \log M_\mathrm{halo}`;
* ``haloPeakHeightNu``: The peak height of the halo, :math:`\nu = \delta_\mathrm{c}/\sigma(M)`;
* ``haloMassFunctionM``: The halo mass function per halo mass, :math:`\mathrm{d}n/\mathrm{d}M_\mathrm{halo}` (in units of :math:`\mathrm{Mpc}^{-3} M_\odot^{-1}`);
* ``haloMassFunctionLnM``: The halo mass function per logarithmic halo mass, :math:`\mathrm{d}n/\mathrm{d}\log M_\mathrm{halo}` (in units of :math:`\mathrm{Mpc}^{-3}`);
* ``haloMassFunctionLnMBinAveraged``: The halo mass function per logarithmic halo mass averaged over the finite width of the bin (in units of :math:`\mathrm{Mpc}^{-3}`);
* ``haloMassFunctionNuFNu``: The halo mass function defined in terms of the peak-height parameter, :math:`\nu F(\nu)`;
* ``haloMassFractionCumulative``: The mass fraction in halos above the current halo mass;
* ``haloMassFunctionCumulative``: The cumulative number of halos per unit volume above the current halo mass (in units of :math:`\mathrm{Mpc}^{-3}`);
* ``subhaloMassFunctionCumulative``: The cumulative number of sub-halos per unit volume above the current halo mass (in units of :math:`\mathrm{Mpc}^{-3}`);
* ``haloBias``: The large scale linear theory bias of the halo;
* ``haloVirialRadius``: The virial radius (in units of Mpc) of the current halo mass;
* ``haloVirialTemperature``: The virial temperature (in units of Kelvin) of the current halo mass;
* ``haloVirialVelocity``: The virial velocity (in units of km/s) of the current halo mass;
* ``haloScaleRadius``: The scale radius (in units of Mpc) of the current halo mass;
* ``haloVelocityMaximum``: The peak of the rotation curve (in units of km/s) of the current halo mass;

Dimensionful datasets have a ``unitsInSI`` attribute that gives their units in the SI system.

Additionally, an attribute giving the critical density of the universe (in units of :math:`M_\odot \mathrm{Mpc}^{-3}`) is written to the ``cosmology`` group and, if such a scale is well-defined, an attribute given the mass corresponding to the scale at which the power spectrum is reduce by half relative to a CDM power spectrum in units of :math:`M_\odot` is written as ``massHalfMode`` to the ``powerSpectrum`` group.
