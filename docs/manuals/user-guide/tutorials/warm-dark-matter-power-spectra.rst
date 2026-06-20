Warm Dark Matter Power Spectra
==============================

In this tutorial we'll use Galacticus to compute and output the power spectrum for a warm dark matter model. You should first follow the `CDM power spectrum tutorial <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/power-spectra.html>`_ which explains how to run a power spectrum calculation, and describes the outputs.

Running the calculation
-----------------------

To run an example of a warm dark matter power spectrum calculation, do the following:

.. code-block:: console

   $ ./Galacticus.exe parameters/tutorials/powerSpectrumWarmDarkMatter.xml

If everything is working a file ``powerSpectrumnWarmDarkMatter.hdf5`` will have been created.

Understanding the input
-----------------------

You can look at the entire parameter file for this tutorial `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/powerSpectrumWarmDarkMatter.xml>`_. Below we'll explore just those sections of the parameter file specific to the warm dark matter calculation - everything else is just as it would be for the equivalent CDM calculation.

.. code-block:: xml

   <!-- Use a thermal WDM particle - mass is in keV -->
   <darkMatterParticle value="WDMThermal">
     <degreesOfFreedomEffective value="1.5" />
     <mass value="3.0" />
   </darkMatterParticle>

In the above we explicitly specify a thermal warm dark matter particle, with effective degrees of freedom :math:`g_\mathrm{X}=1.5` and a mass of :math:`m_\mathrm{X}=3.0` keV. In the CDM case we didn't specify the type of dark matter particle, since CDM is the default.

.. code-block:: xml

   <!-- Use the Bode et al. (2001) transfer function for thermal WDM -->
   <transferFunction value="bode2001">
     <epsilon value="0.359" />
     <eta value="3.810" />
     <nu value="1.100" />
     <!-- Bode2001 transfer function works by modifying a CDM transfer function - so feed it a CDM transfer function here -->
     <transferFunction value="eisensteinHu1999">
       <!-- Feed this transfer function a CDM particle - otherwise it will see the WDM particle defined above and complain that it
            can not compute WDM transfer functions -->
       <darkMatterParticle value="CDM" />
       <neutrinoNumberEffective value="3.046"/>
       <neutrinoMassSummed value="0.000"/>
     </transferFunction>
   </transferFunction>
   <!-- When computing sigma(M) for power spectra with a cut off it's better to use a filter that is sharp in k-space, instead of
        the usual real-space top-hat (which introduces artificial halos below the cut-off scale -->
   <cosmologicalMassVariance value="filteredPower">
     <monotonicInterpolation value="true" />
     <nonMonotonicIsFatal value="false" />
     <powerSpectrumWindowFunction value="sharpKSpace">
       <normalization value="2.5" />
     </powerSpectrumWindowFunction>
     <sigma_8 value="0.8111" />
     <tolerance value="3.0e-4" />
     <toleranceTopHat value="3.0e-4" />
   </cosmologicalMassVariance>

Here we use the `Bode et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...556...93B>`_ fitting function for the warm dark matter transfer function. Note that this works by modifying the CDM transfer function (supplied here using the ``eisensteinHu1999`` method), and that we explicitly pass a CDM dark matter particle to this class (otherwise it would find the WDM particle class and refuse to work, since it doesn't know how to compute the transfer function for WDM). For warm dark matter, we switch to using a sharp-in-:math:`k`-space filter to compute :math:`\sigma(M)`, as it avoids spurious halos below the cut-off scale which result for warm dark matter power spectra using the usual top-hat-in-real-space filter (see `Benson et al. (2013) <http://adsabs.harvard.edu/abs/2013MNRAS.428.1774B>`_).

Understanding the output
------------------------

The output file ``powerSpectrumWarmDarkMatter.hdf5`` has exactly the same format as for the CDM case, described in the power spectrum tutorial.
