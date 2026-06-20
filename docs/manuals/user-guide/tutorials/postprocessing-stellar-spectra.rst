Postprocessing of Stellar Spectra
=================================

Stellar luminosities are computed by convolving a library of simple
stellar populations with the star formation history of each galaxy.
Galacticus allows the spectra of those simple stellar populations to be
postprocessed (after being read from file or internally generated for
example) before they are utilized in the convolution integral. This
postprocessing can modify the spectra in arbitrary ways that depend on
wavelength, redshift, and age of stellar population. Furthermore,
Galacticus allows you to chain together stellar spectra postprocessors into a set
to allow multiple postprocessings to be applied. Furthermore again, you
can define an arbitrary number of sets and apply different sets to
different luminosities.

Typical uses of stellar spectra postprocessors include accounting for
absorption of galaxy light by the intervening , or capturing only the
light from recent star formation [1]_. A full list of the available
postprocessors can be found `here <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_.

If you don't specify a postprocessing set, the "default" set (consisting
of the `inoue2014 <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ postprocessor is applied
to each luminosity calculation. To specify other postprocessing sets add
the following to your parameter file:

.. code-block:: xml

   <luminosityPostprocessSet value="default recent unabsorbed recentUnabsorbed"/>

where one set must be specified for each luminosity specified in the
``luminosityFilter`` parameter. Note that set names can be
reused in order to apply the same postprocessor set to multiple
luminosities.

The chain of postprocessors to apply for each set is then specified as
follows:

.. code-block:: xml

   <stellarPopulationSpectraPostprocessorBuilder value="lookup">
     <names value="default recent unabsorbed recentUnabsorbed"/>
     <stellarPopulationSpectraPostprocessor      value="inoue2014"/>
     <stellarPopulationSpectraPostprocessor      value="sequence"  >
       <stellarPopulationSpectraPostprocessor    value="inoue2014"/>
       <stellarPopulationSpectraPostprocessor    value="recent"    >
   	   <timeLimit value="1.0e-2"/>
       </stellarPopulationSpectraPostprocessor>
     </stellarPopulationSpectraPostprocessor>
     <stellarPopulationSpectraPostprocessor      value="identity" />
     <stellarPopulationSpectraPostprocessor      value="recent"    >
       <timeLimit value="1.0e-2"/>
     </stellarPopulationSpectraPostprocessor>
   </stellarPopulationSpectraPostprocessorBuilder>

In this case we've constructed four postprocessor sets using the ``lookup`` ``stellarPopulationSpectraPostprocessorBuilder``. This builder is responsible for constructing a suitable postprocessor from the names in the ``luminosityPostprocessSet`` parameter. The ``lookup`` implementation simply takes a list of postprocessor set names, and a corresponding list of ``stellarPopulationSpectraPostprocessor``\ s and select the relevant one based on the name.

In the above the ``default`` set applies the ``inoue2014`` absorption postprocessor, while the ``recent`` set applies both the ``inoue2014`` absorption postprocessor, followed by the ``recent`` postprocessor to retain only recently emitted light. The ``unabsorbed`` set ignores absorption entirely—it does this by using the ``identity`` postprocessor which leaves the spectrum unaffected. Finally, the ``recentUnabsorbed`` set applies only the ``recent`` filter while ignoring absorption.

In this way it is relatively easy to extract multiple different measures of luminosity from a Galacticus model. For example, you could construct four postprocessor sets, each corresponding to one of the four different
absorption models (``lycSuppress``, ``madau1995``, ``meiksin2006``, and ``inoue2014``) and apply these to the same luminosity filter to assess how luminosity depends on the model used.

.. [1] Perhaps so that additional dust extinction can be applied to the light of recently formed stars.
