Instantaneous and Non-Instantaneous Recycling
=============================================

Galacticus tracks metals produced by stellar evolution. How this is done is controlled by the ``stellarPopulationProperties`` parameter which selects between different ways of computing properties of stellar populations.

Choosing ``stellarPopulationProperties=instantaneous`` will cause Galacticus to use the
instantaneous recycling approximation for all calculations of stellar
populations. The recycled fraction and metal yield are determined from the ``stellarPopulation`` provided by a ``stellarPopulationSelector``.

Setting ``stellarPopulationProperties=noninstantaneous`` causes Galacticus to use a fully
non-instantaneous, metal-dependent calculation of recycling, metal
production and rates. These rates are determined from the ``stellarPopulation`` provided by a ``stellarPopulationSelector``. However, it is possible to force this method to
operate in the instantaneous recycling approximation limit (which can be
useful for testing and comparison) by setting:

.. code-block:: xml

   <stellarPopulations value="standard">
     <!-- Force the calculation of recycling, yields etc. to   -->
     <!-- be done assuming instantaneous recycling             -->
     <instantaneousRecyclingApproximation value="true"/>
     <!-- Set the recycled fraction and yield -->
     <recycledFraction value="0.35"/>
     <metalYield       value="0.02"/>
   </stellarPopulationProperties>

where the recycled fraction and metal yield are specified directly, or

.. code-block:: xml

   <stellarPopulations value="standard">
     <!-- Force the calculation of recycling to be done       -->
     <!-- assuming the instantaneous recycling approximation  -->
     <instantaneousRecyclingApproximation value="true"/>
     <!-- Set the mass of stars which should be used as the    -->
     <!-- dividing line between long-lived and instantaneously -->
     <!-- evolving in this approximation.                      -->
     <massLongLived value="1.0"/>
     <!-- Set the effective age of populations to use in this -->
     <!-- approximation when computing SNe numbers.           -->
     <ageEffective value="13.8"/>
   </stellarPopulationProperties>

in which case the recycled fraction and metal yield will be computed assuming that all stars with mass greater than ``massLongLived`` have fully evolved, and energy input (from stellar winds and supernovae) will be computed assuming that stellar populations instantaneously reach an age of ``ageEffective``.

Similar options are available to control whether metal yields and energy input from stellar populations are computed using the fully non-instantaneous or instantaneous approximations, e.g.:

.. code-block:: xml

   <stellarPopulations value="standard">
     <!-- Force the calculation of recycling to be done       -->
     <!-- assuming the instantaneous recycling approximation  -->
     <instantaneousRecyclingApproximation value="true"/>
     <!-- Force the calculation of yields to be done          -->
     <!-- assuming the instantaneous recycling approximation  -->
     <instantaneousYieldApproximation value="true"/>
     <!-- Force the calculation of energy input to be done    -->
     <!-- assuming the instantaneous recycling approximation  -->
     <instantaneousEnergyInputApproximation value="true"/>
     <!-- Set the mass of stars which should be used as the    -->
     <!-- dividing line between long-lived and instantaneously -->
     <!-- evolving in this approximation.                      -->
     <massLongLived value="1.0"/>
     <!-- Set the effective age of populations to use in this -->
     <!-- approximation when computing SNe numbers.           -->
     <ageEffective value="13.8"/>
   </stellarPopulationProperties>
