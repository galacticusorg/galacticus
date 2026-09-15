Generating Mock Catalogs with Lightcones
========================================

Suppose that you want to create a catalog of galaxies as would be found
in a survey of an area of the sky out to some redshift. Such a "mock
catalog" can be built by populating with galaxies all of the dark matter
halos which happen to lie within the cone which that area makes as it is
projected from the observer through the Universe.

Generating such a mock catalog using Galacticus involves first extracting the
halos (and their merger trees) within this "lightcone" from a suitable N-body
simulation, and then processing them through Galacticus. In this tutorial, we will assume
that you have merger trees from a cosmological simulation available in
Galacticus' merger tree file format.

Each such merger tree file can then be run through Galacticus in the
usual way (see the tutorial on :doc:`Using N-body Merger Trees <nbody-merger-trees>`). For example:

.. code-block:: xml

   <!-- Set output redshifts to the minimum and maximum for which the lightcone is to be constructed. -->
   <outputTimes value="list">
      <redshifts value="0.0000 6.1968"/>
   </outputTimes>

   <!-- Prune away trees not appearing in the lightcone -->
   <mergerTreeOperator value="pruneLightcone">
     <splitTrees value="true"/>
   </mergerTreeOperator>

   <!-- Set up the new lightcone output using an outputter passed to the lightcone-crossing merger tree evolution timestepper -->
   <mergerTreeEvolveTimestep value="lightconeCrossing">
     <mergerTreeOutputter value="standard">
      <outputsGroupName value="Lightcone"/>
     </mergerTreeOutputter>
   </mergerTreeEvolveTimestep>
   <mergerTreeOutputter value="null"/>

   <!-- Specify lightcone geometry -->
   <geometryLightcone value="square">
     <origin value="0 0 0"/>
     <unitVector1 value=" 1 1  1"/>
     <unitVector2 value=" 0 1 -1"/>
     <unitVector3 value="-2 1  1"/>
     <lengthReplication value="500"/>
     <lengthHubbleExponent value="-1"/>
     <lengthUnitsInSI value="3.08567758135e22"/>
     <angularSize value="0.5"/>
     <timeEvolvesAlongLightcone value="true"/>
     <redshift value=
      "0.0000 0.0199 0.0414 0.0645 0.0893 0.1159 0.1444 0.1749 0.2075 0.2425
       0.2798 0.3197 0.3623 0.4079 0.4566 0.5086 0.5642 0.6236 0.6871 0.7550
       0.8277 0.9055 0.9887 1.0779 1.1734 1.2758 1.3857 1.5036 1.6303 1.7663
       1.9126 2.0700 2.2395 2.4220 2.6189 2.8312 3.0604 3.3081 3.5759 3.8657
       4.1795 4.5196 4.8884 5.2888 5.7239 6.1968"
     />
   </geometryLightcone>

   <!-- Position interpolation -->
   <nodeOperator value="positionInterpolated">
     <wrapPeriodic value="false"/> <!-- Do not wrap interpolated positions back into the box - this is not needed as we replicate the box. -->
     <lengthBox    value="=[geometryLightcone/lengthReplication]*(([cosmologyParameters/HubbleConstant]/100.0)^[geometryLightcone/lengthHubbleExponent])*([geometryLightcone/lengthUnitsInSI]/3.08567758e+22)"/>
   </nodeOperator>

In the above,
:galacticus-class:`mergerTreeEvolveTimestepLightconeCrossing` causes
galaxies to be output at the time at which they cross the
lightcone - this can be combined with any other timesteppers to control
evolution. Note that for this to work, fully time-dependent positions
must be available for galaxies. Typically this is achieved by using the
:galacticus-class:`nodeOperatorPositionInterpolated` operator as shown
above (which should be incorporated into the list of all
``nodeOperator``\ s used in the model).

We use :galacticus-class:`mergerTreeOperatorPruneLightcone` to
cause any merger trees which have no nodes within the lightcone volume
to be pruned away (as there is no need to process them). Finally, the
:galacticus-class:`geometryLightcone` parameter describes the geometry of
the lightcone to be used.
