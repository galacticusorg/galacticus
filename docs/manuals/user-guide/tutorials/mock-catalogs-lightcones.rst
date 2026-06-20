Generating Mock Catalogs with Lightcones
========================================

Suppose that you want to create a catalog of galaxies as would be found
in a survey of an area of the sky out to some redshift. Such a "mock
catalog" can be built by populating with galaxies all of the dark matter
halos which happen to lie within the cone which that area makes as it is
projected from the observer through the Universe.

Generating such a mock catalog using Galacticus involves first extracting the
halos (and their merger trees) within this "lightcone" from a suitable N-body
simulation, and then processing them through. In this tutorial, we will assume
that you have merger trees from a cosmological simulation available in
Galacticus' merger tree file format.

Each such merger tree file can then be run through Galacticus in the
usual way (see the tutorial on "`Using N-body Merger Trees <https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Using-N-body-Merger-Trees>`_"). Outputs should be requested at every snapshot (up to the largest redshift to be considered), and
the ``lightcone`` filter should be used to cause only those
galaxies which intersect the lightcone to be output—for example:

.. code-block:: xml

   <!-- Set output redshifts to the snapshots in the Millennium Simulation. -->
   <outputRedshifts value=
      "0.0000 0.0199 0.0414 0.0645 0.0893 0.1159 0.1444 0.1749 0.2075 0.2425
       0.2798 0.3197 0.3623 0.4079 0.4566 0.5086 0.5642 0.6236 0.6871 0.7550
       0.8277 0.9055 0.9887 1.0779 1.1734 1.2758 1.3857 1.5036 1.6303 1.7663
       1.9126 2.0700 2.2395 2.4220 2.6189 2.8312 3.0604 3.3081 3.5759 3.8657
       4.1795 4.5196 4.8884 5.2888 5.7239 6.1968"
   />

   <!-- Add a lightcone filter with the required geometry -->
   <mergerTreeOutput>
     <galacticFilter value="lightcone"/>
   </mergerTreeOutput>

   <!-- Switch on output of lightcone data -->
   <outputLightconeData value="true"/>

   <!-- Prune away trees not appearing in the lightcone -->
   <mergerTreeOperator value="pruneLightcone"/>

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

In the above causes lightcone coordinate information (i.e. the position
and velocity of each galaxy in a coordinate system with axes aligned
along the line of sight of the lightcone and parallel to the two edges
of the square field of view, along with the redshift) to be output (see the documentation on the ``lightcone`` `nodePropertyExtractor <https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html>`_), and is set to ``pruneLightcone`` to
cause any merger trees which have no nodes within the lightcone volume
to be pruned away (as there is no need to process them). Finally, the
``geometryLightcone`` parameter describes the geometry of
the lightcone to be used—see `the documentation <https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html>`_ for
details.
