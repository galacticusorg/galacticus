Subsampling of Merger Tree Branches
===================================

In this tutorial we will build on the previous tutorial on `generating dark matter only merger trees <https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Dark-matter-only-merger-trees>`_ to explore Galacticus' functionality to modify how trees are built. Specifically, we'll subsample merger tree branches - keeping only a fraction of the lower mass branches. This can be useful when exploring extremely high-resolution trees to keep both memory and CPU requirements manageable.

How subsampling works
---------------------

A merger tree is constructed by starting from a single halo (often, but not necessarily, at :math:`z=0`) and following it backward in time, as it splits apart into smaller, progenitor halos. (This splitting as we follow the tree backward in time is equivalent to merging if we follow the tree forward in time.) Each splitting results in the formation of a new branch of the tree. Normally, each new branch is followed back in time, undergoing further splitting until some mass resolution limit is reached and the branch is no longer followed.

This results in large numbers of low mass branches. We can choose to subsample these - that is, for each new branch we randomly decide whether or not to follow that branch, or whether to instead ignore it (removing it completely from the tree). The probability to subsample a given branch can be a function of the branch mass - allowing us to keep a smaller fraction of lower mass branches. The probability used by Galacticus is

.. math::

   p(M) = \left\{ \begin{array}{ll} p_0 (m/m_0)^\alpha & \hbox{ if } m < m_0, \\ 1 & \hbox{ otherwise.} \end{array} \right.

Running the calculation
-----------------------

To run the same dark matter only merger tree model as in the previous `tutorial <https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Dark-matter-only-merger-trees>`_ do the following:

.. code-block:: console

   $ ./Galacticus.exe parameters/tutorials/darkMatterOnlyMergerTreesSubsampled.xml

This `parameter file <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/darkMatterOnlyMergerTreesSubsampled.xml>`_ has two additions relative to the `original <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/darkMatterOnlyMergerTrees.xml>`_:

.. code-block:: xml

   <mergerTreeMassResolution  value="fixed"    >
     <massResolution             value="1.0e6"/>
   </mergerTreeMassResolution>
   <mergerTreeBuildController value="subsample">
     <massThreshold              value="5.0e9"/>
     <subsamplingRateAtThreshold value="1.0"  />
     <exponent                   value="1.0"  />
   </mergerTreeBuildController>

The first element sets a much higher mass resolution of :math:`10^6\mathrm{M}_\odot`. The second adds a ``mergerTreeBuildController`` which performs the subsampling. In this case the parameters of the subsampling probability model are set to be:

* :math:`p_0` = ``subsamplingRateAtThreshold`` = 1.0
* :math:`m_0` = ``massThreshold`` = :math:`5\times10^9\mathrm{M}_\odot`
* :math:`\alpha` = ``exponent`` = 1.0

So, for example, branches of mass :math:`10^7\mathrm{M}_\odot` are kept with probability 0.002.

Understanding the output
------------------------

The resulting output file has the same structure as in the previous `tutorial <https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Dark-matter-only-merger-trees>`_. However, if we examine the masses of the halos produced we will now find some much lower mass halos:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/nodeData/basicMass
   basicMass                Dataset {30/Inf}
       Data:
           (0) 2759486893.83384, 2237815.30539422, 5983576.01413539, 18132446785.7874, 1550665478.91184, 1195595360.20214, 1716353527.78637, 914573516.240639, 23964746.1711591,
           (9) 46214619814.4762, 174578529.205376, 26294735509.1231, 11008079013.6822, 29413575.8834082, 13737902521.1931, 582630977.780802, 6915643378.85546, 6868563.99249936,
           (18) 2061933039.46783, 2499819982.83892, 4986276.78940057, 322312344.18355, 38948189.6417105, 46609707.019351, 853076206.708892, 20170000.2038353, 171002134.762819,
           (27) 103151283.610772, 5767875.77411493, 154165799703.412

Note that, because of the subsampling, the number of subhalos of each mass is no longer representative of a real halo. To account for this in any analysis, e.g. when computing a subhalo mass function, it is necessary to re-weight the contribution of each subhalo in the output. Galacticus provides a dataset ``nodeSubsamplingWeight`` which gives the weight that should be applied to each subhalo to correct for the effects of subsampling:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/nodeData/nodeSubsamplingWeight
   nodeSubsamplingWeight    Dataset {30/Inf}
       Data:
           (0) 1.81193105543377, 4048.43744491812, 835.620703771152, 1, 3.22442207426239, 4.1820168983883, 2.91315275032449, 5.46702907006594, 208.639806334246, 1, 28.6404062559031, 1, 1,
           (13) 169.989532038518, 1, 8.58176133896042, 1, 2065.50458763648, 2.42490900737031, 2.00014402409959, 13894.4136303381, 15.5129025934936, 128.375671526601, 107.273791657265,
           (24) 5.86114108057198, 247.892907757595, 29.2394010573905, 48.4724942334878, 866.870264862326, 1

Warning!
--------

.. warning::

   Subsampling of merger tree branches is useful to provide rapid evaluation of high-resolution merger trees. However, it is only guaranteed to give physically correct results if the missing branches would not have affected the evolution and properties of any of the remaining branches.

   In the case of dark matter only merger trees this is often true (the mass of the missing branches will simply be treated as smooth accretion) - but not always. Models which compute, for example, halo concentrations based on the halo merger history will likely give incorrect results if subsampling is applied.

   Similarly, models which include baryonic (galaxy) physics will not necessarily give correct results if subsampling is applied - the galaxies which would have formed in the missing branches may have influenced (via their outflows, or by merging) galaxies forming in the remaining branches of the merger tree.
