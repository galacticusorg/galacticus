Dark-Matter-Only Merger Trees
=============================

In this tutorial we'll use Galacticus to simulate a set of dark matter only `merger trees <https://en.wikipedia.org/wiki/Galaxy_merger#Merger_history_trees>`_. This tutorial assumes that you've already installed Galacticus (if you haven't, you can find the installation instructions `here <https://github.com/galacticusorg/galacticus/wiki#how-do-i-install-and-use-galacticus>`_), and have worked through the earlier tutorials so understand the basics of Galacticus parameter and output files.

Running the calculation
-----------------------

Galacticus works by reading a parameter file that describes what you want it to do, and gives the numerical values for all relevant parameters (e.g. cosmological density parameters). For this tutorial we'll use the ``parameters/tutorials/darkMatterOnlyMergerTrees.xml`` parameter file. To run Galacticus on this parameter file issue the following command from the source directory where you installed Galacticus:

.. code-block:: bash

   $ ./Galacticus.exe parameters/tutorials/darkMatterOnlyMergerTrees.xml

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

   MM: Memory:         code +      nodes +       misc =      total
   MM:            22.242Mib +   1.000  b +  13.149kib =  22.255Mib
   MM: -> Begin task: merger tree evolution
    7:     Memory:         code +      nodes +       misc =      total
    7:                22.242Mib +   3.114Mib +   1.954Mib =  27.294Mib
   10:     Memory:         code +      nodes +       misc =      total
   10:                22.242Mib +   9.536Mib +   2.235Mib =  34.001Mib
    8:     Memory:         code +      nodes +       misc =      total
    8:                22.242Mib +  17.697Mib +   2.517Mib =  42.441Mib
   13:     Memory:         code +      nodes +       misc =      total
   13:                22.242Mib +  26.206Mib +   2.798Mib =  51.230Mib
   11:     Memory:         code +      nodes +       misc =      total
   11:                22.242Mib + 356.793Mib +   3.049Mib = 382.062Mib
    9:     -> Evolving tree number 4 {1}
    9:         Output tree data at t=   5.98 Gyr
    9:         Output tree data at t=  13.79 Gyr
    9:     <- Finished tree
    4:     -> Evolving tree number 2 {1}
    4:         Output tree data at t=   5.98 Gyr
    3:     -> Evolving tree number 5 {1}
    4:         Output tree data at t=  13.79 Gyr
    4:     <- Finished tree
    5:     -> Evolving tree number 3 {1}
    3:         Output tree data at t=   5.98 Gyr
    5:         Output tree data at t=   5.98 Gyr
    5:         Output tree data at t=  13.79 Gyr
    5:     <- Finished tree
    3:         Output tree data at t=  13.79 Gyr
    3:     <- Finished tree
   15:     -> Evolving tree number 1 {1}
   15:         Output tree data at t=   5.98 Gyr
   15:         Output tree data at t=  13.79 Gyr
   15:     <- Finished tree
    0:     -> Evolving tree number 6 {1}
    0:         Output tree data at t=   5.98 Gyr
    0:         Output tree data at t=  13.79 Gyr
    0:     <- Finished tree
   MM: <- Done task: merger tree evolution

and a file ``darkMatterOnlyMergerTrees.hdf5`` will have been created.

As Galacticus runs it outputs information about its progress to the terminal. It first tells you what task it is working on:

.. code-block:: text

   MM: -> Begin task: merger tree evolution

and then informs you about some of the sub-tasks within that task:

.. code-block:: text

    2:     -> Evolving tree number 3 {1}
    2:         Output tree data at t=   5.98 Gyr
    2:         Output tree data at t=  13.79 Gyr
    2:     <- Finished tree

Note the prefixes to these lines: "MM:" and " 2:". Galacticus runs in parallel whenever possible. In this case using `OpenMP <https://www.openmp.org/>`_ parallelism - so if you have multiple cores on your computer Galacticus will attempt to use all of them. The "MM:" prefix indicates a message from the "master" thread, while a numerical prefix such as " 2:" indicates a message from one of the parallel threads. Note that messages from parallel threads can be interleaved in the output.

Understanding the input
-----------------------

We'll explore the relevant parts of the input parameter file - you can explore the full file `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/darkMatterOnlyMergerTrees.xml>`_. Sections that we've already explored in previous tutorials and which haven't changed we'll just skip over.

The first block in the parameter file is:

.. code-block:: xml

   <!-- Specify tasks to perform -->
   <task value="evolveForests"/>

As usual, we have to tell Galacticus what task we want it to perform (although, in this case we could omit this parameter as "``evolveForests``" is the default task). The "``evolveForests``" tasks tells Galacticus to generate a set of merger trees for dark matter halos, then evolve them forward in time according to whatever physics you specify and output the results.

A merger tree is made up of a collection of "nodes" - you can think of a node as representing a single dark matter halo (or subhalo) along with everything inside it, such as gas, stars, black holes, etc. Each thing that a node may contain is called a "component" - so we next specify what components we want to use:

.. code-block:: xml

   <!-- Component selection -->
   <componentBasic             value="standard"/>
   <componentBlackHole         value="null"    />
   <componentDarkMatterProfile value="scale"   />
   <componentDisk              value="null"    />
   <componentHotHalo           value="null"    />
   <componentSatellite         value="standard"/>
   <componentSpheroid          value="null"    />
   <componentSpin              value="null"    />

In this tutorial we're going to run a dark matter only model, so most of the components above which describe baryonic content are set to "null" - this deactivates those components. We set three components to something other than "null":

* The "basic" component tracks the total mass of each node and the time at which the node exists - so is (almost) always needed - we choose the "standard" implementation of this component (don't worry yet about what other choices there might be).
* The "darkMatterProfile" component tracks properties of the dark matter halo density profile. We're going to use NFW profiles for our halos, and these are defined by a scale-length. So, we choose the "scale" implementation of this component to allow Galacticus to track those scale-lengths.
* The "satellite" component tracks details about the orbits of subhalos within their host halo. We choose the "standard" implementation which just tracks how long until the subhalo merges with its host halo, ignoring any other details.

We next specify cosmological parameters - we've already seen this in previous tutorials, but there's a slight difference:

.. code-block:: xml

   <!-- Cosmological parameters and options -->
   <cosmologyFunctions  value="matterLambda"/>
   <cosmologyParameters value="simple">
     <HubbleConstant  value="70.20000"/>
     <OmegaMatter     value=" 0.27250"/>
     <OmegaDarkEnergy value=" 0.72750"/>
     <OmegaBaryon     value=" 0.00000"/>
     <temperatureCMB  value=" 2.72548"/>
   </cosmologyParameters>

Note that we've set the "``OmegaBaryon``" parameter to zero since we want to run a dark matter only calculation.

Next we specify the power spectrum. We've already seen this in previous tutorials, but here we again do something slightly different:

.. code-block:: xml

    <!-- Power spectrum options -->
   <transferFunction value="eisensteinHu1999">
     <neutrinoNumberEffective value="3.046"/>
     <neutrinoMassSummed      value="0.000"/>
     <cosmologyParameters value="simple">
       <HubbleConstant  value="70.20000"/>
       <OmegaMatter     value=" 0.27250"/>
       <OmegaDarkEnergy value=" 0.72750"/>
       <OmegaBaryon     value=" 0.04550"/>
       <temperatureCMB  value=" 2.72548"/>
     </cosmologyParameters>
   </transferFunction>

Note that within the transfer function block we include a block of cosmological parameters. We've already specified cosmological parameters above, so why do so again here? Notice that here we've set a non-zero value for "``OmegaBaryon``", while above we set this parameter to zero. What's going on? Often, in dark matter-only N-body simulations the initial conditions (i.e. the density field, derived from the power spectrum) is computed assuming that the universe contains baryons (which is important since the baryons affect the power spectrum) even though the simulation itself is going to be run using only dark matter. This can be useful to capture features such as BAO in the simulation.

Here we're doing the same thing. We want to run our calculation with no baryons, but construct our power spectrum including the effects of baryons. In this case, the transfer function method we've selected (`Eisenstein & Hu 1999 <http://adsabs.harvard.edu/abs/1999ApJ...511....5E>`_) needs to know the cosmological model. It will first look for a cosmological model specified inside its own block in the parameter file. If it finds one there, it will use it. If not, it moves up one level in the parameters structure and looks again, and uses whatever cosmological parameters it finds there.

The task we've selected will generate and evolve a set of merger trees. Our next block of parameters therefore specifies what trees to create:

.. code-block:: xml

   <!-- Merger tree building options -->
   <mergerTreeConstructor value="build"   />
   <mergerTreeBuilder     value="cole2000" >
     <accretionLimit   value="0.1"/>
     <mergeProbability value="0.1"/>
   </mergerTreeBuilder>
   <mergerTreeBranchingProbability value="parkinsonColeHelly">
     <G0                 value="+0.57"/>
     <gamma1             value="+0.38"/>
     <gamma2             value="-0.01"/>
     <accuracyFirstOrder value="+0.10"/>
   </mergerTreeBranchingProbability>
   <mergerTreeBuildMasses value="sampledDistributionUniform">
     <massTreeMinimum value="1.0e10"/>
     <massTreeMaximum value="1.0e13"/>
     <treesPerDecade  value="2"     />
   </mergerTreeBuildMasses>

First we specify how to generate the trees. We select the "``build``" method, which means trees are built internally within Galacticus (other options including reading trees from a file). Then we need to specify the algorithm to use for that tree building - we go with the `Cole et al. (2000) <http://adsabs.harvard.edu/abs/2000MNRAS.319..168C>`_ algorithm and specify two parameter values that control the numerical precision of the algorithm. That method requires some description of the branching probabilities in trees (i.e. the merger rates of dark matter halos) - we select the `Parkinson, Cole, & Helly (2008) <http://adsabs.harvard.edu/abs/2008MNRAS.383..557P>`_ branching probabilities, and choose parameter values for that algorithm taken from their paper (we also specify the "``accuracyFirstOrder``" parameter which controls the numerical precision of the algorithm).

Finally in this section we need to specify the *z*\ =0 masses of the trees and how many trees to create. There are many ways to do this, but we're going to use a uniform distribution of tree masses sampled from a distribution (the default distribution is just the halo mass function). We specify the minimum and maximum tree masses to use (masses are always in units of Solar masses, :math:`M_\odot` - **not** :math:`h^{-1}M_\odot`). We also specify ``treesPerDecade``\ =2 - this is the number of trees per decade of halo mass to generate. In this case we have :math:`\log_{10}(1.0\mathrm{e}13/1.0\mathrm{e}10)=3` decades of halo mass, so we'll get a total of 6 trees.

The next block of parameters specifies how Galacticus should handle the hierarchy of substructures that form as halos merge. We'll choose "``singleLevelHierarchy``" (currently the only option) which means that we don't allow for any levels of the hierarchy deeper than the first - so there are substructures, but not sub-substructures etc.

.. code-block:: xml

   <!-- Substructure hierarchy options -->
   <mergerTreeNodeMerger value="singleLevelHierarchy"/>

We next specify how the structure of our dark matter halos is to be determined. We'll use `NFW <http://adsabs.harvard.edu/abs/1997ApJ...490..493N>`_ density profiles, and choose their concentrations using the `Gao et al. (2008) <http://adsabs.harvard.edu/abs/2008MNRAS.387..536G>`_ fitting function. We'll also set a limit on the concentration so that it can never go below 4 (which is unphysical, but can happen if the fitting function is applied far outside its range of validity):

.. code-block:: xml

   <!-- Dark matter halo structure options -->
   <darkMatterProfileDMO            value="NFW"    />
   <darkMatterProfileConcentration  value="gao2008"/>
   <darkMatterProfileMinimumConcentration value="4"      />

Since we want to run a dark matter-only calculation we need to switch off some baryonic physics that is on by default. Specifically we set the distribution of hot gas in halos to "``null``".

.. code-block:: xml

   <!-- Switch off baryonic physics -->
   <hotHaloMassDistribution value="null"/>

We've chosen how to build our merger trees, next we need to choose how to evolve them. We're just going to choose the standard evolvers for trees and nodes, and choose some reasonable values for the numerical precision parameters:

.. code-block:: xml

   <!-- Tree evolution -->
   <mergerTreeNodeEvolver value="standard">
     <odeToleranceAbsolute value="0.01"/>
     <odeToleranceRelative value="0.01"/>
   </mergerTreeNodeEvolver>
   <mergerTreeEvolver value="standard">
     <timestepHostAbsolute value="1.0"/>
     <timestepHostRelative value="0.1"/>
   </mergerTreeEvolver>

Finally we specify what we want to be output. We'll choose a file name for the output, specify the redshifts at which we want to output the properties of nodes in the merger tree, and finally choose the outputter itself - we'll use the standard outputter which just writes lists of node properties to the file:

.. code-block:: xml

   <!-- Output options -->
   <galacticusOutputFileName value="darkMatterOnlyMergerTrees.hdf5"/>
   <outputTimes value="list">
     <redshifts value="0.0 1.0"/>
   </outputTimes>
   <mergerTreeOutputter value="standard">
     <outputReferences value="false"/>
   </mergerTreeOutputter>

Understanding the output
------------------------

We'll explore the output file, ``darkMatterOnlyMergerTrees.hdf5`` using the command-line tools `h5ls <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Ls>`_ and `h5dump <https://docs.hdfgroup.org/archive/support/HDF5/doc/RM/Tools.html#Tools-Dump>`_.

The structure of the output file is similar to what we've seen in previous tutorials, so let's take a look at the *z*\ =0 output group:

.. code-block:: console

   $ h5ls darkMatterOnlyMergerTrees.hdf5
   $ h5ls darkMatterOnlyMergerTrees.hdf5/Outputs/Output2
   mergerTreeCount          Dataset {6/Inf}
   mergerTreeIndex          Dataset {6/Inf}
   mergerTreeSeed           Dataset {6/Inf}
   mergerTreeStartIndex     Dataset {6/Inf}
   mergerTreeWeight         Dataset {6/Inf}
   nodeData                 Group

This output group now contains five datasets, plus another group. Let's look at the datasets first. The "``{6/Inf}``" after each dataset means that they each contain 6 entries - one for each of our 6 merger trees. Let's take a look at one of the datasets:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/mergerTreeIndex
   mergerTreeIndex          Dataset {6/Inf}
       Data:
           (0) 4, 1, 2, 3, 5, 6

Each merger tree generated by Galacticus is assigned an index - in this case they just run from 1 to 6. The "``mergerTreeIndex``" dataset simply tells you the index of each tree written to the file. Note that they're not in order - trees are output in the order that they finish being computed - and when running in parallel that order can change depending on how fast each thread works.

We can also look at the "``mergerTreeWeight``" dataset:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/mergerTreeWeight
   mergerTreeWeight         Dataset {6/Inf}
       Data:
           (0) 0.0190718657104551, 0.01883711019377, 0.0188628587556674, 0.0189205005044752, 0.0196372577870537, 0.0508381724391334

This dataset tells you the number of each tree per unit volume (in units of :math:`\mathrm{Mpc}^{-3}`) of the universe that you'd need to have to construct a fair sample of halos. So, if you wanted to construct a halo mass function for example you could just bin up the halo masses, weighting by these values. The weights will depend on how you choose to sample halo masses.

The "``mergerTreeSeed``" dataset records the seed value used by each tree when generating random numbers - this can be useful to reproduce the same tree again in a different run.

We'll come back to the remaining two datasets soon, but for now let's look at the "``nodeData``" group:

.. code-block:: console

   $ h5ls darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/nodeData
   basicMass                Dataset {7/Inf}
   basicTimeLastIsolated    Dataset {7/Inf}
   darkMatterProfileScale   Dataset {7/Inf}
   nodeIndex                Dataset {7/Inf}
   nodeIsIsolated           Dataset {7/Inf}
   parentIndex              Dataset {7/Inf}
   satelliteBoundMass       Dataset {7/Inf}
   satelliteIndex           Dataset {7/Inf}
   satelliteMergeTime       Dataset {7/Inf}
   siblingIndex             Dataset {7/Inf}

This group contains several datasets, each of which contains 7 elements. Each of those 7 elements corresponds to a node in a merger tree. Why 7 elements when we have just 6 trees? A tree contains multiple nodes, so the number of nodes here can exceed the number of trees. In this small example, only one tree happens to contain more than 1 node, but in a realistic case you could expect there to be many more nodes than trees.

Nodes from all 6 trees are stored together in these datasets. So, how do you know which nodes belong to which trees? We can figure this out using the "``mergerTreeStartIndex``" and "``mergerTreeCount``" datasets:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/mergerTreeStartIndex
   mergerTreeStartIndex     Dataset {6/Inf}
       Data:
           (0) 0, 1, 2, 3, 4, 5
   [abenson@mies galacticus]$ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/mergerTreeCount
   mergerTreeCount          Dataset {6/Inf}
       Data:
           (0) 1, 1, 1, 1, 1, 2

The values in ``mergerTreeStartIndex`` tell us at what position in the ``nodeData`` datasets the nodes for each tree begin at. The values in the ``mergerTreeCount`` dataset tell us how many nodes are output for each tree. Nodes for a tree are always stored contiguously in the ``nodeData`` datasets, so you can easily figure out which entries correspond to a given tree. For example, suppose we want to know the nodes corresponding to the tree with index 3. From above we know that this tree was the 4th tree output. So, according to ``mergerTreeStartIndex`` its nodes begin at location 3, and there's only one node. So, we just pull out location 3 from the ``nodeData`` datasets (note that they are indexed starting at zero). For tree index 6 we would find that the nodes start at location 5, and there are two nodes. So, we'd pull out locations 5-6 from the ``nodeData`` datasets.

Let's go through the datasets in the ``nodeData`` group:

.. code-block:: console

   $ h5ls -d darkMatterOnlyMergerTrees.hdf5/Outputs/Output2/nodeData
   basicMass                Dataset {7/Inf}
       Data:
           (0) 26294735584.2762, 11008079278.5949, 13737902820.2165, 18132446384.304, 46214619906.483, 7557240227.04057, 154165801427.933
   basicTimeLastIsolated    Dataset {7/Inf}
       Data:
           (0) 13.7915100007606, 13.7915100007606, 13.7915100007606, 13.7915100007606, 13.7915100007606, 9.60061103648177, 13.7915100007606
   darkMatterProfileScale   Dataset {7/Inf}
       Data:
           (0) 0.00451184258191652, 0.00299512345054255, 0.00332410700870255, 0.00378791722462339, 0.00588278500721434, 0.00282708687493467, 0.0103690332192463
   nodeIndex                Dataset {7/Inf}
       Data:
           (0) 1, 1, 1, 1, 1, 3, 1
   nodeIsIsolated           Dataset {7/Inf}
       Data:
           (0) 1, 1, 1, 1, 1, 0, 1
   parentIndex              Dataset {7/Inf}
       Data:
           (0) -1, -1, -1, -1, -1, 1, -1
   satelliteBoundMass       Dataset {7/Inf}
       Data:
           (0) 26294735584.2762, 11008079278.5949, 13737902820.2165, 18132446384.304, 46214619906.483, 7557240227.04057, 154165801427.933
   satelliteIndex           Dataset {7/Inf}
       Data:
           (0) -1, -1, -1, -1, -1, -1, 3
   satelliteMergeTime       Dataset {7/Inf}
       Data:
           (0) -1, -1, -1, -1, -1, 4.09827855235655, -1
   siblingIndex             Dataset {7/Inf}
       Data:
           (0) -1, -1, -1, -1, -1, -1, -1

There are several datasets that end with "``Index``" - these are useful to see how nodes are related to each other:

* ``nodeIndex`` - each node in each merger tree is assigned an index which is unique (at least within the tree). You'll see that all but one entry in this dataset is "1" - corresponding to the *z*\ =0 halo. For the final tree output there's a second node, with index "3" - this is a subhalo within the *z*\ =0 halo.
* ``parentIndex`` - for subhalos this is the index of their host halo, for other halos it corresponds to the index of the halo that they will merge with in the future. Since this is the *z*\ =0 output there is no future for this tree, so values of "-1" are listed instead.
* ``satelliteIndex`` - if a halo has any subhalos then this dataset will give the index of the first subhalo.
* ``siblingIndex`` - for halos which have siblings (other halos with which they will eventually merge) the first sibling is given by this dataset; for subhalos if there is more than one subhalo in the same halo this dataset gives the index of the next subhalo.

From these you can figure out which nodes are subhalos. To make this easier you can also use the ``nodeIsIsolated`` dataset - it has a value of 1 for halos, and 0 for subhalos.

Each component of each node can also output some properties. For example, the ``basic`` component has output

* ``basicMass`` - the total mass of the node
* ``basicTimeLastIsolated`` - the time at which this node was last an isolated node - i.e. when it was last a halo and not a subhalo. For nodes that are still halos at the output time this will correspond to the current time. For subhalos, this is the time at which they first merged into a larger halo.

The ``satellite`` component outputs:

* ``satelliteMergeTime`` - this is the time (in units of Gyr) until the subhalo will merge with its host halo. For nodes that are not subhalos this is set to -1 instead.

There are many more properties that can be output if you request them.
