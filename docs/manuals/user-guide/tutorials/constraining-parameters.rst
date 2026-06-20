Constraining Galacticus Parameters
==================================

This tutorial guides you through constraining the parameters of a Galacticus model to achieve a good fit to an observational dataset. Galacticus has built-in `MCMC <https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo>`_ functionality which can be used for this purpose.

Here we will use a very simple example - we'll constrain a single parameter of a Galacticus model to obtain a match to a single point in the stellar mass-halo mass relation of `Leauthaud et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...744..159L>`_. Much more complex cases are possible (multiple parameters, coupled parameters, multiple target datasets etc.), but this simple example will illustrate the key features.

Running the simulation
----------------------

For this tutorial we require Galacticus to be compiled with `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ parallelism. To do this:

.. code-block:: bash

   make -j8 GALACTICUS_BUILD_OPTION=MPI Galacticus.exe

Note that you must have MPI installed for this to work.

Once compilation is completed, to run the tutorial model:

.. code-block:: bash

   export OMP_NUM_THREADS=1
   mpirun -np 4 Galacticus.exe parameters/tutorials/mcmcConfig.xml

The ``export OMP_NUM_THREADS=1`` effectively switches off `OpenMP <https://en.wikipedia.org/wiki/OpenMP>`_ parallelism (which we don't want to use for this tutorial). The ``mpirun -np 4`` prefix command launches 4 parallel Galacticus processes which will communicate via MPI.

Expect this example to run for around 15 minutes - it should output something like this:

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

   1:M: -> Begin task: posterior sampling
   1:M:     Process 0001 [PID: 13486] is running on host 'node081'
   2:M: -> Begin task: posterior sampling
   2:M:     Process 0002 [PID: 13487] is running on host 'node081'
   3:M: -> Begin task: posterior sampling
   3:M:     Process 0003 [PID: 13488] is running on host 'node081'
   0:M: -> Begin task: posterior sampling
   0:M:     Process 0000 [PID: 13485] is running on host 'node081'
   0:M:     -> Load balancing report
   0:M:         Not performing load balancing - missing work cost data
   0:M:     <- done
   0:M:     Chain 0 has logℒ= -0.8306E+00
   1:M:     Chain 1 has logℒ= -0.4578E+01
   2:M:     Chain 2 has logℒ= -0.3024E+02
   3:M:     Chain 3 has logℒ= -0.1131E+01
   0:M:     -> Node work done vs. expected:
   0:M:         Node 0000: work (actual/estimated) =     10.65 /     -1.00
   0:M:         Node 0001: work (actual/estimated) =      8.99 /     -1.00
   0:M:         Node 0002: work (actual/estimated) =      8.70 /     -1.00
   0:M:         Node 0003: work (actual/estimated) =      8.07 /     -1.00
   0:M:     <- done
   0:M:     -> Load balancing report
   0:M:         -> Chain redistribution:
   0:M:             Chain 0000 -> process/node 0000/0001 (work =     10.65)
   0:M:             Chain 0001 -> process/node 0001/0001 (work =      8.99)
   0:M:             Chain 0002 -> process/node 0002/0001 (work =      8.70)
   0:M:             Chain 0003 -> process/node 0003/0001 (work =      8.07)
   0:M:         <- done
   0:M:         -> Node work loads:
   0:M:             Node 0001: work =     36.40
   0:M:         <- done
   0:M:     <- done
   0:M:     Chain 0 has logℒ= -0.1154E+01
   1:M:     Chain 1 has logℒ= -0.3044E+02
   2:M:     Chain 2 has logℒ= -0.4176E+01
   3:M:     Chain 3 has logℒ= -0.3042E+02
   0:M:     -> Node work done vs. expected:
   0:M:         Node 0000: work (actual/estimated) =      8.01 /     10.65
   0:M:         Node 0001: work (actual/estimated) =      7.93 /      8.99
   0:M:         Node 0002: work (actual/estimated) =      8.71 /      8.70
   0:M:         Node 0003: work (actual/estimated) =      7.70 /      8.07
   0:M:     <- done

   .
   .
   .
   0:M:     Converged after 10 steps
   0:M:     -> Load balancing report
   0:M:         -> Chain redistribution:
   0:M:             Chain 0000 -> process/node 0002/0001 (work =      7.84)
   0:M:             Chain 0001 -> process/node 0003/0001 (work =      7.81)
   0:M:             Chain 0002 -> process/node 0001/0001 (work =      8.14)
   0:M:             Chain 0003 -> process/node 0000/0001 (work =      9.10)
   0:M:         <- done
   0:M:         -> Node work loads:
   0:M:             Node 0001: work =     32.87
   0:M:         <- done
   0:M:     <- done
   0:M:     Chain 3 has logℒ= -0.8918E+00
   1:M:     Chain 2 has logℒ= -0.3917E+00
   2:M:     Chain 0 has logℒ= -0.9008E+00
   3:M:     Chain 1 has logℒ= -0.2899E+01
   0:M:     -> Node work done vs. expected:
   0:M:         Node 0000: work (actual/estimated) =      8.41 /      7.84
   0:M:         Node 0001: work (actual/estimated) =      8.42 /      7.81
   0:M:         Node 0002: work (actual/estimated) =      8.53 /      8.14
   0:M:         Node 0003: work (actual/estimated) =      8.55 /      9.10
   0:M:     <- done
   .
   .
   .
   0:M:     -> Load balancing report
   0:M:         -> Chain redistribution:
   0:M:             Chain 0000 -> process/node 0001/0001 (work =      8.85)
   0:M:             Chain 0001 -> process/node 0000/0001 (work =     13.86)
   0:M:             Chain 0002 -> process/node 0002/0001 (work =      8.22)
   0:M:             Chain 0003 -> process/node 0003/0001 (work =      8.06)
   0:M:         <- done
   0:M:         -> Node work loads:
   0:M:             Node 0001: work =     39.00
   0:M:         <- done
   0:M:     <- done
   0:M:     Chain 1 has logℒ= -0.3850E+01
   1:M:     Chain 0 has logℒ= -0.1265E+02
   2:M:     Chain 2 has logℒ= -0.1516E-01
   3:M:     Chain 3 has logℒ= -0.3263E+02
   0:M:     -> Node work done vs. expected:
   0:M:         Node 0000: work (actual/estimated) =      7.64 /      8.85
   0:M:         Node 0001: work (actual/estimated) =      9.02 /     13.86
   0:M:         Node 0002: work (actual/estimated) =      9.03 /      8.22
   0:M:         Node 0003: work (actual/estimated) =      8.44 /      8.06
   0:M:     <- done
   0:M: <- Done task: posterior sampling
   2:M: <- Done task: posterior sampling
   1:M: <- Done task: posterior sampling
   3:M: <- Done task: posterior sampling

where we have cut out some of the output of chain likelihoods for brevity.

The output starts by reporting on where each MPI process is running and lists the PID for that process (which can be useful for debugging). Since we're running four MPI processes we get 4 MCMC chains, which begin to evaluate our model for different combinations of parameters. After each evaluation you'll see a report on the likelihood of the model for that chain:

.. code-block:: text

   0:M:     Chain 0 has logℒ= -0.8306E+00

The simulation will continue evaluating models until it completes. In between this chain output you'll see various reports, for example a report on run time for each chain:

.. code-block:: text

   0:M:     -> Node work done vs. expected:
   0:M:         Node 0000: work (actual/estimated) =      8.41 /      7.84
   0:M:         Node 0001: work (actual/estimated) =      8.42 /      7.81
   0:M:         Node 0002: work (actual/estimated) =      8.53 /      8.14
   0:M:         Node 0003: work (actual/estimated) =      8.55 /      9.10
   0:M:     <- done

which reports how much CPU time was spent on each process vs. how much was estimated would be used. This is used in load balancing, which is also reported:

.. code-block:: text

   0:M:     -> Load balancing report
   0:M:         -> Chain redistribution:
   0:M:             Chain 0000 -> process/node 0002/0001 (work =      7.84)
   0:M:             Chain 0001 -> process/node 0003/0001 (work =      7.81)
   0:M:             Chain 0002 -> process/node 0001/0001 (work =      8.14)
   0:M:             Chain 0003 -> process/node 0000/0001 (work =      9.10)
   0:M:         <- done
   0:M:         -> Node work loads:
   0:M:             Node 0001: work =     32.87
   0:M:         <- done
   0:M:     <- done

which reports which chain was computed on which node and how much time total was used by each node.

Understanding the input parameter file
---------------------------------------

You can view the complete input parameter `file <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/mcmcConfig.xml>`_ for this tutorial here. Here we'll focus on each section of the parameter file and understand what it does:

.. code-block:: xml

   <task value="posteriorSample">
     <initializeNodeClassHierarchy value="false"/>
   </task>

We begin by specifying the "task" to perform - we choose ``posteriorSample`` which causes Galacticus to run a posterior sampling simulation, which will generate a set of parameters sampled from the posterior distribution given some constraining datasets.

.. code-block:: xml

   <posteriorSampleLikelihood value="galaxyPopulation">
     <baseParametersFileName   value="parameters/tutorials/mcmcBase.xml"/>
     <failedParametersFileName value="./failedParameters.xml"           />
     <randomize                value="false"                            />
     <evolveForestsVerbosity   value="silent"                           />
   </posteriorSampleLikelihood>

In this section we specify how the likelihood function for the model should be computed. Many different likelihood functions can be used, but here we use the ``galaxyPopulation`` likelihood function - which causes Galacticus to run a model generating a population of galaxies, computed some observables from that population, and then find the likelihood of that model given some target datasets. We specify a ``baseParametersFileName`` which gives the path to a parameter file describing the galaxy population model to run - this will be described in more detail below (note that a set of `parameter change files <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/parameter-files.html#changing-parameters>`_ can be applied to this base parameter file before it is used - these can be specified as a space-separated list via the ``changeParametersFileNames`` parameter added to the above section). We also specify a path for the output of failed parameter sets - sometimes while exploring the model parameter space a model will fail to complete - its parameters will be output to this file for further exploration. The ``randomize`` option specifies whether or not to change the random number seed for the model on each evaluation. Finally, we set ``evolveForestsVerbosity`` to ``silent`` to switch off output of logging messages while running the galaxy population calculation (this just keeps the amount of output reasonable).

.. code-block:: xml

   <!-- MCMC -->
   <posteriorSampleSimulation value="differentialEvolution">
     <stepsMaximum           value="1000"      />
     <acceptanceAverageCount value="  10"      />
     <stateSwapCount         value=" 100"      />
     <logFileRoot            value="mcmcChains"/>
     <reportCount            value="  10"      />
     <sampleOutliers         value="false"     />
     <logFlushCount          value="   1"      />

This next section controls the type of posterior sampling simulation that we want to perform. Here we choose to use the differential evolution MCMC algorithm (which runs multiple MCMC chains in parallel). Note that this section is not closed at this point - the sections below are all therefore subsections of this section.

.. code-block:: xml

   <posteriorSampleState value="correlation">
     <acceptedStateCount value="100"/>
   </posteriorSampleState>

Next we specify the type of object used to store the state of each chain - here we choose a type which can store and compute the correlation structure of each chain's state.

.. code-block:: xml

   <posteriorSampleStateInitialize value="latinHypercube">
     <maximinTrialCount value="100"/>
   </posteriorSampleStateInitialize>

Each MCMC chain must be assigned some initial position in the model parameter space. He we choose to assign those positions using a Latin hypercube design, in which we generate 100 trial Latin hypercubes and keep the hypercube with the maximum minimum distance between chains - this ensures a good coverage of the parameter space.

.. code-block:: xml

   <posteriorSampleConvergence value="gelmanRubin">
     <thresholdHatR              value=" 1.30"/>
     <burnCount                  value="10"   />
     <testCount                  value="10"   />
     <outlierCountMaximum        value=" 1"   />
     <outlierSignificance        value=" 0.95"/>
     <outlierLogLikelihoodOffset value="60"   />
     <reportCount                value=" 1"   />
     <logFileName                value="mcmcConvergence.log"/>
   </posteriorSampleConvergence>

To judge the convergence of our simulation we choose here to use the Gelman-Rubin statistic - specifying the threshold value for convergence to be declared, how many steps to burn before checking convergence, and how often to test convergence. There are also options to control detection of outlier chains.

.. code-block:: xml

   <posteriorSampleStoppingCriterion value="stepCount">
     <stopAfterCount value="10"/>
   </posteriorSampleStoppingCriterion>

This option controls when the simulation will be stopped. Here we choose a simple example which stops the simulation after a specified number of steps post-convergence.

.. code-block:: xml

   <posteriorSampleDffrntlEvltnRandomJump   value="adaptive"/>

In generating proposed new states in the differential evolution MCMC method a random component is added to each chain to ensure that the parameter space is fully converged. This option choose to use the ``adaptive`` method for this random component, which scales the size of the random proportion based on the current chain distribution.

.. code-block:: xml

   <posteriorSampleDffrntlEvltnProposalSize value="adaptive" >
     <gammaInitial          value="0.500e+0"/>
     <gammaAdjustFactor     value="1.100e+0"/>
     <gammaMinimum          value="1.000e-4"/>
     <gammaMaximum          value="3.000e+0"/>
     <acceptanceRateMinimum value="0.100e+0"/>
     <acceptanceRateMaximum value="0.900e+0"/>
     <updateCount           value="10"     />
   </posteriorSampleDffrntlEvltnProposalSize>

In the differential evolution MCMC algorithm new proposed states are generated by taking some fraction, :math:`\gamma`, of the difference of two randomly selected states and adding that to the current state. The choice of :math:`\gamma` is controlled by the above section. We chose the ``adaptive`` option which will adjust the size of :math:`\gamma` to ensure that the acceptance rate of proposed states remains between a given minimum and maximum.

.. code-block:: xml

   <!-- Feedback -->
   <modelParameter value="active">
     <name value="nodeOperator/nodeOperator[@value='stellarFeedbackDisks']/stellarFeedbackOutflows/stellarFeedbackOutflows/velocityCharacteristic"/>
     <distributionFunction1DPrior value="uniform">
       <limitLower value="25.0"/>
       <limitUpper value="500.0"/>
     </distributionFunction1DPrior>
     <operatorUnaryMapper value="identity"/>
     <distributionFunction1DPerturber value="cauchy">
       <median value="0.0"/>
       <scale value="1.0e-3"/>
     </distributionFunction1DPerturber>
   </modelParameter>

Most importantly, we have to specify the model parameters that we wish to constrain. Here we use just one, but any number can be included. We choose an ``active`` parameter (one which will be adjusted to obtain the best match to the target data). The ``name`` elements specifies the parameter in the base parameter file to vary, as we will describe below. We must specify a prior for this function - here we choose a uniform prior between given lower and upper limits. We also specify a "mapping" for the parameter - here we use an identity mapping, but we could also use a logarithmic mapping such that the internal calculations would be performed with the log of the parameter value. Finally, we must specify a distribution function for the random perturbation added to each proposal generated for this parameter - here we use a Cauchy distribution which typically produces a small perturbation, but has a tail to arbitrarily large perturbations.

.. code-block:: xml

   </posteriorSampleSimulation>

Finally, we close the ``posteriodSampleSimulation`` section.

.. code-block:: xml

   <!-- Random seed -->
   <randomNumberGenerator value="GSL">
     <seed          value="219" />
     <mpiRankOffset value="true"/>
   </randomNumberGenerator>

One last thing - since the MCMC algorithm is a Monte Carlo algorithm it needs to generate random numbers. Here we choose a random number generator (we use the `GSL <https://www.gnu.org/software/gsl/>`_ standard algorithm), and explicitly request that the random number seeds be offset for each different MPI process rank - this ensures that each of our MPI processes generates a different sequence of random numbers.

The "Base" Parameter File
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the ``mcmcConfig.xml`` parameter file that we examined above, a second parameter file is used in this tutorial. As explained above when discussing the ``posteriorSampleLikelihood`` section, we supply a parameter file which defines the "base" model, which will then have its parameters varied to attempt to find a good fit to the target dataset. In this example that base parameter file is `parameters/tutorials/mcmcBase.xml <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/mcmcBase.xml>`_. This is just a normal Galacticus parameter file, using the ``evolveForests`` task, that you could run directly to generate a model if you wanted to. The only important feature of it is that it contains each parameter which will be varied as specified in the ``mcmcConfig.xml`` file. In this case it contains a section:

.. code-block:: xml

   <!-- Node evolution and physics -->
   <nodeOperator value="multi">
     <!-- Star formation options -->
     <nodeOperator value="starFormationDisks"    >
       <luminositiesStellarInactive value="true"/>
     </nodeOperator>
     <nodeOperator value="starFormationSpheroids">
       <luminositiesStellarInactive value="true"/>
     </nodeOperator>
     <!--Stellar feedback outflows-->
     <nodeOperator value="stellarFeedbackDisks">
       <stellarFeedbackOutflows value="rateLimit">
         <timescaleOutflowFractionalMinimum value="0.001"/>
         <stellarFeedbackOutflows value="powerLaw">
           <velocityCharacteristic value="175.0"/>  <!-- This is the parameter being varied -->
           <exponent               value="  3.5"/>
         </stellarFeedbackOutflows>
       </stellarFeedbackOutflows>
     </nodeOperator>
     <nodeOperator value="stellarFeedbackSpheroids">
       <stellarFeedbackOutflows value="rateLimit">
         <timescaleOutflowFractionalMinimum value="0.001"/>
         <stellarFeedbackOutflows value="powerLaw">
           <velocityCharacteristic value=" 50.0"/>
           <exponent               value="  3.5"/>
         </stellarFeedbackOutflows>
       </stellarFeedbackOutflows>
     </nodeOperator>
   </nodeOperator>

In the ``mcmcConfig.xml`` file we referred to a parameter "``nodeOperator/nodeOperator[@value='stellarFeedbackDisks']/stellarFeedbackOutflows/stellarFeedbackOutflows/velocityCharacteristic``". This refers to the parameter indicated (by a comment) in the above - we first find the ``nodeOperator`` element in ``mcmcBase.xml``. Next we look for the ``nodeOperator`` element within that element which has a ``value`` attribute equal to ``stellarFeedbackDisks`` (we can also reference by index, e.g. ``nodeOperator[3]`` would reference this same element as it is the third ``nodeOperator`` in this list - element as our indexing starts at 1 as is standard for `XPath <http://www.w3schools.com/Xml/xpath_syntax.asp>`_ indexing). Next we find the ``stellarFeedbackOutflows`` element inside that element, then the ``stellarFeedbackOutflows`` inside that element, and finally locate the ``velocityCharacteristic`` element inside that element. It is the value of this parameter which will be varied.

Of course, you can specify any other parameters in the base parameter file also - their values will be fixed throughout the MCMC simulation.

The other key part of the base parameter file is in defining the target datasets to be used to constrain the model. These make use of Galacticus' ability to compute predictions for observables as it runs. The base parameter file contains a section:

.. code-block:: xml

   <!-- Analyses -->
   <outputAnalysis value="stellarVsHaloMassRelationLeauthaud2012" >
     <redshiftInterval                     value="1"      />
     <computeScatter                       value="false"  />
     <systematicErrorPolynomialCoefficient value="0.0 0.0"/>
     <likelihoodBins                       value="9"      />
   </outputAnalysis>

This tells Galacticus to compute its prediction for the stellar mass-halo mass relation of `Leauthaud et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...744..159L>`_ - specifically in the first redshift interval (`Leauthaud et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...744..159L>`_ compute the relation in three redshift intervals), and to compute the likelihood of the model given this dataset using only halo mass bin number 11 (since for this simple tutorial we only compute halos falling within that bin). The ``computeScatter`` option is set to ``false`` such that we compute the mean of the relation - if it were instead set to ``true`` the scatter in the relation would be computed instead. The ``systematicErrorPolynomialCoefficient`` parameter allows for the possibility of including and constraining a model of observational systematic errors here, but we ignore it for now (setting the values of this parameter to zero).

Understanding the output
------------------------

The results of the MCMC simulation are output to files named ``mcmcChains_????.log`` where ``????`` corresponds to the MPI process number and each file holds the output of states from that chain. An example of the content of each file is:

.. code-block:: text

              1           0   8.65868378     F  -7.3171034777503241       -1.1537886737156828        321.98163229607087
              2           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              3           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              4           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              5           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              6           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              7           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              8           0   8.96963501     F  -7.3564114151060496       -1.1930966110714085        322.62448285700106
              9           0   8.95562744     F  -7.3727669523445538       -1.2094521483099128        322.80007092398876
             10           0   8.95562744     F  -7.3727669523445538       -1.2094521483099128        322.80007092398876
             11           0   8.36373901     T  -7.0640684583132103      -0.90075365427856924        318.64017850185974
             12           0   8.06176758     T  -7.1758380406243720       -1.0125232365897303        324.84524606237727
             13           0   8.06176758     T  -7.1758380406243720       -1.0125232365897303        324.84524606237727
             14           0   8.76565552     T  -7.3388034666121991       -1.1754886625775576        322.42077843550794
             15           0   8.08679199     T  -6.7269143477283295      -0.56359954369368803        314.22563572091480
             16           0   8.08679199     T  -6.7269143477283295      -0.56359954369368803        314.22563572091480
             17           0   7.77185059     T  -7.1396858801635883      -0.97637107612894714        319.78935277370073
             18           0   8.35467529     T  -6.5888317014128344      -0.42551689737819348        306.82365204243695
             19           0   7.88580322     T  -6.2438500587941803       -8.0535254759539335E-002   292.61032842104959
             20           0   7.88580322     T  -6.2438500587941803       -8.0535254759539335E-002   292.61032842104959
             21           0   7.88580322     T  -6.2438500587941803       -8.0535254759539335E-002   292.61032842104959

Each line corresponds to one step in the MCMC chain. The columns are as follows:

#. The step number;
#. The chain number;
#. The CPU time taken to evaluate this model;
#. T/F to indicate whether the simulation is converged or not;
#. The logarithm of the posterior probability of this state;
#. The logarithm of the likelihood of this state;
#. (and any further columns) The values of the model parameters in this state.

Analyzing the chains
--------------------

Once the MCMC simulation has completed (or while it is still running) you will likely want to analyze the resulting chains - for example to assess convergence, extract a maximum-likelihood model, or visualize the posterior distribution. This is supported by the `Dendros <https://github.com/galacticusorg/dendros>`_ package, the companion analysis and visualization package for Galacticus, which provides:

* Convergence diagnostics (e.g. Gelman-Rubin statistic, acceptance rate);
* Maximum-likelihood model extraction;
* Posterior predictive checks;
* Kernel density estimates and corner plots of the posterior distribution.

Dendros is available on `PyPI <https://pypi.org/project/dendros/>`_:

.. code-block:: bash

   pip install dendros

See the `Dendros documentation <https://github.com/galacticusorg/dendros>`_ for usage examples.
