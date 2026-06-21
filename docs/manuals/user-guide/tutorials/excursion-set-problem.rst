Solving the Excursion Set Problem
=================================

In this tutorial we'll use Galacticus to compute and output solutions to the excursion set problem, which is an important ingredient in Press-Schechter-like halo mass function calculations. Briefly, the problem is to consider a set of random walk trajectories in overdensity, :math:`\delta`, as a function of variance (of the cosmological density field), :math:`S`, and then compute the fraction, :math:`f(S)\mathrm{d}S`, of these trajectories which make their first upcrossing through a barrier, :math:`B(S)`, between :math:`S` and :math:`S+\mathrm{d}S`.

This tutorial assumes that you've previously completed some of the earlier tutorials, so are familiar with the basics of Galacticus usage. If you haven't, a good place to start is with the `halo mass function <https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Dark-matter-halo-mass-function>`_ tutorial.

A more mathematical and detailed discussion of the excursion set problem, and how it is implemented in Galacticus can be found `here <https://hackmd.io/@galacticus/r1gOkpTat>`_.

Running the calculation
-----------------------

To run this tutorial:

.. code-block:: console

   $ ./Galacticus.exe parameters/tutorials/excursionSets.xml

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
   M: -> Begin task: excursion sets
   M: <- Done task: excursion sets

and a file ``excursionSets.hdf5`` will have been created.

Understanding the input
-----------------------

You can look at the entire parameter file for this tutorial `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/excursionSets.xml>`_. Below we'll explore just those sections which are specific to this tutorial.

.. code-block:: xml

   <!-- Specify tasks to perform -->
   <task value="excursionSets"/>

This block tells Galacticus what task we want it to perform - here we're just asking it to solve the excursion set problem.

.. code-block:: xml

   <excursionSetBarrier       value="criticalOverdensity"/>
   <excursionSetFirstCrossing value="linearBarrier"      />

The first option specifies the form of the barrier function, :math:`B(S)`. Here we choose :math:`B(S)=\delta_\mathrm{c}` where :math:`\delta_\mathrm{c}` is the critical overdensity for collapse of a dark matter halo.

The second option selects the solver to use to compute the distribution of first crossings of the barrier. Here we choose a solver which utilizes the analytic solution for a linear barrier, i.e. a barrier of the form :math:`B(S) = a + b S`.

Understanding the output
------------------------

The output file ``excursionSets.hdf5`` has the excursion set problem solutions in the ``excursionSets`` group:

.. code-block:: console

   $ h5ls excursionSets.hdf5/excursionSets
   barrier                  Dataset {15, 51}
   firstCrossingProbability Dataset {15, 51}
   firstCrossingRate        Dataset {15, 51, 51}
   mass                     Dataset {51}
   massFunction             Dataset {15, 51}
   powerSpectrum            Dataset {15, 51}
   time                     Dataset {15}
   variance                 Dataset {15, 51}
   wavenumber               Dataset {51}

These datasets contain the following information:

* ``mass``: Halo mass [:math:`\mathrm{M}_\odot`];
* ``time``: Cosmic time [Gyr];
* ``wavenumber``: Wavenumber corresponding to this halo mass [:math:`\mathrm{Mpc}^{-1}`];
* ``powerSpectrum``: Power spectrum at this wavenumber [:math:`\mathrm{Mpc}^3`];
* ``variance``: The variance, :math:`S(M) = \sigma^2(M)`, at this halo mass;
* ``barrier``: The excursion set barrier, :math:`B(S)`;
* ``firstCrossingProbability``: The probability of first crossing this barrier between :math:`S` and :math:`S+\mathrm{d}S`;
* ``firstCrossingRate``: The rate of first crossing of the barrier per unit time [:math:`\mathrm{Gyr}^{-1}`] for all pairs of halo mass;
* ``massFunction``: The halo mass function [:math:`\mathrm{M}_\odot^{-1}\,\mathrm{Mpc}^{-3}`].
