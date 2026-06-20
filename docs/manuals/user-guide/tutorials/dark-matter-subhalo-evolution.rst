Dark Matter Subhalo Evolution
=============================

In this tutorial we'll use Galacticus to evolve populations of subhalos within a dark matter halo merger tree and output their properties. No baryonic physics (well, almost no baryonic physics) is included in this model.

If you haven't already installed Galacticus, you can find the installation instructions `here <https://github.com/galacticusorg/galacticus/wiki#how-do-i-install-and-use-galacticus>`_.

Running the calculation
-----------------------

For this tutorial we'll use the ``parameters/tutorials/darkMatterOnlySubHalos.xml`` parameter file:

.. code-block:: bash

   $ ./Galacticus.exe parameters/tutorials/darkMatterOnlySubHalos.xml

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
      2017, 2018, 2019, 2020, 2021
      - Andrew Benson

   M: Memory:         code +      nodes +       misc =      total
   M:            22.246Mib +   1.000  b +   9.000  b =  22.246Mib
   M: -> Begin multiple tasks
   .
   .
   .

and a file ``darkMatterOnlySubHalos.hdf5`` will have been created. Typically this model will take around 30s to run, but it may be (much) slower the first time you run it as various power spectra and linear growth tables are computed and stored.

Understanding the input
-----------------------

You can look at the entire parameter file for this tutorial `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/parameters/tutorials/darkMatterOnlySubHalos.xml>`_. Below we'll explore sections of the file to understand what they do. (Sections that have been explored in previous tutorials, e.g. cosmological parameters, will not be re-examined here.)

.. code-block:: xml

     <!-- Task -->
     <task                   value="multi"  >
       <task value="powerSpectra"      >
         <wavenumberMinimum value="1.0e-3"/>
         <wavenumberMaximum value="1.0e+2"/>
         <pointsPerDecade   value="10"    />
       </task>
       <task value="haloMassFunction"  >
         <haloMassMinimum value="1.0e06"/>
         <haloMassMaximum value="1.0e15"/>
         <pointsPerDecade value="10"    />
       </task>
       <task value="evolveForests"    />
     </task>
     <evolveForestsWorkShare value="cyclic"/>

This block tells Galacticus what task we want it to perform - here we're just asking it to compute the power spectrum, the halo mass function, and then to evolve merger tree forests.

.. code-block:: xml

     <!-- Structure formation options -->
     <linearGrowth          value="baryonsDarkMatter"                    />
     <criticalOverdensity   value="sphericalCollapseBrynsDrkMttrDrkEnrgy"/>
     <virialDensityContrast value="sphericalCollapseBrynsDrkMttrDrkEnrgy"/>
     <haloMassFunction      value="shethTormen"                           >
       <a             value="0.791"/> <!-- Best fit values from Benson, Ludlow, & Cole (2019). -->
       <normalization value="0.302"/>
       <p             value="0.218"/>
     </haloMassFunction>
     <excursionSetBarrier value="remapScale"     >
       <!-- Remap the barrier height by a constant factor to account for the difference in sigma(M) on large scales introduced by our
            choice of using a sharp-in-k-space filter on the power spectrum. -->
       <factor                    value="1.1965"            />
       <applyTo                   value="nonRates"          />
       <!-- Remap the barrier height according to the parameterization of Sheth, Mo, & Tormen (2001) to account for ellipsoidal
            collapse. -->
       <excursionSetBarrier value="remapShethMoTormen" >
         <a                         value="0.707"              />
         <b                         value="0.500"              />
         <c                         value="0.600"              />
         <applyTo                   value="nonRates"           />
         <!-- Use the critical overdensity as the barrier height in the excursion set problem. -->
         <excursionSetBarrier value="criticalOverdensity"/>
       </excursionSetBarrier>
     </excursionSetBarrier>
     <excursionSetFirstCrossing value="linearBarrier"/>

This block specifies various options related to large scale structure growth. The linear growth function is set to the ``baryonsDarkMatter`` model of `Benson (2020) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1268B>`_, in which the suppression of growth of large wavenumber modes due to baryon pressure is accounted for. Options for `solving the excursion set problem <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/excursion-set-problem.html>`_ are also included here, even though it is not used with the `Sheth & Tormen (2002) <http://adsabs.harvard.edu/abs/2002MNRAS.329...61S>`_ mass function which we select. It is included since it is required if the mass function is switched to the `Press & Schechter (1974) <http://adsabs.harvard.edu/abs/1974ApJ...187..425P>`_ form used for some non-CDM dark matter models.

.. code-block:: xml

     <intergalacticMediumState value="instantReionization">
       <reionizationRedshift           value="8.0e0"  />
       <reionizationTemperature        value="1.5e4"  />
       <presentDayTemperature          value="1.0e3"  />
       <intergalacticMediumState value="recFast"/>
     </intergalacticMediumState>

Here we select a simple model of the evolution of the IGM. The thermal history of the IGM is needed to solve for the growth of large wavenumber modes where baryon pressure delays the growth of structure. This IGM model assumes instant reionization at :math:`z`\ =8, heating to a temperature of 15,000K, and cooling to 1,000K by the present day. Prior to reionization the IGM state is computed using the `RecFast <https://www.astro.ubc.ca/people/scott/recfast.html>`_ code.

.. code-block:: xml

     <!-- Dark matter halo structure options -->
     <darkMatterProfileDMO         value="heated"              >
       <darkMatterProfileDMO value="NFW"      />
       <nonAnalyticSolver          value="numerical"/>
     </darkMatterProfileDMO>
     <darkMatterProfileHeating     value="tidal"              />
     <darkMatterProfileScaleRadius value="ludlow2016"          >
       <C                                  value="700.27000"    /> <!-- Best fit values from Johnson, Benson, & Grin (2020). -->
       <f                                  value="  0.07534"    />
       <timeFormationSeekDelta             value="  0.00000"    />
       <darkMatterProfileScaleRadius value="concentration" >
         <correctForConcentrationDefinition    value="true"              />
         <darkMatterProfileConcentration value="diemerKravtsov2014" >
   	<alpha   value="1.12"/>
   	<beta    value="1.69"/>
   	<eta0    value="6.82"/>
   	<eta1    value="1.42"/>
   	<kappa   value="0.69"/>
   	<phi0    value="6.58"/>
   	<phi1    value="1.37"/>
   	<scatter value="0.00"/>
         </darkMatterProfileConcentration>
       </darkMatterProfileScaleRadius>
     </darkMatterProfileScaleRadius>
     <darkMatterProfileMinimumConcentration value="  4.0"/>
     <darkMatterProfileMaximumConcentration value="100.0"/>

     <!-- Concentration model -->
     <darkMatterProfileScaleVirialTheoremUnresolvedEnergy value="0.5500"/> <!-- Best fit value from Johnson, Benson, & Grin (2020) -->
     <darkMatterProfileScaleVirialTheoremMassExponent     value="1.5552"/>
     <darkMatterProfileScaleVirialTheoremEnergyBoost      value="0.6773"/>

Here we select our model for dark matter halo profiles - we use an NFW profile but allowing for "heating" of the profile as described in `Pullen et al. (2016) <http://adsabs.harvard.edu/abs/2014ApJ...792...24P>`_. We also define parameters of the concentration model. In this model scale radii of halos are computed using the random walk approach of `Johnson, Benson & Grin (2021) <https://ui.adsabs.harvard.edu/abs/2020arXiv200615231J>`_. The parameters of that model are given here, along with the parameters of the "fallback" model used for nodes in the merger tree for which the random walk model can not be applied.

.. code-block:: xml

     <!-- Dark matter halo spin -->
     <haloSpinDistribution value="bett2007"> <!-- Values from Benson (2017) -->
       <alpha   value="1.7091800"/>
       <lambda0 value="0.0420190"/>
     </haloSpinDistribution>
     <spinVitvitskaMassExponent value="0.10475"/> <!-- Best fit value from Benson, Behrens, & Lu (2020) -->

Halo spin parameters are also modeled using a random walk approach `Benson, Behrens & Lu (2020) <http://adsabs.harvard.edu/abs/2020MNRAS.496.3371B>`_ - here we specify the parameters of that random walk model and a "fallback" spin distribution to be used for halos for which the random walk model can not be applied.

.. code-block:: xml

     <!-- Halo accretion options -->
     <accretionHalo value="zero"/>

     <!-- Hot halo gas model options -->
     <hotHaloMassDistribution value="null"          />

     <!-- Galactic structure solver options -->
     <galacticStructureSolver value="null"          />
     <darkMatterProfile       value="darkMatterOnly"/>

     <!-- Galaxy mergers -->
     <mergerRemnantSize value="null"/>

The preceeding block acts to switch off various baryonic physics which we want to exclude from this model.

.. code-block:: xml

     <!-- Satellite orbit options -->
     <virialOrbit value="spinCorrelated">
       <alpha             value="0.47263"  /> <!-- Best fit value from Benson, Behrens, & Lu (2020) -->
       <virialOrbit value="jiang2014" >
         <!-- Best fit value from Benson, Behrens, & Lu (2020) -->
         <bRatioHigh             value="+2.88333 +4.06371 +3.86726"/>
         <bRatioIntermediate     value="+1.05361 +1.56868 +2.89027"/>
         <bRatioLow              value="+0.07432 +0.54554 +1.04721"/>
         <gammaRatioHigh         value="+0.07124 +0.04737 -0.01913"/>
         <gammaRatioIntermediate value="+0.10069 +0.07821 +0.04231"/>
         <gammaRatioLow          value="+0.10866 +0.11260 +0.11698"/>
         <muRatioHigh            value="+1.10168 +1.09639 +1.09819"/>
         <muRatioIntermediate    value="+1.18205 +1.19573 +1.24581"/>
         <muRatioLow             value="+1.22053 +1.22992 +1.25528"/>
         <sigmaRatioHigh         value="+0.09244 +0.14335 +0.21079"/>
         <sigmaRatioIntermediate value="+0.07397 +0.09590 +0.10941"/>
         <sigmaRatioLow          value="+0.07458 +0.09040 +0.06981"/>
       </virialOrbit>
     </virialOrbit>
     <satelliteOrbitStoreOrbitalParameters value="true"/>

When halos first become subhalos they are assigned orbital parameters. The above defines the distribution of those parameters, using the `Jiang et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015MNRAS.448.1674J/abstract>`_ model for the distribution of radial and tangential velocities, plus the model of `Benson, Behrens & Lu (2020) <http://adsabs.harvard.edu/abs/2020MNRAS.496.3371B>`_ for the anisotropy of the distribution of orbits.

.. code-block:: xml

     <!-- Orbiting model of satellites -->
     <!-- Values taken from Yang et al. (2020) for their gamma=0 case using the Caterpillar simulations as calibration target -->
     <satelliteDynamicalFriction value="chandrasekhar1943">
       <logarithmCoulomb value="1.53"/>
     </satelliteDynamicalFriction>
     <satelliteTidalHeatingRate  value="gnedin1999"       >
       <epsilon          value="0.33"/>
       <gamma            value="0.00"/>
     </satelliteTidalHeatingRate>
     <satelliteTidalStripping    value="zentner2005"      >
       <efficiency       value="2.86"/>
     </satelliteTidalStripping>

Here we define parameters of the various physics models applied to the evolution of subhalos, as outlined in `Yang et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.498.3902Y>`_.

.. code-block:: xml

     <!-- Node evolution and physics -->
     <nodeOperator value="multi">
       <!-- Subhalo hierarchy -->
       <nodeOperator value="subsubhaloPromotion"   />
       <!-- Subhalo orbits -->
       <nodeOperator value="satelliteOrbit"                   />
       <nodeOperator value="satelliteDynamicalFriction"       />
       <nodeOperator value="satelliteTidalMassLoss"           />
       <nodeOperator value="satelliteTidalHeating"            />
       <nodeOperator value="satelliteMergingRadiusTrigger"     >
         <radiusVirialFraction      value="0.01"/>
       </nodeOperator>
       <nodeOperator value="satelliteDestructionMassThreshold" >
         <massDestruction           value="=[mergerTreeMassResolution::massResolution]"/>
         <massDestructionFractional value="0.0e0"/>
       </nodeOperator>
     </nodeOperator>

This block applies various physical processes to subhalos.

.. code-block:: xml

     <!-- Output options -->
     <galacticusOutputFileName value="darkMatterOnlySubHalos.hdf5"/>
     <mergerTreeOutputter value="multi">
       <mergerTreeOutputter value="standard">
         <outputReferences value="false"/>
       </mergerTreeOutputter>
       <mergerTreeOutputter value="analyzer"/>
     </mergerTreeOutputter>
     <outputTimes value="list">
       <redshifts value="0.5"/>
     </outputTimes>
     <nodePropertyExtractor value="multi">
       <nodePropertyExtractor value="nodeIndices"          />
       <nodePropertyExtractor value="indicesTree"          />
       <nodePropertyExtractor value="redshiftLastIsolated" />
       <nodePropertyExtractor value="radiusTidal"          />
       <nodePropertyExtractor value="radiusBoundMass"      />
       <nodePropertyExtractor value="virialProperties"     />
       <nodePropertyExtractor value="radiusVelocityMaximum"/>
       <nodePropertyExtractor value="velocityMaximum"      />
       <nodePropertyExtractor value="positionOrbital"      />
       <nodePropertyExtractor value="densityProfile"        >
         <includeRadii     value="true"                                                                                                                                                                            />
         <radiusSpecifiers value="darkMatterScaleRadius:all:all:0.3 darkMatterScaleRadius:all:all:1.0 darkMatterScaleRadius:all:all:3.0 virialRadius:all:all:0.1 virialRadius:all:all:0.3 virialRadius:all:all:1.0"/>
       </nodePropertyExtractor>
       <nodePropertyExtractor value="projectedDensity"     >
         <includeRadii     value="true"                                                                                                                                                                            />
         <radiusSpecifiers value="darkMatterScaleRadius:all:all:0.3 darkMatterScaleRadius:all:all:1.0 darkMatterScaleRadius:all:all:3.0 virialRadius:all:all:0.1 virialRadius:all:all:0.3 virialRadius:all:all:1.0"/>
       </nodePropertyExtractor>
     </nodePropertyExtractor>
     <outputAnalysis value="multi">
       <outputAnalysis value="subhaloMassFunction">
         <fileName                          value="%DATASTATICPATH%/darkMatter/subhaloDistributionsCaterpillar.hdf5"/>
         <negativeBinomialScatterFractional value="0.18"                                                            /> <!-- Boylan-Kolchin et al. (2010) -->
         <virialDensityContrast       value="bryanNorman1998"                                                 />
         <redshift                          value="0.0"                                                             />
       </outputAnalysis>
       <outputAnalysis value="subhaloRadialDistribution">
         <fileName                          value="%DATASTATICPATH%/darkMatter/subhaloDistributionsCaterpillar.hdf5"/>
         <negativeBinomialScatterFractional value="0.18"                                                            /> <!-- Boylan-Kolchin et al. (2010) -->
         <virialDensityContrast       value="bryanNorman1998"                                                 />
         <redshift                          value="0.0"                                                             />
       </outputAnalysis>
       <outputAnalysis value="subhaloVMaxVsMass">
         <fileName                          value="%DATASTATICPATH%/darkMatter/subhaloDistributionsCaterpillar.hdf5"/>
         <virialDensityContrast       value="bryanNorman1998"                                                 />
         <redshift                          value="0.0"                                                             />
       </outputAnalysis>
     </outputAnalysis>

Finally we specify the outputs that we want. Here we request a number of additional properties to be output, and also request thee summary analyses (subhalo mass function, :math:`V_\mathrm{max}` function, and radial distribution - each compared to results from the Caterpillar simulations) to be computed and output.

Understanding the output
------------------------

The main output follows the general form described in the tutorial on dark matter only merger trees, but with some additional properties (subhalo orbital positions, tidal radii, etc.)

In addition, there is now an ``analyses`` group:

.. code-block:: console

   $ h5ls galacticus.hdf5/analyses
   subhaloMassFunction      Group
   subhaloRadialDistribution Group
   subhaloVelocityMaximumMean Group

Each group contains the results of one summary analysis. For example:

.. code-block:: console

   $ h5ls darkMatterOnlySubHalos.hdf5/analyses/subhaloMassFunction
   massFunction             Dataset {20}
   massFunctionCovariance   Dataset {20, 20}
   massFunctionCovarianceTarget Dataset {20, 20}
   massFunctionTarget       Dataset {20}
   massRatio                Dataset {20}

gives the results of a subhalo mass function analysis. The mass function is given by the ``massFunction`` dataset in bins of subhalo-to-host-halo mass ratio as given by the ``massRatio`` dataset. Additionally the covariance of the mass function is given in the ``massFunctionCovariance`` dataset. For comparison the subhalo mass function and its covariance from a target dataset (in this case the `Caterpillar simulations <https://ui.adsabs.harvard.edu/abs/2016ApJ...818...10G>`_) are included also.
