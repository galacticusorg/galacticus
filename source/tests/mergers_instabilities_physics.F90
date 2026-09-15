!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{RST
Contains a program which tests galaxy merger and disk instability physics against independently computed reference values.
!!}

!+ Contributions to this file made by: Andrew Benson, Claude.

program Test_Mergers_Instabilities_Physics
  !!{RST
  Tests the :cite:t:`efstathiou_stability_1982` bar instability model and the :cite:t:`cole_hierarchical_2000` merger remnant
  size algorithm against reference values computed independently by the ``mergersInstabilitiesPhysics.py`` script in the
  `galacticusDevTools <https://github.com/galacticusorg/galacticusDevTools>`_ repository. Tolerances are:

  * the stability estimator for a maximally unstable disk (which equals its value for an isolated exponential disk): :math:`10^{-8}`
    relative;
  * the stability estimator otherwise: :math:`10^{-4}` relative, since it depends on :math:`\sqrt{\mathrm{G}\mathrm{M}_\odot}`,
    which differs by :math:`3.5\times 10^{-5}` between the GSL (Galacticus) and CODATA 2018/IAU 2015 (reference) values;
  * bar formation timescales: :math:`10^{-4}` relative for a maximally unstable disk (for which the timescale is the dynamical
    time, and the conversion from Mpc/(km/s) to Gyr differs by :math:`2\times 10^{-5}` because of different definitions of the
    year), and :math:`2\times 10^{-3}` otherwise, since the timescale scales as :math:`(\epsilon_\mathrm{c}-\epsilon)^{-2}` which
    amplifies the difference in the estimator by a factor :math:`2\epsilon/(\epsilon_\mathrm{c}-\epsilon)\approx 19`;
  * merger remnant half mass radii: :math:`10^{-8}` relative, as they depend only on the masses and radii of the progenitors and
    on no physical constants;
  * merger remnant circular velocities: :math:`10^{-4}` relative, since they depend on
    :math:`\mathrm{G}\mathrm{M}_\odot`, which differs by :math:`7\times 10^{-5}` between the GSL (Galacticus) and CODATA
    2018/IAU 2015 (reference) values;
  * the merger remnant specific angular momentum: an absolute tolerance, since the progenitors here are given zero angular
    momentum and so the result must be exactly zero.
  !!}
  use :: Display                                , only : displayVerbositySet                         , &
  &                                                      verbosityLevelStandard
  use :: Events_Hooks                           , only : eventsHooksInitialize
  use :: Functions_Global_Utilities             , only : Functions_Global_Set
  use :: Galactic_Dynamics_Bar_Instabilities    , only : galacticDynamicsBarInstabilityEfstathiou1982
  use :: Galacticus_Nodes                       , only : nodeClassHierarchyInitialize                , &
  &                                                      nodeComponentBasic                          , &
  &                                                      nodeComponentDisk                           , &
  &                                                      nodeComponentSpheroid                       , &
  &                                                      treeNode
  use :: ISO_Varying_String                     , only : assignment(=)                               , &
  &                                                      varying_string
  use :: Input_Parameters                       , only : inputParameters
  use :: Node_Components                        , only : Node_Components_Initialize                  , &
  &                                                      Node_Components_Thread_Initialize           , &
  &                                                      Node_Components_Thread_Uninitialize         , &
  &                                                      Node_Components_Uninitialize
  use :: Satellite_Merging_Mass_Movements       , only : destinationMergerSpheroid                   , &
  &                                                      mergerMassMovementsSimple
  use :: Satellite_Merging_Progenitor_Properties, only : mergerProgenitorPropertiesStandard
  use :: Satellite_Merging_Remnant_Sizes        , only : mergerRemnantSizeCole2000
  use :: Unit_Tests                             , only : Assert                                      , &
  &                                                      Unit_Tests_Begin_Group                      , &
  &                                                      Unit_Tests_End_Group                        , &
  &                                                      Unit_Tests_Finish
  implicit none
  type            (treeNode                                    ), pointer      :: node
  class           (nodeComponentBasic                          ), pointer      :: basic
  class           (nodeComponentDisk                           ), pointer      :: disk
  type            (galacticDynamicsBarInstabilityEfstathiou1982)               :: galacticDynamicsBarInstability_
  type            (varying_string                              )               :: parameterFile
  type            (inputParameters                             )               :: parameters
  ! Disk velocities (at the scale length) giving stable, unstable, and maximally unstable disks.
  double precision                                              , dimension(3) :: velocityDisk                           =[150.0d0,90.0d0,10.0d0]
  double precision                                              , dimension(3) :: estimatorReference                     =[1.4782983544d+0,8.8697901264d-1,6.2212973163d-1] , &
       &                                                                          timescaleReference                     =[-1.0000000000d+0,4.8240818905d-1,2.9333766650d-1]
  double precision                                              , dimension(3) :: estimator                                                                                 , &
       &                                                                          timescale
  double precision                                                             :: externalDrivingSpecificTorque                                                             , &
       &                                                                          fractionAngularMomentumRetainedDisk                                                       , &
       &                                                                          fractionAngularMomentumRetainedSpheroid
  integer                                                                      :: i
  ! Objects and nodes used in testing merger remnant sizes.
  type            (mergerMassMovementsSimple                   )               :: mergerMassMovements_
  type            (mergerProgenitorPropertiesStandard          )               :: mergerProgenitorProperties_
  type            (mergerRemnantSizeCole2000                   )               :: mergerRemnantSize_
  type            (treeNode                                    ), pointer      :: nodeHost                                                                                  , &
       &                                                                          nodeSatellite
  class           (nodeComponentSpheroid                       ), pointer      :: spheroidHost                                                                              , &
       &                                                                          spheroidSatellite
  double precision                                                             :: radiusRemnant                                                                             , &
       &                                                                          velocityRemnant                                                                           , &
       &                                                                          angularMomentumSpecificRemnant
  ! Scale radii of the progenitor spheroids, and the half mass radius of a Hernquist profile in units of its scale radius.
  double precision                                              , parameter    :: radiusScaleHost                        =1.0d-3                                            , &
       &                                                                          radiusScaleSatellite                   =5.0d-4
  double precision                                              , parameter    :: radiusHalfMassHernquist                =2.4142135624d0

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Mergers and instabilities physics")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/diskSpheroidComponents.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Construct the bar instability model using the parameter values of the Galacticus reference models.
  !![
  <referenceConstruct object="galacticDynamicsBarInstability_">
   <constructor>
    galacticDynamicsBarInstabilityEfstathiou1982(stabilityThresholdStellar=1.1d0,stabilityThresholdGaseous=0.7d0,timescaleMinimum=1.0d-3,fractionAngularMomentumRetainedDisk=1.0d0,fractionAngularMomentumRetainedSpheroid=1.0d0)
   </constructor>
  </referenceConstruct>
  !!]
  ! Create a node with a disk of gas fraction 0.3.
  node  => treeNode  (                 )
  basic => node%basic(autoCreate=.true.)
  disk  => node%disk (autoCreate=.true.)
  call disk%massGasSet        (3.0d9 )
  call disk%massStellarSet    (7.0d9 )
  call disk%radiusSet         (3.0d-3)
  call disk%angularMomentumSet(1.0d0 )

  ! Efstathiou et al. (1982) bar instability.
  call Unit_Tests_Begin_Group("Efstathiou et al. (1982) bar instability")
  do i=1,size(velocityDisk)
     call disk%velocitySet(velocityDisk(i))
     estimator(i)=galacticDynamicsBarInstability_%estimator(node)
     call galacticDynamicsBarInstability_%timescale(node,timescale(i),externalDrivingSpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
  end do
  call Assert("stability estimator"                         ,estimator(1:2),estimatorReference(1:2),relTol=1.0d-4)
  call Assert("stability estimator {maximally unstable}"    ,estimator(3  ),estimatorReference(3  ),relTol=1.0d-8)
  call Assert("bar formation timescale {stable}"            ,timescale(1  ),timescaleReference(1  ),relTol=1.0d-8)
  call Assert("bar formation timescale {unstable}"          ,timescale(2  ),timescaleReference(2  ),relTol=2.0d-3)
  call Assert("bar formation timescale {maximally unstable}",timescale(3  ),timescaleReference(3  ),relTol=1.0d-4)
  call Unit_Tests_End_Group()

  ! Cole et al. (2000) merger remnant sizes. Each galaxy is a Hernquist spheroid with no disk mass, and no dark matter profile
  ! component is active, so the dark matter terms in the algorithm vanish and the remnant half mass radius is given exactly by
  ! eqn. (4.19) of that paper. Both mergers tested here are major (mass ratios of 1.0 and 0.3 exceed the 0.25 threshold), so all
  ! mass moves to the spheroid of the remnant. The progenitor angular momenta are zero, so the specific angular momentum of the
  ! remnant must also be zero. A fresh pair of nodes is used for each configuration since the remnant size and mass movements
  ! classes memoize their results by node unique ID.
  !![
  <referenceConstruct object="mergerMassMovements_"       >
   <constructor>
    mergerMassMovementsSimple(massRatioMajorMerger=0.25d0,destinationGasMinorMerger=destinationMergerSpheroid,destinationStarsMinorMerger=destinationMergerSpheroid)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="mergerProgenitorProperties_">
   <constructor>
    mergerProgenitorPropertiesStandard(mergerMassMovements_=mergerMassMovements_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="mergerRemnantSize_"         >
   <constructor>
    mergerRemnantSizeCole2000(energyOrbital=1.0d0,ignoreUnphysicalConditions=.false.,mergerProgenitorProperties_=mergerProgenitorProperties_)
   </constructor>
  </referenceConstruct>
  !!]
  call Unit_Tests_Begin_Group("Cole et al. (2000) merger remnant sizes")
  ! Equal mass merger: the remnant half mass radius must be 4/3 of that of the progenitors.
  nodeHost                  => treeNode              (                 )
  nodeSatellite             => treeNode              (                 )
  basic                     => nodeHost     %basic   (autoCreate=.true.)
  basic                     => nodeSatellite%basic   (autoCreate=.true.)
  spheroidHost              => nodeHost     %spheroid(autoCreate=.true.)
  spheroidSatellite         => nodeSatellite%spheroid(autoCreate=.true.)
  nodeSatellite%mergeTarget => nodeHost
  call spheroidHost     %massStellarSet    (1.0d10         )
  call spheroidHost     %radiusSet         (radiusScaleHost)
  call spheroidHost     %angularMomentumSet(0.0d0          )
  call spheroidSatellite%massStellarSet    (1.0d10         )
  call spheroidSatellite%radiusSet         (radiusScaleHost)
  call spheroidSatellite%angularMomentumSet(0.0d0          )
  call Assert("spheroid half mass radius",spheroidHost%halfMassRadius(),radiusHalfMassHernquist*radiusScaleHost,relTol=1.0d-8)
  call mergerRemnantSize_%get(nodeSatellite,radiusRemnant,velocityRemnant,angularMomentumSpecificRemnant)
  call Assert("remnant radius {equal mass}"          ,radiusRemnant                 ,3.2189514165d-3,relTol=1.0d-8 )
  call Assert("remnant radius is 4/3 of progenitors" ,radiusRemnant                 ,4.0d0*radiusHalfMassHernquist*radiusScaleHost/3.0d0,relTol=1.0d-8)
  call Assert("remnant circular velocity {equal mass}",velocityRemnant              ,1.6347010111d+2,relTol=1.0d-4 )
  call Assert("remnant specific angular momentum"    ,angularMomentumSpecificRemnant,0.0d+0         ,absTol=1.0d-30)
  ! Unequal mass merger, still major.
  nodeHost                  => treeNode              (                 )
  nodeSatellite             => treeNode              (                 )
  basic                     => nodeHost     %basic   (autoCreate=.true.)
  basic                     => nodeSatellite%basic   (autoCreate=.true.)
  spheroidHost              => nodeHost     %spheroid(autoCreate=.true.)
  spheroidSatellite         => nodeSatellite%spheroid(autoCreate=.true.)
  nodeSatellite%mergeTarget => nodeHost
  call spheroidHost     %massStellarSet    (1.0d10              )
  call spheroidHost     %radiusSet         (radiusScaleHost     )
  call spheroidHost     %angularMomentumSet(0.0d0               )
  call spheroidSatellite%massStellarSet    (3.0d09              )
  call spheroidSatellite%radiusSet         (radiusScaleSatellite)
  call spheroidSatellite%angularMomentumSet(0.0d0               )
  call mergerRemnantSize_%get(nodeSatellite,radiusRemnant,velocityRemnant,angularMomentumSpecificRemnant)
  call Assert("remnant radius {unequal mass}"           ,radiusRemnant  ,2.5822917218d-3,relTol=1.0d-8)
  call Assert("remnant circular velocity {unequal mass}",velocityRemnant,1.4714637550d+2,relTol=1.0d-4)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Mergers_Instabilities_Physics
