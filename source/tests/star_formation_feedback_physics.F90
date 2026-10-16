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
Contains a program which tests star formation and stellar feedback physics against independently computed reference values.
!!}

!+ Contributions to this file made by: Andrew Benson, Claude.

program Test_Star_Formation_Feedback_Physics
  !!{RST
  Tests star formation and stellar feedback physics (the Kennicutt-Schmidt star formation rate surface density with and without
  truncation, the disk-integrated star formation rate, the dynamical time star formation timescale, and power-law and
  rate-limited stellar feedback outflows) against reference values computed independently by the
  ``starFormationFeedbackPhysics.py`` script in the `galacticusDevTools <https://github.com/galacticusorg/galacticusDevTools>`_
  repository. Tolerances are:

  * quantities involving no physical constants (star formation rate surface densities where no truncation is applied,
    power-law outflow rates): :math:`10^{-8}` relative, limited only by double precision arithmetic;
  * the truncated star formation rate surface density: :math:`10^{-3}` relative, since the critical surface density depends on
    :math:`\mathrm{G}`, which differs by :math:`7\times 10^{-5}` between the GSL (Galacticus) and CODATA 2018 (reference)
    values, raised to the power of the truncation exponent (6);
  * dynamical times, and quantities derived from them: :math:`10^{-4}` relative, since the conversion from Mpc/(km/s) to Gyr
    differs by :math:`2\times 10^{-5}` because of different definitions of the year;
  * disk-integrated star formation rates: :math:`10^{-5}` (untruncated) and :math:`10^{-4}` (truncated) relative, set by the
    :math:`10^{-6}` relative tolerance of the numerical integration, and the difference in :math:`\mathrm{G}` in the truncated
    case.
  !!}
  use :: Display                                  , only : displayVerbositySet                                 , &
  &                                                        verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize                        , &
  &                                                        nodeComponentBasic                                  , &
  &                                                        nodeComponentDisk                                   , &
  &                                                        nodeComponentSpheroid                               , &
  &                                                        treeNode
  use :: ISO_Varying_String                       , only : assignment(=)                                       , &
  &                                                        varying_string
  use :: Input_Parameters                         , only : inputParameters
  use :: Node_Components                          , only : Node_Components_Initialize                          , &
  &                                                        Node_Components_Thread_Initialize                   , &
  &                                                        Node_Components_Thread_Uninitialize                 , &
  &                                                        Node_Components_Uninitialize
  use :: Star_Formation_Rate_Surface_Density_Disks, only : starFormationRateSurfaceDensityDisksKennicuttSchmidt
  use :: Star_Formation_Rates_Disks               , only : starFormationRateDisksIntgrtdSurfaceDensity
  use :: Star_Formation_Timescales                , only : starFormationTimescaleDynamicalTime
  use :: Stellar_Feedback_Outflows                , only : stellarFeedbackOutflowsPowerLaw                     , &
  &                                                        stellarFeedbackOutflowsRateLimit
  use :: Unit_Tests                               , only : Assert                                              , &
  &                                                        Unit_Tests_Begin_Group                              , &
  &                                                        Unit_Tests_End_Group                                , &
  &                                                        Unit_Tests_Finish
  implicit none
  type            (treeNode                                            ), pointer      :: node
  class           (nodeComponentBasic                                  ), pointer      :: basic
  class           (nodeComponentDisk                                   ), pointer      :: disk
  class           (nodeComponentSpheroid                               ), pointer      :: spheroid
  type            (starFormationRateSurfaceDensityDisksKennicuttSchmidt)               :: starFormationRateSurfaceDensityDisksTruncated_                            , &
       &                                                                                  starFormationRateSurfaceDensityDisksUntruncated_
  type            (starFormationRateDisksIntgrtdSurfaceDensity         )               :: starFormationRateDisksTruncated_                                          , &
       &                                                                                  starFormationRateDisksUntruncated_
  type            (starFormationTimescaleDynamicalTime                 )               :: starFormationTimescale_
  type            (stellarFeedbackOutflowsPowerLaw                     )               :: stellarFeedbackOutflowsPowerLaw_
  type            (stellarFeedbackOutflowsRateLimit                    )               :: stellarFeedbackOutflowsRateLimit_
  type            (varying_string                                      )               :: parameterFile
  type            (inputParameters                                     )               :: parameters
  ! Properties of the test disk (metal-free gas, so that the hydrogen mass fraction is the primordial value).
  double precision                                                      , parameter    :: massGas                                         =1.0d10                   , &
       &                                                                                  radiusDisk                                      =3.0d-3                   , &
       &                                                                                  velocityDisk                                    =1.5d02
  ! Radii (in units of the disk scale length) at which to test star formation rate surface densities.
  double precision                                                      , dimension(4) :: radiusDimensionless                             =[0.5d0,1.0d0,3.0d0,6.0d0]
  double precision                                                      , dimension(4) :: rateSurfaceDensityTruncatedReference      =[6.8570850400d+13,3.4051276577d+13,2.0706602611d+12,1.2213488275d+05], &
       &                                                                                  rateSurfaceDensityUntruncatedReference    =[6.8570850400d+13,3.4051276577d+13,2.0706602611d+12,3.1050745015d+10]
  double precision                                                      , dimension(4) :: rateSurfaceDensityTruncated                                               , &
       &                                                                                  rateSurfaceDensityUntruncated
  double precision                                                                     :: rateOutflowEjective                                                       , &
       &                                                                                  rateOutflowExpulsive
  integer                                                                              :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Star formation and feedback physics")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/diskSpheroidComponents.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Construct the physics objects to test, using the parameter values of the Galacticus reference models where available.
  !![
  <referenceConstruct object="starFormationRateSurfaceDensityDisksTruncated_"  >
   <constructor>
    starFormationRateSurfaceDensityDisksKennicuttSchmidt(normalization=0.147d0,exponent=1.4d0,truncate=.true. ,exponentTruncated=6.0d0,velocityDispersionDiskGas=10.0d0,toomreParameterCritical=0.4d0)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="starFormationRateSurfaceDensityDisksUntruncated_">
   <constructor>
    starFormationRateSurfaceDensityDisksKennicuttSchmidt(normalization=0.147d0,exponent=1.4d0,truncate=.false.,exponentTruncated=6.0d0,velocityDispersionDiskGas=10.0d0,toomreParameterCritical=0.4d0)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="starFormationRateDisksTruncated_"                >
   <constructor>
    starFormationRateDisksIntgrtdSurfaceDensity(tolerance=1.0d-6,starFormationRateSurfaceDensityDisks_=starFormationRateSurfaceDensityDisksTruncated_  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="starFormationRateDisksUntruncated_"              >
   <constructor>
    starFormationRateDisksIntgrtdSurfaceDensity(tolerance=1.0d-6,starFormationRateSurfaceDensityDisks_=starFormationRateSurfaceDensityDisksUntruncated_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="starFormationTimescale_"                         >
   <constructor>
    starFormationTimescaleDynamicalTime(efficiency=0.04d0,exponentVelocity=2.0d0,timescaleMinimum=1.0d-3)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="stellarFeedbackOutflowsPowerLaw_"                >
   <constructor>
    stellarFeedbackOutflowsPowerLaw(velocityCharacteristic_=250.0d0,exponent=2.0d0)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="stellarFeedbackOutflowsRateLimit_"               >
   <constructor>
    stellarFeedbackOutflowsRateLimit(timescaleOutflowFractionalMinimum=1.0d-3,stellarFeedbackOutflows_=stellarFeedbackOutflowsPowerLaw_)
   </constructor>
  </referenceConstruct>
  !!]
  ! Create a node with an exponential disk.
  node     => treeNode      (                 )
  basic    => node%basic    (autoCreate=.true.)
  disk     => node%disk     (autoCreate=.true.)
  spheroid => node%spheroid (autoCreate=.true.)
  call disk%massGasSet (massGas     )
  call disk%radiusSet  (radiusDisk  )
  call disk%velocitySet(velocityDisk)

  ! Kennicutt-Schmidt star formation rate surface density.
  call Unit_Tests_Begin_Group("Kennicutt-Schmidt star formation rate surface density")
  do i=1,size(radiusDimensionless)
     rateSurfaceDensityTruncated  (i)=starFormationRateSurfaceDensityDisksTruncated_  %rate(node,radiusDimensionless(i)*radiusDisk)
     rateSurfaceDensityUntruncated(i)=starFormationRateSurfaceDensityDisksUntruncated_%rate(node,radiusDimensionless(i)*radiusDisk)
  end do
  call Assert("untruncated"                       ,rateSurfaceDensityUntruncated     ,rateSurfaceDensityUntruncatedReference     ,relTol=1.0d-8)
  call Assert("truncated {above critical density}",rateSurfaceDensityTruncated  (1:3),rateSurfaceDensityTruncatedReference  (1:3),relTol=1.0d-8)
  call Assert("truncated {below critical density}",rateSurfaceDensityTruncated  (4  ),rateSurfaceDensityTruncatedReference  (4  ),relTol=1.0d-3)
  call Unit_Tests_End_Group()

  ! Disk-integrated star formation rate.
  call Unit_Tests_Begin_Group("Disk-integrated star formation rate")
  call Assert("untruncated {closed form}"       ,starFormationRateDisksUntruncated_%rate(node),3.9838828798d+09,relTol=1.0d-5)
  call Assert("truncated"                       ,starFormationRateDisksTruncated_  %rate(node),3.7715535711d+09,relTol=1.0d-4)
  call Unit_Tests_End_Group()

  ! Star formation timescale.
  call Unit_Tests_Begin_Group("Dynamical time star formation timescale")
  call Assert("disk"                            ,starFormationTimescale_%timescale(disk    ),2.7500406235d-01,relTol=1.0d-4)
  call spheroid%radiusSet  (1.0d-3)
  call spheroid%velocitySet(3.0d+2)
  call Assert("spheroid"                        ,starFormationTimescale_%timescale(spheroid),1.8333604157d-01,relTol=1.0d-4)
  call spheroid%radiusSet  (1.0d-6)
  call Assert("spheroid {timescale floor}"      ,starFormationTimescale_%timescale(spheroid),1.0000000000d-03,relTol=1.0d-8)
  call Unit_Tests_End_Group()

  ! Stellar feedback outflows.
  call Unit_Tests_Begin_Group("Stellar feedback outflows")
  call stellarFeedbackOutflowsPowerLaw_ %outflowRate(disk,rateStarFormation=0.0d0,rateEnergyInput=1.0d14,rateOutflowEjective=rateOutflowEjective,rateOutflowExpulsive=rateOutflowExpulsive)
  call Assert("power-law"                       ,[rateOutflowEjective,rateOutflowExpulsive],[6.1496076550d+08,0.0d0],relTol=1.0d-8)
  call stellarFeedbackOutflowsRateLimit_%outflowRate(disk,rateStarFormation=0.0d0,rateEnergyInput=1.0d14,rateOutflowEjective=rateOutflowEjective,rateOutflowExpulsive=rateOutflowExpulsive)
  call Assert("rate-limited {limit not reached}",rateOutflowEjective                      ,6.1496076550d+08        ,relTol=1.0d-8)
  call stellarFeedbackOutflowsRateLimit_%outflowRate(disk,rateStarFormation=0.0d0,rateEnergyInput=1.0d20,rateOutflowEjective=rateOutflowEjective,rateOutflowExpulsive=rateOutflowExpulsive)
  call Assert("rate-limited {limit reached}"    ,rateOutflowEjective                      ,5.1135608252d+14        ,relTol=1.0d-4)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Star_Formation_Feedback_Physics
