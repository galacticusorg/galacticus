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

!!{
Contains a code to test the radiative transfer state solver.
!!}

program Test_Radiative_Transfer_State_Solver
  !!{
  Test the radiative transfer state solver.
  !!}
  use :: Atomic_Cross_Sections_Ionization_Photo      , only : atomicCrossSectionIonizationPhotoVerner
  use :: Atomic_Ionization_Potentials                , only : atomicIonizationPotentialVerner
  use :: Atomic_Radiation_Gaunt_Factors              , only : gauntFactorVanHoof2014
  use :: Atomic_Rates_Excitation_Collisional         , only : atomicExcitationRateCollisionalScholzWalters1991
  use :: Atomic_Rates_Ionization_Collisional         , only : atomicIonizationRateCollisionalVerner1996
  use :: Atomic_Rates_Recombination_Dielectronic     , only : atomicRecombinationRateDielectronicArnaud1985   , atomicRecombinationRateDielectronicZero
  use :: Atomic_Rates_Recombination_Radiative        , only : atomicRecombinationRateRadiativeFixed           , atomicRecombinationRateRadiativeVerner1996
  use :: Atomic_Rates_Recombination_Radiative_Cooling, only : atomicRecombinationRateRadiativeCoolingFixed    , atomicRecombinationRateRadiativeCoolingHummer
  use :: Display                                     , only : displayVerbositySet                             , verbosityLevelStandard
  use :: Error                            , only : errorStatusSuccess
  use :: Galacticus_Nodes                            , only : nodeClassHierarchyInitialize
  use :: Functions_Global_Utilities                  , only : Functions_Global_Set
  use :: Node_Components                             , only : Node_Components_Initialize                      , Node_Components_Uninitialize
  use :: ISO_Varying_String                          , only : var_str
  use :: Input_Parameters                            , only : inputParameters
#ifdef USEMPI
  use :: MPI_F08                                     , only : MPI_Thread_Single
#endif
#ifdef USEMPI
  use :: MPI_Utilities                               , only : mpiFinalize                                     , mpiInitialize
#endif
  use :: Mass_Distributions                          , only : massDistributionConstantDensityCloud
  use :: Numerical_Constants_Astronomical            , only : metallicitySolar
  use :: Radiative_Transfer_Matters                  , only : radiativeTransferMatterAtomic                   , radiativeTransferPropertiesMatterAtomic
  use :: Unit_Tests                                  , only : Assert                                          , Unit_Tests_Begin_Group                       , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type   (inputParameters                                 ), target      :: parameters
  type   (atomicCrossSectionIonizationPhotoVerner         ), pointer     :: atomicCrossSectionIonizationPhoto_
  type   (atomicIonizationPotentialVerner                 ), pointer     :: atomicIonizationPotential_
  type   (atomicExcitationRateCollisionalScholzWalters1991), pointer     :: atomicExcitationRateCollisional_
  type   (gauntFactorVanHoof2014                          ), pointer     :: gauntFactor_
  type   (atomicIonizationRateCollisionalVerner1996       ), pointer     :: atomicIonizationRateCollisional_
  type   (atomicRecombinationRateDielectronicZero         ), pointer     :: atomicRecombinationRateDielectronicZero_
  type   (atomicRecombinationRateDielectronicArnaud1985   ), pointer     :: atomicRecombinationRateDielectronicArnaud1985_
  type   (atomicRecombinationRateRadiativeFixed           ), pointer     :: atomicRecombinationRateRadiativeFixed_
  type   (atomicRecombinationRateRadiativeVerner1996      ), pointer     :: atomicRecombinationRateRadiativeVerner1996_
  type   (atomicRecombinationRateRadiativeCoolingFixed    ), pointer     :: atomicRecombinationRateRadiativeCoolingFixed_
  type   (atomicRecombinationRateRadiativeCoolingHummer   ), pointer     :: atomicRecombinationRateRadiativeCoolingHummer_
  type   (massDistributionConstantDensityCloud            ), pointer     :: massDistribution_
  type   (radiativeTransferMatterAtomic                   ), pointer     :: radiativeTransferMatter_
  type   (radiativeTransferPropertiesMatterAtomic         ), allocatable :: properties
  integer                                                                :: status

#ifdef USEMPI
  call mpiInitialize(MPI_Thread_Single  )
#endif
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize parameters.
  parameters=inputParameters(var_str('testSuite/parameters/test-radiativeTransfer-atomicMatterStateSolver.xml'))
  call Functions_Global_Set        (          )
  call nodeClassHierarchyInitialize(parameters)
  call Node_Components_Initialize  (parameters)
  ! Construct atomic matter.
  allocate(atomicCrossSectionIonizationPhoto_            )
  allocate(atomicIonizationPotential_                    )
  allocate(atomicExcitationRateCollisional_              )
  allocate(gauntFactor_                                  )
  allocate(atomicIonizationRateCollisional_              )
  allocate(atomicRecombinationRateDielectronicZero_      )
  allocate(atomicRecombinationRateDielectronicArnaud1985_)
  allocate(atomicRecombinationRateRadiativeFixed_        )
  allocate(atomicRecombinationRateRadiativeVerner1996_   )
  allocate(atomicRecombinationRateRadiativeCoolingFixed_ )
  allocate(atomicRecombinationRateRadiativeCoolingHummer_)
  allocate(massDistribution_                             )
  !![
  <referenceConstruct object="atomicCrossSectionIonizationPhoto_"            >
   <constructor>
    atomicCrossSectionIonizationPhotoVerner         (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicIonizationPotential_"                    >
   <constructor>
    atomicIonizationPotentialVerner                 (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicExcitationRateCollisional_"              >
   <constructor>
    atomicExcitationRateCollisionalScholzWalters1991(                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="gauntFactor_"                                  >
   <constructor>
    gauntFactorVanHoof2014                          (                                                                               &amp;
      &amp;                                          atomicIonizationPotential_       =atomicIonizationPotential_                   &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateDielectronicZero_"      >
   <constructor>
    atomicRecombinationRateDielectronicZero         (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateDielectronicArnaud1985_">
   <constructor>
    atomicRecombinationRateDielectronicArnaud1985   (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeFixed_"        >
   <constructor>
    atomicRecombinationRateRadiativeFixed           (                                                                               &amp;
      &amp;                                          rateCoefficient                  =2.0d-13                                      &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeVerner1996_"   >
   <constructor>
    atomicRecombinationRateRadiativeVerner1996      (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeCoolingFixed_" >
   <constructor>
    atomicRecombinationRateRadiativeCoolingFixed    (                                                                               &amp;
      &amp;                                          atomicRecombinationRateRadiative_=atomicRecombinationRateRadiativeFixed_     , &amp;
      &amp;                                          gamma=0.75d0                                                                   &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeCoolingHummer_">
   <constructor>
    atomicRecombinationRateRadiativeCoolingHummer   (                                                                               &amp;
      &amp;                                          atomicRecombinationRateRadiative_=atomicRecombinationRateRadiativeVerner1996_, &amp;
      &amp;                                          gamma=0.67d0                                                                   &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicIonizationRateCollisional_"             >
   <constructor>
    atomicIonizationRateCollisionalVerner1996       (                                                                               &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="massDistribution_"                            >
   <constructor>
    massDistributionConstantDensityCloud            (                                                                               &amp;
      &amp;                                          mass                             =1.0d0                                      , &amp;
      &amp;                                          radius                           =1.0d0                                        &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  !!]
  ! Tests on atomic matter - H only.
  call Unit_Tests_Begin_Group("Atomic matter - H only")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_)
  allocate(properties              )
  !![
  <referenceConstruct object="radiativeTransferMatter_">
   <constructor>
    radiativeTransferMatterAtomic(                                                                                        &amp;
      &amp;                       abundancePattern                        =var_str('solar')                             , &amp;
      &amp;                       metallicity                             =0.0d0                                        , &amp;
      &amp;                       elements                                =['H ']                                       , &amp;
      &amp;                       iterationAverageCount                   =1                                            , &amp;
      &amp;                       temperatureMinimum                      =3.0d0                                        , &amp;
      &amp;                       outputRates                             =.false.                                      , &amp;
      &amp;                       outputAbsorptionCoefficients            =.false.                                      , &amp;
      &amp;                       convergencePercentile                   =0.9d0                                        , &amp;
      &amp;                       massDistribution_                       =massDistribution_                            , &amp;
      &amp;                       atomicCrossSectionIonizationPhoto_      =atomicCrossSectionIonizationPhoto_           , &amp;
      &amp;                       atomicRecombinationRateRadiative_       =atomicRecombinationRateRadiativeFixed_       , &amp;
      &amp;                       atomicRecombinationRateRadiativeCooling_=atomicRecombinationRateRadiativeCoolingFixed_, &amp;  
      &amp;                       atomicIonizationRateCollisional_        =atomicIonizationRateCollisional_             , &amp;
      &amp;                       atomicRecombinationRateDielectronic_    =atomicRecombinationRateDielectronicZero_     , &amp;
      &amp;                       atomicIonizationPotential_              =atomicIonizationPotential_                   , &amp;
      &amp;                       atomicExcitationRateCollisional_        =atomicExcitationRateCollisional_             , &amp;
      &amp;                       gauntFactor_                            =gauntFactor_                                   &amp;
      &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  allocate(properties%elements(1)                                      )
  allocate(properties%elements(1)%photoIonizationRateHistory    (1,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory       (1,0:0))
  allocate(properties%elements(1)%ionizationStateFractionHistory(1,0:1))
  allocate(properties%elements(1)%photoIonizationRate           (  0:0))
  allocate(properties%elements(1)%photoHeatingRate              (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious   (  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious      (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction       (  0:1))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties            %iterationCount             = 1
  properties            %volume                     = 0.172800000000d-20
  properties            %temperature                = 0.103698470589d+05
  properties%elements(1)%densityNumber              = 0.101064330490d+03
  properties%elements(1)%ionizationStateFraction    =[0.706879980476d+00,0.293120019524d+00]
  properties%elements(1)%photoIonizationRate        = 0.209229241840d-06
  properties%elements(1)%photoHeatingRate           = 0.331582285728d-24
  properties%elements(1)%photoIonizationRatePrevious= 0.228307228801d-09
  properties%elements(1)%photoHeatingRatePrevious   = 0.164662972542d-26
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  !! Test 2.
  properties            %iterationCount             = 1
  properties            %volume                     = 0.172800000000d-20
  properties            %temperature                = 0.100000000000d+05
  properties%elements(1)%densityNumber              = 0.101064330490d+03
  properties%elements(1)%ionizationStateFraction    =[0.100000000000d+01,0.000000000000d+00]
  properties%elements(1)%photoIonizationRate        = 0.245726973297d-04
  properties%elements(1)%photoHeatingRate           = 0.330670408342d-22
  properties%elements(1)%photoIonizationRatePrevious=-huge(0.0d0)
  properties%elements(1)%photoHeatingRatePrevious   =-huge(0.0d0)
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #2',status,errorStatusSuccess)
  !! Test 3.
  properties            %iterationCount             = 1
  properties            %volume                     = 0.17280000000000000d-20
  properties            %temperature                = 0.30000000000000000d+01
  properties%elements(1)%densityNumber              = 0.93764865414000000d-02
  properties%elements(1)%ionizationStateFraction    =[0.10000000000000000d+01,0.000000000000d+00]
  properties%elements(1)%photoIonizationRate        = 0.87630070734000000d-22
  properties%elements(1)%photoHeatingRate           = 0.15147475879338856d-38
  properties%elements(1)%photoIonizationRatePrevious= 0.00000000000000000d+00
  properties%elements(1)%photoHeatingRatePrevious   = 0.00000000000000000d+00
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #3',status,errorStatusSuccess)
  !! Test 4.
  properties            %iterationCount             = 1
  properties            %volume                     = 0.1727999999999987d-20
  properties            %temperature                = 0.8464271217162063d+05
  properties%elements(1)%densityNumber              = 0.1010643304897267d+03
  properties%elements(1)%ionizationStateFraction    =[0.4697602133585263d-09,0.9999999995302397d+00]
  properties%elements(1)%photoIonizationRate        = 0.2219143893216550d-13
  properties%elements(1)%photoHeatingRate           = 0.1516777172160954d-31
  properties%elements(1)%photoIonizationRatePrevious= 0.3471374150826434d-12
  properties%elements(1)%photoHeatingRatePrevious   = 0.2365279824300775d-30
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #4',status,errorStatusSuccess)
  ! Done with atomic matter - H only.
  call Unit_Tests_End_Group()
  deallocate(properties)
  !![
  <objectDestructor name="radiativeTransferMatter_"/>
  !!]
  ! Tests on atomic matter - H, He, O.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, O (fixed recombination coefficient)")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_)
  allocate(properties              )
  !![
  <referenceConstruct object="radiativeTransferMatter_">
   <constructor>
    radiativeTransferMatterAtomic(                                                                                        &amp;
      &amp;                       abundancePattern                        =var_str('solar')                             , &amp;
      &amp;                       metallicity                             =1.0d0                                        , &amp;
      &amp;                       elements                                =['H ','He','O ']                             , &amp;
      &amp;                       iterationAverageCount                   =2                                            , &amp;
      &amp;                       temperatureMinimum                      =3.0d0                                        , &amp;
      &amp;                       outputRates                             =.false.                                      , &amp;
      &amp;                       outputAbsorptionCoefficients            =.false.                                      , &amp;
      &amp;                       convergencePercentile                   =0.9d0                                        , &amp;
      &amp;                       massDistribution_                       =massDistribution_                            , &amp;
      &amp;                       atomicCrossSectionIonizationPhoto_      =atomicCrossSectionIonizationPhoto_           , &amp;
      &amp;                       atomicRecombinationRateRadiative_       =atomicRecombinationRateRadiativeFixed_       , &amp;
      &amp;                       atomicRecombinationRateRadiativeCooling_=atomicRecombinationRateRadiativeCoolingFixed_, &amp;  
      &amp;                       atomicIonizationRateCollisional_        =atomicIonizationRateCollisional_             , &amp;
      &amp;                       atomicRecombinationRateDielectronic_    =atomicRecombinationRateDielectronicZero_     , &amp;
      &amp;                       atomicIonizationPotential_              =atomicIonizationPotential_                   , &amp;
      &amp;                       atomicExcitationRateCollisional_        =atomicExcitationRateCollisional_             , &amp;
      &amp;                       gauntFactor_                            =gauntFactor_                                   &amp;
      &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  allocate(properties%elements(3)                                      )
  allocate(properties%elements(1)%photoIonizationRateHistory    (2,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory       (2,0:0))
  allocate(properties%elements(1)%ionizationStateFractionHistory(2,0:1))
  allocate(properties%elements(1)%photoIonizationRate           (  0:0))
  allocate(properties%elements(1)%photoHeatingRate              (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious   (  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious      (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction       (  0:1))
  allocate(properties%elements(2)%photoIonizationRateHistory    (2,0:1))
  allocate(properties%elements(2)%photoHeatingRateHistory       (2,0:1))
  allocate(properties%elements(2)%ionizationStateFractionHistory(2,0:2))
  allocate(properties%elements(2)%photoIonizationRate           (  0:1))
  allocate(properties%elements(2)%photoHeatingRate              (  0:1))
  allocate(properties%elements(2)%photoIonizationRatePrevious   (  0:1))
  allocate(properties%elements(2)%photoHeatingRatePrevious      (  0:1))
  allocate(properties%elements(2)%ionizationStateFraction       (  0:2))
  allocate(properties%elements(3)%photoIonizationRateHistory    (2,0:7))
  allocate(properties%elements(3)%photoHeatingRateHistory       (2,0:7))
  allocate(properties%elements(3)%ionizationStateFractionHistory(2,0:8))
  allocate(properties%elements(3)%photoIonizationRate           (  0:7))
  allocate(properties%elements(3)%photoHeatingRate              (  0:7))
  allocate(properties%elements(3)%photoIonizationRatePrevious   (  0:7))
  allocate(properties%elements(3)%photoHeatingRatePrevious      (  0:7))
  allocate(properties%elements(3)%ionizationStateFraction       (  0:8))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties            %iterationCount             = 1
  properties            %volume                     = 0.2160000000000002d-18
  properties            %temperature                = 0.1481414020036419d+05
  properties%elements(1)%densityNumber              = 0.1011365525512600d+03
  properties%elements(1)%ionizationStateFraction    =[0.1260106732174199d-01,0.9873989326782581d+00]
  properties%elements(1)%photoIonizationRate        =[0.1721111089433651d-06                       ]
  properties%elements(1)%photoHeatingRate           =[0.1879261275044316d-24                       ]
  properties%elements(1)%photoIonizationRatePrevious=[0.1721003789285937d-06                       ]
  properties%elements(1)%photoHeatingRatePrevious   =[0.1910456649414847d-24                       ]
  properties%elements(2)%densityNumber              = 0.9870346689011134d+01
  properties%elements(2)%ionizationStateFraction    =[0.9019180567862987d-02,0.9909808194321362d+00,0.7780903148005057d-15]
  properties%elements(2)%photoIonizationRate        =[0.2004719134478281d-07,0.3244112744730688d-09                       ]
  properties%elements(2)%photoHeatingRate           =[0.2339803216069218d-25,0.0000000000000000d+00                       ]
  properties%elements(2)%photoIonizationRatePrevious=[0.2379621870382762d-07,0.3178194427298822d-27                       ]
  properties%elements(2)%photoHeatingRatePrevious   =[0.3044954784678252d-25,0.0000000000000000d+00                       ]
  properties%elements(3)%densityNumber              = 0.5082336812053197d-01
  properties%elements(3)%ionizationStateFraction    =[0.2147362946297900d-02,0.9978526353275143d+00,0.1726187852461528d-08,0.3548529101814460d-24,0.3127108511424789d-48,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRate        =[0.3010419326041556d-09,0.2236398465502293d-08,0.9331322409705859d-20,0.4049712122482999d-37,0.7189693608635254d-64,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00                       ]
  properties%elements(3)%photoHeatingRate           =[0.6511272500195937d-27,0.2120363621345119d-26,0.9951339529553514d-38,0.4599199679905641d-55,0.4207242937159813d-82,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00                       ]
  properties%elements(3)%photoIonizationRatePrevious=[0.5176521606053964d-09,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00                       ]
  properties%elements(3)%photoHeatingRatePrevious   =[0.1135594953109664d-26,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00,0.0000000000000000d+00                       ]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  !! Test 2.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.2160000000000002d-18
  properties            %temperature                =  0.4198998414060693d+05
  properties%elements(1)%densityNumber              =  0.1011365525512600d+03
  properties%elements(1)%ionizationStateFraction    =[ 0.5403575500822897d-03, 0.9994596424499178d+00]
  properties%elements(1)%photoIonizationRate        =[ 0.2986551157222103d-07]
  properties%elements(1)%photoHeatingRate           =[ 0.2506301481143166d-25]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.2782220518204434d-07]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.2365973160630620d-25]
  properties%elements(2)%densityNumber              =  0.9870346689011134d+01
  properties%elements(2)%ionizationStateFraction    =[ 0.3139403621081789d-02, 0.9968592824112124d+00, 0.1313967705898324d-05]
  properties%elements(2)%photoIonizationRate        =[ 0.1925530934796398d-08, 0.2397328271031403d-10]
  properties%elements(2)%photoHeatingRate           =[ 0.1960262625478408d-26, 0.2561382715174376d-28]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.1333270819311719d-08, 0.2916001332806510d-15]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.1529378689726634d-26, 0.2265943898326056d-32]
  properties%elements(3)%densityNumber              =  0.5082336812053197d-01
  properties%elements(3)%ionizationStateFraction    =[ 0.1044467850935537d-06, 0.9142007578952905d-02, 0.9908519296681630d+00, 0.5958274853955481d-05, 0.3124499730897977d-10, 0.9208877712985887d-17, 0.4932243755023585d-45, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.2515279620229301d-11, 0.7630000867196537d-11, 0.4010527420431380d-12, 0.6265690908072481d-19, 0.7171298199418202d-27, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRate           =[ 0.6495141615887409d-29, 0.8310307960798297d-29, 0.4649877019176625d-30, 0.7306792163932982d-37, 0.4172737138621947d-45, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.2514595625513513d-11, 0.6071714478015695d-11, 0.8105863123278223d-17, 0.5841335099465117d-22, 0.2189249460811778d-28, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.6494175980698739d-29, 0.6883690263450297d-29, 0.6337570041419400d-34, 0.2453000698019630d-39, 0.1573199018660068d-46, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #2',status,errorStatusSuccess)
  !! Test 3.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.2160000000000002d-18
  properties            %temperature                =  0.1000000000000000d+05
  properties%elements(1)%densityNumber              =  0.3838830883589838d+00
  properties%elements(1)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00]
  properties%elements(1)%photoIonizationRate        =[ 0.0000000000000000d+00]
  properties%elements(1)%photoHeatingRate           =[ 0.0000000000000000d+00]
  properties%elements(1)%photoIonizationRatePrevious=[-0.1000000000000000d+300]
  properties%elements(1)%photoHeatingRatePrevious   =[-0.1000000000000000d+300]
  properties%elements(2)%densityNumber              =  0.3746478473478734d-01
  properties%elements(2)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(2)%photoIonizationRate        =[ 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(2)%photoHeatingRate           =[ 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(2)%photoIonizationRatePrevious=[-0.1000000000000000d+300,-0.1000000000000000d+300]
  properties%elements(2)%photoHeatingRatePrevious   =[-0.1000000000000000d+300,-0.1000000000000000d+300]
  properties%elements(3)%densityNumber              =  0.1929097939642225d-03
  properties%elements(3)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRate           =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRatePrevious=[-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300]
  properties%elements(3)%photoHeatingRatePrevious   =[-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300,-0.1000000000000000d+300]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #3',status,errorStatusSuccess)
  ! Done with atomic matter - H, He, O.
  call Unit_Tests_End_Group()
  deallocate(properties)
  !![
  <objectDestructor name="radiativeTransferMatter_"                     />
  <objectDestructor name="atomicRecombinationRateRadiativeCoolingFixed_"/>
  !!]
  ! Tests on atomic matter - H, He, O.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, O")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_                     )
  allocate(atomicRecombinationRateRadiativeCoolingFixed_)
  allocate(properties                                   )
  !![
  <referenceConstruct object="atomicRecombinationRateRadiativeCoolingFixed_">
   <constructor>
    atomicRecombinationRateRadiativeCoolingFixed(                                                                               &amp;
      &amp;                                      atomicRecombinationRateRadiative_=atomicRecombinationRateRadiativeVerner1996_, &amp;
      &amp;                                      gamma=0.75d0                                                                   &amp;
      &amp;                                     )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="radiativeTransferMatter_">
   <constructor>
    radiativeTransferMatterAtomic(                                                                                        &amp;
      &amp;                       abundancePattern                        =var_str('solar')                             , &amp;
      &amp;                       metallicity                             =1.0d0                                        , &amp;
      &amp;                       elements                                =['H ','He','O ']                             , &amp;
      &amp;                       iterationAverageCount                   =2                                            , &amp;
      &amp;                       temperatureMinimum                      =3.0d0                                        , &amp;
      &amp;                       outputRates                             =.false.                                      , &amp;
      &amp;                       outputAbsorptionCoefficients            =.false.                                      , &amp;
      &amp;                       convergencePercentile                   =0.9d0                                        , &amp;
      &amp;                       massDistribution_                       =massDistribution_                            , &amp;
      &amp;                       atomicCrossSectionIonizationPhoto_      =atomicCrossSectionIonizationPhoto_           , &amp;
      &amp;                       atomicRecombinationRateRadiative_       =atomicRecombinationRateRadiativeVerner1996_  , &amp;
      &amp;                       atomicRecombinationRateRadiativeCooling_=atomicRecombinationRateRadiativeCoolingFixed_, &amp;  
      &amp;                       atomicIonizationRateCollisional_        =atomicIonizationRateCollisional_             , &amp;
      &amp;                       atomicRecombinationRateDielectronic_    =atomicRecombinationRateDielectronicZero_     , &amp;
      &amp;                       atomicIonizationPotential_              =atomicIonizationPotential_                   , &amp;
      &amp;                       atomicExcitationRateCollisional_        =atomicExcitationRateCollisional_             , &amp;
      &amp;                       gauntFactor_                            =gauntFactor_                                   &amp;
      &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  allocate(properties%elements(3)                                      )
  allocate(properties%elements(1)%photoIonizationRateHistory    (2,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory       (2,0:0))
  allocate(properties%elements(1)%ionizationStateFractionHistory(2,0:1))
  allocate(properties%elements(1)%photoIonizationRate           (  0:0))
  allocate(properties%elements(1)%photoHeatingRate              (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious   (  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious      (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction       (  0:1))
  allocate(properties%elements(2)%photoIonizationRateHistory    (2,0:1))
  allocate(properties%elements(2)%photoHeatingRateHistory       (2,0:1))
  allocate(properties%elements(2)%ionizationStateFractionHistory(2,0:2))
  allocate(properties%elements(2)%photoIonizationRate           (  0:1))
  allocate(properties%elements(2)%photoHeatingRate              (  0:1))
  allocate(properties%elements(2)%photoIonizationRatePrevious   (  0:1))
  allocate(properties%elements(2)%photoHeatingRatePrevious      (  0:1))
  allocate(properties%elements(2)%ionizationStateFraction       (  0:2))
  allocate(properties%elements(3)%photoIonizationRateHistory    (2,0:7))
  allocate(properties%elements(3)%photoHeatingRateHistory       (2,0:7))
  allocate(properties%elements(3)%ionizationStateFractionHistory(2,0:8))
  allocate(properties%elements(3)%photoIonizationRate           (  0:7))
  allocate(properties%elements(3)%photoHeatingRate              (  0:7))
  allocate(properties%elements(3)%photoIonizationRatePrevious   (  0:7))
  allocate(properties%elements(3)%photoHeatingRatePrevious      (  0:7))
  allocate(properties%elements(3)%ionizationStateFraction       (  0:8))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.2160000000000002d-18
  properties            %temperature                =  0.3000000000000000d+01
  properties%elements(1)%densityNumber              =  0.1011365525512600d+03
  properties%elements(1)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00]
  properties%elements(1)%photoIonizationRate        =[ 0.3427463768689085d-12]
  properties%elements(1)%photoHeatingRate           =[ 0.1475726281552722d-29]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.0000000000000000d+00]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.0000000000000000d+00]
  properties%elements(2)%densityNumber              =  0.9870346689011134d+01
  properties%elements(2)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(2)%photoIonizationRate        =[ 0.3465361558728980d-12, 0.0000000000000000d+00]
  properties%elements(2)%photoHeatingRate           =[ 0.9739883853483459d-30, 0.0000000000000000d+00]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%densityNumber              =  0.5082336812053197d-01
  properties%elements(3)%ionizationStateFraction    =[ 0.1000000000000000d+01, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.4748522330473645d-14, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRate           =[ 0.2258782769889694d-31, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  ! Done with atomic matter - H, He, O.
  call Unit_Tests_End_Group()
  deallocate(properties)
  !![
  <objectDestructor name="radiativeTransferMatter_"                     />
  <objectDestructor name="atomicRecombinationRateRadiativeCoolingFixed_"/>
  !!]
  ! Tests on atomic matter - H, He, C, O.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, C, O")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_                     )
  allocate(properties                                   )
  !![
  <referenceConstruct object="radiativeTransferMatter_">
   <constructor>
    radiativeTransferMatterAtomic(                                                                                         &amp;
      &amp;                       abundancePattern                        =var_str('solar')                              , &amp;
      &amp;                       metallicity                             =1.0d0                                         , &amp;
      &amp;                       elements                                =['H ','He','C ','O ']                         , &amp;
      &amp;                       iterationAverageCount                   =2                                             , &amp;
      &amp;                       temperatureMinimum                      =3.0d0                                         , &amp;
      &amp;                       outputRates                             =.false.                                       , &amp;
      &amp;                       outputAbsorptionCoefficients            =.false.                                       , &amp;
      &amp;                       convergencePercentile                   =0.9d0                                         , &amp;
      &amp;                       massDistribution_                       =massDistribution_                             , &amp;
      &amp;                       atomicCrossSectionIonizationPhoto_      =atomicCrossSectionIonizationPhoto_            , &amp;
      &amp;                       atomicRecombinationRateRadiative_       =atomicRecombinationRateRadiativeVerner1996_   , &amp;
      &amp;                       atomicRecombinationRateRadiativeCooling_=atomicRecombinationRateRadiativeCoolingHummer_, &amp;  
      &amp;                       atomicIonizationRateCollisional_        =atomicIonizationRateCollisional_              , &amp;
      &amp;                       atomicRecombinationRateDielectronic_    =atomicRecombinationRateDielectronicArnaud1985_, &amp;
      &amp;                       atomicIonizationPotential_              =atomicIonizationPotential_                    , &amp;
      &amp;                       atomicExcitationRateCollisional_        =atomicExcitationRateCollisional_              , &amp;
      &amp;                       gauntFactor_                            =gauntFactor_                                    &amp;
      &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  allocate(properties%elements(4)                                    )
  allocate(properties%elements(1)%photoIonizationRateHistory    (2,0: 0))
  allocate(properties%elements(1)%photoHeatingRateHistory       (2,0: 0))
  allocate(properties%elements(1)%ionizationStateFractionHistory(2,0: 1))
  allocate(properties%elements(1)%photoIonizationRate           (  0: 0))
  allocate(properties%elements(1)%photoHeatingRate              (  0: 0))
  allocate(properties%elements(1)%photoIonizationRatePrevious   (  0: 0))
  allocate(properties%elements(1)%photoHeatingRatePrevious      (  0: 0))
  allocate(properties%elements(1)%ionizationStateFraction       (  0: 1))
  allocate(properties%elements(2)%photoIonizationRateHistory    (2,0: 1))
  allocate(properties%elements(2)%photoHeatingRateHistory       (2,0: 1))
  allocate(properties%elements(2)%ionizationStateFractionHistory(2,0: 2))
  allocate(properties%elements(2)%photoIonizationRate           (  0: 1))
  allocate(properties%elements(2)%photoHeatingRate              (  0: 1))
  allocate(properties%elements(2)%photoIonizationRatePrevious   (  0: 1))
  allocate(properties%elements(2)%photoHeatingRatePrevious      (  0: 1))
  allocate(properties%elements(2)%ionizationStateFraction       (  0: 2))
  allocate(properties%elements(3)%photoIonizationRateHistory    (2,0: 5))
  allocate(properties%elements(3)%photoHeatingRateHistory       (2,0: 5))
  allocate(properties%elements(3)%ionizationStateFractionHistory(2,0: 6))
  allocate(properties%elements(3)%photoIonizationRate           (  0: 5))
  allocate(properties%elements(3)%photoHeatingRate              (  0: 5))
  allocate(properties%elements(3)%photoIonizationRatePrevious   (  0: 5))
  allocate(properties%elements(3)%photoHeatingRatePrevious      (  0: 5))
  allocate(properties%elements(3)%ionizationStateFraction       (  0: 6))
  allocate(properties%elements(4)%photoIonizationRateHistory    (2,0: 7))
  allocate(properties%elements(4)%photoHeatingRateHistory       (2,0: 7))
  allocate(properties%elements(4)%ionizationStateFractionHistory(2,0: 8))
  allocate(properties%elements(4)%photoIonizationRate           (  0: 7))
  allocate(properties%elements(4)%photoHeatingRate              (  0: 7))
  allocate(properties%elements(4)%photoIonizationRatePrevious   (  0: 7))
  allocate(properties%elements(4)%photoHeatingRatePrevious      (  0: 7))
  allocate(properties%elements(4)%ionizationStateFraction       (  0: 8))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.1099557428756440d-09
  properties            %temperature                =  0.5988888361234766d+04
  properties%elements(1)%densityNumber              =  0.2699662043065019d-06
  properties%elements(1)%ionizationStateFraction    =[ 0.9999970041055877d+00, 0.2995894412308993d-05]
  properties%elements(1)%photoIonizationRate        =[ 0.8528711304312218d-33]
  properties%elements(1)%photoHeatingRate           =[ 0.2495627959916713d-49]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.2246656067830424d-34]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.7682933462658509d-51]
  properties%elements(2)%densityNumber              =  0.2634715109031453d-07
  properties%elements(2)%ionizationStateFraction    =[ 0.9999498767086821d+00, 0.5012268169036738d-04, 0.6096276172311671d-09]
  properties%elements(2)%photoIonizationRate        =[ 0.2075258708097711d-32, 0.1795915552913853d-36]
  properties%elements(2)%photoHeatingRate           =[ 0.5708478755126612d-49, 0.4531390353136864d-53]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.5591329836617498d-34, 0.1743588752315702d-36]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.1814322764111481d-50, 0.4889431865125934d-53]
  properties%elements(3)%densityNumber              =  0.6783201243990366d-10
  properties%elements(3)%ionizationStateFraction    =[ 0.2492551483642307d-10, 0.9999780447546435d+00, 0.2195493489120175d-04, 0.2855386069381062d-09, 0.1219671347129668d-14, 0.2880010539242122d-28, 0.1853110394089593d-48]
  properties%elements(3)%photoIonizationRate        =[ 0.9689329732428958d-32, 0.1178163694270943d-34, 0.5991324991077991d-39, 0.1017355051899161d-42, 0.4723970299935798d-54, 0.9113786933031165d-74]
  properties%elements(3)%photoHeatingRate           =[ 0.1808960267003195d-50, 0.3248272677009963d-51, 0.1572053282715295d-55, 0.2579133827496697d-59, 0.5475067136514638d-72, 0.1039732323283453d-91]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.1907746573409413d-31, 0.3560491805000814d-36, 0.5743112054861547d-39, 0.1527952246756383d-42, 0.6927894530173352d-54, 0.9602369022756033d-74]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.3569503397674029d-50, 0.1159386339825431d-52, 0.1674296689676867d-55, 0.4026527998389647d-59, 0.8240384944854353d-72, 0.1106733245937702d-91]
  properties%elements(4)%densityNumber              =  0.1356640248798074d-09
  properties%elements(4)%ionizationStateFraction    =[ 0.9995265025837968d+00, 0.4734501979155830d-03, 0.4721674899750196d-07, 0.1538660975361139d-11, 0.1998492641871421d-16, 0.1190734247247572d-21, 0.2974637051184671d-27, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.7122483784241892d-34, 0.7944419410095338d-37, 0.1826625607269363d-39, 0.5114610738706758d-42, 0.1708773616739430d-44, 0.4579415050848281d-47, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRate           =[ 0.2085789077634737d-50, 0.2248797572064010d-53, 0.4911060706766282d-56, 0.1195228710981356d-58, 0.2996390694964834d-61, 0.6252060571641930d-64, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.2020021886284078d-35, 0.7660649149954055d-37, 0.2755282907326145d-39, 0.6049800776927208d-42, 0.1753600546876111d-44, 0.4586454271507444d-47, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.6917156832344250d-52, 0.2388874982424921d-53, 0.7683942507663893d-56, 0.1438580729006541d-58, 0.3086857184992557d-61, 0.6263533039293486d-64, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  !! Test 2.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.1099557428756440d-09
  properties            %temperature                =  0.5988888361234766d+04
  properties%elements(1)%densityNumber              =  0.2699662043065019d-06
  properties%elements(1)%ionizationStateFraction    =[ 0.9999970041055877d+00, 0.2995894412308993d-05]
  properties%elements(1)%photoIonizationRate        =[ 0.8528711304312218d-33]
  properties%elements(1)%photoHeatingRate           =[ 0.2495627959916713d-49]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.2246656067830424d-34]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.7682933462658509d-51]
  properties%elements(2)%densityNumber              =  0.2634715109031453d-07
  properties%elements(2)%ionizationStateFraction    =[ 0.9999498767086821d+00, 0.5012268169036738d-04, 0.6096276172311671d-09]
  properties%elements(2)%photoIonizationRate        =[ 0.2075258708097711d-32, 0.1795915552913853d-36]
  properties%elements(2)%photoHeatingRate           =[ 0.5708478755126612d-49, 0.4531390353136864d-53]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.5591329836617498d-34, 0.1743588752315702d-36]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.1814322764111481d-50, 0.4889431865125934d-53]
  properties%elements(3)%densityNumber              =  0.6783201243990366d-10
  properties%elements(3)%ionizationStateFraction    =[ 0.2492551483642307d-10, 0.9999780447546435d+00, 0.2195493489120175d-04, 0.2855386069381062d-09, 0.1219671347129668d-14, 0.2880010539242122d-28, 0.1853110394089593d-48]
  properties%elements(3)%photoIonizationRate        =[ 0.9689329732428958d-32, 0.1178163694270943d-34, 0.5991324991077991d-39, 0.1017355051899161d-42, 0.4723970299935798d-54, 0.9113786933031165d-74]
  properties%elements(3)%photoHeatingRate           =[ 0.1808960267003195d-50, 0.3248272677009963d-51, 0.1572053282715295d-55, 0.2579133827496697d-59, 0.5475067136514638d-72, 0.1039732323283453d-91]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.1907746573409413d-31, 0.3560491805000814d-36, 0.5743112054861547d-39, 0.1527952246756383d-42, 0.6927894530173352d-54, 0.9602369022756033d-74]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.3569503397674029d-50, 0.1159386339825431d-52, 0.1674296689676867d-55, 0.4026527998389647d-59, 0.8240384944854353d-72, 0.1106733245937702d-91]
  properties%elements(4)%densityNumber              =  0.1356640248798074d-09
  properties%elements(4)%ionizationStateFraction    =[ 0.9995265025837968d+00, 0.4734501979155830d-03, 0.4721674899750196d-07, 0.1538660975361139d-11, 0.1998492641871421d-16, 0.1190734247247572d-21, 0.2974637051184671d-27, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.7122483784241892d-34, 0.7944419410095338d-37, 0.1826625607269363d-39, 0.5114610738706758d-42, 0.1708773616739430d-44, 0.4579415050848281d-47, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRate           =[ 0.2085789077634737d-50, 0.2248797572064010d-53, 0.4911060706766282d-56, 0.1195228710981356d-58, 0.2996390694964834d-61, 0.6252060571641930d-64, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.2020021886284078d-35, 0.7660649149954055d-37, 0.2755282907326145d-39, 0.6049800776927208d-42, 0.1753600546876111d-44, 0.4586454271507444d-47, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.6917156832344250d-52, 0.2388874982424921d-53, 0.7683942507663893d-56, 0.1438580729006541d-58, 0.3086857184992557d-61, 0.6263533039293486d-64, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #2',status,errorStatusSuccess)
  ! Done with atomic matter - H, He, C, O.
  call Unit_Tests_End_Group()
  deallocate(properties)
  !![
  <objectDestructor name="radiativeTransferMatter_"/>
  !!]
  ! Tests on atomic matter - H, He, C, O, Si, Fe.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, C, O, Si, Fe")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_                     )
  allocate(properties                                   )
  !![
  <referenceConstruct object="radiativeTransferMatter_">
   <constructor>
    radiativeTransferMatterAtomic(                                                                                         &amp;
      &amp;                       abundancePattern                        =var_str('solar')                              , &amp;
      &amp;                       metallicity                             =1.0d0                                         , &amp;
      &amp;                       elements                                =['H ','He','C ','O ','Si','Fe']               , &amp;
      &amp;                       iterationAverageCount                   =2                                             , &amp;
      &amp;                       temperatureMinimum                      =3.0d0                                         , &amp;
      &amp;                       outputRates                             =.false.                                       , &amp;
      &amp;                       outputAbsorptionCoefficients            =.false.                                       , &amp;
      &amp;                       convergencePercentile                   =0.9d0                                         , &amp;
      &amp;                       massDistribution_                       =massDistribution_                             , &amp;
      &amp;                       atomicCrossSectionIonizationPhoto_      =atomicCrossSectionIonizationPhoto_            , &amp;
      &amp;                       atomicRecombinationRateRadiative_       =atomicRecombinationRateRadiativeVerner1996_   , &amp;
      &amp;                       atomicRecombinationRateRadiativeCooling_=atomicRecombinationRateRadiativeCoolingHummer_, &amp;  
      &amp;                       atomicIonizationRateCollisional_        =atomicIonizationRateCollisional_              , &amp;
      &amp;                       atomicRecombinationRateDielectronic_    =atomicRecombinationRateDielectronicArnaud1985_, &amp;
      &amp;                       atomicIonizationPotential_              =atomicIonizationPotential_                    , &amp;
      &amp;                       atomicExcitationRateCollisional_        =atomicExcitationRateCollisional_              , &amp;
      &amp;                       gauntFactor_                            =gauntFactor_                                    &amp;
      &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  allocate(properties%elements(6)                                       )
  allocate(properties%elements(1)%photoIonizationRateHistory    (2,0: 0))
  allocate(properties%elements(1)%photoHeatingRateHistory       (2,0: 0))
  allocate(properties%elements(1)%ionizationStateFractionHistory(2,0: 1))
  allocate(properties%elements(1)%photoIonizationRate           (  0: 0))
  allocate(properties%elements(1)%photoHeatingRate              (  0: 0))
  allocate(properties%elements(1)%photoIonizationRatePrevious   (  0: 0))
  allocate(properties%elements(1)%photoHeatingRatePrevious      (  0: 0))
  allocate(properties%elements(1)%ionizationStateFraction       (  0: 1))
  allocate(properties%elements(2)%photoIonizationRateHistory    (2,0: 1))
  allocate(properties%elements(2)%photoHeatingRateHistory       (2,0: 1))
  allocate(properties%elements(2)%ionizationStateFractionHistory(2,0: 2))
  allocate(properties%elements(2)%photoIonizationRate           (  0: 1))
  allocate(properties%elements(2)%photoHeatingRate              (  0: 1))
  allocate(properties%elements(2)%photoIonizationRatePrevious   (  0: 1))
  allocate(properties%elements(2)%photoHeatingRatePrevious      (  0: 1))
  allocate(properties%elements(2)%ionizationStateFraction       (  0: 2))
  allocate(properties%elements(3)%photoIonizationRateHistory    (2,0: 5))
  allocate(properties%elements(3)%photoHeatingRateHistory       (2,0: 5))
  allocate(properties%elements(3)%ionizationStateFractionHistory(2,0: 6))
  allocate(properties%elements(3)%photoIonizationRate           (  0: 5))
  allocate(properties%elements(3)%photoHeatingRate              (  0: 5))
  allocate(properties%elements(3)%photoIonizationRatePrevious   (  0: 5))
  allocate(properties%elements(3)%photoHeatingRatePrevious      (  0: 5))
  allocate(properties%elements(3)%ionizationStateFraction       (  0: 6))
  allocate(properties%elements(4)%photoIonizationRateHistory    (2,0: 7))
  allocate(properties%elements(4)%photoHeatingRateHistory       (2,0: 7))
  allocate(properties%elements(4)%ionizationStateFractionHistory(2,0: 8))
  allocate(properties%elements(4)%photoIonizationRate           (  0: 7))
  allocate(properties%elements(4)%photoHeatingRate              (  0: 7))
  allocate(properties%elements(4)%photoIonizationRatePrevious   (  0: 7))
  allocate(properties%elements(4)%photoHeatingRatePrevious      (  0: 7))
  allocate(properties%elements(4)%ionizationStateFraction       (  0: 8))
  allocate(properties%elements(5)%photoIonizationRateHistory    (2,0:13))
  allocate(properties%elements(5)%photoHeatingRateHistory       (2,0:13))
  allocate(properties%elements(5)%ionizationStateFractionHistory(2,0:14))
  allocate(properties%elements(5)%photoIonizationRate           (  0:13))
  allocate(properties%elements(5)%photoHeatingRate              (  0:13))
  allocate(properties%elements(5)%photoIonizationRatePrevious   (  0:13))
  allocate(properties%elements(5)%photoHeatingRatePrevious      (  0:13))
  allocate(properties%elements(5)%ionizationStateFraction       (  0:14))
  allocate(properties%elements(6)%photoIonizationRateHistory    (2,0:25))
  allocate(properties%elements(6)%photoHeatingRateHistory       (2,0:25))
  allocate(properties%elements(6)%ionizationStateFractionHistory(2,0:26))
  allocate(properties%elements(6)%photoIonizationRate           (  0:25))
  allocate(properties%elements(6)%photoHeatingRate              (  0:25))
  allocate(properties%elements(6)%photoIonizationRatePrevious   (  0:25))
  allocate(properties%elements(6)%photoHeatingRatePrevious      (  0:25))
  allocate(properties%elements(6)%ionizationStateFraction       (  0:26))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.1149822911213856d-09
  properties            %temperature                =  0.8324950624162928d+04
  properties%elements(1)%densityNumber              =  0.3186589640433595d+00
  properties%elements(1)%ionizationStateFraction    =[ 0.9999098414056881d+00, 0.9015859431191938d-04]
  properties%elements(1)%photoIonizationRate        =[ 0.1747521378101921d-37]
  properties%elements(1)%photoHeatingRate           =[ 0.9575027811767581d-54]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.4723915395620179d-38]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.2683198191199928d-54]
  properties%elements(2)%densityNumber              =  0.3109928479196422d-01
  properties%elements(2)%ionizationStateFraction    =[ 0.9999999999953881d+00, 0.4611924046107046d-11, 0.5090522689287785d-29]
  properties%elements(2)%photoIonizationRate        =[ 0.4640764635443423d-37, 0.1306729874579425d-48]
  properties%elements(2)%photoHeatingRate           =[ 0.2461207546850978d-53, 0.6257411038929028d-65]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.1259026054096125d-37, 0.9499091472697100d-52]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.6930043459978893d-54, 0.4823304551800508d-68]
  properties%elements(3)%densityNumber              =  0.8006661007292382d-04
  properties%elements(3)%ionizationStateFraction    =[ 0.9600119852107292d-04, 0.9999039987962151d+00, 0.5263896227772151d-11, 0.2017650901029773d-28, 0.2554732864600911d-44, 0.1244086572343986d-60, 0.2566427460904302d-76]
  properties%elements(3)%photoIonizationRate        =[ 0.8769229530300367d-20, 0.8132037014306334d-38, 0.3589045369706536d-49, 0.1916900705289919d-67, 0.8902110270102893d-85, 0.2059434741778787d-106]
  properties%elements(3)%photoHeatingRate           =[ 0.1610627167355485d-38, 0.4317582005626209d-54, 0.1756644235797546d-65, 0.9244227108705084d-84, 0.1288033403736800d-102, 0.2378904459386037d-124]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.7698042136217909d-20, 0.2273310096784754d-38, 0.2967054187910934d-52, 0.2374095747488914d-70, 0.1599641087652282d-90, 0.9526203446545196d-112]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.1434851143853008d-38, 0.1252650299389780d-54, 0.1538099239625228d-68, 0.1167592697867909d-86, 0.4115322949586260d-108, 0.1194490307947537d-129]
  properties%elements(4)%densityNumber              =  0.1601332201458477d-03
  properties%elements(4)%ionizationStateFraction    =[ 0.9999188686893512d+00, 0.8113131064881596d-04, 0.3324547763203279d-21, 0.1506770532116384d-36, 0.3613923341975327d-49, 0.4711607011287509d-59, 0.9277669614495633d-63, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.2193028958608482d-38, 0.1582852907220823d-42, 0.6668229313108563d-60, 0.2937667868922605d-75, 0.6657323134397482d-88, 0.4902426493232423d-98, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRate           =[ 0.1201931587024160d-54, 0.8070588740021033d-59, 0.3188269769780404d-76, 0.1298871398741979d-91, 0.2554473286279817d-104, 0.1690833243249769d-114, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.6051638467776794d-39, 0.1056869035713972d-44, 0.5674112062342350d-62, 0.1035055752449272d-79, 0.1009349561669927d-97, 0.1272125847550894d-112, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.3438291550933957d-55, 0.5694637722412547d-61, 0.2877259122577769d-78, 0.4876327254303932d-96, 0.4165501282634390d-114, 0.4756201444608914d-129, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%densityNumber              =  0.1134004640624677d-04
  properties%elements(5)%ionizationStateFraction    =[ 0.1800375213445502d-319, 0.5165649506937757d-315, 0.4644217070907718d-321, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.7773736237071671d-314, 0.1363562107607444d-232, 0.1401156150848983d-114, 0.1000000000000000d+01, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoIonizationRate        =[ 0.5824759609678382d-21, 0.3167508523444822d-39, 0.5425454546046776d-50, 0.2039457961487241d-69, 0.2793679692982451d-89, 0.2422973845983362d-109, 0.1109761848275629d-129, 0.3607786859227427d-150, 0.6244050724856305d-171, 0.6480787392769458d-194, 0.7134245870431799d-221, 0.1903660275121619d-163, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRate           =[ 0.1859272377783804d-39, 0.1777355958184574d-55, 0.2895313954269416d-66, 0.1050291399043164d-85, 0.8943335893189123d-106, 0.6269692509382309d-126, 0.2134970904910397d-146, 0.3662414010324437d-167, 0.1588001335245284d-188, 0.6363165642524713d-212, 0.7683354927441095d-239, 0.2053614361779898d-181, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoIonizationRatePrevious=[ 0.1185145474035515d-20, 0.4895745811852008d-39, 0.1714898168054845d-46, 0.4816359026528189d-63, 0.6309650314958056d-80, 0.5297303627447833d-97, 0.1867543616295396d-110, 0.1102419573791060d-114, 0.6589066646832619d-107, 0.2419164523834549d-89, 0.2503814931699532d-75, 0.8498026707056565d-106, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRatePrevious   =[ 0.3783453180326249d-39, 0.2760619377891972d-55, 0.9285645167426375d-63, 0.2518013997680918d-79, 0.2069221329356089d-96, 0.1412157198908340d-113, 0.3738602256243666d-127, 0.1205346761413901d-131, 0.2165981001509157d-124, 0.4851936721553088d-107, 0.2880845986568805d-93, 0.7461176597172840d-124, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%densityNumber              =  0.9215830220638581d-05
  properties%elements(6)%ionizationStateFraction    =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.3651703705127574d-308, 0.2650634252547611d-155, 0.1000000000000000d+01, 0.1594648386177952d-21, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoIonizationRate        =[ 0.2242262449862456d-39, 0.1628694890384073d-41, 0.4741304431177864d-52, 0.3907115703759017d-70, 0.5649420562021145d-90, 0.3872368145192501d-110, 0.1512877334844691d-130, 0.3689086546218649d-151, 0.6129842296590605d-172, 0.6751378777982560d-193, 0.5057747936770733d-214, 0.1808189605067596d-194, 0.6203508945442314d-40, 0.3452132117154760d-63, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRate           =[ 0.1288939485878879d-55, 0.9145291873825075d-58, 0.2552412524969626d-68, 0.1952223861810273d-86, 0.2639987437579580d-106, 0.1660087802695062d-126, 0.5858170093299154d-147, 0.1274543323688591d-167, 0.1307097275509825d-188, 0.1131603025327740d-209, 0.6192476920464600d-231, 0.6574795795461012d-212, 0.7850865736653437d-58, 0.5405880888700499d-81, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoIonizationRatePrevious=[ 0.2242412227029631d-39, 0.1638203913933925d-41, 0.1054806655396607d-50, 0.2605227116046093d-67, 0.2887573957851829d-84, 0.1663096365338628d-101, 0.5016487935261105d-114, 0.5347036737003211d-115, 0.2868094666188384d-100, 0.3612971107118585d-70, 0.9979947550422008d-40, 0.3977298625412024d-58, 0.9047845364672958d-257, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRatePrevious   =[ 0.1289026750467316d-55, 0.9199426906207138d-58, 0.5756889066281677d-67, 0.1321985457014657d-83, 0.1371856398139031d-100, 0.7259179836216655d-118, 0.1981521726211578d-130, 0.1888929258193318d-131, 0.6338838057357936d-117, 0.6336474928977924d-87, 0.1299375089964609d-56, 0.2598693957284839d-75, 0.8444584341919007d-275, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  !! Test 2.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.1570796326794902d-10
  properties            %temperature                =  0.4331819438256507d+04
  properties%elements(1)%densityNumber              =  0.9822759748911494d-01
  properties%elements(1)%ionizationStateFraction    =[ 0.9999997675722596d+00, 0.2324277404218597d-06]
  properties%elements(1)%photoIonizationRate        =[ 0.3149549812326454d-23]
  properties%elements(1)%photoHeatingRate           =[ 0.7370018538232106d-40]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.1343942595755334d-23]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.3137248227211692d-40]
  properties%elements(2)%densityNumber              =  0.9586449381442010d-02
  properties%elements(2)%ionizationStateFraction    =[ 0.9999963719905798d+00, 0.3627395022383571d-05, 0.6143978330395593d-09]
  properties%elements(2)%photoIonizationRate        =[ 0.7363330008294868d-23, 0.1464064510743349d-28]
  properties%elements(2)%photoHeatingRate           =[ 0.1596927019611142d-39, 0.2445419642128941d-45]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.3167859212418551d-23, 0.1906190017275457d-28]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.6864929320508911d-40, 0.3025908388802434d-45]
  properties%elements(3)%densityNumber              =  0.2468077673625692d-04
  properties%elements(3)%ionizationStateFraction    =[ 0.9665523129616057d-09, 0.3400614349091487d-03, 0.1154340240548943d-04, 0.7168758379430924d-05, 0.9996412254377536d+00, 0.4846450778090257d-18, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.3620736090349337d-21, 0.9546894228896723d-29, 0.2975864812121577d-30, 0.1105075610128530d-30, 0.8341327014182825d-37, 0.1033867872326408d-61]
  properties%elements(3)%photoHeatingRate           =[ 0.6561514646842094d-40, 0.2069736759747428d-45, 0.5317540292006181d-47, 0.1679315684478469d-47, 0.9339690958341944d-55, 0.1153662131412653d-79]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.1722891529674162d-17, 0.5762574176970526d-26, 0.6471870329314789d-37, 0.3745956425240425d-43, 0.4125202373915922d-37, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.3135915832099038d-36, 0.1208328801632835d-42, 0.1214106224565391d-53, 0.6031257345089323d-60, 0.4602652478134562d-55, 0.0000000000000000d+00]
  properties%elements(4)%densityNumber              =  0.4936155347251386d-04
  properties%elements(4)%ionizationStateFraction    =[ 0.9984625278605992d+00, 0.1537468966753481d-02, 0.3172644925912438d-08, 0.2403164076055079d-14, 0.3241787176639670d-35, 0.7778714307784842d-45, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.1693469748404474d-24, 0.3845926002732460d-27, 0.7546885277813425d-33, 0.4759983551487899d-25, 0.3315188382422845d-57, 0.5146996353857849d-70, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRate           =[ 0.3952786421195824d-41, 0.7767283112059360d-44, 0.1284265074617931d-49, 0.6678907529377906d-42, 0.2729781683936908d-74, 0.1772829112409763d-87, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.4173153310524982d-25, 0.1406973951396781d-27, 0.2483040136723652d-33, 0.4759983551487864d-25, 0.3311477653916303d-57, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.9524785711063459d-42, 0.2914760801253564d-44, 0.4375945938241276d-50, 0.6678907529377861d-42, 0.2727064113186357d-74, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%densityNumber              =  0.3495603888767819d-05
  properties%elements(5)%ionizationStateFraction    =[ 0.3767560993503195d-08, 0.2645882773037337d+00, 0.8026474517818991d-01, 0.6132234425947926d-02, 0.6490142514052385d+00, 0.4879193290196600d-06, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoIonizationRate        =[ 0.9453493051505720d-22, 0.2164460997783761d-25, 0.6457903464474168d-26, 0.4284089584323139d-27, 0.2201672798531063d-25, 0.3200800089784840d-33, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRate           =[ 0.2975935572935132d-40, 0.5007087973332078d-42, 0.1313602330428169d-42, 0.8060852097889101d-44, 0.4903161905043047d-43, 0.4686795860685526d-51, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoIonizationRatePrevious=[ 0.2428693262642737d-18, 0.1767015518116780d-25, 0.3566923368401594d-30, 0.8869065528141460d-35, 0.1075659203139376d-25, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRatePrevious   =[ 0.7614685799081697d-37, 0.4026114342627846d-42, 0.7582904504466677d-47, 0.1729271267659154d-51, 0.2385016322390415d-43, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%densityNumber              =  0.2840807771275287d-05
  properties%elements(6)%ionizationStateFraction    =[ 0.4060717363064254d-02, 0.3801928042714697d+00, 0.2458387081479202d-01, 0.3432663940908279d-03, 0.5285681421378351d-03, 0.5209656209845447d-12, 0.5902907730139244d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoIonizationRate        =[ 0.1959785105791393d-27, 0.1841098738229581d-25, 0.1197560632215285d-26, 0.1652759480258016d-28, 0.1706238156163757d-25, 0.2050501038832741d-37, 0.1961928233440531d-25, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRate           =[ 0.4783484865014018d-44, 0.4247237438812811d-42, 0.2484566853292507d-43, 0.2791331893832819d-45, 0.2500730897000279d-42, 0.2029910288970880d-54, 0.1137185634567314d-42, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoIonizationRatePrevious=[ 0.5244593854603088d-26, 0.5240237048245339d-26, 0.2807840447414703d-30, 0.1560227340816332d-34, 0.1703900172920338d-25, 0.6323395977280726d-56, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRatePrevious   =[ 0.1257104878377717d-42, 0.1186735446768677d-42, 0.6082106486983901d-47, 0.2776947298548271d-51, 0.2497519897618752d-42, 0.6856408960141033d-73, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #2',status,errorStatusSuccess)
  !! Test 3.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.9504364974660357d-08
  properties            %temperature                =  0.6238572866739802d+07
  properties%elements(1)%densityNumber              =  0.3723988139810234d-01
  properties%elements(1)%ionizationStateFraction    =[ 0.2427359836500234d-08, 0.9999999975726401d+00]
  properties%elements(1)%photoIonizationRate        =[ 0.1743541138531698d-19]
  properties%elements(1)%photoHeatingRate           =[ 0.1505328858502127d-37]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.6270082804677979d-20]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.5425955207281664d-38]
  properties%elements(2)%densityNumber              =  0.3634398551113628d-02
  properties%elements(2)%ionizationStateFraction    =[ 0.9240078219850975d-07, 0.4999999537996089d+00, 0.4999999537996089d+00]
  properties%elements(2)%photoIonizationRate        =[ 0.7089804421292520d-19, 0.1080272036791481d-14]
  properties%elements(2)%photoHeatingRate           =[ 0.7619789452569050d-37, 0.1420948694478112d-32]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.8182833725594046d-20, 0.7857592020093728d-15]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.8852800591361641d-38, 0.1043002735671610d-32]
  properties%elements(3)%densityNumber              =  0.9356934527214736d-05
  properties%elements(3)%ionizationStateFraction    =[ 0.2623161367143534d-113, 0.2684886997385219d-63, 0.1859522832993539d-28, 0.4030384898842573d-10, 0.2829585236170896d-05, 0.4100236063630426d-02, 0.9958969343108296d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.1699246048824414d-97, 0.1481706695762122d-70, 0.6258527934899935d-45, 0.2703339162853195d-28, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRate           =[ 0.2130574892375046d-115, 0.1729357728059584d-88, 0.3997175185746750d-63, 0.3683025785311067d-46, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.8235319818713464d-54, 0.5562055627737527d-41, 0.2006585397990070d-33, 0.1669931234471270d-29, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.1029235766521625d-71, 0.6494630890787865d-59, 0.1099743984582255d-51, 0.2346155955038406d-47, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%densityNumber              =  0.1871386905442948d-04
  properties%elements(4)%ionizationStateFraction    =[ 0.5037486717670505d-205, 0.2895305213181756d-120, 0.1926870046609294d-60, 0.1644674460478032d-25, 0.1903470222279577d-11, 0.5128009355204521d-07, 0.4729289401679540d-03, 0.4816765320856194d-01, 0.9513593665692731d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.3815227279626898d-137, 0.6504997377635548d-101, 0.9110744758362662d-73, 0.1738588773779941d-43, 0.1330023963039691d-31, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRate           =[ 0.7304192152846555d-155, 0.6702926527089620d-119, 0.1333305896397693d-90, 0.2499319234148136d-61, 0.8350091997729027d-50, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.4191236729660200d-59, 0.5096407573286559d-46, 0.2702717311220421d-40, 0.4889818202878533d-36, 0.9739397051581318d-33, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.8035953468366910d-77, 0.5199784396827263d-64, 0.4035614507570237d-58, 0.7253375468967899d-54, 0.6306067407118138d-51, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%densityNumber              =  0.1325247461609598d-05
  properties%elements(5)%ionizationStateFraction    =[ 0.2364451884614657d-301, 0.1804827597549781d-178, 0.6005452650243013d-91, 0.1046911517365954d-40, 0.1072687092854404d-24, 0.3435557565259838d-20, 0.4224542713699439d-16, 0.1632168770011651d-12, 0.2776406670553139d-09, 0.2343212952636804d-06, 0.8660740957046704d-04, 0.1608597335836803d-01, 0.9113469004732619d+00, 0.7103473394893621d-01, 0.1445550210764179d-02]
  properties%elements(5)%photoIonizationRate        =[ 0.6348087059246838d-196, 0.6486762311990663d-148, 0.6756134237414268d-103, 0.1509827334442576d-58, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRate           =[ 0.3318799356463403d-214, 0.9361916390334283d-166, 0.7805448412920390d-121, 0.1125653366146218d-76, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoIonizationRatePrevious=[ 0.1918191142396492d-79, 0.1157717293852240d-65, 0.3107435624783703d-56, 0.4160377745594081d-51, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(5)%photoHeatingRatePrevious   =[ 0.9874522343524975d-98, 0.1680158312549026d-83, 0.3548652196368760d-74, 0.2815394247549216d-69, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%densityNumber              =  0.1077002259867166d-05
  properties%elements(6)%ionizationStateFraction    =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.6277779114638693d-232, 0.2291042341416442d-113, 0.9735175309625375d-50, 0.2198673006021114d-35, 0.8395701347830732d-31, 0.1527955012554552d-26, 0.8757218338965942d-23, 0.2947297891201133d-19, 0.5560663108009888d-16, 0.5357665077219562d-13, 0.3004026941338127d-10, 0.9393610848182059d-08, 0.1538582243196779d-05, 0.7631497850316645d-04, 0.1677252492227540d-02, 0.1946645599275418d-01, 0.1216270174511344d+00, 0.3302895769527080d+00, 0.3735325448969531d+00, 0.1353213397354907d+00, 0.1742684346592096d-01, 0.5811060245154114d-03, 0.3844591619283513d-11, 0.4427089070696599d-20]
  properties%elements(6)%photoIonizationRate        =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.5662899972933536d-255, 0.1014739777455875d-191, 0.7155212943055485d-132, 0.1814898531515846d-69, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRate           =[ 0.0000000000000000d+00, 0.0000000000000000d+00, 0.8240309535311761d-273, 0.1719296506719612d-209, 0.9606808295073146d-150, 0.1886558277956047d-87, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoIonizationRatePrevious=[ 0.1948737248517191d-104, 0.2495908895529311d-92, 0.5703371289781859d-81, 0.8342928017805746d-76, 0.5813668927397507d-71, 0.7279428559000560d-67, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  properties%elements(6)%photoHeatingRatePrevious   =[ 0.2100822070811749d-121, 0.2746176909876323d-109, 0.8259402434055722d-99, 0.1448306808781059d-93, 0.8173234468847486d-89, 0.7855646361797616d-85, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00, 0.0000000000000000d+00]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #3',status,errorStatusSuccess)
  !! Test 4.
  properties            %iterationCount             = 1
  properties            %volume                     =  0.8143427037125286d-07
  properties            %temperature                =  0.1549276837919467d+10
  properties%elements(1)%densityNumber              =  0.4006318280879603d-07
  properties%elements(1)%ionizationStateFraction    =[ 0.9658068224642371d-12, 0.9999999999990342d+00]
  properties%elements(1)%photoIonizationRate        =[ 0.1995352020950228d-33]
  properties%elements(1)%photoHeatingRate           =[ 0.1974358486496125d-51]
  properties%elements(1)%photoIonizationRatePrevious=[ 0.1850281941656038d-29]
  properties%elements(1)%photoHeatingRatePrevious   =[ 0.1882785355789361d-47]
  properties%elements(2)%densityNumber              =  0.3909936554220832d-08
  properties%elements(2)%ionizationStateFraction    =[ 0.2026117084999547d-16, 0.1985247394529095d-07, 0.9999999801475260d+00]
  properties%elements(2)%photoIonizationRate        =[ 0.5411155491129962d-39, 0.3064876087982742d-32]
  properties%elements(2)%photoHeatingRate           =[ 0.6438386028442484d-57, 0.3936805444852026d-50]
  properties%elements(2)%photoIonizationRatePrevious=[ 0.9510859442068140d-35, 0.7127729393119092d-30]
  properties%elements(2)%photoHeatingRatePrevious   =[ 0.1157283516681535d-52, 0.8980778838247629d-48]
  properties%elements(3)%densityNumber              =  0.1006632041832555d-10
  properties%elements(3)%ionizationStateFraction    =[ 0.1850433732270316d-48, 0.6695388453193269d-38, 0.4297950753316255d-28, 0.2737398116211722d-19, 0.1382520977158208d-11, 0.3737138182054625d-05, 0.9999962628604354d+00]
  properties%elements(3)%photoIonizationRate        =[ 0.1493218462607758d-71, 0.3789174653314264d-63, 0.8837953308680280d-55, 0.2356771587328955d-47, 0.1451622566923310d-41, 0.2139907856059930d-35]
  properties%elements(3)%photoHeatingRate           =[ 0.3426067592422157d-90, 0.4900270577145279d-81, 0.9885426670261587d-73, 0.1806149095165274d-64, 0.6446562409872093d-57, 0.9259746271191389d-51]
  properties%elements(3)%photoIonizationRatePrevious=[ 0.2422306821536439d-50, 0.8319343645764637d-44, 0.6206842264203831d-39, 0.5125323630055196d-35, 0.1038012593310584d-32, 0.4662553624227832d-31]
  properties%elements(3)%photoHeatingRatePrevious   =[ 0.5593793357908285d-69, 0.1101962797456126d-61, 0.7883581994849857d-57, 0.3700140887164601d-52, 0.4575908501287427d-48, 0.1978251538211450d-46]
  properties%elements(4)%densityNumber              =  0.2013264083665111d-10
  properties%elements(4)%ionizationStateFraction    =[ 0.3610597581516900d-55, 0.2393965531311841d-45, 0.3002397502549896d-37, 0.1239930359363018d-29, 0.9342629006722523d-23, 0.1106623891865137d-15, 0.1259683907922584d-09, 0.2229703337667918d-04, 0.9999777028406548d+00]
  properties%elements(4)%photoIonizationRate        =[ 0.1991511901665005d-79, 0.2337508942028020d-70, 0.8116032579875602d-64, 0.1811773218459454d-57, 0.7538936012834717d-52, 0.7360081924449229d-45, 0.7767779298160135d-39, 0.7648691831308045d-34]
  properties%elements(4)%photoHeatingRate           =[ 0.4274686059787210d-97, 0.2863913588602610d-88, 0.2209617999367187d-81, 0.4666073794417764d-74, 0.3332939956171900d-67, 0.3808865539451950d-60, 0.3342219366055016d-54, 0.3110566194099163d-49]
  properties%elements(4)%photoIonizationRatePrevious=[ 0.8203601570111106d-56, 0.9315442806533606d-49, 0.4781139074320958d-44, 0.2050401787679102d-40, 0.2391270031949154d-37, 0.1126043709255960d-33, 0.1916557200637668d-30, 0.9150973747201850d-30]
  properties%elements(4)%photoHeatingRatePrevious   =[ 0.1782079748781783d-73, 0.1118162380238818d-66, 0.1245629222662095d-61, 0.5656432673705159d-57, 0.1053694649886003d-52, 0.5843215172160241d-49, 0.8212721941294146d-46, 0.3725554195297447d-45]
  properties%elements(5)%densityNumber              =  0.1425719667411824d-11
  properties%elements(5)%ionizationStateFraction    =[ 0.6340711790993047d-84, 0.1443671193258700d-72, 0.8742994166396034d-64, 0.1243594109258829d-56, 0.4752967705066366d-48, 0.6113605956848655d-42, 0.2465397184676559d-36, 0.6274336147074745d-31, 0.3966210938230325d-26, 0.1633601671397108d-21, 0.7206544266937757d-17, 0.2168489775304743d-12, 0.1804651113795612d-07, 0.2557434417198595d-03, 0.9997442385115521d+00]
  properties%elements(5)%photoIonizationRate        =[ 0.4302335828196386d-107, 0.3330556734823719d-99, 0.3675624312309984d-91, 0.1037038069923369d-84, 0.1657722324207264d-77, 0.2112184785563256d-71, 0.8474724472048651d-66, 0.2144387533280919d-60, 0.1337242632021735d-55, 0.5424783597111897d-51, 0.2338053787603392d-46, 0.6648965385671366d-42, 0.4901747625213514d-37, 0.3426537710026839d-33]
  properties%elements(5)%photoHeatingRate           =[ 0.1231701353986665d-125, 0.8865368577662103d-117, 0.2306987246876700d-108, 0.2680909269257332d-101, 0.9446404603991582d-93, 0.1195534118613604d-86, 0.4761009281895888d-81, 0.1192163894736258d-75, 0.7389119773295600d-71, 0.2986999503998836d-66, 0.1274775445338553d-61, 0.3626221167967652d-57, 0.1243679316031691d-52, 0.8294296719553294d-49]
  properties%elements(5)%photoIonizationRatePrevious=[ 0.3645336678099931d-70, 0.4237381034815603d-64, 0.4582715005654088d-58, 0.2802335208994288d-53, 0.1741478645363542d-49, 0.9110570811835007d-46, 0.1723958098906676d-42, 0.9302997344138548d-40, 0.2014411758082995d-37, 0.2096045060342091d-35, 0.1036409756039569d-33, 0.3360870592501588d-32, 0.6141850402512858d-31, 0.2998671282534678d-30]
  properties%elements(5)%photoHeatingRatePrevious   =[ 0.1009636387349113d-88, 0.1170495217662300d-81, 0.3068524172999380d-75, 0.8010811564836488d-70, 0.9864403269544908d-65, 0.5124440821476048d-61, 0.9668569534835932d-58, 0.5162078358126769d-55, 0.1110768764903739d-52, 0.1151875516670028d-50, 0.5642946262793183d-49, 0.1833612151381285d-47, 0.1559581600865099d-46, 0.7247909052713006d-46]
  properties%elements(6)%densityNumber              =  0.1158654023660329d-11
  properties%elements(6)%ionizationStateFraction    =[ 0.1517906829866800d-137, 0.1655015158589035d-127, 0.5229142870351517d-119, 0.2512852832737620d-110, 0.1643651371303520d-103, 0.5774828151025926d-97, 0.1018026401202543d-90, 0.1180182783514059d-84, 0.7129860407192726d-79, 0.2800311437027643d-73, 0.8018041783387322d-68, 0.1374181925747399d-62, 0.1381496020623042d-57, 0.9405987759708819d-53, 0.2442679037676192d-47, 0.3691287247919358d-42, 0.4265926688286993d-37, 0.3000124176320092d-32, 0.1304473334675312d-27, 0.3966310740667937d-23, 0.8224793000028341d-19, 0.1063207636862050d-14, 0.6515818528873430d-11, 0.2881647152305536d-07, 0.5464121743744143d-04, 0.1465753747831528d-01, 0.9852877924812589d+00]
  properties%elements(6)%photoIonizationRate        =[ 0.1868999670878200d-165, 0.7325498462801437d-156, 0.5080951135620962d-145, 0.7981981361098975d-138, 0.6789077664319846d-132, 0.4058305432375498d-126, 0.4589053713027022d-120, 0.5258780022018362d-114, 0.3173325379895824d-108, 0.1234724794072276d-102, 0.3499736499853008d-97, 0.5940582403150359d-92, 0.5896937031960514d-87, 0.3965499915323305d-82, 0.1012096026809569d-76, 0.1480261425234796d-71, 0.1668366683785624d-66, 0.1097364040056853d-61, 0.4406285587380056d-57, 0.1230517592034127d-52, 0.2269945309140933d-48, 0.2493187492460336d-44, 0.1263193162552482d-40, 0.3463648353231874d-37, 0.2407745834358789d-34, 0.2763380210959416d-32]
  properties%elements(6)%photoHeatingRate           =[ 0.7058072489927640d-182, 0.6296483378560901d-172, 0.1008146055528818d-162, 0.9665430221517006d-155, 0.5541150353961255d-148, 0.1909596410690022d-141, 0.3345875232527808d-135, 0.3856578049193219d-129, 0.2280922706976698d-123, 0.8828806021268265d-118, 0.2487574821938886d-112, 0.4186378398699912d-107, 0.4129582003866029d-102, 0.2758452491988906d-97, 0.6944717438230551d-92, 0.1012790311070455d-86, 0.9408200064361346d-82, 0.6088588449027779d-77, 0.2395967545912135d-72, 0.6689156687239744d-68, 0.1247811783257415d-63, 0.1409776925197841d-59, 0.7540267348206002d-56, 0.2492226187532916d-52, 0.5438526139092387d-50, 0.5656203105614208d-48]
  properties%elements(6)%photoIonizationRatePrevious=[ 0.3160581500464020d-94, 0.1141424423749116d-86, 0.1409975173620722d-77, 0.1470365339221957d-72, 0.2800958297835990d-68, 0.1669233592997055d-64, 0.2268839777937830d-60, 0.2057176405455841d-56, 0.7229243102915538d-53, 0.6079730386198588d-50, 0.2739685993772693d-47, 0.6686120111687445d-45, 0.8121110104848075d-43, 0.5985089540586486d-41, 0.2511314953586159d-39, 0.6731060586480003d-38, 0.1209996393579578d-36, 0.2435378293607639d-35, 0.3955222039623559d-34, 0.5187185033989526d-33, 0.4808372275582017d-32, 0.3136944719543164d-31, 0.1431387935058248d-30, 0.3450410675056916d-30, 0.1946361825860747d-30, 0.3437609050817295d-31]
  properties%elements(6)%photoHeatingRatePrevious   =[ 0.1146530317157367d-110, 0.1038590353849624d-102, 0.2872392471634188d-95, 0.1659077196016000d-89, 0.2315385119728523d-84, 0.8444916695970308d-80, 0.1658128073080524d-75, 0.1503364182969386d-71, 0.5186191764014226d-68, 0.4340273804164533d-65, 0.1946995458259745d-62, 0.4710969535727516d-60, 0.5689248912762104d-58, 0.4168387665766990d-56, 0.1725308811762099d-54, 0.4611030427058273d-53, 0.6828269506811836d-52, 0.1351983294074198d-50, 0.2155071466349329d-49, 0.2820226300695007d-48, 0.2645941725137840d-47, 0.1774824845065237d-46, 0.8538609609499634d-46, 0.2482422775225330d-45, 0.4402951431747823d-46, 0.7043203683847899d-47]
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #4',status,errorStatusSuccess)
  ! Done with atomic matter - H, He, C, O, Si, Fe.
  call Unit_Tests_End_Group()
  deallocate(properties)
  !![
  <objectDestructor name="radiativeTransferMatter_"                      />
  <objectDestructor name="atomicCrossSectionIonizationPhoto_"            />
  <objectDestructor name="atomicIonizationPotential_"                    />
  <objectDestructor name="atomicExcitationRateCollisional_"              />
  <objectDestructor name="gauntFactor_"                                  />
  <objectDestructor name="atomicRecombinationRateDielectronicZero_"      />
  <objectDestructor name="atomicRecombinationRateDielectronicArnaud1985_"/>
  <objectDestructor name="atomicRecombinationRateRadiativeFixed_"        />
  <objectDestructor name="atomicRecombinationRateRadiativeVerner1996_"   />
  <objectDestructor name="atomicRecombinationRateRadiativeCoolingHummer_"/>
  <objectDestructor name="atomicIonizationRateCollisional_"              />
  <objectDestructor name="massDistribution_"                             />
  !!]
  ! Done with unit tests.
  call Unit_Tests_Finish()
  call Node_Components_Uninitialize()
  call parameters%destroy()
#ifdef USEMPI
  call mpiFinalize()
#endif
end program Test_Radiative_Transfer_State_Solver
