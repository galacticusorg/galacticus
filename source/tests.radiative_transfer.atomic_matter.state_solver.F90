!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a code to test the radiative transfer state solver.

program Test_Radiative_Transfer_State_Solver
  !% Test the radiative transfer state solver.
  use :: Atomic_Cross_Sections_Ionization_Photo , only : atomicCrossSectionIonizationPhotoVerner
  use :: Atomic_Ionization_Potentials           , only : atomicIonizationPotentialVerner
  use :: Atomic_Radiation_Gaunt_Factors         , only : gauntFactorVanHoof2014
  use :: Atomic_Rates_Excitation_Collisional    , only : atomicExcitationRateCollisionalScholzWalters1991
  use :: Atomic_Rates_Ionization_Collisional    , only : atomicIonizationRateCollisionalVerner1996
  use :: Atomic_Rates_Recombination_Dielectronic, only : atomicRecombinationRateDielectronicZero
  use :: Atomic_Rates_Recombination_Radiative   , only : atomicRecombinationRateRadiativeFixed           , atomicRecombinationRateRadiativeVerner1996
  use :: Galacticus_Display                     , only : Galacticus_Verbosity_Level_Set                  , verbosityStandard
  use :: Galacticus_Error                       , only : errorStatusSuccess
  use :: Input_Parameters                       , only : inputParameters
  use :: ISO_Varying_String                     , only : var_str
  use :: Mass_Distributions                     , only : massDistributionConstantDensityCloud
  use :: Numerical_Constants_Astronomical       , only : metallicitySolar
  use :: Radiative_Transfer_Matters             , only : radiativeTransferMatterAtomic                   , radiativeTransferMatterPropertiesAtomic
  use :: Unit_Tests                             , only : Assert                                          , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group, Unit_Tests_Finish
#ifdef USEMPI
  use :: MPI                                    , only : MPI_Thread_Single
  use :: MPI_Utilities                          , only : mpiFinalize                                     , mpiInitialize
#endif
  implicit none
  type   (inputParameters                                 ), target      :: parameters
  type   (atomicCrossSectionIonizationPhotoVerner         ), pointer     :: atomicCrossSectionIonizationPhoto_
  type   (atomicIonizationPotentialVerner                 ), pointer     :: atomicIonizationPotential_
  type   (atomicExcitationRateCollisionalScholzWalters1991), pointer     :: atomicExcitationRateCollisional_
  type   (gauntFactorVanHoof2014                          ), pointer     :: gauntFactor_
  type   (atomicIonizationRateCollisionalVerner1996       ), pointer     :: atomicIonizationRateCollisional_
  type   (atomicRecombinationRateDielectronicZero         ), pointer     :: atomicRecombinationRateDielectronic_
  type   (atomicRecombinationRateRadiativeFixed           ), pointer     :: atomicRecombinationRateRadiativeFixed_
  type   (atomicRecombinationRateRadiativeVerner1996      ), pointer     :: atomicRecombinationRateRadiativeVerner1996_
  type   (massDistributionConstantDensityCloud            ), pointer     :: massDistribution_
  type   (radiativeTransferMatterAtomic                   ), pointer     :: radiativeTransferMatter_
  type   (radiativeTransferMatterPropertiesAtomic         ), allocatable :: properties
  integer                                                                :: status

#ifdef USEMPI
  call mpiInitialize(MPI_Thread_Single  )
#endif
  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Initialize parameters.
  parameters=inputParameters(var_str('testSuite/parameters/test-radiativeTransfer-atomicMatterStateSolver.xml'))
  call parameters%markGlobal()
  ! Construct atomic matter.
  allocate(atomicCrossSectionIonizationPhoto_         )
  allocate(atomicIonizationPotential_                 )
  allocate(atomicExcitationRateCollisional_           )
  allocate(gauntFactor_                               )
  allocate(atomicIonizationRateCollisional_           )
  allocate(atomicRecombinationRateDielectronic_       )
  allocate(atomicRecombinationRateRadiativeFixed_     )
  allocate(atomicRecombinationRateRadiativeVerner1996_)
  allocate(massDistribution_                          )
  !# <referenceConstruct object="atomicCrossSectionIonizationPhoto_"  >
  !#  <constructor>
  !#   atomicCrossSectionIonizationPhotoVerner         (                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicIonizationPotential_"          >
  !#  <constructor>
  !#   atomicIonizationPotentialVerner                 (                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicExcitationRateCollisional_"    >
  !#  <constructor>
  !#   atomicExcitationRateCollisionalScholzWalters1991(                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="gauntFactor_"                        >
  !#  <constructor>
  !#   gauntFactorVanHoof2014                          (                                                       &amp;
  !#     &amp;                                          atomicIonizationPotential_=atomicIonizationPotential_  &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicRecombinationRateDielectronic_">
  !#  <constructor>
  !#   atomicRecombinationRateDielectronicZero         (                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicRecombinationRateRadiativeFixed_">
  !#  <constructor>
  !#   atomicRecombinationRateRadiativeFixed           (                                                       &amp;
  !#     &amp;                                          rateCoefficient           =2.0d-13                     &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicRecombinationRateRadiativeVerner1996_">
  !#  <constructor>
  !#   atomicRecombinationRateRadiativeVerner1996      (                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="atomicIonizationRateCollisional_"    >
  !#  <constructor>
  !#   atomicIonizationRateCollisionalVerner1996       (                                                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="massDistribution_"                   >
  !#  <constructor>
  !#   massDistributionConstantDensityCloud            (                                                       &amp;
  !#     &amp;                                          mass                      =1.0d0                     , &amp;
  !#     &amp;                                          radius                    =1.0d0                       &amp;
  !#     &amp;                                         )
  !#  </constructor>
  !# </referenceConstruct>
  ! Tests on atomic matter - H only.
  call Unit_Tests_Begin_Group("Atomic matter - H only")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_)
  allocate(properties              )
  !# <referenceConstruct object="radiativeTransferMatter_">
  !#  <constructor>
  !#   radiativeTransferMatterAtomic(                                                                             &amp;
  !#     &amp;                       abundancePattern                    =var_str('solar')                      , &amp;
  !#     &amp;                       metallicity                         =0.0d0                                 , &amp;
  !#     &amp;                       elements                            =['H ']                                , &amp;
  !#     &amp;                       iterationAverageCount               =1                                     , &amp;
  !#     &amp;                       temperatureMinimum                  =3.0d0                                 , &amp;
  !#     &amp;                       massDistribution_                   =massDistribution_                     , &amp;
  !#     &amp;                       atomicCrossSectionIonizationPhoto_  =atomicCrossSectionIonizationPhoto_    , &amp;
  !#     &amp;                       atomicRecombinationRateRadiative_   =atomicRecombinationRateRadiativeFixed_, &amp;
  !#     &amp;                       atomicIonizationRateCollisional_    =atomicIonizationRateCollisional_      , &amp;
  !#     &amp;                       atomicRecombinationRateDielectronic_=atomicRecombinationRateDielectronic_  , &amp;
  !#     &amp;                       atomicIonizationPotential_          =atomicIonizationPotential_            , &amp;
  !#     &amp;                       atomicExcitationRateCollisional_    =atomicExcitationRateCollisional_      ,  &amp;
  !#     &amp;                       gauntFactor_                        =gauntFactor_                            &amp;
  !#     &amp;                      )
  !#  </constructor>
  !# </referenceConstruct>
  allocate(properties%elements(1)                                   )
  allocate(properties%elements(1)%photoIonizationRateHistory (1,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory    (1,0:0))
  allocate(properties%elements(1)%photoIonizationRate        (  0:0))
  allocate(properties%elements(1)%photoHeatingRate           (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious(  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious   (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction    (  0:1))
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
  !# <objectDestructor name="radiativeTransferMatter_"/>
  ! Tests on atomic matter - H, He, O.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, O (fixed recombination coefficient)")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_)
  allocate(properties              )
  !# <referenceConstruct object="radiativeTransferMatter_">
  !#  <constructor>
  !#   radiativeTransferMatterAtomic(                                                                             &amp;
  !#     &amp;                       abundancePattern                    =var_str('solar')                      , &amp;
  !#     &amp;                       metallicity                         =1.0d0                                 , &amp;
  !#     &amp;                       elements                            =['H ','He','O ']                      , &amp;
  !#     &amp;                       iterationAverageCount               =2                                     , &amp;
  !#     &amp;                       temperatureMinimum                  =3.0d0                                 , &amp;
  !#     &amp;                       massDistribution_                   =massDistribution_                     , &amp;
  !#     &amp;                       atomicCrossSectionIonizationPhoto_  =atomicCrossSectionIonizationPhoto_    , &amp;
  !#     &amp;                       atomicRecombinationRateRadiative_   =atomicRecombinationRateRadiativeFixed_, &amp;
  !#     &amp;                       atomicIonizationRateCollisional_    =atomicIonizationRateCollisional_      , &amp;
  !#     &amp;                       atomicRecombinationRateDielectronic_=atomicRecombinationRateDielectronic_  , &amp;
  !#     &amp;                       atomicIonizationPotential_          =atomicIonizationPotential_            , &amp;
  !#     &amp;                       atomicExcitationRateCollisional_    =atomicExcitationRateCollisional_      , &amp;
  !#     &amp;                       gauntFactor_                        =gauntFactor_                            &amp;
  !#     &amp;                      )
  !#  </constructor>
  !# </referenceConstruct>
  allocate(properties%elements(3)                                   )
  allocate(properties%elements(1)%photoIonizationRateHistory (1,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory    (1,0:0))
  allocate(properties%elements(1)%photoIonizationRate        (  0:0))
  allocate(properties%elements(1)%photoHeatingRate           (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious(  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious   (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction    (  0:1))
  allocate(properties%elements(2)%photoIonizationRateHistory (1,0:1))
  allocate(properties%elements(2)%photoHeatingRateHistory    (1,0:1))
  allocate(properties%elements(2)%photoIonizationRate        (  0:1))
  allocate(properties%elements(2)%photoHeatingRate           (  0:1))
  allocate(properties%elements(2)%photoIonizationRatePrevious(  0:1))
  allocate(properties%elements(2)%photoHeatingRatePrevious   (  0:1))
  allocate(properties%elements(2)%ionizationStateFraction    (  0:2))
  allocate(properties%elements(3)%photoIonizationRateHistory (1,0:7))
  allocate(properties%elements(3)%photoHeatingRateHistory    (1,0:7))
  allocate(properties%elements(3)%photoIonizationRate        (  0:7))
  allocate(properties%elements(3)%photoHeatingRate           (  0:7))
  allocate(properties%elements(3)%photoIonizationRatePrevious(  0:7))
  allocate(properties%elements(3)%photoHeatingRatePrevious   (  0:7))
  allocate(properties%elements(3)%ionizationStateFraction    (  0:8))
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
  !# <objectDestructor name="radiativeTransferMatter_"/>
  ! Tests on atomic matter - H, He, O.
  call Unit_Tests_Begin_Group("Atomic matter - H, He, O")
  ! Initialize the properties of the solver and the atomic matter.
  allocate(radiativeTransferMatter_)
  allocate(properties              )
  !# <referenceConstruct object="radiativeTransferMatter_">
  !#  <constructor>
  !#   radiativeTransferMatterAtomic(                                                                           &amp;
  !#     &amp;                       abundancePattern                    =var_str('solar')                    , &amp;
  !#     &amp;                       metallicity                         =1.0d0                               , &amp;
  !#     &amp;                       elements                            =['H ','He','O ']                    , &amp;
  !#     &amp;                       iterationAverageCount               =2                                   , &amp;
  !#     &amp;                       temperatureMinimum                  =3.0d0                               , &amp;
  !#     &amp;                       massDistribution_                   =massDistribution_                   , &amp;
  !#     &amp;                       atomicCrossSectionIonizationPhoto_  =atomicCrossSectionIonizationPhoto_  , &amp;
  !#     &amp;                       atomicRecombinationRateRadiative_   =atomicRecombinationRateRadiativeVerner1996_, &amp;
  !#     &amp;                       atomicIonizationRateCollisional_    =atomicIonizationRateCollisional_    , &amp;
  !#     &amp;                       atomicRecombinationRateDielectronic_=atomicRecombinationRateDielectronic_, &amp;
  !#     &amp;                       atomicIonizationPotential_          =atomicIonizationPotential_          , &amp;
  !#     &amp;                       atomicExcitationRateCollisional_    =atomicExcitationRateCollisional_    , &amp;
  !#     &amp;                       gauntFactor_                        =gauntFactor_                          &amp;
  !#     &amp;                      )
  !#  </constructor>
  !# </referenceConstruct>
  allocate(properties%elements(3)                                   )
  allocate(properties%elements(1)%photoIonizationRateHistory (1,0:0))
  allocate(properties%elements(1)%photoHeatingRateHistory    (1,0:0))
  allocate(properties%elements(1)%photoIonizationRate        (  0:0))
  allocate(properties%elements(1)%photoHeatingRate           (  0:0))
  allocate(properties%elements(1)%photoIonizationRatePrevious(  0:0))
  allocate(properties%elements(1)%photoHeatingRatePrevious   (  0:0))
  allocate(properties%elements(1)%ionizationStateFraction    (  0:1))
  allocate(properties%elements(2)%photoIonizationRateHistory (1,0:1))
  allocate(properties%elements(2)%photoHeatingRateHistory    (1,0:1))
  allocate(properties%elements(2)%photoIonizationRate        (  0:1))
  allocate(properties%elements(2)%photoHeatingRate           (  0:1))
  allocate(properties%elements(2)%photoIonizationRatePrevious(  0:1))
  allocate(properties%elements(2)%photoHeatingRatePrevious   (  0:1))
  allocate(properties%elements(2)%ionizationStateFraction    (  0:2))
  allocate(properties%elements(3)%photoIonizationRateHistory (1,0:7))
  allocate(properties%elements(3)%photoHeatingRateHistory    (1,0:7))
  allocate(properties%elements(3)%photoIonizationRate        (  0:7))
  allocate(properties%elements(3)%photoHeatingRate           (  0:7))
  allocate(properties%elements(3)%photoIonizationRatePrevious(  0:7))
  allocate(properties%elements(3)%photoHeatingRatePrevious   (  0:7))
  allocate(properties%elements(3)%ionizationStateFraction    (  0:8))
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
  !# <objectDestructor name="radiativeTransferMatter_"/>
  ! Done with unit tests.
  call Unit_Tests_Finish()
  call parameters%destroy()
#ifdef USEMPI
  call mpiFinalize()
#endif
end program Test_Radiative_Transfer_State_Solver
