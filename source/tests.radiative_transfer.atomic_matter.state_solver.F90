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
  use :: Atomic_Rates_Recombination_Radiative   , only : atomicRecombinationRateRadiativeFixed
  use :: Galacticus_Display                     , only : Galacticus_Verbosity_Level_Set                  , verbosityStandard
  use :: Galacticus_Error                       , only : errorStatusSuccess
  use :: Mass_Distributions                     , only : massDistributionConstantDensityCloud
  use :: Radiative_Transfer_Matters             , only : radiativeTransferMatterAtomic                   , radiativeTransferMatterPropertiesAtomic
  use :: Unit_Tests                             , only : Assert                                          , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group, Unit_Tests_Finish
#ifdef USEMPI
  use :: MPI                                    , only : MPI_Thread_Single
  use :: MPI_Utilities                          , only : mpiFinalize                                     , mpiInitialize
#endif
  implicit none
  type   (atomicCrossSectionIonizationPhotoVerner         ) :: atomicCrossSectionIonizationPhoto_
  type   (atomicIonizationPotentialVerner                 ) :: atomicIonizationPotential_
  type   (atomicExcitationRateCollisionalScholzWalters1991) :: atomicExcitationRateCollisional_
  type   (gauntFactorVanHoof2014                          ) :: gauntFactor_
  type   (atomicIonizationRateCollisionalVerner1996       ) :: atomicIonizationRateCollisional_
  type   (atomicRecombinationRateDielectronicZero         ) :: atomicRecombinationRateDielectronic_
  type   (atomicRecombinationRateRadiativeFixed           ) :: atomicRecombinationRateRadiative_
  type   (massDistributionConstantDensityCloud            ) :: massDistribution_
  type   (radiativeTransferMatterAtomic                   ) :: radiativeTransferMatter_
  type   (radiativeTransferMatterPropertiesAtomic         ) :: properties
  integer                                                   :: status

#ifdef USEMPI
  call mpiInitialize(MPI_Thread_Single  )
#endif
  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Construct atomic matter.
  atomicCrossSectionIonizationPhoto_  =atomicCrossSectionIonizationPhotoVerner         (                                                                           &
       &                                                                               )
  atomicIonizationPotential_          =atomicIonizationPotentialVerner                 (                                                                           &
       &                                                                               )
  atomicExcitationRateCollisional_    =atomicExcitationRateCollisionalScholzWalters1991(                                                                           &
       &                                                                               )
  gauntFactor_                        =gauntFactorVanHoof2014                          (                                                                           &
       &                                                                                atomicIonizationPotential_          =atomicIonizationPotential_            &
       &                                                                               )
  atomicRecombinationRateDielectronic_=atomicRecombinationRateDielectronicZero         (                                                                           &
       &                                                                               )
  atomicRecombinationRateRadiative_   =atomicRecombinationRateRadiativeFixed           (                                                                           &
       &                                                                                rateCoefficient                     =2.0d-13                               &
       &                                                                               )
  atomicIonizationRateCollisional_    =atomicIonizationRateCollisionalVerner1996       (                                                                           &
       &                                                                               )
  massDistribution_                   =massDistributionConstantDensityCloud            (                                                                           &
       &                                                                                mass                                =1.0d0                               , &
       &                                                                                radius                              =1.0d0                                 &
       &                                                                               )
  radiativeTransferMatter_            =radiativeTransferMatterAtomic                   (                                                                           &
       &                                                                                iterationAverageCount               =5                                   , &
       &                                                                                temperatureMinimum                  =3.0d0                               , &
       &                                                                                massDistribution_                   =massDistribution_                   , &
       &                                                                                atomicCrossSectionIonizationPhoto_  =atomicCrossSectionIonizationPhoto_  , &
       &                                                                                atomicRecombinationRateRadiative_   =atomicRecombinationRateRadiative_   , &
       &                                                                                atomicIonizationRateCollisional_    =atomicIonizationRateCollisional_    , &
       &                                                                                atomicRecombinationRateDielectronic_=atomicRecombinationRateDielectronic_, &
       &                                                                                atomicIonizationPotential_          =atomicIonizationPotential_          , &
       &                                                                                atomicExcitationRateCollisional_    =atomicExcitationRateCollisional_    , &
       &                                                                                gauntFactor_                        =gauntFactor_                          &
       &                                                                               )
  ! Tests on atomic matter.
  call Unit_Tests_Begin_Group("Atomic matter")
  ! Initialize the properties of the atomic matter.
  allocate(properties%photoIonizationRateHistory(5,0:0))
  allocate(properties%photoHeatingRateHistory   (5,0:0))
  ! Perform tests. These are cases which have previously failed while running radiative transfer models due to inadequacies of the
  ! solver.  
  !! Test 1.
  properties%iterationCount             = 1
  properties%densityNumber              = 0.101064330490d+03
  properties%volume                     = 0.172800000000d-20
  properties%temperature                = 0.103698470589d+05
  properties%ionizationStateFraction    =[0.706879980476d+00,0.293120019524d+00]
  properties%photoIonizationRate        = 0.209229241840d-06
  properties%photoHeatingRate           = 0.331582285728d-24
  properties%photoIonizationRatePrevious= 0.228307228801d-09
  properties%photoHeatingRatePrevious   = 0.164662972542d-26
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #1',status,errorStatusSuccess)
  !! Test 2.
  properties%iterationCount             = 1
  properties%densityNumber              = 0.101064330490d+03
  properties%volume                     = 0.172800000000d-20
  properties%temperature                = 0.100000000000d+05
  properties%ionizationStateFraction    =[0.100000000000d+01,0.000000000000d+00]
  properties%photoIonizationRate        = 0.245726973297d-04
  properties%photoHeatingRate           = 0.330670408342d-22
  properties%photoIonizationRatePrevious=-huge(0.0d0)
  properties%photoHeatingRatePrevious   =-huge(0.0d0)
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #2',status,errorStatusSuccess)
  !! Test 3.
  properties%iterationCount             = 1
  properties%densityNumber              = 0.93764865414000000d-02
  properties%volume                     = 0.17280000000000000d-20
  properties%temperature                = 0.30000000000000000d+01
  properties%ionizationStateFraction    =[0.10000000000000000d+01,0.000000000000d+00]
  properties%photoIonizationRate        = 0.87630070734000000d-22
  properties%photoHeatingRate           = 0.15147475879338856d-38
  properties%photoIonizationRatePrevious= 0.00000000000000000d+00
  properties%photoHeatingRatePrevious   = 0.00000000000000000d+00
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #3',status,errorStatusSuccess)
  !! Test 4.
  properties%iterationCount             = 1
  properties%densityNumber              = 0.1010643304897267d+03
  properties%volume                     = 0.1727999999999987d-20
  properties%temperature                = 0.8464271217162063d+05
  properties%ionizationStateFraction    =[0.4697602133585263d-09,0.9999999995302397d+00]
  properties%photoIonizationRate        = 0.2219143893216550d-13
  properties%photoHeatingRate           = 0.1516777172160954d-31
  properties%photoIonizationRatePrevious= 0.3471374150826434d-12
  properties%photoHeatingRatePrevious   = 0.2365279824300775d-30
  call radiativeTransferMatter_%stateSolve(properties,status)
  call Assert('case #4',status,errorStatusSuccess)
  ! Done with atomic matter.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
#ifdef USEMPI
  call mpiFinalize()
#endif
end program Test_Radiative_Transfer_State_Solver
