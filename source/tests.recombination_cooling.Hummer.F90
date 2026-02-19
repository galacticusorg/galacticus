!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a code to test the ``Hummer'' recombination cooling class.
!!}

program Test_Recombination_Cooling_Hummer
  !!{
  Test the radiative transfer state solver.
  !!}
  use :: Atomic_Cross_Sections_Ionization_Photo      , only : atomicCrossSectionIonizationPhotoVerner
  use :: Atomic_Ionization_Potentials                , only : atomicIonizationPotentialVerner
  use :: Atomic_Rates_Recombination_Radiative        , only : atomicRecombinationRateRadiativeVerner1996     , recombinationCase1
  use :: Atomic_Rates_Recombination_Radiative_Cooling, only : atomicRecombinationRateRadiativeCoolingComputed, atomicRecombinationRateRadiativeCoolingHummer
  use :: Display                                     , only : displayVerbositySet                            , verbosityLevelStandard
  use :: Error                                       , only : errorStatusSuccess
  use :: Unit_Tests                                  , only : Assert                                         , Unit_Tests_Begin_Group                       , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type(atomicCrossSectionIonizationPhotoVerner        ), pointer :: atomicCrossSectionIonizationPhoto_
  type(atomicIonizationPotentialVerner                ), pointer :: atomicIonizationPotential_
  type(atomicRecombinationRateRadiativeCoolingComputed), pointer :: atomicRecombinationRateRadiativeCoolingComputed_
  type(atomicRecombinationRateRadiativeCoolingHummer  ), pointer :: atomicRecombinationRateRadiativeCoolingHummer_
  type(atomicRecombinationRateRadiativeVerner1996     ), pointer :: atomicRecombinationRateRadiative_
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Construct atomic matter.
  allocate(atomicCrossSectionIonizationPhoto_              )
  allocate(atomicIonizationPotential_                      )
  allocate(atomicRecombinationRateRadiativeCoolingComputed_)
  allocate(atomicRecombinationRateRadiativeCoolingHummer_  )
  allocate(atomicRecombinationRateRadiative_               )
  !![
  <referenceConstruct object="atomicCrossSectionIonizationPhoto_"              >
   <constructor>
    atomicCrossSectionIonizationPhotoVerner         (                                   &amp;
      &amp;                                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicIonizationPotential_"                      >
   <constructor>
    atomicIonizationPotentialVerner                (                                    &amp;
      &amp;                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiative_"               >
   <constructor>
    atomicRecombinationRateRadiativeVerner1996     (                                    &amp;
      &amp;                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeCoolingComputed_">
   <constructor>
    atomicRecombinationRateRadiativeCoolingComputed(                                    &amp;
      &amp;                                         atomicCrossSectionIonizationPhoto_, &amp;
      &amp;                                         atomicIonizationPotential_          &amp;
      &amp;                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="atomicRecombinationRateRadiativeCoolingHummer_"  >
   <constructor>
    atomicRecombinationRateRadiativeCoolingHummer  (                                    &amp;
      &amp;                                         0.0d0                             , &amp;
      &amp;                                         atomicRecombinationRateRadiative_   &amp;
      &amp;                                        )
   </constructor>
  </referenceConstruct>
  !!]
  call Unit_Tests_Begin_Group("Recombination cooling coefficient (Hummer)")
  ! 1²2S level of hydrogen, target values from Table 3.2 of "Astrophysics of Gaseous Nebulae and Active Galactic Nuclei", 2nd
  ! edition, by Osterbrock and Ferland, 2006, University Science Books.
  call Assert("H; 1²S; 2.5 × 10³K",atomicRecombinationRateRadiativeCoolingHummer_%rate(1,1,2.5d3,level=recombinationCase1),3.22d-13,relTol=1.0d-2)
  call Assert("H; 1²S; 5.0 × 10³K",atomicRecombinationRateRadiativeCoolingHummer_%rate(1,1,5.0d3,level=recombinationCase1),2.23d-13,relTol=1.0d-2)
  call Assert("H; 1²S; 1.0 × 10⁴K",atomicRecombinationRateRadiativeCoolingHummer_%rate(1,1,1.0d4,level=recombinationCase1),1.52d-13,relTol=1.0d-2)
  call Assert("H; 1²S; 2.0 × 10⁴K",atomicRecombinationRateRadiativeCoolingHummer_%rate(1,1,2.0d4,level=recombinationCase1),1.00d-13,relTol=1.0d-2)
  call Unit_Tests_End_Group()
  ! Done with unit tests.
  call Unit_Tests_Finish()
end program Test_Recombination_Cooling_Hummer
