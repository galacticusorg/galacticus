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
Contains custom functions for the standard black hole component.
!!}

double precision function Node_Component_Black_Hole_Standard_Seed_Spin(self)
  !!{
  Return the seed spin for standard black holes.
  !!}
  implicit none
  class(nodeComponentBlackHoleStandard), intent(inout) :: self
  !$GLC attributes unused :: Self

  Node_Component_Black_Hole_Standard_Seed_Spin=1.0d-3
  return
end function Node_Component_Black_Hole_Standard_Seed_Spin

double precision function Node_Component_Black_Hole_Standard_Spin(self)
  !!{
  Return the spin of a standard black hole.
  !!}
  implicit none
  class           (nodeComponentBlackHoleStandard), intent(inout) :: self
  double precision                                , parameter     :: spinMaximumBlackHole=0.9999d0 !  Maximum allowed spin (useful to avoid infinities).

  ! Keep the spin in range.
  Node_Component_Black_Hole_Standard_Spin=min(max(self%spinData,0.0d0),spinMaximumBlackHole)
  return
end function Node_Component_Black_Hole_Standard_Spin

double precision function Node_Component_Black_Hole_Standard_Mass_Baryonic(self) result(massBaryonic)
  !!{
  Return the baryonic mass for the standard black hole component.
  !!}
  implicit none
  class(nodeComponentBlackHoleStandard), intent(inout) :: self

  massBaryonic=max(0.0d0,self%mass())
  return
end function Node_Component_Black_Hole_Standard_Mass_Baryonic

function Node_Component_Black_Hole_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
  !!{
  Return the mass distribution for the standard black hole component.
  !!}
  use :: Galactic_Structure_Options             , only : componentTypeBlackHole   , massTypeBlackHole           , weightByMass            , weightByLuminosity         , &
       &                                                 enumerationWeightByType  , enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions                     , only : massDistributionBlackHole, massDistributionClass       , massDistributionMatches_, kinematicsDistributionLocal
  use :: Node_Component_Black_Hole_Standard_Data, only : massDistributionPool
  implicit none
  class  (massDistributionClass         ), pointer                 :: massDistribution_
  type   (kinematicsDistributionLocal   ), pointer                 :: kinematicsDistribution_
  class  (nodeComponentBlackHoleStandard), intent(inout)           :: self
  type   (enumerationComponentTypeType  ), intent(in   ), optional :: componentType
  type   (enumerationMassTypeType       ), intent(in   ), optional :: massType
  type   (enumerationWeightByType       ), intent(in   ), optional :: weightBy
  integer                                , intent(in   ), optional :: weightIndex
  logical                                                          :: reused
  integer                                                          :: i
  !![
  <optionalArgument name="weightBy" defaultsTo="weightByMass" />
  !!]
  !$GLC attributes unused :: weightIndex, componentType, massType

  if     (                                                                                           &
       &   massDistributionMatches_(componentTypeBlackHole,massTypeBlackHole,componentType,massType) &
       &  .and.                                                                                      &
       &   weightBy_ == weightByMass                                                                 &
       &  .and.                                                                                      &
       &   self%mass() > 0.0d0                                                                       &
       & ) then
     ! Acquire a black hole mass distribution from the per-thread pool, creating a new one only if no
     ! re-usable object (reference count of unity) is available.
     call massDistributionPool%acquire(i,reused)
     if (.not.reused) allocate(massDistributionBlackHole :: massDistributionPool%slots(i)%object_)
     select type (massDistribution__ => massDistributionPool%slots(i)%object_)
     type is (massDistributionBlackHole)
        if (reused) then
           ! An existing distribution is available for re-use - update its mass for this node. Its
           ! attached kinematics distribution is unchanged and so is retained.
           call massDistribution__%initialize(mass=self%mass())
        else
           ! No pooled object was available - construct a new mass distribution and attach its
           ! kinematics distribution.
           !![
           <referenceConstruct object="massDistribution__" constructor="massDistributionBlackHole(mass=self%mass(),componentType=componentTypeBlackHole,massType=massTypeBlackHole)"/>
           !!]
           allocate(kinematicsDistribution_)
           !![
           <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
           !!]
           call massDistribution__%setKinematicsDistribution(kinematicsDistribution_)
           !![
           <objectDestructor name="kinematicsDistribution_"/>
           !!]
        end if
        ! Increment the reference count since we are returning the object to a calling function which
        ! will therefore hold a reference to it.
        !![
        <referenceCountIncrement object="massDistribution__"/>
        !!]
        massDistribution_ => massDistribution__
     end select
  else
     massDistribution_ => null()
  end if
  return
end function Node_Component_Black_Hole_Standard_Mass_Distribution
