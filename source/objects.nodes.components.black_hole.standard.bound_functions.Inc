!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

function Node_Component_Black_Hole_Standard_Mass_Distribution(self,weightBy,weightIndex) result(massDistribution_)
  !!{
  Return the mass distribution for the standard black hole component.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeBlackHole   , massTypeBlackHole    , weightByMass, weightByLuminosity, &
       &                                    enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionBlackHole, massDistributionClass
  implicit none
  class  (massDistributionClass         ), pointer                 :: massDistribution_
  class  (nodeComponentBlackHoleStandard), intent(inout)           :: self
  type   (enumerationWeightByType       ), intent(in   ), optional :: weightBy
  integer                                , intent(in   ), optional :: weightIndex
  !![
  <optionalArgument name="weightBy" defaultsTo="weightByMass" />
  !!]
  !$GLC attributes unused :: weightIndex

  if (weightBy_ /= weightByMass) then
     massDistribution_ => null()
     return
  end if
  ! Create a black hole mass distribution.
  allocate(massDistributionBlackHole :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBlackHole)
     !![
     <referenceConstruct object="massDistribution_" constructor="massDistributionBlackHole(mass=max(0.0d0,self%mass()),componentType=componentTypeBlackHole,massType=massTypeBlackHole)"/>
     !!]
  end select
  return
end function Node_Component_Black_Hole_Standard_Mass_Distribution
