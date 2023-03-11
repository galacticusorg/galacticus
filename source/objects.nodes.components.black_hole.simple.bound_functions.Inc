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
Contains custom functions for the simple black hole component.
!!}

function Node_Component_Black_Hole_Simple_Mass_Distribution(self) result(massDistribution_)
  !!{
  Return the mass distribution for the simple black hole component.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeBlackHole   , massTypeBlackHole
  use :: Mass_Distributions        , only : massDistributionBlackHole, massDistributionClass
  implicit none
  class(massDistributionClass       ), pointer       :: massDistribution_
  class(nodeComponentBlackHoleSimple), intent(inout) :: self

  ! Create a black hole mass distribution.
  allocate(massDistributionBlackHole :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBlackHole)
     !![
     <referenceConstruct object="massDistribution_" constructor="massDistributionBlackHole(mass=self%mass(),componentType=componentTypeBlackHole,massType=massTypeBlackHole)"/>
     !!]
  end select
  return
end function Node_Component_Black_Hole_Simple_Mass_Distribution

double precision function Node_Component_Black_Hole_Simple_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
  !!{
  Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
  !!}
  use :: Galactic_Structure_Options, only : weightByMass         , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentBlackHoleSimple), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  class           (massDistributionClass       ), pointer       :: massDistribution_
  !$GLC attributes unused :: weightIndex

  if (weightBy == weightByMass) then
     massDistribution_ => self             %massDistribution    (                             )
     massEnclosed      =  massDistribution_%massEnclosedBySphere(radius,componentType,massType)
     !![
     <objectDestructor name="massDistribution_"/>
     !!]
  else
     massEnclosed      =  0.0d0
  end if
  return
end function Node_Component_Black_Hole_Simple_Enclosed_Mass
