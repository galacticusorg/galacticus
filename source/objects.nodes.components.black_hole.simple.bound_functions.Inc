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

function Node_Component_Black_Hole_Simple_Mass_Distribution(self,weightBy,weightIndex) result(massDistribution_)
  !!{
  Return the mass distribution for the simple black hole component.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeBlackHole   , massTypeBlackHole    , weightByMass, weightByLuminosity, &
       &                                    enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionBlackHole, massDistributionClass
  implicit none
  class  (massDistributionClass       ), pointer                 :: massDistribution_
  class  (nodeComponentBlackHoleSimple), intent(inout)           :: self
  type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy
  integer                              , intent(in   ), optional :: weightIndex
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
     <referenceConstruct object="massDistribution_" constructor="massDistributionBlackHole(mass=self%mass(),componentType=componentTypeBlackHole,massType=massTypeBlackHole)"/>
     !!]
  end select
  return
end function Node_Component_Black_Hole_Simple_Mass_Distribution
