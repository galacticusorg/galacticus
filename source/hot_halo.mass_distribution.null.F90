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
A null implementation of the hot halo mass distribution class.
!!}

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionNull">
   <description>
    A hot halo mass distribution class that assumes no hot halo mass distribution. It is useful, for example, when performing
    dark matter-only calculations.
   </description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionNull
     !!{
     A null implementation of the hot halo mass distribution class.
     !!}
     private
   contains
     procedure :: get => nullGet
  end type hotHaloMassDistributionNull

  interface hotHaloMassDistributionNull
     !!{
     Constructors for the null hot halo mass distribution class.
     !!}
     module procedure nullConstructorParameters
  end interface hotHaloMassDistributionNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the null hot halo mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(hotHaloMassDistributionNull)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=hotHaloMassDistributionNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  function nullGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return a null hot halo mass distribution.
    !!}
    implicit none
    class  (massDistributionClass      ), pointer                 :: massDistribution_
    class  (hotHaloMassDistributionNull), intent(inout)           :: self
    type   (treeNode                   ), intent(inout)           :: node
    type   (enumerationWeightByType    ), intent(in   ), optional :: weightBy
    integer                             , intent(in   ), optional :: weightIndex

    massDistribution_ => null()
    return
  end function nullGet
