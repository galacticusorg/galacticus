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
An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.
!!}

  !![
  <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationZero">
   <description>An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.</description>
  </hotHaloOutflowReincorporation>
  !!]
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationZero
     !!{
     An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.
     !!}
     private
   contains
     procedure :: rate => zeroRate
  end type hotHaloOutflowReincorporationZero

  interface hotHaloOutflowReincorporationZero
     !!{
     Constructors for the \refClass{hotHaloOutflowReincorporationZero} hot halo outflow reincorporation class.
     !!}
     module procedure zeroConstructorParameters
  end interface hotHaloOutflowReincorporationZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily zero} hot halo outflow reincorporation class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(hotHaloOutflowReincorporationZero)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self=hotHaloOutflowReincorporationZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroRate(self,node)
    !!{
    Return the rate of mass reincorporation for outflowed gas in the hot halo.
    !!}
    implicit none
    class(hotHaloOutflowReincorporationZero), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    zeroRate=0.0d0
    return
  end function zeroRate
