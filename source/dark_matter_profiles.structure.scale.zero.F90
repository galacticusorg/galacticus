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
  An implementation of dark matter halo profile scale radii which returns zero radii---useful when scale radii are not
  relevant.
  !!}

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusZero">
   <description>Dark matter halo scale radii class in which are assumed to be zero.</description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusZero
     !!{
     A dark matter halo profile scale radius class in which are assumed to be zero.
     !!}
     private
   contains
     procedure :: radius => zeroRadius
  end type darkMatterProfileScaleRadiusZero

  interface darkMatterProfileScaleRadiusZero
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusZero} dark matter halo profile scale radius class.
     !!}
     module procedure zeroConstructorParameters
  end interface darkMatterProfileScaleRadiusZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusZero} dark matter halo profile scale radius class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileScaleRadiusZero)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=darkMatterProfileScaleRadiusZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroRadius(self,node)
    !!{
    Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileScaleRadiusZero), intent(inout), target :: self
    type (treeNode                        ), intent(inout), target :: node
    !$GLC attributes unused :: self, node

    zeroRadius=0.0d0
    return
  end function zeroRadius
