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
  Implementation of zero mass loss rate from dark matter halos.
  !!}

  !![
  <darkMatterHaloMassLossRate name="darkMatterHaloMassLossRateZero">
   <description>
    A dark matter halo mass loss rate class which assumes a zero rate of mass loss from dark matter halos.
   </description>
  </darkMatterHaloMassLossRate>
  !!]
  type, extends(darkMatterHaloMassLossRateClass) :: darkMatterHaloMassLossRateZero
     !!{
     Implementation of a dark matter halo mass loss rate class which assumes a zero rate of mass loss.
     !!}
     private
   contains
     procedure :: rate => zeroRate
  end type darkMatterHaloMassLossRateZero

  interface darkMatterHaloMassLossRateZero
     !!{
     Constructors for the zero dark matter halo mass loss rate class.
     !!}
     module procedure zeroConstructorParameters
  end interface darkMatterHaloMassLossRateZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the zero dark matter halo mass loss rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterHaloMassLossRateZero)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=darkMatterHaloMassLossRateZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroRate(self,node)
    !!{
    Returns the mass loss rate from the dark matter halo of the given \gls{node} in units of $M_\odot$/Gyr.
    !!}
    implicit none
    class(darkMatterHaloMassLossRateZero), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    !$GLC attributes unused :: self,node

    zeroRate=0.0d0
    return
  end function zeroRate
