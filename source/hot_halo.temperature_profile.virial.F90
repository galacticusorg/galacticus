!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
An implementation of the hot halo temperature class which uses an isothermal virial temperature.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloTemperatureProfile name="hotHaloTemperatureProfileVirial">
   <description>
    A hot halo temperture profile class which assumes an isothermal halo with a temperature equal to the virial temperature of
    the halo.
   </description>
  </hotHaloTemperatureProfile>
  !!]
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileVirial
     !!{
     An implementation of the hot halo temperature profile class which uses an isothermal virial temperature.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                        virialDestructor
     procedure :: temperature         => virialTemperature
     procedure :: temperatureLogSlope => virialTemperatureLogSlope
  end type hotHaloTemperatureProfileVirial

  interface hotHaloTemperatureProfileVirial
     !!{
     Constructors for the virial hot halo temperature profile class.
     !!}
     module procedure virialConstructorParameters
     module procedure virialConstructorInternal
  end interface hotHaloTemperatureProfileVirial

contains

  function virialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the virial cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloTemperatureProfileVirial)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloTemperatureProfileVirial(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function virialConstructorParameters

  function virialConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the virial cooling rate class.
    !!}
    implicit none
    type (hotHaloTemperatureProfileVirial)                        :: self
    class(darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function virialConstructorInternal

  subroutine virialDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily virial} hot halo temperature profile class.
    !!}
    implicit none
    type(hotHaloTemperatureProfileVirial), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine virialDestructor

  double precision function virialTemperature(self,node,radius)
    !!{
    Return the density in a {\normalfont \ttfamily virial} hot halo mass distribution.
    !!}
    implicit none
    class           (hotHaloTemperatureProfileVirial), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    !$GLC attributes unused :: radius

    virialTemperature=self%darkMatterHaloScale_%temperatureVirial(node)
    return
  end function virialTemperature

  double precision function virialTemperatureLogSlope(self,node,radius)
    !!{
    Return the logarithmic slope of the density profile in a {\normalfont \ttfamily virial} hot halo mass
    distribution.
    !!}
    implicit none
    class           (hotHaloTemperatureProfileVirial), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius

    virialTemperatureLogSlope=0.0d0
    return
  end function virialTemperatureLogSlope

