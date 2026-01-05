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
Implements a cold mode hot halo mass distribution core radius class which sets the core radius to a fraction of the virial
radius.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloColdModeCoreRadii name="hotHaloColdModeCoreRadiiVirialFraction">
   <description>Provides an implementation of the cold mode hot halo mass distribution core radius class which sets the core radius to a fraction of the virial radius.</description>
  </hotHaloColdModeCoreRadii>
  !!]
  type, extends(hotHaloColdModeCoreRadiiClass) :: hotHaloColdModeCoreRadiiVirialFraction
     !!{
     An implementation of the cold mode hot halo mass distribution core radius class which sets the core radius to a fraction of the virial radius.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: coreRadiusOverVirialRadius
   contains
     final     ::           virialRadiusFractionDestructor
     procedure :: radius => virialRadiusFractionRadius
  end type hotHaloColdModeCoreRadiiVirialFraction

  interface hotHaloColdModeCoreRadiiVirialFraction
     !!{
     Constructors for the \refClass{hotHaloColdModeCoreRadiiVirialFraction} hot halo mass distribution core radius class.
     !!}
     module procedure virialRadiusFractionConstructorParameters
     module procedure virialRadiusFractionConstructorInternal
  end interface hotHaloColdModeCoreRadiiVirialFraction

contains

  function virialRadiusFractionConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily virialRadiusFraction} cold mode hot halo mass distribution core radius class which builds the
    object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloColdModeCoreRadiiVirialFraction)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    double precision                                                        :: coreRadiusOverVirialRadius

    !![
    <inputParameter>
      <name>coreRadiusOverVirialRadius</name>
      <defaultValue>0.3d0</defaultValue>
      <description>The core radius in the hot halo density profile in units of the virial radius.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloColdModeCoreRadiiVirialFraction(coreRadiusOverVirialRadius,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function virialRadiusFractionConstructorParameters

  function virialRadiusFractionConstructorInternal(coreRadiusOverVirialRadius,darkMatterHaloScale_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily virialRadiusFraction} cold mode hot halo mass distribution core radius class.
    !!}
    implicit none
    type            (hotHaloColdModeCoreRadiiVirialFraction)                        :: self
    double precision                                        , intent(in   )         :: coreRadiusOverVirialRadius
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="coreRadiusOverVirialRadius, *darkMatterHaloScale_"/>
    !!]

    return
  end function virialRadiusFractionConstructorInternal

  subroutine virialRadiusFractionDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloColdModeCoreRadiiVirialFraction} cold mode hot halo mass distribution core radius class.
    !!}
    implicit none
    type(hotHaloColdModeCoreRadiiVirialFraction), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine virialRadiusFractionDestructor

  double precision function virialRadiusFractionRadius(self,node)
    !!{
    Return the core radius of the hot halo mass distribution.
    !!}
    implicit none
    class(hotHaloColdModeCoreRadiiVirialFraction), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node

    virialRadiusFractionRadius=+self                     %coreRadiusOverVirialRadius       &
         &                     *self%darkMatterHaloScale_%radiusVirial              (node)
    return
  end function virialRadiusFractionRadius
