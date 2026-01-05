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
  Implements a class for ram pressure stripping which simply returns the virial radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloRamPressureStripping name="hotHaloRamPressureStrippingVirialRadius">
   <description>
    A hot halo ram pressure stripping class which sets the ram pressure stripping radius equal to the virial radius of the
    halo. The effectively results in no ram pressure stripping.
   </description>
  </hotHaloRamPressureStripping>
  !!]
  type, extends(hotHaloRamPressureStrippingClass) :: hotHaloRamPressureStrippingVirialRadius
     !!{
     Implementation of a hot halo ram pressure stripping class which simply returns the virial radius.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                   virialRadiusDestructor
     procedure :: radiusStripped => virialRadiusRadiusStripped
  end type hotHaloRamPressureStrippingVirialRadius

  interface hotHaloRamPressureStrippingVirialRadius
     !!{
     Constructors for the \refClass{hotHaloRamPressureStrippingVirialRadius} hot halo ram pressure stripping class.
     !!}
     module procedure virialRadiusConstructorParameters
     module procedure virialRadiusConstructorInternal
  end interface hotHaloRamPressureStrippingVirialRadius

contains

  function virialRadiusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloRamPressureStrippingVirialRadius} hot halo ram pressure stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureStrippingVirialRadius)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloRamPressureStrippingVirialRadius(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function virialRadiusConstructorParameters

  function virialRadiusConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloRamPressureStrippingVirialRadius} hot halo ram pressure stripping class.
    !!}
    implicit none
    type (hotHaloRamPressureStrippingVirialRadius)                        :: self
    class(darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function virialRadiusConstructorInternal

  subroutine virialRadiusDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloRamPressureStrippingVirialRadius} hot halo ram pressure stripping class.
    !!}
    implicit none
    type(hotHaloRamPressureStrippingVirialRadius), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine virialRadiusDestructor

  double precision function virialRadiusRadiusStripped(self,node)
    !!{
    Return the ram pressure stripping radius which is assumed to be equal to the virial radius.
    !!}
    implicit none
    class(hotHaloRamPressureStrippingVirialRadius), intent(inout), target :: self
    type (treeNode                               ), intent(inout), target :: node

    virialRadiusRadiusStripped=self%darkMatterHaloScale_%radiusVirial(node)
    return
  end function virialRadiusRadiusStripped
