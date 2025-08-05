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
  Implements a class for the timescale of ram pressure stripping of hot halos in which the timescale is equal to the halo
  dynamical timescale.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloRamPressureTimescale name="hotHaloRamPressureTimescaleHaloDynamicalTime">
   <description>
    A hot halo ram pressure timescale class in which the timescale is equal to the halo dynamical time of the associated halo.
   </description>
  </hotHaloRamPressureTimescale>
  !!]
  type, extends(hotHaloRamPressureTimescaleClass) :: hotHaloRamPressureTimescaleHaloDynamicalTime
     !!{
     Implementation of a hot halo ram pressure timescale class in which the timescale is equal to the halo dynamical time.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::              haloDynamicalTimeDestructor
     procedure :: timescale => haloDynamicalTimeTimescale
  end type hotHaloRamPressureTimescaleHaloDynamicalTime

  interface hotHaloRamPressureTimescaleHaloDynamicalTime
     !!{
     Constructors for the \refClass{hotHaloRamPressureTimescaleHaloDynamicalTime} hot halo ram pressure timescale class.
     !!}
     module procedure haloDynamicalTimeConstructorParameters
     module procedure haloDynamicalTimeConstructorInternal
  end interface hotHaloRamPressureTimescaleHaloDynamicalTime

contains

  function haloDynamicalTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloRamPressureTimescaleHaloDynamicalTime} hot halo ram pressure timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureTimescaleHaloDynamicalTime)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                    ), pointer       :: darkMatterHaloScale_


    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloRamPressureTimescaleHaloDynamicalTime(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function haloDynamicalTimeConstructorParameters

  function haloDynamicalTimeConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloRamPressureTimescaleHaloDynamicalTime} hot halo ram pressure timescale class.
    !!}
    implicit none
    type (hotHaloRamPressureTimescaleHaloDynamicalTime)                        :: self
    class(darkMatterHaloScaleClass                    ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function haloDynamicalTimeConstructorInternal

  subroutine haloDynamicalTimeDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloRamPressureTimescaleHaloDynamicalTime} hot halo ram pressure timescale class.
    !!}
    implicit none
    type(hotHaloRamPressureTimescaleHaloDynamicalTime), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine haloDynamicalTimeDestructor

  double precision function haloDynamicalTimeTimescale(self,node)
    !!{
    Return a ram pressure timescale due to the hot halo assuming that it equals the halo dynamical time.
    !!}
    implicit none
    class(hotHaloRamPressureTimescaleHaloDynamicalTime), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node

    haloDynamicalTimeTimescale=self%darkMatterHaloScale_%timescaleDynamical(node)
    return
  end function haloDynamicalTimeTimescale
