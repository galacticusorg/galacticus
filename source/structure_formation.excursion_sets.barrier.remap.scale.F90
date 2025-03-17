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
Implements an excursion set barrier class which remaps another class by multiplying by a constant.
!!}

  !![
  <excursionSetBarrier name="excursionSetBarrierRemapScale">
   <description>
    An excursion set barrier class which remaps another class by multiplying by a constant given by {\normalfont \ttfamily
    [factor]}.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierRemapScale
     !!{
     An excursion set barrier class which remaps another class by multiplying by a constant.
     !!}
     private
     class           (excursionSetBarrierClass        ), pointer :: excursionSetBarrier_ => null()
     double precision                                            :: factor
     type            (enumerationExcursionSetRemapType)          :: applyTo
   contains
     final     ::                    remapScaleDestructor
     procedure :: barrier         => remapScaleBarrier
     procedure :: barrierGradient => remapScaleBarrierGradient
  end type excursionSetBarrierRemapScale

  interface excursionSetBarrierRemapScale
     !!{
     Constructors for the remap scale excursion set barrier class.
     !!}
     module procedure remapScaleConstructorParameters
     module procedure remapScaleConstructorInternal
  end interface excursionSetBarrierRemapScale

contains

  function remapScaleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the critical overdensity excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetBarrierRemapScale   )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (excursionSetBarrierClass        ), pointer       :: excursionSetBarrier_
    double precision                                                  :: factor
    type            (varying_string                  )                :: applyTo

    !![
    <inputParameter>
      <name>factor</name>
      <source>parameters</source>
      <variable>factor</variable>
      <defaultValue>1.0d0</defaultValue>
      <description>The factor by which to rescale the excursion set barrier.</description>
    </inputParameter>
    <inputParameter>
      <name>applyTo</name>
      <source>parameters</source>
      <variable>applyTo</variable>
      <defaultValue>var_str('nonRates')</defaultValue>
      <description>Specifies whether rescaling is to be applied to the barrier when used for rate calculation, for other calculations, or both.</description>
    </inputParameter>
    <objectBuilder class="excursionSetBarrier" name="excursionSetBarrier_" source="parameters"/>
    !!]
    self=excursionSetBarrierRemapScale(factor,enumerationExcursionSetRemapEncode(applyTo,includesPrefix=.false.),excursionSetBarrier_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="excursionSetBarrier_"/>
    !!]
    return
  end function remapScaleConstructorParameters

  function remapScaleConstructorInternal(factor,applyTo,excursionSetBarrier_) result(self)
    !!{
    Internal constructor for the critical overdensity excursion set class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (excursionSetBarrierRemapScale   )                        :: self
    class           (excursionSetBarrierClass        ), intent(in   ), target :: excursionSetBarrier_
    double precision                                  , intent(in   )         :: factor
    type            (enumerationExcursionSetRemapType), intent(in   )         :: applyTo
    !![
    <constructorAssign variables="factor, applyTo, *excursionSetBarrier_"/>
    !!]

    if (.not.enumerationExcursionSetRemapIsValid(applyTo)) call Error_Report('applyTo is invalid'//{introspection:location})
    return
  end function remapScaleConstructorInternal

  subroutine remapScaleDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetBarrierRemapScale), intent(inout) :: self

    !![
    <objectDestructor name="self%excursionSetBarrier_"/>
    !!]
    return
  end subroutine remapScaleDestructor

  double precision function remapScaleBarrier(self,variance,time,node,rateCompute)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierRemapScale), intent(inout) :: self
    double precision                               , intent(in   ) :: variance   , time
    type            (treeNode                     ), intent(inout) :: node
    logical                                        , intent(in   ) :: rateCompute

    remapScaleBarrier=self%excursionSetBarrier_%barrier(variance,time,node,rateCompute)
    if     (                                                                    &
         &    self%applyTo == excursionSetRemapBoth                             &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapRates    .and.      rateCompute) &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapNonRates .and. .not.rateCompute) &
         & )                                                                    &
         & remapScaleBarrier=+remapScaleBarrier                                 &
         &                   *self%factor
    return
  end function remapScaleBarrier

  double precision function remapScaleBarrierGradient(self,variance,time,node,rateCompute)
    !!{
    Return the gradient with respect to variance of the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierRemapScale), intent(inout) :: self
    double precision                               , intent(in   ) :: variance   , time
    type            (treeNode                     ), intent(inout) :: node
    logical                                        , intent(in   ) :: rateCompute

    remapScaleBarrierGradient=self%excursionSetBarrier_%barrierGradient(variance,time,node,rateCompute)
    if     (                                                                    &
         &    self%applyTo == excursionSetRemapBoth                             &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapRates    .and.      rateCompute) &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapNonRates .and. .not.rateCompute) &
         & )                                                                    &
         & remapScaleBarrierGradient=+remapScaleBarrierGradient                 &
         &                           *self%factor
    return
  end function remapScaleBarrierGradient
