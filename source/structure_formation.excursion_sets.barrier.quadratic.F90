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
Implements a quadratic excursion set barrier class.
!!}

  !![
  <excursionSetBarrier name="excursionSetBarrierQuadratic">
   <description>
  !!]

  !![
    A quadratic excursion set barrier class. The barrier is given by:
    \begin{equation}
    B(S) = B_0 + B_1 S + B_2 S^2,
    \end{equation}
    where $B_0=${\normalfont \ttfamily [coefficientConstant]}, $B_0=${\normalfont \ttfamily
    [coefficientLinear]}, and $B_2=${\normalfont \ttfamily [coefficientQuadratic]}.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierQuadratic
     !!{
     A quadratic excursion set barrier class.
     !!}
     private
     double precision :: coefficientConstant, coefficientLinear, &
          &              coefficientQuadratic
    contains
     procedure :: barrier         => quadraticBarrier
     procedure :: barrierGradient => quadraticBarrierGradient
  end type excursionSetBarrierQuadratic

  interface excursionSetBarrierQuadratic
     !!{
     Constructors for the quadratic excursion set barrier class.
     !!}
     module procedure quadraticConstructorParameters
     module procedure quadraticConstructorInternal
  end interface excursionSetBarrierQuadratic

contains

  function quadraticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the quadratic excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetBarrierQuadratic)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: coefficientConstant , coefficientLinear, &
         &                                                           coefficientQuadratic

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>coefficientConstant</name>
      <source>parameters</source>
      <variable>coefficientConstant</variable>
      <defaultValue>1.67d0</defaultValue>
      <description>The constant coefficient in the quadratic excursion set barrier.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientLinear</name>
      <source>parameters</source>
      <variable>coefficientLinear</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>The linear coefficient in the quadratic excursion set barrier.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientQuadratic</name>
      <source>parameters</source>
      <variable>coefficientQuadratic</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>The quadratic coefficient in the quadratic excursion set barrier.</description>
    </inputParameter>
    !!]
    self=excursionSetBarrierQuadratic(coefficientConstant,coefficientLinear,coefficientQuadratic)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function quadraticConstructorParameters

  function quadraticConstructorInternal(coefficientConstant,coefficientLinear,coefficientQuadratic) result(self)
    !!{
    Internal constructor for the quadratic excursion set class.
    !!}
    implicit none
    type            (excursionSetBarrierQuadratic)                :: self
    double precision                              , intent(in   ) :: coefficientConstant , coefficientLinear, &
         &                                                           coefficientQuadratic
    !![
    <constructorAssign variables="coefficientConstant, coefficientLinear, coefficientQuadratic"/>
    !!]

    return
  end function quadraticConstructorInternal

  double precision function quadraticBarrier(self,variance,time,node,rateCompute)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierQuadratic), intent(inout) :: self
    double precision                              , intent(in   ) :: variance   , time
    type            (treeNode                    ), intent(inout) :: node
    logical                                       , intent(in   ) :: rateCompute
    !$GLC attributes unused :: time, rateCompute, node

    quadraticBarrier=+self%coefficientConstant              &
         &           +self%coefficientLinear   *variance    &
         &           +self%coefficientQuadratic*variance**2
    return
  end function quadraticBarrier

  double precision function quadraticBarrierGradient(self,variance,time,node,rateCompute)
    !!{
    Return the gradient with respect to variance of the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierQuadratic), intent(inout) :: self
    double precision                              , intent(in   ) :: variance   , time
    type            (treeNode                    ), intent(inout) :: node
    logical                                       , intent(in   ) :: rateCompute
    !$GLC attributes unused :: variance, time, rateCompute, node

    quadraticBarrierGradient=+      self%coefficientLinear             &
         &                   +2.0d0*self%coefficientQuadratic*variance
    return
  end function quadraticBarrierGradient
