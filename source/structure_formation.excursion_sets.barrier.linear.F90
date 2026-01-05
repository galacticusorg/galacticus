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
Implements a linear excursion set barrier class.
!!}

  !![
  <excursionSetBarrier name="excursionSetBarrierLinear">
   <description>
    A linear excursion set barrier class. The barrier is given by:
    \begin{equation}
    B(S) = B_0 + B_1 S,
    \end{equation}
    where $B_0=${\normalfont \ttfamily [coefficientConstant]}, and $B_0=${\normalfont \ttfamily [coefficientLinear]}.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierLinear
     !!{
     A linear excursion set barrier class.
     !!}
     private
     double precision :: coefficientConstant, coefficientLinear
    contains
     procedure :: barrier         => linearBarrier
     procedure :: barrierGradient => linearBarrierGradient
  end type excursionSetBarrierLinear

  interface excursionSetBarrierLinear
     !!{
     Constructors for the linear excursion set barrier class.
     !!}
     module procedure linearConstructorParameters
     module procedure linearConstructorInternal
  end interface excursionSetBarrierLinear

contains

  function linearConstructorParameters(parameters) result(self)
    !!{
    Constructor for the linear excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetBarrierLinear)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: coefficientConstant, coefficientLinear

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>coefficientConstant</name>
      <source>parameters</source>
      <variable>coefficientConstant</variable>
      <defaultValue>1.67d0</defaultValue>
      <description>The constant coefficient in the linear excursion set barrier.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientLinear</name>
      <source>parameters</source>
      <variable>coefficientLinear</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>The linear coefficient in the linear excursion set barrier.</description>
    </inputParameter>
    !!]
    self=excursionSetBarrierLinear(coefficientConstant,coefficientLinear)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function linearConstructorParameters

  function linearConstructorInternal(coefficientConstant,coefficientLinear) result(self)
    !!{
    Internal constructor for the linear excursion set class.
    !!}
    implicit none
    type            (excursionSetBarrierLinear)                :: self
    double precision                           , intent(in   ) :: coefficientConstant, coefficientLinear
    !![
    <constructorAssign variables="coefficientConstant, coefficientLinear"/>
    !!]

    return
  end function linearConstructorInternal

  double precision function linearBarrier(self,variance,time,node,rateCompute)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierLinear), intent(inout) :: self
    double precision                           , intent(in   ) :: variance, time
    type            (treeNode                 ), intent(inout) :: node
    logical                                    , intent(in   ) :: rateCompute
    !$GLC attributes unused :: time, rateCompute, node

    linearBarrier=+self%coefficientConstant          &
         &        +self%coefficientLinear  *variance
    return
  end function linearBarrier

  double precision function linearBarrierGradient(self,variance,time,node,rateCompute)
    !!{
    Return the gradient with respect to variance of the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierLinear), intent(inout) :: self
    double precision                           , intent(in   ) :: variance   , time
    type            (treeNode                 ), intent(inout) :: node
    logical                                    , intent(in   ) :: rateCompute
    !$GLC attributes unused :: variance, time, rateCompute, node

    linearBarrierGradient=+self%coefficientLinear
    return
  end function linearBarrierGradient
