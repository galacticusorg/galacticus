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

!!{RST
Implements a linear excursion set barrier class.
!!}

  !![
  <excursionSetBarrier name="excursionSetBarrierLinear" docformat="rst">
   <description>
   A linear excursion set barrier class. The barrier is given by:

   .. math::

      B(S) = B_0 + B_1 S,

   where :math:`B_0=`\ ``[coefficientConstant]``, and :math:`B_0=`\ ``[coefficientLinear]``.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierLinear
     !!{RST
     A linear excursion set barrier class.
     !!}
     private
     double precision :: coefficientConstant, coefficientLinear
    contains
     procedure :: barrier         => linearBarrier
     procedure :: barrierGradient => linearBarrierGradient
  end type excursionSetBarrierLinear

  interface excursionSetBarrierLinear
     !!{RST
     Constructors for the linear excursion set barrier class.
     !!}
     module procedure linearConstructorParameters
     module procedure linearConstructorInternal
  end interface excursionSetBarrierLinear

contains

  function linearConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the linear excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetBarrierLinear)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: coefficientConstant, coefficientLinear

    ! Check and read parameters.
    !![
    <inputParameter docformat="rst">
      <name>coefficientConstant</name>
      <source>parameters</source>
      <variable>coefficientConstant</variable>
      <defaultValue>1.67d0</defaultValue>
      <description>
      The constant (zero-order) coefficient :math:`B_0` in the linear excursion-set barrier :math:`B(\sigma^2) = B_0 + B_1\,\sigma^2`; corresponds to the spherical collapse threshold :math:`\delta_\mathrm{c} \approx 1.686` in the simplest case.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>coefficientLinear</name>
      <source>parameters</source>
      <variable>coefficientLinear</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The linear (first-order in :math:`\sigma^2`) coefficient :math:`B_1` in the excursion-set barrier :math:`B(\sigma^2) = B_0 + B_1\,\sigma^2`; a non-zero value produces a moving barrier that mimics ellipsoidal collapse corrections to the halo mass function.
      </description>
    </inputParameter>
    !!]
    self=excursionSetBarrierLinear(coefficientConstant,coefficientLinear)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function linearConstructorParameters

  function linearConstructorInternal(coefficientConstant,coefficientLinear) result(self)
    !!{RST
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
    !!{RST
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
    !!{RST
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
