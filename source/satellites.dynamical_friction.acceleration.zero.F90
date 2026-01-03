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
  Implementation of a zero acceleration satellite dynamical friction class.
  !!}

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionZero">
   <description>A satellite dynamical friction class in which the acceleration is always zero.</description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionZero
     !!{
     Implementation of a satellite dynamical fiction class in which the acceleration is always zero.
     !!}
     private
   contains
     procedure :: acceleration => zeroAcceleration
  end type satelliteDynamicalFrictionZero

  interface satelliteDynamicalFrictionZero
     !!{
     Constructors for the zero satellite dynamical friction class.
     !!}
     module procedure zeroConstructorParameters
  end interface satelliteDynamicalFrictionZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the zero satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteDynamicalFrictionZero)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=satelliteDynamicalFrictionZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  function zeroAcceleration(self,node)
    !!{
    Return a zero acceleration for satellites due to dynamical friction.
    !!}
    implicit none
    double precision                                , dimension(3)          :: zeroAcceleration
    class           (satelliteDynamicalFrictionZero), intent(inout), target :: self
    type            (treeNode                      ), intent(inout)         :: node
    !$GLC attributes unused :: self, node

    zeroAcceleration=[0.0d0,0.0d0,0.0d0]
    return
  end function zeroAcceleration
