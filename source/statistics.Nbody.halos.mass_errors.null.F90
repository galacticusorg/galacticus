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
Implements a null N-body dark matter halo mass error class.
!!}

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorNull">
   <description>A null N-body dark matter halo mass error class. Errors are always zero.</description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorNull
     !!{
     A null N-body halo mass error class.
     !!}
     private
    contains
     procedure :: errorFractional => nullErrorFractional
     procedure :: correlation     => nullCorrelation
     procedure :: errorZeroAlways => nullErrorZeroAlways
  end type nbodyHaloMassErrorNull

  interface nbodyHaloMassErrorNull
     !!{
     Constructors for the \refClass{nbodyHaloMassErrorNull} N-body halo mass error class.
     !!}
     module procedure nbodyHaloMassErrorNullParameters
  end interface nbodyHaloMassErrorNull

contains

  function nbodyHaloMassErrorNullParameters(parameters)
    !!{
    Constructor for the \refClass{nbodyHaloMassErrorNull} N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyHaloMassErrorNull)                :: nbodyHaloMassErrorNullParameters
    type(inputParameters       ), intent(inout) :: parameters

    ! Check and read parameters.
    nbodyHaloMassErrorNullParameters=nbodyHaloMassErrorNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nbodyHaloMassErrorNullParameters

  double precision function nullErrorFractional(self,node)
    !!{
    Return the fractional error on the mass of an N-body halo.
    !!}
    implicit none
    class(nbodyHaloMassErrorNull), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullErrorFractional=0.0d0
    return
  end function nullErrorFractional

  double precision function nullCorrelation(self,node1,node2)
    !!{
    Return the correlation of the masses of a pair of N-body halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nbodyHaloMassErrorNull), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node1 , node2
    class(nodeComponentBasic    ), pointer       :: basic1, basic2
    !$GLC attributes unused :: self

    basic1 => node1%basic()
    basic2 => node2%basic()
    if     (                                &
         &   basic1%mass() == basic2%mass() &
         &  .and.                           &
         &   basic1%time() == basic2%time() &
         & ) then
       nullCorrelation=1.0d0
    else
       nullCorrelation=0.0d0
    end if
    return
  end function nullCorrelation

  logical function nullErrorZeroAlways(self)
    !!{
    Return true since errors are always zero in this model.
    !!}
    implicit none
    class(nbodyHaloMassErrorNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullErrorZeroAlways=.true.
    return
  end function nullErrorZeroAlways

