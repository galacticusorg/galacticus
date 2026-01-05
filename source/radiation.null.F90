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
  Implements a class for null radiation fields.
  !!}

  !![
  <radiationField name="radiationFieldNull">
   <description>A radiation field class for null fields.</description>
  </radiationField>
  !!]
  type, extends(radiationFieldClass) :: radiationFieldNull
     !!{
     A radiation field class for null fields.
     !!}
     private
     double precision :: time_
   contains
     procedure :: flux              => nullFlux
     procedure :: time              => nullTime
     procedure :: timeSet           => nullTimeSet
     procedure :: timeDependentOnly => nullTimeDependentOnly
  end type radiationFieldNull

  interface radiationFieldNull
     !!{
     Constructors for the \refClass{radiationFieldNull} radiation field class.
     !!}
     module procedure nullConstructorParameters
     module procedure nullConstructorInternal
  end interface radiationFieldNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiationFieldNull} radiation field class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(radiationFieldNull)                :: self
    type(inputParameters   ), intent(inout) :: parameters

    self=radiationFieldNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  function nullConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{radiationFieldNull} radiation field class.
    !!}
    implicit none
    type(radiationFieldNull) :: self

    self%time_=-huge(0.0d0)
    return
  end function nullConstructorInternal

  double precision function nullFlux(self,wavelength,node)
    !!{
    Return the flux of a null radiation field.
    !!}
    implicit none
    class           (radiationFieldNull), intent(inout) :: self
    double precision                    , intent(in   ) :: wavelength
    type            (treeNode          ), intent(inout) :: node
    !$GLC attributes unused :: self, wavelength, node

    nullFlux=0.0d0
    return
  end function nullFlux

  double precision function nullTime(self)
    !!{
    Return the time for which this radiation field is set.
    !!}
    implicit none
    class(radiationFieldNull), intent(inout) :: self

    nullTime=self%time_
    return
  end function nullTime

  subroutine nullTimeSet(self,time)
    !!{
    Set the time for this radiation field.
    !!}
    implicit none
    class           (radiationFieldNull), intent(inout) :: self
    double precision                    , intent(in   ) :: time

    self%time_=time
    return
  end subroutine nullTimeSet

  logical function nullTimeDependentOnly(self)
    !!{
    Return true as this radiation field depends on time only.
    !!}
    implicit none
    class(radiationFieldNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullTimeDependentOnly=.true.
    return
  end function nullTimeDependentOnly
