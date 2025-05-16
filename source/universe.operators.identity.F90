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

  !![
  <universeOperator name="universeOperatorIdentity">
   <description>An identity operator on universes.</description>
  </universeOperator>
  !!]
  type, extends(universeOperatorClass) :: universeOperatorIdentity
     !!{
     Implementation of an identity operator on universes.
     !!}
     private
   contains
     procedure :: operate => identityOperate
  end type universeOperatorIdentity

  interface universeOperatorIdentity
     !!{
     Constructors for the \refClass{universeOperatorIdentity} universe operator.
     !!}
     module procedure identityConstructorParameters
  end interface universeOperatorIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{universeOperatorIdentity} universe operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(universeOperatorIdentity)                :: self
    type(inputParameters         ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=universeOperatorIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

   subroutine identityOperate(self,universe_)
     !!{
     Perform an identity operation on a universe.
     !!}
     implicit none
     class(universeOperatorIdentity), intent(inout), target :: self
     type (universe                ), intent(inout)         :: universe_
     !$GLC attributes unused :: self, universe_

     return
   end subroutine identityOperate

