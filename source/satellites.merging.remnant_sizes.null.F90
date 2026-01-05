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
  Implements a merger remnant size class which takes no action.
  !!}

  !![
  <mergerRemnantSize name="mergerRemnantSizeNull">
   <description>
    A merger remnant size class which does nothing at all. It is useful, for example, when running \glc\ to study dark matter
    only (i.e. when no galaxy properties are computed).
   </description>
  </mergerRemnantSize>
  !!]
  type, extends(mergerRemnantSizeClass) :: mergerRemnantSizeNull
     !!{
     A merger remnant size class which uses takes no action.
     !!}
     private
   contains
     procedure :: get => nullGet
  end type mergerRemnantSizeNull

  interface mergerRemnantSizeNull
     !!{
     Constructors for the \refClass{mergerRemnantSizeNull} merger remnant size class.
     !!}
     module procedure nullConstructorParameters
  end interface mergerRemnantSizeNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerRemnantSizeNull} merger remnant size class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerRemnantSizeNull)                :: self
    type(inputParameters      ), intent(inout) :: parameters

    self=mergerRemnantSizeNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  subroutine nullGet(self,node,radius,velocityCircular,angularMomentumSpecific)
    !!{
    Do not compute the size of the merger remnant for {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (mergerRemnantSizeNull), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(  out) :: radius                 , velocityCircular, &
         &                                                    angularMomentumSpecific
    !$GLC attributes unused :: self,node,radius,velocityCircular,angularMomentumSpecific

    return
  end subroutine nullGet
