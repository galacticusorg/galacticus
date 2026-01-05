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
  Implementation of a cooling rate class which modifies another cooling rate by cutting off cooling side branches.
  !!}


  !![
  <coolingRate name="coolingRateMainBranchOnly">
   <description>A cooling rate class which modifies another cooling rate by cutting off cooling in side branches.</description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateMainBranchOnly
     !!{
     Implementation of cooling rate class which modifies another cooling rate by cutting off cooling in side branches.
     !!}
     private
     class(coolingRateClass), pointer :: coolingRate_ => null()
   contains
     final     ::         mainBranchOnlyDestructor
     procedure :: rate => mainBranchOnlyRate
  end type coolingRateMainBranchOnly

  interface coolingRateMainBranchOnly
     !!{
     Constructors for the cut off cooling rate class.
     !!}
     module procedure mainBranchOnlyConstructorParameters
     module procedure mainBranchOnlyConstructorInternal
  end interface coolingRateMainBranchOnly

contains

  function mainBranchOnlyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cut off cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRateMainBranchOnly)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(coolingRateClass         ), pointer       :: coolingRate_

    !![
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=coolingRateMainBranchOnly(coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function mainBranchOnlyConstructorParameters

  function mainBranchOnlyConstructorInternal(coolingRate_) result(self)
    !!{
    Internal constructor for the cut off cooling rate class.
    !!}
    type (coolingRateMainBranchOnly)                        :: self
    class(coolingRateClass         ), intent(in   ), target :: coolingRate_

    !![
    <constructorAssign variables="*coolingRate_"/>
    !!]
    return
  end function mainBranchOnlyConstructorInternal

  subroutine mainBranchOnlyDestructor(self)
    !!{
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(coolingRateMainBranchOnly), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine mainBranchOnlyDestructor

  double precision function mainBranchOnlyRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is cut off
    in side branches.
    !!}
    implicit none
    class(coolingRateMainBranchOnly), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node

    if (node%isOnMainBranch()) then
    else
       mainBranchOnlyRate=self%coolingRate_%rate(node)
       mainBranchOnlyRate=0.0d0
    end if
    return
  end function mainBranchOnlyRate

