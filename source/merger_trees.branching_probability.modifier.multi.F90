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
  Implements a merger tree branching probability rate modifier which chains multiple other modifiers.
  !!}

  type, public :: multiModifierList
     class(mergerTreeBranchingProbabilityModifierClass), pointer :: modifier_ => null()
     type (multiModifierList                          ), pointer :: next      => null()
  end type multiModifierList

  !![
  <mergerTreeBranchingProbabilityModifier name="mergerTreeBranchingProbabilityModifierMulti" docformat="rst">
   <description>
   Chains multiple other merger tree branch probability modifiers by taking their product.
   </description>
   <linkedList type="multiModifierList" variable="modifiers" next="next" object="modifier_" objectType="mergerTreeBranchingProbabilityModifierClass"/>
  </mergerTreeBranchingProbabilityModifier>
  !!]
  type, extends(mergerTreeBranchingProbabilityModifierClass) :: mergerTreeBranchingProbabilityModifierMulti
     !!{RST
     A merger tree branching probability rate modifier which chains multiple other modifiers.
     !!}
     private
     type(multiModifierList), pointer :: modifiers => null()
   contains
     final     ::                 multiDestructor
     procedure :: rateModifier => multiRateModifier
  end type mergerTreeBranchingProbabilityModifierMulti

  interface mergerTreeBranchingProbabilityModifierMulti
     !!{RST
     Constructors for the ``mergerTreeBranchingProbabilityModifierMulti`` merger tree branching probability rate class.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface mergerTreeBranchingProbabilityModifierMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{RST
    A constructor for the ``multi`` merger tree branching probability rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeBranchingProbabilityModifierMulti)                :: self
    type   (inputParameters                            ), intent(inout) :: parameters
    type   (multiModifierList                          ), pointer       :: modifier_
    integer                                                             :: i

    self     %modifiers => null()
    modifier_           => null()
    do i=1,parameters%copiesCount('mergerTreeBranchingProbabilityModifier',zeroIfNotPresent=.true.)
       if (associated(modifier_)) then
          allocate(modifier_%next)
          modifier_ => modifier_%next
       else
          allocate(self%modifiers)
          modifier_ => self%modifiers
       end if
       !![
       <objectBuilder class="mergerTreeBranchingProbabilityModifier" name="modifier_%modifier_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeBranchingProbabilityModifier"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(modifiers) result(self)
    !!{RST
    Default constructor for the ``multi`` merger tree branching probability rate class.
    !!}
    implicit none
    type(mergerTreeBranchingProbabilityModifierMulti)                        :: self
    type(multiModifierList                          ), target, intent(in   ) :: modifiers
    type(multiModifierList                          ), pointer               :: modifier_

    self     %modifiers => modifiers
    modifier_           => modifiers
    do while (associated(modifier_))
       !![
       <referenceCountIncrement owner="modifier_" object="modifier_"/>
       !!]
       modifier_ => modifier_%next
    end do
     return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{RST
    Destructor for the ``mergerTreeBranchingProbabilityModifierMulti`` class.
    !!}
    implicit none
    type(mergerTreeBranchingProbabilityModifierMulti), intent(inout) :: self
    type(multiModifierList                          ), pointer       :: modifier_, modifierNext

    if (associated(self%modifiers)) then
       modifier_ => self%modifiers
       do while (associated(modifier_))
          modifierNext => modifier_%next
          !![
          <objectDestructor name="modifier_%modifier_"/>
          !!]
          deallocate(modifier_)
          modifier_ => modifierNext
       end do
    end if
    return
  end subroutine multiDestructor

  double precision function multiRateModifier(self,nodeParent,massParent,sigmaParent,sigmaChild,timeParent) result(modifier)
    !!{RST
    Returns a modifier for merger tree branching rates by taking the product over other modifiers.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityModifierMulti), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: nodeParent
    double precision                                             , intent(in   ) :: sigmaChild , timeParent, &
         &                                                                          sigmaParent, massParent
    type            (multiModifierList                          ), pointer       :: modifier_

    modifier  = 1.0d0
    modifier_ => self%modifiers
    do while (associated(modifier_))
       modifier  =  +modifier                                                                                  &
            &       *modifier_%modifier_%rateModifier(nodeParent,massParent,sigmaParent,sigmaChild,timeParent)
       modifier_ =>  modifier_%next
    end do
     return
  end function multiRateModifier
