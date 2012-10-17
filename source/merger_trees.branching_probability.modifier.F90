!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements modifiers for merger tree branching probabilities.

module Merger_Tree_Branching_Modifiers
  !% Implements modifiers for merger tree branching probabilities.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Branching_Modifier

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: treeBranchingModifierInitialized=.false.

  ! Name of branching method used.
  type(varying_string)                           :: treeBranchingModifierMethod

  ! Pointer to the functions that return branching probability modifiers.
  procedure(Merger_Tree_Branching_Modifier_Template), pointer :: Merger_Tree_Branching_Modifier_Get => null()
  abstract interface
     double precision function Merger_Tree_Branching_Modifier_Template(parentDelta,childSigma,parentSigma)
       double precision, intent(in) :: parentDelta,childSigma,parentSigma
     end function Merger_Tree_Branching_Modifier_Template
  end interface
 
contains

  double precision function Merger_Tree_Branching_Modifier(parentDelta,childSigma,parentSigma)
    !% Return a modifier for merger tree branching probabilities.
    implicit none
    double precision, intent(in) :: parentDelta,childSigma,parentSigma

    ! Initialize if necessary.
    call Tree_Branching_Modifiers_Initialize

    ! Call the function to complete the calculation.
    Merger_Tree_Branching_Modifier=Merger_Tree_Branching_Modifier_Get(parentDelta,childSigma,parentSigma)
    return
  end function Merger_Tree_Branching_Modifier

  subroutine Tree_Branching_Modifiers_Initialize
    !% Initializes the tree branching modifier module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="treeBranchingModifierMethod" type="moduleUse">
    include 'merger_trees.branching_probability.modifier.modules.inc'
    !# </include>
    implicit none
 
    ! Initialize if necessary.
    !$omp critical(Tree_Branching_Modifiers_Initialization) 
    if (.not.treeBranchingModifierInitialized) then
       ! Get the tree branching method parameter.
       !@ <inputParameter>
       !@   <name>treeBranchingModifierMethod</name>
       !@   <defaultValue>null</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing modifiers to merger tree branching probabilitie.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('treeBranchingModifierMethod',treeBranchingModifierMethod,defaultValue='null')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="treeBranchingModifierMethod" type="code" action="subroutine">
       !#  <subroutineArgs>treeBranchingModifierMethod,Merger_Tree_Branching_Modifier_Get</subroutineArgs>
       include 'merger_trees.branching_probability.modifier.inc'
       !# </include>
       if (.not.associated(Merger_Tree_Branching_Modifier_Get))  &
            & call Galacticus_Error_Report('Tree_Branching_Modifiers_Initialize','method '//char(treeBranchingModifierMethod)//' is unrecognized')
       
       treeBranchingModifierInitialized=.true.
    end if
    !$omp end critical(Tree_Branching_Modifiers_Initialization)
    return
  end subroutine Tree_Branching_Modifiers_Initialize
  
end module Merger_Tree_Branching_Modifiers
