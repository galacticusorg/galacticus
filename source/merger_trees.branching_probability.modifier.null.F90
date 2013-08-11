!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of branching probabilties in modified Press-Schechter theory.

module Merger_Tree_Branching_Modifiers_Null
  !% Implements a null modifier of merger tree branching rates.
  implicit none
  private
  public :: Merger_Tree_Branching_Modifiers_Null_Initialize

contains

  !# <treeBranchingModifierMethod>
  !#  <unitName>Merger_Tree_Branching_Modifiers_Null_Initialize</unitName>
  !# </treeBranchingModifierMethod>
  subroutine Merger_Tree_Branching_Modifiers_Null_Initialize(treeBranchingModifierMethod,Merger_Tree_Branching_Modifier_Get)
    !% Initialize the null modifier method for merger tree branching rates.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: treeBranchingModifierMethod
    procedure(double precision), intent(inout), pointer :: Merger_Tree_Branching_Modifier_Get

    if (treeBranchingModifierMethod == 'null') Merger_Tree_Branching_Modifier_Get => Merger_Tree_Branching_Modifier_Null
    return
  end subroutine Merger_Tree_Branching_Modifiers_Null_Initialize

  double precision function Merger_Tree_Branching_Modifier_Null(parentDelta,childSigma,parentSigma)
    !% Returns a null (multiplicative) modifier for merger tree branching rates.
    implicit none
    double precision, intent(in   ) :: childSigma, parentDelta, parentSigma

    Merger_Tree_Branching_Modifier_Null=1.0d0
    return
  end function Merger_Tree_Branching_Modifier_Null

end module Merger_Tree_Branching_Modifiers_Null
