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

!% Contains a module which handles resetting of calculations before a new or updated node is processed.

module Galacticus_Calculations_Resets
  !% Handles resetting of calculations before a new or updated node is processed.
  implicit none
  private
  public :: Galacticus_Calculations_Reset
  
contains
  
  subroutine Galacticus_Calculations_Reset(thisNode)
    !% Calls any routines required to reset all calculation for a new or updated node.
    use Galacticus_Nodes
    !# <include directive="calculationResetTask" type="moduleUse">
    include 'galacticus.calculation_reset.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    
    !# <include directive="calculationResetTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'galacticus.calculation_reset.tasks.inc'
    !# </include>
    
    return
  end subroutine Galacticus_Calculations_Reset
  
end module Galacticus_Calculations_Resets
