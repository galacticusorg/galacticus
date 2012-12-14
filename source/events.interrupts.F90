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

!% Contains a module which specifies the template for interrupt procedures.

module Events_Interrupts
  !% Specifies the template for interrupt procedures.
  use Galacticus_Nodes
  implicit none
  private
  public :: Interrupt_Procedure_Template

  ! Procedure template for interrupt routines.
  abstract interface
     subroutine Interrupt_Procedure_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
     end subroutine Interrupt_Procedure_Template
  end interface
  
end module Events_Interrupts
