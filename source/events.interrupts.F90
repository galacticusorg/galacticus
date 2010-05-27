!% Contains a module which specifies the template for interrupt procedures.

module Events_Interrupts
  !% Specifies the template for interrupt procedures.
  use Tree_Nodes
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
