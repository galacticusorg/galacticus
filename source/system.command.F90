!% Contains a module which executes system commands.

module System_Command
  !% Executes system commands.
  private
  public :: System_Command_Do
  
contains
  
  subroutine System_Command_Do(command,iStatus)
    !% Executes the system command {\tt command}, optionally returning the resulting status in {\tt iStatus}.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string), intent(in)            :: command
    integer,              intent(out), optional :: iStatus
    integer                                     :: iStatusActual

    call System(char(command),iStatusActual)
    if (present(iStatus)) then
       iStatus=iStatusActual
    else
       if (iStatusActual.ne.0) call Galacticus_Error_Report('System_Command_Do','failed to execute system command')
    end if
    return
  end subroutine System_Command_Do
  
end module System_Command
