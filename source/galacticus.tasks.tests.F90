!% Contains a module which performs self tests of \glc.

module Galacticus_Tasks_Tests
!% Performs self tests of \glc.
  private
  public :: Galacticus_Task_Test
  
contains

  !# <galacticusTask>
  !#  <unitName>Galacticus_Task_Test</unitName>
  !#  <after>Galacticus_Task_Start</after>
  !#  <before>Galacticus_Task_End</before>
  !# </galacticusTask>
  logical function Galacticus_Task_Test()
    use Unit_Tests
    ! Include modules needed for pre-evolution tasks.
    !# <include directive="galacticusSelfTest" type="moduleUse">
    include 'galacticus.tasks.tests.moduleUse.inc'
    !# </include>
    implicit none
    logical, save :: doneTests=.false.
    
    if (.not.doneTests) then
       doneTests=.true.
       
       ! Begin a unit testing group.
       call Unit_Tests_Begin_Group("Galacticus self tests")

       ! Perform any test required.
       !# <include directive="galacticusSelfTest" type="code" action="subroutine">
       include 'galacticus.tasks.tests.inc'
       !# </include>

       ! End the unit testing group.
       call Unit_Tests_End_Group
       
       ! Finish testing and report results.
       call Unit_Tests_Finish

    end if
    Galacticus_Task_Test=.false.
    return
  end function Galacticus_Task_Test

end module Galacticus_Tasks_Tests
