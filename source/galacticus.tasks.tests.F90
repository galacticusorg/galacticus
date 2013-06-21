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

!% Contains a module which performs self tests of \glc.

module Galacticus_Tasks_Tests
!% Performs self tests of \glc.
  implicit none
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
       !# <include directive="galacticusSelfTest" type="functionCall" functionType="void">
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
