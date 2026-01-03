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

  !![
  <task name="taskReport">
   <description>A task which reports on version and build information.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskReport
     !!{
     Implementation of a task which reports on version and build information.
     !!}
     private
   contains
     procedure :: perform            => reportPerform
     procedure :: requiresOutputFile => reportRequiresOutputFile
  end type taskReport

  interface taskReport
     !!{
     Constructors for the \refClass{taskReport} task.
     !!}
     module procedure reportParameters
  end interface taskReport

contains

  function reportParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskReport} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskReport     )                :: self
    type(inputParameters), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskReport()
    return
  end function reportParameters

  subroutine reportPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display          , only : displayIndent      , displayMessage, displayUnindent
    use :: Output_Build     , only : Output_Build_String
    use :: Error            , only : errorStatusSuccess
    use :: Output_Versioning, only : Version_String
    implicit none
    class  (taskReport), intent(inout), target   :: self
    integer            , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call displayIndent  ('Begin task: report'                         )
    call displayMessage ('This is Galacticus: '//Version_String     ())
    call displayMessage ('Built with: '        //Output_Build_String())
    call displayUnindent('Done task: report'                          )
    if (present(status)) status=errorStatusSuccess
    return
  end subroutine reportPerform

  logical function reportRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskReport), intent(inout) :: self
    !$GLC attributes unused :: self

    reportRequiresOutputFile=.false.
    return
  end function reportRequiresOutputFile
