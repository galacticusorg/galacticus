!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  <task name="taskBuildToolRecFast">
   <description>A task which builds the RecFast tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolRecFast
     !!{
     Implementation of a task which builds the RecFast tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolRecFastPerform
     procedure :: requiresOutputFile => buildToolRecFastRequiresOutputFile
  end type taskBuildToolRecFast

  interface taskBuildToolRecFast
     !!{
     Constructors for the \refClass{taskBuildToolRecFast} task.
     !!}
     module procedure buildToolRecFastParameters
  end interface taskBuildToolRecFast

contains

  function buildToolRecFastParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolRecFast} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolRecFast)                :: self
    type(inputParameters     ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolRecFast()
    return
  end function buildToolRecFastParameters

  subroutine buildToolRecFastPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display           , only : displayIndent               , displayMessage, displayUnindent
    use :: Error  , only : errorStatusSuccess
    use :: Interfaces_RecFast, only : Interface_RecFast_Initialize
    implicit none
    class  (taskBuildToolRecFast), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status
    type   (varying_string      )                          :: recfastPath, recfastVersion
    !$GLC attributes unused :: self
#include "os.inc"

    call displayIndent ('Begin task: RecFast tool build')
    call Interface_RecFast_Initialize(                       &
         &                                   recfastPath   , &
         &                                   recfastVersion, &
#ifdef __APPLE__
         &                            static=.false.         &
#else
         &                            static=.true.          &
#endif
         &                           )
    call displayMessage('RecFast version '//recfastVersion//' successfully built in: '//recfastPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: RecFast tool build')
    return
  end subroutine buildToolRecFastPerform

  logical function buildToolRecFastRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolRecFast), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolRecFastRequiresOutputFile=.false.
    return
  end function buildToolRecFastRequiresOutputFile
