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
  <task name="taskBuildToolCAMB">
   <description>A task which builds the CAMB tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolCAMB
     !!{
     Implementation of a task which builds the CAMB tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolCAMBPerform
     procedure :: requiresOutputFile => buildToolCAMBRequiresOutputFile
  end type taskBuildToolCAMB

  interface taskBuildToolCAMB
     !!{
     Constructors for the \refClass{taskBuildToolCAMB} task.
     !!}
     module procedure buildToolCAMBParameters
  end interface taskBuildToolCAMB

contains

  function buildToolCAMBParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolCAMB} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolCAMB)                :: self
    type(inputParameters  ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolCAMB()
    return
  end function buildToolCAMBParameters

  subroutine buildToolCAMBPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display         , only : displayIndent            , displayMessage, displayUnindent
    use :: Error           , only : errorStatusSuccess
    use :: Interfaces_CAMB , only : Interface_CAMB_Initialize
    implicit none
    class  (taskBuildToolCAMB), intent(inout), target   :: self
    integer                   , intent(  out), optional :: status
    type   (varying_string   )                          :: cambPath, cambVersion
    !$GLC attributes unused :: self
#include "os.inc"
    
    call displayIndent  ('Begin task: CAMB tool build')
    call Interface_CAMB_Initialize(                    &
         &                                cambPath   , &
         &                                cambVersion, &
#ifdef __APPLE__
         &                         static=.false.      &
#else
         &                         static=.true.       &
#endif
         &                        )
    call displayMessage('CAMB version '//cambVersion//' successfully built in: '//cambPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: CAMB tool build')
    return
  end subroutine buildToolCAMBPerform

  logical function buildToolCAMBRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolCAMB), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolCAMBRequiresOutputFile=.false.
    return
  end function buildToolCAMBRequiresOutputFile
