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
  <task name="taskBuildToolCLASS">
   <description>A task which builds the CLASS tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolCLASS
     !!{
     Implementation of a task which builds the CLASS tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolCLASSPerform
     procedure :: requiresOutputFile => buildToolCLASSRequiresOutputFile
  end type taskBuildToolCLASS

  interface taskBuildToolCLASS
     !!{
     Constructors for the \refClass{taskBuildToolCLASS} task.
     !!}
     module procedure buildToolCLASSParameters
  end interface taskBuildToolCLASS

contains

  function buildToolCLASSParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolCLASS} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolCLASS)                :: self
    type(inputParameters   ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolCLASS()
    return
  end function buildToolCLASSParameters

  subroutine buildToolCLASSPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display         , only : displayIndent             , displayMessage, displayUnindent
    use :: Error, only : errorStatusSuccess
    use :: Interfaces_CLASS, only : Interface_CLASS_Initialize
    implicit none
    class  (taskBuildToolCLASS), intent(inout), target   :: self
    integer                    , intent(  out), optional :: status
    type   (varying_string    )                          :: classPath, classVersion
    !$GLC attributes unused :: self
#include "os.inc"

    call displayIndent  ('Begin task: CLASS tool build')
    call Interface_CLASS_Initialize(                     &
         &                                 classPath   , &
         &                                 classVersion, &
#ifdef __APPLE__
         &                          static=.false.       &
#else
         &                          static=.true.        &
#endif
         &                         )
    call displayMessage('CLASS version '//classVersion//' successfully built in: '//classPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: CLASS tool build')
    return
  end subroutine buildToolCLASSPerform

  logical function buildToolCLASSRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolCLASS), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolCLASSRequiresOutputFile=.false.
    return
  end function buildToolCLASSRequiresOutputFile
