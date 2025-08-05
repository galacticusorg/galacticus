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
  <task name="taskBuildToolFSPS">
   <description>A task which builds the FSPS tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolFSPS
     !!{
     Implementation of a task which builds the FSPS tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolFSPSPerform
     procedure :: requiresOutputFile => buildToolFSPSRequiresOutputFile
  end type taskBuildToolFSPS

  interface taskBuildToolFSPS
     !!{
     Constructors for the \refClass{taskBuildToolFSPS} task.
     !!}
     module procedure buildToolFSPSParameters
  end interface taskBuildToolFSPS

contains

  function buildToolFSPSParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolFSPS} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolFSPS)                :: self
    type(inputParameters  ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolFSPS()
    return
  end function buildToolFSPSParameters

  subroutine buildToolFSPSPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display         , only : displayIndent            , displayMessage, displayUnindent
    use :: Error, only : errorStatusSuccess
    use :: Interfaces_FSPS , only : Interface_FSPS_Initialize
    implicit none
    class  (taskBuildToolFSPS), intent(inout), target   :: self
    integer                   , intent(  out), optional :: status
    type   (varying_string   )                          :: fspsPath, fspsVersion
    !$GLC attributes unused :: self
#include "os.inc"

    call displayIndent  ('Begin task: FSPS tool build')
    call Interface_FSPS_Initialize(                    &
         &                                fspsPath   , &
         &                                fspsVersion, &
#ifdef __APPLE__
         &                         static=.false.      &
#else
         &                         static=.true.       &
#endif
         &                        )
    call displayMessage('FSPS version '//fspsVersion//' successfully built in: '//fspsPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: FSPS tool build')
    return
  end subroutine buildToolFSPSPerform

  logical function buildToolFSPSRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolFSPS), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolFSPSRequiresOutputFile=.false.
    return
  end function buildToolFSPSRequiresOutputFile
