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
  <task name="taskBuildToolCloudy">
   <description>A task which builds the Cloudy tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolCloudy
     !!{
     Implementation of a task which builds the Cloudy tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolCloudyPerform
     procedure :: requiresOutputFile => buildToolCloudyRequiresOutputFile
  end type taskBuildToolCloudy

  interface taskBuildToolCloudy
     !!{
     Constructors for the \refClass{taskBuildToolCloudy} task.
     !!}
     module procedure buildToolCloudyParameters
  end interface taskBuildToolCloudy

contains

  function buildToolCloudyParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolCloudy} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolCloudy)                :: self
    type(inputParameters    ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolCloudy()
    return
  end function buildToolCloudyParameters

  subroutine buildToolCloudyPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display          , only : displayIndent              , displayMessage, displayUnindent
    use :: Error , only : errorStatusSuccess
    use :: Interfaces_Cloudy, only : Interface_Cloudy_Initialize
    implicit none
    class  (taskBuildToolCloudy), intent(inout), target   :: self
    integer                     , intent(  out), optional :: status
    type   (varying_string     )                          :: cloudyPath, cloudyVersion
    !$GLC attributes unused :: self
#include "os.inc"

    call displayIndent  ('Begin task: Cloudy tool build')
    call Interface_Cloudy_Initialize(                      &
         &                                  cloudyPath   , &
         &                                  cloudyVersion, &
#ifdef __APPLE__
         &                           static=.false.        &
#else
         &                           static=.true.         &
#endif
         &                          )
    call displayMessage('Cloudy version '//cloudyVersion//' successfully built in: '//cloudyPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: Cloudy tool build')
    return
  end subroutine buildToolCloudyPerform

  logical function buildToolCloudyRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolCloudy), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolCloudyRequiresOutputFile=.false.
    return
  end function buildToolCloudyRequiresOutputFile
