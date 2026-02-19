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
  <task name="taskBuildToolAxionCAMB">
   <description>A task which builds the AxionCAMB tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolAxionCAMB
     !!{
     Implementation of a task which builds the AxionCAMB tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolAxionCAMBPerform
     procedure :: requiresOutputFile => buildToolAxionCAMBRequiresOutputFile
  end type taskBuildToolAxionCAMB

  interface taskBuildToolAxionCAMB
     !!{
     Constructors for the \refClass{taskBuildToolAxionCAMB} task.
     !!}
     module procedure buildToolAxionCAMBParameters
  end interface taskBuildToolAxionCAMB

contains

  function buildToolAxionCAMBParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolAxionCAMB} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolAxionCAMB)                :: self
    type(inputParameters  ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolAxionCAMB()
    return
  end function buildToolAxionCAMBParameters

  subroutine buildToolAxionCAMBPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display             , only : displayIndent                 , displayMessage, displayUnindent
    use :: Error               , only : errorStatusSuccess
    use :: Interfaces_AxionCAMB, only : Interface_AxionCAMB_Initialize
    implicit none
    class  (taskBuildToolAxionCAMB), intent(inout), target   :: self
    integer                        , intent(  out), optional :: status
    type   (varying_string        )                          :: axionCambPath, axionCambVersion
    !$GLC attributes unused :: self
#include "os.inc"

    call displayIndent  ('Begin task: AxionCAMB tool build')
    call Interface_AxionCAMB_Initialize(                         &
         &                                     axionCambPath   , &
         &                                     axionCambVersion, &
#ifdef __APPLE__
         &                              static=.false.           &
#else
         &                              static=.true.            &
#endif
         &                             )
    call displayMessage('AxionCAMB version '//axionCambVersion//' successfully built in: '//axionCambPath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: AxionCAMB tool build')
    return
  end subroutine buildToolAxionCAMBPerform

  logical function buildToolAxionCAMBRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolAxionCAMB), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolAxionCAMBRequiresOutputFile=.false.
    return
  end function buildToolAxionCAMBRequiresOutputFile
