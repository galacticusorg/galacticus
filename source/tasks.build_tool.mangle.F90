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
  <task name="taskBuildToolMangle">
   <description>A task which builds the {\normalfont \ttfamily mangle} tool.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildToolMangle
     !!{
     Implementation of a task which builds the {\normalfont \ttfamily mangle} tool.
     !!}
     private
   contains
     procedure :: perform            => buildToolManglePerform
     procedure :: requiresOutputFile => buildToolMangleRequiresOutputFile
  end type taskBuildToolMangle

  interface taskBuildToolMangle
     !!{
     Constructors for the \refClass{taskBuildToolMangle} task.
     !!}
     module procedure buildToolMangleParameters
  end interface taskBuildToolMangle

contains

  function buildToolMangleParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildToolMangle} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolMangle)                :: self
    type(inputParameters    ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolMangle()
    return
  end function buildToolMangleParameters

  subroutine buildToolManglePerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display        , only : displayIndent      , displayMessage, displayUnindent
    use :: Error          , only : errorStatusSuccess
    use :: Geometry_Mangle, only : geometryMangleBuild
    implicit none
    class  (taskBuildToolMangle), intent(inout), target   :: self
    integer                     , intent(  out), optional :: status
    type   (varying_string     )                          :: manglePath, mangleVersion
    !$GLC attributes unused :: self
#include "os.inc"
    
    call displayIndent('Begin task: mangle tool build')
    call geometryMangleBuild(                      &
         &                          manglePath   , &
         &                          mangleVersion, &
#ifdef __APPLE__
         &                   static=.false.        &
#else
         &                   static=.true.         &
#endif
         &                  )
    call displayMessage('mangle version '//mangleVersion//' successfully built in: '//manglePath)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: mangle tool build')
    return
  end subroutine buildToolManglePerform

  logical function buildToolMangleRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildToolMangle), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolMangleRequiresOutputFile=.false.
    return
  end function buildToolMangleRequiresOutputFile
