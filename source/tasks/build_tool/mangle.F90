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
  <task name="taskBuildToolMangle" docformat="rst">
   <description>
   A task which downloads, compiles, and installs the ``mangle`` angular mask software, making it available for defining and manipulating survey angular geometries used in galaxy clustering and correlation function analyses.
   </description>
  </task>
  !!]
  type, extends(taskBuildTool) :: taskBuildToolMangle
     !!{RST
     Implementation of a task which builds the ``mangle`` tool.
     !!}
     private
   contains
     procedure :: perform => buildToolManglePerform
  end type taskBuildToolMangle

  interface taskBuildToolMangle
     !!{RST
     Constructors for the :galacticus-class:`taskBuildToolMangle` task.
     !!}
     module procedure buildToolMangleParameters
  end interface taskBuildToolMangle

contains

  function buildToolMangleParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskBuildToolMangle` task class which takes a parameter set as input.
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
    !!{RST
    Builds the tool.
    !!}
    use :: Geometry_Mangle           , only : geometryMangleBuild
    use :: Tasks_Build_Tool_Utilities, only : Task_Build_Tool
    implicit none
    class  (taskBuildToolMangle), intent(inout), target   :: self
    integer                     , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call Task_Build_Tool('mangle',geometryMangleBuild,status)
    return
  end subroutine buildToolManglePerform
