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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements an abstract task class for tasks which build external tools.
  !!}

  !![
  <task name="taskBuildTool" abstract="yes" docformat="rst">
   <description>
   An abstract task class for tasks which download, compile, and install an external tool. This is an abstract class---each concrete class must provide the ``perform`` method which builds its tool.
   </description>
  </task>
  !!]
  type, abstract, extends(taskClass) :: taskBuildTool
     !!{RST
     An abstract task class for tasks which build external tools.
     !!}
     private
   contains
     procedure :: requiresOutputFile => buildToolRequiresOutputFile
  end type taskBuildTool

contains

  logical function buildToolRequiresOutputFile(self)
    !!{RST
    Specifies that this task does not require the main output file.
    !!}
    implicit none
    class(taskBuildTool), intent(inout) :: self
    !$GLC attributes unused :: self

    buildToolRequiresOutputFile=.false.
    return
  end function buildToolRequiresOutputFile
