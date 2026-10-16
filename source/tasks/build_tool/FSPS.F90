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
  <task name="taskBuildToolFSPS" docformat="rst">
   <description>
   A task which downloads, compiles, and installs the FSPS (Flexible Stellar Population Synthesis) code, making it available for computing stellar population spectral energy distributions and broadband luminosities from star formation histories.
   </description>
  </task>
  !!]
  type, extends(taskBuildTool) :: taskBuildToolFSPS
     !!{RST
     Implementation of a task which builds the FSPS tool.
     !!}
     private
   contains
     procedure :: perform => buildToolFSPSPerform
  end type taskBuildToolFSPS

  interface taskBuildToolFSPS
     !!{RST
     Constructors for the :galacticus-class:`taskBuildToolFSPS` task.
     !!}
     module procedure buildToolFSPSParameters
  end interface taskBuildToolFSPS

contains

  function buildToolFSPSParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskBuildToolFSPS` task class which takes a parameter set as input.
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
    !!{RST
    Builds the tool.
    !!}
    use :: Interfaces_FSPS           , only : Interface_FSPS_Initialize
    use :: Tasks_Build_Tool_Utilities, only : Task_Build_Tool
    implicit none
    class  (taskBuildToolFSPS), intent(inout), target   :: self
    integer                   , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call Task_Build_Tool('FSPS',Interface_FSPS_Initialize,status)
    return
  end subroutine buildToolFSPSPerform
