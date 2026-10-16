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
  <task name="taskBuildToolCloudy" docformat="rst">
   <description>
   A task which downloads, compiles, and installs the Cloudy photoionization and spectral synthesis code, making it available for computing cooling functions, chemical state tables, and emission line strengths of astrophysical plasmas.
   </description>
  </task>
  !!]
  type, extends(taskBuildTool) :: taskBuildToolCloudy
     !!{RST
     Implementation of a task which builds the Cloudy tool.
     !!}
     private
   contains
     procedure :: perform => buildToolCloudyPerform
  end type taskBuildToolCloudy

  interface taskBuildToolCloudy
     !!{RST
     Constructors for the :galacticus-class:`taskBuildToolCloudy` task.
     !!}
     module procedure buildToolCloudyParameters
  end interface taskBuildToolCloudy

contains

  function buildToolCloudyParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskBuildToolCloudy` task class which takes a parameter set as input.
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
    !!{RST
    Builds the tool.
    !!}
    use :: Interfaces_Cloudy         , only : Interface_Cloudy_Initialize
    use :: Tasks_Build_Tool_Utilities, only : Task_Build_Tool
    implicit none
    class  (taskBuildToolCloudy), intent(inout), target   :: self
    integer                     , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call Task_Build_Tool('Cloudy',Interface_Cloudy_Initialize,status)
    return
  end subroutine buildToolCloudyPerform
