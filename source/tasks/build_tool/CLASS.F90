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
  <task name="taskBuildToolCLASS" docformat="rst">
   <description>
   A task which downloads, compiles, and installs the CLASS (Cosmic Linear Anisotropy Solving System) Boltzmann code, making it available as an alternative to CAMB for computing CMB anisotropies, matter transfer functions, and linear power spectra.
   </description>
  </task>
  !!]
  type, extends(taskBuildTool) :: taskBuildToolCLASS
     !!{RST
     Implementation of a task which builds the CLASS tool.
     !!}
     private
   contains
     procedure :: perform => buildToolCLASSPerform
  end type taskBuildToolCLASS

  interface taskBuildToolCLASS
     !!{RST
     Constructors for the :galacticus-class:`taskBuildToolCLASS` task.
     !!}
     module procedure buildToolCLASSParameters
  end interface taskBuildToolCLASS

contains

  function buildToolCLASSParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskBuildToolCLASS` task class which takes a parameter set as input.
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
    !!{RST
    Builds the tool.
    !!}
    use :: Interfaces_CLASS          , only : Interface_CLASS_Initialize
    use :: Tasks_Build_Tool_Utilities, only : Task_Build_Tool
    implicit none
    class  (taskBuildToolCLASS), intent(inout), target   :: self
    integer                    , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call Task_Build_Tool('CLASS',Interface_CLASS_Initialize,status)
    return
  end subroutine buildToolCLASSPerform
