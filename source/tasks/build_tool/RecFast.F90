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
  <task name="taskBuildToolRecFast" docformat="rst">
   <description>
   A task which downloads, compiles, and installs the RecFast recombination code, making it available for computing the ionization history of the universe during hydrogen and helium recombination for use in Boltzmann solvers and CMB calculations.
   </description>
  </task>
  !!]
  type, extends(taskBuildTool) :: taskBuildToolRecFast
     !!{RST
     Implementation of a task which builds the RecFast tool.
     !!}
     private
   contains
     procedure :: perform => buildToolRecFastPerform
  end type taskBuildToolRecFast

  interface taskBuildToolRecFast
     !!{RST
     Constructors for the :galacticus-class:`taskBuildToolRecFast` task.
     !!}
     module procedure buildToolRecFastParameters
  end interface taskBuildToolRecFast

contains

  function buildToolRecFastParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskBuildToolRecFast` task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildToolRecFast)                :: self
    type(inputParameters     ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=taskBuildToolRecFast()
    return
  end function buildToolRecFastParameters

  subroutine buildToolRecFastPerform(self,status)
    !!{RST
    Builds the tool.
    !!}
    use :: Interfaces_RecFast        , only : Interface_RecFast_Initialize
    use :: Tasks_Build_Tool_Utilities, only : Task_Build_Tool
    implicit none
    class  (taskBuildToolRecFast), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status
    !$GLC attributes unused :: self

    call Task_Build_Tool('RecFast',Interface_RecFast_Initialize,status)
    return
  end subroutine buildToolRecFastPerform
