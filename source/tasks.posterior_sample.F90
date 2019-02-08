!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  use Posterior_Sampling_Simulation, only : posteriorSampleSimulationClass, posteriorSampleSimulation

  !# <task name="taskPosteriorSample">
  !#  <description>A task which performs sampling from a posterior distribution.</description>
  !# </task>
  type, extends(taskClass) :: taskPosteriorSample
     !% Implementation of a task which performs sampling from a posterior distribution.
     private
     class(posteriorSampleSimulationClass), pointer :: posteriorSampleSimulation_
   contains
     final     ::            posteriorSampleDestructor
     procedure :: perform => posteriorSamplePerform
  end type taskPosteriorSample

  interface taskPosteriorSample
     !% Constructors for the {\normalfont \ttfamily posteriorSample} task.
     module procedure posteriorSampleConstructorParameters
     module procedure posteriorSampleConstructorInternal
  end interface taskPosteriorSample

contains

  function posteriorSampleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily posteriorSample} task class which takes a parameter set as input.
    use Input_Parameters
    use Node_Components , only : Node_Components_Initialize
    use Galacticus_Nodes, only : nodeClassHierarchyInitialize
    implicit none
    type (taskPosteriorSample           )                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(posteriorSampleSimulationClass), pointer       :: posteriorSampleSimulation_

    call nodeClassHierarchyInitialize(parameters)
    call Node_Components_Initialize  (parameters)
    !# <objectBuilder class="posteriorSampleSimulation" name="posteriorSampleSimulation_" source="parameters"/>
    self=taskPosteriorSample(posteriorSampleSimulation_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="posteriorSampleSimulation_"/>
    return
  end function posteriorSampleConstructorParameters

  function posteriorSampleConstructorInternal(posteriorSampleSimulation_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily posteriorSample} task class.
    implicit none
    type (taskPosteriorSample           )                        :: self
    class(posteriorSampleSimulationClass), intent(in   ), target :: posteriorSampleSimulation_
    !# <constructorAssign variables="*posteriorSampleSimulation_"/>

    return
  end function posteriorSampleConstructorInternal

  subroutine posteriorSampleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily posteriorSample} task class.
    use Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskPosteriorSample), intent(inout) :: self

    !# <objectDestructor name="self%posteriorSampleSimulation_"/>
    call Node_Components_Uninitialize()
    return
  end subroutine posteriorSampleDestructor

  subroutine posteriorSamplePerform(self)
    !% Perform the posterior sampling.
    use Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent
    implicit none
    class(taskPosteriorSample), intent(inout) :: self

    call Galacticus_Display_Indent('Begin task: posterior sampling')
    call self%posteriorSampleSimulation_%simulate()
    call Galacticus_Display_Unindent('Done task: posterior sampling' )
    return
  end subroutine posteriorSamplePerform
