!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements stopping criteria objects for use when constraining \glc.

module Constraints_Stopping_Criteria
  !% Implements stopping criteria objects for use when constraining \glc.
  use Constraints_State
  use Constraints_Convergence
  private
  public :: stoppingCriterionNew

  ! Define the basic stoppingCriterion class.
  type, abstract, public :: stoppingCriterion
   contains
     !@ <objectMethods>
     !@   <object>stoppingCriterion</object>
     !@   <objectMethod>
     !@     <method>stop</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if the simulation should stop.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(stoppingCriterionStop), deferred :: stop
  end type stoppingCriterion

  abstract interface
     logical function stoppingCriterionStop(self,simulationState,simulationConvergence)
       import :: stoppingCriterion, state, convergence
       class(stoppingCriterion), intent(in   ) :: self
       class(state            ), intent(inout) :: simulationState
       class(convergence      ), intent(inout) :: simulationConvergence
     end function stoppingCriterionStop
  end interface

  ! Include all stoppingCriterion types.
  include 'constraints.stopping_criteria.never.type.inc'
  include 'constraints.stopping_criteria.step_count.type.inc'
  include 'constraints.stopping_criteria.correlation_length_count.type.inc'

contains

  function stoppingCriterionNew(definition) result (newStoppingCriterion)
    !% Create a new {\tt stoppingCriterion} from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class(stoppingCriterion), pointer                :: newStoppingCriterion
    type (node             ), pointer, intent(in   ) :: definition
    type (node             ), pointer                :: stopAfterCountDefinition
    integer                                          :: stopAfterCount
   
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("never")
       allocate(stoppingCriterionNever :: newStoppingCriterion)
       select type (newStoppingCriterion)
       type is (stoppingCriterionNever)
          newStoppingCriterion=stoppingCriterionNever()
       end select
   case ("stepCount")
       allocate(stoppingCriterionStepCount :: newStoppingCriterion)
       select type (newStoppingCriterion)
       type is (stoppingCriterionStepCount)
          stopAfterCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"stopAfterCount")
          call extractDataContent(stopAfterCountDefinition,stopAfterCount)
          newStoppingCriterion=stoppingCriterionStepCount(stopAfterCount)
       end select
    case ("correlationLengthCount")
       allocate(stoppingCriterionCorrelationLengthCount :: newStoppingCriterion)
       select type (newStoppingCriterion)
       type is (stoppingCriterionCorrelationLengthCount)
          stopAfterCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"stopAfterCount")
          call extractDataContent(stopAfterCountDefinition,stopAfterCount)
          newStoppingCriterion=stoppingCriterionCorrelationLengthCount(stopAfterCount)
       end select
    case default
       call Galacticus_Error_Report('stoppingCriterionNew','stoppingCriterion type is unrecognized')
    end select
    return
  end function stoppingCriterionNew

  ! Include all stopping criterion methods.
  include 'constraints.stopping_criteria.never.methods.inc'
  include 'constraints.stopping_criteria.step_count.methods.inc'
  include 'constraints.stopping_criteria.correlation_length_count.methods.inc'

end module Constraints_Stopping_Criteria
