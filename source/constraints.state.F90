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

!% Contains a module which implements state objects for use when constraining \glc.

module Constraints_State
  !% Implements state objects for use when constraining \glc.
  private
  public :: stateNew

  ! Define the basic state class.
  type, abstract, public :: state
     integer :: parameterCount, stepCount
   contains
     !@ <objectMethods>
     !@   <object>state</object>
     !@   <objectMethod>
     !@     <method>count</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return a count of the number of logged state steps.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dimension</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the dimension of the state.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>get</method>
     !@     <type>\doubleone</type>
     !@     <arguments></arguments>
     !@     <description>Return the current state.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>update</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ stateNew\argin, \logicalzero\ logState\argin, \logicalzero isConverged\argin, \logicalone outlierMask\argin</arguments>
     !@     <description>Update the current state to the specified new state. Only log this state (i.e. increase step count etc.) if the {\tt logState} argument is true.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>mean</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return the mean of the state over logged steps.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>variance</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return the variance of the state over logged steps.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>acceptanceRate</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return the acceptance rate of the state over logged steps.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset the state.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                :: count          => stateCount
     procedure                                :: dimension      => stateDimension
     procedure                                :: reset          => stateReset
     procedure(stateGet           ), deferred :: get
     procedure(stateUpdate        ), deferred :: update
     procedure(stateMean          ), deferred :: mean
     procedure(stateVariance      ), deferred :: variance
     procedure(stateAcceptanceRate), deferred :: acceptanceRate
  end type state

  abstract interface
     function stateGet(self)
       import :: state
       class           (state), intent(in   )                  :: self
       double precision       , dimension(self%parameterCount) :: stateGet
     end function stateGet
  end interface

  abstract interface
     subroutine stateUpdate(self,stateNew,logState,isConverged,outlierMask)
       import :: state
       class           (state), intent(inout)                         :: self
       double precision       , intent(in   ), dimension(:)           :: stateNew
       logical                , intent(in   )                         :: logState
       logical                , intent(in   )                         :: isConverged
       logical                , intent(in   ), dimension(:), optional :: outlierMask
     end subroutine stateUpdate
  end interface

  abstract interface
     function stateMean(self)
       import :: state
       class           (state), intent(in   )                  :: self
       double precision       , dimension(self%parameterCount) :: stateMean
     end function stateMean
  end interface

  abstract interface
     function stateVariance(self)
       import :: state
       class           (state), intent(in   )                  :: self
       double precision       , dimension(self%parameterCount) :: stateVariance
     end function stateVariance
  end interface

  abstract interface
     double precisionfunction stateAcceptanceRate(self)
       import :: state
       class(state), intent(in   ) :: self
     end function stateAcceptanceRate
  end interface

  ! Include all state types.
  include 'constraints.state.simple.type.inc'
  include 'constraints.state.history.type.inc'
  include 'constraints.state.correlation.type.inc'

contains

  function stateNew(definition,parameterCount) result (newState)
    !% Create a new state from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class  (state), pointer                :: newState
    type   (node ), pointer, intent(in   ) :: definition
    integer                , intent(in   ) :: parameterCount
    type   (node ), pointer                :: stateAcceptedCountDefinition
    integer                                :: stateAcceptedCount
   
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("simple")
       allocate(stateSimple :: newState)
       select type (newState)
       type is (stateSimple)
          stateAcceptedCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptedStateCount")
          call extractDataContent(stateAcceptedCountDefinition,stateAcceptedCount)
          newState=stateSimple(parameterCount,stateAcceptedCount)
       end select
    case ("history")
       allocate(stateHistory :: newState)
       select type (newState)
       type is (stateHistory)
          stateAcceptedCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptedStateCount")
          call extractDataContent(stateAcceptedCountDefinition,stateAcceptedCount)
          newState=stateHistory(parameterCount,stateAcceptedCount)
       end select
    case ("correlation")
       allocate(stateCorrelation :: newState)
       select type (newState)
       type is (stateCorrelation)
          stateAcceptedCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptedStateCount")
          call extractDataContent(stateAcceptedCountDefinition,stateAcceptedCount)
          newState=stateCorrelation(parameterCount,stateAcceptedCount)
       end select
    case default
       call Galacticus_Error_Report('stateNew','state type is unrecognized')
    end select
    return
  end function stateNew

  integer function stateCount(self)
    !% Returns the number of steps in the current state.
    implicit none
    class(state), intent(in   ) :: self

    stateCount=self%stepCount
    return
  end function stateCount

  integer function stateDimension(self)
    !% Returns the dimension of the state.
    implicit none
    class(state), intent(in   ) :: self

    stateDimension=self%parameterCount
    return
  end function stateDimension

  subroutine stateReset(self)
    !% Reset the state object.
    implicit none
    class(state), intent(inout) :: self

    self%stepCount=0
    return
  end subroutine stateReset

  ! Include all state methods.
  include 'constraints.state.simple.methods.inc'
  include 'constraints.state.history.methods.inc'
  include 'constraints.state.correlation.methods.inc'

end module Constraints_State
