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

  !!{
  Implements a node operator class that updates positions of nodes in discrete steps.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorPositionDiscrete">
   <description>
    A node operator class that interpolates positions of nodes in discrete steps.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPositionDiscrete
     !!{
     A node operator class that interpolates positions of nodes in discrete steps.
     !!}
     private
   contains
     procedure :: differentialEvolutionStepFinalState => positionDiscreteDifferentialEvolutionStepFinalState
  end type nodeOperatorPositionDiscrete
  
  interface nodeOperatorPositionDiscrete
     !!{
     Constructors for the \refClass{nodeOperatorPositionDiscrete} node operator class.
     !!}
     module procedure positionDiscreteConstructorParameters
  end interface nodeOperatorPositionDiscrete

contains
  
  function positionDiscreteConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorPositionDiscrete} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorPositionDiscrete)                :: self
    type(inputParameters             ), intent(inout) :: parameters
  
    self=nodeOperatorPositionDiscrete()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function positionDiscreteConstructorParameters

  subroutine positionDiscreteDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class(nodeOperatorPositionDiscrete), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentPosition       ), pointer       :: position
    !$GLC attributes unused :: self

    position => node%position()
    call position%positionAnalytic()
    call position%velocityAnalytic()
    return
  end subroutine positionDiscreteDifferentialEvolutionAnalytics

  subroutine positionDiscreteDifferentialEvolutionStepFinalState(self,node)
    !!{
    Compute the discrete position and velocity of the node.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Galacticus_Nodes       , only : nodeComponentPosition, nodeComponentBasic
    use            :: Numerical_Interpolation, only : interpolator
    use            :: Histories              , only : history
    implicit none
    class           (nodeOperatorPositionDiscrete), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (nodeComponentPosition       ), pointer       :: position
    integer                                       , parameter     :: positionBegin=1, positionEnd=3, &
         &                                                           velocityBegin=4, velocityEnd=6
    type            (interpolator                )                :: interpolator_
    type            (history                     )                :: history_
    integer         (c_size_t                    )                :: iTime
    !$GLC attributes unused :: self
    
    ! Check if this node has a position history attached to it.
    position => node    %position       ()
    history_ =  position%positionHistory()
    if (history_%exists()) then
       basic => node%basic()
       if (history_%time(1) <= basic%time()) then
          interpolator_=interpolator        (history_%time                 )
          iTime        =interpolator_%locate(basic   %time(),closest=.true.)
          call position%positionSet(history_%data(iTime,positionBegin:positionEnd))
          call position%velocitySet(history_%data(iTime,velocityBegin:velocityEnd))
       end if
    end if
    return
  end subroutine positionDiscreteDifferentialEvolutionStepFinalState
  
