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
  Implements a node operator class that sets the positions of orphaned subhalos.
  !!}

  !![
  <nodeOperator name="nodeOperatorPositionOrphans">
   <description>A node operator class that sets the positions of orphaned subhalos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorPositionDiscrete) :: nodeOperatorPositionOrphans
     !!{
     A node operator class that sets the positions of orphaned subhalos.
     !!}
     private
     class(satelliteOrphanDistributionClass), pointer :: satelliteOrphanDistribution_ => null()
   contains
     procedure :: differentialEvolutionStepFinalState => positionOrphansDifferentialEvolutionStepFinalState
  end type nodeOperatorPositionOrphans
  
  interface nodeOperatorPositionOrphans
     !!{
     Constructors for the \refClass{nodeOperatorPositionOrphans} node operator class.
     !!}
     module procedure positionOrphansConstructorParameters
     module procedure positionOrphansConstructorInternal
  end interface nodeOperatorPositionOrphans
  
contains
  
  function positionOrphansConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorPositionOrphans} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorPositionOrphans     )                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(satelliteOrphanDistributionClass), pointer       :: satelliteOrphanDistribution_

    !![
    <objectBuilder class="satelliteOrphanDistribution" name="satelliteOrphanDistribution_" source="parameters"/>
    !!]
    self=nodeOperatorPositionOrphans(satelliteOrphanDistribution_)
    !![
    <objectDestructor name="satelliteOrphanDistribution_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function positionOrphansConstructorParameters

  function positionOrphansConstructorInternal(satelliteOrphanDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorPositionOrphans} node operator class.
    !!}
    implicit none
    type (nodeOperatorPositionOrphans     )                        :: self
    class(satelliteOrphanDistributionClass), intent(in   ), target :: satelliteOrphanDistribution_
    !![
    <constructorAssign variables="*satelliteOrphanDistribution_"/>
    !!]

    return
  end function positionOrphansConstructorInternal

  subroutine positionOrphansDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorPositionOrphans} node operator class.
    !!}
    implicit none
    type(nodeOperatorPositionOrphans), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteOrphanDistribution_"/>
    !!]
    return
  end subroutine positionOrphansDestructor

  subroutine positionOrphansDifferentialEvolutionStepFinalState(self,node)
    !!{
    Assign a position to orphaned galaxies after each successful differential evolution step.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentBasic
    use :: Histories       , only : history
    implicit none
    class  (nodeOperatorPositionOrphans), intent(inout) :: self
    type   (treeNode                   ), intent(inout) :: node
    class  (nodeComponentBasic         ), pointer       :: basic
    class  (nodeComponentPosition      ), pointer       :: position
    type   (history                    )                :: positionHistory
    logical                                             :: isOrphan
    
    if (node%isSatellite()) then
       position        => node    %position       ()
       positionHistory =  position%positionHistory()
       if (positionHistory%exists()) then
          ! Subhalo has a position history - it is an orphan if this history ends prior to the current time.
          basic    => node           %basic(                          )
          isOrphan =  positionHistory%time (size(positionHistory%time)) < basic%time()
       else
          ! Subhalo has no position history - must be an orphan.
          isOrphan=.true.
       end if
    else
       ! Not a subhalo - can not be an orphan.
       isOrphan=.false.
    end if
    if (isOrphan) then
       position => node%position()
       call position%positionSet(self%satelliteOrphanDistribution_%position(node))
       call position%velocitySet(self%satelliteOrphanDistribution_%velocity(node))
    else
       call self%nodeOperatorPositionDiscrete%differentialEvolutionStepFinalState(node)
    end if
    return
  end subroutine positionOrphansDifferentialEvolutionStepFinalState
