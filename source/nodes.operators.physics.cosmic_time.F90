!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node operator class that simply evolves the cosmic time of a node.
  !!}

  !![
  <nodeOperator name="nodeOperatorCosmicTime">
   <description>
    A node operator class that simply evolves the cosmic time of a node.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCosmicTime
     !!{
     A node operator class that simply evolves the cosmic time of a node.
     !!}
     private
   contains
     procedure :: differentialEvolutionAnalytics      => cosmicTimeDifferentialEvolutionAnalytics
     procedure :: preDeterminedSolveAnalytics         => cosmicTimeSolveAnalytics
     procedure :: differentialEvolutionSolveAnalytics => cosmicTimeSolveAnalytics
     procedure :: nodesMerge                          => cosmicTimeNodesMerge
  end type nodeOperatorCosmicTime
  
  interface nodeOperatorCosmicTime
     !!{
     Constructors for the \refClass{nodeOperatorCosmicTime} node operator class.
     !!}
     module procedure cosmicTimeConstructorParameters
  end interface nodeOperatorCosmicTime
  
contains
  
  function cosmicTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCosmicTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorCosmicTime)                :: self
    type (inputParameters       ), intent(inout) :: parameters

    self=nodeOperatorCosmicTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cosmicTimeConstructorParameters

  subroutine cosmicTimeDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorCosmicTime), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentBasic    ), pointer       :: basic

    basic => node%basic()
    call basic%timeAnalytic()
    return
  end subroutine cosmicTimeDifferentialEvolutionAnalytics

  subroutine cosmicTimeSolveAnalytics(self,node,time)
    !!{
    Set values of analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorCosmicTime), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    double precision                        , intent(in   ) :: time
    class           (nodeComponentBasic    ), pointer       :: basic
    !$GLC attributes unused :: self

    basic => node%basic()
    call basic%timeSet(time)
    return
  end subroutine cosmicTimeSolveAnalytics

  subroutine cosmicTimeNodesMerge(self,node)
    !!{
    Act on a merger between nodes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorCosmicTime), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentBasic    ), pointer       :: basic

    ! Record the time at which the node became a satellite - used for computing halo scales etc.
    basic => node%basic()
    call basic%timeLastIsolatedSet(basic%time())
    return
  end subroutine cosmicTimeNodesMerge
  
