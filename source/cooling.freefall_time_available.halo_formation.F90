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
  Implementation of the \cite{cole_hierarchical_2000} method for computing the time available for freefall in cooling
  calculations in hot halos.
  !!}

  !![
  <freefallTimeAvailable name="freefallTimeAvailableHaloFormation">
   <description>
    A freefall time available class in which the time available for freefall is equal to
    \begin{equation}
     t_\mathrm{available} = t - t_\mathrm{form},
    \end{equation}
    where $t_\mathrm{form}$ is the time at which the halo formed (see \S\ref{sec:ComponentFormationTimes}).
   </description>
  </freefallTimeAvailable>
  !!]
  type, extends(freefallTimeAvailableClass) :: freefallTimeAvailableHaloFormation
     !!{
     Implementation of freefall time available class in which the time available is determined by the halo formation time.
     !!}
     private
     integer :: nodeFormationTimeID
   contains
     procedure :: timeAvailable             => haloFormationTimeAvailable
     procedure :: timeAvailableIncreaseRate => haloFormationTimeAvailableIncreaseRate
  end type freefallTimeAvailableHaloFormation

  interface freefallTimeAvailableHaloFormation
     !!{
     Constructors for the haloFormation freefall time available class.
     !!}
     module procedure haloFormationConstructorParameters
     module procedure haloFormationConstructorInternal
  end interface freefallTimeAvailableHaloFormation

contains

  function haloFormationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{freefallTimeAvailableHaloFormation} freefall time available class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters  , only : inputParameters
    implicit none
    type(freefallTimeAvailableHaloFormation)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=freefallTimeAvailableHaloFormation()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloFormationConstructorParameters

  function haloFormationConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{freefallTimeAvailableHaloFormation} freefall time available class.
    !!}
    implicit none
    type(freefallTimeAvailableHaloFormation) :: self

    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function haloFormationConstructorInternal
  
  double precision function haloFormationTimeAvailableIncreaseRate(self,node)
    !!{
    Compute the rate of increase of the time available for freefall using the \cite{cole_hierarchical_2000} method. We return a rate
    of 1.
    !!}
    implicit none
    class(freefallTimeAvailableHaloFormation), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Simply return unit rate.
    haloFormationTimeAvailableIncreaseRate=1.0d0
    return
  end function haloFormationTimeAvailableIncreaseRate

  double precision function haloFormationTimeAvailable(self,node)
    !!{
    Compute the time available for freefall using the \cite{cole_hierarchical_2000} method. Specifically, the time available is
    assumed to be the time since the halo formation event.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(freefallTimeAvailableHaloFormation), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    class(nodeComponentBasic                ), pointer       :: basic

    basic                      =>  node %basic                    (                        )
    haloFormationTimeAvailable =  +basic%time                     (                        ) &
         &                        -basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID)
    return
  end function haloFormationTimeAvailable
