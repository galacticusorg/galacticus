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
  Implementation of the \cite{cole_hierarchical_2000} time available for cooling class.
  !!}

  !![
  <coolingTimeAvailable name="coolingTimeAvailableFormationTime">
   <description>
    A time available for cooling class which implements the algorithm of \cite{cole_hierarchical_2000}, that is, the time
    available is equal to
    \begin{equation}
     t_\mathrm{available} = t - t_\mathrm{form},
    \end{equation}
    where $t_\mathrm{form}$ is the time at which the halo formed (see \S\ref{sec:ComponentFormationTimes}).
   </description>
  </coolingTimeAvailable>
  !!]
  type, extends(coolingTimeAvailableClass) :: coolingTimeAvailableFormationTime
     !!{
     Implementation of a time available for cooling class which implements the algorithm of \cite{cole_hierarchical_2000}.
     !!}
     private
     integer :: nodeFormationTimeID
   contains
     procedure :: timeAvailable             => formationTimeTimeAvailable
     procedure :: timeAvailableIncreaseRate => formationTimeTimeAvailableIncreaseRate
  end type coolingTimeAvailableFormationTime

  interface coolingTimeAvailableFormationTime
     !!{
     Constructors for the \cite{cole_hierarchical_2000} time available for cooling class.
     !!}
     module procedure formationTimeConstructorParameters
     module procedure formationTimeConstructorInternal
  end interface coolingTimeAvailableFormationTime

contains

  function formationTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} time available for cooling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(coolingTimeAvailableFormationTime)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self=coolingTimeAvailableFormationTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function formationTimeConstructorParameters

  function formationTimeConstructorInternal() result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} time available for cooling class.
    !!}
    implicit none
    type(coolingTimeAvailableFormationTime) :: self

    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function formationTimeConstructorInternal
  
  double precision function formationTimeTimeAvailable(self,node)
    !!{
    Returns the time available for cooling (in units of Gyr) in the hot atmosphere for the \cite{cole_hierarchical_2000} model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(coolingTimeAvailableFormationTime), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(nodeComponentBasic               ), pointer       :: basic

    basic                      =>  node %basic                    (                        )
    formationTimeTimeAvailable =  +basic%time                     (                        ) &
         &                        -basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID)
    return
  end function formationTimeTimeAvailable

  double precision function formationTimeTimeAvailableIncreaseRate(self,node)
    !!{
    Compute the rate of increase of the time available for cooling using the \cite{cole_hierarchical_2000} method. We return a rate
    of 1, even though technically it can depend on halo properties.
    !!}
    implicit none
    class(coolingTimeAvailableFormationTime), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Simply return unit rate.
    formationTimeTimeAvailableIncreaseRate=1.0d0
    return
  end function formationTimeTimeAvailableIncreaseRate

