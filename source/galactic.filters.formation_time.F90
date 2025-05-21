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
Implements a galactic filter which removes recently-formed halos.
!!}

  !![
  <galacticFilter name="galacticFilterFormationTime">
   <description>
   A filter which removes recently-formed halos. Halos with a formation time greater than the current time minus $\Delta
   t=${\normalfont \ttfamily [timeRecent]} are removed.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterFormationTime
     !!{
     A galactic filter which implements a galactic filter which removes recently-formed halos.
     !!}
     private
     integer          :: nodeFormationTimeID
     double precision :: timeRecent
   contains
     procedure :: passes => formationTimePasses
  end type galacticFilterFormationTime

  interface galacticFilterFormationTime
     !!{
     Constructors for the \refClass{galacticFilterFormationTime} galactic filter class.
     !!}
     module procedure formationTimeConstructorParameters
     module procedure formationTimeConstructorInternal
  end interface galacticFilterFormationTime

contains

  function formationTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterFormationTime} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterFormationTime)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: timeRecent

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>timeRecent</name>
      <source>parameters</source>
      <variable>timeRecent</variable>
      <description>The parameter $\Delta t$ (in units of Gyr) appearing in the formation time galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterFormationTime(timeRecent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function formationTimeConstructorParameters

  function formationTimeConstructorInternal(timeRecent) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterFormationTime} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterFormationTime)                :: self
    double precision                             , intent(in   ) :: timeRecent
    !![
    <constructorAssign variables="timeRecent"/>
    !!]

    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function formationTimeConstructorInternal

  logical function formationTimePasses(self,node)
    !!{
    Implement a filter which rejects halos that formed too recently.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticFilterFormationTime), intent(inout)         :: self
    type (treeNode                   ), intent(inout), target :: node
    class(nodeComponentBasic         ), pointer               :: basic

    basic               =>   node %basic                    (                        )
    formationTimePasses =   +basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID) &
         &                 <=                                                          &
         &                  +basic%time                     (                        ) &
         &                  -self %timeRecent
    return
  end function formationTimePasses
