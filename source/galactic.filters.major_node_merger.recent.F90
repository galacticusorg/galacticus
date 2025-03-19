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
Implements a galactic low-pass filter for time since the last major node merger.
!!}

  !![
  <galacticFilter name="galacticFilterNodeMajorMergerRecent">
   <description>
   A low-pass filter for time since the last major node merger. Halos with a time of the last major node merger greater than or equal to the current time minus
   $\Delta t=${\normalfont \ttfamily [timeRecent]} are passed.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterNodeMajorMergerRecent
     !!{
     A galactic low-pass filter for time since the last major node merger.
     !!}
     private
     integer          :: nodeMajorMergerTimeID
     double precision :: timeRecent
   contains
     procedure :: passes => nodeMajorMergerRecentPasses
  end type galacticFilterNodeMajorMergerRecent

  interface galacticFilterNodeMajorMergerRecent
     !!{
     Constructors for the ``nodeMajorMergerRecent'' galactic filter class.
     !!}
     module procedure nodeMajorMergerRecentConstructorParameters
     module procedure nodeMajorMergerRecentConstructorInternal
  end interface galacticFilterNodeMajorMergerRecent

contains

  function nodeMajorMergerRecentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``nodeMajorMergerRecent'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterNodeMajorMergerRecent)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: timeRecent

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>timeRecent</name>
      <source>parameters</source>
      <variable>timeRecent</variable>
      <description>The parameter $\Delta t$ (in units of Gyr) appearing in the recent node major merger galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterNodeMajorMergerRecent(timeRecent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeMajorMergerRecentConstructorParameters

  function nodeMajorMergerRecentConstructorInternal(timeRecent) result(self)
    !!{
    Internal constructor for the ``nodeMajorMergerRecent'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterNodeMajorMergerRecent)                :: self
    double precision                                     , intent(in   ) :: timeRecent
    !![
    <constructorAssign variables="timeRecent"/>
    !!]

    !![
    <addMetaProperty component="basic" name="nodeMajorMergerTime" id="self%nodeMajorMergerTimeID"/>
    !!]
    return
  end function nodeMajorMergerRecentConstructorInternal

  logical function nodeMajorMergerRecentPasses(self,node)
    !!{
    Implement a low-pass filter for time since the last major node merger.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(galacticFilterNodeMajorMergerRecent), intent(inout)         :: self
    type (treeNode                           ), intent(inout), target :: node
    class(nodeComponentBasic                 ), pointer               :: basic

    basic                       =>   node %basic                    (                          )
    nodeMajorMergerRecentPasses =   +basic%floatRank0MetaPropertyGet(self%nodeMajorMergerTimeID) &
         &                         >=                                                            &
         &                          +basic%time                     (                          ) &
         &                          -self %timeRecent
    return
  end function nodeMajorMergerRecentPasses
