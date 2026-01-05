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
Implements a galactic filter which tests whether the given node has a specified label.
!!}
  
  !![
  <galacticFilter name="galacticFilterLabelled">
   <description>Tests whether the given node has a specified label.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterLabelled
     !!{
     Tests whether the given node has a specified label.
     !!}
     private
     type   (varying_string) :: label
     integer                 :: labelID
   contains
     procedure :: passes => labelledPasses
  end type galacticFilterLabelled

  interface galacticFilterLabelled
     !!{
     Constructors for the \refClass{galacticFilterLabelled} galactic filter class.
     !!}
     module procedure labelledConstructorParameters
     module procedure labelledConstructorInternal
  end interface galacticFilterLabelled

contains

  function labelledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterLabelled} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(galacticFilterLabelled)                :: self
    type(inputParameters       ), intent(inout) :: parameters
    type(varying_string        )                :: label

    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>The label upon which to filter.</description>
    </inputParameter>
    !!]
    self=galacticFilterLabelled(label)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function labelledConstructorParameters
  
  function labelledConstructorInternal(label) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterLabelled} galactic filter class.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type(galacticFilterLabelled)                :: self
    type(varying_string        ), intent(inout) :: label
    !![
    <constructorAssign variables="label"/>
    !!]

    self%labelID=nodeLabelRegister(char(label))
    return
  end function labelledConstructorInternal

  logical function labelledPasses(self,node)
    !!{
    Implement a filter on node labels.
    !!}
    use :: Nodes_Labels, only : nodeLabelIsPresent
    implicit none
    class(galacticFilterLabelled), intent(inout)         :: self
    type (treeNode              ), intent(inout), target :: node

    labelledPasses=nodeLabelIsPresent(self%labelID,node)
    return
  end function labelledPasses
