!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Implements a node operator class that applies labels to nodes during tree initialization based on a \refClass{galacticFilterClass}.
!!}

  use :: Galactic_Filters, only : galacticFilterClass

  !![
  <nodeOperator name="nodeOperatorLabel">
   <description>A node operator class that applies labels to nodes during tree initialization based on a \refClass{galacticFilterClass}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorLabel
     !!{
     A node operator class that applies labels to nodes during tree initialization based on a \refClass{galacticFilterClass}.
     !!}
     private
     type   (varying_string     )          :: label
     integer                               :: labelID
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::                       labelDestructor
     procedure :: nodeTreeInitialize => labelTreeInitialize
  end type nodeOperatorLabel

  interface nodeOperatorLabel
     !!{
     Constructors for the {\normalfont \ttfamily label} node operator class.
     !!}
     module procedure labelConstructorParameters
     module procedure labelConstructorInternal
  end interface nodeOperatorLabel

contains

  function labelConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily label} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorLabel  )                :: self
    type (inputParameters    ), intent(inout) :: parameters
    class(galacticFilterClass), pointer       :: galacticFilter_
    type (varying_string     )                :: label

    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>The label to apply.</description>
    </inputParameter>
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]    
    self=nodeOperatorLabel(label,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function labelConstructorParameters

  function labelConstructorInternal(label,galacticFilter_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily label} node operator class which takes a parameter set as input.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type (nodeOperatorLabel  )                        :: self
    type (varying_string     ), intent(inout)         :: label
    class(galacticFilterClass), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="label, *galacticFilter_"/>
    !!]

    self%labelID=nodeLabelRegister(char(label))
    return
  end function labelConstructorInternal
  
  subroutine labelDestructor(self)
    !!{
    Destructor for  the ``label'' galactic filter class.
    !!}[<0;149;44m
    implicit none
    type(nodeOperatorLabel), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine labelDestructor

  subroutine labelTreeInitialize(self,node)
    !!{
    Initialize node branch tip indices.
    !!}
    use :: Nodes_Labels, only : nodeLabelSet
    implicit none
    class(nodeOperatorLabel), intent(inout), target  :: self
    type (treeNode         ), intent(inout), target  :: node

    if (self%galacticFilter_%passes(node)) &
         & call nodeLabelSet(self%labelID,node)
    return
  end subroutine labelTreeInitialize
