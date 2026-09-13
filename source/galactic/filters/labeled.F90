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

!!{RST
Implements a galactic filter which tests whether the given node has a specified label.
!!}
  
  !![
  <galacticFilter name="galacticFilterLabeled" docformat="rst">
   <description>
   Tests whether the given node has been assigned the label specified by ``[label]``. This filter passes only nodes that carry the designated label, enabling targeted selection of nodes based on categorical metadata attached during tree construction or post-processing.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterLabeled
     !!{RST
     Tests whether the given node has a specified label.
     !!}
     private
     type   (varying_string) :: label
     integer                 :: labelID
   contains
     procedure :: passes => labeledPasses
  end type galacticFilterLabeled

  interface galacticFilterLabeled
     !!{RST
     Constructors for the :galacticus-class:`galacticFilterLabeled` galactic filter class.
     !!}
     module procedure labeledConstructorParameters
     module procedure labeledConstructorInternal
  end interface galacticFilterLabeled

contains

  function labeledConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticFilterLabeled` galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(galacticFilterLabeled)                :: self
    type(inputParameters      ), intent(inout) :: parameters
    type(varying_string       )                :: label

    !![
    <inputParameter docformat="rst">
      <name>label</name>
      <source>parameters</source>
      <description>
      The label string that a node must carry in order to pass this filter; only nodes assigned this exact label during tree construction or post-processing will be selected.
      </description>
    </inputParameter>
    !!]
    self=galacticFilterLabeled(label)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function labeledConstructorParameters
  
  function labeledConstructorInternal(label) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticFilterLabeled` galactic filter class.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type(galacticFilterLabeled)                :: self
    type(varying_string       ), intent(in   ) :: label
    !![
    <constructorAssign variables="label"/>
    !!]

    self%labelID=nodeLabelRegister(char(label))
    return
  end function labeledConstructorInternal

  logical function labeledPasses(self,node)
    !!{RST
    Implement a filter on node labels.
    !!}
    use :: Nodes_Labels, only : nodeLabelIsPresent
    implicit none
    class(galacticFilterLabeled), intent(inout)         :: self
    type (treeNode             ), intent(inout), target :: node

    labeledPasses=nodeLabelIsPresent(self%labelID,node)
    return
  end function labeledPasses
