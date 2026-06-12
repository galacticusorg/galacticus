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
  Implements a filtered node operator class.
  !!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <nodeOperator name="nodeOperatorFiltered" docformat="rst">
   <description>
   A node operator class that applies a collection of child :galacticus-class:`nodeOperatorClass` objects only to nodes that pass a :galacticus-class:`galacticFilterClass` test, enabling conditional application of physical processes based on node properties.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorMulti) :: nodeOperatorFiltered
     !!{RST
     A filtered node operator class, which applies multiple node operators but only if the node passes the provided filter.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::             filteredDestructor
     procedure :: isActive => filteredIsActive
  end type nodeOperatorFiltered

  interface nodeOperatorFiltered
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorFiltered` node operator class.
     !!}
     module procedure filteredConstructorParameters
     module procedure filteredConstructorInternal
  end interface nodeOperatorFiltered

contains

  function filteredConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorFiltered` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Error           , only : Error_Report
    implicit none
    type(nodeOperatorFiltered)                :: self
    type(inputParameters     ), intent(inout) :: parameters


    if (.not.parameters%isPresent('galacticFilter')) call Error_Report('A [galacticFilter] must be given explicitly for this class'//{introspection:location})
    self%nodeOperatorMulti=nodeOperatorMulti(parameters)
    !![
    <objectBuilder class="galacticFilter" name="self%galacticFilter_" source="parameters"/>
    <inputParametersValidate source="parameters" multiParameters="nodeOperator"/>
    !!]
    return
  end function filteredConstructorParameters

  function filteredConstructorInternal(processes,galacticFilter_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorFiltered` node operator class.
    !!}
    implicit none
    type (nodeOperatorFiltered)                        :: self
    type (multiProcessList    ), intent(in   ), target :: processes
    class(galacticFilterClass ), intent(inout), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]
    
    self%nodeOperatorMulti=nodeOperatorMulti(processes)
    return
  end function filteredConstructorInternal

  subroutine filteredDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorFiltered` node operator class.
    !!}
    implicit none
    type(nodeOperatorFiltered), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine filteredDestructor

  logical function filteredIsActive(self,node) result(isActive)
    !!{RST
    Return true if the operators are active for the given ``node``.
    !!}
    implicit none
    class(nodeOperatorFiltered), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    
    isActive=self%galacticFilter_%passes(node)
    return
  end function filteredIsActive
