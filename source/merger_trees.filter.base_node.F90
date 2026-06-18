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
Implements a merger tree filter which passes if the base node in the tree passes the given galactic filter.
!!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <mergerTreeFilter name="mergerTreeFilterBaseNode" docformat="rst">
   <description>
   A merger tree filter which passes if the base node in the tree passes the given galactic filter.
   </description>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterBaseNode
     !!{RST
     A merger tree filter class which passes if the base node in the tree passes the given galactic filter.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           baseNodeDestructor
     procedure :: passes => baseNodePasses
  end type mergerTreeFilterBaseNode

  interface mergerTreeFilterBaseNode
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeFilterBaseNode` merger tree filter class.
     !!}
     module procedure baseNodeConstructorParameters
     module procedure baseNodeConstructorInternal
  end interface mergerTreeFilterBaseNode

contains
  
  function baseNodeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeFilterBaseNode` merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (mergerTreeFilterBaseNode)                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(galacticFilterClass     ), pointer       :: galacticFilter_
    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=mergerTreeFilterBaseNode(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function baseNodeConstructorParameters

  function baseNodeConstructorInternal(galacticFilter_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerTreeFilterBaseNode` merger tree filter class.
    !!}
    implicit none
    type (mergerTreeFilterBaseNode)                        :: self
    class(galacticFilterClass     ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function baseNodeConstructorInternal

  subroutine baseNodeDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerTreeFilterBaseNode` merger tree filter class.
    !!}
    implicit none
    type(mergerTreeFilterBaseNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine baseNodeDestructor

  logical function baseNodePasses(self,tree) result(passes)
    !!{RST
    Implement a merger tree filter which passes if the base node in the tree passes the given merger tree filter.
    !!}
    implicit none
    class(mergerTreeFilterBaseNode), intent(inout) :: self
    type (mergerTree              ), intent(in   ) :: tree
    
    passes=self%galacticFilter_%passes(tree%nodeBase)
    return
  end function baseNodePasses
