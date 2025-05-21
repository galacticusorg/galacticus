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
Implements a merger tree filter which passes if the base node in the tree passes the given galactic filter.
!!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <mergerTreeFilter name="mergerTreeFilterBaseNode">
   <description>A merger tree filter which passes if the base node in the tree passes the given galactic filter.</description>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterBaseNode
     !!{
     A merger tree filter class which passes if the base node in the tree passes the given galactic filter.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           baseNodeDestructor
     procedure :: passes => baseNodePasses
  end type mergerTreeFilterBaseNode

  interface mergerTreeFilterBaseNode
     !!{
     Constructors for the \refClass{mergerTreeFilterBaseNode} merger tree filter class.
     !!}
     module procedure baseNodeConstructorParameters
     module procedure baseNodeConstructorInternal
  end interface mergerTreeFilterBaseNode

contains
  
  function baseNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeFilterBaseNode} merger tree filter class which takes a parameter set as input.
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
    !!{
    Internal constructor for the \refClass{mergerTreeFilterBaseNode} merger tree filter class.
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
    !!{
    Destructor for the \refClass{mergerTreeFilterBaseNode} merger tree filter class.
    !!}
    implicit none
    type(mergerTreeFilterBaseNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine baseNodeDestructor

  logical function baseNodePasses(self,tree) result(passes)
    !!{
    Implement a merger tree filter which passes if the base node in the tree passes the given merger tree filter.
    !!}
    implicit none
    class(mergerTreeFilterBaseNode), intent(inout) :: self
    type (mergerTree              ), intent(in   ) :: tree
    
    passes=self%galacticFilter_%passes(tree%nodeBase)
    return
  end function baseNodePasses
