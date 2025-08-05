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
  Implements a merger tree constructor class which filters out trees from other constructors.
  !!}

  use :: Merger_Tree_Filters  , only : mergerTreeFilterClass
  use :: Merger_Tree_Operators, only : mergerTreeOperatorClass

  !![
  <mergerTreeConstructor name="mergerTreeConstructorFilter">
   <description>
    A merger tree constructor class filters trees from another constructor. Trees which do not pass the filter are
    dropped. Those that do pass the filter are returned.
   </description>
  </mergerTreeConstructor>
  !!]
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorFilter
     !!{
     A class implementing merger tree construction by filtering out trees from another constructor.
     !!}
     private
     class(mergerTreeConstructorClass), pointer :: mergerTreeConstructor_ => null()
     class(mergerTreeFilterClass     ), pointer :: mergerTreeFilter_      => null()
     class(nodeOperatorClass         ), pointer :: nodeOperator_          => null()
     class(mergerTreeOperatorClass   ), pointer :: mergerTreeOperator_    => null()
   contains
     final     ::              filterDestructor
     procedure :: construct => filterConstruct
  end type mergerTreeConstructorFilter

  interface mergerTreeConstructorFilter
     !!{
     Constructors for the \refClass{mergerTreeConstructorFilter} merger tree constructor class.
     !!}
     module procedure filterConstructorParameters
     module procedure filterConstructorInternal
  end interface mergerTreeConstructorFilter

contains

  function filterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeConstructorFilter} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeConstructorFilter)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(mergerTreeConstructorClass ), pointer       :: mergerTreeConstructor_
    class(mergerTreeFilterClass      ), pointer       :: mergerTreeFilter_
    class(nodeOperatorClass          ), pointer       :: nodeOperator_
    class(mergerTreeOperatorClass    ), pointer       :: mergerTreeOperator_
 
    !![
    <objectBuilder class="mergerTreeConstructor" name="mergerTreeConstructor_" source="parameters"/>
    <objectBuilder class="mergerTreeFilter"      name="mergerTreeFilter_"      source="parameters"/>
    <objectBuilder class="nodeOperator"          name="nodeOperator_"          source="parameters"/>
    <objectBuilder class="mergerTreeOperator"    name="mergerTreeOperator_"    source="parameters"/>
    !!]
    self=mergerTreeConstructorFilter(mergerTreeConstructor_,mergerTreeFilter_,mergerTreeOperator_,nodeOperator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeConstructor_"/>
    <objectDestructor name="mergerTreeFilter_"     />
    <objectDestructor name="nodeOperator_"         />
    <objectDestructor name="mergerTreeOperator_"   />
    !!]
    return
  end function filterConstructorParameters

  function filterConstructorInternal(mergerTreeConstructor_,mergerTreeFilter_,mergerTreeOperator_,nodeOperator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeConstructorFilter} merger tree operator class.
    !!}
    implicit none
    type (mergerTreeConstructorFilter)                        :: self
    class(mergerTreeConstructorClass ), intent(in   ), target :: mergerTreeConstructor_
    class(mergerTreeFilterClass      ), intent(in   ), target :: mergerTreeFilter_
    class(nodeOperatorClass          ), intent(in   ), target :: nodeOperator_
    class(mergerTreeOperatorClass    ), intent(in   ), target :: mergerTreeOperator_
    !![
    <constructorAssign variables="*mergerTreeConstructor_, *mergerTreeFilter_, *mergerTreeOperator_, *nodeOperator_"/>
    !!]

    return
  end function filterConstructorInternal

  subroutine filterDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeConstructorFilter} merger tree constructor class.
    !!}
    implicit none
    type(mergerTreeConstructorFilter), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeConstructor_"/>
    <objectDestructor name="self%mergerTreeFilter_"     />
    <objectDestructor name="self%mergerTreeOperator_"   />
    <objectDestructor name="self%nodeOperator_"         />
    !!]
    return
  end subroutine filterDestructor

  function filterConstruct(self,treeNumber,finished) result(tree)
    !!{
    Construct a merger tree by filtering out trees from another constructor.
    !!}
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    use :: Galacticus_Nodes   , only : treeNode
    implicit none
    type   (mergerTree                 ), pointer       :: tree
    class  (mergerTreeConstructorFilter), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: treeNumber
    logical                             , intent(  out) :: finished
    type   (treeNode                   ), pointer       :: nodeWork
    type   (mergerTreeWalkerAllNodes   )                :: treeWalkerAll

    tree => self%mergerTreeConstructor_%construct(treeNumber,finished)
    if (associated(tree)) then
       ! Apply any tree initialization operators.
       call self%mergerTreeOperator_%operatePreInitialization(tree)
       treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
       do while (treeWalkerAll%next(nodeWork))
          call self%nodeOperator_%nodeTreeInitialize(nodeWork)
       end do
       if (.not.self%mergerTreeFilter_%passes(tree)) then
          call tree%destroy()
          deallocate(tree)
          nullify   (tree)
       end if
    end if
    return
  end function filterConstruct
