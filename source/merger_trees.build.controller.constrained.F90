!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a merger tree build controller class which builds constrained trees.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerConstrained">
   <description>A merger tree build controller class which builds constrained trees.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerConstrained
     !!{     
     A merger tree build controller class which builds constrained trees.
     !!}
     private
  contains
     procedure :: branchingProbabilityObject => branchingProbabilityObjectConstrained
     procedure :: is_contrained => is_constrained
  end type mergerTreeBuildControllerConstrained

  interface mergerTreeBuildControllerConstrained
     !!{
     Constructors for the ``constrained'' merger tree build controller class.
     !!}
     module procedure constrainedConstructorParameters
  end interface mergerTreeBuildControllerConstrained

contains

  function constrainedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``constrained'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeBuildControllerConstrained)               :: self
    type(inputParameters                    ), intent(inout) :: parameters
  
    self=mergerTreeBuildControllerConstrained()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function constrainedConstructorParameters

  function branchingProbabilityObjectConstrained(self,mergerTreeBranchingProbability_)
    !!{
    Return a pointer to a constrained branching probability object.
    }
    implicit none
    class(mergerTreeBuildControllerConstrained), intent(inout)         :: self
    class(mergerTreeBranchingProbabilityClass),  intent(inout)         :: mergerTreeBranchingProbability_
    !$GLC attributes unused :: self

    ! Do something to associate mergerTreeBranchingProbability_ with a constrained branching probability object?

  logical function is_constrained(self,node,treeWalker_)
    !!{
    Mark a node as "constrained" or "unconstrained."
    !!}
    implicit none
    class(mergerTreeBuildControllerConstrained), intent(inout)         :: self    
    type (treeNode                           ), intent(inout), pointer :: node
    class(mergerTreeWalkerClass              ), intent(inout)          :: treeWalker_
    !$GLC attributes unused :: self

    is_constrained=.true.
    ! Returned "true" if this it the primary node, and "false" otherwise (but do primary and secondary always map to nodeNew1 and nodeNew2)?
    if (is_constrained.and.node%isPrimaryProgenitor())
       is_constrained=.true.
    else
       is_constrained=.false.
    return
  end function is_constrained
