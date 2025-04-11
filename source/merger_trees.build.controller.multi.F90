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
Implements a merger tree build controller class which applies multiple other controllers.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerMulti">
   <description>A merger tree build controller class which applies multiple other controllers.</description>
   <linkedList type="controllerList" variable="controllers" next="next" object="controller_" objectType="mergerTreeBuildControllerClass"/>
  </mergerTreeBuildController>
  !!]

  type, public :: controllerList
     class(mergerTreeBuildControllerClass), pointer :: controller_ => null()
     type (controllerList                ), pointer :: next        => null()
  end type controllerList

  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerMulti
     !!{     
     A merger tree build controller class which applies multiple other controllers.
     !!}
     private
     type (controllerList                ), pointer :: controllers            => null()
     class(mergerTreeBuildControllerClass), pointer :: timeMaximumController_ => null()
  contains
     final     ::                               multiDestructor
     procedure :: control                    => multiControl
     procedure :: timeMinimum                => multiTimeMinimum
     procedure :: timeMaximum                => multiTimeMaximum
     procedure :: controlTimeMaximum         => multiControlTimeMaximum
     procedure :: branchingProbabilityObject => multiBranchingProbabilityObject
     procedure :: nodesInserted              => multiNodesInserted
  end type mergerTreeBuildControllerMulti

  interface mergerTreeBuildControllerMulti
     !!{
     Constructors for the ``multi'' merger tree build controller class.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface mergerTreeBuildControllerMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``multi'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeBuildControllerMulti)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    type   (controllerList                ), pointer       :: controller_
    integer                                                :: i

    self       %controllers => null()
    controller_             => null()
    do i=1,parameters%copiesCount('mergerTreeBuildController',zeroIfNotPresent=.true.)
       if (associated(controller_)) then
          allocate(controller_%next)
          controller_ => controller_%next
       else
          allocate(self%controllers)
          controller_ => self%controllers
       end if
       !![
       <objectBuilder class="mergerTreeBuildController" name="controller_%controller_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeBuildController"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(controllers) result(self)
    !!{
    Internal constructor for the ``multi'' merger tree build controller class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(mergerTreeBuildControllerMulti)                        :: self
    type(controllerList                ), target, intent(in   ) :: controllers
    type(controllerList                ), pointer               :: controller_

    self     %controllers => controllers
    controller_           => controllers
    do while (associated(controller_))
       !![
       <referenceCountIncrement owner="controller_" object="controller_"/>
       !!]
       controller_ => controller_%next
    end do

    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily multi} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerMulti), intent(inout) :: self
    type(controllerList                ), pointer       :: controller_, controllerNext

    if (associated(self%controllers)) then
       controller_ => self%controllers
       do while (associated(controller_))
          controllerNext => controller_%next
          !![
	  <objectDestructor name="controller_%controller_"/>
          !!]
          deallocate(controller_)
          controller_ => controllerNext
       end do
    end if
    return
  end subroutine multiDestructor

  logical function multiControl(self,node,treeWalker_) result(control)
    !!{
    Apply control from multiple controllers to a tree under construction.
    !!}
    implicit none
    class(mergerTreeBuildControllerMulti), intent(inout)           :: self    
    type (treeNode                      ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass         ), intent(inout), optional :: treeWalker_
    type (controllerList                ), pointer                 :: controller_

    control     =  .true.
    controller_ => self%controllers
    do while (associated(controller_))
       control     =  controller_%controller_%control(node,treeWalker_)
       controller_ => controller_%next
       if (.not.control) return
    end do
    return
  end function multiControl

  double precision function multiTimeMinimum(self,node,massBranch,criticalOverdensityBranch) result(timeMinimum)
    !!{
    Return the maximum allowed time for this node.
    !!}
    implicit none
    class           (mergerTreeBuildControllerMulti), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: massBranch  , criticalOverdensityBranch,
    type            (controllerList                ), pointer       :: controller_
    double precision                                                :: timeMinimum_

    timeMinimum =  -huge(0.0d0)
    controller_ => self%controllers
    do while (associated(controller_))
       timeMinimum_ =  controller_%controller_%timeMinimum(node,massBranch,criticalOverdensityBranch)
       timeMinimum  =  max(timeMinimum,timeMinimum_)
       controller_  => controller_%next
    end do
    return
  end function multiTimeMinimum
  
  double precision function multiTimeMaximum(self,node,massBranch,criticalOverdensityBranch,timeReference,insertNode) result(timeMaximum)
    !!{
    Return the maximum allowed time for this node.
    !!}
    implicit none
    class           (mergerTreeBuildControllerMulti), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: massBranch   , criticalOverdensityBranch, &
         &                                                             timeReference
    logical                                         , intent(  out) :: insertNode
    type            (controllerList                ), pointer       :: controller_
    double precision                                                :: timeMaximum_
    logical                                                         :: insertNode_

    timeMaximum =  +huge(0.0d0)
    insertNode  =  .false.
    controller_ => self%controllers
    do while (associated(controller_))
       timeMaximum_ =  controller_%controller_%timeMaximum(node,massBranch,criticalOverdensityBranch,timeReference,insertNode_)
       if (timeMaximum_ < timeMaximum) then
          timeMaximum                 =  timeMaximum_
          insertNode                  =  insertNode_
          self%timeMaximumController_ => controller_%controller_
       end if
       controller_  => controller_%next
    end do
    return
  end function multiTimeMaximum
  
  logical function multiControlTimeMaximum(self,node,massBranch,criticalOverdensityBranch,nodeIndex) result(control)
    !!{
    Control when the maximum time is reached.
    !!}
    implicit none
    class           (mergerTreeBuildControllerMulti), intent(inout)         :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: massBranch, criticalOverdensityBranch
    integer         (kind=kind_int8                ), intent(inout)         :: nodeIndex

    control                     =  self%timeMaximumController_%controlTimeMaximum(node,massBranch,criticalOverdensityBranch,nodeIndex)
    self%timeMaximumController_ => null()
    return
  end function multiControlTimeMaximum

  function multiBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerMulti     ), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    type (controllerList                     ), pointer       :: controller_
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability__

    mergerTreeBranchingProbability_ => null()
    controller_                     => self%controllers
    do while (associated(controller_))
       mergerTreeBranchingProbability__ => controller_%controller_%branchingProbabilityObject(node)
       if (associated(mergerTreeBranchingProbability_)) then
          if (.not.associated(mergerTreeBranchingProbability_,mergerTreeBranchingProbability__)) &
               & call Error_Report('ambiguous branching probability'//{introspection:location})
       else
          mergerTreeBranchingProbability_ => mergerTreeBranchingProbability__
       end if
       controller_ => controller_%next
    end do
    return
  end function multiBranchingProbabilityObject

  subroutine multiNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerMulti), intent(inout)           :: self
    type   (treeNode                      ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                      ), intent(inout), optional :: nodeProgenitor2
    logical                                , intent(in   ), optional :: didBranch

    return
  end subroutine multiNodesInserted
