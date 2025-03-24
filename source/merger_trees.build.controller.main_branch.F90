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

!+    Contributions to this file made by: Xiaolong Du
  
!!{
Implements a merger tree build controller class which builds trees containing only the main branch and
progenitors of the node on the main branch above a certain mass fraction.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerMainBranch">
    <description>
      A merger tree build controller class which builds trees containing only the main branch and progenitors of the node on the
      main branch above a certain mass fraction. Specifically, if a progenitor node of the node on the main branch has a mass
      below a certain mass fraction relative to the main branch node, the branch will not grow any further.
    </description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerMainBranch
     !!{     
     A merger tree build controller class which builds trees containing only the main branch and progenitors of the node on the
     main branch above a certain mass fraction. Specifically, if a progenitor node of the node on the main branch has a mass
     below a certain mass fraction relative to the main branch node, the branch will not grow any further.
     !!}
     private
     class           (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     double precision                                               :: massFraction
   contains
     final     ::                               mainBranchDestructor
     procedure :: control                    => mainBranchControl
     procedure :: branchingProbabilityObject => mainBranchBranchingProbabilityObject
     procedure :: nodesInserted              => mainBranchNodesInserted
  end type mergerTreeBuildControllerMainBranch

  interface mergerTreeBuildControllerMainBranch
     !!{
     Constructors for the {\normalfont \ttfamily mainBranch} merger tree build controller class.
     !!}
     module procedure mainBranchConstructorParameters
     module procedure mainBranchConstructorInternal
  end interface mergerTreeBuildControllerMainBranch

contains

  function mainBranchConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily mainBranch} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildControllerMainBranch)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    double precision                                                     :: massFraction

    !![
    <inputParameter>
      <name>massFraction</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
	Mass fraction relative to the descendant node on the main branch below which the progenitor branch does not grow any further.
      </description>
    </inputParameter>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=mergerTreeBuildControllerMainBranch(massFraction,mergerTreeBranchingProbability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function mainBranchConstructorParameters

  function mainBranchConstructorInternal(massFraction,mergerTreeBranchingProbability_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily mainBranch} merger tree build controller class.
    !!}
    implicit none
    type            (mergerTreeBuildControllerMainBranch)                        :: self
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    double precision                                     , intent(in   )         :: massFraction
    !![
    <constructorAssign variables="massFraction, *mergerTreeBranchingProbability_"/>
    !!]
    
    return
  end function mainBranchConstructorInternal

  subroutine mainBranchDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mainBranch} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerMainBranch), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    return
  end subroutine mainBranchDestructor

  logical function mainBranchControl(self,node,treeWalker_) result(control)
    !!{
    Skip side branches of a tree under construction if the mass of the node is below a certain fraction relative to its descendant
    on the main branch.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(mergerTreeBuildControllerMainBranch), intent(inout)           :: self    
    type (treeNode                           ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass              ), intent(inout), optional :: treeWalker_
    class(nodeComponentBasic                 )               , pointer  :: basic

    control=.true.
    ! Move to the next node in the tree while such exists, and the mass of current node is below a certain fraction
    ! relative to its descendant on the main branch.
    do while (control.and.associated(node%parent))
       basic => node%basic()
       if (basic%mass() < self%massFraction*massDescendantOnMainBranch(node)) then
          if (present(treeWalker_)) then
             control=treeWalker_%next(node)
          else
             control=.false.
          end if
       else
          return
       end if
    end do
    return
  end function mainBranchControl

  double precision function massDescendantOnMainBranch(node) result(mass)
    !!{
    Find the mass of the descendant node on the main branch.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    type (treeNode          ), intent(inout), pointer  :: node
    type (treeNode          )               , pointer  :: nodeWork
    class(nodeComponentBasic)               , pointer  :: basic

    nodeWork => node
    do while (.not.nodeWork%isOnMainBranch().and.associated(nodeWork%parent))
       nodeWork => nodeWork%parent
    end do
    basic => nodeWork%basic()
    mass  =  basic   %mass ()
    return
  end function massDescendantOnMainBranch

  function mainBranchBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerMainBranch), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function mainBranchBranchingProbabilityObject

  subroutine mainBranchNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerMainBranch), intent(inout)           :: self
    type   (treeNode                           ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                           ), intent(inout), optional :: nodeProgenitor2
    logical                                     , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine mainBranchNodesInserted
