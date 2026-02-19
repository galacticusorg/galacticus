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
  A merger tree builder class that enforces constraints on the merger tree built by some other builder.
  !!}

  use, intrinsic :: ISO_C_Binding      , only : c_size_t             , c_long
  use            :: Merger_Tree_Filters, only : mergerTreeFilterClass
  
  type, public :: constrainedBuilderList
     class(mergerTreeBuilderClass), pointer :: mergerTreeBuilder_ => null()
     class(mergerTreeFilterClass ), pointer :: mergerTreeFilter_  => null()
     type (constrainedBuilderList), pointer :: next               => null()
  end type constrainedBuilderList

  !![
  <mergerTreeBuilder name="mergerTreeBuilderConstrained">
   <description>
    A merger tree builder class that enforces constraints on the merger tree built by some other builder.
   </description>
   <linkedList type="constrainedBuilderList" variable="mergerTreeBuilders" next="next" object="mergerTreeBuilder_ mergerTreeFilter_" objectType="mergerTreeBuilderClass mergerTreeFilterClass"/>
  </mergerTreeBuilder>
  !!]
  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderConstrained
     !!{
     A merger tree builder class that enforces constraints on the merger tree built by some other builder.
     !!}
     private
     type  (constrainedBuilderList), pointer :: mergerTreeBuilders => null()
     integer(c_size_t             )          :: trialCountMaximum
   contains
     final     ::          constrainedDestructor
     procedure :: build => constrainedBuild
  end type mergerTreeBuilderConstrained

  interface mergerTreeBuilderConstrained
     !!{
     Constructors for the \refClass{mergerTreeBuilderConstrained} merger tree builder class.
     !!}
     module procedure constrainedConstructorParameters
     module procedure constrainedConstructorInternal
  end interface mergerTreeBuilderConstrained

  integer(c_long) :: constrainedSeed=0_c_long
  
contains

  function constrainedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the constrained merger tree building class which reads parameters from a provided parameter list.
    !!}
    implicit none
    type   (mergerTreeBuilderConstrained)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    type   (constrainedBuilderList      ), pointer       :: mergerTreeBuilder_
    integer                                              :: i

    if (parameters%copiesCount('mergerTreeBuilder') == 0                                         ) &
         & call Error_Report('at least one [mergerTreeBuilder] must be specified'                        //{introspection:location})
    if (parameters%copiesCount('mergerTreeBuilder') /= parameters%copiesCount('mergerTreeFilter')) &
         & call Error_Report('number of [mergerTreeBuilder] and [mergerTreeFilter] parameters must match'//{introspection:location})
    mergerTreeBuilder_ => null()
    do i=1,parameters%copiesCount('mergerTreeBuilder')
       if (associated(mergerTreeBuilder_)) then
          allocate(mergerTreeBuilder_%next)
          mergerTreeBuilder_ => mergerTreeBuilder_%next
       else
          allocate(self%mergerTreeBuilders)
          mergerTreeBuilder_ => self%mergerTreeBuilders
       end if
       !![
       <objectBuilder class="mergerTreeBuilder" name="mergerTreeBuilder_%mergerTreeBuilder_" source="parameters" copy="i" />
       !!]
    end do
    mergerTreeBuilder_ => self%mergerTreeBuilders
    do i=1,parameters%copiesCount('mergerTreeFilter')
       !![
       <objectBuilder class="mergerTreeFilter" name="mergerTreeBuilder_%mergerTreeFilter_" source="parameters" copy="i" />
       !!]
       mergerTreeBuilder_ => mergerTreeBuilder_%next
    end do
    !![
    <inputParameter>
      <name>trialCountMaximum</name>
      <variable>self%trialCountMaximum</variable>
      <source>parameters</source>
      <description>The maximum number of trials to attempt before failing.</description>
      <defaultValue>huge(1_c_size_t)</defaultValue>
    </inputParameter>
    <inputParametersValidate source="parameters" multiParameters="mergerTreeBuilder, mergerTreeFilter"/>
    !!]
    return
  end function constrainedConstructorParameters

  function constrainedConstructorInternal(mergerTreeBuilders,trialCountMaximum) result(self)
    !!{
    Internal constructor for the constrained merger tree building class.
    !!}
    implicit none
    type   (mergerTreeBuilderConstrained)                         :: self
    type   (constrainedBuilderList      ), intent(in   ), target  :: mergerTreeBuilders
    integer(c_size_t                    ), intent(in   )          :: trialCountMaximum
    type   (constrainedBuilderList      )               , pointer :: mergerTreeBuilder_
    !![
    <constructorAssign variables="trialCountMaximum"/>
    !!]
    
    self %mergerTreeBuilders => mergerTreeBuilders
    mergerTreeBuilder_       => mergerTreeBuilders
    do while (associated(mergerTreeBuilder_))
       !![
       <referenceCountIncrement owner="mergerTreeBuilder_" object="mergerTreeBuilder_"/>
       <referenceCountIncrement owner="mergerTreeBuilder_" object="mergerTreeFilter_" />
       !!]
       mergerTreeBuilder_ => mergerTreeBuilder_%next
    end do 
    return
  end function constrainedConstructorInternal

  subroutine constrainedDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuilderConstrained} merger tree builder class.
    !!}
    implicit none
    type(mergerTreeBuilderConstrained), intent(inout) :: self
    type(constrainedBuilderList      ), pointer       :: mergerTreeBuilder_, mergerTreeBuilderNext

    if (associated(self%mergerTreeBuilders)) then
       mergerTreeBuilder_ => self%mergerTreeBuilders
       do while (associated(mergerTreeBuilder_))
          mergerTreeBuilderNext => mergerTreeBuilder_%next
          !![
          <objectDestructor name="mergerTreeBuilder_%mergerTreeBuilder_"/>
          <objectDestructor name="mergerTreeBuilder_%mergerTreeFilter_" />
          !!]
          deallocate(mergerTreeBuilder_)
          mergerTreeBuilder_ => mergerTreeBuilderNext
       end do
    end if

    return
  end subroutine constrainedDestructor

  subroutine constrainedBuild(self,tree)
    !!{
    Build a constrained merger tree.
    !!}
    use            :: Display            , only : displayMessage          , displayIndent, displayUnindent, verbosityLevelWorking
    use            :: Kind_Numbers       , only : kind_int8
    use            :: Galacticus_Nodes   , only : treeNodeLinkedList
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    use            :: String_Handling    , only : operator(//)
    use            :: ISO_Varying_String , only : operator(//)            , var_str
    implicit none
    class  (mergerTreeBuilderConstrained), intent(inout), target  :: self
    type   (mergerTree                  ), intent(inout), target  :: tree
    type   (treeNode                    )               , pointer :: nodeChild         , nodeNext
    type   (treeNodeLinkedList          )               , pointer :: nodeLeavesHead    , nodeLeavesCurrent, &
         &                                                           nodeLeavesNext
    type   (constrainedBuilderList      )               , pointer :: mergerTreeBuilder_
    logical                                                       :: passes
    integer(c_size_t                    )                         :: trialCount        , trialCountTotal
    integer                                                       :: stage
    integer(kind_int8                   )                         :: treeIndex
    type   (mergerTreeWalkerAllNodes    )                         :: treeWalker
    
    call displayIndent('Begin constrained tree build',verbosityLevelWorking)
    ! Store the tree index so that we can restore this after building our tree.    
    treeIndex=tree%index
    ! Iterate over builders.
    trialCountTotal    =  1_c_size_t
    stage              =  0
    mergerTreeBuilder_ => self%mergerTreeBuilders
    do while (associated(mergerTreeBuilder_))
       stage=stage+1
       call displayIndent(var_str('Begin stage ')//stage,verbosityLevelWorking)
       ! Initialize a counter for the number of trials used.
       trialCount=0_c_size_t
       ! Find leaf nodes of the current (partial) tree. This will allow us to restore state back to this current tree on a failed
       ! build attempt.
       nodeLeavesHead    => null()
       nodeLeavesCurrent => null()
       treeWalker        =  mergerTreeWalkerAllNodes(tree)
       do while (treeWalker%next(nodeChild))
          if (.not.associated(nodeChild%firstChild)) then
             if (.not.associated(nodeLeavesHead)) then
                allocate(nodeLeavesHead        )
                nodeLeavesCurrent => nodeLeavesHead
             else
                allocate(nodeLeavesCurrent%next)
                nodeLeavesCurrent => nodeLeavesCurrent%next
             end if
             nodeLeavesCurrent%node => nodeChild
          end if
       end do
       ! Iterate until we find an acceptable tree.
       passes=.false.
       do while (.not.passes)
          ! Increment counters.
          trialCount=trialCount+1_c_size_t
          if (mod(trialCount,100_c_size_t) == 0_c_size_t) &
               & call displayMessage(var_str('trial number ')//trialCount,verbosityLevelWorking)
          if (trialCount > self%trialCountMaximum) &
               & call Error_Report('maximum number of trials exceeded'//{introspection:location})
          ! Set tree seed and random number generator to ensure that a different random realization is made each time.
          !$omp atomic
          constrainedSeed      =constrainedSeed+1_c_long
          tree           %index=constrainedSeed
          call tree%randomNumberGenerator_%seedSet(constrainedSeed,offset=.true.)
          ! Destroy any branches from the previous build attempt.
          nodeLeavesCurrent => nodeLeavesHead
          do while (associated(nodeLeavesCurrent))
             nodeChild => nodeLeavesCurrent%node%firstChild
             do while (associated(nodeChild))
                nodeNext => nodeChild%sibling
                call nodeChild%destroyBranch()
                deallocate(nodeChild)
                nodeChild => nodeNext
             end do
             nodeLeavesCurrent%node%firstChild => null()
             nodeLeavesCurrent => nodeLeavesCurrent%next
          end do
          ! Build the tree.
          call mergerTreeBuilder_%mergerTreeBuilder_%build(tree)
          ! Determine if the tree is acceptable.
          passes=mergerTreeBuilder_%mergerTreeFilter_%passes(tree)
          if (passes) call displayMessage(var_str('success after ')//trialCount//' trials',verbosityLevelWorking)
       end do
       ! Destroy the list of leaf nodes.
       nodeLeavesCurrent => nodeLeavesHead
       do while (associated(nodeLeavesCurrent))
          nodeLeavesNext => nodeLeavesCurrent%next
          deallocate(nodeLeavesCurrent)
          nodeLeavesCurrent => nodeLeavesNext
       end do
       ! Accumulate the number of trials.
       trialCountTotal=+trialCountTotal &
            &          *trialCount
       ! Move to the next builder
       mergerTreeBuilder_ => mergerTreeBuilder_%next
       call displayUnindent('done',verbosityLevelWorking)
    end do
    ! Adjust the tree weight to account for the number of trials taken to construct this tree.
    tree%volumeWeight=tree%volumeWeight/dble(trialCountTotal)
    ! Restore the tree index.
    tree%index=treeIndex
    ! Report on the effective number of trials used.
    call displayMessage (var_str('Effective number of trails = ')//trialCountTotal,verbosityLevelWorking)
    call displayUnindent('done'                                                   ,verbosityLevelWorking)
    return
  end subroutine constrainedBuild

