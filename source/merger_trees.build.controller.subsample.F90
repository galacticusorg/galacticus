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
Implements a merger tree build controller class which performs subsampling of branches.
!!}

  ! Options controlling when to destroy stub branches.
  !![
  <enumeration>
    <name>destroyStubs</name>
    <description>Enumeration of options controlling when to destroy stub branches.</description>
    <encodeFunction>yes</encodeFunction>
    <visibility>private</visibility>
    <entry label="always"          />
    <entry label="never"           />
    <entry label="sideBranchesOnly"/>
  </enumeration>
  !!]

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerSubsample">
   <description>A merger tree build controller class which performs subsampling of branches.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerSubsample
     !!{     
     A merger tree build controller class which performs subsampling of branches. A branch of mass $M$ will be retained with
     probability
     \begin{equation}
       P(M) = \left\{ \begin{array}{ll} 1 &amp; \hbox{if } M \ge M_0 \\ P_0 (M/M_0)^\alpha &amp; \hbox{if } M &lt; M_0, \end{array} \right.
     \end{equation}
     where $M_0=${\normalfont \ttfamily [massThreshold]}, $P_0=${\normalfont \ttfamily [subsamplingRateAtThreshold]} and
     $\alpha=${\normalfont \ttfamily [exponent]}, otherwise being pruned. Node weights are adjusted to account for this pruning.

     If, after pruning, a section of tree is branchless, the nodes along that branch can be consolidated into fewer nodes with the
     constraint that the mass of the node increases by a fractional amount {\normalfont \ttfamily [factorMassGrowthConsolidate]}
     relative to its child node. This avoids having very long, non-branching runs of nodes with only tiny mass differences between
     each parent and child. If {\normalfont \ttfamily [factorMassGrowthConsolidate]}$\le 0$ no consolidation will be performed.     
     !!}
     private
     class           (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     double precision                                               :: massThreshold                            , subsamplingRateAtThreshold , &
          &                                                            exponent                                 , factorMassGrowthConsolidate
     type            (enumerationDestroyStubsType        )          :: destroyStubs
     integer         (kind_int8                          )          :: uniqueIDKnownMainBranchNode
  contains
     final     ::                               subsampleDestructor
     procedure :: control                    => subsampleControl
     procedure :: branchingProbabilityObject => subsampleBranchingProbabilityObject
     procedure :: nodesInserted              => subsampleNodesInserted
  end type mergerTreeBuildControllerSubsample

  interface mergerTreeBuildControllerSubsample
     !!{
     Constructors for the {\normalfont \ttfamily subsample} merger tree build controller class.
     !!}
     module procedure subsampleConstructorParameters
     module procedure subsampleConstructorInternal
  end interface mergerTreeBuildControllerSubsample

contains

  function subsampleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily subsample} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildControllerSubsample )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    double precision                                                     :: massThreshold                  , subsamplingRateAtThreshold , &
         &                                                                  exponent                       , factorMassGrowthConsolidate
    type            (varying_string                     )                :: destroyStubs

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The mass threshold, $M_0$, below which subsampling is applied.</description>
    </inputParameter>
    <inputParameter>
      <name>subsamplingRateAtThreshold</name>
      <source>parameters</source>
      <description>The subsampling rate at the mass threshold, $P_0$.</description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <description>The exponent, $\alpha$, of mass in the subsampling probability, i.e. $P(M) = P_0 (M/M_0)^\alpha$ for $M &lt; M_0$.</description>
    </inputParameter>
    <inputParameter>
      <name>factorMassGrowthConsolidate</name>
      <source>parameters</source>
      <description>The maximum factor by which the mass is allowed to grow between child and parent when consolidating nodes. A non-positive value prevents consolidation.</description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>destroyStubs</name>
      <source>parameters</source>
      <defaultValue>var_str('always')</defaultValue>
      <description>Parameter controlling when to destroy stub branches. Options are `always`, `never`, and `sideBranchesOnly`.</description>
    </inputParameter>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=mergerTreeBuildControllerSubsample(massThreshold,subsamplingRateAtThreshold,exponent,factorMassGrowthConsolidate,enumerationDestroyStubsEncode(char(destroyStubs),includesPrefix=.false.),mergerTreeBranchingProbability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function subsampleConstructorParameters

  function subsampleConstructorInternal(massThreshold,subsamplingRateAtThreshold,exponent,factorMassGrowthConsolidate,destroyStubs,mergerTreeBranchingProbability_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily subsample} merger tree build controller class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeBuildControllerSubsample )                        :: self
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    double precision                                     , intent(in   )         :: massThreshold                  , subsamplingRateAtThreshold , &
         &                                                                          exponent                       , factorMassGrowthConsolidate
    type            (enumerationDestroyStubsType        ), intent(in   )         :: destroyStubs
    !![
    <constructorAssign variables="massThreshold, subsamplingRateAtThreshold, exponent, factorMassGrowthConsolidate, destroyStubs, *mergerTreeBranchingProbability_"/>
    !!]

    if (self%destroyStubs == destroyStubsNever .and. self%factorMassGrowthConsolidate > 0.0d0) &
         & call Error_Report("branch consolidation is not supported when branch stubs are not to be detroyed"//{introspection:location})
    ! Initialize the known main branch unique ID to an impossible value.
    self%uniqueIDKnownMainBranchNode=-1_kind_int8
    return
  end function subsampleConstructorInternal

  subroutine subsampleDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily subsample} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerSubsample), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    return
  end subroutine subsampleDestructor

  logical function subsampleControl(self,node,treeWalker_)
    !!{
    Subsample branches of a tree under construction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: Error           , only : Error_Report
    implicit none
    class           (mergerTreeBuildControllerSubsample), intent(inout)           :: self    
    type            (treeNode                          ), intent(inout), pointer  :: node
    class           (mergerTreeWalkerClass             ), intent(inout), optional :: treeWalker_
    type            (treeNode                          )               , pointer  :: nodeNext       , nodeChild     , &
         &                                                                           nodeParent     , nodeGrandchild
    class           (nodeComponentBasic                )               , pointer  :: basic          , basicParent
    double precision                                                              :: rateSubsampling
    integer         (c_size_t                          )                          :: countNodes
    logical                                                                       :: finished       , destroyStub,isMainBranch

    ! If the node has been left in place as a stub, never process it.
    ! The node which we return to the tree builder must be one that we have determined will not be pruned, since this node will be
    ! fully-processed by the tree builder. Therefore, if we prune a node we must check for pruning of the next node, and so on
    ! until we reach a node that is not pruned.
    finished=.false.
    do while (.not.finished)
       finished        =.true.
       subsampleControl=.true.
       ! Root node is not eligible for pruning.
       if (.not.associated(node%parent)) return
       ! Set the subsampling weight for this node to equal that of its parent.
       call node%subsamplingWeightSet(node%parent%subsamplingWeight())
       ! Primary progenitors are not eligible for pruning.
       if (node%isPrimaryProgenitor()) return
       ! Nodes above the mass threshold are not eligible for pruning.
       basic => node%basic()
       if (basic%mass() >= self%massThreshold) return
       ! Compute subsampling rate, perform sampling.
       rateSubsampling=+self%subsamplingRateAtThreshold &
            &          *(                               &
            &            +basic%mass         ()         &
            &            /self %massThreshold           &
            &           )**self%exponent
       if (node%hostTree%randomNumberGenerator_%uniformSample() < rateSubsampling) then
          ! Node is to be kept - increase its weight to account for corresponding branches which will have been lost.
          call node%subsamplingWeightSet(node%subsamplingWeight()/rateSubsampling)
       else
          ! Prune the node.
          !! Get the next node to walk to in the tree.
          if (present(treeWalker_)) then
             subsampleControl =  treeWalker_%next(nodeNext)
          else
             subsampleControl =  .false.
             nodeNext         => null()
          end if
          !! Determine if the stub should be destroyed.
          select case (self%destroyStubs%ID)
          case (destroyStubsNever           %ID)
             destroyStub=.false.
          case (destroyStubsAlways          %ID)
             destroyStub=.true.
          case (destroyStubsSideBranchesOnly%ID)
             ! Determine if the node's parent is on the main branch. We implement this test directly here as we can often exploit
             ! the order in which the tree is built to speed up this test.
             nodeParent   => node%parent
             isMainBranch =  .true.
             do while (associated(nodeParent%parent))
                ! If we have reached another node known to be on the main branch, then our node must be on the main branch also.
                if (nodeParent%uniqueID() == self%uniqueIDKnownMainBranchNode) exit
                ! Test if we are the primary progenitor - if we are not then we are not on the main branch.
                if (.not.nodeParent%isPrimaryProgenitor()) then
                   isMainBranch=.false.
                   exit
                end if
                ! Move to the parent node and test again.
                nodeParent => nodeParent%parent
             end do
             ! If our parent node was on the main branch, then its parent is also. Record the unique ID of that node - if we find
             ! it again we know we are on the main branch.
             if (isMainBranch.and.associated(node%parent%parent)) self%uniqueIDKnownMainBranchNode=node%parent%parent%uniqueID()
             ! Destroy the stub only if it is not on the main branch.
             destroyStub=.not.isMainBranch
          case default
             destroyStub=.false.
             call Error_Report('unknown `destroyStubs` option'//{introspection:location})
          end select
          !! Decouple the node from the tree.
          if (destroyStub) then
             nodeParent => node      %parent
             nodeChild  => nodeParent%firstChild
             do while (.not.associated(nodeChild%sibling,node))
                nodeChild => nodeChild%sibling
             end do
             nodeChild%sibling => node%sibling
             ! Destroy and deallocate the node.
             call node%destroy()
             deallocate(node)
             ! Determine if we can consolidate any nodes down the parent branch.
             if (self%factorMassGrowthConsolidate > 0.0d0) then
                ! Seek down through the branch until which find a node which either has a sibling (so can't be consolidated), or
                ! which has a mass sufficiently different from that of the starting node. Count how many such nodes we find.
                nodeChild   => nodeParent
                basic       => nodeChild %basic()
                basicParent => nodeParent%basic()
                countNodes  =  0_c_size_t
                do while (                                                              &
                     &         associated(nodeChild%firstChild)                         &
                     &    .and.                                                         &
                     &    .not.associated(nodeChild%sibling   )                         &
                     &    .and.                                                         &
                     &      basic      %mass()*(1.0d0+self%factorMassGrowthConsolidate) &
                     &     >                                                            &
                     &      basicParent%mass()                                          &
                     &   )
                   nodeChild  => nodeChild%firstChild
                   basic      => nodeChild%basic     ()
                   countNodes =  countNodes+1_c_size_t
                end do
                ! If we have found nodes that can be consolidated, remove the intervening nodes.
                if (countNodes > 1_c_size_t) then
                   nodeChild => nodeParent%firstChild
                   do while (countNodes > 1_c_size_t)
                      nodeGrandchild => nodeChild%firstChild
                      call nodeChild%destroy()
                      deallocate(nodeChild)
                      countNodes = countNodes-1_c_size_t
                      nodeChild => nodeGrandchild
                   end do
                   nodeParent%firstChild => nodeChild
                   nodeChild %parent     => nodeParent
                   do while (associated(nodeChild%sibling))
                      nodeChild        => nodeChild %sibling
                      nodeChild%parent => nodeParent
                   end do
                end if
             end if
          else
             ! Stubs are not being destroyed. Mark the stub by assigning a negative subsampling weight to it.
             call node%subsamplingWeightSet(-1.0d0)
          end if
          ! Set the current node to the next node in the tree walk.
          node => nodeNext
          ! We pruned a node. Therefore, if a next node was found we must now check whether we want to prune it too.
          finished=.not.subsampleControl
       end if
    end do
    return
  end function subsampleControl

  function subsampleBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerSubsample ), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function subsampleBranchingProbabilityObject

  subroutine subsampleNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerSubsample), intent(inout)           :: self
    type   (treeNode                          ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                          ), intent(inout), optional :: nodeProgenitor2
    logical                                    , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine subsampleNodesInserted
