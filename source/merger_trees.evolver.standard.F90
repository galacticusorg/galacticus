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
  Implements the standard class for evolving merger trees.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions        , cosmologyFunctionsClass
  use :: Galactic_Structure_Solvers, only : galacticStructureSolver   , galacticStructureSolverClass
  use :: Galacticus_Nodes          , only : treeNode
  use :: Kind_Numbers              , only : kind_int8
  use :: Merger_Tree_Initialization, only : mergerTreeInitializorClass
  use :: Merger_Tree_Timesteps     , only : mergerTreeEvolveTimestep  , mergerTreeEvolveTimestepClass
  use :: Merger_Trees_Evolve_Node  , only : mergerTreeNodeEvolver     , mergerTreeNodeEvolverClass
  use :: Galactic_Structure        , only : galacticStructureClass

  ! Structure used to store list of nodes for deadlock reporting.
  type :: deadlockList
     type   (deadlockList  ), pointer :: next      => null()
     type   (treeNode      ), pointer :: nodeLock           , node
     integer(kind=kind_int8)          :: treeIndex
     type   (varying_string)          :: lockType
  end type deadlockList

  !![
  <mergerTreeEvolver name="mergerTreeEvolverStandard">
   <description>
    The standard merger tree evolver. Each merger tree forest is evolved by repeatedly walking the
    trees and evolving each node forward in time by some timestep $\Delta t$. Nodes are evolved
    individually such that nodes in different branches of a tree may have reached different cosmic
    times at any given point in the execution of \glc. Each node is evolved over the interval
    $\Delta t$ using an adaptive \gls{ode} solver, which adjusts the smaller timesteps, $\delta
    t$, taken in evolving the system of \glspl{ode} to maintain a specified precision.
    
    The choice of $\Delta t$ then depends on other considerations. For example, a node should not
    be evolved beyond the time at which it is due to merge with another galaxy. Also, we typically
    don't want satellite nodes to evolve too far ahead of their host node, such that any
    interactions between satellite and host occur (near) synchronously.
    
    The following timestep criteria ensure that tree evolution occurs in a way which correctly
    preserves tree structure and ordering of interactions between \glspl{node}. All criteria are
    considered and the largest $\Delta t$ consistent with all criteria is selected.
    
    \begin{description}
     \item[Branch segment criterion] For \glspl{node} which are the \gls{primary progenitor} of
     their \gls{parent}, the ``branch segment'' criterion asserts that
      \begin{equation}
       \Delta t \le t_\mathrm{parent} - t
      \end{equation}
     where $t$ is current time in the \gls{node} and $t_\mathrm{parent}$ is the time of the
     \gls{parent} \gls{node}. This ensures that \gls{primary progenitor} \glspl{node} to not evolve
     beyond the time at which their \gls{parent} (which they will replace) exists.  If this
     criterion is the limiting criteria for $\Delta t$ then the \gls{node} will be promoted to
     replace its \gls{parent} at the end of the timestep.
    
     \item[Parent criterion] For \glspl{node} which are satellites in a hosting \gls{node} the
     ``\gls{parent}'' timestep criterion asserts that
      \begin{eqnarray}
       \Delta t &amp;\le&amp; t_\mathrm{host}, \\
       \Delta t &amp;\le&amp; \epsilon_\mathrm{host} (a/\dot{a}),
      \end{eqnarray}
     where $t_\mathrm{host}=${\normalfont \ttfamily [timestepHostAbsolute]},
     $\epsilon_\mathrm{host}=${\normalfont \ttfamily [timestepHostRelative]}, and $a$ is expansion
     factor. These criteria are intended to prevent a satellite for evolving too far ahead of the
     host node before the host is allowed to ``catch up''.
    
     \item[Satellite criterion] For \glspl{node} which host satellite \glspl{node}, the
     ``satellite'' criterion asserts that
      \begin{equation}
       \Delta t \le \hbox{min}(t_\mathrm{satellite}) - t,
      \end{equation}
     where $t$ is the time of the host \gls{node} and $t_\mathrm{satellite}$ are the times of all
     satellite \glspl{node} in the host. This criterion prevents a host from evolving ahead of any
     satellites.
    
     \item[Sibling criterion] For \glspl{node} which are \glspl{primary progenitor}, the
     ``sibling'' criterion asserts that
      \begin{equation}
       \Delta t \le \hbox{min}(t_\mathrm{sibling}) - t,
      \end{equation}
     where $t$ is the time of the host \gls{node} and $t_\mathrm{sibling}$ are the times of all
     siblings of the \gls{node}. This criterion prevents a \gls{node} from reaching its
     \gls{parent} (and being promoted to replace it) before all of its siblings have reach the
     \gls{parent} and have become satellites within it.
    
     \item[Mergee criterion] For \glspl{node} with \gls{mergee} \glspl{node}, the ``\gls{mergee}''
     criterion asserts that
      \begin{equation}
       \Delta t \le \hbox{min}(t_\mathrm{merge}) - t,
      \end{equation}
     where $t$ is the time of the host \gls{node} and $t_\mathrm{merge}$ are the times at which
     the \glspl{mergee} will merge. This criterion prevents a \gls{node} from evolving past the
     time at which a merger event takes place.
    \end{description}
    </description>
  </mergerTreeEvolver>
  !!]
  type, extends(mergerTreeEvolverClass) :: mergerTreeEvolverStandard
     !!{
     Implementation of the standars merger tree evolver.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (mergerTreeEvolveTimestepClass), pointer :: mergerTreeEvolveTimestep_        => null()
     class           (galacticStructureSolverClass ), pointer :: galacticStructureSolver_         => null()
     class           (mergerTreeNodeEvolverClass   ), pointer :: mergerTreeNodeEvolver_           => null()
     class           (mergerTreeInitializorClass   ), pointer :: mergerTreeInitializor_           => null()
     class           (galacticStructureClass       ), pointer :: galacticStructure_               => null()
     logical                                                  :: allTreesExistAtFinalTime                  , dumpTreeStructure    , &
          &                                                      backtrackToSatellites
     double precision                                         :: timestepHostAbsolute                      , timestepHostRelative , &
          &                                                      fractionTimestepSatelliteMinimum
     type            (deadlockList                 ), pointer :: deadlockHeadNode                 => null()
   contains
     !![
     <methods>
       <method description="Find the time to which a node can be evolved." method="timeEvolveTo" />
       <method description="Add a node to the deadlock list." method="deadlockAddNode" />
       <method description="Output a description of a deadlocked tree." method="deadlockOutputTree" />
     </methods>
     !!]
     final     ::                       standardDestructor
     procedure :: evolve             => standardEvolve
     procedure :: timeEvolveTo       => standardTimeEvolveTo
     procedure :: deadlockAddNode    => standardDeadlockAddNode
     procedure :: deadlockOutputTree => standardDeadlockOutputTree
  end type mergerTreeEvolverStandard

  interface mergerTreeEvolverStandard
     !!{
     Constructors for the {\normalfont \ttfamily standard} merger tree evolver.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeEvolverStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily standard} merger tree evolver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolverStandard    )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass), pointer       :: mergerTreeEvolveTimestep_
    class           (galacticStructureSolverClass ), pointer       :: galacticStructureSolver_
    class           (mergerTreeNodeEvolverClass   ), pointer       :: mergerTreeNodeEvolver_
    class           (mergerTreeInitializorClass   ), pointer       :: mergerTreeInitializor_
    class           (galacticStructureClass       ), pointer       :: galacticStructure_
    logical                                                        :: allTreesExistAtFinalTime        , dumpTreeStructure         , &
         &                                                            backtrackToSatellites
    double precision                                               :: timestepHostRelative            , timestepHostAbsolute      , &
         &                                                            fractionTimestepSatelliteMinimum

    !![
    <inputParameter>
      <name>allTreesExistAtFinalTime</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not all merger trees are expected to exist at the final requested output time. If set to false,
         then trees which finish before a given output time will be ignored.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dumpTreeStructure</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether merger tree structure should be dumped to a \href{http://www.graphviz.org/}{\normalfont \scshape dot} file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timestepHostRelative</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The maximum allowed relative timestep for node evolution relative to the time of the host halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timestepHostAbsolute</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The maximum allowed absolute timestep (in Gyr) for node evolution relative to the time of the host halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionTimestepSatelliteMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum fraction of the timestep imposed by the ``satellite in host'' criterion to evolve over. If the timestep allowed is smaller than this fraction, the actual timestep will be reduced to zero. This avoids forcing satellites to take a large number of very small timesteps, and instead defers evolving a satellite until a large timestep can be taken.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>backtrackToSatellites</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, after successfully evolving a node with satellites, revisit the satellites and attempt to evolve them again.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    ! A galacticStructureSolver is built here. Even though this is not called explicitly by this mergerTreeEvolver, the
    ! galacticStructureSolver is expected to hook itself to any events which will trigger a change in galactic structure.
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="mergerTreeEvolveTimestep" name="mergerTreeEvolveTimestep_" source="parameters"/>
    <objectBuilder class="galacticStructureSolver"  name="galacticStructureSolver_"  source="parameters"/>
    <objectBuilder class="mergerTreeNodeEvolver"    name="mergerTreeNodeEvolver_"    source="parameters"/>
    <objectBuilder class="mergerTreeInitializor"    name="mergerTreeInitializor_"    source="parameters"/>
    <objectBuilder class="galacticStructure"        name="galacticStructure_"        source="parameters"/>
    !!]
    self=mergerTreeEvolverStandard(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,galacticStructureSolver_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="mergerTreeEvolveTimestep_"/>
    <objectDestructor name="mergerTreeNodeEvolver_"   />
    <objectDestructor name="galacticStructureSolver_" />
    <objectDestructor name="mergerTreeInitializor_"   />
    <objectDestructor name="galacticStructure_"       />
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,galacticStructureSolver_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily standard} merger tree evolver class.
    !!}
    implicit none
    type            (mergerTreeEvolverStandard    )                        :: self
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass), intent(in   ), target :: mergerTreeEvolveTimestep_
    class           (galacticStructureSolverClass ), intent(in   ), target :: galacticStructureSolver_
    class           (mergerTreeNodeEvolverClass   ), intent(in   ), target :: mergerTreeNodeEvolver_
    class           (mergerTreeInitializorClass   ), intent(in   ), target :: mergerTreeInitializor_
    class           (galacticStructureClass       ), intent(in   ), target :: galacticStructure_
    logical                                        , intent(in   )         :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                    backtrackToSatellites
    double precision                               , intent(in   )         :: timestepHostRelative            , timestepHostAbsolute, &
         &                                                                    fractionTimestepSatelliteMinimum
    !![
    <constructorAssign variables="allTreesExistAtFinalTime, dumpTreeStructure, timestepHostRelative, timestepHostAbsolute, fractionTimestepSatelliteMinimum, backtrackToSatellites, *cosmologyFunctions_, *mergerTreeNodeEvolver_, *mergerTreeEvolveTimestep_, *mergerTreeInitializor_, *galacticStructureSolver_, *galacticStructure_"/>
    !!]

    self%deadlockHeadNode => null()
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily standard}m erger tree evolver class.
    !!}
    implicit none
    type(mergerTreeEvolverStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%mergerTreeEvolveTimestep_"/>
    <objectDestructor name="self%galacticStructureSolver_" />
    <objectDestructor name="self%mergerTreeNodeEvolver_"   />
    <objectDestructor name="self%mergerTreeInitializor_"   />
    <objectDestructor name="self%galacticStructure_"       />
    !!]
    return
  end subroutine standardDestructor

  subroutine standardEvolve(self,tree,timeEnd,treeDidEvolve,suspendTree,deadlockReporting,systemClockMaximum,initializationLock,status)
    !!{
    Evolves all properties of a merger tree to the specified time.
    !!}
    use    :: Display                            , only : displayIndent               , displayMessage                    , displayUnindent          , displayVerbosity           , &
         &                                                displayGreen                , displayReset
    use    :: Galacticus_Error                   , only : Galacticus_Error_Report     , errorStatusSuccess
    use    :: Galacticus_Nodes                   , only : interruptTask               , mergerTree                        , nodeComponentBasic       , nodeEvent                  , &
          &                                               nodeEventBranchJumpInterTree, nodeEventSubhaloPromotionInterTree, treeNode
    use    :: Merger_Tree_Timesteps              , only : timestepTask
    use    :: Merger_Tree_Walkers                , only : mergerTreeWalkerAllNodes
    use    :: Merger_Trees_Dump                  , only : Merger_Tree_Dump
    use    :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsDeadlocked  , deadlockStatusIsNotDeadlocked     , deadlockStatusIsReporting, deadlockStatusIsSuspendable
    !$ use :: OMP_Lib                            , only : OMP_Set_Lock                , OMP_Unset_Lock                    , omp_lock_kind
    use    :: String_Handling                    , only : operator(//)
    implicit none
    class           (mergerTreeEvolverStandard )                    , intent(inout) :: self
    integer                                     , optional          , intent(  out) :: status
    type            (mergerTree                )           , target , intent(inout) :: tree
    double precision                                                , intent(in   ) :: timeEnd
    logical                                                         , intent(  out) :: treeDidEvolve                     , suspendTree
    logical                                                         , intent(in   ) :: deadlockReporting
    integer         (kind_int8                 ), optional          , intent(in   ) :: systemClockMaximum
    integer         (omp_lock_kind             ), optional          , intent(inout) :: initializationLock
    type            (treeNode                  )           , pointer                :: nodeLock                          , nodeNext         , &
         &                                                                             nodeParent                        , node             , &
         &                                                                             nodeSiblingNext                   , nodeParentNext
    class           (nodeEvent                 )           , pointer                :: event
    double precision                            , parameter                         :: timeTolerance              =1.0d-5
    double precision                            , parameter                         :: largeTime                  =1.0d10
    procedure       (interruptTask             )           , pointer                :: interruptProcedure
    procedure       (timestepTask              )           , pointer                :: timestepTask_
    class           (*                         )           , pointer                :: timestepSelf
    integer                                     , parameter                         :: verbosityLevel             =3
    class           (nodeComponentBasic        )           , pointer                :: basicBase                         , basicParent      , &
         &                                                                             basic
    type            (mergerTree                )           , pointer                :: currentTree
    type            (mergerTreeWalkerAllNodes  )                                    :: treeWalker
    integer                                                                         :: statusDeadlock                    , nodesEvolvedCount, &
         &                                                                             nodesTotalCount                   , treeWalkCount    , &
         &                                                                             treeWalkCountPreviousOutput
    double precision                                                                :: earliestTimeInTree                , timeEndThisNode  , &
         &                                                                             finalTimeInTree
    character       (len=24                    )                                    :: label
    character       (len=35                    )                                    :: message
    type            (varying_string            )                                    :: lockType                          , vMessage
    logical                                                                         :: anyTreeExistsAtOutputTime         , hasIntertreeEvent, &
         &                                                                             hasParent                         , treeLimited      , &
         &                                                                             nodeProgressed                    , nextNodeFound    , &
         &                                                                             didEvolve                         , interrupted      , &
         &                                                                             nodesRemain

    ! Iterate through all trees.
    suspendTree               =  .false.
    anyTreeExistsAtOutputTime =  .false.
    treeDidEvolve             =  .false.
    currentTree               => tree
    do while (associated(currentTree))
       ! Skip empty trees.
       if (associated(currentTree%nodeBase)) then
          ! Initialize the tree if necessary.
          !$ if (present(initializationLock)) call OMP_Set_Lock  (initializationLock)
          call self%mergerTreeInitializor_%initialize(currentTree,timeEnd)
          !$ if (present(initializationLock)) call OMP_Unset_Lock(initializationLock)
          ! Check that the output time is not after the end time of this tree.
          basicBase => currentTree%nodeBase%basic()
          if (timeEnd > basicBase%time()) then
             ! Final time is exceeded. Check if by a significant factor.
             if (timeEnd > basicBase%time()*(1.0d0+timeTolerance)) then
                ! Exceeded by a significant factor - report an error. Check if such behavior is expected.
                if (self%allTreesExistAtFinalTime) then
                   ! It is not, write an error and exit.
                   vMessage='requested time exceeds the final time in the tree'//char(10)
                   vMessage=vMessage//displayGreen()//' HELP:'//displayReset()//' If you expect that not all trees will exist at the latest requested'//char(10)
                   vMessage=vMessage//                                         '       output time (this can happen when using trees extracted from N-body'//char(10)
                   vMessage=vMessage//                                         '       simulations for example) set the following in your input parameter file:'//char(10)//char(10)
                   vMessage=vMessage//                                         '         <allTreesExistAtFinalTime value="false" />'//char(10)
                   call Galacticus_Error_Report(vMessage//{introspection:location})
                end if
             else
                ! Not exceeded by a significant factor (can happen due to approximation errors). Unless there is an event
                ! associated with this node at the current time, simply reset to actual time requested.
                event => currentTree%nodeBase%event
                do while (associated(event))
                   if (event%time == basicBase%time()) then
                      vMessage=          'requested time exceeds the final time in the tree by a small factor'  //char(10)
                      vMessage=vMessage//'refusing to adjust the final time in the tree due to associated event'//char(10)
                      write (label,'(e24.16)') timeEnd
                      vMessage=vMessage//'  requested time: '//trim(label)//' Gyr'//char(10)
                      write (label,'(e24.16)') basicBase%time()
                      vMessage=vMessage//'      final time: '//trim(label)//' Gyr'//char(10)
                      write (label,'(e24.16)') event%time
                      vMessage=vMessage//'      event time: '//trim(label)//' Gyr'//char(10)
                      vMessage=vMessage//'      event ID  : '//event%ID           //char(10)
                      vMessage=vMessage//displayGreen()//' HELP:'//displayReset()//' if you are reading merger trees from file and are attempting to'//char(10)
                      vMessage=vMessage//                                          '       output at a "snapshot time" consider setting:'                  //char(10)
                      vMessage=vMessage//                                          '           <mergerTreeReadOutputTimeSnapTolerance value="1.0e-3"/>'    //char(10)
                      vMessage=vMessage//                                          '       or similar in your parameter file to ensure that nodes exist'   //char(10)
                      vMessage=vMessage//                                          '       precisely at the output times you request'
                      call Galacticus_Error_Report(vMessage//{introspection:location})
                   end if
                   event => event%next
                end do
                call basicBase%timeSet(timeEnd)
                anyTreeExistsAtOutputTime=.true.
             end if
          else
             anyTreeExistsAtOutputTime=.true.
          end if
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    ! If none of these trees exist at the output time, check if they contain any inter-tree events. If they do, we need to evolve
    ! the tree anyway, as it interacts with another tree that may exist at the output time. Otherwise, we can ignore this tree.
    if (.not.anyTreeExistsAtOutputTime) then
       ! Walk over all trees in the forest.
       treeWalker       =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
       hasInterTreeEvent=.false.
       do while (treeWalker%next(node).and..not.hasIntertreeEvent)
          ! Iterate over events.
          event => node%event
          do while (associated(event).and..not.hasIntertreeEvent)
             select type (event)
             type is (nodeEventSubhaloPromotionInterTree)
                hasIntertreeEvent=.true.
             type is (nodeEventBranchJumpInterTree      )
                hasIntertreeEvent=.true.
             end select
             event => event%next
          end do
       end do
       if (.not.hasInterTreeEvent) then
          ! Mark the tree as evolved here, as the only reason that we did not evolve it was the given time target.
          treeDidEvolve=.true.
          return
       end if
    end if
    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible.
    didEvolve                  =.true.
    treeWalkCount              =0
    treeWalkCountPreviousOutput=0
    outerLoop: do while (didEvolve) ! Keep looping through the tree until we make a pass during which no nodes were evolved.
       ! Flag that no nodes have been evolved yet.
       didEvolve=.false.
       ! Increment tree walks counter.
       treeWalkCount=treeWalkCount+1
       ! Reset tree progress variables.
       nodesEvolvedCount =0
       nodesTotalCount   =0
       earliestTimeInTree=largeTime
       ! Set the deadlock status to deadlocked initially, unless we have been specifically asked for deadlock reporting via a
       ! function argument, in which case set to reporting status.
       statusDeadlock=deadlockStatusIsDeadlocked
       if (deadlockReporting) statusDeadlock=deadlockStatusIsReporting
       ! Enter loop for deadlock reporting.
       deadlock : do while (statusDeadlock /= deadlockStatusIsNotDeadlocked)
          ! Post a deadlocking message.
          if (statusDeadlock == deadlockStatusIsReporting) call displayIndent("Deadlock report follows")
          ! Iterate through all trees.
          currentTree => tree
          treesLoop: do while (associated(currentTree))
             ! Skip empty trees.
             if (associated(currentTree%nodeBase)) then
                ! Find the final time in this tree.
                node            => currentTree%nodeBase
                basic           => node       %basic   ()
                finalTimeInTree =  basic      %time    ()
                ! Report on current tree if deadlocked.
                if (statusDeadlock == deadlockStatusIsReporting) then
                   vMessage="tree "
                   vMessage=vMessage//currentTree%index
                   call displayIndent(vMessage)
                end if
                ! Point to the base of the tree.
                node => currentTree%nodeBase
                ! Get the basic component of the node.
                basic => node%basic()
                ! Tree walk loop: Walk to each node in the tree and consider whether or not to evolve it.
                treeWalker =mergerTreeWalkerAllNodes(currentTree,spanForest=.false.)
                nodesRemain=treeWalker%next(node)
                treeWalkLoop: do while (associated(node))
                   ! Get the basic component of the node.
                   basic => node%basic()
                   ! Count nodes in the tree.
                   nodesTotalCount=nodesTotalCount+1
                   ! Find the next node that we will process.
                   call treeWalker%setNode(node)
                   ! Store pointers to the parent and sibling of this node, to be possibly used when deciding the next node to
                   ! process. Initialize the next node pointer to null.
                   nodeParentNext  => node%parent
                   nodeSiblingNext => node%sibling
                   nodeNext        => null()
                   ! Evolve this node if it has a parent (or will transfer to another tree where it will have a parent), exists
                   ! before the output time, has no children (i.e. they've already all been processed), and either exists before
                   ! the final time in its tree, or exists precisely at that time and has some attached event yet to occur.
                   event       =>            node%event
                   hasParent   =  associated(node%parent)
                   treeLimited =  .true.
                   do while (associated(event).and.treeLimited)
                      ! Skip events which occur after the current evolution end time.
                      if (event%time <= timeEnd) then
                         ! Detect inter-tree events.
                         select type (event)
                         type is (nodeEventSubhaloPromotionInterTree)
                            hasParent  =.true.
                            treeLimited=.false.
                         type is (nodeEventBranchJumpInterTree      )
                            hasParent  =.true.
                            treeLimited=.false.
                         end select
                      end if
                      event => event%next
                   end do
                   evolveCondition: if (                                     &
                        &                     hasParent                      &
                        &               .and.                                &
                        &                .not.associated(node%firstChild  )  &
                        &               .and.                                &
                        &                (                                   &
                        &                  basic%time() <  timeEnd           &
                        &                 .or.                               &
                        &                  (                                 &
                        &                    .not.treeLimited                & ! For nodes that are not tree limited (i.e. have a node which
                        &                   .and.                            & ! will jump to another tree), allow them to evolve if the node
                        &                    basic%time() == timeEnd         & ! is at the end time also, since the jump may occur at that time.
                        &                  )                                 &
                        &                )                                   &
                        &               .and.                                &
                        &                (                                   &
                        &                   .not.treeLimited                 &
                        &                .or.                                &
                        &                    basic%time() <  finalTimeInTree &
                        &                .or.                                &
                        &                 (                                  &
                        &                  (                                 &
                        &                     associated(node%event      )   &
                        &                   .or.                             &
                        &                     associated(node%mergeTarget)   &
                        &                  )                                 &
                        &                  .and.                             &
                        &                    basic%time() <= finalTimeInTree &
                        &                 )                                  &
                        &                )                                   &
                        &              ) then
                      ! Flag that a node was evolved.
                      didEvolve=.true.
                      ! Flag that this node has not yet progressed in time, and that we have not yet determined which node we will
                      ! evolve next.
                      nodeProgressed=.false.
                      nextNodeFound =.false.
                      ! Update tree progress counter.
                      nodesEvolvedCount=nodesEvolvedCount+1
                      ! Dump the merger tree structure for later plotting.
                      if (self%dumpTreeStructure) call Merger_Tree_Dump(currentTree,[node%index()])
                      ! Evolve the node, handling interrupt events. We keep on evolving it until no interrupt is returned (in which case
                      ! the node has reached the requested end time) or the node no longer exists (e.g. if it was destroyed).
                      interrupted=.true.
                      do while (interrupted.and.associated(node))
                         interrupted=.false.
                         ! Find maximum allowed end time for this particular node.
                         if (statusDeadlock == deadlockStatusIsReporting) then
                            vMessage="node "
                            write (label,'(e12.6)') basic%time()
                            vMessage=vMessage//node%index()//" (current:target times = "//label
                            write (label,'(e12.6)') timeEnd
                            vMessage=vMessage//":"//label//")"
                            call displayIndent(vMessage)
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.true. ,nodeLock=nodeLock,lockType=lockType)
                            call displayUnindent("end node")
                            call self%deadlockAddNode(node,currentTree%index,nodeLock,lockType)
                         else
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.false.                                    )
                         end if
                         ! If this node is able to evolve by a finite amount, the tree is not deadlocked.
                         if (timeEndThisNode > basic%time()) then
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                            ! Record that the node will be progressed forward in time.
                            nodeProgressed=.true.
                            ! Update record of earliest time in the tree.
                            earliestTimeInTree=min(earliestTimeInTree,timeEndThisNode)
                            ! Evolve the node to the next interrupt event, or the end time.
                            call self%mergerTreeNodeEvolver_%evolve(currentTree,node,timeEndThisNode,interrupted,interruptProcedure,self%galacticStructureSolver_,systemClockMaximum,status)
                            if (present(status)) then
                               if (status /= errorStatusSuccess) return
                            end if
                         end if
                         ! Check for interrupt.
                         if (interrupted) then
                            ! If an interrupt occured call the specified procedure to handle it.
                            call interruptProcedure(node)
                            ! Something happened so the tree is not deadlocked.
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                         else
                            ! Call routine to handle end of timestep processing.
                            if (associated(timestepTask_).and.associated(node)) call timestepTask_(timestepSelf,currentTree,node,statusDeadlock)
                         end if
                      end do
                      ! If this halo has reached its parent halo, decide how to handle it.
                      if (associated(node)) then
                         if (associated(node%parent)) then
                            nodeParent  => node      %parent
                            basicParent => nodeParent%basic ()
                            if (basic%time() >= basicParent%time()) then
                               ! Parent halo has been reached. Check if the node is the primary (major) progenitor of the parent node.
                               select case (node%isPrimaryProgenitor())
                               case (.false.)
                                  ! It is not the major progenitor, so this could be a halo merger event unless the halo is already a
                                  ! satellite. Check for satellite status and, if it's not a satellite, process this halo merging
                                  ! event. Also record that the tree is not deadlocked, as we are changing the tree state.
                                  if (.not.node%isSatellite()) then
                                     statusDeadlock=deadlockStatusIsNotDeadlocked
                                     call self%mergerTreeNodeEvolver_%merge(node)
                                  end if
                               case (.true.)
                                  ! This is the major progenitor, so promote the node to its parent providing that the node has no
                                  ! siblings - this ensures that any siblings have already been evolved and become satellites of the
                                  ! parent halo. Also record that the tree is not deadlocked, as we are changing the tree state.
                                  if (.not.associated(node%sibling).and..not.associated(node%event)) then
                                     statusDeadlock=deadlockStatusIsNotDeadlocked
                                     call self%mergerTreeNodeEvolver_%promote(node)
                                     ! As this is a node promotion, we want to attempt to continue evolving the same node. Mark the
                                     ! parent as the next node to evolve, and flag that our next node has been identified.
                                     nodeNext      => nodeParent
                                     nextNodeFound =  .true.
                                  end if
                               end select
                            end if
                         end if
                      end if
                      ! Determine the next node to process.
                      if (.not.nextNodeFound) then
                         ! If the node still exists and advanced forward in time, attempt to evolve it, or its satellites,
                         ! again. Set the next node pointer back to the current node so that we will evolve it again, or (if
                         ! selected) to any satellites which may now be evolvable.
                         if (associated(node).and.nodeProgressed) then
                            if (associated(node%firstSatellite).and.self%backtrackToSatellites) then
                               nodeNext => node%firstSatellite
                            else
                               nodeNext => node
                            end if
                         end if
                      end if
                   else
                      if (statusDeadlock == deadlockStatusIsReporting) then
                         vMessage="node "
                         write (label,'(e12.6)') basic%time()
                         vMessage=vMessage//node%index()//" (current:target times = "//label
                         write (label,'(e12.6)') timeEnd
                         vMessage=vMessage//":"//label//")"
                         call displayIndent(vMessage)
                         call displayUnindent("end node")
                         ! Determine why this node could not be evolved. We check the "has child" condition first as it's the only
                         ! one that provides additional connection between nodes, so leads to the most informative deadlock graph.
                         if      (associated(node%firstChild)) then
                            call self%deadlockAddNode(node,currentTree%index,node%firstChild,var_str("has child"          ))
                         else if (.not.associated(node%parent)) then
                            call self%deadlockAddNode(node,currentTree%index,node           ,var_str("no parent"          ))
                         else if (basic%time() >= timeEnd) then
                            call self%deadlockAddNode(node,currentTree%index,node           ,var_str("in future of output"))
                         else
                            call self%deadlockAddNode(node,currentTree%index,node           ,var_str("in future of tree"  ))
                         end if
                      end if
                   end if evolveCondition
                   ! Step to the next node to consider.
                   if (associated(nodeNext)) then
                      ! The next node to consider has already been determined.
                      nodesRemain=.true.
                   else
                      ! The next node to consider has not yet been determined.
                      if (associated(node)) then
                         ! The current node still exists - simply walk to the next node in the tree using our tree walker.
                         nodesRemain=treeWalker%next(nodeNext)
                      else
                         ! The current node no longer exists. Attempt to used the sibling or parent of the node (saved from before
                         ! it was destroyed) as our next node. If either is available step down to the deepest satellite or child
                         ! of that node.
                         if (associated(nodeSiblingNext).or.associated(nodeParentNext)) then
                            if (associated(nodeSiblingNext)) then
                               nodeNext => nodeSiblingNext
                            else
                               nodeNext => nodeParentNext
                            end if
                            do while (associated(nodeNext%firstSatellite).or.associated(nodeNext%firstChild))
                               if (associated(nodeNext%firstSatellite)) then
                                  nodeNext => nodeNext%firstSatellite
                               else
                                  nodeNext => nodeNext%firstChild
                               end if
                            end do
                         else
                            ! The current node no longer exists, and had no sibling or parent - we having nothing left to process
                            ! in this tree so set the next node to null to force the loop to exit.
                            nodeNext => null()
                         end if
                      end if
                   end if
                   node => nodeNext
                end do treeWalkLoop
                ! Output tree progress information.
                if (treeWalkCount > int(treeWalkCountPreviousOutput*1.1d0)+1) then
                   if (displayVerbosity() >= verbosityLevel) then
                      !![
                      <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
                       <description>Internal file I/O in gfortran can be non-thread safe.</description>
                      </workaround>
                      !!]
#ifdef THREADSAFEIO
                      !$omp critical(gfortranInternalIO)
#endif
                      write (message,'(a,i9,a )') 'Evolving tree [',treeWalkCount,']'
#ifdef THREADSAFEIO
                      !$omp end critical(gfortranInternalIO)
#endif
                      call displayIndent(message,verbosityLevel)
#ifdef THREADSAFEIO
                      !$omp critical(gfortranInternalIO)
#endif
                      write (message,'(a,i9   )') 'Nodes in tree:         ',nodesTotalCount
#ifdef THREADSAFEIO
                      !$omp end critical(gfortranInternalIO)
#endif
                      call displayMessage(message,verbosityLevel)
#ifdef THREADSAFEIO
                      !$omp critical(gfortranInternalIO)
#endif
                      write (message,'(a,i9   )') 'Nodes evolved:         ',nodesEvolvedCount
#ifdef THREADSAFEIO
                      !$omp end critical(gfortranInternalIO)
#endif
                      call displayMessage(message,verbosityLevel)
#ifdef THREADSAFEIO
                      !$omp critical(gfortranInternalIO)
#endif
                      write (message,'(a,e10.4)') 'Earliest time in tree: ',earliestTimeInTree
#ifdef THREADSAFEIO
                      !$omp end critical(gfortranInternalIO)
#endif
                      call displayMessage(message,verbosityLevel)
                      call displayUnindent('done',verbosityLevel)
                      treeWalkCountPreviousOutput=treeWalkCount
                   end if
                end if
                ! Report on current tree if deadlocked.
                if (statusDeadlock == deadlockStatusIsReporting) call displayUnindent('end tree')
             end if
             ! Move to the next tree.
             currentTree => currentTree%nextTree
          end do treesLoop
          ! Perform any tree events.
          call standardTreeEventsPerform(tree,statusDeadlock)
          ! Check deadlocking.
          if (didEvolve .and. statusDeadlock /= deadlockStatusIsNotDeadlocked) then
             if (statusDeadlock == deadlockStatusIsReporting) then
                call displayUnindent("report done")
                call self%deadlockOutputTree(timeEnd)
                if (.not.deadlockReporting) then
                   call Galacticus_Error_Report('merger tree appears to be deadlocked (see preceding report) - check timestep criteria'//{introspection:location})
                else
                   return
                end if
             else
                ! Tree appears to be deadlocked. Check if it is suspendable.
                if (statusDeadlock == deadlockStatusIsSuspendable) then
                   ! Tree is suspendable, so do not attempt to process further, but simply return and flag it for suspension.
                   suspendTree  =.true.
                   return
                else
                   ! Tree is truly deadlocked. Switch to reporting mode and do one more pass through the tree.
                   statusDeadlock=deadlockStatusIsReporting
                end if
             end if
          else
             ! No evolution could occur, so the tree is not deadlocked.
             statusDeadlock=deadlockStatusIsNotDeadlocked
          end if
       end do deadlock
       ! Record tree evolution status.
       if (didEvolve) treeDidEvolve=.true.
    end do outerLoop
    return
  end subroutine standardEvolve

  recursive function standardTimeEvolveTo(self,node,timeEnd,timestepTask_,timestepSelf,report,nodeLock,lockType) result(evolveToTime)
    !!{
    Determine the time to which {\normalfont \ttfamily node} should be evolved.
    !!}
    use :: Display               , only : displayIndent                     , displayMessage        , displayUnindent, verbosityLevelInfo
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Error      , only : Galacticus_Error_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic                , nodeComponentSatellite, nodeEvent      , nodeEventBranchJumpInterTree, &
          &                               nodeEventSubhaloPromotionInterTree, treeEvent             , treeNode
    use :: Merger_Tree_Timesteps , only : timestepTask
    use :: String_Handling       , only : operator(//)
    implicit none
    class           (mergerTreeEvolverStandard    ), intent(inout)                    :: self
    double precision                                                                  :: evolveToTime
    type            (treeNode                     ), intent(inout)          , pointer :: node
    double precision                               , intent(in   )                    :: timeEnd
    type            (treeNode                     )                         , pointer :: nodeSatellite       , nodeSibling
    procedure       (timestepTask                 ), intent(  out)          , pointer :: timestepTask_
    class           (*                            ), intent(  out)          , pointer :: timestepSelf
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(  out), optional, pointer :: nodeLock
    type            (varying_string               ), intent(  out), optional          :: lockType
    procedure       (timestepTask                 )                         , pointer :: timestepTaskInternal
    class           (nodeComponentBasic           )                         , pointer :: basicParent         , basicSatellite    , &
         &                                                                               basicSibling        , basic
    class           (nodeComponentSatellite       )                         , pointer :: satelliteSatellite
    class           (nodeEvent                    )                         , pointer :: event
    class           (treeEvent                    )                         , pointer :: treeEvent_
    double precision                                                                  :: expansionFactor     , expansionTimescale, &
         &                                                                               hostTimeLimit       , time              , &
         &                                                                               timeEarliest        , evolveToTimeStep  , &
         &                                                                               hostTimeStep        , timeNode          , &
         &                                                                               timeSatellite
    character       (len=9                        )                                   :: timeFormatted
    type            (varying_string               ), save                             :: message
    !$omp threadprivate(message)

    ! Initially set to the global end time.
    evolveToTime=timeEnd
    if (report) call Evolve_To_Time_Report("start (target): ",evolveToTime)
    ! Initialize the lock node if present.
    if (present(nodeLock)) nodeLock => null()
    if (present(lockType)) lockType = 'null'
    ! Get the basic component of the node.
    basic    => node %basic()
    timeNode =  basic%time ()
    ! Ensure that the timestep doesn't exceed any event attached to the tree.
    treeEvent_ => node%hostTree%event
    do while (associated(treeEvent_))
       if (max(treeEvent_%time,timeNode) <= evolveToTime) then
          if (present(nodeLock)) nodeLock => node%hostTree%nodeBase
          if (present(lockType)) then
             lockType =  "tree event ("
             lockType=lockType//treeEvent_%ID//")"
          end if
          evolveToTime=max(treeEvent_%time,timeNode)
       end if
       if (report) then
          message="tree event ("
          message=message//treeEvent_%ID//"): "
          call Evolve_To_Time_Report(char(message),evolveToTime,node%hostTree%nodeBase%index())
       end if
       treeEvent_ => treeEvent_%next
    end do
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Ensure that the timestep doesn't exceed any event attached to the node.
    event => node%event
    do while (associated(event))
       if (max(event%time,timeNode) <= evolveToTime) then
          if (present(nodeLock)) nodeLock => event%node
          if (present(lockType)) then
             lockType =  "event ("
             select type (event)
             type is (nodeEventSubhaloPromotionInterTree)
                lockType=lockType//event%splitForestUniqueID//":"//event%pairedNodeID
             type is (nodeEventBranchJumpInterTree      )
                lockType=lockType//event%splitForestUniqueID//":"//event%pairedNodeID
             class default
                lockType=lockType//event%ID
             end select
             lockType=lockType//")"
          end if
          evolveToTime=max(event%time,timeNode)
          timestepTask_ => standardNodeEventsPerform
       end if
       if (report) then
          message="event ("
          message=message//event%ID//"): "
          call Evolve_To_Time_Report(char(message),evolveToTime,event%node%index())
       end if
       event => event%next
    end do
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Also ensure that the timestep taken does not exceed the allowed timestep for this specific node.
    if (report) call displayIndent("timestepping criteria")
    evolveToTimeStep=self%mergerTreeEvolveTimestep_%timeEvolveTo(evolveToTime,node,timestepTaskInternal,timestepSelf,report,nodeLock,lockType)
    if (evolveToTimeStep <= evolveToTime) then
       evolveToTime  =  evolveToTimeStep
       timestepTask_ => timestepTaskInternal
    else
       timestepTask_ => null()
       timestepSelf  => null()
    end if
    if (report) call displayUnindent("done")
    if (evolveToTime == timeNode) return
    ! Ensure that this node is not evolved beyond the time of any of its current satellites.
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       basicSatellite => nodeSatellite %basic()
       timeSatellite  =  basicSatellite%time ()
       if (max(timeSatellite,timeNode) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) lockType =  "hosted satellite"
          evolveToTime=max(timeSatellite,timeNode)
       end if
       if (report) call Evolve_To_Time_Report("hosted satellite: ",evolveToTime,nodeSatellite%index())
       if (evolveToTime == timeNode) exit
       nodeSatellite => nodeSatellite%sibling
    end do
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Also ensure that this node is not evolved beyond the time at which any of its mergees merge. In some cases, the node may
    ! already be in the future of a mergee. In such cases, simply freeze it at the current time.
    nodeSatellite => node%firstMergee
    do while (associated(nodeSatellite))
       satelliteSatellite => nodeSatellite%satellite()
       if (max(satelliteSatellite%timeOfMerging(),timeNode) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) then
             write (timeFormatted,'(f7.4)') max(satelliteSatellite%timeOfMerging(),timeNode)
             lockType =  "mergee ("//trim(timeFormatted)//")"
          end if
          evolveToTime=max(satelliteSatellite%timeOfMerging(),timeNode)
       end if
       if (report) call Evolve_To_Time_Report("mergee limit: ",evolveToTime,nodeSatellite%index())
       nodeSatellite => nodeSatellite%siblingMergee
    end do
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Also ensure that a primary progenitor does not evolve in advance of siblings. This is important since we can not promote a
    ! primary progenitor into its parent until all siblings have become satellites in that parent.
    if (node%isPrimaryProgenitor()) then
       nodeSibling => node%sibling
       do while (associated(nodeSibling))
          basicSibling => nodeSibling%basic()
          if (max(timeNode,basicSibling%time()) < evolveToTime) then
             if (present(nodeLock)) nodeLock => nodeSibling
             if (present(lockType)) lockType =  "sibling"
             evolveToTime=max(timeNode,basicSibling%time())
          end if
          if (report) call Evolve_To_Time_Report("sibling: ",evolveToTime,nodeSibling%index())
          nodeSibling => nodeSibling%sibling
       end do
    end if
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Limit time based on satellite status.
    select case (node%isSatellite())
    case (.false.)
       ! Limit to the time of its parent node if this node is not a satellite.
       if (associated(node%parent)) then
          basicParent => node%parent%basic()
          if (basicParent%time() < evolveToTime) then
             if (present(nodeLock)) nodeLock => node%parent
             if (present(lockType)) lockType =  "promotion"
             evolveToTime=basicParent%time()
          end if
       end if
       if (report) call Evolve_To_Time_Report("promotion limit: ",evolveToTime)
    case (.true.)
       ! Do not let satellite evolve too far beyond parent.
       if (associated(node%parent%parent)) then
          ! The host halo has a parent, so use the host halo time to limit satellite evolution.
          basicParent => node%parent%basic()
          time=basicParent%time()
       else
          ! The host halo has no parent. The satellite must therefore be evolving to some event (e.g. a merger). We have to allow
          ! it to evolve ahead of the host halo in this case to avoid deadlocks.
          basicParent => node%parent%basic()
          time=max(basicParent%time(),timeNode)
       end if
       ! Check if the host has a child.
       select case (associated(node%parent%firstChild))
       case (.true. )
          ! Host still has a child - do not let the satellite evolve beyond the host.
          hostTimeLimit=max(time,timeNode)
          ! Check for any merge targets directed at this node.
          nodeSatellite => node%firstMergee
          timeEarliest=huge(1.0d0)
          do while (associated(nodeSatellite))
             satelliteSatellite => nodeSatellite%satellite()
             if (nodeSatellite%isSatellite().and.nodeSatellite%parent%isProgenitorOf(node%parent)) &
                  & timeEarliest=min(timeEarliest,satelliteSatellite%timeOfMerging())
             nodeSatellite => nodeSatellite%siblingMergee
          end do
          if (timeEarliest < huge(1.0d0)) hostTimeLimit=max(hostTimeLimit,timeEarliest)
       case (.false.)
          ! Find current expansion timescale.
          if (self%timestepHostRelative > 0.0d0) then
             expansionFactor   =      self%cosmologyFunctions_%expansionFactor(time           )
             expansionTimescale=1.0d0/self%cosmologyFunctions_%expansionRate  (expansionFactor)
             hostTimeStep      =min(self%timestepHostRelative*expansionTimescale,self%timestepHostAbsolute)
          else
             ! Avoid use of expansion timescale if host absolute timestep is non-positive. This allows static universe cases to be handled.
             hostTimeStep      =                                                 self%timestepHostAbsolute
          end if
          hostTimeLimit=max(time+hostTimeStep,timeNode)
          ! Check if this criterion will actually limit the evolution time.
          if (hostTimeLimit < evolveToTime) then
             ! Satellite evolution will be limited by being required to not advance too far ahead of the host halo. If the
             ! timestep we can advance is less than a specified fraction of what is possible, then skip this for now.
             if (hostTimeLimit-timeNode < self%fractionTimestepSatelliteMinimum*hostTimeStep) then
                hostTimeLimit=timeNode
             end if
          end if
       end select
       ! Limit to this time.
       if (hostTimeLimit < evolveToTime) then
          if (present(nodeLock)) nodeLock => node%parent
          if (present(lockType)) lockType =  "satellite in host"
          evolveToTime=hostTimeLimit
       end if
       if (report) call Evolve_To_Time_Report("satellite in host limit: ",evolveToTime,node%parent%index())
    end select
    ! Check that end time exceeds current time.
    if (evolveToTime < timeNode) then
       message='end time ('
       write (timeFormatted,'(f7.4)') evolveToTime
       message=message//trim(timeFormatted)//' Gyr) is before current time ('
       write (timeFormatted,'(f7.4)') timeNode
       message=message//trim(timeFormatted)//' Gyr) of node '
       message=message//node%index()
       message=message//' (time difference is '
       write (timeFormatted,'(e8.2)') timeNode-evolveToTime
       message=message//trim(timeFormatted)
       if (.not.self%mergerTreeNodeEvolver_%isAccurate(timeNode,evolveToTime)) then
          ! End time is well before current time. This is an error. Call ourself with reporting switched on to generate a report
          ! on the time limits.
          message=message//' Gyr)'
          if (.not.report) time=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.true.)
          call Galacticus_Error_Report(message//{introspection:location})
       else
          ! End time is before current time, but only by a small amount, simply reset the current time to the end time.
          message=message//' Gyr) - this should happen infrequently'
          call displayMessage(message,verbosityLevelInfo)
          call basic%timeSet(evolveToTime)
       end if
    end if
    return
  end function standardTimeEvolveTo

  subroutine standardDeadlockAddNode(self,node,treeIndex,nodeLock,lockType)
    !!{
    Add a node to the deadlocked nodes list.
    !!}
    implicit none
    class  (mergerTreeEvolverStandard), intent(inout)          :: self
    type   (treeNode                 ), intent(in   ), target  :: nodeLock        , node
    integer(kind=kind_int8           ), intent(in   )          :: treeIndex
    type   (varying_string           ), intent(in   )          :: lockType
    type   (deadlockList             )               , pointer :: deadlockThisNode

    ! Add a node to the deadlock linked list.
    if (associated(self%deadlockHeadNode)) then
       deadlockThisNode => self%deadlockHeadNode
       do while (associated(deadlockThisNode%next))
          deadlockThisNode => deadlockThisNode%next
       end do
       allocate(deadlockThisNode%next)
       deadlockThisNode => deadlockThisNode%next
    else
       allocate(self%deadlockHeadNode)
       deadlockThisNode => self%deadlockHeadNode
    end if
    ! Set properties.
    deadlockThisNode%node      => node
    deadlockThisNode%treeIndex =  treeIndex
    deadlockThisNode%nodeLock  => nodeLock
    deadlockThisNode%lockType  =  lockType
    return
  end subroutine standardDeadlockAddNode

  subroutine standardDeadlockOutputTree(self,timeEnd)
    !!{
    Output the deadlocked nodes in {\normalfont \ttfamily dot} format.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: String_Handling , only : operator(//)
    implicit none
    class           (mergerTreeEvolverStandard), intent(inout) :: self
    double precision                           , intent(in   ) :: timeEnd
    type            (deadlockList             ), pointer       :: nodeLock        , testNode  , &
         &                                                        node
    class           (nodeComponentBasic       ), pointer       :: basic
    type            (treeNode                 ), pointer       :: nodeParent
    logical                                                    :: foundLockNode
    integer                                                    :: treeUnit
    integer         (kind=kind_int8           )                :: uniqueID
    logical                                                    :: inCycle         , nodesAdded
    character       (len=20                   )                :: color           , style
    type            (varying_string           )                :: deadlockFileName

    ! If no deadlock list exists, simply return.
    if (.not.associated(self%deadlockHeadNode)) return
    ! Begin tree.
    deadlockFileName=var_str('galacticusDeadlockTree_')//self%deadlockHeadNode%node%hostTree%nodeBase%uniqueID()//'.gv'
    open(newUnit=treeUnit,file=char(deadlockFileName),status='unknown',form='formatted')
    write (treeUnit,*) 'digraph Tree {'
    ! Find any nodes that cause a lock but which are not in our list.
    nodesAdded=.true.
    do while (nodesAdded)
       nodesAdded=.false.
       node => self%deadlockHeadNode
       do while (associated(node))
          if (associated(node%nodeLock)) then
             testNode => self%deadlockHeadNode
             foundLockNode=.false.
             do while (associated(testNode).and..not.foundLockNode)
                foundLockNode=(associated(node%nodeLock,testNode%node))
                testNode => testNode%next
             end do
             if (.not.foundLockNode) then
                nodesAdded =  .true.
                testNode   => self%deadlockHeadNode
                do while (associated(testNode%next))
                   testNode => testNode%next
                end do
                allocate(testNode%next)
                testNode => testNode%next
                ! Find root node.
                nodeParent => node%nodeLock
                do while (associated(nodeParent%parent))
                   nodeParent => nodeParent%parent
                end do
                ! Set properties.
                testNode%node      => node%nodeLock
                testNode%treeIndex =  nodeParent%index()
                testNode%nodeLock  => null()
                testNode%lockType  =  "unknown"
                basic => node%nodeLock%basic()
                if (associated(node%nodeLock%firstChild)) then
                   testNode%lockType = "child"
                   nodeLock        => self%deadlockHeadNode
                   do while (associated(nodeLock))
                      if (associated(node%nodeLock%firstChild,nodeLock%node)) then
                         testNode%nodeLock => node%nodeLock%firstChild
                         exit
                      end if
                      nodeLock => nodeLock%next
                   end do
                end if
                if (basic%time() >= timeEnd) testNode%lockType = "end time"
             end if
          end if
          node => node%next
       end do
    end do
    ! Iterate over all nodes visited.
    node => self%deadlockHeadNode
    do while (associated(node))
       ! Detect cycles.
       inCycle=.false.
       uniqueID=node%node%uniqueID()
       testNode => node%next
       do while (associated(testNode))
          if (testNode%node%uniqueID() == uniqueID) then
             inCycle=.true.
             exit
          end if
          testNode => testNode%next
       end do
       ! Output node.
       basic => node%node%basic()
       if (inCycle) then
          color="green"
          style="filled"
       else
          color="black"
          style="solid"
       end if
       write (treeUnit,'(a,i16.16,a,a,a,a,a,i16.16,a,i16.16,a,i16.16,a,f13.10,a,a,a)') '"',node%node%uniqueID(),'" [shape=circle, color=',trim(color),', style=',trim(style),' label="',node%node%uniqueid(),':',node%node%index(),'\ntree: ',node%treeIndex,'\ntime: ',basic%time(),'\n',char(node%lockType),'"];'
       if (associated(node%nodeLock)) write (treeUnit,'(a,i16.16,a,i16.16,a)') '"',node%node%uniqueID(),'" -> "',node%nodeLock%uniqueID(),'"'
       node => node%next
    end do
    ! Close the tree.
    write (treeUnit,*) '}'
    close(treeUnit)
    ! Clean up the deadlock node list.
    node => self%deadlockHeadNode
    do while (associated(node))
       testNode => node%next
       nullify(node%node    )
       nullify(node%nodeLock)
       deallocate(node)
       node => testNode
    end do
    self%deadlockHeadNode => null()
    return
  end subroutine standardDeadlockOutputTree

  subroutine standardNodeEventsPerform(self,tree,node,statusDeadlock)
    !!{
    Perform any events associated with {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, nodeComponentBasic, nodeEvent, treeNode
    implicit none
    class           (*                 ), intent(inout)          :: self
    type            (mergerTree        ), intent(in   )          :: tree
    type            (treeNode          ), intent(inout), pointer :: node
    integer                             , intent(inout)          :: statusDeadlock
    class           (nodeEvent         )               , pointer :: eventLast            , eventNext        , &
         &                                                          event
    class           (nodeComponentBasic)               , pointer :: basic
    double precision                                             :: timeNode             , timeEventEarliest
    logical                                                      :: mergerTreeEvolverDone
    !$GLC attributes unused :: self, tree

    ! Get the current time.
    basic    => node %basic()
    timeNode =  basic%time ()
    ! Find the current earliest event.
    event => node%event
    timeEventEarliest=huge(1.0d0)
    do while (associated(event))
       timeEventEarliest=min(timeEventEarliest,event%time)
       event => event%next
    end do
    ! Get the first event.
    event     => node%event
    eventLast => node%event
    ! Iterate over all events.
    do while (associated(event))
       ! Process the event if it occurs at the present time.
       if (event%time <= timeNode .and. event%time == timeEventEarliest .and. associated(event%task)) then
          mergerTreeEvolverDone=event%task(node,statusDeadlock)
          ! If the node is no longer associated, simply exit (as any events associated with it must have been processed already).
          if (.not.associated(node)) exit
          ! Move to the next event.
          if (mergerTreeEvolverDone) then
             ! The mergerTreeEvolver was performed successfully, so remove it and move to the next event.
             if (associated(event,node%event)) then
                node     %event => event%next
                eventLast       => node %event
             else
                eventLast%next  => event%next
             end if
             eventNext => event%next
             deallocate(event)
             event => eventNext
          else
             ! The mergerTreeEvolver was not performed, so simply move to the next event.
             eventLast => event
             event     => event%next
          end if
       else
          eventLast => event
          event     => event%next
       end if
    end do
    return
  end subroutine standardNodeEventsPerform

  subroutine standardTreeEventsPerform(tree,statusDeadlock)
    !!{
    Perform any events associated with {\normalfont \ttfamily tree}.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, treeEvent
    implicit none
    type            (mergerTree), intent(inout), target  :: tree
    integer                     , intent(inout)          :: statusDeadlock
    type            (treeEvent )               , pointer :: eventLast            , eventNext, event
    double precision                                     :: treeTimeEarliest
    logical                                              :: mergerTreeEvolverDone

    ! Find the earliest time in the tree.
    treeTimeEarliest=tree%earliestTime()
    ! Get the first event.
    event     => tree%event
    eventLast => tree%event
    ! Iterate over all events.
    do while (associated(event))
       ! Process the event if it occurs at the present time.
       if (event%time <= treeTimeEarliest .and. associated(event%task)) then
          mergerTreeEvolverDone=event%task(tree,statusDeadlock)
          ! Move to the next event.
          if (mergerTreeEvolverDone) then
             ! The mergerTreeEvolver was performed successfully, so remove it and move to the next event.
             if (associated(event,tree%event)) then
                tree%event => event%next
                eventLast      => tree %event
             else
                eventLast%next => event%next
             end if
             eventNext => event%next
             deallocate(event)
             event => eventNext
          else
             ! The mergerTreeEvolver was not performed, so simply move to the next event.
             event => event%next
          end if
       else
          eventLast => event
          event => event%next
       end if
    end do
    return
  end subroutine standardTreeEventsPerform
