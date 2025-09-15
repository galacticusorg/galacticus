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
  Implements the standard class for evolving merger trees.
  !!}

  use :: Cosmology_Functions         , only : cosmologyFunctions           , cosmologyFunctionsClass
  use :: Galactic_Structure_Solvers  , only : galacticStructureSolver      , galacticStructureSolverClass
  use :: Galacticus_Nodes            , only : treeNode
  use :: Kind_Numbers                , only : kind_int8
  use :: Merger_Tree_Evolve_Profilers, only : mergerTreeEvolveProfilerClass
  use :: Merger_Tree_Initialization  , only : mergerTreeInitializorClass
  use :: Meta_Tree_Compute_Times     , only : metaTreeProcessingTimeClass
  use :: Merger_Tree_Timesteps       , only : mergerTreeEvolveTimestep  , mergerTreeEvolveTimestepClass
  use :: Merger_Trees_Evolve_Node    , only : mergerTreeNodeEvolver     , mergerTreeNodeEvolverClass

  ! Structure used to store list of nodes for deadlock reporting.
  type :: deadlockList
     type   (deadlockList  ), pointer :: next      => null()
     type   (treeNode      ), pointer :: nodeLock  => null(), node => null()
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
     Implementation of the standard merger tree evolver.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (mergerTreeEvolveTimestepClass), pointer :: mergerTreeEvolveTimestep_        => null()
     class           (galacticStructureSolverClass ), pointer :: galacticStructureSolver_         => null()
     class           (mergerTreeNodeEvolverClass   ), pointer :: mergerTreeNodeEvolver_           => null()
     class           (mergerTreeInitializorClass   ), pointer :: mergerTreeInitializor_           => null()
     class           (mergerTreeEvolveProfilerClass), pointer :: mergerTreeEvolveProfiler_        => null()
     class           (metaTreeProcessingTimeClass  ), pointer :: metaTreeProcessingTime_          => null()
     logical                                                  :: allTreesExistAtFinalTime                  , dumpTreeStructure    , &
          &                                                      backtrackToSatellites                     , profileSteps
     double precision                                         :: timestepHostAbsolute                      , timestepHostRelative , &
          &                                                      fractionTimestepSatelliteMinimum          , timeHostPrevious     , &
          &                                                      timeStepHost
     type            (deadlockList                 ), pointer :: deadlockHeadNode                 => null()
   contains
     !![
     <methods>
       <method method="initializeTree"     description="Initialize the tree(s)."                      />
       <method method="nodeIsEvolvable"    description="Determine if a node is evolvable."            />
       <method method="timeEvolveTo"       description="Find the time to which a node can be evolved."/>
       <method method="deadlockAddNode"    description="Add a node to the deadlock list."             />
       <method method="deadlockOutputTree" description="Output a description of a deadlocked tree."   />
     </methods>
     !!]
     final     ::                       standardDestructor
     procedure :: evolve             => standardEvolve
     procedure :: initializeTree     => standardInitializeTree
     procedure :: nodeIsEvolvable    => standardNodeIsEvolvable
     procedure :: timeEvolveTo       => standardTimeEvolveTo
     procedure :: deadlockAddNode    => standardDeadlockAddNode
     procedure :: deadlockOutputTree => standardDeadlockOutputTree
  end type mergerTreeEvolverStandard

  interface mergerTreeEvolverStandard
     !!{
     Constructors for the \refClass{mergerTreeEvolverStandard} merger tree evolver.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeEvolverStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolverStandard} merger tree evolver class which takes a parameter set as input.
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
    class           (mergerTreeEvolveProfilerClass), pointer       :: mergerTreeEvolveProfiler_
    class           (metaTreeProcessingTimeClass  ), pointer       :: metaTreeProcessingTime_
    logical                                                        :: allTreesExistAtFinalTime        , dumpTreeStructure         , &
         &                                                            backtrackToSatellites           , profileSteps
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
    <inputParameter>
      <name>profileSteps</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not to profile the ODE evolver.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="mergerTreeEvolveTimestep" name="mergerTreeEvolveTimestep_" source="parameters"/>
    <objectBuilder class="galacticStructureSolver"  name="galacticStructureSolver_"  source="parameters"/>
    <objectBuilder class="mergerTreeNodeEvolver"    name="mergerTreeNodeEvolver_"    source="parameters"/>
    <objectBuilder class="mergerTreeInitializor"    name="mergerTreeInitializor_"    source="parameters"/>
    <objectBuilder class="mergerTreeEvolveProfiler" name="mergerTreeEvolveProfiler_" source="parameters"/>
    <objectBuilder class="metaTreeProcessingTime"   name="metaTreeProcessingTime_"   source="parameters"/>
    !!]
    self=mergerTreeEvolverStandard(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,profileSteps,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,galacticStructureSolver_,mergerTreeEvolveProfiler_,metaTreeProcessingTime_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="mergerTreeEvolveTimestep_"/>
    <objectDestructor name="mergerTreeNodeEvolver_"   />
    <objectDestructor name="galacticStructureSolver_" />
    <objectDestructor name="mergerTreeInitializor_"   />
    <objectDestructor name="mergerTreeEvolveProfiler_"/>
    <objectDestructor name="metaTreeProcessingTime_"  />
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,profileSteps,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,galacticStructureSolver_,mergerTreeEvolveProfiler_,metaTreeProcessingTime_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolverStandard} merger tree evolver class.
    !!}
    implicit none
    type            (mergerTreeEvolverStandard    )                        :: self
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass), intent(in   ), target :: mergerTreeEvolveTimestep_
    class           (galacticStructureSolverClass ), intent(in   ), target :: galacticStructureSolver_
    class           (mergerTreeNodeEvolverClass   ), intent(in   ), target :: mergerTreeNodeEvolver_
    class           (mergerTreeInitializorClass   ), intent(in   ), target :: mergerTreeInitializor_
    class           (mergerTreeEvolveProfilerClass), intent(in   ), target :: mergerTreeEvolveProfiler_
    class           (metaTreeProcessingTimeClass  ), intent(in   ), target :: metaTreeProcessingTime_
    logical                                        , intent(in   )         :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                    backtrackToSatellites           , profileSteps
    double precision                               , intent(in   )         :: timestepHostRelative            , timestepHostAbsolute, &
         &                                                                    fractionTimestepSatelliteMinimum
    !![
    <constructorAssign variables="allTreesExistAtFinalTime, dumpTreeStructure, timestepHostRelative, timestepHostAbsolute, fractionTimestepSatelliteMinimum, backtrackToSatellites, profileSteps, *cosmologyFunctions_, *mergerTreeNodeEvolver_, *mergerTreeEvolveTimestep_, *mergerTreeInitializor_, *galacticStructureSolver_, *mergerTreeEvolveProfiler_, *metaTreeProcessingTime_"/>
    !!]

    self%deadlockHeadNode =>  null(     )
    self%timeHostPrevious =  -huge(0.0d0)
    self%timeStepHost     =  -huge(0.0d0)
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolverStandard} merger tree evolver class.
    !!}
    implicit none
    type(mergerTreeEvolverStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%mergerTreeEvolveTimestep_"/>
    <objectDestructor name="self%galacticStructureSolver_" />
    <objectDestructor name="self%mergerTreeNodeEvolver_"   />
    <objectDestructor name="self%mergerTreeInitializor_"   />
    <objectDestructor name="self%mergerTreeEvolveProfiler_"/>
    <objectDestructor name="self%metaTreeProcessingTime_"  />
    !!]
    return
  end subroutine standardDestructor

  subroutine standardEvolve(self,tree,timeEnd,treeDidEvolve,suspendTree,deadlockReporting,systemClockMaximum,initializationLock,status)
    !!{
    Evolves all properties of a merger tree to the specified time.
    !!}
    use    :: Display                            , only : displayIndent                , displayMessage                    , displayUnindent              , displayVerbosity           , &
         &                                                enumerationVerbosityLevelType, verbosityLevelWarn
    use    :: Error                              , only : Error_Report                 , errorStatusSuccess
    use    :: Galacticus_Nodes                   , only : interruptTask                , mergerTree                        , nodeComponentBasic
    use    :: Merger_Tree_Timesteps              , only : timestepTask
    use    :: Merger_Tree_Walkers                , only : mergerTreeWalkerAllNodes
    use    :: Merger_Trees_Dump                  , only : Merger_Tree_Dump
    use    :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsDeadlocked   , deadlockStatusIsNotDeadlocked     , deadlockStatusIsReporting, deadlockStatusIsSuspendable, &
         &                                                enumerationDeadlockStatusType
    use :: Numerical_Constants_Prefixes          , only : siFormat
    !$ use :: OMP_Lib                            , only : omp_lock_kind
    use    :: Locks                              , only : ompLockClass
    use    :: String_Handling                    , only : operator(//)
    implicit none
    class           (mergerTreeEvolverStandard    )                    , intent(inout) :: self
    integer                                        , optional          , intent(  out) :: status
    type            (mergerTree                   )           , target , intent(inout) :: tree
    double precision                                                   , intent(in   ) :: timeEnd
    logical                                                            , intent(  out) :: treeDidEvolve                                  , suspendTree
    logical                                                            , intent(in   ) :: deadlockReporting
    integer         (kind_int8                    ), optional          , intent(in   ) :: systemClockMaximum
    integer         (omp_lock_kind                ), optional          , intent(inout) :: initializationLock
    type            (treeNode                     )           , pointer                :: nodeLock                                       , nodeNext         , &
         &                                                                                nodeParent                                     , node             , &
         &                                                                                nodeSiblingNext                                , nodeParentNext
    double precision                               , parameter                         :: largeTime                  =1.0d10
    procedure       (interruptTask                )           , pointer                :: interruptProcedure
    procedure       (timestepTask                 )           , pointer                :: timestepTask_
    class           (*                            )           , pointer                :: timestepSelf
    type            (enumerationVerbosityLevelType), parameter                         :: verbosityLevel             =verbosityLevelWarn
    class           (nodeComponentBasic           )           , pointer                :: basic                                          , basicParent
    type            (mergerTree                   )           , pointer                :: currentTree
    type            (ompLockClass                 )                                    :: lockTree
    type            (mergerTreeWalkerAllNodes     )                                    :: treeWalker
    type            (enumerationDeadlockStatusType)                                    :: statusDeadlock
    integer                                                                            :: treeWalkCountPreviousOutput                    , nodesEvolvedCount, &
         &                                                                                nodesTotalCount                                , treeWalkCount
    double precision                                                                   :: earliestTimeInTree                             , timeEndThisNode  , &
         &                                                                                finalTimeInTree                                , timeRemaining
    character       (len=24                       )                                    :: label
    character       (len=35                       )                                    :: message
    type            (varying_string               )                                    :: lockType                                       , vMessage
    logical                                                                            :: anyTreeExistsAtOutputTime                      , hasIntertreeEvent, &
         &                                                                                nodeProgressed                                 , nextNodeFound    , &
         &                                                                                didEvolve                                      , interrupted      , &
         &                                                                                nodesRemain

    ! Initialize trees.
    suspendTree               =  .false.
    anyTreeExistsAtOutputTime =  .false.
    treeDidEvolve             =  .false.
    call self%initializeTree(tree,timeEnd,treeDidEvolve,anyTreeExistsAtOutputTime,hasInterTreeEvent,initializationLock)
    if (.not.anyTreeExistsAtOutputTime.and..not.hasInterTreeEvent) then
       ! Mark the tree as evolved here, as the only reason that we did not evolve it was the given time target.
       treeDidEvolve=.true.
       return
    end if
    ! Initialize remaining time calculations.
    timeRemaining=self%metaTreeProcessingTime_%timeRemaining(tree,timeEnd)
    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible.
    didEvolve                  =.true.
    treeWalkCount              =0
    treeWalkCountPreviousOutput=0
    lockTree                   =ompLockClass()
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
                   ! Determine if the node is evolvable.
                   evolveCondition: if (self%nodeIsEvolvable(node,timeEnd,finalTimeInTree)) then
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
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%cosmologyFunctions_,self%mergerTreeEvolveTimestep_,self%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.true. ,nodeLock=nodeLock,lockType=lockType)
                            call displayUnindent("end node")
                            call self%deadlockAddNode(node,currentTree%index,nodeLock,lockType)
                         else if (self%profileSteps) then
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%cosmologyFunctions_,self%mergerTreeEvolveTimestep_,self%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.false.                  ,lockType=lockType)
                            call self%mergerTreeEvolveProfiler_%stepDescriptor(lockType)
                         else
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%cosmologyFunctions_,self%mergerTreeEvolveTimestep_,self%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.false.                                    )
                         end if
                         ! If this node is able to evolve by a finite amount, the tree is not deadlocked.
                         if (timeEndThisNode > basic%time()) then
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                            ! Record that the node will be progressed forward in time.
                            nodeProgressed=.true.
                            ! Update record of earliest time in the tree.
                            earliestTimeInTree=min(earliestTimeInTree,timeEndThisNode)
                            ! Evolve the node to the next interrupt event, or the end time.
                            call self%mergerTreeNodeEvolver_%evolve(currentTree,node,timeEndThisNode,interrupted,interruptProcedure,self%galacticStructureSolver_,lockTree,systemClockMaximum,status)
                            if (present(status)) then
                               if (status /= errorStatusSuccess) return
                            end if
                         end if
                         ! Check for interrupt.
                         if     (                                                         &
                              &   interrupted                                             & ! An interrupt occured.
                              &  .and.                                                    &
                              &   (                                                       &
                              &    basic%time() < timeEndThisNode                         & ! The end of the timestep was not reached.
                              &    .or.                                                   &
                              &     .not.(associated(timestepTask_).and.associated(node)) & ! No end of timestep task is possible.
                              &   )                                                       &
                              & ) then
                              ! If an interrupt occurred call the specified procedure to handle it.
                            call interruptProcedure(node,timeEnd)
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
                ! Estimate remaining time to process the tree.
                timeRemaining=self%metaTreeProcessingTime_%timeRemaining(tree,timeEnd)
                if (timeRemaining > 0.0d0) &
                     & call displayMessage("Estimated time remaining to process tree: "//trim(adjustl(siFormat(timeRemaining,'f7.2,1x')))//"s")
                ! Output tree progress information.
                if (treeWalkCount > int(treeWalkCountPreviousOutput*1.1d0)+1) then
                   if (displayVerbosity() >= verbosityLevel) then
                      write (message,'(a,i9,a )') 'Evolving tree [',treeWalkCount,']'
                      call displayIndent(message,verbosityLevel)
                      write (message,'(a,i9   )') 'Nodes in tree:         ',nodesTotalCount
                      call displayMessage(message,verbosityLevel)
                      write (message,'(a,i9   )') 'Nodes evolved:         ',nodesEvolvedCount
                      call displayMessage(message,verbosityLevel)
                      write (message,'(a,e10.4)') 'Earliest time in tree: ',earliestTimeInTree
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
                   call Error_Report('merger tree appears to be deadlocked (see preceding report) - check timestep criteria'//{introspection:location})
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

  subroutine standardInitializeTree(self,tree,timeEnd,treeDidEvolve,anyTreeExistsAtOutputTime,hasInterTreeEvent,initializationLock)
    !!{
    Initialize trees prior to evolution.
    !!}
    use    :: Display            , only : displayBlue             , displayYellow , displayGreen                , displayBold                       , &
         &                                displayReset
    use    :: Galacticus_Nodes   , only : nodeComponentBasic      , nodeEvent     , nodeEventBranchJumpInterTree, nodeEventSubhaloPromotionInterTree
    use    :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    use    :: String_Handling    , only : operator(//)            , stringXMLFormat
    !$ use :: OMP_Lib            , only : OMP_Set_Lock            , OMP_Unset_Lock , omp_lock_kind
    implicit none
    class           (mergerTreeEvolverStandard)           , intent(inout) :: self
    type            (mergerTree               ), target   , intent(inout) :: tree
    double precision                                      , intent(in   ) :: timeEnd
    logical                                               , intent(inout) :: treeDidEvolve
    logical                                               , intent(inout) :: anyTreeExistsAtOutputTime
    logical                                               , intent(  out) :: hasInterTreeEvent
    integer         (omp_lock_kind            ), optional , intent(inout) :: initializationLock
    double precision                           , parameter                :: timeTolerance            =1.0d-5
    type            (mergerTree               ), pointer                  :: currentTree
    type            (treeNode                 ), pointer                  :: node
    class           (nodeComponentBasic       ), pointer                  :: basicBase
    class           (nodeEvent                ), pointer                  :: event
    type            (varying_string           )                           :: vMessage
    type            (mergerTreeWalkerAllNodes )                           :: treeWalker
    character       (len=24                   )                           :: label

    currentTree => tree
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
                   vMessage=vMessage//displayGreen()//' HELP:'//displayReset()//' If you expect that not all trees will exist at the latest requested output'//char(10)
                   vMessage=vMessage//                                         '    time (this can happen when using trees extracted from N-body simulations for'//char(10)
                   vMessage=vMessage//                                         '    example) set the highlighted option in your input parameter file as shown below:'//char(10)//char(10)
                   vMessage=vMessage//stringXMLFormat('<mergerTreeEvolver value="doop">**B<allTreesExistAtFinalTime value="false" />**C</mergerTreeEvolver>',indentInitial=6)//char(10)
                  call Error_Report(vMessage//{introspection:location})
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
                      call Error_Report(vMessage//{introspection:location})
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
    hasInterTreeEvent=.false.
    if (.not.anyTreeExistsAtOutputTime) then
       ! Walk over all trees in the forest.
       treeWalker       =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
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
       if (.not.hasInterTreeEvent) &
            & treeDidEvolve=.true.   ! Mark the tree as evolved here, as the only reason that we did not evolve it was the given time target.
    end if
    return
  end subroutine standardInitializeTree
  
  logical function standardNodeIsEvolvable(self,node,timeEnd,finalTimeInTree)
    !!{
    Return true if the given {\normalfont \ttfamily node} is evolvable.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeEventBranchJumpInterTree , nodeEventSubhaloPromotionInterTree, nodeEvent
    implicit none
    class           (mergerTreeEvolverStandard), intent(inout)          :: self
    type            (treeNode                 ), intent(inout)          :: node
    double precision                           , intent(in   )          :: timeEnd  , finalTimeInTree
    class           (nodeEvent                )               , pointer :: event
    class           (nodeComponentBasic       )               , pointer :: basic
    logical                                                             :: hasParent, treeLimited
    
    ! Evolve this node if it has a parent (or will transfer to another tree where it will have a parent), exists
    ! before the output time, has no children (i.e. they've already all been processed), and either exists before
    ! the final time in its tree, or exists precisely at that time and has some attached event yet to occur.
    basic       =>            node%basic ()
    event       =>            node%event
    hasParent   =  associated(node%parent  )
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
    standardNodeIsEvolvable= hasParent                           &
         &                  .and.                                &
         &                   .not.associated(node%firstChild  )  &
         &                  .and.                                &
         &                   (                                   &
         &                     basic%time() <  timeEnd           &
         &                    .or.                               &
         &                     (                                 &
         &                       .not.treeLimited                & ! For nodes that are not tree limited (i.e. have a node which
         &                      .and.                            & ! will jump to another tree), allow them to evolve if the node
         &                       basic%time() == timeEnd         & ! is at the end time also, since the jump may occur at that time.
         &                     )                                 &
         &                   )                                   &
         &                  .and.                                &
         &                   (                                   &
         &                      .not.treeLimited                 &
         &                   .or.                                &
         &                       basic%time() <  finalTimeInTree &
         &                   .or.                                &
         &                    (                                  &
         &                     (                                 &
         &                        associated(node%event      )   &
         &                      .or.                             &
         &                        associated(node%mergeTarget)   &
         &                     )                                 &
         &                     .and.                             &
         &                       basic%time() <= finalTimeInTree &
         &                    )                                  &
         &                  )
    return
  end function standardNodeIsEvolvable
  
  recursive function standardTimeEvolveTo(self,node,timeEnd,cosmologyFunctions_,mergerTreeEvolveTimestep_,mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report,nodeLock,lockType) result(evolveToTime)
    !!{
    Determine the time to which {\normalfont \ttfamily node} should be evolved.
    !!}
    use :: Display               , only : displayIndent                     , displayMessage        , displayUnindent, verbosityLevelInfo
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Error                 , only : Error_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic                , nodeComponentSatellite, nodeEvent      , nodeEventBranchJumpInterTree, &
          &                               nodeEventSubhaloPromotionInterTree, treeEvent             , treeNode
    use :: Merger_Tree_Timesteps , only : timestepTask
    use :: String_Handling       , only : operator(//)
    implicit none
    class           (mergerTreeEvolverStandard    ), intent(inout)                    :: self
    double precision                                                                  :: evolveToTime
    type            (treeNode                     ), intent(inout)          , pointer :: node
    double precision                               , intent(in   )                    :: timeEnd
    class           (cosmologyFunctionsClass      ), intent(inout)                    :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass), intent(inout)                    :: mergerTreeEvolveTimestep_
    class           (mergerTreeNodeEvolverClass   ), intent(inout)                    :: mergerTreeNodeEvolver_
    type            (treeNode                     )                         , pointer :: nodeSatellite                       , nodeSibling
    procedure       (timestepTask                 ), intent(  out)          , pointer :: timestepTask_
    class           (*                            ), intent(  out)          , pointer :: timestepSelf
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(  out), optional, pointer :: nodeLock
    type            (varying_string               ), intent(  out), optional          :: lockType
    procedure       (timestepTask                 )                         , pointer :: timestepTaskInternal
    class           (nodeComponentBasic           )                         , pointer :: basicParent                         , basicSatellite    , &
         &                                                                               basicSibling                        , basic
    class           (nodeComponentSatellite       )                         , pointer :: satelliteSatellite
    class           (nodeEvent                    )                         , pointer :: event
    class           (treeEvent                    )                         , pointer :: treeEvent_
    type            (varying_string               ), save                             :: message
    !$omp threadprivate(message)
    double precision                                                                  :: expansionFactor                     , expansionTimescale, &
         &                                                                               hostTimeLimit                       , time              , &
         &                                                                               timeEarliest                        , evolveToTimeStep  , &
         &                                                                               timeSatellite                       , timeNode
    logical                                                                           :: isLimitedByTimestepper
    character       (len=9                        )                                   :: timeFormatted

    ! Initially set to the global end time.
    evolveToTime=timeEnd
    if (report) call Evolve_To_Time_Report("start (target): ",evolveToTime)
    ! Initialize the lock node if present.
    if (present(nodeLock)) nodeLock => null()
    if (present(lockType)) lockType = 'null'
    ! Get the basic component of the node.
    basic    => node %basic()
    timeNode =  basic%time ()
    ! Ensure that the timestep does not exceed any event attached to the tree.
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
    ! Ensure that the timestep does not exceed any event attached to the node.
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
    evolveToTimeStep=mergerTreeEvolveTimestep_%timeEvolveTo(evolveToTime,node,timestepTaskInternal,timestepSelf,report,nodeLock,lockType)
    if (evolveToTimeStep <= evolveToTime) then
       evolveToTime           =  evolveToTimeStep
       timestepTask_          => timestepTaskInternal
       isLimitedByTimestepper =  .true.
    else
       timestepTask_          => null()
       timestepSelf           => null()
       isLimitedByTimestepper =  .false.
    end if
    if (report) call displayUnindent("done")
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
             evolveToTime          =basicParent%time()
             isLimitedByTimestepper=.false.
          end if
       end if
       if (report) call Evolve_To_Time_Report("promotion limit: ",evolveToTime)
    case (.true.)
       ! Do not let satellite evolve too far beyond parent.
       if (associated(node%parent%parent)) then
          ! The host halo has a parent, so use the host halo time to limit satellite evolution.
          basicParent =>     node       %parent%basic()
          time        =      basicParent       %time ()
       else
          ! The host halo has no parent. The satellite must therefore be evolving to some event (e.g. a merger). We have to allow
          ! it to evolve ahead of the host halo in this case to avoid deadlocks.
          basicParent =>     node       %parent%basic()
          time        =  max(basicParent       %time (),timeNode)
       end if
       ! Check if the host has a child.
       select case (associated(node%parent%firstChild))
       case (.true. )
          ! Host still has a child - do not let the satellite evolve beyond the host.
          hostTimeLimit=max(time,timeNode)
          ! Check for any merge targets directed at this node.
          nodeSatellite => node%firstMergee
          timeEarliest  =  huge(1.0d0)
          do while (associated(nodeSatellite))
             satelliteSatellite => nodeSatellite%satellite()
             if (nodeSatellite%isSatellite().and.nodeSatellite%parent%isProgenitorOf(node%parent)) &
                  & timeEarliest=min(timeEarliest,satelliteSatellite%timeOfMerging())
             nodeSatellite => nodeSatellite%siblingMergee
          end do
          if (timeEarliest < huge(1.0d0)) hostTimeLimit=max(hostTimeLimit,timeEarliest)
       case (.false.)
          ! Find current expansion timescale.
          if (time /= self%timeHostPrevious) then
             ! Memoize the host timestep as it depends only on the host time and can be re-used for all satellites in a host.
             if (self%timestepHostRelative > 0.0d0) then
                expansionFactor   =      cosmologyFunctions_%expansionFactor(time           )
                expansionTimescale=1.0d0/cosmologyFunctions_%expansionRate  (expansionFactor)
                self%timeStepHost =min(self%timestepHostRelative*expansionTimescale,self%timestepHostAbsolute)
             else
                ! Avoid use of expansion timescale if host absolute timestep is non-positive. This allows static universe cases to be handled.
                self%timeStepHost =                                                 self%timestepHostAbsolute
             end if
             ! Update the memoization time.
             self%timeHostPrevious=time
          end if
          hostTimeLimit=max(time+self%timeStepHost,timeNode)
          ! Check if this criterion will actually limit the evolution time.
          if (hostTimeLimit < evolveToTime) then
             ! Satellite evolution will be limited by being required to not advance too far ahead of the host halo. If the
             ! timestep we can advance is less than a specified fraction of what is possible, then skip this for now.
             if (hostTimeLimit-timeNode < self%fractionTimestepSatelliteMinimum*self%timeStepHost) then
                hostTimeLimit=timeNode
             end if
          end if
       end select
       ! Limit to this time.
       if (hostTimeLimit < evolveToTime) then
          if (present(nodeLock)) nodeLock => node%parent
          if (present(lockType)) lockType =  "satellite in host"
          evolveToTime          =hostTimeLimit
          isLimitedByTimestepper=.false.
       end if
       if (report) call Evolve_To_Time_Report("satellite in host limit: ",evolveToTime,node%parent%index())
    end select
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
             evolveToTime          =max(timeNode,basicSibling%time())
             isLimitedByTimestepper=.false.
          end if
          if (report) call Evolve_To_Time_Report("sibling: ",evolveToTime,nodeSibling%index())
          nodeSibling => nodeSibling%sibling
       end do
    end if
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Ensure that this node is not evolved beyond the time of any of its current satellites.
    if (timeNode < evolveToTime) then
       nodeSatellite => node%firstSatellite
       do while (associated(nodeSatellite))
          basicSatellite => nodeSatellite %basic()
          timeSatellite  =  basicSatellite%time ()
          if (timeSatellite < evolveToTime) then
             if (present(nodeLock)) nodeLock => nodeSatellite
             if (present(lockType)) lockType =  "hosted satellite"
             evolveToTime          =max(timeSatellite,timeNode)
             isLimitedByTimestepper=.false.
             if (evolveToTime == timeNode) exit
          end if
          if (report) then
             if (node%isPrimaryProgenitor()) then
                call Evolve_To_Time_Report("promotion limit: " ,evolveToTime,node%parent%index())
             else
                call Evolve_To_Time_Report("node merge limit: ",evolveToTime,node%parent%index())
             end if
          end if
          nodeSatellite => nodeSatellite%sibling
       end do
    end if
    ! Return early if the timestep is already zero.
    if (evolveToTime == timeNode) return
    ! Also ensure that this node is not evolved beyond the time at which any of its mergees merge. In some cases, the node may
    ! already be in the future of a mergee. In such cases, simply freeze it at the current time.
    ! then need to also figure out how to speed up satellite evolve checks
    nodeSatellite => node%firstMergee
    do while (associated(nodeSatellite))
       satelliteSatellite => nodeSatellite%satellite()
       if (max(satelliteSatellite%timeOfMerging(),timeNode) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) then
             write (timeFormatted,'(f7.4)') max(satelliteSatellite%timeOfMerging(),timeNode)
             lockType =  "mergee ("//trim(timeFormatted)//")"
          end if
          evolveToTime          =max(satelliteSatellite%timeOfMerging(),timeNode)
          isLimitedByTimestepper=.false.
       end if
       if (report) call Evolve_To_Time_Report("mergee limit: ",evolveToTime,nodeSatellite%index())
       nodeSatellite => nodeSatellite%siblingMergee
    end do
    ! If the timestepper class provided the limit, allow it to optionally refuse to evolve (e.g. if the step is too small to be
    ! efficient).
    if (isLimitedByTimestepper) then
       if (mergerTreeEvolveTimestep_%refuseToEvolve(node)) evolveToTime=timeNode
    end if
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
       if (.not.mergerTreeNodeEvolver_%isAccurate(timeNode,evolveToTime)) then
          ! End time is well before current time. This is an error. Call ourself with reporting switched on to generate a report
          ! on the time limits.
          message=message//' Gyr)'
          if (.not.report) time=self%timeEvolveTo(node,timeEnd,cosmologyFunctions_,mergerTreeEvolveTimestep_,mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.true.)
          call Error_Report(message//{introspection:location})
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
    use :: Galacticus_Nodes                   , only : mergerTree                   , nodeComponentBasic, nodeEvent, treeNode
    use :: Merger_Trees_Evolve_Deadlock_Status, only : enumerationDeadlockStatusType
    implicit none
    class           (*                            ), intent(inout)          :: self
    type            (mergerTree                   ), intent(in   )          :: tree
    type            (treeNode                     ), intent(inout), pointer :: node
    type            (enumerationDeadlockStatusType), intent(inout)          :: statusDeadlock
    class           (nodeEvent                    )               , pointer :: eventLast            , eventNext        , &
         &                                                                     event
    class           (nodeComponentBasic           )               , pointer :: basic
    double precision                                                        :: timeNode             , timeEventEarliest
    logical                                                                 :: mergerTreeEvolverDone
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
    use :: Galacticus_Nodes                   , only : mergerTree                   , treeEvent
    use :: Merger_Trees_Evolve_Deadlock_Status, only : enumerationDeadlockStatusType
    implicit none
    type            (mergerTree                   ), intent(inout), target  :: tree
    type            (enumerationDeadlockStatusType), intent(inout)          :: statusDeadlock
    type            (treeEvent                    )               , pointer :: eventLast            , eventNext, event
    double precision                                                        :: treeTimeEarliest
    logical                                                                 :: mergerTreeEvolverDone

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
