!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implements an OpenMP threaded merger tree evolver.
  !!}

  use :: Cosmology_Functions            , only : cosmologyFunctions              , cosmologyFunctionsClass
  use :: Galactic_Structure_Solvers     , only : galacticStructureSolver         , galacticStructureSolverClass
  use :: Galacticus_Nodes               , only : treeNode
  use :: Input_Parameters               , only : inputParameters
  use :: Kind_Numbers                   , only : kind_int8
  use :: Merger_Tree_Evolve_Profilers   , only : mergerTreeEvolveProfilerClass
  use :: Merger_Tree_Initialization     , only : mergerTreeInitializorClass
  use :: Merger_Tree_Timesteps          , only : mergerTreeEvolveTimestepClass   , timestepTask
  use :: Merger_Trees_Evolve_Node       , only : mergerTreeNodeEvolverClass
  use :: Merger_Trees_Evolve_Concurrency, only : mergerTreeEvolveConcurrencyClass

  ! Structure used to hold worker copies of objects.
  type :: worker
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(mergerTreeEvolveTimestepClass), pointer :: mergerTreeEvolveTimestep_ => null()
     class(mergerTreeNodeEvolverClass   ), pointer :: mergerTreeNodeEvolver_    => null()
     class(mergerTreeEvolveProfilerClass), pointer :: mergerTreeEvolveProfiler_ => null()
     class(galacticStructureSolverClass ), pointer :: galacticStructureSolver_  => null()
  end type worker

  ! Linked list used for storing nodes to be evolved/post-processed.
  type :: nodeTaskList
     type     (nodeTaskList)        , pointer :: next          => null()
     type     (treeNode    )        , pointer :: node          => null()
     procedure(timestepTask), nopass, pointer :: timestepTask_ => null()
     class    (*           )        , pointer :: timestepSelf  => null()
  end type nodeTaskList
  
  !![
  <mergerTreeEvolver name="mergerTreeEvolverThreaded">
   <description>
    An OpenMP threaded merger tree evolver. This class extends the \refClass{mergerTreeEvolverStandard} class. To evolve the tree,
    a list of evolvable nodes is constructed and then a set of parallel threads is spawned which take nodes from that list, evolve
    them, and add them to a second list for postprocessing. This repeats until tree evolution is completed.
    </description>
  </mergerTreeEvolver>
  !!]
  type, extends(mergerTreeEvolverClass) :: mergerTreeEvolverThreaded
     !!{
     Implementation of an OpenMP threaded merger tree evolver.
     !!}
     private
     class           (cosmologyFunctionsClass         ), pointer                   :: cosmologyFunctions_              => null()
     class           (mergerTreeEvolveTimestepClass   ), pointer                   :: mergerTreeEvolveTimestep_        => null()
     class           (galacticStructureSolverClass    ), pointer                   :: galacticStructureSolver_         => null()
     class           (mergerTreeNodeEvolverClass      ), pointer                   :: mergerTreeNodeEvolver_           => null()
     class           (mergerTreeInitializorClass      ), pointer                   :: mergerTreeInitializor_           => null()
     class           (mergerTreeEvolveProfilerClass   ), pointer                   :: mergerTreeEvolveProfiler_        => null()
     class           (mergerTreeEvolveConcurrencyClass), pointer                   :: mergerTreeEvolveConcurrency_     => null()
     type            (worker                          ), allocatable, dimension(:) :: workers
     logical                                                                       :: allTreesExistAtFinalTime                   , dumpTreeStructure    , &
          &                                                                           backtrackToSatellites                      , profileSteps         , &
          &                                                                           workersInitialized               =  .false.
     double precision                                                              :: timestepHostAbsolute                       , timestepHostRelative , &
          &                                                                           fractionTimestepSatelliteMinimum
     type            (deadlockList                    ), pointer                   :: deadlockHeadNode                 => null()
     type            (inputParameters                 ), pointer                   :: parameters                       => null()
   contains
     !![
     <methods>
       <method description="Find the time to which a node can be evolved." method="timeEvolveTo"      />
       <method description="Add a node to the deadlock list."              method="deadlockAddNode"   />
       <method description="Output a description of a deadlocked tree."    method="deadlockOutputTree"/>
     </methods>
     !!]
     final     ::                       threadedDestructor
     procedure :: evolve             => threadedEvolve
     procedure :: timeEvolveTo       => threadedTimeEvolveTo
     procedure :: deadlockAddNode    => threadedDeadlockAddNode
     procedure :: deadlockOutputTree => threadedDeadlockOutputTree
  end type mergerTreeEvolverThreaded

  interface mergerTreeEvolverThreaded
     !!{
     Constructors for the {\normalfont \ttfamily threaded} merger tree evolver.
     !!}
     module procedure threadedConstructorParameters
     module procedure threadedConstructorInternal
  end interface mergerTreeEvolverThreaded

  ! OpenMP thread worker number.
  integer :: numberWorker
  !$omp threadprivate(numberWorker)

contains

  function threadedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily threaded} merger tree evolver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolverThreaded       )                        :: self
    type            (inputParameters                 ), intent(inout), target :: parameters
    class           (cosmologyFunctionsClass         ), pointer               :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass   ), pointer               :: mergerTreeEvolveTimestep_
    class           (galacticStructureSolverClass    ), pointer               :: galacticStructureSolver_
    class           (mergerTreeNodeEvolverClass      ), pointer               :: mergerTreeNodeEvolver_
    class           (mergerTreeInitializorClass      ), pointer               :: mergerTreeInitializor_
    class           (mergerTreeEvolveProfilerClass   ), pointer               :: mergerTreeEvolveProfiler_
    class           (mergerTreeEvolveConcurrencyClass), pointer               :: mergerTreeEvolveConcurrency_
    type            (inputParameters                 ), pointer               :: parameters_
    logical                                                                   :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                       backtrackToSatellites           , profileSteps
    double precision                                                          :: timestepHostRelative            , timestepHostAbsolute, &
         &                                                                       fractionTimestepSatelliteMinimum

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
    !!]
    ! A galacticStructureSolver is built here. Even though this is not called explicitly by this mergerTreeEvolver, the
    ! galacticStructureSolver is expected to hook itself to any events which will trigger a change in galactic structure.
    !![
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="mergerTreeEvolveTimestep"    name="mergerTreeEvolveTimestep_"    source="parameters"/>
    <objectBuilder class="galacticStructureSolver"     name="galacticStructureSolver_"     source="parameters"/>
    <objectBuilder class="mergerTreeNodeEvolver"       name="mergerTreeNodeEvolver_"       source="parameters"/>
    <objectBuilder class="mergerTreeInitializor"       name="mergerTreeInitializor_"       source="parameters"/>
    <objectBuilder class="mergerTreeEvolveProfiler"    name="mergerTreeEvolveProfiler_"    source="parameters"/>
    <objectBuilder class="mergerTreeEvolveConcurrency" name="mergerTreeEvolveConcurrency_" source="parameters"/>
    !!]
    if (associated(parameters%parent)) then
       parameters_ => parameters%parent
    else
       parameters_ => parameters
    end if
    self=mergerTreeEvolverThreaded(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,profileSteps,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,mergerTreeEvolveConcurrency_,galacticStructureSolver_,mergerTreeEvolveProfiler_,parameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="mergerTreeEvolveTimestep_"/>
    <objectDestructor name="mergerTreeNodeEvolver_"   />
    <objectDestructor name="galacticStructureSolver_" />
    <objectDestructor name="mergerTreeInitializor_"   />
    <objectDestructor name="mergerTreeEvolveProfiler_"/>
    !!]
    return
  end function threadedConstructorParameters

  function threadedConstructorInternal(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,backtrackToSatellites,profileSteps,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,mergerTreeEvolveConcurrency_,galacticStructureSolver_,mergerTreeEvolveProfiler_,parameters) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily threaded} merger tree evolver class.
    !!}
    implicit none
    type            (mergerTreeEvolverThreaded       )                        :: self
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (mergerTreeEvolveTimestepClass   ), intent(in   ), target :: mergerTreeEvolveTimestep_
    class           (galacticStructureSolverClass    ), intent(in   ), target :: galacticStructureSolver_
    class           (mergerTreeNodeEvolverClass      ), intent(in   ), target :: mergerTreeNodeEvolver_
    class           (mergerTreeInitializorClass      ), intent(in   ), target :: mergerTreeInitializor_
    class           (mergerTreeEvolveProfilerClass   ), intent(in   ), target :: mergerTreeEvolveProfiler_
    class           (mergerTreeEvolveConcurrencyClass), intent(in   ), target :: mergerTreeEvolveConcurrency_
    logical                                           , intent(in   )         :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                       backtrackToSatellites           , profileSteps
    double precision                                  , intent(in   )         :: timestepHostRelative            , timestepHostAbsolute, &
         &                                                                       fractionTimestepSatelliteMinimum
    type            (inputParameters                 ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="allTreesExistAtFinalTime, dumpTreeStructure, timestepHostRelative, timestepHostAbsolute, fractionTimestepSatelliteMinimum, backtrackToSatellites, profileSteps, *cosmologyFunctions_, *mergerTreeNodeEvolver_, *mergerTreeEvolveTimestep_, *mergerTreeInitializor_, *mergerTreeEvolveConcurrency_, *galacticStructureSolver_, *mergerTreeEvolveProfiler_, *parameters"/>
    !!]

    self%deadlockHeadNode   => null()
    self%workersInitialized =  .false.
    return
  end function threadedConstructorInternal

  subroutine threadedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily threaded} merger tree evolver class.
    !!}
    implicit none
    type   (mergerTreeEvolverThreaded), intent(inout) :: self
    integer                                           :: i

    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%mergerTreeEvolveTimestep_"   />
    <objectDestructor name="self%galacticStructureSolver_"    />
    <objectDestructor name="self%mergerTreeNodeEvolver_"      />
    <objectDestructor name="self%mergerTreeInitializor_"      />
    <objectDestructor name="self%mergerTreeEvolveProfiler_"   />
    <objectDestructor name="self%mergerTreeEvolveConcurrency_"/>
    !!]
    if (self%workersInitialized) then
       do i=lbound(self%workers,dim=1),ubound(self%workers,dim=1)
          !![
	  <objectDestructor name="self%workers(i)%cosmologyFunctions_"      />
	  <objectDestructor name="self%workers(i)%mergerTreeEvolveTimestep_"/>
	  <objectDestructor name="self%workers(i)%mergerTreeNodeEvolver_"   />
	  <objectDestructor name="self%workers(i)%mergerTreeEvolveProfiler_"/>
	  <objectDestructor name="self%workers(i)%galacticStructureSolver_" />
          !!]
       end do
    end if
    return
  end subroutine threadedDestructor

  subroutine threadedEvolve(self,tree,timeEnd,treeDidEvolve,suspendTree,deadlockReporting,systemClockMaximum,initializationLock,status)
    !!{
    Evolves all properties of a merger tree to the specified time.
    !!}
    use    :: Display                            , only : displayIndent                    , displayMessage                     , displayUnindent              , displayVerbosity           , &
         &                                                displayGreen                     , displayReset                       , enumerationVerbosityLevelType, verbosityLevelWarn
    use    :: Events_Hooks                       , only : eventsHooksFutureThread
    use    :: Error                              , only : Error_Report                     , errorStatusSuccess
    use    :: Galacticus_Nodes                   , only : interruptTask                    , mergerTree                         , nodeComponentBasic           , nodeEvent                  , &
          &                                               nodeEventBranchJumpInterTree     , nodeEventSubhaloPromotionInterTree , treeNode                     , treeNodeLinkedList
    use    :: Merger_Tree_Walkers                , only : mergerTreeWalkerAllNodes
    use    :: Merger_Trees_Dump                  , only : Merger_Tree_Dump
    use    :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsDeadlocked       , deadlockStatusIsNotDeadlocked      , deadlockStatusIsReporting    , deadlockStatusIsSuspendable, &
         &                                                enumerationDeadlockStatusType
    use    :: Node_Components                    , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    !$ use :: OMP_Lib                            , only : OMP_Set_Lock                     , OMP_Unset_Lock                     , omp_lock_kind                , OMP_Get_Max_Threads        , &
    !$   &                                                OMP_Get_Thread_Num
    use    :: Locks                              , only : ompLock
    use    :: String_Handling                    , only : operator(//)
    implicit none
    class           (mergerTreeEvolverThreaded    )                    , intent(inout) :: self
    integer                                        , optional          , intent(  out) :: status
    type            (mergerTree                   )           , target , intent(inout) :: tree
    double precision                                                   , intent(in   ) :: timeEnd
    logical                                                            , intent(  out) :: treeDidEvolve                                  , suspendTree
    logical                                                            , intent(in   ) :: deadlockReporting
    integer         (kind_int8                    ), optional          , intent(in   ) :: systemClockMaximum
    integer         (omp_lock_kind                ), optional          , intent(inout) :: initializationLock
    type            (treeNode                     )           , pointer                :: nodeParent
    type            (treeNode                     ), save     , pointer                :: node                                           , nodeLock
    !$omp threadprivate(node,nodeLock)
    double precision                               , save                              :: timeEndThisNode
    !$omp threadprivate(timeEndThisNode)
    type            (varying_string               ), save                              :: vMessage
    !$omp threadprivate(vMessage)
    class           (nodeEvent                    )           , pointer                :: event
    double precision                               , parameter                         :: timeTolerance              =1.0d-5
    double precision                               , parameter                         :: largeTime                  =1.0d10
    procedure       (interruptTask                ), save     , pointer                :: interruptProcedure
    !$omp threadprivate(interruptProcedure)
    procedure       (timestepTask                 ), save     , pointer                :: timestepTask_
    !$omp threadprivate(timestepTask_)
    class           (*                            ), save     , pointer                :: timestepSelf
    !$omp threadprivate(timestepSelf)
    type            (enumerationVerbosityLevelType), parameter                         :: verbosityLevel             =verbosityLevelWarn
    class           (nodeComponentBasic           )           , pointer                :: basicBase                                     , basicParent
    class           (nodeComponentBasic           ), save     , pointer                :: basic
    !$omp threadprivate(basic)
    type            (mergerTree                   ), save     , pointer                :: currentTree
    !$omp threadprivate(currentTree)
    type            (nodeTaskList                 )           , pointer                :: evolveStackHead                               , evolveStackTip        , &
         &                                                                                postProcessStackHead                          , postProcessStackTip
    type            (nodeTaskList                 ), save     , pointer                :: stackNext
    !$omp threadprivate(stackNext)
    type            (inputParameters              ), save     , allocatable            :: parameters
    !$omp threadprivate(parameters)
    type            (ompLock                      )                                    :: lockEvolveStack                               , lockPostProcessStack , &
         &                                                                                lockTree
    type            (mergerTreeWalkerAllNodes     )                                    :: treeWalker
    type            (enumerationDeadlockStatusType)                                    :: statusDeadlock
    integer                                                                            :: treeWalkCountPreviousOutput                   , nodesEvolvedCount    , &
         &                                                                                nodesTotalCount                               , treeWalkCount        , &
         &                                                                                countWorkers                                  , status_
    double precision                                                                   :: earliestTimeInTree                            , finalTimeInTree
    character       (len=24                       )                                    :: label
    character       (len=35                       )                                    :: message
    type            (varying_string               ), save                              :: lockType       
    !$omp threadprivate(lockType)                               
    logical                                                                            :: anyTreeExistsAtOutputTime                     , hasIntertreeEvent    , &
         &                                                                                hasParent                                     , treeLimited          , &
         &                                                                                didEvolve                                     , evolutionFailed      , &
         &                                                                                evolutionExiting
    logical                                        , save                              :: interrupted
    !$omp threadprivate(interrupted)
    integer                                        , save                              :: evolutionPhase                                , evolutionPhaseMaximum
    !$omp threadprivate(evolutionPhase)
    
    ! Iterate through all trees.
    suspendTree               =  .false.
    anyTreeExistsAtOutputTime =  .false.
    treeDidEvolve             =  .false.
    evolutionExiting          =  .false.
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
    ! Create copies of objects needed for evolution.
    if (.not.self%workersInitialized) then
       countWorkers   =1
       !$ countWorkers=OMP_Get_Max_Threads()
       allocate(self%workers(0:countWorkers-1))
       self%workersInitialized=.true.
       do numberWorker=0,countWorkers-1
          call eventsHooksFutureThread(numberWorker)
          allocate(self%workers(numberWorker)%cosmologyFunctions_      ,mold=self%cosmologyFunctions_      )
          allocate(self%workers(numberWorker)%mergerTreeEvolveTimestep_,mold=self%mergerTreeEvolveTimestep_)
          allocate(self%workers(numberWorker)%mergerTreeNodeEvolver_   ,mold=self%mergerTreeNodeEvolver_   )
          allocate(self%workers(numberWorker)%mergerTreeEvolveProfiler_,mold=self%mergerTreeEvolveProfiler_)
          allocate(self%workers(numberWorker)%galacticStructureSolver_ ,mold=self%galacticStructureSolver_ )
          !$omp critical(mergerTreeEvolverThreadedDeepCopy)
          !![
	  <deepCopyReset variables="self%cosmologyFunctions_ self%mergerTreeEvolveTimestep_ self%mergerTreeNodeEvolver_ self%mergerTreeEvolveProfiler_ self%galacticStructureSolver_"/>
	  <deepCopy source="self%cosmologyFunctions_"       destination="self%workers(numberWorker)%cosmologyFunctions_"      />
	  <deepCopy source="self%mergerTreeEvolveTimestep_" destination="self%workers(numberWorker)%mergerTreeEvolveTimestep_"/>
	  <deepCopy source="self%mergerTreeNodeEvolver_"    destination="self%workers(numberWorker)%mergerTreeNodeEvolver_"   />
	  <deepCopy source="self%mergerTreeEvolveProfiler_" destination="self%workers(numberWorker)%mergerTreeEvolveProfiler_"/>
	  <deepCopy source="self%galacticStructureSolver_"  destination="self%workers(numberWorker)%galacticStructureSolver_" />
	  <deepCopyFinalize variables="self%workers(numberWorker)%cosmologyFunctions_ self%workers(numberWorker)%mergerTreeEvolveTimestep_ self%workers(numberWorker)%mergerTreeNodeEvolver_ self%workers(numberWorker)%mergerTreeEvolveProfiler_ self%workers(numberWorker)%galacticStructureSolver_"/>  
          !!]
          !$omp end critical(mergerTreeEvolverThreadedDeepCopy)
          call eventsHooksFutureThread()
        end do
    end if
    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible. Begin a parallel region
    ! here - this allows us to retain a team of threads for the entirety of this trees' evolution. We will use OpenMP `single`
    ! sections for work that must be done in serial.
    !$omp parallel
    ! Initialize nodes for these threads.
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call parameters%parametersGroupCopy(self%parameters)
    call Node_Components_Thread_Initialize(parameters)   
    ! Determine out worker number.
    numberWorker=OMP_Get_Thread_Num()
    ! Initialize evolution state and locks.
    !$omp single
    didEvolve                  =.true.
    treeWalkCount              =0
    treeWalkCountPreviousOutput=0
    lockEvolveStack            =ompLock()
    lockPostProcessStack       =ompLock()
    lockTree                   =ompLock()
    !$omp end single
    outerLoop: do while (didEvolve) ! Keep looping through the tree until we make a pass during which no nodes were evolved.
       !$omp barrier
       !$omp single
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
       !$omp end single
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
                !$omp single
                node            => currentTree%nodeBase
                basic           => node       %basic   ()
                finalTimeInTree =  basic      %time    ()
                ! Report on current tree if deadlocked.
                if (statusDeadlock == deadlockStatusIsReporting) then
                   vMessage="tree "
                   vMessage=vMessage//currentTree%index
                   call displayIndent(vMessage)
                end if
                call self%mergerTreeEvolveConcurrency_%initializeTree()
                evolutionPhaseMaximum=self%mergerTreeEvolveConcurrency_%countPhases()
                !$omp end single
                ! Build a stack of all evolvable nodes.
                evolutionPhase=0
                concurrencyLoop : do while(evolutionPhase < evolutionPhaseMaximum)
                   evolutionPhase=evolutionPhase+1
                   !$omp single
                   evolveStackHead => null()
                   evolveStackTip  => null()
                   treeWalker      =  mergerTreeWalkerAllNodes(currentTree,spanForest=.false.)
                   do while (treeWalker%next(node))
                      ! Skip nodes that are not to be evolved during this phase/
                      if (.not.self%mergerTreeEvolveConcurrency_%includeInEvolution(evolutionPhase,node)) cycle
                      ! Get the basic component of the node.
                      basic => node%basic()
                      ! Count nodes in the tree.
                      nodesTotalCount=nodesTotalCount+1
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
                      if     (                                     &
                           &        hasParent                      &
                           &  .and.                                &
                           &   .not.associated(node%firstChild  )  &
                           &  .and.                                &
                           &   (                                   &
                           &     basic%time() <  timeEnd           &
                           &    .or.                               &
                           &     (                                 &
                           &       .not.treeLimited                & ! For nodes that are not tree limited (i.e. have a node which
                           &      .and.                            & ! will jump to another tree), allow them to evolve if the node
                           &       basic%time() == timeEnd         & ! is at the end time also, since the jump may occur at that time.
                           &     )                                 &
                           &   )                                   &
                           &  .and.                                &
                           &   (                                   &
                           &      .not.treeLimited                 &
                           &   .or.                                &
                           &       basic%time() <  finalTimeInTree &
                           &   .or.                                &
                           &    (                                  &
                           &     (                                 &
                           &        associated(node%event      )   &
                           &      .or.                             &
                           &        associated(node%mergeTarget)   &
                           &     )                                 &
                           &     .and.                             &
                           &       basic%time() <= finalTimeInTree &
                           &    )                                  &
                           &   )                                   &
                           & ) then
                         if (associated(evolveStackHead)) then
                            allocate(evolveStackTip%next)
                            evolveStackTip => evolveStackTip %next
                         else
                            allocate(evolveStackHead)
                            evolveStackTip => evolveStackHead
                         end if
                         evolveStackTip%node => node
                         evolveStackTip%timestepTask_ => null()
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
                      end if
                   end do
                   evolveStackTip        => evolveStackHead
                   postProcessStackHead  => null()
                   postProcessStackTip   => null()
                   evolutionFailed       =  .false.
                   evolutionPhaseMaximum = self%mergerTreeEvolveConcurrency_%countPhases()
                   !$omp end single
                   ! Evolve all nodes on the stack.
                   do while (.not.evolutionFailed)
                      call lockEvolveStack%set()
                      if (associated(evolveStackTip)) then
                         node           => evolveStackTip%node
                         basic          => node          %basic()
                         stackNext      => evolveStackTip%next
                         deallocate(evolveStackTip)
                         evolveStackTip =>                stackNext
                         call lockEvolveStack%unset()
                      else
                         call lockEvolveStack%unset()
                         exit
                      end if
                      ! Flag that a node was evolved.
                      didEvolve        =.true.
                      ! Update tree progress counter.
                      !$omp atomic
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
                         else if (self%profileSteps) then
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.false.                  ,lockType=lockType)
                            call self%workers(numberWorker)%mergerTreeEvolveProfiler_%stepDescriptor(lockType)
                         else
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.false.                                    )
                         end if
                         ! If this node is able to evolve by a finite amount, the tree is not deadlocked.
                         if (timeEndThisNode > basic%time()) then
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                            ! Update record of earliest time in the tree.
                            earliestTimeInTree=min(earliestTimeInTree,timeEndThisNode)
                            ! Evolve the node to the next interrupt event, or the end time.
                            !![
			    <conditionalCall>
			      <call>
				call self%workers(numberWorker)%mergerTreeNodeEvolver_%evolve(currentTree,node,timeEndThisNode,interrupted,interruptProcedure,self%workers(numberWorker)%galacticStructureSolver_,lockTree,systemClockMaximum{conditions})
			      </call>
			      <argument name="status" value="status_" condition="present(status)"/>
			    </conditionalCall>
                            !!]
                            if (present(status)) then
                               if (status_ /= errorStatusSuccess) then
                                  status         =status_
                                  evolutionFailed=.true.
                                  cycle
                               end if
                            end if
                         end if
                         ! Check for interrupt. If no interrupt occurred, we will exit the evolution loop. Any end of timestep task
                         ! will be handled during the post-processing phase.
                         if (interrupted) then
                            ! If an interrupt occurred call the specified procedure to handle it.
                            call interruptProcedure(node)
                            ! Something happened so the tree is not deadlocked.
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                         end if
                      end do
                      ! Push the node onto a postprocessing stack, if it still exists.
                      if (associated(node)) then
                         call lockPostProcessStack%set()
                         if (associated(postProcessStackHead)) then
                            allocate(postProcessStackTip%next)
                            postProcessStackTip => postProcessStackTip %next
                         else
                            allocate(postProcessStackHead)
                            postProcessStackTip => postProcessStackHead
                         end if
                         postProcessStackTip%node          => node
                         postProcessStackTip%timestepTask_ => timestepTask_
                         postProcessStackTip%timestepSelf  => timestepSelf
                         call lockPostProcessStack%unset()
                      end if
                   end do
                   !$omp barrier
                   !$omp single
                   if (evolutionFailed) then
                      ! Evolution failed - deallocate any remaining stack and then exit and return.
                      evolveStackTip      => evolveStackHead
                      do while (associated(evolveStackTip     ))
                         node                => evolveStackTip     %node
                         stackNext           => evolveStackTip     %next
                         deallocate(evolveStackTip     )
                         evolveStackTip      =>                     stackNext
                      end do
                      postProcessStackTip => postProcessStackHead
                      do while (associated(postProcessStackTip))
                         node                => postProcessStackTip%node
                         stackNext           => postProcessStackTip%next
                         deallocate(postProcessStackTip)
                         postProcessStackTip =>                     stackNext
                      end do
                      evolutionExiting=.true.
                   else
                      ! Perform post-processing of nodes.
                      postProcessStackTip => postProcessStackHead
                      do while (associated(postProcessStackTip))
                         ! Pop a node from the stack.
                         node                => postProcessStackTip%node
                         timestepTask_       => postProcessStackTip%timestepTask_
                         timestepSelf        => postProcessStackTip%timestepSelf
                         stackNext           => postProcessStackTip%next
                         deallocate(postProcessStackTip)
                         postProcessStackTip =>                     stackNext
                         ! Handle any end-of-timestep task associated with this node.
                         if (associated(timestepTask_))                                          &
                              & call timestepTask_(timestepSelf,currentTree,node,statusDeadlock)
                         ! Handle node promotion or merging.
                         if (associated(node%parent)) then
                            nodeParent  => node      %parent
                            basic       => node      %basic ()
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
                                     call self%workers(numberWorker)%mergerTreeNodeEvolver_%merge(node)
                                  end if
                               case (.true.)
                                  ! This is the major progenitor, so promote the node to its parent providing that the node has no
                                  ! siblings - this ensures that any siblings have already been evolved and become satellites of the
                                  ! parent halo. Also record that the tree is not deadlocked, as we are changing the tree state.
                                  if (.not.associated(node%sibling).and..not.associated(node%event)) then
                                     statusDeadlock=deadlockStatusIsNotDeadlocked
                                     call self%workers(numberWorker)%mergerTreeNodeEvolver_%promote(node)
                                  end if
                               end select
                            end if
                         end if
                      end do
                   end if
                   !$omp end single                
                   if (evolutionExiting) exit
                end do concurrencyLoop
                !$omp barrier
                if (evolutionExiting) exit
                !$omp single
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
                !$omp end single
             end if
             ! Move to the next tree.
             currentTree => currentTree%nextTree
          end do treesLoop
          !$omp barrier
          if (evolutionExiting) exit
          !$omp single
          ! Perform any tree events.
          call threadedTreeEventsPerform(tree,statusDeadlock)
          ! Check deadlocking.
          if (didEvolve .and. statusDeadlock /= deadlockStatusIsNotDeadlocked) then
             if (statusDeadlock == deadlockStatusIsReporting) then
                call displayUnindent("report done")
                call self%deadlockOutputTree(timeEnd)
                if (.not.deadlockReporting) then
                   call Error_Report('merger tree appears to be deadlocked (see preceding report) - check timestep criteria'//{introspection:location})
                else
                   ! Exit evolution and return
                   evolutionExiting=.true.
                end if
             else
                ! Tree appears to be deadlocked. Check if it is suspendable.
                if (statusDeadlock == deadlockStatusIsSuspendable) then
                   ! Tree is suspendable, so do not attempt to process further, but simply return and flag it for suspension.
                   suspendTree     =.true.
                   evolutionExiting=.true.
                else
                   ! Tree is truly deadlocked. Switch to reporting mode and do one more pass through the tree.
                   statusDeadlock=deadlockStatusIsReporting
                end if
             end if
          else
             ! No evolution could occur, so the tree is not deadlocked.
             statusDeadlock=deadlockStatusIsNotDeadlocked
          end if
          !$omp end single
          if (evolutionExiting) exit
       end do deadlock
       !$omp barrier
       if (evolutionExiting) exit
       ! Record tree evolution status.
       !$omp single
       if (didEvolve) treeDidEvolve=.true.
       !$omp end single
    end do outerLoop
    !$omp barrier
    call Node_Components_Thread_Uninitialize()
    !$omp barrier
    !$omp critical(evolverReset)
    call parameters%reset()
    !$omp end critical(evolverReset)
    deallocate(parameters)
    !$omp end parallel
    return
  end subroutine threadedEvolve

  recursive function threadedTimeEvolveTo(self,node,timeEnd,timestepTask_,timestepSelf,report,nodeLock,lockType) result(evolveToTime)
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
    class           (mergerTreeEvolverThreaded    ), intent(inout)                    :: self
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
    class           (nodeComponentBasic           )                         , pointer :: basicParent           , basicSatellite    , &
         &                                                                               basicSibling          , basic
    class           (nodeComponentSatellite       )                         , pointer :: satelliteSatellite
    class           (nodeEvent                    )                         , pointer :: event
    class           (treeEvent                    )                         , pointer :: treeEvent_
    double precision                                                                  :: expansionFactor       , expansionTimescale, &
         &                                                                               hostTimeLimit         , time              , &
         &                                                                               timeEarliest          , evolveToTimeStep  , &
         &                                                                               hostTimeStep          , timeNode          , &
         &                                                                               timeSatellite
    logical                                                                           :: isLimitedByTimestepper
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
          timestepTask_ => threadedNodeEventsPerform
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
    evolveToTimeStep=self%workers(numberWorker)%mergerTreeEvolveTimestep_%timeEvolveTo(evolveToTime,node,timestepTaskInternal,timestepSelf,report,nodeLock,lockType)
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
    ! Ensure that this node is not evolved beyond the time of any of its current satellites.
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       basicSatellite => nodeSatellite %basic()
       timeSatellite  =  basicSatellite%time ()
       if (max(timeSatellite,timeNode) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) lockType =  "hosted satellite"
          evolveToTime          =max(timeSatellite,timeNode)
          isLimitedByTimestepper=.false.
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
          evolveToTime          =max(satelliteSatellite%timeOfMerging(),timeNode)
          isLimitedByTimestepper=.false.
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
             evolveToTime          =max(timeNode,basicSibling%time())
             isLimitedByTimestepper=.false.
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
             evolveToTime          =basicParent%time()
             isLimitedByTimestepper=.false.
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
             expansionFactor   =      self%workers(numberWorker)%cosmologyFunctions_%expansionFactor(time           )
             expansionTimescale=1.0d0/self%workers(numberWorker)%cosmologyFunctions_%expansionRate  (expansionFactor)
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
          evolveToTime          =hostTimeLimit
          isLimitedByTimestepper=.false.
       end if
       if (report) call Evolve_To_Time_Report("satellite in host limit: ",evolveToTime,node%parent%index())
    end select
    ! If the timestepper class provided the limit, allow it to optionally refuse to evolve (e.g. if the step is too small to be
    ! efficient).
    if (isLimitedByTimestepper) then
       if (self%workers(numberWorker)%mergerTreeEvolveTimestep_%refuseToEvolve(node)) evolveToTime=timeNode
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
       if (.not.self%workers(numberWorker)%mergerTreeNodeEvolver_%isAccurate(timeNode,evolveToTime)) then
          ! End time is well before current time. This is an error. Call ourself with reporting switched on to generate a report
          ! on the time limits.
          message=message//' Gyr)'
          if (.not.report) time=self%timeEvolveTo(node,timeEnd,timestepTask_,timestepSelf,report=.true.)
          call Error_Report(message//{introspection:location})
       else
          ! End time is before current time, but only by a small amount, simply reset the current time to the end time.
          message=message//' Gyr) - this should happen infrequently'
          call displayMessage(message,verbosityLevelInfo)
          call basic%timeSet(evolveToTime)
       end if
    end if
    return
  end function threadedTimeEvolveTo

  subroutine threadedDeadlockAddNode(self,node,treeIndex,nodeLock,lockType)
    !!{
    Add a node to the deadlocked nodes list.
    !!}
    implicit none
    class  (mergerTreeEvolverThreaded), intent(inout)          :: self
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
  end subroutine threadedDeadlockAddNode

  subroutine threadedDeadlockOutputTree(self,timeEnd)
    !!{
    Output the deadlocked nodes in {\normalfont \ttfamily dot} format.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: String_Handling , only : operator(//)
    implicit none
    class           (mergerTreeEvolverThreaded), intent(inout) :: self
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
  end subroutine threadedDeadlockOutputTree

  subroutine threadedNodeEventsPerform(self,tree,node,statusDeadlock)
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
  end subroutine threadedNodeEventsPerform

  subroutine threadedTreeEventsPerform(tree,statusDeadlock)
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
  end subroutine threadedTreeEventsPerform
