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
  Implements an OpenMP threaded merger tree evolver.
  !!}

  use :: Merger_Trees_Evolve_Concurrency, only : mergerTreeEvolveConcurrencyClass

  ! Structure used to hold worker copies of objects.
  type :: worker
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(mergerTreeEvolveTimestepClass), pointer :: mergerTreeEvolveTimestep_ => null()
     class(mergerTreeNodeEvolverClass   ), pointer :: mergerTreeNodeEvolver_    => null()
     class(mergerTreeEvolveProfilerClass), pointer :: mergerTreeEvolveProfiler_ => null()
     class(galacticStructureSolverClass ), pointer :: galacticStructureSolver_  => null()
     class(metaTreeProcessingTimeClass  ), pointer :: metaTreeProcessingTime_   => null()
  end type worker

  ! Linked list used for storing nodes to be evolved/post-processed.
  type :: nodeTask
     type     (treeNode    )        , pointer :: node          => null()
     procedure(timestepTask), nopass, pointer :: timestepTask_ => null()
     class    (*           )        , pointer :: timestepSelf  => null()
  end type nodeTask
  
  !![
  <mergerTreeEvolver name="mergerTreeEvolverThreaded">
   <description>
    An OpenMP threaded merger tree evolver. This class extends the \refClass{mergerTreeEvolverStandard} class. To evolve the tree,
    a list of evolvable nodes is constructed and then a set of parallel threads is spawned which take nodes from that list, evolve
    them, and add them to a second list for postprocessing. This repeats until tree evolution is completed.
    </description>
  </mergerTreeEvolver>
  !!]
  type, extends(mergerTreeEvolverStandard) :: mergerTreeEvolverThreaded
     !!{
     Implementation of an OpenMP threaded merger tree evolver.
     !!}
     private
     type   (worker                          ), allocatable, dimension(:) :: workers
     class  (mergerTreeEvolveConcurrencyClass), pointer                   :: mergerTreeEvolveConcurrency_ => null()
     type   (inputParameters                 ), pointer                   :: parameters                   => null()
     logical                                                              :: workersInitialized           =  .false., reportTiming
     integer                                                              :: timeEvolveID
   contains
     final     ::           threadedDestructor
     procedure :: evolve => threadedEvolve
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
    class           (metaTreeProcessingTimeClass     ), pointer               :: metaTreeProcessingTime_
    type            (inputParameters                 ), pointer               :: parameters_
    logical                                                                   :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                       profileSteps                    , reportTiming
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
      <name>profileSteps</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not to profile the ODE evolver.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportTiming</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, report on timing of serial and parallel sections.</description>
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
    <objectBuilder class="metaTreeProcessingTime"      name="metaTreeProcessingTime_"      source="parameters"/>
    !!]
    if (associated(parameters%parent)) then
       parameters_ => parameters%parent
    else
       parameters_ => parameters
    end if
    self=mergerTreeEvolverThreaded(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,profileSteps,reportTiming,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,mergerTreeEvolveConcurrency_,metaTreeProcessingTime_,galacticStructureSolver_,mergerTreeEvolveProfiler_,parameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="mergerTreeEvolveTimestep_"   />
    <objectDestructor name="mergerTreeNodeEvolver_"      />
    <objectDestructor name="galacticStructureSolver_"    />
    <objectDestructor name="mergerTreeInitializor_"      />
    <objectDestructor name="mergerTreeEvolveProfiler_"   />
    <objectDestructor name="mergerTreeEvolveConcurrency_"/>
    <objectDestructor name="metaTreeProcessingTime_"     />
    !!]
    return
  end function threadedConstructorParameters

  function threadedConstructorInternal(allTreesExistAtFinalTime,dumpTreeStructure,timestepHostRelative,timestepHostAbsolute,fractionTimestepSatelliteMinimum,profileSteps,reportTiming,cosmologyFunctions_,mergerTreeNodeEvolver_,mergerTreeEvolveTimestep_,mergerTreeInitializor_,mergerTreeEvolveConcurrency_,metaTreeProcessingTime_,galacticStructureSolver_,mergerTreeEvolveProfiler_,parameters) result(self)
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
    class           (metaTreeProcessingTimeClass     ), intent(in   ), target :: metaTreeProcessingTime_
    logical                                           , intent(in   )         :: allTreesExistAtFinalTime        , dumpTreeStructure   , &
         &                                                                       profileSteps                    , reportTiming
    double precision                                  , intent(in   )         :: timestepHostRelative            , timestepHostAbsolute, &
         &                                                                       fractionTimestepSatelliteMinimum
    type            (inputParameters                 ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="allTreesExistAtFinalTime, dumpTreeStructure, timestepHostRelative, timestepHostAbsolute, fractionTimestepSatelliteMinimum, profileSteps, reportTiming, *cosmologyFunctions_, *mergerTreeNodeEvolver_, *mergerTreeEvolveTimestep_, *mergerTreeInitializor_, *mergerTreeEvolveConcurrency_, *metaTreeProcessingTime_, *galacticStructureSolver_, *mergerTreeEvolveProfiler_, *parameters"/>
    !!]

    self%deadlockHeadNode   => null()
    self%workersInitialized =  .false.
    self%timeHostPrevious   =  -huge(0.0d0)
    self%timeStepHost       =  -huge(0.0d0)
    !![
    <addMetaProperty component="basic" name="timeEvolve" id="self%timeEvolveID" isCreator="yes"/>
    !!]
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
    <objectDestructor name="self%mergerTreeEvolveConcurrency_"/>
    <objectDestructor name="self%metaTreeProcessingTime_"     />
    !!]
    if (self%workersInitialized) then
       do i=lbound(self%workers,dim=1),ubound(self%workers,dim=1)
          !![
	  <objectDestructor name="self%workers(i)%cosmologyFunctions_"      />
	  <objectDestructor name="self%workers(i)%mergerTreeEvolveTimestep_"/>
	  <objectDestructor name="self%workers(i)%mergerTreeNodeEvolver_"   />
	  <objectDestructor name="self%workers(i)%mergerTreeEvolveProfiler_"/>
	  <objectDestructor name="self%workers(i)%galacticStructureSolver_" />
          <objectDestructor name="self%workers(i)%metaTreeProcessingTime_"  />
          !!]
       end do
    end if
    return
  end subroutine threadedDestructor

  subroutine threadedEvolve(self,tree,timeEnd,treeDidEvolve,suspendTree,deadlockReporting,systemClockMaximum,initializationLock,status)
    !!{
    Evolves all properties of a merger tree to the specified time.
    !!}
    use, intrinsic :: ISO_C_Binding                      , only : c_size_t
    use            :: Display                            , only : displayIndent                    , displayMessage                     , displayUnindent              , displayVerbosity           , &
         &                                                        displayGreen                     , displayReset                       , enumerationVerbosityLevelType, verbosityLevelWarn
    use            :: Events_Hooks                       , only : eventsHooksFutureThread
    use            :: Error                              , only : Error_Report                     , errorStatusSuccess
    use            :: Galacticus_Nodes                   , only : interruptTask                    , mergerTree                         , nodeComponentBasic
    use            :: Merger_Tree_Walkers                , only : mergerTreeWalkerAllNodes
    use            :: Merger_Trees_Dump                  , only : Merger_Tree_Dump
    use            :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsDeadlocked       , deadlockStatusIsNotDeadlocked      , deadlockStatusIsReporting    , deadlockStatusIsSuspendable, &
         &                                                 enumerationDeadlockStatusType
    use            :: Node_Components                    , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    !$ use         :: OMP_Lib                            , only : OMP_Set_Lock                     , OMP_Unset_Lock                     , omp_lock_kind                , OMP_Get_Max_Threads        , &
    !$   &                                                        OMP_Get_Thread_Num
    use            :: Locks                              , only : ompLock
    use            :: String_Handling                    , only : operator(//)
    use            :: Timers                             , only : timer
    use            :: Numerical_Constants_Prefixes       , only : siFormat
    use            :: Sorting                            , only : sortIndex
    implicit none
    class           (mergerTreeEvolverThreaded    )                       , intent(inout) :: self
    integer                                        , optional             , intent(  out) :: status
    type            (mergerTree                   )              , target , intent(inout) :: tree
    double precision                                                      , intent(in   ) :: timeEnd
    logical                                                               , intent(  out) :: treeDidEvolve                                  , suspendTree
    logical                                                               , intent(in   ) :: deadlockReporting
    integer         (kind_int8                    ), optional             , intent(in   ) :: systemClockMaximum
    integer         (omp_lock_kind                ), optional             , intent(inout) :: initializationLock
    type            (treeNode                     )              , pointer                :: nodeParent
    type            (treeNode                     ), save        , pointer                :: node                                           , nodeLock
    !$omp threadprivate(node,nodeLock)
    double precision                               , save                                 :: timeEndThisNode
    !$omp threadprivate(timeEndThisNode)
    type            (varying_string               ), save                                 :: vMessage
    !$omp threadprivate(vMessage)
    double precision                               , parameter                            :: largeTime                  =1.0d10
    procedure       (interruptTask                ), save        , pointer                :: interruptProcedure
    !$omp threadprivate(interruptProcedure)
    procedure       (timestepTask                 ), save        , pointer                :: timestepTask_
    !$omp threadprivate(timestepTask_)
    class           (*                            ), save        , pointer                :: timestepSelf
    !$omp threadprivate(timestepSelf)
    type            (enumerationVerbosityLevelType), parameter                            :: verbosityLevel             =verbosityLevelWarn
    class           (nodeComponentBasic           )              , pointer                :: basicParent
    class           (nodeComponentBasic           ), save        , pointer                :: basic
    !$omp threadprivate(basic)
    type            (mergerTree                   ), save        , pointer                :: currentTree
    !$omp threadprivate(currentTree)
    type            (nodeTask                     ), dimension(:), allocatable            :: evolveList                                    , evolveListTmp           , &
         &                                                                                   postProcessList                               , postProcessListTmp
    double precision                               , dimension(:), allocatable            :: evolveListCost                                , evolveListCostTmp
    integer         (c_size_t                     ), dimension(:), allocatable            :: evolveListOrder
    integer                                                                               :: countEvolve                                   , iEvolve                 , &
         &                                                                                   countPostProcess                              , iPostProcess
    type            (inputParameters              ), save        , allocatable            :: parameters
    !$omp threadprivate(parameters)
    type            (ompLock                      )                                       :: lockEvolveList                                , lockPostProcessList     , &
         &                                                                                   lockTree
    type            (mergerTreeWalkerAllNodes     )                                       :: treeWalker
    type            (enumerationDeadlockStatusType)                                       :: statusDeadlock
    integer                                                                               :: treeWalkCountPreviousOutput                   , nodesEvolvedCount       , &
         &                                                                                   nodesTotalCount                               , treeWalkCount           , &
         &                                                                                   countWorkers                                  , status_
    double precision                                                                      :: earliestTimeInTree                            , finalTimeInTree
    character       (len=24                       )                                       :: label
    character       (len=35                       )                                       :: message
    type            (varying_string               ), save                                 :: lockType       
    !$omp threadprivate(lockType)                               
    logical                                                                               :: anyTreeExistsAtOutputTime                     , hasIntertreeEvent       , &
         &                                                                                   didEvolve                                     , evolutionFailed         , &
         &                                                                                   evolutionExiting
    logical                                        , save                                 :: interrupted
    !$omp threadprivate(interrupted)
    integer                                        , save                                 :: evolutionPhase                                , evolutionPhaseMaximum
    !$omp threadprivate(evolutionPhase)
    type            (timer                        ), save                                 :: timer_
    !$omp threadprivate(timer_)
    double precision                                                                      :: timeBuildEvolveList                           , timeWaitEvolve          , &
         &                                                                                   timePostProcess                               , timeEvolve              , &
         &                                                                                   timeTotal                                     , timeInitialize          , &
         &                                                                                   timeCopyObjects                               , timeThreadInitialize    , &
         &                                                                                   timeThreadUninitialize                        , timeEvolveListLock      , &
         &                                                                                   timePostProcessListLock                       , timeBuildEvolveList__   , &
         &                                                                                   timePostProcess__
    double precision                               , save                                 :: timeWaitEvolve_                               , timeEvolve_             , &
         &                                                                                   timeThreadInitialize_                         , timeThreadUninitialize_ , &
         &                                                                                   timeEvolveListLock_                           , timePostProcessListLock_, &
         &                                                                                   timeBuildEvolveList_                          , timePostProcess_        , &
         &                                                                                   timeRemaining
    !$omp threadprivate(timeWaitEvolve_,timeEvolve_,timeThreadInitialize_,timeThreadUninitialize_,timeEvolveListLock_,timePostProcessListLock_,timeBuildEvolveList_,timePostProcess_,timeRemaining)
    character       (len=5                        )                                       :: percentBuildEvolveList                        , percentEvolve           , &
         &                                                                                   percentWaitEvolve                             , percentPostProcess      , &
         &                                                                                   percentPostProcessListLock                    , percentInitialize       , &
         &                                                                                   percentCopyObjects                            , percentThreadInitialize , &
         &                                                                                   percentThreadUninitialize                     , percentEvolveListLock

    ! Begin timing.
    if (self%reportTiming) then
       timeBuildEvolveList    =0.0d0
       timePostProcess        =0.0d0
       timeCopyObjects        =0.0d0
       timeThreadInitialize   =0.0d0
       timeThreadUninitialize =0.0d0
       timeEvolve             =0.0d0
       timeWaitEvolve         =0.0d0
       timeEvolveListLock     =0.0d0
       timePostProcessListLock=0.0d0
    end if
    ! Initialize trees.
    if (self%reportTiming) then
       timer_=timer()
       call timer_%start()
    end if
    suspendTree               =  .false.
    anyTreeExistsAtOutputTime =  .false.
    treeDidEvolve             =  .false.
    evolutionExiting          =  .false.
    call self%initializeTree(tree,timeEnd,treeDidEvolve,anyTreeExistsAtOutputTime,hasInterTreeEvent,initializationLock)
    if (.not.anyTreeExistsAtOutputTime.and..not.hasInterTreeEvent) then
       ! Mark the tree as evolved here, as the only reason that we did not evolve it was the given time target.
       treeDidEvolve=.true.
       return
    end if
    if (self%reportTiming) then
       call timer_%stop()
       timeInitialize=timer_%report()
       call timer_%start()
    end if
    ! Create copies of objects needed for evolution.
    countWorkers   =1
    !$ countWorkers=OMP_Get_Max_Threads()
    if (.not.self%workersInitialized) then
       allocate(self%workers(0:countWorkers-1))
       self%workersInitialized=.true.
       do numberWorker=0,countWorkers-1
          call eventsHooksFutureThread(numberWorker)
          allocate(self%workers(numberWorker)%cosmologyFunctions_      ,mold=self%cosmologyFunctions_      )
          allocate(self%workers(numberWorker)%mergerTreeEvolveTimestep_,mold=self%mergerTreeEvolveTimestep_)
          allocate(self%workers(numberWorker)%mergerTreeNodeEvolver_   ,mold=self%mergerTreeNodeEvolver_   )
          allocate(self%workers(numberWorker)%mergerTreeEvolveProfiler_,mold=self%mergerTreeEvolveProfiler_)
          allocate(self%workers(numberWorker)%galacticStructureSolver_ ,mold=self%galacticStructureSolver_ )
          allocate(self%workers(numberWorker)%metaTreeProcessingTime_  ,mold=self%metaTreeProcessingTime_  )
          !$omp critical(mergerTreeEvolverThreadedDeepCopy)
          !![
	  <deepCopyReset variables="self%cosmologyFunctions_ self%mergerTreeEvolveTimestep_ self%mergerTreeNodeEvolver_ self%mergerTreeEvolveProfiler_ self%galacticStructureSolver_ self%metaTreeProcessingTime_"/>
	  <deepCopy source="self%cosmologyFunctions_"       destination="self%workers(numberWorker)%cosmologyFunctions_"      />
	  <deepCopy source="self%mergerTreeEvolveTimestep_" destination="self%workers(numberWorker)%mergerTreeEvolveTimestep_"/>
	  <deepCopy source="self%mergerTreeNodeEvolver_"    destination="self%workers(numberWorker)%mergerTreeNodeEvolver_"   />
	  <deepCopy source="self%mergerTreeEvolveProfiler_" destination="self%workers(numberWorker)%mergerTreeEvolveProfiler_"/>
	  <deepCopy source="self%galacticStructureSolver_"  destination="self%workers(numberWorker)%galacticStructureSolver_" />
	  <deepCopy source="self%metaTreeProcessingTime_"   destination="self%workers(numberWorker)%metaTreeProcessingTime_"  />
	  <deepCopyFinalize variables="self%workers(numberWorker)%cosmologyFunctions_ self%workers(numberWorker)%mergerTreeEvolveTimestep_ self%workers(numberWorker)%mergerTreeNodeEvolver_ self%workers(numberWorker)%mergerTreeEvolveProfiler_ self%workers(numberWorker)%galacticStructureSolver_ self%workers(numberWorker)%metaTreeProcessingTime_"/>  
          !!]
          !$omp end critical(mergerTreeEvolverThreadedDeepCopy)
          call eventsHooksFutureThread()
        end do
    end if
    if (self%reportTiming) then
       call timer_%stop()
       timeCopyObjects=timer_%report()
    end if
    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible. Begin a parallel region
    ! here - this allows us to retain a team of threads for the entirety of this trees' evolution. We will use OpenMP `single`
    ! sections for work that must be done in serial.
    !$omp parallel
    ! Initialize nodes for these threads.
    if (self%reportTiming) then
       timer_=timer()
       call timer_%start()
    end if
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call parameters%parametersGroupCopy(self%parameters)
    call Node_Components_Thread_Initialize(parameters)   
    ! Determine our worker number.
    numberWorker=OMP_Get_Thread_Num()
    if (self%reportTiming) then
       call timer_%stop()
       timeThreadInitialize_=timer_%report()
       !$omp atomic
       timeThreadInitialize =max(                       &
            &                    timeThreadInitialize , &
            &                    timeThreadInitialize_  &
            &                   )
    end if
    !$omp single
    ! Initialize remaining time calculations.
    timeRemaining=self%metaTreeProcessingTime_%timeRemaining(tree,timeEnd)
    ! Initialize evolution state and locks.
    didEvolve                  =.true.
    treeWalkCount              =0
    treeWalkCountPreviousOutput=0
    lockEvolveList             =ompLock()
    lockPostProcessList        =ompLock()
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
          !$omp single
          if (statusDeadlock == deadlockStatusIsReporting) call displayIndent("Deadlock report follows")
          !$omp end single
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
                   if (self%reportTiming) call timer_%start()
                   evolutionPhase=evolutionPhase+1
                   !$omp single
                   countEvolve          =0
                   timeBuildEvolveList__=0.0d0
                   timePostProcess__    =0.0d0
                   treeWalker           =mergerTreeWalkerAllNodes(currentTree,spanForest=.false.)
                   do while (treeWalker%next(node))
                      ! Count nodes in the tree.
                      if (evolutionPhase == 1) nodesTotalCount=nodesTotalCount+1
                      ! Skip nodes that are not to be evolved during this phase.
                      if (.not.self%mergerTreeEvolveConcurrency_%includeInEvolution(evolutionPhase,node)) cycle
                      ! Get the basic component of the node.
                      basic => node%basic()
                      ! Determine if the node is evolvable.
                      if (self%nodeIsEvolvable(node,timeEnd,finalTimeInTree)) then                      
                         if (allocated(evolveList)) then
                            if (countEvolve == size(evolveList)) then
                               call move_alloc(evolveList    ,evolveListTmp    )
                               call move_alloc(evolveListCost,evolveListCostTmp)
                               deallocate(evolveListOrder               )
                               allocate  (evolveList     (countEvolve*2))
                               allocate  (evolveListCost (countEvolve*2))
                               allocate  (evolveListOrder(countEvolve*2))
                               evolveList    (1:countEvolve)=evolveListTmp
                               evolveListCost(1:countEvolve)=evolveListCostTmp
                               deallocate(evolveListTmp    )
                               deallocate(evolveListCostTmp)
                            end if
                         else
                            allocate(evolveList     (1))
                            allocate(evolveListCost (1))
                            allocate(evolveListOrder(1))
                         end if
                         countEvolve                               =  countEvolve+1
                         evolveList    (countEvolve)%node          => node
                         evolveList    (countEvolve)%timestepTask_ => null()
                         evolveListCost(countEvolve)               =  basic%floatRank0MetaPropertyGet(self%timeEvolveID)
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
                   evolveListOrder(1:countEvolve)=sortIndex(evolveListCost(1:countEvolve))
                   iEvolve                       =0
                   countPostProcess              =0
                   evolutionFailed               =.false.
                   evolutionPhaseMaximum         =self%mergerTreeEvolveConcurrency_%countPhases()
                   !$omp end single
                   if (self%reportTiming) then
                      call timer_%stop()
                      timeBuildEvolveList_ =timer_%report()
                      !$omp atomic
                      timeBuildEvolveList__=max(timeBuildEvolveList__,timeBuildEvolveList_)
                      !$omp barrier
                      !$omp single
                      timeBuildEvolveList=timeBuildEvolveList+timeBuildEvolveList__
                      !$omp end single
                   end if
                   ! Evolve all nodes on the stack.
                   do while (.not.evolutionFailed)
                      if (self%reportTiming) call timer_%start()
                      call lockEvolveList%set()
                      iEvolve=iEvolve+1
                      if (iEvolve <= countEvolve) then
                         node  => evolveList(evolveListOrder(countEvolve+1-iEvolve))%node
                         basic => node                                              %basic()
                         call lockEvolveList%unset()
                      else
                         call lockEvolveList%unset()
                         exit
                      end if
                      if (self%reportTiming) then
                         call timer_%stop ()
                         timeEvolveListLock_=+timer_             %report()
                         !$omp atomic
                         timeEvolveListLock = timeEvolveListLock           &
                              &              +timeEvolveListLock_
                      end if
                      call timer_%start()
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
                            !$omp critical (mergerTreeEvolverThreadNodeReport)
                            vMessage="node "
                            write (label,'(e12.6)') basic%time()
                            vMessage=vMessage//node%index()//" (current:target times = "//label
                            write (label,'(e12.6)') timeEnd
                            vMessage=vMessage//":"//label//")"
                            call displayIndent(vMessage)
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%workers(numberWorker)%cosmologyFunctions_,self%workers(numberWorker)%mergerTreeEvolveTimestep_,self%workers(numberWorker)%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.true. ,nodeLock=nodeLock,lockType=lockType)
                            call displayUnindent("end node")
                            !$omp end critical (mergerTreeEvolverThreadNodeReport)
                            call self%deadlockAddNode(node,currentTree%index,nodeLock,lockType)
                         else if (self%profileSteps) then
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%workers(numberWorker)%cosmologyFunctions_,self%workers(numberWorker)%mergerTreeEvolveTimestep_,self%workers(numberWorker)%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.false.                  ,lockType=lockType)
                            call self%workers(numberWorker)%mergerTreeEvolveProfiler_%stepDescriptor(lockType)
                         else
                            timeEndThisNode=self%timeEvolveTo(node,timeEnd,self%workers(numberWorker)%cosmologyFunctions_,self%workers(numberWorker)%mergerTreeEvolveTimestep_,self%workers(numberWorker)%mergerTreeNodeEvolver_,timestepTask_,timestepSelf,report=.false.                                    )
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
                            call lockTree%set  ()
                            call interruptProcedure(node)
                            call lockTree%unset()
                            ! Something happened so the tree is not deadlocked.
                            statusDeadlock=deadlockStatusIsNotDeadlocked
                         end if
                      end do
                      call timer_%stop ()
                      timeEvolve_=+timer_%report()
                      if (self%reportTiming) then                         
                         !$omp atomic
                         timeEvolve = timeEvolve          &
                              &      +timeEvolve_
                         call timer_%start()
                      end if
                      ! Add the node onto a postprocessing list, if it still exists.
                      if (associated(node)) then
                         call basic%floatRank0MetaPropertySet(self%timeEvolveID,timeEvolve_)
                         call lockPostProcessList%set()
                         if (allocated(postProcessList)) then
                            if (countPostProcess == size(postProcessList)) then
                               call move_alloc(postProcessList,postProcessListTmp)
                               allocate(postProcessList(size(postProcessListTmp)*2))
                               postProcessList(1:size(postProcessListTmp))=postProcessListTmp
                               deallocate(postProcessListTmp)
                            end if
                         else
                            allocate(postProcessList(1))
                         end if
                         countPostProcess                                =  countPostProcess+1
                         postProcessList(countPostProcess)%node          => node
                         postProcessList(countPostProcess)%timestepTask_ => timestepTask_
                         postProcessList(countPostProcess)%timestepSelf  => timestepSelf
                         call lockPostProcessList%unset()
                      end if
                      if (self%reportTiming) then
                         call timer_%stop ()
                         timePostProcessListLock_=+timer_                  %report()
                         !$omp atomic
                         timePostProcessListLock = timePostProcessListLock          &
                              &                   +timePostProcessListLock_
                      end if
                   end do
                   if (self%reportTiming) call timer_%start()
                   !$omp barrier
                   if (self%reportTiming) then
                      call timer_%stop ()
                      timeWaitEvolve_=+timer_          %report()
                      !$omp atomic
                      timeWaitEvolve = timeWaitEvolve            &
                           &          +timeWaitEvolve_
                      call timer_%start()
                   end if
                   !$omp single
                   if (evolutionFailed) then
                      ! Evolution failed - exit and return.
                      evolutionExiting=.true.
                   else
                      ! Perform post-processing of nodes.
                      do iPostProcess=1,countPostProcess
                         node                => postProcessList(iPostProcess)%node
                         timestepTask_       => postProcessList(iPostProcess)%timestepTask_
                         timestepSelf        => postProcessList(iPostProcess)%timestepSelf
                         ! Handle any end-of-timestep task associated with this node.
                         if (associated(timestepTask_))                                          &
                              & call timestepTask_(timestepSelf,currentTree,node,statusDeadlock)
                         ! Handle node promotion or merging.
                         if (associated(node)) then
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
                         end if
                      end do
                   end if
                   !$omp end single
                   if (self%reportTiming) then
                      call timer_%stop()
                      timePostProcess_ =timer_%report()
                      !$omp atomic
                      timePostProcess__=max(timePostProcess__,timePostProcess_)
                      !$omp barrier
                      !$omp single
                      timePostProcess=timePostProcess+timePostProcess__
                      !$omp end single
                   end if
                   if (evolutionExiting) exit
                end do concurrencyLoop
                !$omp barrier
                if (evolutionExiting) exit
                !$omp single
                ! Estimate remaining time to process the tree.
                timeRemaining=self%metaTreeProcessingTime_%timeRemaining(tree,timeEnd)
                if (timeRemaining > 0.0d0) then
                   write (label,'(i16)') int(timeRemaining)
                   call displayMessage("Estimated time remaining to process tree: "//trim(adjustl(label))//"s")
                end if
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
          call standardTreeEventsPerform(tree,statusDeadlock)
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
    if (self%reportTiming) call timer_%start()
    call Node_Components_Thread_Uninitialize()
    !$omp barrier
    !$omp critical(evolverReset)
    call parameters%reset()
    !$omp end critical(evolverReset)
    deallocate(parameters)
    if (self%reportTiming) then
       call timer_%stop()
       timeThreadUninitialize_=timer_%report()
       !$omp atomic
       timeThreadUninitialize =max(                         &
            &                      timeThreadUninitialize , &
            &                      timeThreadUninitialize_  &
            &                     )
    end if
    !$omp end parallel
    if (self%reportTiming) then
       timeEvolve             =+timeEvolve             /dble(countWorkers)
       timeEvolveListLock     =+timeEvolveListLock     /dble(countWorkers)
       timePostProcessListLock=+timePostProcessListLock/dble(countWorkers)
       timeWaitEvolve         =+timeWaitEvolve         /dble(countWorkers)
       timeTotal              =+timeInitialize          &
            &                  +timeCopyObjects         &
            &                  +timeThreadInitialize    &
            &                  +timeBuildEvolveList     &
            &                  +timeEvolveListLock      &
            &                  +timeEvolve              &
            &                  +timePostProcessListLock &
            &                  +timeWaitEvolve          &
            &                  +timePostProcess         &
            &                  +timeThreadUninitialize  
       write (percentInitialize         ,'(f5.1)') 100.0d0*timeInitialize         /timeTotal
       write (percentCopyObjects        ,'(f5.1)') 100.0d0*timeCopyObjects        /timeTotal
       write (percentThreadInitialize   ,'(f5.1)') 100.0d0*timeThreadInitialize   /timeTotal
       write (percentBuildEvolveList    ,'(f5.1)') 100.0d0*timeBuildEvolveList    /timeTotal
       write (percentEvolveListLock     ,'(f5.1)') 100.0d0*timeEvolveListLock     /timeTotal
       write (percentEvolve             ,'(f5.1)') 100.0d0*timeEvolve             /timeTotal
       write (percentPostProcessListLock,'(f5.1)') 100.0d0*timePostProcessListLock/timeTotal
       write (percentWaitEvolve         ,'(f5.1)') 100.0d0*timeWaitEvolve         /timeTotal
       write (percentPostProcess        ,'(f5.1)') 100.0d0*timePostProcess        /timeTotal
       write (percentThreadUninitialize ,'(f5.1)') 100.0d0*timeThreadUninitialize /timeTotal
       call displayIndent ('Timing report')
       call displayMessage('Total                       = '//trim(siFormat(timeTotal             ,'f10.6,1x'))//'s')
       call displayMessage('  ∦ Initialize tree        = '//trim(siFormat(timeInitialize         ,'f10.6,1x'))//'s ('//percentInitialize         //'%)')
       call displayMessage('  ∦ Copy objects           = '//trim(siFormat(timeCopyObjects        ,'f10.6,1x'))//'s ('//percentCopyObjects        //'%)')
       call displayMessage('  ∥ Initialize threads     = '//trim(siFormat(timeThreadInitialize   ,'f10.6,1x'))//'s ('//percentThreadInitialize   //'%)')
       call displayMessage('  ∦ Build evolve lists     = '//trim(siFormat(timeBuildEvolveList    ,'f10.6,1x'))//'s ('//percentBuildEvolveList    //'%)')
       call displayMessage('  ∦ Evolve array lock      = '//trim(siFormat(timeEvolveListLock     ,'f10.6,1x'))//'s ('//percentEvolveListLock     //'%)')
       call displayMessage('  ∥ Evolve nodes           = '//trim(siFormat(timeEvolve             ,'f10.6,1x'))//'s ('//percentEvolve             //'%)')
       call displayMessage('  ∦ Postprocess array lock = '//trim(siFormat(timePostProcessListLock,'f10.6,1x'))//'s ('//percentPostProcessListLock//'%)')
       call displayMessage('  ∥ Wait on threads        = '//trim(siFormat(timeWaitEvolve         ,'f10.6,1x'))//'s ('//percentWaitEvolve         //'%)')
       call displayMessage('  ∦ Postprocess nodes      = '//trim(siFormat(timePostProcess        ,'f10.6,1x'))//'s ('//percentPostProcess        //'%)')
       call displayMessage('  ∥ Uninitialize threads   = '//trim(siFormat(timeThreadUninitialize ,'f10.6,1x'))//'s ('//percentThreadUninitialize //'%)')
       call displayUnindent('done')
    end if
    return
  end subroutine threadedEvolve
