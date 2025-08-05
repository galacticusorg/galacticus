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

  use, intrinsic :: ISO_C_Binding                  , only : c_size_t
  use            :: Galacticus_Nodes               , only : mergerTree                 , treeNode                , universe
  use            :: Input_Parameters               , only : inputParameters
  use            :: Kind_Numbers                   , only : kind_int8
  use            :: Merger_Tree_Construction       , only : mergerTreeConstructorClass
  use            :: Merger_Tree_Initialization     , only : mergerTreeInitializorClass
  use            :: Merger_Tree_Operators          , only : mergerTreeOperatorClass
  use            :: Merger_Tree_Outputters         , only : mergerTreeOutputter        , mergerTreeOutputterClass
  use            :: Merger_Trees_Evolve            , only : mergerTreeEvolver          , mergerTreeEvolverClass
  use            :: Merger_Tree_Seeds              , only : mergerTreeSeedsClass
  use            :: Nodes_Operators                , only : nodeOperatorClass
  use            :: Numerical_Random_Numbers       , only : randomNumberGeneratorClass
  use            :: Output_Times                   , only : outputTimesClass
  use            :: Task_Evolve_Forests_Work_Shares, only : evolveForestsWorkShareClass
  use            :: Timers                         , only : timer
  use            :: Universe_Operators             , only : universeOperator           , universeOperatorClass

  !![
  <task name="taskEvolveForests">
   <description>A task which evolves galaxies within a set of merger tree forests.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskEvolveForests
     !!{
     Implementation of a task which evolves galaxies within a set of merger tree forests.
     !!}
     private
     ! Parameter controlling maximum number of forests to evolve.
     integer         (c_size_t )                            :: countForestsMaximum
     ! Parameter controlling maximum wall time for which forest evolution can run.
     integer         (kind_int8)                            :: walltimeMaximum
     ! Parameters controlling tree suspension.
     logical                                                :: suspendToRAM
     type            (varying_string             )          :: suspendPath
     ! Parameter controlling how threads process forests.
     logical                                                :: evolveForestsInParallel
     ! Tree universes used while processing all trees.
     type            (universe                   )          :: universeWaiting                         , universeProcessed
     ! Objects used in tree processing.
     class           (mergerTreeConstructorClass ), pointer :: mergerTreeConstructor_        => null()
     class           (mergerTreeOperatorClass    ), pointer :: mergerTreeOperator_           => null()
     class           (mergerTreeEvolverClass     ), pointer :: mergerTreeEvolver_            => null()
     class           (mergerTreeOutputterClass   ), pointer :: mergerTreeOutputter_          => null()
     class           (mergerTreeInitializorClass ), pointer :: mergerTreeInitializor_        => null()
     class           (nodeOperatorClass          ), pointer :: nodeOperator_                 => null()
     class           (evolveForestsWorkShareClass), pointer :: evolveForestsWorkShare_       => null()
     class           (outputTimesClass           ), pointer :: outputTimes_                  => null()
     class           (universeOperatorClass      ), pointer :: universeOperator_             => null()
     class           (randomNumberGeneratorClass ), pointer :: randomNumberGenerator_        => null()
     class           (mergerTreeSeedsClass       ), pointer :: mergerTreeSeeds_              => null()
     ! Pointer to the parameters for this task.
     type            (inputParameters            ), pointer :: parameters                    => null()
     logical                                                :: initialized                   =  .false., nodeComponentsInitialized=.false.
     ! Checkpointing.
     integer         (kind_int8                  )          :: timeIntervalCheckpoint
     type            (varying_string             )          :: fileNameCheckpoint
     type            (timer                      )          :: timer_
     ! Output time display format.
     integer                                                :: outputTimePrecision
     character       (len=9                      )          :: outputTimeFormat
   contains
     !![
     <methods>
       <method description="Suspend a tree (to memory or to file)." method="suspendTree"/>
       <method description="Restore a suspended tree."              method="resumeTree" />
     </methods>
     !!]
     final     ::                evolveForestsDestructor
     procedure :: perform     => evolveForestsPerform
     procedure :: suspendTree => evolveForestsSuspendTree
     procedure :: resumeTree  => evolveForestsResumeTree
     procedure :: autoHook    => evolveForestsAutoHook
  end type taskEvolveForests

  interface taskEvolveForests
     !!{
     Constructors for the \refClass{taskEvolveForests} task.
     !!}
     module procedure evolveForestsConstructorParameters
     module procedure evolveForestsConstructorInternal
  end interface taskEvolveForests

  ! Copies of objects used by each thread.
  class(mergerTreeOutputterClass  ), pointer :: mergerTreeOutputter_   => null()
  class(mergerTreeInitializorClass), pointer :: mergerTreeInitializor_ => null()
  class(mergerTreeEvolverClass    ), pointer :: mergerTreeEvolver_     => null()
  class(mergerTreeConstructorClass), pointer :: mergerTreeConstructor_ => null()
  class(mergerTreeOperatorClass   ), pointer :: mergerTreeOperator_    => null()
  class(nodeOperatorClass         ), pointer :: nodeOperator_          => null()
  !$omp threadprivate(mergerTreeOutputter_,mergerTreeInitializor_,mergerTreeEvolver_,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_)

contains

  function evolveForestsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskEvolveForests} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskEvolveForests          )                        :: self
    type            (inputParameters            ), intent(inout), target :: parameters
    class           (mergerTreeOperatorClass    ), pointer               :: mergerTreeOperator_
    class           (nodeOperatorClass          ), pointer               :: nodeOperator_
    class           (evolveForestsWorkShareClass), pointer               :: evolveForestsWorkShare_
    class           (mergerTreeConstructorClass ), pointer               :: mergerTreeConstructor_
    class           (outputTimesClass           ), pointer               :: outputTimes_
    class           (universeOperatorClass      ), pointer               :: universeOperator_
    class           (mergerTreeEvolverClass     ), pointer               :: mergerTreeEvolver_
    class           (mergerTreeOutputterClass   ), pointer               :: mergerTreeOutputter_
    class           (mergerTreeInitializorClass ), pointer               :: mergerTreeInitializor_
    class           (randomNumberGeneratorClass ), pointer               :: randomNumberGenerator_
    class           (mergerTreeSeedsClass       ), pointer               :: mergerTreeSeeds_
    type            (inputParameters            ), pointer               :: parametersRoot
    logical                                                              :: evolveForestsInParallel, suspendToRAM
    integer         (kind_int8                  )                        :: walltimeMaximum        , timeIntervalCheckpoint
    integer         (c_size_t                   )                        :: countForestsMaximum
    type            (varying_string             )                        :: suspendPath            , fileNameCheckpoint

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => null()
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    !![
    <inputParameter>
      <name>countForestsMaximum</name>
      <defaultValue>-1_c_size_t</defaultValue>
      <description>If set to a positive number, this is the maximum number of forests that will be evolved.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>walltimeMaximum</name>
      <defaultValue>-1_kind_int8</defaultValue>
      <description>If set to a positive number, this is the maximum wall time for which forest evolution is allowed to proceed before the task gives up.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>evolveForestsInParallel</name>
      <defaultValue>.true.</defaultValue>
      <description>If true then each forest is evolved by a separate OpenMP thread. Otherwise, a single thread evolves all forests.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>suspendToRAM</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether trees should be suspended to RAM (otherwise they are suspend to file).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (.not.suspendToRAM) then
       !![
       <inputParameter>
         <name>suspendPath</name>
         <description>The path to which tree suspension files will be stored.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>timeIntervalCheckpoint</name>
      <defaultValue>-1_kind_int8</defaultValue>
      <description>If positive, gives the time in seconds between storing of checkpoint files. If zero or negative, no checkpointing is performed..</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (timeIntervalCheckpoint> 0_kind_int8) then
       !![
       <inputParameter>
         <name>fileNameCheckpoint</name>
         <description>The path to which checkpoint data will be stored.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="mergerTreeConstructor"  name="mergerTreeConstructor_"  source="parameters"/>
    <objectBuilder class="mergerTreeOperator"     name="mergerTreeOperator_"     source="parameters"/>
    <objectBuilder class="nodeOperator"           name="nodeOperator_"           source="parameters"/>
    <objectBuilder class="evolveForestsWorkShare" name="evolveForestsWorkShare_" source="parameters"/>
    <objectBuilder class="outputTimes"            name="outputTimes_"            source="parameters"/>
    <objectBuilder class="universeOperator"       name="universeOperator_"       source="parameters"/>
    <objectBuilder class="mergerTreeEvolver"      name="mergerTreeEvolver_"      source="parameters"/>
    <objectBuilder class="mergerTreeOutputter"    name="mergerTreeOutputter_"    source="parameters"/>
    <objectBuilder class="mergerTreeInitializor"  name="mergerTreeInitializor_"  source="parameters"/>
    <objectBuilder class="randomNumberGenerator"  name="randomNumberGenerator_"  source="parameters"/>
    <objectBuilder class="mergerTreeSeeds"        name="mergerTreeSeeds_"        source="parameters"/>
    !!]
    if (associated(parametersRoot)) then
       self=taskEvolveForests(evolveForestsInParallel,countForestsMaximum,walltimeMaximum,suspendToRAM,suspendPath,timeIntervalCheckpoint,fileNameCheckpoint,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parametersRoot)
    else
       self=taskEvolveForests(evolveForestsInParallel,countForestsMaximum,walltimeMaximum,suspendToRAM,suspendPath,timeIntervalCheckpoint,fileNameCheckpoint,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parameters    )
    end if
    self%nodeComponentsInitialized=.true.
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeConstructor_" />
    <objectDestructor name="mergerTreeOperator_"    />
    <objectDestructor name="nodeOperator_"          />
    <objectDestructor name="evolveForestsWorkShare_"/>
    <objectDestructor name="outputTimes_"           />
    <objectDestructor name="universeOperator_"      />
    <objectDestructor name="mergerTreeEvolver_"     />
    <objectDestructor name="mergerTreeOutputter_"   />
    <objectDestructor name="mergerTreeInitializor_" />
    <objectDestructor name="randomNumberGenerator_" />
    <objectDestructor name="mergerTreeSeeds_"       />
    !!]
    return
  end function evolveForestsConstructorParameters

  function evolveForestsConstructorInternal(evolveForestsInParallel,countForestsMaximum,walltimeMaximum,suspendToRAM,suspendPath,timeIntervalCheckpoint,fileNameCheckpoint,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parameters) result(self)
    !!{
    Internal constructor for the \refClass{taskEvolveForests} task class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    type            (taskEvolveForests          )                        :: self
    logical                                      , intent(in   )         :: evolveForestsInParallel, suspendToRAM
    integer         (kind_int8                  ), intent(in   )         :: walltimeMaximum        , timeIntervalCheckpoint
    integer         (c_size_t                   ), intent(in   )         :: countForestsMaximum
    type            (varying_string             ), intent(in   )         :: suspendPath            , fileNameCheckpoint
    class           (mergerTreeConstructorClass ), intent(in   ), target :: mergerTreeConstructor_
    class           (mergerTreeOperatorClass    ), intent(in   ), target :: mergerTreeOperator_
    class           (nodeOperatorClass          ), intent(in   ), target :: nodeOperator_
    class           (evolveForestsWorkShareClass), intent(in   ), target :: evolveForestsWorkShare_
    class           (outputTimesClass           ), intent(in   ), target :: outputTimes_
    class           (universeOperatorClass      ), intent(in   ), target :: universeOperator_
    class           (mergerTreeEvolverClass     ), intent(in   ), target :: mergerTreeEvolver_
    class           (mergerTreeOutputterClass   ), intent(in   ), target :: mergerTreeOutputter_
    class           (mergerTreeInitializorClass ), intent(in   ), target :: mergerTreeInitializor_
    class           (randomNumberGeneratorClass ), intent(in   ), target :: randomNumberGenerator_
    class           (mergerTreeSeedsClass       ), intent(in   ), target :: mergerTreeSeeds_
    type            (inputParameters            ), intent(in   ), target :: parameters
    integer         (c_size_t                   )                        :: i
    double precision                                                     :: timeStepMinimum
    !![
    <constructorAssign variables="evolveForestsInParallel, countForestsMaximum, walltimeMaximum, suspendToRAM, suspendPath, timeIntervalCheckpoint, fileNameCheckpoint, *mergerTreeConstructor_, *mergerTreeOperator_, *nodeOperator_, *evolveForestsWorkShare_, *outputTimes_, *universeOperator_, *mergerTreeEvolver_, *mergerTreeOutputter_, *mergerTreeInitializor_, *randomNumberGenerator_, *mergerTreeSeeds_"/>
    !!]

    self%parameters  => parameters
    self%initialized =  .true.
    self%timer_      = timer()
    ! Validate.
    if (evolveForestsInParallel .and. timeIntervalCheckpoint > 0_kind_int8) call Error_Report('Checkpointing is not possible when evolving forests in parallel'//{introspection:location})
    ! Find the minimum step in output times and compute the precision for outputting times such that all are distinct.
    if (self%outputTimes_%count() > 0) then
       timeStepMinimum=self%outputTimes_%time(1_c_size_t)
       if (self%outputTimes_%count() > 1) then
          do i=2,self%outputTimes_%count()
             timeStepMinimum=min(self%outputTimes_%time(i)-self%outputTimes_%time(i-1),timeStepMinimum)
          end do
       end if
       self%outputTimePrecision=max(2,-floor(log10(timeStepMinimum)))
       write (self%outputTimeFormat,'(a2,i2.2,a1,i2.2,a1)') "(f",self%outputTimePrecision+2+max(0,floor(log10(self%outputTimes_%time(self%outputTimes_%count())))),".",self%outputTimePrecision,")"
    else
       self%outputTimeFormat="(f)"
    end if
    return 
  end function evolveForestsConstructorInternal

  subroutine evolveForestsAutoHook(self)
    !!{
    Attach to state store/restore event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingNone, stateRestoreEventGlobal, stateStoreEventGlobal
    implicit none
    class(taskEvolveForests), intent(inout) :: self

    call stateStoreEventGlobal  %attach(self,evolveForestsStateStore  ,openMPThreadBindingNone,label='evolveForests')
    call stateRestoreEventGlobal%attach(self,evolveForestsStateRestore,openMPThreadBindingNone,label='evolveForests')
    return
  end subroutine evolveForestsAutoHook

  subroutine evolveForestsStateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Store the internal state of this object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_ptr, c_size_t
    implicit none
    class  (taskEvolveForests), intent(inout) :: self
    integer                   , intent(in   ) :: stateFile
    type   (c_ptr            ), intent(in   ) :: gslStateFile
    integer(c_size_t         ), intent(in   ) :: stateOperationID
    !$GLC attributes unused :: self

    call mergerTreeEvolver_    %stateStore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeOutputter_  %stateStore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeInitializor_%stateStore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeConstructor_%stateStore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeOperator_   %stateStore(stateFile,gslStateFile,stateOperationID)
    call nodeOperator_         %stateStore(stateFile,gslStateFile,stateOperationID)
    return
  end subroutine evolveForestsStateStore

  subroutine evolveForestsStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Store the internal state of this object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_ptr, c_size_t
    implicit none
    class  (taskEvolveForests), intent(inout) :: self
    integer                   , intent(in   ) :: stateFile
    type   (c_ptr            ), intent(in   ) :: gslStateFile
    integer(c_size_t         ), intent(in   ) :: stateOperationID
    !$GLC attributes unused :: self

    call mergerTreeEvolver_    %stateRestore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeOutputter_  %stateRestore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeInitializor_%stateRestore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeConstructor_%stateRestore(stateFile,gslStateFile,stateOperationID)
    call mergerTreeOperator_   %stateRestore(stateFile,gslStateFile,stateOperationID)
    call nodeOperator_         %stateRestore(stateFile,gslStateFile,stateOperationID)
    return
  end subroutine evolveForestsStateRestore

  subroutine evolveForestsCheckpoint(self,node)
    !!{
    Checkpoint the current tree.
    !!}
    use :: File_Utilities          , only : File_Rename
    use :: Dates_and_Times         , only : Formatted_Date_and_Time
    use :: Display                 , only : displayMessage         , displayIndent, displayUnindent, verbosityLevelWorking
    use :: Galacticus_Nodes        , only : treeNode
    use :: Merger_Tree_Construction, only : mergerTreeStateStore
    use :: ISO_Varying_String      , only : operator(//)
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    class is (taskEvolveForests)
       call self%timer_%stop()
       if (int(self%timer_%report(),kind=kind_int8) >= self%timeIntervalCheckpoint) then
          call displayIndent("Checkpointing",verbosityLevelWorking)
          call displayMessage("Begin at "    //Formatted_Date_and_Time(),verbosityLevelWorking)
          call self%timer_%start()
          call mergerTreeStateStore(node%hostTree,char(self%fileNameCheckpoint)//'.tmp',snapshot=.false.,append=.false.)
          call File_Rename(self%fileNameCheckpoint//'.tmp',self%fileNameCheckpoint,overwrite=.true.)
          call displayMessage("Completed at "//Formatted_Date_and_Time(),verbosityLevelWorking)
          call displayUnindent("done",verbosityLevelWorking)
       end if
    end select
    return
  end subroutine evolveForestsCheckpoint

  subroutine evolveForestsDestructor(self)
    !!{
    Destructor for the \refClass{taskEvolveForests} task class.
    !!}
    use :: Events_Hooks    , only : stateRestoreEventGlobal     , stateStoreEventGlobal
    use :: Node_Components , only : Node_Components_Uninitialize
    use :: Galacticus_Nodes, only : nodeClassHierarchyFinalize
    implicit none
    type(taskEvolveForests), intent(inout) :: self

    if (.not.self%initialized) return
    !![
    <objectDestructor name="self%mergerTreeConstructor_" />
    <objectDestructor name="self%mergerTreeOperator_"    />
    <objectDestructor name="self%nodeOperator_"          />
    <objectDestructor name="self%evolveForestsWorkShare_"/>
    <objectDestructor name="self%outputTimes_"           />
    <objectDestructor name="self%universeOperator_"      />
    <objectDestructor name="self%mergerTreeEvolver_"     />
    <objectDestructor name="self%mergerTreeOutputter_"   />
    <objectDestructor name="self%mergerTreeInitializor_" />
    <objectDestructor name="self%randomNumberGenerator_" />
    <objectDestructor name="self%mergerTreeSeeds_"       />
    !!]
    if (stateStoreEventGlobal  %isAttached(self,evolveForestsStateStore  )) call stateStoreEventGlobal  %detach(self,evolveForestsStateStore  )
    if (stateRestoreEventGlobal%isAttached(self,evolveForestsStateRestore)) call stateRestoreEventGlobal%detach(self,evolveForestsStateRestore)
    if (self%nodeComponentsInitialized                                    ) then
       call Node_Components_Uninitialize()
       call nodeClassHierarchyFinalize  ()
    end if
    return
  end subroutine evolveForestsDestructor

  subroutine evolveForestsPerform(self,status)
    !!{
    Evolves the complete set of merger trees as specified.
    !!}
    use               :: Display                 , only : displayIndent                    , displayMessage                     , displayUnindent, verbosityLevelInfo
    use               :: Error                   , only : Error_Report                     , errorStatusSuccess
    use               :: Events_Hooks            , only : openMPThreadBindingAllLevels     , postEvolveEvent
    use               :: File_Utilities          , only : File_Exists
    use               :: Galacticus_Nodes        , only : mergerTree                       , nodeComponentBasic                 , treeNode       , universe          , &
          &                                               universeEvent
    use   , intrinsic :: ISO_C_Binding           , only : c_size_t
    use               :: Memory_Reporting        , only : reportMemoryUsage
    use               :: Merger_Tree_Construction, only : mergerTreeStateFromFile
    use               :: Merger_Tree_Walkers     , only : mergerTreeWalkerAllNodes
    use               :: Node_Components         , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use               :: Node_Events_Inter_Tree  , only : Inter_Tree_Event_Post_Evolve
    !$ use            :: OMP_Lib                 , only : OMP_Destroy_Lock                 , OMP_Get_Thread_Num                 , OMP_Init_Lock  , omp_lock_kind
    use               :: Sorting                 , only : sortIndex
    use               :: String_Handling         , only : operator(//)
    ! Include modules needed for tasks.
    !![
    <include directive="universePostEvolveTask" type="moduleUse" functionType="void">
    !!]
    include 'tasks.evolve_tree.universePostEvolveTask.moduleUse.inc'
    !![
    </include>
    !!]
    implicit none
    class           (taskEvolveForests       ), intent(inout), target           :: self
    integer                                   , intent(  out), optional         :: status
    type            (mergerTree              ), pointer                  , save :: tree
    logical                                                              , save :: finished             , treeIsNew
    integer         (c_size_t                )                           , save :: iOutput
    double precision                                                     , save :: evolveToTime         , treeTimeEarliest       , &
         &                                                                         universalEvolveToTime, treeTimeLatest         , &
         &                                                                         outputTimeNext
    type            (varying_string          )                           , save :: message
    character       (len=20                  )                           , save :: label
    !$omp threadprivate(tree,finished,iOutput,evolveToTime,message,label,treeIsNew,treeTimeEarliest,treeTimeLatest,outputTimeNext)
    logical                                                              , save :: treeIsFinished       , evolutionIsEventLimited, &
         &                                                                         success              , removeTree             , &
         &                                                                         suspendTree          , treesDidEvolve         , &
         &                                                                         treeDidEvolve        , treesCouldEvolve       , &
         &                                                                         deadlockReport
    type            (mergerTree              ), pointer                  , save :: currentTree          , previousTree           , &
         &                                                                         nextTree
    type            (mergerTreeWalkerAllNodes)                           , save :: treeWalkerAll
    !$omp threadprivate(currentTree,previousTree,nextTree,treeWalkerAll)
    type            (treeNode                ), pointer                  , save :: satelliteNode
    class           (nodeComponentBasic      ), pointer                  , save :: basicNodeBase
    !$omp threadprivate(satelliteNode,basicNodeBase,treeIsFinished,evolutionIsEventLimited,success,removeTree,suspendTree,treeDidEvolve)
    type            (universeEvent           ), pointer                  , save :: event_
    !$omp threadprivate(event_)
    type            (treeNode                ), pointer                  , save :: node
    class           (nodeComponentBasic      ), pointer                  , save :: basic
    logical                                                              , save :: treesFinished
    integer         (c_size_t                )                           , save :: treeNumber
    type            (inputParameters         ), allocatable              , save :: parameters
    integer         (c_size_t                )                                  :: treeCount
    integer         (omp_lock_kind           )                                  :: initializationLock
    integer         (kind_int8               )                                  :: systemClockRate      , systemClockMaximum
    !$omp threadprivate(node,basic,treeNumber,treesFinished,parameters)
    logical                                                                     :: checkpointRestored   , checkpointing         , &
         &                                                                         universeUpdated

    ! The following processes merger trees, one at a time, to each successive output time, then dumps their contents to file. It
    ! allows for the possibility of "universal events" - events which require all merger trees to reach the same cosmic time. If
    ! such an event exists, each tree is processed up to that time and then pushed onto a stack where it waits to be
    ! processed. Once all trees reach the event time the stack of trees is passed to the event task. Tree processing then
    ! continues by popping trees off of the stack and processing them further (possibly to the next universal event).

    call displayIndent('Begin task: merger tree evolution')

    ! Set status to success by default.
    if (present(status)) status=errorStatusSuccess

    ! Initialize a lock used for controlling tree initialization.
    !$ call OMP_Init_Lock(initializationLock)
    ! Initialize checkpoint restoration state.
    checkpointing               =.true.
    checkpointRestored          =.false.
    ! Initialize tree counter and record that we are not finished processing trees.
    deadlockReport              =.false.
    finished                    =.false.
    treeCount                   =0_c_size_t
    ! Initialize universes which will act as tree stacks. We use two stacks: one for trees waiting to be processed, one for trees
    ! that have already been processed.
    self%universeWaiting  =universe()
    self%universeProcessed=universe()

    ! Set record of whether any trees were evolved to false initially.
    treesDidEvolve  =.false.
    treesCouldEvolve=.false.

    ! Set the maximum allowed wall time.
    if (self%wallTimeMaximum > 0_kind_int8) then
       call System_Clock(systemClockMaximum,systemClockRate)
       systemClockMaximum=systemClockMaximum+self%wallTimeMaximum*systemClockRate
    else
       systemClockMaximum=0_kind_int8
    end if

    ! Begin parallel processing of trees until all work is done.
    !$omp parallel copyin(finished) if (self%evolveForestsInParallel)
    allocate(mergerTreeOutputter_  ,mold=self%mergerTreeOutputter_  )
    allocate(mergerTreeInitializor_,mold=self%mergerTreeInitializor_)
    allocate(mergerTreeEvolver_    ,mold=self%mergerTreeEvolver_    )
    allocate(mergerTreeConstructor_,mold=self%mergerTreeConstructor_)
    allocate(mergerTreeOperator_   ,mold=self%mergerTreeOperator_   )
    allocate(nodeOperator_         ,mold=self%nodeOperator_         )
    !$omp critical(evolveForestsDeepCopy)
    !![
    <deepCopyReset variables="self%mergerTreeEvolver_ self%mergerTreeOutputter_ self%mergerTreeInitializor_ self%mergerTreeConstructor_ self%mergerTreeOperator_ self%nodeOperator_"/>
    <deepCopy source="self%mergerTreeEvolver_"     destination="mergerTreeEvolver_"    />
    <deepCopy source="self%mergerTreeOutputter_"   destination="mergerTreeOutputter_"  />
    <deepCopy source="self%mergerTreeInitializor_" destination="mergerTreeInitializor_"/>
    <deepCopy source="self%mergerTreeConstructor_" destination="mergerTreeConstructor_"/>
    <deepCopy source="self%mergerTreeOperator_"    destination="mergerTreeOperator_"   />
    <deepCopy source="self%nodeOperator_"          destination="nodeOperator_"         />
    <deepCopyFinalize variables="mergerTreeEvolver_ mergerTreeOutputter_ mergerTreeInitializor_ mergerTreeConstructor_ mergerTreeOperator_ nodeOperator_"/>
    !!]
    !$omp end critical(evolveForestsDeepCopy)
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call Node_Components_Thread_Initialize(parameters)
    ! Allow events to be attached to the universe.
    !$omp master
    self%universeWaiting%event => null()
    !![
    <eventHook name="universePreEvolve">
     <callWith>self%universeWaiting</callWith>
    </eventHook>
    !!]
    call self%universeOperator_%operate(self%universeWaiting)
    if (self%timeIntervalCheckpoint > 0_kind_int8) then
       call postEvolveEvent%attach(self,evolveForestsCheckpoint,openMPThreadBindingAllLevels,label='evolveForests')
       call self%timer_    %start (                                                                               )
    end if
    !$omp end master
    !$omp barrier
    ! Begin processing trees.
    treeProcess : do while (.not.finished)
       ! Attempt to get a new tree to process. We first try to get a new tree. If no new trees exist, we will look for a tree on
       ! the stack waiting to be processed.
       ! Perform any pre-tree construction tasks.
       call mergerTreeOperator_%operatePreConstruction()
       ! Get a tree.
       if (checkpointing .and. self%timeIntervalCheckpoint > 0 .and. File_Exists(self%fileNameCheckpoint)) then
          ! Resume from a checkpointed tree.
          if (checkpointRestored) then
             tree => null()
          else
             allocate(tree)
             call mergerTreeStateFromFile(tree,char(self%fileNameCheckpoint),self%randomNumberGenerator_,self%mergerTreeSeeds_,deleteAfterRead=.false.)
             checkpointRestored=.true.
          end if
       else
          checkpointing =  .false.
          treesFinished =  .false.
          tree          => null()
          do while (.not.associated(tree).and..not.treesFinished)
             ! Get the number of the next tree to process.
             treeNumber =  self                  %evolveForestsWorkShare_%forestNumber(utilizeOpenMPThreads=self%evolveForestsInParallel)
             tree       => mergerTreeConstructor_                        %construct   (treeNumber,treesFinished)
          end do
       end if
       if (associated(tree)) then
          tree%hostUniverse => self%universeWaiting
          ! Limit to the maximum number of forests allowed to be run.
          if (self%countForestsMaximum >= 0_c_size_t .and. treeNumber > self%countForestsMaximum) then
             call tree%destroy()
             deallocate(tree)
             tree => null()
          end if
       end if
       finished                                =  finished.or..not.associated(tree)
       treeIsNew                               =  .not.finished.and..not.checkpointRestored       
       ! If no new tree was available, attempt to pop one off the universe stack.
       if (finished) then
          call self%resumeTree(tree)
          treeIsNew=.false.
          finished =.not.associated(tree)
       end if
       ! Report on memory utilization.
       call reportMemoryUsage()
       ! If we got a tree (i.e. we are not "finished") process it.
       if (.not.finished) then
          treeIsFinished=.false.
          ! Count trees.
          !$omp atomic
          treeCount=treeCount+1_c_size_t
          if (treeCount > 1_c_size_t .and. self%timeIntervalCheckpoint > 0_kind_int8) call Error_Report('more than 1 tree not permitted when checkpointing is enabled'//{introspection:location})
          ! If this is a new tree, perform any initialization and pre-evolution tasks on it.
          if (treeIsNew) then
             ! Walk over all nodes and perform "node tree" initialization. This typically includes initialization related to
             ! the static structure of the tree (e.g. assign scale radii, merging orbits, etc.). Initialization related to
             ! evolution of the tree (e.g. growth rates of scale radii, baryonic component initialization) are typically handled
             ! by the mergerTreeInitializor class which is called later.
             call    mergerTreeOperator_%operatePreInitialization(tree)
             treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
             do while (treeWalkerAll%next(node))
                call nodeOperator_      %nodeTreeInitialize (node)
             end do
             call    mergerTreeOperator_%operatePreEvolution(tree)
             message="Evolving tree number "
          else
             message="Resuming tree number "
          end if
          ! Display a message.
          message=message//tree%index//" {"//tree%nodeBase%index()//"}"
          call displayIndent(message)
          ! Get the next time to which the tree should be evolved.
          treeTimeEarliest=tree%earliestTime()
          outputTimeNext  =self%outputTimes_%timeNext(treeTimeEarliest,indexOutput=iOutput)
          ! For new trees, if the earliest time in the tree exactly coincides with an output
          ! time, then process the tree to that output. This ensures that we include all
          ! halos from this time in the output, even though they will be devoid of any
          ! galaxies.
          if (treeIsNew .and. iOutput > 1) then
             if (treeTimeEarliest == self%outputTimes_%time(iOutput-1)) iOutput=iOutput-1
          end if
          ! Resumed trees must always be allowed to evolve - since they by definition were not finished (otherwise they would
          ! not have been suspended).
          if (.not.treeIsNew .and. iOutput > self%outputTimes_%count()) iOutput=self%outputTimes_%count()
          ! Catch cases where there is no next output time.
          if (outputTimeNext > 0.0d0) then
             treeIsFinished=.false.
          else
             iOutput       =self%outputTimes_%count()+1
             treeIsFinished=.true.
          end if
          ! Iterate evolving the tree until no more outputs are required.
          treeEvolveLoop : do while (iOutput <= self%outputTimes_%count() .and. associated(tree))
             ! We want to find the maximum time to which we can evolve this tree. This will be the minimum of the next output
             ! time (at which we must stop and output the tree) and the next universal event time (at which we must stop and
             ! perform the event task). Find the next output time.
             evolveToTime=self%outputTimes_%time(iOutput)
             ! Find the earliest universe event.
             call self%universeWaiting%lock%set()
             event_                  => self%universeWaiting%event
             evolutionIsEventLimited =  .false.
             do while (associated(event_))
                if (event_%time < evolveToTime) then
                   evolveToTime           =event_%time
                   evolutionIsEventLimited=.true.
                   universalEvolveToTime  =evolveToTime
                end if
                event_ => event_%next
             end do
             call self%universeWaiting%lock%unset()
             if (tree%earliestTime() <= evolveToTime) treesCouldEvolve=.true.
             ! Evolve the tree to the computed time.
             call mergerTreeEvolver_%evolve(tree,evolveToTime,treeDidEvolve,suspendTree,deadlockReport,systemClockMaximum,status=status)
             if (present(status)) then
                if (status /= errorStatusSuccess) then
                   ! Tree evolution failed - abort further evolution and return the failure code.
                   treeIsFinished=.true.
                   finished      =.true.
                   exit
                end if
             end if
             !$omp critical (universeStatus)
             ! Record that evolution of the universe of trees occurred if this tree evolved and was not suspended.
             if (treeDidEvolve) treesDidEvolve=.true.
             !$omp end critical (universeStatus)
             ! If tree was marked to be suspended, record that evolution is limited by this suspension event.
             if (suspendTree) evolutionIsEventLimited=.true.
             ! Locate trees which consist of only a base node with no progenitors. These have reached
             ! the end of their evolution, and can be removed from the forest of trees.
             previousTree => null()
             currentTree  => tree
             do while (associated(currentTree))
                ! Skip empty trees.
                removeTree=.false.
                if (associated(currentTree%nodeBase)) then
                   basicNodeBase => currentTree%nodeBase%basic()
                   removeTree    =   .not.associated(currentTree%nodeBase%firstChild) &
                        &           .and.                                             &
                        &            (basicNodeBase%time() < evolveToTime)
                   if (removeTree) then
                      ! Does the node have attached satellites which are about to merge?
                      satelliteNode => currentTree%nodeBase%firstSatellite
                      do while (associated(satelliteNode))
                         if (associated(satelliteNode%mergeTarget)) then
                            removeTree=.false.
                            exit
                         end if
                         satelliteNode => satelliteNode%sibling
                      end do
                      ! Does the node have attached events?
                      if (associated(currentTree%nodeBase%event)) removeTree=.false.
                   end if
                else
                   ! No need to remove already empty trees.
                   removeTree=.false.
                end if
                if (removeTree) then
                   message="Removing remnant tree "
                   message=message//currentTree%index//" {"//currentTree%nodeBase%index()//"}"
                   call displayMessage(message,verbosityLevelInfo)
                   if (.not.associated(previousTree)) then
                      ! No previous tree is set, so the current tree is the first tree in the forest. Destroy it, but reset the
                      ! `tree` pointer to the next tree in the forest (otherwise `tree` will become a dangling pointer).
                      nextTree    => currentTree%nextTree
                      call currentTree%destroy()
                      deallocate(currentTree)
                      currentTree => nextTree
                      tree        => currentTree
                   else
                      ! A previous tree exists - simply connect its `next` pointer to the `next` pointer of the current tree, and
                      ! then destroy the current tree.
                      previousTree%nextTree => currentTree%nextTree
                      call currentTree%destroy()
                      deallocate(currentTree)
                      currentTree => previousTree%nextTree
                   end if
                else
                   previousTree => currentTree
                   currentTree  => currentTree%nextTree
                end if
             end do
             ! Check that tree reached required time. If it did not, we can evolve it no further.
             if (associated(tree)) then
                treeTimeEarliest=tree%earliestTimeEvolving()
                treeTimeLatest  =tree%  latestTime        ()
                if     (                                 &
                     &   treeTimeLatest   > evolveToTime &
                     &  .and.                            &
                     &   treeTimeEarliest < evolveToTime &
                     &  .and.                            &
                     &   .not.evolutionIsEventLimited    &
                     & ) then
                   if (deadlockReport) exit
                   message='failed to evolve tree to required time'//char(10)
                   write (label,'(f7.2)') evolveToTime
                   message=message//"            target time = "//trim(label)//" Gyr"//char(10)
                   write (label,'(f7.2)') treeTimeEarliest
                   message=message//"  earliest time in tree = "//trim(label)//" Gyr"//char(10)
                   write (label,'(f7.2)') treeTimeLatest
                   message=message//"    latest time in tree = "//trim(label)//" Gyr"
                   call Error_Report(message//{introspection:location})
                end if
             end if
             ! Determine what limited evolution.
             if (evolutionIsEventLimited) then
                ! Tree evolution was limited by a universal event. Therefore it can evolve no further
                ! until that event's task is performed.
                exit
             else
                ! Report on memory utilization.
                call reportMemoryUsage()
                ! Tree reached an output time, so output it. We can then continue evolving.
                write (label,self%outputTimeFormat) evolveToTime
                message="Output tree data at t="//trim(label)//" Gyr"
                call displayMessage(message)
                if (associated(tree)) then
                   call mergerTreeOutputter_%outputTree(tree,iOutput,evolveToTime)
                   ! Perform any extra output and post-output processing on nodes.
                   treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
                   do while (treeWalkerAll%next(node))
                      basic => node%basic()
                      if (basic%time() == evolveToTime) call node%postOutput(evolveToTime)
                   end do
                end if
                iOutput=iOutput+1
                ! If all output times have been reached, or all trees have been destroyed, we're finished.
                if (iOutput > self%outputTimes_%count() .or. .not.associated(tree)) then
                   treeIsFinished=.true.
                   exit
                end if
             end if
          end do treeEvolveLoop
          ! If tree could not evolve further, but is not finished, push it to the universe stack.
          if (.not.treeIsFinished) then
             ! Suspend the tree.
             call self%suspendTree(tree)
             ! Unindent messages.
             call displayUnindent('Suspending tree')
          else
             ! Unindent messages.
             call displayUnindent('Finished tree'  )
          end if
       end if
       ! Destroy the tree.
       if (associated(tree)) then
          currentTree => tree
          do while (associated(currentTree))
             previousTree => currentTree
             currentTree  => currentTree%nextTree
             call previousTree%destroy()
             ! Deallocate the tree.
             deallocate(previousTree)
          end do
          nullify(tree)
       end if
       ! Perform any post-evolution operations on the tree.
       if (treeIsFinished) call mergerTreeOperator_%operatePostEvolution()
       ! If any trees were pushed onto the processed stack, then there must be an event to process.
       if (finished) then
          !$omp barrier
          ! The following section - which transfers processed trees back to the active universe - is done within an OpenMP
          ! master section. We choose not to use an OpenMP single section here (which would be marginally more efficient) since
          ! single sections have an implicit barrier at the end, which would desynchronize threads when multiple threads are
          ! used to process a single tree. Using a master section avoids that implicit barrier - we instead handle the barrier
          ! (and sharing of the "finished" status back to other threads) explicitly if needed.
          !$omp master
          ! Check whether any tree evolution occurred. If it did not, we have a universe-level deadlock.
          if ((.not.treesDidEvolve.and.treesCouldEvolve).or.deadlockReport) then
             ! If we already did the deadlock reporting pass it's now time to finish that report and exit. Otherwise, set deadlock
             ! reporting status to true and continue for one more pass through the universe.
             if (deadlockReport) then
                message="Universe appears to be deadlocked"//char(10)
                call self%universeProcessed%lock%set()
                treeCount=0_c_size_t
                if (associated(self%universeProcessed%trees)) then
                   do while (associated(self%universeProcessed%trees))
                      tree => self%universeProcessed%popTree()
                      treeCount=treeCount+1_c_size_t
                   end do
                   message=message//" --> There are "//treeCount//" trees pending further processing"
                else
                   message=message//" --> There are no trees pending further processing"
                end if
                call self%universeProcessed%lock%unset()
                call displayMessage(message)
                call Inter_Tree_Event_Post_Evolve()
                call Error_Report('exiting'//{introspection:location})
             else
                deadlockReport=.true.
             end if
          end if
          treesDidEvolve=.false.
          call self%universeWaiting  %lock%set()
          call self%universeProcessed%lock%set()
          universeUpdated=.false.
          if (associated(self%universeProcessed%trees)) then
             ! Transfer processed trees back to the waiting universe.
             self%universeWaiting  %trees => self%universeProcessed%trees
             self%universeProcessed%trees => null()
             ! Find the event to process.
             event_ => self%universeWaiting%event
             do while (associated(event_))
                if (event_%time < universalEvolveToTime) then
                   call Error_Report('a universal event exists in the past - this should not happen'//{introspection:location})
                else if (event_%time == universalEvolveToTime) then
                   universeUpdated=.true.
                   success        =event_%task(self%universeWaiting)
                   if (success) call self%universeWaiting%removeEvent(event_)
                   exit
                end if
                event_ => event_%next
             end do
             call displayMessage('Finished universe evolution pass')
          end if
          call self%universeWaiting  %lock%unset()
          call self%universeProcessed%lock%unset()
          !$omp end master
          !$omp barrier
          if (universeUpdated) finished=.false.
          !$omp barrier
       end if
    end do treeProcess
    ! Finalize any merger tree operator.
    call mergerTreeOperator_ %finalize(                         )
    ! Reduce outputs back into the original outputter object.
    call mergerTreeOutputter_%reduce  (self%mergerTreeOutputter_)
    ! Explicitly deallocate objects.
    !![
    <objectDestructor name="mergerTreeOutputter_"  />
    <objectDestructor name="mergerTreeInitializor_"/>
    <objectDestructor name="mergerTreeEvolver_"    />
    <objectDestructor name="mergerTreeConstructor_"/>
    <objectDestructor name="mergerTreeOperator_"   />
    <objectDestructor name="nodeOperator_"         />
    !!]
    call Node_Components_Thread_Uninitialize()
    !$omp barrier
    !$omp critical(evolveForestReset)
    call parameters%reset()
    !$omp end critical(evolveForestReset)
    !$omp barrier
    deallocate(parameters)
    !$omp master
    if (postEvolveEvent%isAttached(self,evolveForestsCheckpoint)) call postEvolveEvent%detach(self,evolveForestsCheckpoint)
    !$omp end master
    !$omp end parallel

    ! Finalize outputs.
    call self%mergerTreeOutputter_%finalize()

    ! Destroy tree initialization lock.
    !$ call OMP_Destroy_Lock(initializationLock)
    ! Perform any post universe evolve tasks
    !![
    <include directive="universePostEvolveTask" type="functionCall" functionType="void">
    !!]
    include 'tasks.evolve_tree.universePostEvolveTask.inc'
    !![
    </include>
    !!]

    call displayUnindent('Done task: merger tree evolution')

    return
  end subroutine evolveForestsPerform

  subroutine evolveForestsSuspendTree(self,tree)
    !!{
    Suspend processing of a tree.
    !!}
#ifdef USEMPI
    use :: Error                   , only : Error_Report
#endif
    use :: ISO_Varying_String      , only : operator(//)        , varying_string
    use :: Kind_Numbers            , only : kind_int8
    use :: Merger_Tree_Construction, only : mergerTreeStateStore
    use :: String_Handling         , only : operator(//)
    implicit none
    class  (taskEvolveForests)         , intent(inout) :: self
    type   (mergerTree       ), pointer, intent(inout) :: tree
    type   (mergerTree       ), pointer                :: treeCurrent     , branchNext
    integer(kind_int8        )                         :: uniqueIDNodeBase
    type   (varying_string   )                         :: fileName

#ifdef USEMPI
    call Error_Report('suspending trees is not supported under MPI'//{introspection:location})
#endif
    ! If the tree is to be suspended to file do so now.
    if (.not.self%suspendToRAM) then
       ! Make a copy of the unique ID of the base node.
       uniqueIDNodeBase=tree%nodeBase%uniqueID()
       ! Generate a suitable file name.
       fileName=self%suspendPath//'/suspendedTree_'//uniqueIDNodeBase
       ! Store the tree to file.
       call mergerTreeStateStore(tree,char(fileName),snapshot=.false.,append=.false.)
       ! Destroy the tree(s).
       treeCurrent => tree
       do while (associated(treeCurrent))
          branchNext => treeCurrent%nextTree
          call treeCurrent%destroy()
          treeCurrent => branchNext
       end do
       ! Set the tree index to the base node unique ID so that we can resume from the correct file.
       tree%index=uniqueIDNodeBase
    end if
    call self%universeProcessed%lock%set  (    )
    call self%universeProcessed%pushTree  (tree)
    call self%universeProcessed%lock%unset(    )
    tree => null()
    return
  end subroutine evolveForestsSuspendTree

  subroutine evolveForestsResumeTree(self,tree)
    !!{
    Resume processing of a tree.
    !!}
    use :: ISO_Varying_String      , only : operator(//)           , varying_string
    use :: Merger_Tree_Construction, only : mergerTreeStateFromFile
    use :: String_Handling         , only : operator(//)
    implicit none
    class(taskEvolveForests)         , intent(inout) :: self
    type (mergerTree       ), pointer, intent(  out) :: tree
    type (varying_string   )                         :: fileName

    call self%universeWaiting%lock%set  ()
    tree => self%universeWaiting%popTree()
    call self%universeWaiting%lock%unset()
    ! If the tree was suspended to file, restore it now.
    if (.not.self%suspendToRAM.and.associated(tree)) then
       ! Generate the file name.
       fileName=self%suspendPath//'/suspendedTree_'//tree%index
       ! Read the tree from file.
       call mergerTreeStateFromFile(tree,char(fileName),self%randomNumberGenerator_,self%mergerTreeSeeds_,deleteAfterRead=.true.)
    end if
    return
  end subroutine evolveForestsResumeTree
