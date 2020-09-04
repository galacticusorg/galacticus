!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  use :: Galactic_Filters               , only : galacticFilter             , galacticFilterClass
  use :: Galacticus_Nodes               , only : mergerTree                 , treeNode                , universe
  use :: Input_Parameters               , only : inputParameters
  use :: Kind_Numbers                   , only : kind_int8
  use :: Merger_Tree_Construction       , only : mergerTreeConstructorClass
  use :: Merger_Tree_Operators          , only : mergerTreeOperatorClass
  use :: Merger_Tree_Outputters         , only : mergerTreeOutputter        , mergerTreeOutputterClass
  use :: Merger_Trees_Evolve            , only : mergerTreeEvolver          , mergerTreeEvolverClass
  use :: Output_Times                   , only : outputTimesClass
  use :: Task_Evolve_Forests_Work_Shares, only : evolveForestsWorkShareClass
  use :: Universe_Operators             , only : universeOperator           , universeOperatorClass

  !# <task name="taskEvolveForests">
  !#  <description>A task which evolves galaxies within a set of merger tree forests.</description>
  !# </task>
  type, extends(taskClass) :: taskEvolveForests
     !% Implementation of a task which evolves galaxies within a set of merger tree forests.
     private
     ! Parameter controlling maximum walltime for which forest evolution can run.
     integer         (kind_int8)                            :: walltimeMaximum
     ! Parameters controlling tree suspension.
     logical                                                :: suspendToRAM
     type            (varying_string             )          :: suspendPath
     ! Parameters controlling how threads process forests.
     logical                                                :: evolveSingleForest
     integer                                                :: evolveSingleForestSections
     double precision                                       :: evolveSingleForestMassMinimum
     ! Tree universes used while processing all trees.
     type            (universe                   ), pointer :: universeWaiting               => null(), universeProcessed       => null()
     ! Objects used in tree processing.
     class           (mergerTreeConstructorClass ), pointer :: mergerTreeConstructor_        => null()
     class           (mergerTreeOperatorClass    ), pointer :: mergerTreeOperator_           => null()
     class           (mergerTreeEvolverClass     ), pointer :: mergerTreeEvolver_            => null()
     class           (mergerTreeOutputterClass   ), pointer :: mergerTreeOutputter_          => null()
     class           (galacticFilterClass        ), pointer :: galacticFilter_               => null()
     class           (evolveForestsWorkShareClass), pointer :: evolveForestsWorkShare_       => null()
     class           (outputTimesClass           ), pointer :: outputTimes_                  => null()
     class           (universeOperatorClass      ), pointer :: universeOperator_             => null()
     ! Pointer to the parameters for this task.
     type            (inputParameters            )          :: parameters
   contains
     !@ <objectMethods>
     !@   <object>taskEvolveForests</object>
     !@   <objectMethod>
     !@     <method>suspendTree</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(mergerTree)\textgreater} tree\argout</arguments>
     !@     <description>Suspend a tree (to memory or to file).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>resumeTree</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(mergerTree)\textgreater} tree\argout</arguments>
     !@     <description>Restore a suspended tree.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                evolveForestsDestructor
     procedure :: perform     => evolveForestsPerform
     procedure :: suspendTree => evolveForestsSuspendTree
     procedure :: resumeTree  => evolveForestsResumeTree
     procedure :: autoHook    => evolveForestsAutoHook
  end type taskEvolveForests

  interface taskEvolveForests
     !% Constructors for the {\normalfont \ttfamily evolveForests} task.
     module procedure evolveForestsConstructorParameters
     module procedure evolveForestsConstructorInternal
  end interface taskEvolveForests

  ! Class used to build lists of tree branches which get processed independently.
  type :: evolveForestsBranchList
     type(treeNode               ), pointer :: nodeParent
     type(mergerTree             ), pointer :: branch
     type(evolveForestsBranchList), pointer :: next
  end type evolveForestsBranchList

  ! Copies of objects used by each thread.
  class(mergerTreeOutputterClass  ), pointer :: evolveForestsMergerTreeOutputter_   => null()
  class(mergerTreeEvolverClass    ), pointer :: evolveForestsMergerTreeEvolver_     => null()
  class(mergerTreeConstructorClass), pointer :: evolveForestsMergerTreeConstructor_ => null()
  class(mergerTreeOperatorClass   ), pointer :: evolveForestsMergerTreeOperator_    => null()
  class(galacticFilterClass       ), pointer :: evolveForestsGalacticFilter_        => null()
  !$omp threadprivate(evolveForestsMergerTreeOutputter_,evolveForestsMergerTreeEvolver_,evolveForestsMergerTreeConstructor_,evolveForestsMergerTreeOperator_,evolveForestsGalacticFilter_)

contains

  function evolveForestsConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily evolveForests} task class which takes a parameter set as input.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskEvolveForests          )                        :: self
    type            (inputParameters            ), intent(inout), target :: parameters
    class           (mergerTreeOperatorClass    ), pointer               :: mergerTreeOperator_
    class           (evolveForestsWorkShareClass), pointer               :: evolveForestsWorkShare_
    class           (mergerTreeConstructorClass ), pointer               :: mergerTreeConstructor_
    class           (outputTimesClass           ), pointer               :: outputTimes_
    class           (universeOperatorClass      ), pointer               :: universeOperator_
    class           (mergerTreeEvolverClass     ), pointer               :: mergerTreeEvolver_
    class           (mergerTreeOutputterClass   ), pointer               :: mergerTreeOutputter_
    class           (galacticFilterClass        ), pointer               :: galacticFilter_
    type            (inputParameters            ), pointer               :: parametersRoot
    logical                                                              :: evolveSingleForest           , suspendToRAM
    integer                                                              :: evolveSingleForestSections
    double precision                                                     :: evolveSingleForestMassMinimum
    integer         (kind_int8                  )                        :: walltimeMaximum
    type            (varying_string             )                        :: suspendPath

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
    !# <inputParameter>
    !#   <name>walltimeMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-1_kind_int8</defaultValue>
    !#   <description>If set to a positive number, this is the maximum wall time for which forest evolution is allowed to proceed before the task gives up.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>evolveSingleForest</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true then each forest is processed sequentially, with multiple parallel threads (if available) working on the same forest. If false, multiple forests are processed simultaneously, with a single parallel thread (if available) working on each.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>evolveSingleForestSections</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>100</defaultValue>
    !#   <description>The number of timesteps into which forests should be split when processing single forests in parallel.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>evolveSingleForestMassMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The minimum tree mass for which forests should be processed in parallel.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>suspendToRAM</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Specifies whether trees should be suspended to RAM (otherwise they are suspend to file).</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    if (.not.suspendToRAM) then
       !# <inputParameter>
       !#   <name>suspendPath</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The path to which tree suspension files will be stored.</description>
       !#   <source>parameters</source>
       !#   <type>string</type>
       !# </inputParameter>
    end if
    !# <objectBuilder class="mergerTreeConstructor"  name="mergerTreeConstructor_"  source="parameters"/>
    !# <objectBuilder class="mergerTreeOperator"     name="mergerTreeOperator_"     source="parameters"/>
    !# <objectBuilder class="evolveForestsWorkShare" name="evolveForestsWorkShare_" source="parameters"/>
    !# <objectBuilder class="outputTimes"            name="outputTimes_"            source="parameters"/>
    !# <objectBuilder class="universeOperator"       name="universeOperator_"       source="parameters"/>
    !# <objectBuilder class="mergerTreeEvolver"      name="mergerTreeEvolver_"      source="parameters"/>
    !# <objectBuilder class="mergerTreeOutputter"    name="mergerTreeOutputter_"    source="parameters"/>
    !# <objectBuilder class="galacticFilter"         name="galacticFilter_"         source="parameters"/>
    if (associated(parametersRoot)) then
       self=taskEvolveForests(evolveSingleForest,evolveSingleForestSections,evolveSingleForestMassMinimum,walltimeMaximum,suspendToRAM,suspendPath,mergerTreeConstructor_,mergerTreeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,galacticFilter_,parametersRoot)
    else
       self=taskEvolveForests(evolveSingleForest,evolveSingleForestSections,evolveSingleForestMassMinimum,walltimeMaximum,suspendToRAM,suspendPath,mergerTreeConstructor_,mergerTreeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,galacticFilter_,parameters    )
    end if
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="mergerTreeConstructor_" />
    !# <objectDestructor name="mergerTreeOperator_"    />
    !# <objectDestructor name="evolveForestsWorkShare_"/>
    !# <objectDestructor name="outputTimes_"           />
    !# <objectDestructor name="universeOperator_"      />
    !# <objectDestructor name="mergerTreeEvolver_"     />
    !# <objectDestructor name="mergerTreeOutputter_"   />
    !# <objectDestructor name="galacticFilter_"        />
    return
  end function evolveForestsConstructorParameters

  function evolveForestsConstructorInternal(evolveSingleForest,evolveSingleForestSections,evolveSingleForestMassMinimum,walltimeMaximum,suspendToRAM,suspendPath,mergerTreeConstructor_,mergerTreeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,galacticFilter_,parameters) result(self)
    !% Internal constructor for the {\normalfont \ttfamily evolveForests} task class.
    implicit none
    type            (taskEvolveForests          )                        :: self
    logical                                      , intent(in   )         :: evolveSingleForest           , suspendToRAM
    integer                                      , intent(in   )         :: evolveSingleForestSections
    double precision                             , intent(in   )         :: evolveSingleForestMassMinimum
    integer         (kind_int8                  ), intent(in   )         :: walltimeMaximum
    type            (varying_string             ), intent(in   )         :: suspendPath
    class           (mergerTreeConstructorClass ), intent(in   ), target :: mergerTreeConstructor_
    class           (mergerTreeOperatorClass    ), intent(in   ), target :: mergerTreeOperator_
    class           (evolveForestsWorkShareClass), intent(in   ), target :: evolveForestsWorkShare_
    class           (outputTimesClass           ), intent(in   ), target :: outputTimes_
    class           (universeOperatorClass      ), intent(in   ), target :: universeOperator_
    class           (mergerTreeEvolverClass     ), intent(in   ), target :: mergerTreeEvolver_
    class           (mergerTreeOutputterClass   ), intent(in   ), target :: mergerTreeOutputter_
    class           (galacticFilterClass        ), intent(in   ), target :: galacticFilter_
    type            (inputParameters            ), intent(in   ), target :: parameters
    !# <constructorAssign variables="evolveSingleForest, evolveSingleForestSections, evolveSingleForestMassMinimum, walltimeMaximum, suspendToRAM, suspendPath, *mergerTreeConstructor_, *mergerTreeOperator_, *evolveForestsWorkShare_, *outputTimes_, *universeOperator_, *mergerTreeEvolver_, *mergerTreeOutputter_, *galacticFilter_"/>

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    return
  end function evolveForestsConstructorInternal

  subroutine evolveForestsAutoHook(self)
    !% Attach to state store/restore event hooks.
    use :: Events_Hooks, only : openMPThreadBindingNone, stateRestoreEvent, stateStoreEvent
    implicit none
    class(taskEvolveForests), intent(inout) :: self

    call stateStoreEvent  %attach(self,evolveForestsStateStore  ,openMPThreadBindingNone)
    call stateRestoreEvent%attach(self,evolveForestsStateRestore,openMPThreadBindingNone)
    return
  end subroutine evolveForestsAutoHook

  subroutine evolveForestsStateStore(self,stateFile,gslStateFile,stateOperationID)
    !% Store the internal state of this object.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_ptr
    implicit none
    class  (taskEvolveForests), intent(inout) :: self
    integer                   , intent(in   ) :: stateFile
    type   (c_ptr            ), intent(in   ) :: gslStateFile
    integer(c_size_t         ), intent(in   ) :: stateOperationID
    !$GLC attributes unused :: self

    call evolveForestsMergerTreeConstructor_%stateStore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeOperator_   %stateStore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeEvolver_    %stateStore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeOutputter_  %stateStore(stateFile,gslStateFile,stateOperationID)
    return
  end subroutine evolveForestsStateStore

  subroutine evolveForestsStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !% Store the internal state of this object.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_ptr
    implicit none
    class  (taskEvolveForests), intent(inout) :: self
    integer                   , intent(in   ) :: stateFile
    type   (c_ptr            ), intent(in   ) :: gslStateFile
    integer(c_size_t         ), intent(in   ) :: stateOperationID
    !$GLC attributes unused :: self

    call evolveForestsMergerTreeConstructor_%stateRestore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeOperator_   %stateRestore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeEvolver_    %stateRestore(stateFile,gslStateFile,stateOperationID)
    call evolveForestsMergerTreeOutputter_  %stateRestore(stateFile,gslStateFile,stateOperationID)
    return
  end subroutine evolveForestsStateRestore

  subroutine evolveForestsDestructor(self)
    !% Destructor for the {\normalfont \ttfamily evolveForests} task class.
    use :: Events_Hooks   , only : stateRestoreEvent           , stateStoreEvent
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskEvolveForests), intent(inout) :: self

    !# <objectDestructor name="self%mergerTreeConstructor_" />
    !# <objectDestructor name="self%mergerTreeOperator_"    />
    !# <objectDestructor name="self%evolveForestsWorkShare_"/>
    !# <objectDestructor name="self%outputTimes_"           />
    !# <objectDestructor name="self%universeOperator_"      />
    !# <objectDestructor name="self%mergerTreeEvolver_"     />
    !# <objectDestructor name="self%mergerTreeOutputter_"   />
    !# <objectDestructor name="self%galacticFilter_"        />
    if (associated(self%universeWaiting  )) deallocate(self%universeWaiting  )
    if (associated(self%universeProcessed)) deallocate(self%universeProcessed)
    call stateStoreEvent  %detach(self,evolveForestsStateStore  )
    call stateRestoreEvent%detach(self,evolveForestsStateRestore)
    call Node_Components_Uninitialize()
    return
  end subroutine evolveForestsDestructor

  subroutine evolveForestsPerform(self,status)
    !% Evolves the complete set of merger trees as specified.
    use               :: Galacticus_Display                  , only : Galacticus_Display_Indent          , Galacticus_Display_Message         , Galacticus_Display_Unindent, verbosityInfo
    use               :: Galacticus_Error                    , only : Galacticus_Error_Report            , errorStatusSuccess
    use               :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
    use               :: Galacticus_Nodes                    , only : mergerTree                         , nodeComponentBasic                 , treeNode                   , universe        , &
          &                                                           universeEvent
    use   , intrinsic :: ISO_C_Binding                       , only : c_size_t
    
    use               :: Memory_Management                   , only : Memory_Usage_Record                , memoryTypeNodes
    use               :: Merger_Tree_Walkers                 , only : mergerTreeWalkerAllNodes           , mergerTreeWalkerIsolatedNodes
    use               :: Merger_Trees_Initialize             , only : Merger_Tree_Initialize
    use               :: Node_Components                     , only : Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize
    use               :: Node_Events_Inter_Tree              , only : Inter_Tree_Event_Post_Evolve
    !$ use            :: OMP_Lib                             , only : omp_lock_kind                      , OMP_Init_Lock                      , OMP_Get_Thread_Num         , OMP_Destroy_Lock
    use               :: Sorting                             , only : sortIndex
    use               :: String_Handling                     , only : operator(//)
    ! Include modules needed for tasks.
    !# <include directive="universePostEvolveTask" type="moduleUse" functionType="void">
    include 'galacticus.tasks.evolve_tree.universePostEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="mergerTreeExtraOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.extra.modules.inc'
    !# </include>
    implicit none
    class           (taskEvolveForests            ), intent(inout), target           :: self
    integer                                        , intent(  out), optional         :: status
    type            (mergerTree                   ), pointer                  , save :: tree
    logical                                                                   , save :: finished                                  , treeIsNew
    integer         (c_size_t                     )                           , save :: iOutput
    double precision                                                          , save :: evolveToTime                              , treeTimeEarliest            , &
         &                                                                              universalEvolveToTime                     , treeTimeLatest              , &
         &                                                                              outputTimeNext
    type            (varying_string               )                           , save :: message
    character       (len=20                       )                           , save :: label
    !$omp threadprivate(tree,finished,iOutput,evolveToTime,message,label,treeIsNew,treeTimeEarliest,outputTimeNext)
    logical                                                                   , save :: treeIsFinished                            , evolutionIsEventLimited     , &
         &                                                                              success                                   , removeTree                  , &
         &                                                                              suspendTree                               , treesDidEvolve              , &
         &                                                                              treeDidEvolve                             , treesCouldEvolve            , &
         &                                                                              deadlockReport
    type            (mergerTree                   ), pointer                  , save :: currentTree                               , previousTree                , &
         &                                                                              nextTree
    type            (mergerTreeWalkerIsolatedNodes)                                  :: treeWalkerIsolated
    type            (mergerTreeWalkerAllNodes     )                           , save :: treeWalkerAll
    !$omp threadprivate(currentTree,previousTree,treeWalkerAll)
    type            (treeNode                     ), pointer                  , save :: satelliteNode
    class           (nodeComponentBasic           ), pointer                  , save :: baseNodeBasic
    !$omp threadprivate(satelliteNode,baseNodeBasic,treeIsFinished,evolutionIsEventLimited,success,removeTree,suspendTree,treeDidEvolve)
    type            (universeEvent                ), pointer                  , save :: event_
    !$omp threadprivate(event_)
    ! Variables used in processing individual forests in parallel.
    double precision                                                          , save :: timeBranchSplit
    type            (treeNode                     ), pointer                  , save :: node
    class           (nodeComponentBasic           ), pointer                  , save :: basic                                     , basicChild
    type            (evolveForestsBranchList      ), pointer                  , save :: branchList_                               , branchNew                   , &
         &                                                                              branchNext
    logical                                                                   , save :: branchAccept
    integer         (c_size_t                     )                           , save :: iBranch                                   , i                           , &
         &                                                                              countBranch                               , treeNumber
    integer         (c_size_t                     ), allocatable, dimension(:), save :: rankBranch
    double precision                               , allocatable, dimension(:), save :: massBranch
    integer                                                                   , save :: forestSection
    double precision                                                          , save :: timeSectionForestBegin
    logical                                                                          :: triggerExit                               , triggerFinishUniverse       , &
         &                                                                              disableSingleForestEvolution              , triggerFinishFinal          , &
         &                                                                              triggerFinish
    integer         (c_size_t                     )                                  :: iBranchAcceptedLast                       , treeCount
    integer         (omp_lock_kind                )                                  :: initializationLock
    integer         (kind_int8                    )                                  :: systemClockRate                           , systemClockMaximum
    double precision                                                                 :: evolveToTimeForest
    !$omp threadprivate(node,basic,basicChild,timeBranchSplit,branchNew,branchNext,i,iBranch,branchAccept,massBranch,timeSectionForestBegin,forestSection,treeNumber)

    ! The following processes merger trees, one at a time, to each successive output time, then dumps their contents to file. It
    ! allows for the possibility of "universal events" - events which require all merger trees to reach the same cosmic time. If
    ! such an event exists, each tree is processed up to that time and then pushed onto a stack where it waits to be
    ! processed. Once all trees reach the event time the stack of trees is passed to the event task. Tree processing then
    ! continues by popping trees off of the stack and processing them further (possibly to the next universal event).

    call Galacticus_Display_Indent('Begin task: merger tree evolution')

    ! Set status to success by default.
    if (present(status)) status=errorStatusSuccess

    ! Initialize a lock used for controling tree initialization.
    !$ call OMP_Init_Lock(initializationLock)

    ! Initialize tree counter and record that we are not finished processing trees.
    deadlockReport              =.false.
    finished                    =.false.
    triggerFinish               =.false.
    triggerFinishFinal          =.false.
    triggerFinishUniverse       =.false.
    disableSingleForestEvolution=.false.
    ! Initialize universes which will act as tree stacks. We use two stacks: one for trees waiting to be processed, one for trees
    ! that have already been processed.
    allocate(self%universeWaiting  )
    allocate(self%universeProcessed)
    self%universeWaiting  %trees => null()
    self%universeProcessed%trees => null()
    call self%universeWaiting%attributes%initialize()

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
    !$omp parallel copyin(finished)
    allocate(evolveForestsMergerTreeOutputter_  ,mold=self%mergerTreeOutputter_  )
    allocate(evolveForestsMergerTreeEvolver_    ,mold=self%mergerTreeEvolver_    )
    allocate(evolveForestsMergerTreeConstructor_,mold=self%mergerTreeConstructor_)
    allocate(evolveForestsMergerTreeOperator_   ,mold=self%mergerTreeOperator_   )
    allocate(evolveForestsGalacticFilter_       ,mold=self%galacticFilter_       )
    !$omp critical(evolveForestsDeepCopy)
    !# <deepCopyReset variables="self%mergerTreeEvolver_ self%mergerTreeOutputter_ self%mergerTreeConstructor_ self%mergerTreeOperator_ self%galacticFilter_"/>
    !# <deepCopy source="self%mergerTreeEvolver_"     destination="evolveForestsMergerTreeEvolver_"    />
    !# <deepCopy source="self%mergerTreeOutputter_"   destination="evolveForestsMergerTreeOutputter_"  />
    !# <deepCopy source="self%mergerTreeConstructor_" destination="evolveForestsMergerTreeConstructor_"/>
    !# <deepCopy source="self%mergerTreeOperator_"    destination="evolveForestsMergerTreeOperator_"   />
    !# <deepCopy source="self%galacticFilter_"        destination="evolveForestsGalacticFilter_"       />
    !$omp end critical(evolveForestsDeepCopy)
    ! Call routines to perform initializations which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Allow events to be attached to the universe.
    !$omp master
    self%universeWaiting%event => null()
    !# <eventHook name="universePreEvolve">
    !#  <callWith>self%universeWaiting</callWith>
    !# </eventHook>
    call self%universeOperator_%operate(self%universeWaiting)
    !$omp end master
    !$omp barrier
    ! Begin processing trees.
    treeProcess : do while (.not.finished)
       ! For single forest evolution, only the master thread should retrieve a merger tree.
       singleForestTreeFetch : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
          ! Attempt to get a new tree to process. We first try to get a new tree. If no new trees exist, we will look for a tree on
          ! the stack waiting to be processed.
          ! Perform any pre-tree construction tasks.
          call evolveForestsMergerTreeOperator_%operatePreConstruction()
          ! Get the number of the next tree to process.
          treeNumber=self%evolveForestsWorkShare_%forestNumber(utilizeOpenMPThreads=.not.self%evolveSingleForest)
          ! Get a tree.
          tree                                    => evolveForestsMergerTreeConstructor_%construct(treeNumber)
          if (associated(tree)) tree%hostUniverse => self%universeWaiting
          finished                                =  finished.or..not.associated(tree)
          treeIsNew                               =  .not.finished
          ! If no new tree was available, attempt to pop one off the universe stack.
          if (finished) then
             call self%resumeTree(tree)
             treeIsNew=.false.
             finished =.not.associated(tree)
          end if
          if (self%evolveSingleForest .and. finished) triggerFinish=.true.
       end if singleForestTreeFetch
       ! For single forest evolution, block threads until master has retrieved a merger tree.
       if (self%evolveSingleForest) then
          !$omp barrier
          finished=triggerFinish
          !$omp barrier
       end if
       ! If we got a tree (i.e. we are not "finished") process it.
       if (.not.finished) then
          treeIsFinished=.false.
          ! For single forest evolution, only the master thread should determine evolution time.
          singleForestEvolveTime : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
             ! If this is a new tree, perform any pre-evolution tasks on it.
             if (treeIsNew) then
                call evolveForestsMergerTreeOperator_%operatePreEvolution(tree)
                message="Evolving tree number "
             else
                message="Resuming tree number "
             end if
             ! Display a message.
             message=message//tree%index//" {"//tree%baseNode%index()//"}"
             call Galacticus_Display_Indent(message)
             ! Determine if this will be the final forest evolved by multiple threads.
             if (self%evolveSingleForest) then
                basic =>tree%baseNode%basic()
                if (basic%mass() < self%evolveSingleForestMassMinimum) disableSingleForestEvolution=.true.
             end if
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
                iOutput      =self%outputTimes_%count()+1
                treeIsFinished=.true.
             end if
          end if singleForestEvolveTime
          ! For single forest evolution, block threads until master thread has determined evolution time.
          if (self%evolveSingleForest) then
             if (OMP_Get_Thread_Num() == 0) triggerExit=.false.
             !$omp barrier
          end if
          ! Iterate evolving the tree until no more outputs are required.
          treeEvolveLoop : do while (iOutput <= self%outputTimes_%count())
             ! Ping the work-share object.
             call self%evolveForestsWorkShare_%ping()
             ! For single forest evolution, maximum evolution time is determined by the master thread only.
             singleForestMaximumTime : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
                ! We want to find the maximum time to which we can evolve this tree. This will be the minimum of the next output
                ! time (at which we must stop and output the tree) and the next universal event time (at which we must stop and
                ! perform the event task). Find the next output time.
                evolveToTime=self%outputTimes_%time(iOutput)
                ! Find the earliest universe event.
                !$omp critical(universeTransform)
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
                !$omp end critical(universeTransform)
                if (tree%earliestTime() <= evolveToTime) treesCouldEvolve=.true.
             end if singleForestMaximumTime
             ! For single forest evolution, block all threads until master thread has determined the maximum evolution time.
             if (self%evolveSingleForest) then
                !$omp barrier
             end if
             ! If single trees are to be broken between multiple threads and processed in parallel, then do that here. Begin by
             ! finding evolvable branches at some earlier time.
             if (self%evolveSingleForest) then
                timeSectionForestBegin=tree%earliestTimeEvolving()
                do forestSection=1,self%evolveSingleForestSections
                   ! Master thread alone performs the splitting of the tree into branches.
                   singleForestTreeSplit : if (OMP_Get_Thread_Num() == 0) then
                      timeBranchSplit    =  timeSectionForestBegin+(evolveToTime-timeSectionForestBegin)*dble(forestSection)/dble(self%evolveSingleForestSections)
                      countBranch        =  0
                      evolveToTimeForest =  evolveToTime
                      treeWalkerIsolated=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
                      do while (treeWalkerIsolated%next(node))
                         if (associated(node%firstChild).and.associated(node%parent)) then
                            basic      => node           %basic()
                            basicChild => node%firstChild%basic()
                            if (basic%time() >= timeBranchSplit .and. basicChild%time() < timeBranchSplit) then
                               countBranch=countBranch+1
                               allocate(branchNew      )
                               allocate(branchNew%branch)
                               ! Build the new branch. Ensure that our node, which forms the base node of the new branch, has
                               ! its hostTree pointer set to point to the new branch so that merger tree walkers can correctly
                               ! navigate the branch without straying into a neighboring branch.
                               branchNew%branch            =  node%hostTree
                               !# <referenceCountIncrement owner="branchNew%branch" object="randomNumberGenerator_"/>
                               branchNew%branch%baseNode   => node
                               branchNew       %nodeParent => node%parent
                               node            %hostTree   => branchNew%branch
                               if (associated(branchList_)) then
                                  branchNew%next => branchList_
                               else
                                  branchNew%next => null()
                               end if
                               branchList_ => branchNew
                            end if
                         end if
                      end do
                      ! If no more branches were found to process, exit.
                      if (countBranch > 0) then
                         ! Rank trees by mass.
                         allocate(massBranch(countBranch))
                         allocate(rankBranch(countBranch))
                         branchNew    => branchList_
                         countBranch =  0
                         do while (associated(branchNew))
                            countBranch             =  countBranch                        +1
                            basic                   => branchNew  %branch%baseNode%basic()
                            massBranch(countBranch) =  basic                      %mass ()
                            branchNew               => branchNew                  %next
                         end do
                         rankBranch=sortIndex(massBranch)
                         ! Process trees.
                         currentTree => tree
                         if (associated(currentTree%baseNode)) then
                            do while (associated(currentTree))
                               basic => currentTree%baseNode%basic()
                               call Merger_Tree_Initialize(currentTree,basic%time())
                               currentTree => currentTree%nextTree
                            end do
                         end if
                         iBranchAcceptedLast=0
                      end if
                   end if singleForestTreeSplit
                   ! Block threads until master thread has completed splitting forest into branches.
                   !$omp barrier
                   ! Threads now process branches of the forest.
                   iBranch=0
                   do while (iBranch < countBranch)
                      iBranch     =iBranch+1
                      branchAccept=.false.
                      !$omp critical (branchList_Share)
                      if (iBranch > iBranchAcceptedLast) then
                         branchAccept       =.true.
                         iBranchAcceptedLast=iBranch
                      end if
                      !$omp end critical (branchList_Share)
                      if (branchAccept) then
                         branchNew => branchList_
                         do i=1,rankBranch(countBranch+1-iBranch)-1
                            branchNew => branchNew%next
                         end do
                         branchNew%branch%baseNode%parent =>                           null ()
                         basic                            => branchNew%branch%baseNode%basic()
                         call evolveForestsMergerTreeEvolver_%evolve(                                        &
                              &                                          branchNew    %branch              , &
                              &                                      min(                                    &
                              &                                          basic        %time              (), &
                              &                                                        evolveToTimeForest    &
                              &                                         )                                  , &
                              &                                          treeDidEvolve                     , &
                              &                                          suspendTree                       , &
                              &                                          deadlockReport                    , &
                              &                                          systemClockMaximum                  &
                           !$ &                                         ,initializationLock                  &
                              &                        )
                         !$omp critical (universeStatus)
                         if (treeDidEvolve) treesDidEvolve         =.true.
                         if (suspendTree  ) evolutionIsEventLimited=.true.
                         !$omp end critical (universeStatus)
                      end if
                   end do
                   ! Block all threads until all branches are finished processing.
                   !$omp barrier
                   ! Only the master thread performs clean-up of the branches, and re-linking into the forest.
                   singleForestTreeMerge : if (OMP_Get_Thread_Num() == 0 .and. countBranch > 0) then
                      ! Clean up.
                      branchNew => branchList_
                      do while (associated(branchNew))
                         branchNew%branch%baseNode%hostTree => branchNew%nodeParent%hostTree
                         branchNew%branch%baseNode%parent   => branchNew%nodeParent
                         branchNew%branch%baseNode          => null()
                         call branchNew%branch%destroy()
                         branchNext => branchNew%next
                         deallocate(branchNew)
                         branchNew => branchNext
                      end do
                      nullify   (branchList_)
                      deallocate(massBranch )
                      deallocate(rankBranch )
                   end if singleForestTreeMerge
                   ! Block all threads until branches have been re-merged into the forest by the master thread.
                   !$omp barrier
                end do
             end if
             ! For single forest evolution, only the master thread should finalize evolution of the merger tree.
             singleForestFinalizeEvolution : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
                ! Evolve the tree to the computed time.
                call evolveForestsMergerTreeEvolver_%evolve(tree,evolveToTime,treeDidEvolve,suspendTree,deadlockReport,systemClockMaximum,status=status)
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
                   if (associated(currentTree%baseNode)) then
                      baseNodeBasic => currentTree%baseNode%basic()
                      removeTree    =   .not.associated(currentTree%baseNode%firstChild) &
                           &           .and.                                             &
                           &            (baseNodeBasic%time() < evolveToTime)
                      if (removeTree) then
                         ! Does the node have attached satellites which are about to merge?
                         satelliteNode => currentTree%baseNode%firstSatellite
                         do while (associated(satelliteNode))
                            if (associated(satelliteNode%mergeTarget)) then
                               removeTree=.false.
                               exit
                            end if
                            satelliteNode => satelliteNode%sibling
                         end do
                         ! Does the node have attached events?
                         if (associated(currentTree%baseNode%event)) removeTree=.false.
                      end if
                   else
                      ! No need to remove already empty trees.
                      removeTree=.false.
                   end if
                   if (removeTree) then
                      message="Removing remnant tree "
                      message=message//currentTree%index//" {"//currentTree%baseNode%index()//"}"
                      call Galacticus_Display_Message(message,verbosityInfo)
                      if (.not.associated(previousTree)) then
                         nextTree    => currentTree%nextTree
                         call currentTree%destroy()
                         currentTree => nextTree
                      else
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
                   call Galacticus_Error_Report(                                           &
                        &                       'failed to evolve tree to required time'// &
                        &                       {introspection:location}                   &
                        &                      )
                end if
                ! Determine what limited evolution.
                if (evolutionIsEventLimited) then
                   ! Tree evolution was limited by a universal event. Therefore it can evolve no further
                   ! until that event's task is performed.
                   exit
                else
                   ! Tree reached an output time, so output it. We can then continue evolving.
                   write (label,'(f7.2)') evolveToTime
                   message="Output tree data at t="//trim(label)//" Gyr"
                   call Galacticus_Display_Message(message)
                   call evolveForestsMergerTreeOutputter_%output(tree,iOutput,evolveToTime,.false.)
                   ! Perform any extra output and post-output processing on nodes.
                   treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
                   do while (treeWalkerAll%next(node))
                      basic => node%basic()
                      if (basic%time() == evolveToTime) then
                         !# <include directive="mergerTreeExtraOutputTask" type="functionCall" functionType="void">
                         !#  <functionArgs>node,iOutput,node%hostTree%index,evolveForestsGalacticFilter_%passes(node)</functionArgs>
                         include 'galacticus.output.merger_tree.tasks.extra.inc'
                         !# </include>
                         !# <eventHook name="mergerTreeExtraOutput">
                         !#  <callWith>node,iOutput,node%hostTree,evolveForestsGalacticFilter_%passes(node)</callWith>
                         !# </eventHook>
                         call node%postOutput(evolveToTime)
                      end if
                   end do
                   iOutput=iOutput+1
                   ! If all output times have been reached, we're finished.
                   if (iOutput > self%outputTimes_%count()) then
                      treeIsFinished=.true.
                      ! For single forest evolution, record that we are exiting the evolution loop. This will be used to inform
                      ! all non-master threads to exit also.  Note that this barrier corresponds with the one labelled
                      ! "singleForestExitEvolutionOther" below (which will be seen by all non-master threads).
                      singleForestExitEvolutionMaster : if (self%evolveSingleForest) then
                         triggerExit=.true.
                         !$omp barrier
                      end if singleForestExitEvolutionMaster
                      exit
                   end if
                end if
             end if singleForestFinalizeEvolution
             ! For single forest evolution, block until the master thread has finished finalizing evolution and, if no more
             ! evolution is required, exit the evolution loop. Note that this barrier corresponds with the one labelled
             ! "evolveSingleForest" above (which will be seen by the master thread alone).
             singleForestExitEvolutionOther : if (self%evolveSingleForest) then
                !$omp barrier
                if (triggerExit) exit
             end if singleForestExitEvolutionOther
          end do treeEvolveLoop
          ! For single forest evolution, block threads until all have completed evolution of the forest.
          if (self%evolveSingleForest) then
             !$omp barrier
          end if
          ! For single forest evolution, only the master thread handles suspension of the tree.
          singleForstSuspend : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
             ! If tree could not evolve further, but is not finished, push it to the universe stack.
             if (.not.treeIsFinished) then
                ! Suspend the tree.
                call self%suspendTree(tree)
                ! Unindent messages.
                call Galacticus_Display_Unindent('Suspending tree')
             else
                ! Unindent messages.
                call Galacticus_Display_Unindent('Finished tree'  )
             end if
          end if singleForstSuspend
          ! For single forest evolution, block all threads until the master thread has completed tree suspension.
          if (self%evolveSingleForest) then
             !$omp barrier
          end if
       end if
       ! For single forest evolution, only the master thread performs tree destruction.
       singleForestTreeDestroy : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
          ! Destroy the tree.
          if (associated(tree)) then
             currentTree => tree
             do while (associated(currentTree))
                previousTree => currentTree
                currentTree  => currentTree%nextTree
                call previousTree%destroy()
                ! Deallocate the tree.
                call Memory_Usage_Record(sizeof(previousTree),addRemove=-1,memoryType=memoryTypeNodes)
                deallocate(previousTree)
             end do
             nullify(tree)
          end if
          ! Perform any post-evolution operations on the tree.
          if (treeIsFinished) call evolveForestsMergerTreeOperator_%operatePostEvolution()
       end if singleForestTreeDestroy
       ! For single forest evolution, block all threads until tree destruction is completed by the master thread.
       if (self%evolveSingleForest) then
          !$omp barrier
       end if
       ! For single forest evolution, only the master thread handles universe events.
       singleForestUniverseEvents : if (OMP_Get_Thread_Num() == 0 .or. .not.self%evolveSingleForest) then
          ! If any trees were pushed onto the processed stack, then there must be an event to process.
          if (finished) then
             if (.not.self%evolveSingleForest) then
                !$omp barrier
             end if
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
                   !$omp critical(universeTransform)
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
                   !$omp end critical(universeTransform)
                   call Galacticus_Display_Message(message)
                   call Inter_Tree_Event_Post_Evolve()
                   call Galacticus_Error_Report('exiting'//{introspection:location})
                else
                   deadlockReport=.true.
                end if
             end if
             treesDidEvolve=.false.
             !$omp critical(universeTransform)
             if (associated(self%universeProcessed%trees)) then
                ! Transfer processed trees back to the waiting universe.
                self%universeWaiting  %trees => self%universeProcessed%trees
                self%universeProcessed%trees => null()
                ! Find the event to process.
                event_ => self%universeWaiting%event
                do while (associated(event_))
                   if (event_%time < universalEvolveToTime) then
                      call Galacticus_Error_Report('a universal event exists in the past - this should not happen'//{introspection:location})
                   else if (event_%time == universalEvolveToTime) then
                      success=event_%task(self%universeWaiting)
                      if (success) call self%universeWaiting%removeEvent(event_)
                      exit
                   end if
                   event_ => event_%next
                end do
                ! Mark that there is more work to do.
                triggerFinishUniverse=.true.
                call Galacticus_Display_Message('Finished universe evolution pass')
             end if
             !$omp end critical(universeTransform)
             !$omp end master
             if (.not.self%evolveSingleForest) then
                !$omp barrier
             end if
             if (triggerFinishUniverse) finished=.false.
             if (.not.self%evolveSingleForest) then
                !$omp barrier
                !$omp master
                triggerFinishUniverse=.false.
                !$omp end master
             end if
          end if
       end if singleForestUniverseEvents
       ! For single forest evolution, block all threads until universe events are processed by the master thread.
       if (self%evolveSingleForest) then
          !$omp barrier
          if (OMP_Get_Thread_Num() == 0 .and. finished) triggerFinishFinal=.true.
          !$omp barrier
          finished=triggerFinishFinal
       end if
       ! Decide whether to switch off single forest evolution.
       if (self%evolveSingleForest) then
          !$omp barrier
          if (OMP_Get_Thread_Num() == 0 .and. disableSingleForestEvolution) &
               & self%evolveSingleForest=.false.
          !$omp barrier
       end if
    end do treeProcess
    ! Finalize any merger tree operator.
    call evolveForestsMergerTreeOperator_ %finalize(                         )
    ! Reduce outputs back into the original outputter object.
    call evolveForestsMergerTreeOutputter_%reduce  (self%mergerTreeOutputter_)
    ! Explicitly deallocate objects.
    call Node_Components_Thread_Uninitialize()
    call Galacticus_Function_Classes_Destroy()
    !# <objectDestructor name="evolveForestsMergerTreeOutputter_"  />
    !# <objectDestructor name="evolveForestsMergerTreeEvolver_"    />
    !# <objectDestructor name="evolveForestsMergerTreeConstructor_"/>
    !# <objectDestructor name="evolveForestsMergerTreeOperator_"   />
    !# <objectDestructor name="evolveForestsGalacticFilter_"       />
    !$omp end parallel

    ! Finalize outputs.
    call self%mergerTreeOutputter_%finalize()

    ! Destroy tree initialization lock.
    !$ call OMP_Destroy_Lock(initializationLock)
    ! Perform any post universe evolve tasks
    !# <include directive="universePostEvolveTask" type="functionCall" functionType="void">
    include 'galacticus.tasks.evolve_tree.universePostEvolveTask.inc'
    !# </include>

    call Galacticus_Display_Unindent('Done task: merger tree evolution')

    return
  end subroutine evolveForestsPerform

  subroutine evolveForestsSuspendTree(self,tree)
    !% Suspend processing of a tree.
#ifdef USEMPI
    use :: Galacticus_Error        , only : Galacticus_Error_Report
#endif
    use :: ISO_Varying_String      , only : varying_string         , operator(//)
    use :: Kind_Numbers            , only : kind_int8
    use :: Merger_Tree_Construction, only : mergerTreeStateStore
    use :: String_Handling         , only : operator(//)
    implicit none
    class  (taskEvolveForests)         , intent(inout) :: self
    type   (mergerTree       ), pointer, intent(inout) :: tree
    type   (mergerTree       ), pointer                :: treeCurrent     , branchNext
    integer(kind_int8        )                         :: baseNodeUniqueID
    type   (varying_string   )                         :: fileName

#ifdef USEMPI
    call Galacticus_Error_Report('suspending trees is not supported under MPI'//{introspection:location})
#endif
    ! If the tree is to be suspended to file do so now.
    if (.not.self%suspendToRAM) then
       ! Make a copy of the unique ID of the base node.
       baseNodeUniqueID=tree%baseNode%uniqueID()
       ! Generate a suitable file name.
       fileName=self%suspendPath//'/suspendedTree_'//baseNodeUniqueID
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
       tree%index=baseNodeUniqueID
    end if
    !$omp critical(universeTransform)
    call self%universeProcessed%pushTree(tree)
    !$omp end critical(universeTransform)
    tree => null()
    return
  end subroutine evolveForestsSuspendTree

  subroutine evolveForestsResumeTree(self,tree)
    !% Resume processing of a tree.
    use :: ISO_Varying_String      , only : varying_string         , operator(//)
    use :: Merger_Tree_Construction, only : mergerTreeStateFromFile
    use :: String_Handling         , only : operator(//)
    implicit none
    class(taskEvolveForests)         , intent(inout) :: self
    type (mergerTree       ), pointer, intent(  out) :: tree
    type (varying_string   )                         :: fileName

    !$omp critical(universeTransform)
    tree => self%universeWaiting%popTree()
    !$omp end critical(universeTransform)
    ! If the tree was suspended to file, restore it now.
    if (.not.self%suspendToRAM.and.associated(tree)) then
       ! Generate the file name.
       fileName=self%suspendPath//'/suspendedTree_'//tree%index
       ! Read the tree from file.
       call mergerTreeStateFromFile(tree,char(fileName),deleteAfterRead=.true.)
    end if
    return
  end subroutine evolveForestsResumeTree
