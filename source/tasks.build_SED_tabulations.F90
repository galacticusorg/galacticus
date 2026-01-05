!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Contains a class which implements a task to pre-build tabulations needed for SED calculations.
  !!}

  use :: Node_Property_Extractors, only : nodePropertyExtractorClass
  use :: Output_Times            , only : outputTimesClass
  use :: Star_Formation_Histories, only : starFormationHistoryClass

  !![
  <task name="taskBuildSEDTabulations">
   <description>A task which pre-builds tabulations needed for SED calculations.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildSEDTabulations
     !!{
     Implementation of a task which pre-builds tabulations needed for SED calculations.
     !!}
     private
     class  (nodePropertyExtractorClass), pointer :: nodePropertyExtractor_    => null()
     class  (outputTimesClass          ), pointer :: outputTimes_              => null()
     class  (starFormationHistoryClass ), pointer :: starFormationHistory_     => null()
     type   (inputParameters           )          :: parameters
     logical                                      :: nodeComponentsInitialized =  .false.
   contains
     final     ::                       buildSEDDestructor
     procedure :: perform            => buildSEDTabulationsPerform
     procedure :: requiresOutputFile => buildSEDTabulationsRequiresOutputFile
  end type taskBuildSEDTabulations

  interface taskBuildSEDTabulations
     !!{
     Constructors for the \refClass{taskBuildSEDTabulations} task.
     !!}
     module procedure buildSEDTabulationsParameters
     module procedure buildSEDTabulationsInternal
  end interface taskBuildSEDTabulations

contains

  function buildSEDTabulationsParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildSEDTabulations} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type (taskBuildSEDTabulations   )                        :: self
    type (inputParameters           ), intent(inout), target :: parameters
    class(nodePropertyExtractorClass), pointer               :: nodePropertyExtractor_
    class(outputTimesClass          ), pointer               :: outputTimes_
    class(starFormationHistoryClass ), pointer               :: starFormationHistory_
    type (inputParameters           ), pointer               :: parametersRoot

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => parameters
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    self%nodeComponentsInitialized=.true.
    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    <objectBuilder class="starFormationHistory"  name="starFormationHistory_"  source="parameters"/>
    !!]
    self=taskBuildSEDTabulations(nodePropertyExtractor_,starFormationHistory_,outputTimes_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    <objectDestructor name="outputTimes_"          />
    <objectDestructor name="starFormationHistory_" />
    !!]
    return
  end function buildSEDTabulationsParameters

  function buildSEDTabulationsInternal(nodePropertyExtractor_,starFormationHistory_,outputTimes_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildSEDTabulations} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildSEDTabulations    )                        :: self
    class(nodePropertyExtractorClass), intent(in   ), target :: nodePropertyExtractor_
    class(outputTimesClass          ), intent(in   ), target :: outputTimes_
    class(starFormationHistoryClass ), intent(in   ), target :: starFormationHistory_
    type (inputParameters           ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="*nodePropertyExtractor_, *outputTimes_, *starFormationHistory_"/>
    !!]

    self%parameters=inputParameters(parameters)
    return
  end function buildSEDTabulationsInternal

  subroutine buildSEDDestructor(self)
    !!{
    Destructor for the \refClass{taskBuildSEDTabulations} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskBuildSEDTabulations), intent(inout) :: self
    
    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%starFormationHistory_" />
    <objectDestructor name="self%outputTimes_"          />
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine buildSEDDestructor

  subroutine buildSEDTabulationsPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display                 , only : displayIndent                    , displayUnindent
    use :: Error                   , only : Error_Report                     , errorStatusSuccess
    use :: Histories               , only : history
    use :: Galacticus_Nodes        , only : mergerTree                       , nodeComponentBasic                 , nodeComponentDisk, nodeComponentSpheroid, &
          &                                 nodeComponentNSC                 , treeNode
    use :: Poly_Ranks              , only : assignment(=)                    , polyRankDouble
    use :: Multi_Counters          , only : multiCounter
    use :: Node_Property_Extractors, only : nodePropertyExtractorMulti       , nodePropertyExtractorSED
    use :: Node_Components         , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use :: Locks                   , only : ompLock
#ifdef USEMPI
    use :: MPI_Utilities           , only : mpiBarrier                       , mpiSelf
#endif
    implicit none
    class           (taskBuildSEDTabulations), intent(inout), target         :: self
    integer                                  , intent(  out), optional       :: status
    double precision                         , allocatable  , dimension(:,:) :: doubleArray
    type            (polyRankDouble         ), allocatable  , dimension(:  ) :: doubleProperties
    type            (mergerTree             ), target                        :: tree
    type            (treeNode               ), pointer                       :: node
    class           (nodeComponentBasic     ), pointer                       :: basic
    class           (nodeComponentDisk      ), pointer                       :: disk
    class           (nodeComponentSpheroid  ), pointer                       :: spheroid
    class           (nodeComponentNSC       ), pointer                       :: nuclearStarCluster
    double precision                         , parameter                     :: mass                                  =1.0d+12
    double precision                         , parameter                     :: epsilon                               =1.0d-06
    type            (history                )                                :: starFormationHistoryDisk                      , starFormationHistorySpheroid, &
         &                                                                      starFormationHistorynuclearStarCluster
    integer         (c_size_t               )                                :: i
    double precision                                                         :: time
    type            (multiCounter           )                                :: instance
    type            (ompLock                )                                :: treeLock
    !$GLC attributes unused :: self

    call displayIndent ('Begin task: build SED tabulations')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Build a node and components.
    node                        => treeNode         (                 )
    tree              %nodeBase => node
    node              %hostTree => tree
    disk                        => node    %disk    (autoCreate=.true.)
    spheroid                    => node    %spheroid(autoCreate=.true.)
    nuclearStarCluster          => node    %NSC     (autoCreate=.true.)
    basic                       => node    %basic   (autoCreate=.true.)
    call tree%properties%initialize()
    ! Initialize a single instance.
    instance=multiCounter([1_c_size_t])
    ! Initialize a tree output lock - this is required by the star formation history output functions, even though we don't
    ! actually need any locking here.
    treeLock=ompLock()
    ! Choose a beginning time. Set this to slightly before the first output time to ensure that the initial set of times assigned
    ! to each history corresponds to the first (and not the second) output time).
    time   =+self%outputTimes_%time(1_c_size_t) &
         &  *(                                  &
         &    +1.0d0                            &
         &    -epsilon                          &
         &   )
    call basic%timeSet            (time)
    call basic%timeLastIsolatedSet(time)
    call basic%massSet            (mass)
    ! Create the star formation histories.
    call self              %starFormationHistory_  %create(node,starFormationHistoryDisk              ,time)
    call self              %starFormationHistory_  %create(node,starFormationHistorySpheroid          ,time)
    call self              %starFormationHistory_  %create(node,starFormationHistorynuclearStarCluster,time)

    starFormationHistoryDisk    %data=1.0d0
    starFormationHistorySpheroid%data=1.0d0
    call disk              %starFormationHistorySet       (     starFormationHistoryDisk                   )
    call spheroid          %starFormationHistorySet       (     starFormationHistorySpheroid               )
    call nuclearStarCluster%starFormationHistorySet       (     starFormationHistorynuclearStarCluster     )

    ! Iterate over output times.
    do i=1_c_size_t,self%outputTimes_%count()
       time=self%outputTimes_%time(i)
       call basic   %timeSet            (time)
       call basic   %timeLastIsolatedSet(time)
       call instance%reset              (    )
       if (.not.instance%increment()) call Error_Report('failed to increment instance'//{introspection:location})
#ifdef USEMPI
       if (mod(i,mpiSelf%count()) == mpiSelf%rank()) then
#endif
          select type (extractor_ => self%nodePropertyExtractor_)
          class is (nodePropertyExtractorSED  )
             ! SED property extractor - extract and store the values.
             doubleArray     =extractor_%extract      (node,time,instance)
             deallocate(doubleArray     )
          class is (nodePropertyExtractorMulti)
             ! Multi property extractor - extract and store the values.
             doubleProperties=extractor_%extractDouble(node,time,instance)
             deallocate(doubleProperties)
          end select
#ifdef USEMPI
       end if
#endif
       ! Output star formation history, which also triggers update of the history.
       call starFormationHistoryDisk              %destroy      (                                                   &
            &                                                   )
       call starFormationHistorySpheroid          %destroy      (                                                   &
            &                                                   )
       call starFormationHistorynuclearStarCluster%destroy      (                                                   &
            &                                                   )
       starFormationHistoryDisk              =disk              %starFormationHistory()
       starFormationHistorySpheroid          =spheroid          %starFormationHistory()
       starFormationHistorynuclearStarCluster=nuclearStarCluster%starFormationHistory()

       call self%starFormationHistory_%update                   (                                                             &
            &                                                    node                =node                                  , &
            &                                                    indexOutput         =i                                     , &
            &                                                    historyStarFormation=starFormationHistoryDisk                &
            &                                                   )
       call self%starFormationHistory_%update                   (                                                             &
            &                                                    node                =node                                  , &
            &                                                    indexOutput         =i                                     , &
            &                                                    historyStarFormation=starFormationHistorySpheroid            &
            &                                                   )
       call self%starFormationHistory_%update                   (                                                             &
            &                                                    node                =node                                  , &
            &                                                    indexOutput         =i                                     , &
            &                                                    historyStarFormation=starFormationHistorynuclearStarCluster  &
            &                                                   )
       starFormationHistoryDisk    %data=1.0d0
       starFormationHistorySpheroid%data=1.0d0
       starFormationHistorynuclearStarCluster     %data=1.0d0

       call disk                      %starFormationHistorySet(                                                               &
            &                                                  starFormationHistoryDisk                                       &
            &                                                 )
       call spheroid                  %starFormationHistorySet(                                                               &
            &                                                  starFormationHistorySpheroid                                   &
            &                                                 )
       call nuclearStarCluster        %starFormationHistorySet(                                                               &
            &                                                  starFormationHistorynuclearStarCluster                         &
            &                                                 )
    end do
    call node%destroy()
    deallocate(node)
    if (present(status)) status=errorStatusSuccess
    call Node_Components_Thread_Uninitialize()
    call displayUnindent('Done task: build SED tabulations')
    return
  end subroutine buildSEDTabulationsPerform

  logical function buildSEDTabulationsRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildSEDTabulations ), intent(inout) :: self
    !$GLC attributes unused :: self

    buildSEDTabulationsRequiresOutputFile=.false.
    return
  end function buildSEDTabulationsRequiresOutputFile
