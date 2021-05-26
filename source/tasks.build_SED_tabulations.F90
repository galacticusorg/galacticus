!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Contains a class which implements a task to pre-build tabulations needed for SED calculations.

  use :: Node_Property_Extractors, only : nodePropertyExtractorClass
  use :: Output_Times            , only : outputTimesClass
  use :: Star_Formation_Histories, only : starFormationHistoryClass

  !# <task name="taskBuildSEDTabulations">
  !#  <description>A task which pre-builds tabulations needed for SED calculations.</description>
  !# </task>
  type, extends(taskClass) :: taskBuildSEDTabulations
     !% Implementation of a task which pre-builds tabulations needed for SED calculations.
     private
     class(nodePropertyExtractorClass), pointer :: nodePropertyExtractor_ => null()
     class(outputTimesClass          ), pointer :: outputTimes_           => null()
     class(starFormationHistoryClass ), pointer :: starFormationHistory_  => null()
     type (inputParameters           )          :: parameters
   contains
     final     ::                       buildSEDDestructor
     procedure :: perform            => buildSEDTabulationsPerform
     procedure :: requiresOutputFile => buildSEDTabulationsRequiresOutputFile
  end type taskBuildSEDTabulations

  interface taskBuildSEDTabulations
     !% Constructors for the {\normalfont \ttfamily buildSEDTabulations} task.
     module procedure buildSEDTabulationsParameters
     module procedure buildSEDTabulationsInternal
  end interface taskBuildSEDTabulations

contains

  function buildSEDTabulationsParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily buildSEDTabulations} task class which takes a parameter set as input.
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
    !# <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !# <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    !# <objectBuilder class="starFormationHistory"  name="starFormationHistory_"  source="parameters"/>
    self=taskBuildSEDTabulations(nodePropertyExtractor_,starFormationHistory_,outputTimes_,parametersRoot)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="nodePropertyExtractor_"/>
    !# <objectDestructor name="outputTimes_"          />
    !# <objectDestructor name="starFormationHistory_" />
    return
  end function buildSEDTabulationsParameters

  function buildSEDTabulationsInternal(nodePropertyExtractor_,starFormationHistory_,outputTimes_,parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily buildSEDTabulations} task class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(taskBuildSEDTabulations    )                        :: self
    class(nodePropertyExtractorClass), intent(in   ), target :: nodePropertyExtractor_
    class(outputTimesClass          ), intent(in   ), target :: outputTimes_
    class(starFormationHistoryClass ), intent(in   ), target :: starFormationHistory_
    type (inputParameters           ), intent(in   ), target :: parameters
    !# <constructorAssign variables="*nodePropertyExtractor_, *outputTimes_, *starFormationHistory_"/>

    self%parameters=inputParameters(parameters)
    return
  end function buildSEDTabulationsInternal

  subroutine buildSEDDestructor(self)
    !% Destructor for the {\normalfont \ttfamily buildSEDTabulations} task class.
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskBuildSEDTabulations), intent(inout) :: self
    
    !# <objectDestructor name="self%nodePropertyExtractor_"/>
    !# <objectDestructor name="self%starFormationHistory_" />
    !# <objectDestructor name="self%outputTimes_"          />
    call Node_Components_Uninitialize()
    return
  end subroutine buildSEDDestructor

  subroutine buildSEDTabulationsPerform(self,status)
    !% Builds the tabulation.
    use :: Display                   , only : displayIndent                    , displayUnindent
    use :: Galacticus_Error          , only : errorStatusSuccess               , Galacticus_Error_Report
    use :: Histories                 , only : history
    use :: Galacticus_Nodes          , only : treeNode                         , nodeComponentBasic                 , nodeComponentDisk, nodeComponentSpheroid, &
         &                                    mergerTree
    use :: Poly_Ranks                , only : polyRankDouble                   , assignment(=)
    use :: Multi_Counters            , only : multiCounter
    use :: Node_Property_Extractors  , only : nodePropertyExtractorSED         , nodePropertyExtractorMulti
    use :: Galactic_Structure_Options, only : componentTypeDisk                , componentTypeSpheroid
    use :: Node_Components           , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
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
    double precision                         , parameter                     :: mass                    =1.0d12
    type            (history                )                                :: starFormationHistoryDisk       , starFormationHistorySpheroid
    integer         (c_size_t               )                                :: i
    double precision                                                         :: time
    type            (multiCounter           )                                :: instance
    !$GLC attributes unused :: self

    call displayIndent ('Begin task: build SED tabulations')
    ! Call routines to perform initializations which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Build a node and components.
    node     => treeNode         (                 )
    tree    %baseNode => node
    node%hostTree => tree
    disk     => node    %disk    (autoCreate=.true.)
    spheroid => node    %spheroid(autoCreate=.true.)
    basic    => node    %basic   (autoCreate=.true.)
    ! Initialize a single instance.
    instance=multiCounter([1_c_size_t])
    ! Choose a beginning time.
    time=self%outputTimes_%time(1_c_size_t)
    call basic%timeSet            (time)
    call basic%timeLastIsolatedSet(time)
    call basic%massSet            (mass)
    ! Create the star formation histories.
    call self    %starFormationHistory_  %create(node,starFormationHistoryDisk    ,time)
    call self    %starFormationHistory_  %create(node,starFormationHistorySpheroid,time)
    starFormationHistoryDisk    %data=1.0d0
    starFormationHistorySpheroid%data=1.0d0
    call disk    %starFormationHistorySet       (     starFormationHistoryDisk         )
    call spheroid%starFormationHistorySet       (     starFormationHistorySpheroid     )
    ! Iterate over output times.
    do i=1_c_size_t,self%outputTimes_%count()
       time=self%outputTimes_%time(i)
       call basic   %timeSet            (time)
       call basic   %timeLastIsolatedSet(time)
       call instance%reset              (    )
       if (.not.instance%increment()) call Galacticus_Error_Report('failed to increment instance'//{introspection:location})
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
       ! Output star formation history, which also triggers update of the history.
       call starFormationHistoryDisk    %destroy                (                                                   &
            &                                                   )
       call starFormationHistorySpheroid%destroy                (                                                   &
            &                                                   )
       starFormationHistoryDisk    =disk    %starFormationHistory()
       starFormationHistorySpheroid=spheroid%starFormationHistory()
       call self%starFormationHistory_%output                   (                                                   &
            &                                                    node                =node                        , &
            &                                                    nodePassesFilter    =.false.                     , &
            &                                                    historyStarFormation=starFormationHistoryDisk    , &
            &                                                    indexOutput         =i                           , &
            &                                                    indexTree           =1_c_size_t                  , &
            &                                                    componentType       =componentTypeDisk             &
            &                                                   )
       call self%starFormationHistory_%output                   (                                                   &
            &                                                    node                =node                        , &
            &                                                    nodePassesFilter    =.false.                     , &
            &                                                    historyStarFormation=starFormationHistorySpheroid, &
            &                                                    indexOutput         =i                           , &
            &                                                    indexTree           =1_c_size_t                  , &
            &                                                    componentType       =componentTypeSpheroid         &
            &                                                   )
       starFormationHistoryDisk    %data=1.0d0
       starFormationHistorySpheroid%data=1.0d0
       call disk                        %starFormationHistorySet(                                                   &
            &                                                    starFormationHistoryDisk                           &
            &                                                   )
       call spheroid                    %starFormationHistorySet(                                                   &
            &                                                    starFormationHistorySpheroid                       &
            &                                                   )
    end do
    call node%destroy()
    deallocate(node)
    if (present(status)) status=errorStatusSuccess
    call Node_Components_Thread_Uninitialize()
    call displayUnindent('Done task: build SED tabulations')
    return
  end subroutine buildSEDTabulationsPerform

  logical function buildSEDTabulationsRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskBuildSEDTabulations ), intent(inout) :: self
    !$GLC attributes unused :: self

    buildSEDTabulationsRequiresOutputFile=.false.
    return
  end function buildSEDTabulationsRequiresOutputFile
