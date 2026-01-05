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

  use :: Input_Parameters          , only : inputParameters
  use :: Merger_Tree_Outputters    , only : mergerTreeOutputterClass
  use :: Merger_Tree_Seeds         , only : mergerTreeSeedsClass
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Output_Times              , only : outputTimesClass
  use :: Merger_Tree_Construction  , only : mergerTreeConstructorClass
  use :: Merger_Tree_Initialization, only : mergerTreeInitializorClass
  use :: Merger_Tree_Operators     , only : mergerTreeOperatorClass
  use :: Merger_Trees_Evolve       , only : mergerTreeEvolverClass
  use :: Nodes_Operators           , only : nodeOperatorClass
  use :: Universe_Operators        , only : universeOperatorClass

  !![
  <task name="taskPostprocessForests">
    <description>
      A task which postprocesses galaxies within a set of merger tree forests. This task assumes that a prior model was run, with
      raw forest data written to file using the \refClass{mergerTreeOutputterFullState} class. The name of that file is specified
      via the {\normalfont \ttfamily fileName} parameter. Forests data will be re-read, and re-output. Note that you should use
      the \emph{exact same} parameter file (other than changing the {\normalfont \ttfamily task}, and possibly removing the use of
      the \refClass{mergerTreeOutputterFullState} outputter) as was used to run the original model. This ensures that the raw data
      structures read from the file follow the same format as was used to write them. Also note that forests are not guaranteed to
      be output in the same order as in the original model if OpenMP parallelism is used. If the same order is required, it is
      recommend to run the postprocessing after setting the environment variable {\normalfont \ttfamily OMP\_NUM\_THREADS=1}.
    </description>
  </task>
  !!]
  type, extends(taskClass) :: taskPostprocessForests
     !!{
     Implementation of a task which postprocesses galaxies within a set of merger tree forests.
     !!}
     private
     type   (inputParameters            ), pointer :: parameters                => null()
     class  (mergerTreeConstructorClass ), pointer :: mergerTreeConstructor_    => null()
     class  (mergerTreeOperatorClass    ), pointer :: mergerTreeOperator_       => null()
     class  (mergerTreeEvolverClass     ), pointer :: mergerTreeEvolver_        => null()
     class  (mergerTreeOutputterClass   ), pointer :: mergerTreeOutputter_      => null()
     class  (mergerTreeInitializorClass ), pointer :: mergerTreeInitializor_    => null()
     class  (nodeOperatorClass          ), pointer :: nodeOperator_             => null()
     class  (evolveForestsWorkShareClass), pointer :: evolveForestsWorkShare_   => null()
     class  (outputTimesClass           ), pointer :: outputTimes_              => null()
     class  (universeOperatorClass      ), pointer :: universeOperator_         => null()
     class  (randomNumberGeneratorClass ), pointer :: randomNumberGenerator_    => null()
     class  (mergerTreeSeedsClass       ), pointer :: mergerTreeSeeds_          => null()
     type   (varying_string             )          :: fileName
     logical                                       :: nodeComponentsInitialized =  .false.
 contains
     final     ::            postprocessForestsDestructor
     procedure :: perform => postprocessForestsPerform
  end type taskPostprocessForests

  interface taskPostprocessForests
     !!{
     Constructors for the \refClass{taskPostprocessForests} task.
     !!}
     module procedure postprocessForestsConstructorParameters
     module procedure postprocessForestsConstructorInternal
  end interface taskPostprocessForests

  ! Copies of objects used by each thread.
  class(mergerTreeOutputterClass), pointer :: mergerTreeOutputter_ => null()
  !$omp threadprivate(mergerTreeOutputter_)

contains

  function postprocessForestsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskPostprocessForests} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type (taskPostprocessForests     )                        :: self
    type (inputParameters            ), intent(inout), target :: parameters
    class(mergerTreeOperatorClass    ), pointer               :: mergerTreeOperator_
    class(nodeOperatorClass          ), pointer               :: nodeOperator_
    class(evolveForestsWorkShareClass), pointer               :: evolveForestsWorkShare_
    class(mergerTreeConstructorClass ), pointer               :: mergerTreeConstructor_
    class(outputTimesClass           ), pointer               :: outputTimes_
    class(universeOperatorClass      ), pointer               :: universeOperator_
    class(mergerTreeEvolverClass     ), pointer               :: mergerTreeEvolver_
    class(mergerTreeOutputterClass   ), pointer               :: mergerTreeOutputter_
    class(mergerTreeInitializorClass ), pointer               :: mergerTreeInitializor_
    class(randomNumberGeneratorClass ), pointer               :: randomNumberGenerator_
    class(mergerTreeSeedsClass       ), pointer               :: mergerTreeSeeds_
    type (inputParameters            ), pointer               :: parametersRoot
    type (varying_string             )                        :: fileName

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
    ! Note that we build the same set of objects as used by the `evolveForests` tasks, even though many of them are not used. This
    ! ensures that the same set of metaProperties are added, which ensures that the raw format deserialization is consistent with
    ! the raw format serialization used in the original model.
    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file from which forests should be read.</description>
      <source>parameters</source>
    </inputParameter>
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
       self=taskPostprocessForests(fileName,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parametersRoot)
    else
       self=taskPostprocessForests(fileName,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parameters    )
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
  end function postprocessForestsConstructorParameters

  function postprocessForestsConstructorInternal(fileName,mergerTreeConstructor_,mergerTreeOperator_,nodeOperator_,evolveForestsWorkShare_,outputTimes_,universeOperator_,mergerTreeEvolver_,mergerTreeOutputter_,mergerTreeInitializor_,randomNumberGenerator_,mergerTreeSeeds_,parameters) result(self)
    !!{
    Internal constructor for the \refClass{taskPostprocessForests} task class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    type (taskPostprocessForests     )                        :: self
    type (varying_string             ), intent(in   )         :: fileName
    type (inputParameters            ), intent(in   ), target :: parameters
    class(mergerTreeConstructorClass ), intent(in   ), target :: mergerTreeConstructor_
    class(mergerTreeOperatorClass    ), intent(in   ), target :: mergerTreeOperator_
    class(nodeOperatorClass          ), intent(in   ), target :: nodeOperator_
    class(evolveForestsWorkShareClass), intent(in   ), target :: evolveForestsWorkShare_
    class(outputTimesClass           ), intent(in   ), target :: outputTimes_
    class(universeOperatorClass      ), intent(in   ), target :: universeOperator_
    class(mergerTreeEvolverClass     ), intent(in   ), target :: mergerTreeEvolver_
    class(mergerTreeOutputterClass   ), intent(in   ), target :: mergerTreeOutputter_
    class(mergerTreeInitializorClass ), intent(in   ), target :: mergerTreeInitializor_
    class(randomNumberGeneratorClass ), intent(in   ), target :: randomNumberGenerator_
    class(mergerTreeSeedsClass       ), intent(in   ), target :: mergerTreeSeeds_
    !![
    <constructorAssign variables="fileName, *mergerTreeConstructor_, *mergerTreeOperator_, *nodeOperator_, *evolveForestsWorkShare_, *outputTimes_, *universeOperator_, *mergerTreeEvolver_, *mergerTreeOutputter_, *mergerTreeInitializor_, *randomNumberGenerator_, *mergerTreeSeeds_, *parameters"/>
    !!]

    return 
  end function postprocessForestsConstructorInternal

  subroutine postprocessForestsDestructor(self)
    !!{
    Destructor for the \refClass{taskPostprocessForests} task class.
    !!}
    use :: Node_Components , only : Node_Components_Uninitialize
    use :: Galacticus_Nodes, only : nodeClassHierarchyFinalize
    implicit none
    type(taskPostprocessForests), intent(inout) :: self

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
    if (self%nodeComponentsInitialized) then
       call Node_Components_Uninitialize()
       call nodeClassHierarchyFinalize  ()
    end if
    return
  end subroutine postprocessForestsDestructor

  subroutine postprocessForestsPerform(self,status)
    !!{
    Postprocesses the complete set of merger trees as specified.
    !!}
    use, intrinsic :: ISO_C_Binding           , only : c_size_t
    use            :: Display                 , only : displayIndent                    , displayMessage                     , displayUnindent
    use            :: Galacticus_Nodes        , only : mergerTree                       , nodeComponentBasic                 , treeNode
    use            :: Error                   , only : errorStatusSuccess
    use            :: Node_Components         , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use            :: Merger_Tree_Construction, only : mergerTreeStateFromFile
    use            :: File_Utilities          , only : File_Lock                        , File_Unlock                        , lockDescriptor
    implicit none
    class           (taskPostprocessForests), intent(inout), target   :: self
    integer                                 , intent(  out), optional :: status
    type            (mergerTree            ), pointer      , save     :: tree
    !$omp threadprivate(tree)
    type            (inputParameters       ), allocatable  , save     :: parameters
    !$omp threadprivate(parameters)
    type            (lockDescriptor        )               , save     :: fileLock
    !$omp threadprivate(fileLock)
    integer                                                , save     :: statusRead
    !$omp threadprivate(statusRead)
    integer         (c_size_t              )               , save     :: indexOutput
    !$omp threadprivate(indexOutput)
    integer                                                           :: fileUnit
    logical                                                           :: finished
    
    call displayIndent('Begin task: merger tree postprocess')
    ! Set status to success by default.
    if (present(status)) status=errorStatusSuccess
    ! Open the file.
    open(newUnit=fileUnit,file=char(self%fileName),status='old',form='unformatted')
    ! Begin parallel post-processing of trees until all work is done.
    finished=.false.
    !$omp parallel
    allocate(mergerTreeOutputter_,mold=self%mergerTreeOutputter_)
    !$omp critical(postprocessForestsDeepCopy)
    !![
    <deepCopyReset variables="self%mergerTreeOutputter_"/>
    <deepCopy source="self%mergerTreeOutputter_" destination="mergerTreeOutputter_"/>
    <deepCopyFinalize variables="mergerTreeOutputter_"/>
    !!]
    !$omp end critical(postprocessForestsDeepCopy)
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call parameters%parametersGroupCopy(self%parameters)
    call Node_Components_Thread_Initialize(parameters)
    !$omp barrier
    ! Begin loop to read and post-process trees.
    do while (.not.finished)
       ! Tree read is done in a critical section to avoid thread contentions on the file.
       !$omp critical(taskPostprocessForests)
       if (.not.finished) then
          allocate(tree)
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
          read (fileUnit,iostat=statusRead) indexOutput
          if (statusRead == errorStatusSuccess) call mergerTreeStateFromFile(tree,fileUnit,self%randomNumberGenerator_,self%mergerTreeSeeds_,statusRead)
          call File_Unlock(fileLock)
          if (statusRead /= errorStatusSuccess) then
             deallocate(tree)
             finished=.true.
          end if
       end if
       !$omp end critical(taskPostprocessForests)
       if (finished) cycle
       ! Postprocess the tree by calling the outputter.
       call mergerTreeOutputter_%outputTree(tree,indexOutput,self%outputTimes_%time(indexOutput))
       ! Destroy the tree.
       call tree%destroy()
       deallocate(tree)
       tree => null()
    end do
    !$omp barrier
    ! Reduce outputs back into the original outputter object.
    call mergerTreeOutputter_%reduce(self%mergerTreeOutputter_)
    ! Explicitly deallocate objects.
    !![
    <objectDestructor name="mergerTreeOutputter_"/>
    !!]
    call Node_Components_Thread_Uninitialize()
    !$omp barrier
    !$omp critical(evolveForestReset)
    call parameters%reset()
    !$omp end critical(evolveForestReset)
    !$omp barrier
    deallocate(parameters)
    !$omp end parallel
    ! Finalize outputs.
    call self%mergerTreeOutputter_%finalize()
    ! Close the file.
    close(fileUnit)
    call displayUnindent('Done task: merger tree postprocessing')
    return
  end subroutine postprocessForestsPerform
