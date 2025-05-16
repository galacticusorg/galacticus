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

  use :: NBody_Importers, only : nbodyImporter, nbodyImporterClass
  use :: NBody_Operators, only : nbodyOperator, nbodyOperatorClass

  !![
  <task name="taskNBodyAnalyze">
   <description>A task which analyzes N-body simulation data.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskNBodyAnalyze
     !!{
     Implementation of a task which analyzes N-body simulation data.
     !!}
     private
     class  (nbodyImporterClass), pointer :: nbodyImporter_      => null()
     class  (nbodyOperatorClass), pointer :: nbodyOperator_      => null()
     logical                              :: storeBackToImported          , nodeComponentsInitialized=.false.
     ! Pointer to the parameters for this task.
     type   (inputParameters   )          :: parameters
   contains
     final     ::                       nbodyAnalyzeDestructor
     procedure :: perform            => nbodyAnalyzePerform
     procedure :: requiresOutputFile => nbodyAnalyzeRequiresOutputFile
  end type taskNBodyAnalyze

  interface taskNBodyAnalyze
     !!{
     Constructors for the \refClass{taskNBodyAnalyze} task.
     !!}
     module procedure nbodyAnalyzeConstructorParameters
     module procedure nbodyAnalyzeConstructorInternal
  end interface taskNBodyAnalyze

contains

  function nbodyAnalyzeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskNBodyAnalyze} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type   (taskNBodyAnalyze  )                        :: self
    type   (inputParameters   ), intent(inout), target :: parameters
    class  (nbodyImporterClass), pointer               :: nbodyImporter_
    class  (nbodyOperatorClass), pointer               :: nbodyOperator_
    logical                                            :: storeBackToImported
    type   (inputParameters   ), pointer               :: parametersRoot

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
    <inputParameter>
      <name>storeBackToImported</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, computed properties and results will be stored back to the file from which a simulation was imported (assuming it is of HDF5 type).</description>
    </inputParameter>
    <objectBuilder class="nbodyImporter" name="nbodyImporter_" source="parameters"/>
    <objectBuilder class="nbodyOperator" name="nbodyOperator_" source="parameters"/>
    !!]
    self=taskNBodyAnalyze(storeBackToImported,nbodyImporter_,nbodyOperator_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyImporter_"/>
    <objectDestructor name="nbodyOperator_"/>
    !!]
    return
  end function nbodyAnalyzeConstructorParameters

  function nbodyAnalyzeConstructorInternal(storeBackToImported,nbodyImporter_,nbodyOperator_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskNBodyAnalyze} task class which takes a parameter set as input.
    !!}
    implicit none
    type   (taskNBodyAnalyze  )                        :: self
    logical                    , intent(in   )         :: storeBackToImported
    class  (nbodyImporterClass), intent(in   ), target :: nbodyImporter_
    class  (nbodyOperatorClass), intent(in   ), target :: nbodyOperator_
    type   (inputParameters   ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="storeBackToImported, *nbodyImporter_, *nbodyOperator_"/>
    !!]

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    if (.not.self%nbodyImporter_%isHDF5()) self%storeBackToImported=.false.
    return
  end function nbodyAnalyzeConstructorInternal

  subroutine nbodyAnalyzeDestructor(self)
    !!{
    Destructor for the \refClass{taskNBodyAnalyze} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskNBodyAnalyze), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyOperator_"/>
    <objectDestructor name="self%nbodyImporter_"/>
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine nbodyAnalyzeDestructor

  subroutine nbodyAnalyzePerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use :: Display              , only : displayIndent                    , displayUnindent
    use :: Error                , only : errorStatusSuccess
    use :: Output_HDF5          , only : outputFile
    use :: HDF5_Access          , only : hdf5Access
    use :: NBody_Simulation_Data, only : nBodyData
    use :: Node_Components      , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    implicit none
    class    (taskNBodyAnalyze), intent(inout), target       :: self
    integer                    , intent(  out), optional     :: status
    type     (nBodyData       ), allocatable  , dimension(:) :: simulations
    integer                                                  :: i
    character(len=32          )                              :: label
    
    call displayIndent('Begin task: N-body analyze')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Import N-body data.
    call self%nbodyImporter_%import (simulations)
    ! Open analysis groups if necessary.
    !$ call hdf5Access%set()
    do i=1,size(simulations)
       if (.not.self%storeBackToImported) then
          if (simulations(i)%analysis%isOpen()) call simulations(i)%analysis%close()
          write (label,'(a,i4.4)') 'simulation',i
          simulations(i)%analysis=outputFile%openGroup(label)
          call simulations(i)%analysis%writeAttribute(simulations(i)%label,'label')
       end if
    end do
    !$ call hdf5Access%unset()
    ! Operate on the N-body data.
    call self%nbodyOperator_%operate(simulations)
    ! Close the analysis group.
    !$ call hdf5Access%set()
    do i=1,size(simulations)
       if (simulations(i)%analysis%isOpen()) call simulations(i)%analysis%close()
    end do
    !$ call hdf5Access%unset()
    ! Done.
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: N-body analyze' )
    return
  end subroutine nbodyAnalyzePerform

  logical function nbodyAnalyzeRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskNBodyAnalyze), intent(inout) :: self

    nbodyAnalyzeRequiresOutputFile=.not.self%storeBackToImported
    return
  end function nbodyAnalyzeRequiresOutputFile
