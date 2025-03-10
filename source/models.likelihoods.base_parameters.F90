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
  Implementation of a posterior sampling likelihood class which allows arbitrary modification of a base parameter object.
  !!}

  use :: Input_Parameters, only : inputParameter, inputParameters

  type :: parameterList
     !!{
     Type used to maintain a list of pointers to parameters to be modified.
     !!}
     type   (inputParameter), pointer :: parameter_   => null()
     integer                          :: indexElement
     type   (varying_string)          :: definition
     logical                          :: resolved
  end type parameterList

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodBaseParameters" abstract="yes">
   <description>A posterior sampling likelihood class which allows arbitrary modification of a base parameter object.</description>
  </posteriorSampleLikelihood>
  !!]
  type, abstract, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodBaseParameters
     !!{
     Implementation of a posterior sampling likelihood class which allows arbitrary modification of a base parameter object.
     !!}
     private
     type   (varying_string )                            :: baseParametersFileName
     type   (varying_string ), dimension(:), allocatable :: changeParametersFileNames
     type   (inputParameters), pointer                   :: parametersModel           => null()
     type   (parameterList  ), dimension(:), allocatable :: modelParametersActive_              , modelParametersInactive_
     logical                                             :: reportFileName            =  .false., reportState             =.false.
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize pointers into the base parameters object."/>
       <method method="update"     description="Update values in the base parameters object."        />
     </methods>
     !!]
     procedure :: initialize => baseParametersInitialize
     procedure :: update     => baseParametersUpdate
  end type posteriorSampleLikelihoodBaseParameters

contains

  subroutine baseParametersInitialize(self,modelParametersActive_,modelParametersInactive_)
    !!{
    Initialize pointers into the base parameters object.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char              , operator(//), var_str
    use :: Kind_Numbers      , only : kind_int8
    use :: String_Handling   , only : String_Count_Words, String_Join , String_Split_Words, operator(//)
    implicit none
    class    (posteriorSampleLikelihoodBaseParameters), intent(inout)                 :: self
    type     (modelParameterList                     ), intent(in   ), dimension(:  ) :: modelParametersActive_, modelParametersInactive_
    type     (varying_string                         ), allocatable  , dimension(:  ) :: parameterNames
    integer                                                                           :: i                     , j                       , &
         &                                                                               instance              , indexElement            , &
         &                                                                               parameterCount        , ioStatus                , &
         &                                                                               k
    type     (inputParameters                        ), pointer                       :: parameters_           , subParameters_
    character(len=1024                               )                                :: labelIndex            , labelValue              , &
         &                                                                               parameterValue

    ! On first call we must build pointers to all parameter nodes which will be modified as a function of chain state.
    if (.not.allocated(self%modelParametersActive_)) then
       allocate(self%modelParametersActive_(size(modelParametersActive_)))
       do i=1,size(modelParametersActive_)
          ! Check for duplicated parameters.
          do j=1,size(modelParametersActive_)
             if     (                                                                                                                                                &
                  &   modelParametersActive_(i)%modelParameter_%name() == modelParametersActive_(j)%modelParameter_%name()                                           &
                  &  .and.                                                                                                                                           &
                  &                          i                         /=                        j                                                                   &
                  & ) call Error_Report("duplicated active parameter name '"//char(modelParametersActive_(i)%modelParameter_%name())//"'"//{introspection:location})
          end do
          parameterCount=String_Count_Words(char(modelParametersActive_(i)%modelParameter_%name()),"/")
          allocate(parameterNames(parameterCount))
          call String_Split_Words(parameterNames,char(modelParametersActive_(i)%modelParameter_%name()),"/")
          allocate(parameters_)
          parameters_=inputParameters(self%parametersModel)
          do j=1,parameterCount
             instance    =1
             indexElement=0
             if (index(parameterNames(j),"[") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"[")+1,index(parameterNames(j),"]")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"[")-1)
                read (labelIndex,*,iostat=ioStatus) instance
                if (ioStatus == 0) then
                   ! Integer indexing - nothing further to do.
                else
                   ! Non-integer indexing, expect an XPath expression.
                   if (index(labelIndex,"@value='") == 1) then
                      labelValue=labelIndex(9:index(labelIndex(9:),"'")+8-1)
                      do k=1,parameters_%copiesCount(char(parameterNames(j)),requireValue=.true.)
                         call parameters_%value(char(parameterNames(j)),parameterValue,copyInstance=k)
                         if (trim(parameterValue) == trim(labelValue)) then
                            instance=k
                            exit
                         end if
                      end do
                      if (instance == 0) call Error_Report("unable to find matching parameter `"//parameterNames(j)//"` with value `"//trim(parameterValue)//"`"//{introspection:location})                      
                   else
                      call Error_Report("expected a XPath attribute selector of the form `@value='...'`, but got `"//trim(labelIndex)//"`"//{introspection:location})
                   end if
                end if
             else if (index(parameterNames(j),"{") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"{")+1,index(parameterNames(j),"}")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"{")-1)
                read (labelIndex,*) indexElement
             end if
             if (j == parameterCount) then
                ! This is the final parameter - so get and store a pointer to its node.
                self%modelParametersActive_(i)%parameter_   => parameters_%node         (char(parameterNames(j)),requireValue=.true. ,copyInstance=instance)
                self%modelParametersActive_(i)%indexElement =  indexElement
                self%modelParametersActive_(i)%definition   =  modelParametersActive_(i)%modelParameter_%name()
             else
                ! This is an intermediate parameter, get the appropriate sub-parameters.
                allocate  (subParameters_)
                subParameters_                              =  parameters_%subParameters(char(parameterNames(j)),requireValue=.false.,copyInstance=instance)
                deallocate(   parameters_)
                allocate  (   parameters_)
                parameters_                                 =  inputParameters(subParameters_)
                deallocate(subParameters_)
             end if
          end do
          deallocate(parameters_   )
          deallocate(parameterNames)
       end do
    end if
    if (.not.allocated(self%modelParametersInactive_)) then
       allocate(self%modelParametersInactive_(size(modelParametersInactive_)))
       do i=1,size(modelParametersInactive_)
          ! Check for duplicated parameters.
          do j=1,size(modelParametersInactive_)
             if     (                                                                                                                                                  &
                  &   modelParametersInactive_(i)%modelParameter_%name() == modelParametersInactive_(j)%modelParameter_%name()                                         &
                  &  .and.                                                                                                                                             &
                  &                            i                         /=                          j                                                                 &
                  & ) call Error_Report("duplicated active parameter name '"//char(modelParametersInactive_(i)%modelParameter_%name())//"'"//{introspection:location})
          end do
          parameterCount=String_Count_Words(char(modelParametersInactive_(i)%modelParameter_%name()),"/")
          allocate(parameterNames(parameterCount))
          call String_Split_Words(parameterNames,char(modelParametersInactive_(i)%modelParameter_%name()),"/")
          allocate(parameters_)
          parameters_=inputParameters(self%parametersModel)
          do j=1,parameterCount
             instance    =1
             indexElement=0
              if (index(parameterNames(j),"[") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"[")+1,index(parameterNames(j),"]")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"[")-1)
                read (labelIndex,*,iostat=ioStatus) instance
                if (ioStatus == 0) then
                   ! Integer indexing - nothing further to do.
                else
                   ! Non-integer indexing, expect an XPath expression.
                   if (index(labelIndex,"@value='") == 1) then
                      labelValue=labelIndex(9:index(labelIndex(9:),"'")+8-1)
                      do k=1,parameters_%copiesCount(char(parameterNames(j)),requireValue=.true.)
                         call parameters_%value(char(parameterNames(j)),parameterValue,copyInstance=k)
                         if (trim(parameterValue) == trim(labelValue)) then
                            instance=k
                            exit
                         end if
                      end do
                      if (instance == 0) call Error_Report("unable to find matching parameter `"//parameterNames(j)//"` with value `"//trim(parameterValue)//"`"//{introspection:location})
                   else
                      call Error_Report("expected a XPath attribute selector of the form `@value='...'`, but got `"//trim(labelIndex)//"`"//{introspection:location})
                   end if
                end if
             else if (index(parameterNames(j),"{") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"{")+1,index(parameterNames(j),"}")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"{")-1)
                read (labelIndex,*) indexElement
             end if
             if (j == parameterCount) then
                ! This is the final parameter - so get and store a pointer to its node.
                self%modelParametersInactive_(i)%parameter_   => parameters_   %node         (char(parameterNames(j)),requireValue=.true. ,copyInstance=instance)
                self%modelParametersInactive_(i)%indexElement =  indexElement
             else
                ! This is an intermediate parameter, get the appropriate sub-parameters.
                allocate  (subParameters_)
                subParameters_                                =  parameters_   %subParameters(char(parameterNames(j)),requireValue=.false.,copyInstance=instance)
                deallocate(   parameters_)
                allocate  (   parameters_)
                parameters_                                   =  inputParameters(subParameters_)
                deallocate(subParameters_)
             end if
          end do
          deallocate(parameters_   )
          deallocate(parameterNames)
       end do
    end if
    return
  end subroutine baseParametersInitialize
  
  subroutine baseParametersUpdate(self,simulationState,modelParametersActive_,modelParametersInactive_,stateVector,report)
    use :: Display           , only : displayIndent         , displayMessage    , displayUnindent              , displayVerbositySet, &
         &                            verbosityLevelStandard, displayVerbosity  , enumerationVerbosityLevelType
    use :: Error             , only : Error_Report          , errorStatusSuccess
    use :: Model_Parameters  , only : modelParameterDerived
    use :: ISO_Varying_String, only : operator(//)
    use :: Kind_Numbers      , only : kind_int8
    use :: String_Handling   , only : String_Count_Words    , String_Join       , String_Split_Words, operator(//)
    implicit none
    class           (posteriorSampleLikelihoodBaseParameters), intent(inout)               :: self
    class           (posteriorSampleStateClass              ), intent(inout)               :: simulationState
    type            (modelParameterList                     ), intent(in   ), dimension(:) :: modelParametersActive_, modelParametersInactive_
    double precision                                         , intent(in   ), dimension(:) :: stateVector
    logical                                                  , intent(in   ), optional     :: report
    type            (varying_string                         ), allocatable  , dimension(:) :: parameterNames
    integer                                                                                :: i                     , j                       , &
         &                                                                                    parameterCount
    logical                                                                                :: firstIteration        , dependenciesResolved    , &
         &                                                                                    dependenciesUpdated
    character       (len=24                                 )                              :: valueText
    double precision                                                                       :: valueDerived
    type            (varying_string                         )                              :: parameterText
    type            (enumerationVerbosityLevelType          )                              :: verbosityLevel
    ! Declarations of GNU libmatheval procedures used.
#ifdef MATHEVALAVAIL
    integer         (kind_int8                              )                              :: evaluator
#endif
    integer         (kind_int8                              ), external                    :: Evaluator_Create_
    double precision                                         , external                    :: Evaluator_Evaluate_
    external                                                                               :: Evaluator_Destroy_
    !![
    <optionalArgument name="report" defaultsTo=".false."/>
    !!]
    
    ! Set verbosity for our reports.
    if (report_ .and. (self%reportFileName .or. self%reportState)) then
       verbosityLevel=displayVerbosity()
       call displayVerbositySet(verbosityLevelStandard)
    end if
    ! Report name of parameter file being evaluated if requested.
    if (report_ .and. self%reportFileName) call displayMessage("Evaluating likelihood for file `"//char(self%baseParametersFileName)//"`")
    ! Update parameter values.
    if (report_ .and. self%reportState) call displayIndent("State:")
    do i=1,size(modelParametersActive_)
       if (report_ .and. self%reportState) &
            & call displayMessage(char(modelParametersActive_(i)%modelParameter_%name())//" = "//char(self%modelParametersActive_(i)%parameter_%get()))
       if (self%modelParametersActive_(i)%indexElement == 0) then
          ! Simply overwrite the parameter.
          call self%modelParametersActive_(i)%parameter_%set(modelParametersActive_(i)%modelParameter_%unmap(stateVector(i)))
       else
          ! Overwrite only the indexed parameter in the list.
          parameterText =self%modelParametersActive_(i)%parameter_%get()
          parameterCount=String_Count_Words(char(parameterText))
          if (self%modelParametersActive_(i)%indexElement > parameterCount)                           &
               & call Error_Report(                                                                   &
               &                   var_str('attempt to access non-existant element {')             // &
               &                          (self%modelParametersActive_(i)%indexElement          -1)// &
               &                           '} of parameter "'                                      // &
               &                   char   (     modelParametersActive_(i)%modelParameter_%name()  )// &
               &                           '"'                                                     // &
               &                   {introspection:location}                                           &
               &                  )
          allocate(parameterNames(parameterCount))
          call String_Split_Words(parameterNames,char(parameterText))
          write (valueText,'(e24.16)') modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
          parameterNames(self%modelParametersActive_(i)%indexElement)=trim(valueText)
          call self%modelParametersActive_(i)%parameter_%set(String_Join(parameterNames," "))
          deallocate(parameterNames)
       end if
    end do
    ! Resolve dependencies in derived parameters.
    if (size(modelParametersInactive_) > 0) then
       do i=1,size(modelParametersInactive_)
          select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
          class is (modelParameterDerived)
             self%modelParametersInactive_(i)%definition=modelParameter_%definition()
             self%modelParametersInactive_(i)%resolved  =.false.
          end select
       end do
       firstIteration      =.true.
       dependenciesResolved=.false.
       do while (.not.dependenciesResolved)
          dependenciesResolved=.true.
          dependenciesUpdated =.false.
          do i=1,size(modelParametersInactive_)
             select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
                class is (modelParameterDerived)
                if (index(self%modelParametersInactive_(i)%definition,"%[") /= 0) then
                   ! The expression contains dependencies on other variables. Substitute the actual values where possible.
                   !! For active parameters we only need to consider substitution on the first iteration (since they are fully defined immediately).
                   !! Also handle special parameters here:
                   !!  * %[posteriorSimulationStep] - this is the current step number in the simulation.
                   if (firstIteration) then
                      if (index(self%modelParametersInactive_(i)%definition,"%[posteriorSimulationStep]") /= 0) then
                         write (valueText,'(i10)') simulationState%count()
                         self%modelParametersInactive_(i)%definition=replace(                                                                     &
                              &                                                    self% modelParametersInactive_(i)%definition                 , &
                              &                                                    "%[posteriorSimulationStep]"                                 , &
                              &                                                    valueText                                                    , &
                              &                                              every=.true.                                                         &
                              &                                             )
                      end if
                      do j=1,size(modelParametersActive_)
                         if (index(self%modelParametersInactive_(i)%definition,"%["//modelParametersActive_  (j)%modelParameter_%name()//"]") /= 0) dependenciesUpdated=.true.
                         self%modelParametersInactive_(i)%definition=replace(                                                                     &
                              &                                                    self% modelParametersInactive_(i)%definition                 , &
                              &                                                    "%["//modelParametersActive_  (j)%modelParameter_%name()//"]", &
                              &                                                    self% modelParametersActive_  (j)%parameter_     %get ()     , &
                              &                                              every=.true.                                                         &
                              &                                             )
                      end do
                   end if
                   !! For inactive parameters we must consider them each iteration as they become resolved.
                   do j=1,size(modelParametersInactive_)
                      if (i /= j .and. self%modelParametersInactive_(j)%resolved) then
                         if (index(self%modelParametersInactive_(i)%definition,"%["//modelParametersInactive_(j)%modelParameter_%name()//"]") /= 0) dependenciesUpdated=.true.
                         self%modelParametersInactive_(i)%definition=replace(                                                                     &
                              &                                                    self% modelParametersInactive_(i)%definition                 , &
                              &                                                    "%["//modelParametersInactive_(j)%modelParameter_%name()//"]", &
                              &                                                    self% modelParametersInactive_(j)%parameter_     %get ()     , &
                              &                                              every=.true.                                                         &
                              &                                             )
                      end if
                   end do
                end if
                if (index(self%modelParametersInactive_(i)%definition,"%[") == 0) then
                   ! No dependencies remain, the expression can be evaluated.
                   self%modelParametersInactive_(i)%resolved=.true.
                   dependenciesUpdated                      =.true.
#ifdef MATHEVALAVAIL
                   evaluator   =Evaluator_Create_(char(self%modelParametersInactive_(i)%definition))
                   valueDerived=Evaluator_Evaluate_(evaluator,0,"",0.0d0)
                   call Evaluator_Destroy_(evaluator)
#else
                   call Error_Report('derived parameters require libmatheval, but it is not installed'//{introspection:location})
#endif
                   if (self%modelParametersInactive_(i)%indexElement == 0) then
                      ! Simply overwrite the parameter.
                      call self%modelParametersInactive_(i)%parameter_%set(valueDerived)
                   else
                      ! Overwrite only the indexed parameter in the list.
                      parameterText =self%modelParametersInactive_(i)%parameter_%get()
                      parameterCount=String_Count_Words(char(parameterText))
                      allocate(parameterNames(parameterCount))
                      call String_Split_Words(parameterNames,char(parameterText))
                      write (valueText,'(e24.16)') valueDerived
                      parameterNames(self%modelParametersInactive_(i)%indexElement)=trim(valueText)
                      call self%modelParametersInactive_(i)%parameter_%set(String_Join(parameterNames," "))
                      deallocate(parameterNames)
                   end if
                else
                   dependenciesResolved=.false.
                end if
                class default
                call Error_Report('support for this parameter type is not implemented'//{introspection:location})
             end select
          end do
          if (.not.dependenciesUpdated) then
             call displayVerbositySet(verbosityLevelStandard)
             call displayIndent('unresolved parameters')
             do i=1,size(modelParametersInactive_)
                select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
                   class is (modelParameterDerived)
                   if (index(self%modelParametersInactive_(i)%definition,"%[") /= 0) call displayMessage(modelParametersInactive_(i)%modelParameter_%name()//" : "//self%modelParametersInactive_(i)%definition)
                end select
             end do
             call displayUnindent('unresolved parameters')
             call Error_Report('can not resolve parameter dependencies'//{introspection:location})
          end if
          firstIteration=.false.
       end do
       if (report_ .and. self%reportState) then
          do i=1,size(modelParametersInactive_)
             call displayMessage(char(modelParametersInactive_(i)%modelParameter_%name())//" = "//char(self%modelParametersInactive_(i)%parameter_%get()))
          end do
       end if
    end if
    if (report_ .and. self%reportState) call displayUnindent("")
    call self%parametersModel%reset()
    ! Restore verbosity level.
    if (report_ .and. (self%reportFileName .or. self%reportState)) call displayVerbositySet(verbosityLevel)
    return
  end subroutine baseParametersUpdate
