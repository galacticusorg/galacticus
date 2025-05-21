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
  Implementation of a posterior sampling state initializer class which initializes all chains to values read from a parameter file.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializeParameterFile">
    <description>
      This class initializes all chains to values read from a parameter file.
    </description>
    <runTimeFileDependencies paths="fileName"/>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializeParameterFile
     !!{
     Implementation of a posterior sampling state initialization class that initializes state to values read from a parameter
     file.
     !!}
     private
     type(varying_string) :: fileName
   contains
     procedure :: initialize  => parameterFileInitialize
  end type posteriorSampleStateInitializeParameterFile

  interface posteriorSampleStateInitializeParameterFile
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializeParameterFile} posterior sampling state initialization class.
     !!}
     module procedure parameterFileConstructorParameters
     module procedure parameterFileConstructorInternal
  end interface posteriorSampleStateInitializeParameterFile

contains

  function parameterFileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeParameterFile} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(posteriorSampleStateInitializeParameterFile)                :: self
    type(inputParameters                            ), intent(inout) :: parameters
    type(varying_string                             )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the parameter file from which to read initial state.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleStateInitializeParameterFile(fileName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function parameterFileConstructorParameters

  function parameterFileConstructorInternal(fileName) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeParameterFile} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSampleStateInitializeParameterFile)                :: self
    type(varying_string                             ), intent(in   ) :: fileName
    !![
    <constructorAssign variables="fileName"/>
    !!]

    return
  end function parameterFileConstructorInternal

  subroutine parameterFileInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by reading parameter values from a parameter file.
    !!}
    use :: Display                     , only : displayMessage
    use :: MPI_Utilities               , only : mpiSelf
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateClass
    use :: String_Handling             , only : String_Count_Words       , String_Split_Words, operator(//)
    implicit none
    class           (posteriorSampleStateInitializeParameterFile), intent(inout)               :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    class           (posteriorSampleLikelihoodClass             ), intent(inout)               :: modelLikelihood
    type            (modelParameterList                         ), intent(inout), dimension(:) :: modelParameters_
    double precision                                             , intent(  out)               :: timeEvaluatePrevious, logLikelihood    , &
         &                                                                                        logPosterior
    type            (varying_string                             ), allocatable  , dimension(:) :: parameterNames
    double precision                                             , allocatable  , dimension(:) :: stateVector         , stateVectorMapped, &
         &                                                                                        parameterValues
    type            (inputParameters                            )                              :: parameters          , parameters_      , &
         &                                                                                        subParameters_
    integer                                                                                    :: i                   , j                , &
         &                                                                                        instance            , indexElement     , &
         &                                                                                        parameterCount
    character       (len=  10                                   )                              :: labelIndex
    character       (len=12                                     )                              :: labelValue          , labelMinimum     , &
         &                                                                                        labelMaximum
    character       (len=1024                                   )                              :: parameterContent
    type            (varying_string                             )                              :: message

    ! We have no information about the likelihood of the state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    ! Allocate the state vector.
    allocate(stateVector      (simulationState%dimension()))
    allocate(stateVectorMapped(simulationState%dimension()))
    ! Parse the parameter file.
    parameters=inputParameters(self%fileName,noOutput=.true.)
    ! Read value for each parameter.
    do i=1,size(modelParameters_)
       parameterCount=String_Count_Words(char(modelParameters_(i)%modelParameter_%name()),"::")
       allocate(parameterNames(parameterCount))
       call String_Split_Words(parameterNames,char(modelParameters_(i)%modelParameter_%name()),"::")
       parameters_=inputParameters(parameters)
       do j=1,parameterCount
          instance    =1
          indexElement=0
          if (index(parameterNames(j),"[") /= 0) then
             labelIndex       =extract(parameterNames(j),index(parameterNames(j),"[")+1,index(parameterNames(j),"]")-1)
             parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"[")-1)
             read (labelIndex,*) instance
             instance=instance+1
          else if (index(parameterNames(j),"{") /= 0) then
             labelIndex       =extract(parameterNames(j),index(parameterNames(j),"{")+1,index(parameterNames(j),"}")-1)
             parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"{")-1)
             read (labelIndex,*) indexElement
             indexElement=indexElement+1
          end if
          if (j == parameterCount) then
             ! This is the final parameter - so extract its value.
             call parameters_%value(char(parameterNames(j)),parameterContent,copyInstance=instance)
             if (indexElement == 0) then
                read (parameterContent,*) stateVector(i)
             else
                allocate(parameterValues(indexElement))
                read (parameterContent,*) parameterValues
                stateVector(i)=parameterValues(indexElement)
                deallocate(parameterValues)
             end if
             stateVectorMapped(i)=modelParameters_(i)%modelParameter_%map(stateVector(i))
          else
             ! This is an intermediate parameter, get the appropriate sub-parameters.
             if (parameters_%isPresent(char(parameterNames(j)),requireValue=.false.)) then
                subParameters_=parameters_   %subParameters(char(parameterNames(j)),requireValue=.false.,copyInstance=instance)
                parameters_   =inputParameters(subParameters_)
             else
                call Error_Report('parameter "'//char(parameterNames(j))//'" not found in "'//char(modelParameters_(i)%modelParameter_%name())//'"'//{introspection:location})
             end if
          end if
       end do
       deallocate(parameterNames)
    end do
    call parameters%destroy()
    ! Set the simulation state.
    call simulationState%update(stateVectorMapped,.false.,.false.)
    ! Check for out of range state.
    do i=1,simulationState%dimension()
       if (modelParameters_(i)%modelParameter_%logPrior(stateVector(i)) <= logImpossible) then
          write (labelValue  ,'(e12.6)') modelParameters_(i)%modelParameter_%unmap       (stateVectorMapped(i))
          write (labelMinimum,'(e12.6)') modelParameters_(i)%modelParameter_%priorMinimum(                    )
          write (labelMaximum,'(e12.6)') modelParameters_(i)%modelParameter_%priorMaximum(                    )
          message='Out of range state for parameter '
          message=message                                                                     // &
               &                     modelParameters_(i)%modelParameter_%name()     //char(10)// &
               &  ' -> on chain '  //mpiSelf            %rankLabel           ()     //char(10)// &
               &  ' -> from file "'//self               %fileName              //'"'//char(10)// &
               &  ' -> value   = ' //trim(labelValue  )                             //char(10)// &
               &  ' -> minimum = ' //trim(labelMinimum)                             //char(10)// &
               &  ' -> maximum = ' //trim(labelMaximum)
          call displayMessage(message)
       end if
    end do
    return
  end subroutine parameterFileInitialize
