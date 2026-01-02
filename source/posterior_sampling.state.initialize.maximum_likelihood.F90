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
  Implementation of a posterior sampling state initializer class which initializes to the maximum likelihood state of a previous simulation.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializeMaximumLikelihood">
   <description>
    This class initializes state to the maximum likelihood state from a previous simulation---the state will be the same for all
    chains. The {\normalfont \ttfamily [logFileRoot]} parameter is used to specify the log-file root name used in the previous
    simulation.
  </description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializeMaximumLikelihood
     !!{
     Implementation of a posterior sampling state initialization class that initializes state to the maximum likelihood state
     from a previous simulation---the state will be the same for all chains.
     !!}
     private
     type   (varying_string) :: logFileRoot
   contains
     procedure :: initialize  => maximumLikelihoodInitialize
  end type posteriorSampleStateInitializeMaximumLikelihood

  interface posteriorSampleStateInitializeMaximumLikelihood
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializeMaximumLikelihood} posterior sampling state initialization class.
     !!}
     module procedure maximumLikelihoodConstructorParameters
     module procedure maximumLikelihoodConstructorInternal
  end interface posteriorSampleStateInitializeMaximumLikelihood

contains

  function maximumLikelihoodConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeMaximumLikelihood} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(posteriorSampleStateInitializeMaximumLikelihood)                :: self
    type(inputParameters                                ), intent(inout) :: parameters
    type(varying_string                                 )                :: logFileRoot

    !![
    <inputParameter>
      <name>logFileRoot</name>
      <description>The root file name of the state files from which to maximum likelihood state.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleStateInitializeMaximumLikelihood(logFileRoot)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function maximumLikelihoodConstructorParameters

  function maximumLikelihoodConstructorInternal(logFileRoot) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeMaximumLikelihood} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSampleStateInitializeMaximumLikelihood)                :: self
    type(varying_string                                 ), intent(in   ) :: logFileRoot
    !![
    <constructorAssign variables="logFileRoot"/>
    !!]

    return
  end function maximumLikelihoodConstructorInternal

  subroutine maximumLikelihoodInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by drawing at random from the parameter priors.
    !!}
    use :: Display                     , only : displayMessage
    use :: Error                       , only : Error_Report
    use :: MPI_Utilities               , only : mpiSelf
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateClass
    use :: String_Handling             , only : operator(//)
    implicit none
    class           (posteriorSampleStateInitializeMaximumLikelihood), intent(inout)               :: self
    class           (posteriorSampleStateClass                      ), intent(inout)               :: simulationState
    class           (posteriorSampleLikelihoodClass                 ), intent(inout)               :: modelLikelihood
    type            (modelParameterList                             ), intent(inout), dimension(:) :: modelParameters_
    double precision                                                 , intent(  out)               :: timeEvaluatePrevious, logLikelihood    , &
         &                                                                                            logPosterior
    double precision                                                 , allocatable  , dimension(:) :: stateVector         , stateVectorMapped
    type            (varying_string                                 )                              :: logFileName         , message
    integer                                                                                        :: stateCount          , mpiRank          , &
         &                                                                                            logFileUnit         , ioStatus         , &
         &                                                                                            i                   , j
    double precision                                                                               :: logPosterior_       , logLikelihood_
    logical                                                                                        :: converged
    character       (len=12                                         )                              :: labelValue          , labelMinimum     , &
         &                                                                                            labelMaximum

    ! Assume we have no information about the likelihood of the state by default.
    logLikelihood=logImpossible
    logPosterior =logImpossible
     ! Allocate the state vector.
    allocate(stateVector      (simulationState%dimension()))
    allocate(stateVectorMapped(simulationState%dimension()))
    ! Read state from the log files.
    do i=0,mpiSelf%count()-1
       logFileName=self%logFileRoot//'_'//i//'.log'
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
       ioStatus=0
       do while (ioStatus == 0)
          read (logFileUnit,*,iostat=ioStatus) stateCount          , &
               &                               mpiRank             , &
               &                               timeEvaluatePrevious, &
               &                               converged           , &
               &                               logPosterior_       , &
               &                               logLikelihood_      , &
               &                               stateVector
          if (ioStatus == 0) then
             if (logLikelihood_ > logLikelihood) then
                logLikelihood=logLikelihood_
                logPosterior =logPosterior_
                ! Map the state.
                do j=1,size(stateVector)
                   stateVectorMapped(j)=modelParameters_(j)%modelParameter_%map(stateVector(j))
                end do
             end if
          end if
       end do
       close(logFileUnit)
    end do
    ! Set the simulation state.
    call simulationState%update(stateVectorMapped,.false.,.false.)
    ! Check for out of range state.
    do i=1,simulationState%dimension()
       if (modelParameters_(i)%modelParameter_%logPrior(stateVector(i)) <= logImpossible) then
          write (labelValue  ,'(e12.6)') modelParameters_(i)%modelParameter_%unmap       (stateVectorMapped(i))
          write (labelMinimum,'(e12.6)') modelParameters_(i)%modelParameter_%priorMinimum(                    )
          write (labelMaximum,'(e12.6)') modelParameters_(i)%modelParameter_%priorMaximum(                    )
          message='Out of range state for parameter '
          message=message                                                                           // &
               &                           modelParameters_(i)%modelParameter_%name()     //char(10)// &
               &  ' -> on chain '        //mpiSelf            %rankLabel           ()     //char(10)// &
               &  ' -> from state file "'//logFileName                               //'"'//char(10)// &
               &  ' -> value   = '       //trim(labelValue  )                             //char(10)// &
               &  ' -> minimum = '       //trim(labelMinimum)                             //char(10)// &
               &  ' -> maximum = '       //trim(labelMaximum)
          call displayMessage(message)
       end if
    end do
    ! Deallocate the state vector.
    deallocate(stateVector      )
    deallocate(stateVectorMapped)
    return
  end subroutine maximumLikelihoodInitialize
