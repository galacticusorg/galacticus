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
  Implementation of a posterior sampling state class which stores history.
  !!}

  !![
  <posteriorSampleState name="posteriorSampleStateHistory">
   <description>
    An extension of the {\normalfont \ttfamily simple} state, this class also records the mean and variance of each parameter over the
    history of the simulation.
   </description>
  </posteriorSampleState>
  !!]
  type, extends(posteriorSampleStateSimple) :: posteriorSampleStateHistory
     !!{
     Implementation of a posterior sampling state class which stores history.
     !!}
     private
     double precision, allocatable, dimension(:) :: stateSum, stateSquaredSum
   contains
     procedure :: parameterCountSet => historyParameterCountSet
     procedure :: update            => historyUpdate
     procedure :: mean              => historyMean
     procedure :: variance          => historyVariance
     procedure :: reset             => historyReset
     procedure :: restore           => historyRestore
  end type posteriorSampleStateHistory

  interface posteriorSampleStateHistory
     !!{
     Constructors for the \refClass{posteriorSampleStateHistory} posterior sampling state class.
     !!}
     module procedure historyConstructorParameters
     module procedure historyConstructorInternal
  end interface posteriorSampleStateHistory

contains

  function historyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateHistory} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleStateHistory)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    integer                                             :: acceptedStateCount

    !![
    <inputParameter>
      <name>acceptedStateCount</name>
      <description>The number of states to use in acceptance rate statistics.</description>
      <defaultValue>100</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleStateHistory(acceptedStateCount)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function historyConstructorParameters

  function historyConstructorInternal(acceptedStateCount) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateHistory} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    type   (posteriorSampleStateHistory)                :: self
    integer                             , intent(in   ) :: acceptedStateCount
    !![
    <constructorAssign variables="acceptedStateCount"/>
    !!]

    allocate(self%accepted(acceptedStateCount))
    self%stepCount      =0
    self%accepted       =0
    self%chainIndexValue=mpiSelf%rank()
    return
  end function historyConstructorInternal

  subroutine historyParameterCountSet(self,parameterCount)
    !!{
    Set the number of parameters in this state.
    !!}
    implicit none
    class  (posteriorSampleStateHistory), intent(inout) :: self
    integer                             , intent(in   ) :: parameterCount

    allocate(self%stateSum       (parameterCount))
    allocate(self%stateSquaredSum(parameterCount))
    self%stateSum       =0.0d0
    self%stateSquaredSum=0.0d0
    return
  end subroutine historyParameterCountSet

  subroutine historyUpdate(self,stateNew,logState,isConverged,outlierMask)
    !!{
    Update the current state.
    !!}
    implicit none
    class           (posteriorSampleStateHistory), intent(inout)                         :: self
    double precision                             , intent(in   ), dimension(:)           :: stateNew
    logical                                      , intent(in   )                         :: logState
    logical                                      , intent(in   )                         :: isConverged
    logical                                      , intent(in   ), dimension(:), optional :: outlierMask

    call self%posteriorSampleStateSimple%update(stateNew,logState,isConverged,outlierMask)
    if (logState) then
       self%stateSum       =self%stateSum       +stateNew
       self%stateSquaredSum=self%stateSquaredSum+stateNew**2
    end if
    return
  end subroutine historyUpdate

  function historyMean(self)
    !!{
    Return the mean over state history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleStateHistory), intent(inout)                  :: self
    double precision                             , dimension(self%parameterCount) :: historyMean

    historyMean=self%stateSum/dble(self%stepCount)
    return
  end function historyMean

  function historyVariance(self)
    !!{
    Return the mean over state history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleStateHistory), intent(inout)                  :: self
    double precision                             , dimension(self%parameterCount) :: historyVariance

    historyVariance=(self%stateSquaredSum/dble(self%stepCount)-self%mean()**2)*dble(self%stepCount)/dble(self%stepCount-1)
    return
  end function historyVariance

  subroutine historyReset(self)
    !!{
    Reset the state object.
    !!}
    implicit none
    class(posteriorSampleStateHistory), intent(inout) :: self

    call self%posteriorSampleStateSimple%reset()
    self%stateSum       =0.0d0
    self%stateSquaredSum=0.0d0
    return
  end subroutine historyReset

  subroutine historyRestore(self,stateVector,first)
    !!{
    Restore the state object from file.
    !!}
    implicit none
    class           (posteriorSampleStateHistory), intent(inout)               :: self
    double precision                             , intent(in   ), dimension(:) :: stateVector
    logical                                      , intent(in   )               :: first

    ! On first restore state, reset.
    if (first) call self%reset()
    ! Perform restore of parent class.
    call self%posteriorSampleStateSimple%restore(stateVector,first)
    ! Accumulate state and state-squared.
    self%stateSum       =self%stateSum       +stateVector
    self%stateSquaredSum=self%stateSquaredSum+stateVector**2
    return
  end subroutine historyRestore
