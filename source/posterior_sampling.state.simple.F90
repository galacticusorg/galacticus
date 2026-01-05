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
  Implementation of a simple posterior sampling state class.
  !!}

  !![
  <posteriorSampleState name="posteriorSampleStateSimple">
   <description>
   This class stores the current state but makes no attempt to record a history of the state and so cannot provide measures of the
   mean or variance of state over the simulation history. It does, however, maintain a running average of the state acceptance
   rate. The number of steps over which the acceptance rate should be computed is specified by the {\normalfont \ttfamily
   acceptedStateCount}.
   </description>
  </posteriorSampleState>
  !!]
  type, extends(posteriorSampleStateClass) :: posteriorSampleStateSimple
     !!{
     Implementation of a simple posterior sampling state class.
     !!}
     private
     double precision, allocatable, dimension(:) :: current
     integer         , allocatable, dimension(:) :: accepted
     integer                                     :: acceptedStateCount
   contains
     !![
     <methods>
       <method description="Set the state count." method="countSet" />
     </methods>
     !!]
     procedure :: parameterCountSet => simpleParameterCountSet
     procedure :: get               => simpleGet
     procedure :: update            => simpleUpdate
     procedure :: mean              => simpleMean
     procedure :: variance          => simpleVariance
     procedure :: acceptanceRate    => simpleAcceptanceRate
     procedure :: reset             => simpleReset
     procedure :: restore           => simpleRestore
     procedure :: countSet          => simpleCountSet
  end type posteriorSampleStateSimple

  interface posteriorSampleStateSimple
     !!{
     Constructors for the \refClass{posteriorSampleStateSimple} posterior sampling state class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface posteriorSampleStateSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateSimple} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleStateSimple)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    integer                                            :: acceptedStateCount

    !![
    <inputParameter>
      <name>acceptedStateCount</name>
      <description>The number of states to use in acceptance rate statistics.</description>
      <defaultValue>100</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleStateSimple(acceptedStateCount)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function simpleConstructorParameters

  function simpleConstructorInternal(acceptedStateCount) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateSimple} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    type   (posteriorSampleStateSimple)                :: self
    integer                            , intent(in   ) :: acceptedStateCount
    !![
    <constructorAssign variables="acceptedStateCount"/>
    !!]

    allocate(self%accepted(acceptedStateCount))
    self%stepCount      =0
    self%accepted       =0
    self%chainIndexValue=mpiSelf%rank()
    return
  end function simpleConstructorInternal

  subroutine simpleParameterCountSet(self,parameterCount)
    !!{
    Set the number of parameters in this state.
    !!}
    implicit none
    class  (posteriorSampleStateSimple), intent(inout) :: self
    integer                            , intent(in   ) :: parameterCount

    self%parameterCount=parameterCount
    allocate(self%current(parameterCount))
    self%current=-huge(0.0d0)
    return
  end subroutine simpleParameterCountSet

  function simpleGet(self)
    !!{
    Return the current state.
    !!}
    implicit none
    class           (posteriorSampleStateSimple), intent(inout)                  :: self
    double precision                            , dimension(self%parameterCount) :: simpleGet

    simpleGet=self%current
    return
  end function simpleGet

  subroutine simpleUpdate(self,stateNew,logState,isConverged,outlierMask)
    !!{
    Update the current state.
    !!}
    implicit none
    class           (posteriorSampleStateSimple), intent(inout)                         :: self
    double precision                            , intent(in   ), dimension(:)           :: stateNew
    logical                                     , intent(in   )                         :: logState
    logical                                     , intent(in   )                         :: isConverged
    logical                                     , intent(in   ), dimension(:), optional :: outlierMask
    integer                                                                             :: i
    !$GLC attributes unused :: isConverged, outlierMask

    if (logState) then
       i=mod(self%stepCount,size(self%accepted))+1
       if (any(stateNew /= self%current)) then
          self%accepted(i)=1
       else
          self%accepted(i)=0
       end if
       self%stepCount=self%stepCount+1
    end if
    self%current=stateNew
    return
  end subroutine simpleUpdate

  function simpleMean(self)
    !!{
    Return the mean over state history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleStateSimple), intent(inout)                  :: self
    double precision                            , dimension(self%parameterCount) :: simpleMean

    simpleMean=0.0d0
    call Error_Report('the "simple" state class does not store history'//{introspection:location})
    return
  end function simpleMean

  function simpleVariance(self)
    !!{
    Return the mean over state history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleStateSimple), intent(inout)                  :: self
    double precision                            , dimension(self%parameterCount) :: simpleVariance

    simpleVariance=0.0d0
    call Error_Report('the "simple" state class does not store history'//{introspection:location})
    return
  end function simpleVariance

  double precision function simpleAcceptanceRate(self)
    !!{
    Return the acceptance rate for the state.
    !!}
    implicit none
    class  (posteriorSampleStateSimple), intent(inout) :: self
    integer                                            :: stepCount

    stepCount=min(self%stepCount,size(self%accepted))
    if (stepCount > 0) then
       simpleAcceptanceRate=dble(sum(self%accepted))/dble(stepCount)
    else
       simpleAcceptanceRate=0.0d0
    end if
    return
  end function simpleAcceptanceRate

  subroutine simpleReset(self)
    !!{
    Reset the state object.
    !!}
    implicit none
    class(posteriorSampleStateSimple), intent(inout) :: self

    self%stepCount=0
    self%accepted =0
    return
  end subroutine simpleReset

  subroutine simpleRestore(self,stateVector,first)
    !!{
    Restore the state object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleStateSimple), intent(inout)               :: self
    double precision                            , intent(in   ), dimension(:) :: stateVector
    logical                                     , intent(in   )               :: first
    integer                                                                   :: i

    ! On first restore state, reset.
    if (first) call self%reset()
    ! Increment state count.
    self%stepCount=self%stepCount+1
    ! Increment the accepted count.
    i=mod(self%stepCount,size(self%accepted))+1
    if (first.or.any(stateVector /= self%current)) then
       self%accepted(i)=1
    else
       self%accepted(i)=0
    end if
    ! Update the state.
    self%current=stateVector
    return
  end subroutine simpleRestore

  subroutine simpleCountSet(self,stateCount)
    !!{
    Set the state count.
    !!}
    implicit none
    class  (posteriorSampleStateSimple), intent(inout) :: self
    integer                            , intent(in   ) :: stateCount
    
    self%stepCount=stateCount
    return
  end subroutine simpleCountSet
