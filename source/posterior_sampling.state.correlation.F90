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
  Implementation of a posterior sampling state class which computes correlation lengths.
  !!}

  !![
  <posteriorSampleState name="posteriorSampleStateCorrelation">
   <description>
    An extension of the {\normalfont \ttfamily history} state, this class also computes and stores the correlation length in each
    parameter (which is taken to be the median correlation length over all non-outlier chains).
   </description>
  </posteriorSampleState>
  !!]
  type, extends(posteriorSampleStateHistory) :: posteriorSampleStateCorrelation
     !!{
     Implementation of a correlation posterior sampling state class.
     !!}
     private
     integer                                       :: storedStateCount               , stepComputePrevious, &
          &                                           convergedCorrelationLengthCount
     double precision, allocatable, dimension(:,:) :: states
     integer         , allocatable, dimension(:  ) :: correlationLengths
   contains
     !![
     <methods>
       <method description="Return the current correlation length in the chains." method="correlationLength" />
       <method description="Compute correlation lengths in the chains." method="correlationLengthCompute" />
       <method description="Return the number of post-convergence correlation lengths that have accrued." method="postConvergenceCorrelationCount" />
     </methods>
     !!]
     procedure :: parameterCountSet               => correlationParameterCountSet
     procedure :: update                          => correlationUpdate
     procedure :: reset                           => correlationReset
     procedure :: restore                         => correlationRestore
     procedure :: postConvergenceCorrelationCount => correlationPostConvergenceCorrelationCount
     procedure :: correlationLength               => correlationCorrelationLength
     procedure :: correlationLengthCompute        => correlationCorrelationLengthCompute
  end type posteriorSampleStateCorrelation

  interface posteriorSampleStateCorrelation
     !!{
     Constructors for the \refClass{posteriorSampleStateCorrelation} posterior sampling state class.
     !!}
     module procedure correlationConstructorParameters
     module procedure correlationConstructorInternal
  end interface posteriorSampleStateCorrelation

contains

  function correlationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateCorrelation} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleStateCorrelation)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    integer                                                 :: acceptedStateCount

    !![
    <inputParameter>
      <name>acceptedStateCount</name>
      <description>The number of states to use in acceptance rate statistics.</description>
      <defaultValue>100</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleStateCorrelation(acceptedStateCount)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function correlationConstructorParameters

  function correlationConstructorInternal(acceptedStateCount) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateCorrelation} posterior sampling state class which builds the object from a
    parameter set.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    type   (posteriorSampleStateCorrelation)                :: self
    integer                                 , intent(in   ) :: acceptedStateCount
    !![
    <constructorAssign variables="acceptedStateCount"/>
    !!]

    allocate(self%accepted(acceptedStateCount))
    self%stepCount                      = 0
    self%storedStateCount               = 0
    self%accepted                       = 0
    self%stepComputePrevious            = 0
    self%convergedCorrelationLengthCount= 0
    self%chainIndexValue                =mpiSelf%rank()
    return
  end function correlationConstructorInternal

  subroutine correlationParameterCountSet(self,parameterCount)
    !!{
    Set the number of parameters in this state.
    !!}
    implicit none
    class  (posteriorSampleStateCorrelation), intent(inout) :: self
    integer                                 , intent(in   ) :: parameterCount
    integer                                 , parameter     :: correlationLengthInitial=100

    call self%posteriorSampleStateSimple%parameterCountSet(parameterCount)
    allocate(self%stateSum          (parameterCount                         ))
    allocate(self%stateSquaredSum   (parameterCount                         ))
    allocate(self%states            (parameterCount,correlationLengthInitial))
    allocate(self%correlationLengths(parameterCount                         ))
    self%stateSum          = 0.0d0
    self%stateSquaredSum   = 0.0d0
    self%correlationLengths=-1
    return
  end subroutine correlationParameterCountSet

  subroutine correlationUpdate(self,stateNew,logState,isConverged,outlierMask)
    !!{
    Update the current state.
    !!}
    implicit none
    class           (posteriorSampleStateCorrelation), intent(inout)                           :: self
    double precision                                 , intent(in   ), dimension(:  )           :: stateNew
    logical                                          , intent(in   )                           :: logState
    logical                                          , intent(in   )                           :: isConverged
    logical                                          , intent(in   ), dimension(:  ), optional :: outlierMask
    double precision                                 , allocatable  , dimension(:,:)           :: statesTmp
    integer                                                                                    :: storedStateCountMaximum
    integer                                                                                    :: stepComputeInterval

    call self%posteriorSampleStateHistory%update(stateNew,logState,isConverged,outlierMask)
    ! Return if not logging states.
    if (.not.logState) return
    ! Determine if we should compute the correlation length.
    storedStateCountMaximum=size(self%states,dim=2)
    if (any(self%correlationLengths == -1)) then
       stepComputeInterval=storedStateCountMaximum
    else
       stepComputeInterval=maxval(self%correlationLengths)
    end if
    ! Store full state.
    self%storedStateCount=self%storedStateCount+1
    if (self%storedStateCount > storedStateCountMaximum) then
       self%storedStateCount=storedStateCountMaximum
       self%states(:,1:storedStateCountMaximum-1)=self%states(:,2:storedStateCountMaximum)
    end if
    self%states(:,self%storedStateCount)=stateNew
    if (self%stepCount >= self%stepComputePrevious+stepComputeInterval) then
       ! Update count of number of correlation lengths accrued post-convergence.
       if (isConverged.and..not.any(self%correlationLengths == -1)) &
            & self%convergedCorrelationLengthCount=self%convergedCorrelationLengthCount+1
       call self%correlationLengthCompute(outlierMask)
       self%stepComputePrevious=self%stepCount
    end if
    ! If any correlation length is not found, extend the number of states stored.
    if (any(self%correlationLengths == -1) .and. self%storedStateCount == storedStateCountMaximum) then
       call Move_Alloc(self%states,statesTmp)
       allocate(self%states(size(stateNew),2*size(statesTmp,dim=2)))
       self%states(:,1:size(statesTmp,dim=2))=statesTmp
       deallocate(statesTmp)
    end if
    return
  end subroutine correlationUpdate

  subroutine correlationReset(self)
    !!{
    Reset the state object.
    !!}
    implicit none
    class(posteriorSampleStateCorrelation), intent(inout) :: self

    call self%posteriorSampleStateHistory%reset()
    self%states                         =0.0d0
    self%storedStateCount               =0
    self%convergedCorrelationLengthCount=0
    return
  end subroutine correlationReset

  integer function correlationPostConvergenceCorrelationCount(self)
    !!{
    Return the number of post-convergence correlation lengths that have accrued.
    !!}
    implicit none
    class(posteriorSampleStateCorrelation), intent(inout) :: self

    correlationPostConvergenceCorrelationCount=self%convergedCorrelationLengthCount
    return
  end function correlationPostConvergenceCorrelationCount

  integer function correlationCorrelationLength(self)
    !!{
    Return the correlation length.
    !!}
    implicit none
    class(posteriorSampleStateCorrelation), intent(inout) :: self

    if (any(self%correlationLengths <= 0)) call self%correlationLengthCompute()
    correlationCorrelationLength=maxval(self%correlationLengths)
    return
  end function correlationCorrelationLength

  subroutine correlationCorrelationLengthCompute(self,outlierMask)
    !!{
    Compute correlation lengths.
    !!}
    use :: Display           , only : displayIndent, displayMessage, displayUnindent
    use :: ISO_Varying_String, only : assignment(=), operator(//)  , varying_string
    use :: MPI_Utilities     , only : mpiSelf
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (posteriorSampleStateCorrelation), intent(inout)                            :: self
    logical                                          , intent(in   ), dimension(:), optional    :: outlierMask
    double precision                                                                            :: correlation
    double precision                                                , dimension(:), allocatable :: stateMean
    integer                                                                                     :: i          , n
    type            (varying_string                 )                                           :: message

    ! Compute correlation lengths.
    if (mpiSelf%isMaster()) call displayIndent("Computing correlation lengths")
    allocate(stateMean(size(self%states,dim=1)))
    stateMean              =sum(self%states(:,1:self%storedStateCount),dim=2)/dble(self%storedStateCount)
    self%correlationLengths=-1
    do i=1,size(stateMean)
       do n=1,self%storedStateCount-1
          correlation=sum((self%states(i,1:self%storedStateCount-n)-stateMean(i))*(self%states(i,1+n:self%storedStateCount)-stateMean(i)))
          if (correlation <= 0.0d0) then
             self%correlationLengths(i)=n
             exit
          end if
       end do
    end do
    self%correlationLengths=mpiSelf%median(self%correlationLengths,outlierMask)
    if (mpiSelf%isMaster()) then
       message="Correlation length at step "
       message=message//self%stepCount//" = "
       do i=1,size(stateMean)
          message=message//self%correlationLengths(i)
          if (i < size(stateMean)) message=message//", "
       end do
       call displayMessage (message)
       call displayUnindent("done" )
    end if
    return
  end subroutine correlationCorrelationLengthCompute

  subroutine correlationRestore(self,stateVector,first)
    !!{
    Restore the state object from file.
    !!}
    implicit none
    class           (posteriorSampleStateCorrelation), intent(inout)               :: self
    double precision                                 , intent(in   ), dimension(:) :: stateVector
    logical                                          , intent(in   )               :: first
    integer                                                                        :: storedStateCountMaximum

    ! On first restore state, reset.
    if (first) call self%reset()
    ! Perform restore of parent class.
    call self%posteriorSampleStateHistory%restore(stateVector,first)
    ! Store state.
    self%storedStateCount=self%storedStateCount+1
    storedStateCountMaximum=size(self%states,dim=2)
    if (self%storedStateCount > storedStateCountMaximum) then
       self%storedStateCount=storedStateCountMaximum
       self%states(:,1:storedStateCountMaximum-1)=self%states(:,2:storedStateCountMaximum)
    end if
    self%states(:,self%storedStateCount)=stateVector
    return
  end subroutine correlationRestore
