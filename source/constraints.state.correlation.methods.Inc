!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for the ``correlation'' state class.

function stateCorrelationConstructor(parameterCount,acceptedStateCount)
  !% Constructor for ``correlation'' state class.
  implicit none
  type   (stateCorrelation)                :: stateCorrelationConstructor
  integer                  , intent(in   ) :: parameterCount                 , acceptedStateCount
  integer                  , parameter     :: correlationLengthInitial   =100

  allocate(stateCorrelationConstructor%current          (    parameterCount                         ))
  allocate(stateCorrelationConstructor%accepted         (acceptedStateCount                         ))
  allocate(stateCorrelationConstructor%stateSum         (    parameterCount                         ))
  allocate(stateCorrelationConstructor%stateSquaredSum  (    parameterCount                         ))
  allocate(stateCorrelationConstructor%states           (    parameterCount,correlationLengthInitial))
  allocate(stateCorrelationConstructor%correlationLength(    parameterCount                         ))
  stateCorrelationConstructor%parameterCount                 = parameterCount
  stateCorrelationConstructor%stepCount                      = 0
  stateCorrelationConstructor%storedStateCount               = 0
  stateCorrelationConstructor%stateSum                       = 0.0d0
  stateCorrelationConstructor%stateSquaredSum                = 0.0d0
  stateCorrelationConstructor%accepted                       = 0
  stateCorrelationConstructor%stepComputePrevious            = 0
  stateCorrelationConstructor%correlationLength              =-1
  stateCorrelationConstructor%convergedCorrelationLengthCount= 0
  return
end function stateCorrelationConstructor

subroutine stateCorrelationDestroy(self)
  !% Destructor for ``correlation'' state class.
  implicit none
  type(stateCorrelation), intent(inout) :: self

  if (allocated(self%states)) deallocate(self%states)
  return
end subroutine stateCorrelationDestroy

subroutine stateCorrelationUpdate(self,stateNew,logState,isConverged,outlierMask)
  !% Update the current state.
  use MPI_Utilities
  use Galacticus_Display
  use String_Handling
  use ISO_Varying_String
  implicit none
  class           (stateCorrelation), intent(inout)                                      :: self
  double precision                  , intent(in   ), dimension(:             )           :: stateNew
  logical                           , intent(in   )                                      :: logState
  logical                           , intent(in   )                                      :: isConverged
  logical                           , intent(in   ), dimension(:)             , optional :: outlierMask
  double precision                                 , dimension(size(stateNew))           :: stateMean
  double precision                  , allocatable  , dimension(:,:           )           :: statesTmp
  integer                                                                                :: storedStateCountMaximum
  integer                                                                                :: i                      , n, &
       &                                                                                    stepComputeInterval
  double precision                                                                       :: correlation
  type            (varying_string  )                                                     :: message

  call self%stateHistory%update(stateNew,logState,isConverged,outlierMask)
  ! Store full state.
  self%storedStateCount=self%storedStateCount+1
  storedStateCountMaximum=size(self%states,dim=2)
  if (self%storedStateCount > storedStateCountMaximum) then
     self%storedStateCount=storedStateCountMaximum
     self%states(:,1:storedStateCountMaximum-1)=self%states(:,2:storedStateCountMaximum)
  end if
  self%states(:,self%storedStateCount)=stateNew
  ! Determine if we should compute the correlation length.
  if (any(self%correlationLength == -1)) then
     stepComputeInterval=storedStateCountMaximum
  else
     stepComputeInterval=maxval(self%correlationLength)
  end if
  if (self%stepCount >= self%stepComputePrevious+stepComputeInterval) then
     ! Update count of number of correlation lengths accrued post-convergence.
     if (isConverged.and..not.any(self%correlationLength == -1)) &
          & self%convergedCorrelationLengthCount=self%convergedCorrelationLengthCount+1
     ! Compute correlation lengths.
     if (mpiSelf%isMaster()) call Galacticus_Display_Indent("Computing correlation lengths")
     stateMean             =sum(self%states(:,1:self%storedStateCount),dim=2)/dble(self%storedStateCount)
     self%correlationLength=-1
     do i=1,size(stateMean)
        do n=1,self%storedStateCount-1
           correlation=sum((self%states(i,1:self%storedStateCount-n)-stateMean(i))*(self%states(i,1+n:self%storedStateCount)-stateMean(i)))
           if (correlation <= 0.0d0) then
              self%correlationLength(i)=n
              exit
           end if
        end do
     end do
     self%correlationLength=mpiSelf%median(self%correlationLength,outlierMask)
     self%stepComputePrevious=self%stepCount
     if (mpiSelf%isMaster()) then
        message="Correlation length at step "
        message=message//self%stepCount//" = "
        do i=1,size(stateNew)
           message=message//self%correlationLength(i)
           if (i < size(stateNew)) message=message//", "
        end do
        call Galacticus_Display_Message (message)
        call Galacticus_Display_Unindent("done" )
     end if
  end if
  ! If any correlation length is not found, extend the number of states stored.
  if (any(self%correlationLength == -1) .and. self%storedStateCount == storedStateCountMaximum) then
     call Move_Alloc(self%states,statesTmp)
     allocate(self%states(size(stateNew),2*size(statesTmp,dim=2)))
     self%states(:,1:size(statesTmp,dim=2))=statesTmp
     deallocate(statesTmp)
  end if
  return
end subroutine stateCorrelationUpdate

subroutine stateCorrelationReset(self)
  !% Reset the state object.
  implicit none
  class(stateCorrelation), intent(inout) :: self
  
  call self%stateHistory%reset()
  self%states                         =0.0d0
  self%storedStateCount               =0
  self%convergedCorrelationLengthCount=0
  return
end subroutine stateCorrelationReset

integer function stateCorrelationPostConvergenceCorrelationCount(self)
  !% Return the number of post-convergence correlation lengths that have accrued.
  implicit none
  class(stateCorrelation), intent(inout) :: self
  
  stateCorrelationPostConvergenceCorrelationCount=self%convergedCorrelationLengthCount
  return
end function stateCorrelationPostConvergenceCorrelationCount
