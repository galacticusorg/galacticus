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

!% Defines methods for the ``simple'' state class.

function stateSimpleConstructor(parameterCount,acceptedStateCount)
  !% Constructor for ``simple'' state class.
  implicit none
  type   (stateSimple)                :: stateSimpleConstructor
  integer             , intent(in   ) :: parameterCount        , acceptedStateCount

  allocate(stateSimpleConstructor%current (    parameterCount))
  allocate(stateSimpleConstructor%accepted(acceptedStateCount))
  stateSimpleConstructor%parameterCount    =parameterCount
  stateSimpleConstructor%stepCount         =0
  stateSimpleConstructor%accepted          =0
  return
end function stateSimpleConstructor

subroutine stateSimpleDestroy(self)
  !% Destructor for ``simple'' state class.
  implicit none
  type(stateSimple), intent(inout) :: self

  if (allocated(self%current)) deallocate(self%current)
  return
end subroutine stateSimpleDestroy

function stateSimpleGet(self)
  !% Return the current state.
  implicit none
  class           (stateSimple), intent(in   )                  :: self
  double precision             , dimension(self%parameterCount) :: stateSimpleGet
  
  stateSimpleGet=self%current
  return
end function stateSimpleGet

subroutine stateSimpleUpdate(self,stateNew,logState,isConverged,outlierMask)
  !% Update the current state.
  implicit none
  class           (stateSimple), intent(inout)                         :: self
  double precision             , intent(in   ), dimension(:)           :: stateNew
  logical                      , intent(in   )                         :: logState
  logical                      , intent(in   )                         :: isConverged
  logical                      , intent(in   ), dimension(:), optional :: outlierMask
  integer                                                              :: i

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
end subroutine stateSimpleUpdate

function stateSimpleMean(self)
  !% Return the mean over state history.
  use Galacticus_Error
  implicit none
  class           (stateSimple), intent(in   )                  :: self
  double precision             , dimension(self%parameterCount) :: stateSimpleMean
  
  call Galacticus_Error_Report('stateSimpleMean','the "simple" state class does not store history')
  return
end function stateSimpleMean

function stateSimpleVariance(self)
  !% Return the mean over state history.
  use Galacticus_Error
  implicit none
  class           (stateSimple), intent(in   )                  :: self
  double precision             , dimension(self%parameterCount) :: stateSimpleVariance

  call Galacticus_Error_Report('stateSimpleVariance','the "simple" state class does not store history')
  return
end function stateSimpleVariance

double precision function stateSimpleAcceptanceRate(self)
  !% Return the acceptance rate for the state.
  implicit none
  class  (stateSimple), intent(in   ) :: self
  integer                            :: stepCount

  stepCount=min(self%stepCount,size(self%accepted))
  if (stepCount > 0) then
     stateSimpleAcceptanceRate=dble(sum(self%accepted))/dble(stepCount)
  else
     stateSimpleAcceptanceRate=0.0d0
  end if
  return
end function stateSimpleAcceptanceRate

subroutine stateSimpleReset(self)
  !% Reset the state object.
  implicit none
  class(stateSimple), intent(inout) :: self
  
  self%stepCount=0
  self%accepted =0
  return
end subroutine stateSimpleReset
