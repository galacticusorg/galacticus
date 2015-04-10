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

!% Defines methods for the ``adaptive'' differential evolution proposal size class.

function deProposalSizeAdaptiveConstructor(gamma,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount)
  !% Constructor for the ``adaptive'' differential evolution proposal size class.
  implicit none
  type            (deProposalSizeAdaptive)                :: deProposalSizeAdaptiveConstructor
  double precision                        , intent(in   ) :: gamma                , gammaAdjustFactor    , &
       &                                                     gammaMinimum         , gammaMaximum         , &
       &                                                     acceptanceRateMinimum, acceptanceRateMaximum
  integer                                 , intent(in   ) :: updateCount

  deProposalSizeAdaptiveConstructor%gammaCurrent         =gamma
  deProposalSizeAdaptiveConstructor%gammaAdjustFactor    =gammaAdjustFactor
  deProposalSizeAdaptiveConstructor%gammaMinimum         =gammaMinimum
  deProposalSizeAdaptiveConstructor%gammaMaximum         =gammaMaximum
  deProposalSizeAdaptiveConstructor%acceptanceRateMinimum=acceptanceRateMinimum
  deProposalSizeAdaptiveConstructor%acceptanceRateMaximum=acceptanceRateMaximum
  deProposalSizeAdaptiveConstructor%updateCount          =updateCount
  deProposalSizeAdaptiveConstructor%lastUpdateCount      =0
  return
end function deProposalSizeAdaptiveConstructor

double precision function deProposalSizeAdaptiveGamma(self,simulationState,simulationConvergence)
  !% Return the current state.
  use MPI_Utilities
  use ISO_Varying_String
  use String_Handling
  use Galacticus_Display
  implicit none
  class           (deProposalSizeAdaptive), intent(inout) :: self
  class           (state                 ), intent(in   ) :: simulationState
  class           (convergence           ), intent(inout) :: simulationConvergence
  double precision                        , dimension(1)  :: acceptanceRate
  character       (len=8                 )                :: label
  type            (varying_string        )                :: message
 
  ! Should we consider updating gamma?
  if     (                                                                                   &
       &        simulationState      %count      () >= self%lastUpdateCount+self%updateCount &
       &  .and.                                                                              &
       &   .not.simulationConvergence%isConverged()                                          &
       & ) then
     ! Reset the number of steps remaining.
     self%lastUpdateCount=simulationState%count()
     ! Find the mean acceptance rate across all chains.
     acceptanceRate=mpiSelf%average([simulationState%acceptanceRate()])
     if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
        write (label,'(f5.3)') acceptanceRate(1)
        message='After '
        message=message//simulationState%count()//' steps, acceptance rate is '//trim(label)
        call Galacticus_Display_Message(message)
     end if
     ! If the acceptance rate is out of range, adjust gamma.
     if      (acceptanceRate(1) > self%acceptanceRateMaximum .and. self%gammaCurrent < self%gammaMaximum) then
        self%gammaCurrent=min(self%gammaCurrent*self%gammaAdjustFactor,self%gammaMaximum)
        if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
           write (label,'(f8.5)') self%gammaCurrent
           call Galacticus_Display_Message('Adjusting gamma up to '//label)
        end if
     else if (acceptanceRate(1) < self%acceptanceRateMinimum .and. self%gammaCurrent > self%gammaMinimum) then
        self%gammaCurrent=max(self%gammaCurrent/self%gammaAdjustFactor,self%gammaMinimum)
        if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
           write (label,'(f8.5)') self%gammaCurrent
           call Galacticus_Display_Message('Adjusting gamma down to '//label)
        end if
     end if
  end if
  ! Return the current adaptive size. 
  deProposalSizeAdaptiveGamma=self%gammaCurrent
  return
end function deProposalSizeAdaptiveGamma
