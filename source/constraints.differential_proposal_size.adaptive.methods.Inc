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

function deProposalSizeAdaptiveConstructor(gamma,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount)
  !% Constructor for the ``adaptive'' differential evolution proposal size class.
  implicit none
  type            (deProposalSizeAdaptive)                :: deProposalSizeAdaptiveConstructor
  double precision                        , intent(in   ) :: gamma                , gammaAdjustFactor    , &
       &                                                     acceptanceRateMinimum, acceptanceRateMaximum
  integer                                 , intent(in   ) :: updateCount

  deProposalSizeAdaptiveConstructor%gammaCurrent         =gamma
  deProposalSizeAdaptiveConstructor%gammaAdjustFactor    =gammaAdjustFactor
  deProposalSizeAdaptiveConstructor%acceptanceRateMinimum=acceptanceRateMinimum
  deProposalSizeAdaptiveConstructor%acceptanceRateMaximum=acceptanceRateMaximum
  deProposalSizeAdaptiveConstructor%updateCount          =updateCount
  return
end function deProposalSizeAdaptiveConstructor

double precision function deProposalSizeAdaptiveGamma(self,simulationState,simulationConvergence)
  !% Return the current state.
  use MPI_Utilities
  implicit none
  class           (deProposalSizeAdaptive), intent(inout) :: self
  class           (state                 ), intent(in   ) :: simulationState
  class           (convergence           ), intent(inout) :: simulationConvergence
  double precision                        , dimension(1)  :: acceptanceRate
 
  ! Should we consider updating gamma?
  if     (                                                    &
       &       simulationState%count()                    > 0 &
       &  .and.                                               &
       &   mod(simulationState%count(),self%updateCount) == 0 &
       &  .and.                                               &
       &   .not.simulationConvergence%isConverged()           &
       & ) then
     ! Find the mean acceptance rate across all chains.
     acceptanceRate=mpiSelf%average([simulationState%acceptanceRate()])
     if (mpiSelf%rank() == 0) write (0,*) 'After ',simulationState%count(),' steps, acceptance rate is ',acceptanceRate(1)
     ! If the acceptance rate is out of range, adjust gamma.
     if      (acceptanceRate(1) > self%acceptanceRateMaximum) then
        self%gammaCurrent=self%gammaCurrent*self%gammaAdjustFactor
        if (mpiSelf%rank() == 0) write (0,*) 'Adjusting gamma up to ',self%gammaCurrent
     else if (acceptanceRate(1) < self%acceptanceRateMinimum) then
        self%gammaCurrent=self%gammaCurrent/self%gammaAdjustFactor
        if (mpiSelf%rank() == 0) write (0,*) 'Adjusting gamma down to ',self%gammaCurrent
     end if
  end if
  ! Return the current adaptive size. 
  deProposalSizeAdaptiveGamma=self%gammaCurrent
  return
end function deProposalSizeAdaptiveGamma
