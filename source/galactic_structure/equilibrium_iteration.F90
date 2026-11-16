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

!!{RST
Tracks radius iterations and applies the equilibrium solver's oscillation correction.
!!}

!+    Contributions to this file made by: Andrew Robertson, Codex.

module Galactic_Structure_Equilibrium_Iteration
  !!{RST
  The numerical update used by the equilibrium galactic structure solver. Keeping this independent of the gravitational
  calculation lets its damping and convergence checks be tested with smooth equations of known solution.
  !!}
  implicit none
  private
  public :: equilibriumIterationState

  type equilibriumIterationState
     !!{RST
     Two previous input radii per component. Reset at the start of each solve; allocation is retained for reuse.
     !!}
     private
     double precision, allocatable, dimension(:,:) :: history
   contains
     procedure :: reset  => iterationReset
     procedure :: update => iterationUpdate
  end type equilibriumIterationState

contains

  subroutine iterationReset(self)
    !!{RST
    Discard the previous solve's iteration history without changing the physical radii used as initial guesses.
    !!}
    implicit none
    class(equilibriumIterationState), intent(inout) :: self

    if (allocated(self%history)) self%history=-1.0d0
    return
  end subroutine iterationReset

  subroutine iterationUpdate(self,iteration,component,radius,radiusNew,fitContribution)
    !!{RST
    Apply the existing alternating arithmetic/geometric oscillation correction to the proposed next radius.
    The fit contribution includes the undamped residual so damping alone cannot signal convergence.
    !!}
    implicit none
    class           (equilibriumIterationState), intent(inout)                 :: self
    integer                                    , intent(in   )                 :: iteration                       , &
         &                                                                        component
    double precision                           , intent(in   )                 :: radius
    double precision                           , intent(inout)                 :: radiusNew
    double precision                           , intent(  out)                 :: fitContribution
    integer                                    , parameter                     :: iterationsForBisectionMinimum=10, &
         &                                                                        componentIncrement           =2
    integer                                                                    :: previousSize
    double precision                           , allocatable  , dimension(:,:) :: historyTemporary

    ! Allocate space for newly active components, preserving the history of existing components.
    if (.not.allocated(self%history)) then
       allocate(self%history(2,component+componentIncrement))
       self%history=-1.0d0
    else if (size(self%history,dim=2) < component) then
       previousSize=size(self%history,dim=2)
       call move_alloc(self%history,historyTemporary)
       allocate(self%history(2,component+componentIncrement))
       self%history(:,1:previousSize)=historyTemporary
       self%history(:,previousSize+1:)=-1.0d0
    end if
    ! Preserve the undamped residual before modifying the proposed step.
    fitContribution=0.0d0
    if (radius > 0.0d0 .and. radiusNew > 0.0d0) fitContribution=abs(log(radiusNew/radius))
    if (iteration > iterationsForBisectionMinimum .and. all(self%history(:,component) >= 0.0d0)) then
       if ((self%history(2,component)-self%history(1,component))*(self%history(1,component)-radius) < 0.0d0) then
          ! Alternate the same four averages as the original solver. The correction must change radiusNew, which is accepted
          ! by the caller, rather than radius, which is the location where the potential was evaluated.
          select case (mod(iteration,4))
          case (0)
             radiusNew=sqrt(radius*self%history(1,component))
          case (1)
             radiusNew=0.5d0*(radius+self%history(1,component))
          case (2)
             radiusNew=sqrt(self%history(1,component)*self%history(2,component))
          case (3)
             radiusNew=0.5d0*(self%history(1,component)+self%history(2,component))
          end select
          self%history(:,component)=-1.0d0
       end if
    end if
    ! Store the actual input radius, not an average that was never used to evaluate the gravitational potential.
    self%history(2,component)=self%history(1,component)
    self%history(1,component)=radius
    if (radius > 0.0d0 .and. radiusNew > 0.0d0) fitContribution=max(fitContribution,abs(log(radiusNew/radius)))
    return
  end subroutine iterationUpdate

end module Galactic_Structure_Equilibrium_Iteration
