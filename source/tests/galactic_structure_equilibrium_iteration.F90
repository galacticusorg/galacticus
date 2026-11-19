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
Tests the production equilibrium radius update independently of halo profiles and merger-tree evolution.
!!}

!+    Contributions to this file made by: Andrew Robertson, Codex.

program Test_Galactic_Structure_Equilibrium_Iteration
  !!{RST
  A smooth map with a known root exercises the oscillation correction. Additional cases check false convergence, history
  reset, component isolation and ordinary monotonic convergence. All corrected updates call the production implementation.
  !!}
  use :: Display                                 , only : displayVerbositySet      , verbosityLevelStandard
  use :: Galactic_Structure_Equilibrium_Iteration, only : equilibriumIterationState
  use :: Unit_Tests                              , only : Assert                   , Unit_Tests_Begin_Group, Unit_Tests_End_Group, &
  &                                                       Unit_Tests_Finish
  implicit none
  type            (equilibriumIterationState)            :: state            , fresh
  double precision                                       :: radius           , proposal     , &
       &                                                    fit              , expected     , &
       &                                                    initial          , freshProposal, &
       &                                                    freshFit
  double precision                           , parameter :: tolerance=1.0d-10
  integer                                                :: iteration        , method       , &
       &                                                    count            , repeat
  character       (len=80                   )            :: label
  logical                                                :: converged

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group('Equilibrium radius iteration')

  ! The smooth map P(R)=4/R has a unique positive fixed point R=2 and an undamped two-cycle 1 ↔ 4.
  ! It can be written in the solver's form P=√(j R/V) with j=1 and V(R)=R^3/16.
  ! This is an artificial numerical fixture, not a galaxy model or a contraction-boundary discontinuity.
  radius=1.0d0
  do iteration=2,100
     proposal=4.0d0/radius
     fit=abs(log(proposal/radius))
     radius=proposal
  end do
  call Assert('undamped two-cycle does not converge',fit > tolerance,.true.)
  do repeat=1,2
     if (repeat == 1) then
        initial=1.0d0
     else
        initial=4.0d0
     end if
     call state%reset()
     call solveCycle(state,initial,radius,count,converged)
     call Assert('damping converges within the unchanged 100-iteration limit',converged,.true.)
     call Assert('damping reaches the analytic root',radius,2.0d0,relTol=tolerance)
     call Assert('final equation residual',abs(log(radius**4/16.0d0)) < 4.0d0*tolerance,.true.)
     write (*,'(a,f6.2,a,i0,a,es12.4)') 'Smooth cycle: initial radius ',initial,', iterations ',count,', final residual ',abs(log(radius**4/16.0d0))
  end do

  ! Independently specified averages for history [1,4] and current radius 2.
  ! In particular, the geometric average of the two history entries is exactly the current radius, despite the unmet equation.
  do method=0,3
     call state%reset()
     call seedHistory(state,1)
     radius=2.0d0
     proposal=2.2d0
     select case (method)
     case (0)
        expected=sqrt(2.0d0)
     case (1)
        expected=1.5d0
     case (2)
        expected=2.0d0
     case (3)
        expected=2.5d0
     end select
     call state%update(12+method,1,radius,proposal,fit)
     write (label,'(a,i0)') 'accepted oscillation correction, method ',method
     call Assert(trim(label),proposal,expected,relTol=1.0d-14)
     call Assert('input radius remains the evaluation radius',radius,2.0d0,relTol=1.0d-14)
     call Assert('fit retains the undamped equation residual',fit >= log(1.1d0)-1.0d-14,.true.)
     call Assert('fit includes the accepted step',fit,max(log(1.1d0),abs(log(expected/2.0d0))),absTol=1.0d-14)
  end do

  ! A reset must make a previously used state behave like a fresh state. This also checks the component whose index caused
  ! the history array to grow. Using iteration 14 exposes stale history directly rather than overwriting it during burn-in.
  call state%reset()
  call seedHistory(state,1)
  call seedHistory(state,8)
  proposal=2.2d0
  call state%update(14,1,2.0d0,proposal,fit)
  call Assert('growing the history preserves earlier components',proposal,2.0d0,relTol=1.0d-14)
  call seedHistory(state,1)
  call state%reset()
  call state%reset()
  do method=1,8,7
     proposal=2.2d0
     freshProposal=proposal
     call state%update(14,method,2.0d0,proposal,fit)
     call fresh%update(14,method,2.0d0,freshProposal,freshFit)
     call Assert('reset removes history for every component',proposal,freshProposal,relTol=1.0d-14)
     call Assert('reset restores the undamped fit'          ,fit     ,freshFit     ,absTol=1.0d-14)
     call Assert('first step after reset is not damped'     ,proposal,2.2d0        ,relTol=1.0d-14)
  end do

  ! A populated neighboring component must not damp a component with no history.
  call state%reset()
  call seedHistory(state,1)
  proposal=2.2d0
  call state%update(14,2,2.0d0,proposal,fit)
  call Assert('component histories are independent',proposal,2.2d0,relTol=1.0d-14)
  call state%reset()
  call seedHistory(state,1)
  proposal=2.2d0
  call state%update(10,1,2.0d0,proposal,fit)
  call Assert('burn-in does not damp the proposal',proposal,2.2d0,relTol=1.0d-14)

  ! With V=1 and j=4, the ordinary update √(4R) approaches R=4 monotonically.
  call state%reset()
  radius=1.0d0
  converged=.false.
  do iteration=2,100
     expected=sqrt(4.0d0*radius)
     proposal=expected
     call state%update(iteration,1,radius,proposal,fit)
     call Assert('monotonic proposal is unchanged',proposal,expected,relTol=1.0d-14)
     radius=proposal
     if (fit < tolerance) then
        converged=.true.
        exit
     end if
  end do
  call Assert('monotonic solve converges',converged,.true.)
  call Assert('monotonic analytic root',radius,4.0d0,relTol=2.0d0*tolerance)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine seedHistory(state,component)
    !!{RST
    Supply two actual input radii to the production update without triggering the oscillation handler.
    !!}
    implicit none
    type            (equilibriumIterationState), intent(inout) :: state
    integer                                    , intent(in   ) :: component
    double precision                                           :: proposal , fit

    proposal=4.0d0
    call state%update(2,component,4.0d0,proposal,fit)
    proposal=1.0d0
    call state%update(3,component,1.0d0,proposal,fit)
    return
  end subroutine seedHistory

  subroutine solveCycle(state,initial,radius,count,converged)
    !!{RST
    Solve the smooth two-cycle fixture using the production correction and convergence measure.
    !!}
    implicit none
    type            (equilibriumIterationState), intent(inout) :: state
    double precision                           , intent(in   ) :: initial
    double precision                           , intent(  out) :: radius
    integer                                    , intent(  out) :: count
    logical                                    , intent(  out) :: converged
    double precision                                           :: proposal , fit

    radius=initial
    converged=.false.
    do count=2,100
       proposal=4.0d0/radius
       call state%update(count,1,radius,proposal,fit)
       radius=proposal
       if (fit < tolerance) then
          converged=.true.
          return
       end if
    end do
    return
  end subroutine solveCycle

end program Test_Galactic_Structure_Equilibrium_Iteration
