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

!% Contains a module which provides an interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.

module ODEIV2_Solver
  !% Contains an interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.
  use ODE_Solver_Error_Codes
  use FGSL
  use FODEIV2
  use, intrinsic :: ISO_C_Binding
  private
  public :: ODEIV2_Solve, ODEIV2_Solver_Free

contains

  subroutine ODEIV2_Solve(odeDriver,odeSystem,x0,x1,yCount,y,odeFunction,parameterPointer,toleranceAbsolute,toleranceRelative&
#ifdef PROFILE
       &,Error_Analyzer &
#endif
       &,yScale,errorHandler,algorithm,reset,odeStatus)
    !% Interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    use ISO_Varying_String
    use String_Handling
    implicit none
    double precision                   , intent(in   )                    :: toleranceAbsolute        , toleranceRelative        , x1
    type            (c_ptr            ), intent(in   )                    :: parameterPointer
    integer                            , intent(in   )                    :: yCount
    double precision                   , intent(inout)                    :: x0                       , y                (yCount)
    double precision                   , intent(in   ), optional          :: yScale           (yCount)
    type            (fodeiv2_driver   ), intent(inout)                    :: odeDriver
    type            (fodeiv2_system   ), intent(inout)                    :: odeSystem
    logical                            , intent(inout), optional          :: reset
    procedure       (                 )               , optional, pointer :: errorHandler
    type            (fodeiv2_step_type), intent(in   ), optional          :: algorithm
    integer                            , intent(  out), optional          :: odeStatus
#ifdef PROFILE
    type(c_funptr), intent(in   ) :: Error_Analyzer
#endif
    integer         (kind=4           ), external  :: odeFunction
    integer                            , parameter :: genericFailureCountMaximum=10
    double precision                   , parameter :: dydtScaleUniform          =0.0d0, yScaleUniform=1.0d0
    integer                                        :: status
    integer         (kind=c_size_t    )            :: odeNumber
    double precision                               :: h                               , x                   , &
         &                                            x1Internal
    logical                                        :: forwardEvolve                   , resetActual
    type            (fodeiv2_step_type)            :: algorithmActual
    type            (varying_string   )            :: message

    ! Number of ODEs to solve.
    odeNumber=yCount

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    if (resetActual.or..not.FODEIV2_Driver_Status(odeDriver)) then
       ! Make initial guess for timestep.
       h=(x1-x0)
       odeSystem=FODEIV2_System_Init(odeFunction,odeNumber,parameterPointer)
       ! Select the algorithm to use.
       if (present(algorithm)) then
          algorithmActual=algorithm
       else
          algorithmActual=Fodeiv2_Step_RKCK
       end if
       if (present(yScale)) then
          ! Scales for the absolute tolerance have been given, so use them.
          odeDriver=FODEIV2_Driver_Alloc_Scaled_New(odeSystem,algorithmActual,h,toleranceAbsolute,toleranceRelative&
               &,yScaleUniform,dydtScaleUniform,yScale)
       else
          ! No scales given, assume they are all unity.
          odeDriver=FODEIV2_Driver_Alloc_y_New     (odeSystem,algorithmActual,h,toleranceAbsolute,toleranceRelative)
       end if
    end if
    ! Keep a local copy of the end point as we may reset it.
    x1Internal=x1
    ! Set initial value of x variable.
    x=x0
    ! Determine if we want forward or backward evolution.
    forwardEvolve=x1>x0
    ! Reset the driver.
    status=FODEIV2_Driver_Reset(odeDriver)
    ! Evolve the system until the final time is reached.
    do while ((forwardEvolve.and.x<x1Internal).or.(.not.forwardEvolve.and.x>x1Internal))
      status=FODEIV2_Driver_Apply(odeDriver,x,x1Internal,y&
#ifdef PROFILE
           &,Error_Analyzer&
#endif
           &)
       select case (status)
       case (FGSL_Success)
          ! Successful completion of the step - do nothing except resetting failure count.
       case (FGSL_Failure)
          ! Generic failure - most likely a stepsize underflow.
          if (present(errorHandler)) call errorHandler()
          ! If ODE status was requested, then return it instead of aborting.
          if (present(odeStatus)) then
             x0=x
             odeStatus=status
             return
          end if
          message='ODE integration failed with status '
          message=message//status//' [generic failure]'//char(10)
          message=message//' => most likely a stepsize underflow'
          call Galacticus_Error_Report('ODEIV2_Solve',message)
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1Internal=interruptedAtX
       case default
          ! Some other error condition.
          if (present(errorHandler)) call errorHandler()
          ! If ODE status was requested, then return it instead of aborting.
          if (present(odeStatus)) then
             x0=x
             odeStatus=status
             return
          end if
          message='ODE integration failed with status '
          message=message//status
          call Galacticus_Error_Report('ODEIV2_Solve',message)
       end select
    end do
    ! Return the new value of x.
    x0=x
    if (present(odeStatus)) odeStatus=status
    return
  end subroutine ODEIV2_Solve

  subroutine ODEIV2_Solver_Free(odeDriver,odeSystem)
    !% Free up workspace allocated to ODE solving.
    implicit none
    type(fodeiv2_driver), intent(inout) :: odeDriver
    type(fodeiv2_system), intent(inout) :: odeSystem

    call Fodeiv2_Driver_Free(odeDriver)
    call Fodeiv2_System_Free(odeSystem)
    return
  end subroutine ODEIV2_Solver_Free

end module ODEIV2_Solver
