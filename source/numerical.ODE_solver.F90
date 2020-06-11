!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which provides an interface to the GNU Scientific Library ODEIV differential equation solvers.

module ODE_Solver
  !% Contains an interface to the GNU Scientific Library ODEIV differential equation solvers.
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  implicit none
  private
  public :: ODE_Solve, ODE_Solver_Free

  ! Integrand interface.
  abstract interface
     integer function odesTemplate(x,y,dydx)
       double precision, intent(in   )               :: x
       double precision, intent(in   ), dimension(:) :: y
       double precision, intent(  out), dimension(:) :: dydx
     end function odesTemplate
  end interface

  ! Integrand function.
  procedure(odesTemplate), pointer :: currentODEs
  integer  (c_size_t    )          :: currentODENumber
  !$omp threadprivate(currentODEs,currentODENumber)

contains

  subroutine ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,x0,x1,yCount,y,odes,toleranceAbsolute,toleranceRelative,yScale,reset)
    !% Interface to the GNU Scientific Library ODEIV differential equation solvers.
    use            :: FGSL                  , only : FGSL_ODEiv_Control_Scaled_New, FGSL_ODEiv_Control_y_New, FGSL_ODEiv_Evolve_Alloc, FGSL_ODEiv_Evolve_Apply, &
          &                                          FGSL_ODEiv_Step_Alloc        , FGSL_ODEiv_Step_RKCK    , FGSL_ODEiv_System_Init , fgsl_odeiv_system      , &
          &                                          FGSL_Well_Defined            , fgsl_odeiv_control      , fgsl_odeiv_evolve      , fgsl_odeiv_step
    use            :: Galacticus_Error      , only : Galacticus_Error_Report
    use            :: Interface_GSL         , only : GSL_Success
    use, intrinsic :: ISO_C_Binding         , only : c_ptr                        , c_size_t
    use            :: ODE_Solver_Error_Codes, only : interruptedAtX               , odeSolverInterrupt
    implicit none
    double precision                    , intent(in   )           :: toleranceAbsolute              , toleranceRelative              , x1
    integer                             , intent(in   )           :: yCount
    double precision                    , intent(inout)           :: x0                             , y                (yCount)
    double precision                    , intent(in   ), optional :: yScale           (yCount)
    type            (fgsl_odeiv_step   ), intent(inout)           :: odeStepper
    type            (fgsl_odeiv_control), intent(inout)           :: odeController
    type            (fgsl_odeiv_evolve ), intent(inout)           :: odeEvolver
    type            (fgsl_odeiv_system ), intent(inout)           :: odeSystem
    logical                             , intent(inout), optional :: reset
    procedure       (odesTemplate      )                          :: odes
    procedure       (odesTemplate      ), pointer                 :: previousODEs
    double precision                    , parameter               :: dydtScaleUniform         =0.0d0, yScaleUniform            =1.0d0
    integer                                                       :: status
    integer         (kind=c_size_t     )                          :: previousODENumber
    double precision                                              :: h                              , x                              , x1Internal
    logical                                                       :: forwardEvolve                  , resetActual
    type            (c_ptr             )                          :: parameterPointer

    ! Store the current ODE function (and system size) so that we can restore it on exit. This allows the ODE function to be
    ! called recursively.
    previousODEs      => currentODEs
    previousODENumber =  currentODENumber
    currentODEs       => ODEs
    currentODENumber  =  yCount
    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    if (resetActual.or..not.FGSL_Well_Defined(odeSystem)) then
       odeStepper      =FGSL_ODEiv_Step_Alloc        (FGSL_ODEiv_Step_RKCK,currentODENumber)
       if (present(yScale)) then
          odeController=FGSL_ODEiv_Control_Scaled_New(toleranceAbsolute,toleranceRelative,yScaleUniform,dydtScaleUniform,yScale,currentODENumber)
       else
          odeController=FGSL_ODEiv_Control_y_New     (toleranceAbsolute,toleranceRelative                                                )
       end if
       odeEvolver      =FGSL_ODEiv_Evolve_Alloc      (currentODENumber)
       odeSystem       =FGSL_ODEiv_System_Init       (odesWrapper,currentODENumber,parameterPointer)
    end if
    ! Keep a local copy of the end point as we may reset it.
    x1Internal=x1
    ! Set initial value of x variable.
    x=x0
    ! Make initial guess for timestep.
    h=(x1-x0)
    ! Determine if we want forward or backward evolution.
    forwardEvolve=x1>x0
    ! Evolve the system until the final time is reached.
    do while ((forwardEvolve.and.x<x1Internal).or.(.not.forwardEvolve.and.x>x1Internal))
       status=FGSL_ODEiv_Evolve_Apply(odeEvolver,odeController,odeStepper,odeSystem,x,x1Internal,h,y)
       select case (status)
       case (GSL_Success)
          ! Successful completion of the step - do nothing.
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1Internal=interruptedAtX
       case default
          ! Some other error condition.
          call Galacticus_Error_Report('ODE integration failed'//{introspection:location})
       end select
    end do
    ! Return the new value of x.
    x0=x
    ! Restore the previous OODEs.
    currentODEs      => previousODEs
    currentODENumber =  previousODENumber
    return
  end subroutine ODE_Solve

  function odesWrapper(x,y,dydx,parameterPointer) bind(c)
    !% Wrapper function used for \gls{gsl} ODE functions.
    use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_ptr
    implicit none
    integer(kind=c_int   )                              :: odesWrapper
    real   (kind=c_double), value                       :: x
    real   (kind=c_double), dimension(*), intent(in   ) :: y
    real   (kind=c_double), dimension(*)                :: dydx
    type   (     c_ptr   ), value                       :: parameterPointer
    !$GLC attributes unused :: parameterPointer

    odesWrapper=currentODEs(x,y(1:currentODENumber),dydx(1:currentODENumber))
    return
  end function odesWrapper

  subroutine ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
    !% Free up workspace allocated to ODE solving.
    use :: FGSL, only : FGSL_ODEiv_Control_Free, FGSL_ODEiv_Evolve_Free, FGSL_ODEiv_Step_Free, FGSL_ODEiv_System_Free, &
          &             fgsl_odeiv_control     , fgsl_odeiv_evolve     , fgsl_odeiv_step     , fgsl_odeiv_system
    implicit none
    type(fgsl_odeiv_step   ), intent(inout) :: odeStepper
    type(fgsl_odeiv_control), intent(inout) :: odeController
    type(fgsl_odeiv_evolve ), intent(inout) :: odeEvolver
    type(fgsl_odeiv_system ), intent(inout) :: odeSystem

    call FGSL_ODEiv_Evolve_Free (odeEvolver   )
    call FGSL_ODEiv_Control_Free(odeController)
    call FGSL_ODEiv_Step_Free   (odeStepper   )
    call FGSL_ODEiv_System_Free (odeSystem    )
    return
  end subroutine ODE_Solver_Free

end module ODE_Solver
