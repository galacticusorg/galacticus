module ODE_Solver
  use ODE_Solver_Error_Codes
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: ODE_Solve, ODE_Solver_Free
  
contains
  
  subroutine ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,x0,x1,yCount,y,odeFunction,parameterPointer&
       &,toleranceAbsolute,toleranceRelative,reset)
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    double precision,         intent(in)              :: x1,toleranceAbsolute,toleranceRelative
    type(c_ptr),              intent(in)              :: parameterPointer
    integer,                  intent(in)              :: yCount
    double precision,         intent(inout)           :: x0,y(yCount)
    type(fgsl_odeiv_step),    intent(inout)           :: odeStepper
    type(fgsl_odeiv_control), intent(inout)           :: odeController
    type(fgsl_odeiv_evolve),  intent(inout)           :: odeEvolver
    type(fgsl_odeiv_system),  intent(inout)           :: odeSystem
    logical,                  intent(inout), optional :: reset
    integer(kind=4),          external                :: odeFunction
    integer                                           :: status
    integer(c_size_t)                                 :: odeNumber
    double precision                                  :: x,h,x1Internal
    logical                                           :: resetActual,forwardEvolve

    ! Number of ODEs to solve.
    odeNumber=yCount

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    if (resetActual.or..not.FGSL_Well_Defined(odeSystem)) then
       odeStepper   =FGSL_ODEiv_Step_Alloc   (FGSL_ODEiv_Step_RKF45,odeNumber)
       odeController=FGSL_ODEiv_Control_y_New(toleranceAbsolute,toleranceRelative)
       odeEvolver   =FGSL_ODEiv_Evolve_Alloc (odeNumber)
       odeSystem    =FGSL_ODEiv_System_Init  (odeFunction,odeNumber,parameterPointer)
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
       case (FGSL_Success)
          ! Successful completion of the step - do nothing.
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1Internal=interruptedAtX
       case default
          ! Some other error condition.
          call Galacticus_Error_Report('ODE_Solve','ODE integration failed')
       end select
    end do
    ! Return the new value of x.
    x0=x
    return
  end subroutine ODE_Solve
  
  subroutine ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
    !% Free up workspace allocated to ODE solving.
    implicit none
    type(fgsl_odeiv_step),    intent(inout) :: odeStepper
    type(fgsl_odeiv_control), intent(inout) :: odeController
    type(fgsl_odeiv_evolve),  intent(inout) :: odeEvolver
    type(fgsl_odeiv_system),  intent(inout) :: odeSystem
  
    call FGSL_ODEiv_Evolve_Free (odeEvolver   )
    call FGSL_ODEiv_Control_Free(odeController)
    call FGSL_ODEiv_Step_Free   (odeStepper   )
    call FGSL_ODEiv_System_Free (odeSystem    )
    return
  end subroutine ODE_Solver_Free

end module ODE_Solver
