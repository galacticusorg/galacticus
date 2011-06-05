!% Contains a module which defines internal error codes for the \glc\ ODE solver.

module ODE_Solver_Error_Codes
  !% Defines internal error codes for the \glc\ ODE solver.
  public

  ! An interrupt has been triggered.
  integer,         parameter :: odeSolverInterrupt=1001

  ! Point during the integration at which an interrupt occurred.
  double precision           :: interruptedAtX
  !$omp threadprivate(interruptedAtX)

end module ODE_Solver_Error_Codes
