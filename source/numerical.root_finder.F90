!% Contains a module which does root finding.

module Root_Finder
  !% Implements root finding.
  use FGSL
  private
  public :: Root_Find, Root_Find_Done

contains

  recursive double precision function Root_Find(lowerLimit,upperLimit,root,parameterPointer,rootFunction,rootFunctionSolver &
       &,toleranceAbsolute,toleranceRelative,rootSolver)
    !% Finds the root of the supplied {\tt rootFunction} function.
    use, intrinsic :: ISO_C_Binding                             
    use Galacticus_Error
    implicit none
    double precision,             external               :: root
    type(c_ptr),                  intent(in)             :: parameterPointer
    type(fgsl_function),          intent(inout)          :: rootFunction
    type(fgsl_root_fsolver),      intent(inout)          :: rootFunctionSolver
    double precision,             intent(in)             :: lowerLimit,upperLimit
    double precision,             intent(in),   optional :: toleranceAbsolute,toleranceRelative
    type(fgsl_root_fsolver_type), intent(in),   optional :: rootSolver
    double precision,             parameter              :: toleranceAbsoluteDefault=1.0d-10,toleranceRelativeDefault=1.0d-10
    integer,                      parameter              :: iterationMaximum=1000
    integer                                              :: status,iteration
    type(fgsl_root_fsolver_type)                         :: rootSolverActual
    double precision                                     :: toleranceAbsoluteActual,toleranceRelativeActual,xRoot,xLow,xHigh
    
    ! Set optional parameters if present, otherwise use defaults.
    if (present(toleranceAbsolute)) then
       toleranceAbsoluteActual=toleranceAbsolute
    else
       toleranceAbsoluteActual=toleranceAbsoluteDefault
    end if
    if (present(toleranceRelative)) then
       toleranceRelativeActual=toleranceRelative
    else
       toleranceRelativeActual=toleranceRelativeDefault
    end if
    if (present(rootSolver)) then
       rootSolverActual=rootSolver
    else
       rootSolverActual=FGSL_Root_fSolver_Brent
    end if
 
    ! Initialize the integration variables if necessary.
    if (.not.FGSL_Well_Defined(rootFunctionSolver)) then
       rootFunction      =FGSL_Function_Init(root,parameterPointer)       
       rootFunctionSolver=FGSL_Root_fSolver_Alloc(rootSolverActual)
    end if
    
    ! Do the integration
    status=FGSL_Root_fSolver_Set(rootFunctionSolver,rootFunction,lowerLimit,upperLimit)
    if (status /= FGSL_Success) call Galacticus_Error_Report('Root_Find','failed to initialize solver')
    iteration=0
    do
       iteration=iteration+1
       status=FGSL_Root_fSolver_Iterate(rootFunctionSolver)
       if (status /= FGSL_Success .or. iteration > iterationMaximum) exit
       xRoot=FGSL_Root_fSolver_Root(rootFunctionSolver)
       xLow =FGSL_Root_fSolver_x_Lower(rootFunctionSolver)
       xHigh=FGSL_Root_fSolver_x_Upper(rootFunctionSolver)
       status=FGSL_Root_Test_Interval(xLow,xHigh,toleranceAbsoluteActual,toleranceRelativeActual)
       if (status == FGSL_Success) exit
    end do
    if (status /= FGSL_Success) call Galacticus_Error_Report('Root_Find','failed to find root')
    Root_Find=xRoot
    return
  end function Root_Find
  
  subroutine Find_Root_Done(rootFunction,rootFunctionSolver)
    !% Frees up integration objects that are no longer required.
    implicit none
    type(fgsl_function),     intent(inout) :: rootFunction
    type(fgsl_root_fsolver), intent(inout) :: rootFunctionSolver

    call FGSL_Root_fSolver_Free(rootFunctionSolver)
    call FGSL_Function_Free(rootFunction)
    return
  end subroutine Find_Root_Done

end module Root_Finder
