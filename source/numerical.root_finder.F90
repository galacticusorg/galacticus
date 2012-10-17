!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which does root finding.

module Root_Finder
  !% Implements root finding.
  use FGSL
  implicit none
  private
  public :: Root_Find, Root_Find_Done

contains

  recursive double precision function Root_Find(lowerLimit,upperLimit,root,parameterPointer,rootFunction,rootFunctionSolver &
       &,toleranceAbsolute,toleranceRelative,rootSolver,reset)
    !% Finds the root of the supplied {\tt rootFunction} function.
    use, intrinsic :: ISO_C_Binding                             
    use Galacticus_Error
    implicit none
    double precision,             external                :: root
    type(c_ptr),                  intent(in)              :: parameterPointer
    type(fgsl_function),          intent(inout)           :: rootFunction
    type(fgsl_root_fsolver),      intent(inout)           :: rootFunctionSolver
    double precision,             intent(in)              :: lowerLimit,upperLimit
    double precision,             intent(in),    optional :: toleranceAbsolute,toleranceRelative
    type(fgsl_root_fsolver_type), intent(in),    optional :: rootSolver
    logical,                      intent(inout), optional :: reset
    double precision,             parameter               :: toleranceAbsoluteDefault=1.0d-10,toleranceRelativeDefault=1.0d-10
    integer,                      parameter               :: iterationMaximum=1000
    integer                                               :: status,iteration
    type(fgsl_root_fsolver_type)                          :: rootSolverActual
    double precision                                      :: toleranceAbsoluteActual,toleranceRelativeActual,xRoot,xLow,xHigh
    logical                                               :: resetActual
    
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
 
    ! Determine whether to reset or not.
    if (present(reset)) then
       resetActual=reset
       reset      =.false.
    else
       resetActual=.false.
    end if

    ! Initialize the integration variables if necessary.
    if (.not.FGSL_Well_Defined(rootFunctionSolver).or.resetActual) then
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
  
  subroutine Root_Find_Done(rootFunction,rootFunctionSolver)
    !% Frees up integration objects that are no longer required.
    implicit none
    type(fgsl_function),     intent(inout) :: rootFunction
    type(fgsl_root_fsolver), intent(inout) :: rootFunctionSolver

    call FGSL_Root_fSolver_Free(rootFunctionSolver)
    call FGSL_Function_Free(rootFunction)
    return
  end subroutine Root_Find_Done

end module Root_Finder
