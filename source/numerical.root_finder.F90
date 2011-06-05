!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
