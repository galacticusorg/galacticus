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
       &,yScale,Error_Analyzer,reset)
    !% Interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    use ISO_Varying_String
    use String_Handling
    implicit none
    double precision,         intent(in   )           :: x1,toleranceAbsolute,toleranceRelative
    type(c_ptr),              intent(in   )           :: parameterPointer
    integer,                  intent(in   )           :: yCount
    double precision,         intent(inout)           :: x0,y(yCount)
    double precision,         intent(in   ), optional :: yScale(yCount)
    type(fodeiv2_driver),     intent(inout)           :: odeDriver
    type(fodeiv2_system),     intent(inout)           :: odeSystem
    logical,                  intent(inout), optional :: reset
    procedure(),              pointer,       optional :: Error_Analyzer
    integer(kind=4),          external                :: odeFunction
    integer,                  parameter               :: genericFailureCountMaximum=10
    double precision,         parameter               :: yScaleUniform=1.0d0, dydtScaleUniform=0.0d0
    integer                                           :: status
    integer(c_size_t)                                 :: odeNumber
    double precision                                  :: x,h,x1Internal
    logical                                           :: resetActual,forwardEvolve
    type(varying_string)                              :: message
#ifdef PROFILE
    real(fgsl_double),        dimension(yCount)       :: yError
#endif

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
       if (present(yScale)) then
          ! Scales for the absolute tolerance have been given, so use them.
          odeDriver=FODEIV2_Driver_Alloc_Scaled_New(odeSystem,Fodeiv2_Step_RKCK,h,toleranceAbsolute,toleranceRelative&
               &,yScaleUniform,dydtScaleUniform,yScale)
       else
          ! No scales given, assume they are all unity.
          odeDriver=FODEIV2_Driver_Alloc_y_New     (odeSystem,Fodeiv2_Step_RKCK,h,toleranceAbsolute,toleranceRelative)
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
      status=FODEIV2_Driver_Apply(odeDriver,x,x1Internal,y)

#ifdef PROFILE
       ! If profiling is being performed, extract errors and send them to the specified error analysis function.
       if (present(Error_Analyzer) .and. x /= x0) then
          call FODEIV2_Driver_Error(odeDriver,yError)
          call Error_Analyzer(y,yError,h,status)
       end if
#endif

       select case (status)
       case (FGSL_Success)
          ! Successful completion of the step - do nothing except resetting failure count.
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1Internal=interruptedAtX
       case default
          ! Some other error condition.
          message='ODE integration failed with status '
          message=message//status
          call Galacticus_Error_Report('ODEIV2_Solve',message)
       end select
    end do
    ! Return the new value of x.
    x0=x
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
