!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a program to test ODE solver routines.

program Test_ODE_Solver
  !% Tests that ODE solver routines work.
  use Unit_Tests
  use ISO_Varying_String
  use FGSL
    use ODE_Solver
  use Test_ODE_Solver_Functions
  use, intrinsic :: ISO_C_Binding
  implicit none
  double precision,         dimension(10) :: xEnd
  double precision,         dimension( 2) :: y
  type(fgsl_odeiv_step)                   :: odeStepper
  type(fgsl_odeiv_control)                :: odeController
  type(fgsl_odeiv_evolve)                 :: odeEvolver
  type(fgsl_odeiv_system)                 :: odeSystem
  logical                                 :: odeReset
  type(c_ptr)                             :: parameterPointer
  integer                                 :: i
  double precision                        :: xStart
  character(len=32)                       :: message

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("ODE solver")

  ! Sinusoid.
  call Unit_Tests_Begin_Group("y'=sin(x)")
  xEnd=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  do i=1,size(xEnd)
     odeReset=.true.
     y(1:1)=[0.0d0]
     xStart=0.0d0
     call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,xStart,xEnd(i),1,y(1:1)&
          &,ODE_Set_1,parameterPointer,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-9,reset=odeReset)
     call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
     write (message,'(a,f4.1)') "x=0 to ",xEnd(i)
     call Assert(trim(message),y(1:1),1.0d0-cos(xEnd(i:i)),relTol=1.0d-6)
  end do
  call Unit_Tests_End_Group()

  ! Harmonic oscillator.
  call Unit_Tests_Begin_Group("y''=-y")
  xEnd=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  do i=1,size(xEnd)
     odeReset=.true.
     y(1:2)=[1.0d0,0.0d0]
     xStart=0.0d0
     call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,xStart,xEnd(i),2,y(1:2)&
          &,ODE_Set_2,parameterPointer,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-9,reset=odeReset)
     call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
     write (message,'(a,f4.1)') "x=0 to ",xEnd(i)
     call Assert(trim(message),y(1:2),[cos(xEnd(i)),-sin(xEnd(i))],relTol=1.0d-6)
  end do
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_ODE_Solver
