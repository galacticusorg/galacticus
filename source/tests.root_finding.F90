!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a program to test root finding routines.

program Test_Root_Finding
  !% Tests that routine finding routines work.
  use Unit_Tests
  use ISO_Varying_String
  use FGSL
  use Root_Finder
  use Test_Root_Finding_Functions
  use, intrinsic :: ISO_C_Binding
  implicit none
  type(fgsl_function)     :: rootFunction
  type(fgsl_root_fsolver) :: rootFunctionSolver
  logical                 :: rootFunctionReset
  type(c_ptr)             :: parameterPointer
  double precision        :: xRoot,xMinimum,xMaximum

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Root finding")

  ! Test root finding.
  xMinimum=-1.0d0
  xMaximum= 1.0d0
  rootFunctionReset=.true.
  xRoot=Root_Find(xMinimum,xMaximum,Root_Function_1,parameterPointer &
       &,rootFunction,rootFunctionSolver,reset=rootFunctionReset,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-6)
  call Assert('root of f(x)=x',xRoot,0.0d0,absTol=1.0d-6,relTol=1.0d-6)
  call Root_Find_Done(rootFunction,rootFunctionSolver)

  xMinimum=-1.0d0
  xMaximum= 1.0d0
  rootFunctionReset=.true.
  xRoot=Root_Find(xMinimum,xMaximum,Root_Function_2,parameterPointer &
       &,rootFunction,rootFunctionSolver,reset=rootFunctionReset,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-6)
  call Assert('root of f(x)=x²-5x+1 in range -1<x< 1',xRoot,0.5d0*(5.0d0-sqrt(21.0d0)),absTol=1.0d-6,relTol=1.0d-6)
  call Root_Find_Done(rootFunction,rootFunctionSolver)

  xMinimum= 2.0d0
  xMaximum=10.0d0
  rootFunctionReset=.true.
  xRoot=Root_Find(xMinimum,xMaximum,Root_Function_2,parameterPointer &
       &,rootFunction,rootFunctionSolver,reset=rootFunctionReset,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-6)
  call Assert('root of f(x)=x²-5x+1 in range  2<x<10',xRoot,0.5d0*(5.0d0+sqrt(21.0d0)),absTol=1.0d-6,relTol=1.0d-6)
  call Root_Find_Done(rootFunction,rootFunctionSolver)

  xMinimum=-1.0d0
  xMaximum= 1.0d0
  rootFunctionReset=.true.
  xRoot=Root_Find(xMinimum,xMaximum,Root_Function_3,parameterPointer &
       &,rootFunction,rootFunctionSolver,reset=rootFunctionReset,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-6)
  call Assert('root of f(x)=x*exp(-x)+1',xRoot,-0.567143d0,absTol=1.0d-6,relTol=1.0d-6)
  call Root_Find_Done(rootFunction,rootFunctionSolver)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Root_Finding
