!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module of integrands for unit tests.

module Test_Integration2_Functions
  !% Contains integrands for unit tests.
  use Numerical_Integration2
  use ISO_Varying_String
  implicit none
  private
  public :: testIntegrator , testFunctionsInitialize,                               &
       &    function1Scalar, function1Vector        , function1YEPPP, function1GSL, &
       &    function2Scalar, function2Vector        , function2YEPPP, function2GSL, &
       &    function3Scalar, function3Vector        , function3YEPPP, function3GSL

  type :: testIntegrator
     !% Type used for testing numerical integrators.
     class  (integrator    ), allocatable :: integrator_
     type   (varying_string)              :: description
     logical                              :: useYEPPP
     integer                              :: order
  end type testIntegrator

  type :: testFunction
     !% Type used for referencing functions.
     character       (len=16         )                  :: description
     double precision                                   :: rangeLow   , rangeHigh, solution
     procedure       (function1Scalar), pointer, nopass :: scalar
     procedure       (function1Vector), pointer, nopass :: vector
     procedure       (function1YEPPP ), pointer, nopass :: yeppp
     procedure       (function1GSL   ), pointer, nopass :: gsl
  end type testFunction

  ! Array of functions for integration tests.
  type(testFunction), public, dimension(3) :: testFunctions

contains

  subroutine testFunctionsInitialize()
    !% Initalize an array of test functions for integration tests.
    implicit none
    testFunctions=[                                                                                                                                   &
         &         testFunction('log(x) sin(x)   ',1.0d0,10.0d0, 1.549173238901735869d0,function1Scalar,function1Vector,function1YEPPP,function1GSL), &
         &         testFunction('1/sqrt(x)       ',1.0d0,10.0d0, 4.324555320336759000d0,function2Scalar,function2Vector,function2YEPPP,function2GSL), &
         &         testFunction('1/(10⁻³+[x-3]²) ',1.0d0,10.0d0,98.703068147327100000d0,function3Scalar,function3Vector,function3YEPPP,function3GSL)  &
         &        ]
    return
  end subroutine testFunctionsInitialize

  double precision function function1Scalar(x)
    !% Test function number 1 for numerical integration tests: scalar version.
    implicit none
    double precision, intent(in   ) :: x

    function1Scalar=log(x)*sin(x)
    return
  end function function1Scalar

  function function1Vector(x)
    !% Test function number 1 for numerical integration tests: vector version.
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function1Vector

    function1Vector=log(x)*sin(x)
    return
  end function function1Vector

 function function1YEPPP(x)
    !% Test function number 1 for numerical integration tests: YEPPP! vector version.
    use, intrinsic :: ISO_C_Binding
    use yepCore
    use yepMath
    implicit none
    real   (c_double), intent(in   ), dimension(     : ) :: x
    real   (c_double)               , dimension(size(x)) :: function1YEPPP
    real   (c_double)               , dimension(size(x)) :: y,z
    integer(c_int   )                                    :: s
    integer(c_size_t)                                    :: n
 
    n=size(x)
    s=yepMath_Log_V64f_V64f         (x,y      ,n)
    s=yepMath_Sin_V64f_V64f         (x,z      ,n)
    s=yepCore_Multiply_V64fV64f_V64f(y,z,function1YEPPP,n)
    return
  end function function1YEPPP

  function function1GSL(x,parameterPointer) bind(c)
    !% Test function number 1 for numerical integration tests: \gls{gsl} version.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: function1GSL
    real(kind=c_double), value :: x
    type(c_ptr        ), value :: parameterPointer

    function1GSL=log(x)*sin(x)
    return
  end function function1GSL

  double precision function function2Scalar(x)
    !% Test function number 1 for numerical integration tests: scalar version.
    implicit none
    double precision, intent(in   ) :: x

    function2Scalar=1.0d0/sqrt(x)
    return
  end function function2Scalar

  function function2Vector(x)
    !% Test function number 2 for numerical integration tests: vector version.
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function2Vector

    function2Vector=1.0d0/sqrt(x)
    return
  end function function2Vector

 function function2YEPPP(x)
    !% Test function number 2 for numerical integration tests: YEPPP! vector version.
    use, intrinsic :: ISO_C_Binding
    use yepCore
    use yepMath
    implicit none
    real   (c_double), intent(in   ), dimension(     : ) :: x
    real   (c_double)               , dimension(size(x)) :: function2YEPPP
    real   (c_double)               , dimension(size(x)) :: y,z
    integer(c_int   )                                    :: s
    integer(c_size_t)                                    :: n
 
    n=size(x)
    s=yepMath_Log_V64f_V64f         (x        ,y            ,n)
    s=yepCore_Multiply_V64fS64f_V64f(y,-0.5d0,z             ,n)
    s=yepMath_Exp_V64f_V64f         (z       ,function2YEPPP,n)
    return
  end function function2YEPPP

  function function2GSL(x,parameterPointer) bind(c)
    !% Test function number 2 for numerical integration tests: \gls{gsl} version.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: function2GSL
    real(kind=c_double), value :: x
    type(c_ptr        ), value :: parameterPointer

    function2GSL=1.0d0/sqrt(x)
    return
  end function function2GSL

  double precision function function3Scalar(x)
    !% Test function number 3 for numerical integration tests: scalar version.
    implicit none
    double precision, intent(in   ) :: x

    function3Scalar=1.0d0/(1.0d-3+(x-3.0d0)**2)
    return
  end function function3Scalar

  function function3Vector(x)
    !% Test function number 3 for numerical integration tests: vector version.
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function3Vector

    function3Vector=1.0d0/(1.0d-3+(x-3.0d0)**2)
    return
  end function function3Vector

 function function3YEPPP(x)
    !% Test function number 3 for numerical integration tests: YEPPP! vector version.
    use, intrinsic :: ISO_C_Binding
    use yepCore
    use yepMath
    implicit none
    real   (c_double), intent(in   ), dimension(     : ) :: x
    real   (c_double)               , dimension(size(x)) :: function3YEPPP
    real   (c_double)               , dimension(size(x)) :: w1,w2,w3
    integer(c_int   )                                    :: s
    integer(c_size_t)                                    :: n
 
    n=size(x)
    s=yepCore_Subtract_V64fS64f_V64f(x,3.0d0,w1,n)
    w2=w1**2
    s=yepCore_Add_V64fS64f_V64f(w2,1.0d-3,w3,n)
    function3YEPPP=1.0d0/w3
    return
  end function function3YEPPP

  function function3GSL(x,parameterPointer) bind(c)
    !% Test function number 3 for numerical integration tests: \gls{gsl} version.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: function3GSL
    real(kind=c_double), value :: x
    type(c_ptr        ), value :: parameterPointer

    function3GSL=1.0d0/(1.0d-3+(x-3.0d0)**2)
    return
  end function function3GSL


end module Test_Integration2_Functions
