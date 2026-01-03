!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module of integrands for unit tests.
!!}

module Test_Integration2_Functions
  !!{
  Contains integrands for unit tests.
  !!}
  use :: ISO_Varying_String    , only : varying_string
  use :: Numerical_Integration2, only : integrator2   , integratorMulti
  implicit none
  private
  public :: testIntegrator         , testIntegratorMulti, &
       &    function1Scalar        , function1Vector    , &
       &    function2Scalar        , function2Vector    , &
       &    function3Scalar        , function3Vector    , &
       &    function4Vector        , function14Vector   , &
       &    testFunctionsInitialize

  type :: testIntegrator
     !!{
     Type used for testing numerical integrators.
     !!}
     class  (integrator2   ), allocatable :: integrator_
     type   (varying_string)              :: description
     integer                              :: order
  end type testIntegrator

  type :: testFunction
     !!{
     Type used for referencing functions.
     !!}
     character       (len=22         )                  :: description
     double precision                                   :: rangeLow   , rangeHigh, solution
     procedure       (function1Scalar), pointer, nopass :: scalar
     procedure       (function1Vector), pointer, nopass :: vector
  end type testFunction

  type :: testIntegratorMulti
     !!{
     Type used for testing multi-integrand numerical integrators.
     !!}
     class  (integratorMulti), allocatable :: integrator_
     type   (varying_string )              :: description
     integer                               :: order
  end type testIntegratorMulti

  type :: testFunctionMulti
     !!{
     Type used for referencing functions.
     !!}
     character       (len=22         )                             :: description
     double precision                                              :: rangeLow   , rangeHigh
     double precision                  , allocatable, dimension(:) :: solution
     procedure       (function14Vector), pointer    , nopass       :: vector
  end type testFunctionMulti

  ! Arrays of functions for integration tests.
  type(testFunction     ), public, dimension(4) :: testFunctions
  type(testFunctionMulti), public, dimension(1) :: testFunctionsMulti

contains

  subroutine testFunctionsInitialize()
    !!{
    Initialize an array of test functions for integration tests.
    !!}
    implicit none
    testFunctions     =[                                                                                                                                                      &
         &              testFunction     ('log(x) sin( x)    ',1.0d0,10.0d0,  1.549173238901735869d0                        ,function1Scalar,function1Vector               ), &
         &              testFunction     ('1/sqrt(x)         ',1.0d0,10.0d0,  4.324555320336759000d0                        ,function2Scalar,function2Vector               ), &
         &              testFunction     ('1/(10⁻³+[x-3]²)   ',1.0d0,10.0d0, 98.703068147327100000d0                        ,function3Scalar,function3Vector               ), &
         &              testFunction     ('log(x) sin(2x)    ',1.0d0,10.0d0, -0.659102340089651400d0                        ,function4Scalar,function4Vector               )  &
         &             ]
    testFunctionsMulti=[                                                                                                                                                      &
         &              testFunctionMulti('log(x) sin({1,2}x)',1.0d0,10.0d0,[ 1.549173238901735869d0, -0.6591023400896514d0]                ,function14Vector              )  &
         &             ]
    return
  end subroutine testFunctionsInitialize

  double precision function function1Scalar(x)
    !!{
    Test function number 1 for numerical integration tests: scalar version.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    function1Scalar=log(x)*sin(x)
    return
  end function function1Scalar

  function function1Vector(x)
    !!{
    Test function number 1 for numerical integration tests: vector version.
    !!}
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function1Vector

    function1Vector=log(x)*sin(x)
    return
  end function function1Vector

  double precision function function2Scalar(x)
    !!{
    Test function number 1 for numerical integration tests: scalar version.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    function2Scalar=1.0d0/sqrt(x)
    return
  end function function2Scalar

  function function2Vector(x)
    !!{
    Test function number 2 for numerical integration tests: vector version.
    !!}
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function2Vector

    function2Vector=1.0d0/sqrt(x)
    return
  end function function2Vector

  double precision function function3Scalar(x)
    !!{
    Test function number 3 for numerical integration tests: scalar version.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    function3Scalar=1.0d0/(1.0d-3+(x-3.0d0)**2)
    return
  end function function3Scalar

  function function3Vector(x)
    !!{
    Test function number 3 for numerical integration tests: vector version.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function3Vector

    function3Vector=1.0d0/(1.0d-3+(x-3.0d0)**2)
    return
  end function function3Vector

  double precision function function4Scalar(x)
    !!{
    Test function number 4 for numerical integration tests: scalar version.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    function4Scalar=log(x)*sin(2.0d0*x)
    return
  end function function4Scalar

  function function4Vector(x)
    !!{
    Test function number 4 for numerical integration tests: vector version.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: function4Vector

    function4Vector=log(x)*sin(2.0d0*x)
    return
  end function function4Vector

  subroutine function14Vector(n,x,e,integrand)
    !!{
    Combined functions number 1 and 4 for numerical integration tests: vector version.
    !!}
    implicit none
    integer         , intent(in   )                       :: n
    double precision, intent(in   ), dimension(       : ) :: x
    logical         , intent(inout), dimension(       : ) :: e
    double precision, intent(  out), dimension(n,size(x)) :: integrand
    double precision               , dimension(  size(x)) :: logx

    logx                    =     log(      x)
    if (e(1)) integrand(1,:)=logx*sin(      x)
    if (e(2)) integrand(2,:)=logx*sin(2.0d0*x)
    return
  end subroutine function14Vector

end module Test_Integration2_Functions
