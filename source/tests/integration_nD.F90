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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests multi-dimensional integration.
!!}

program Test_Integration_ND
  !!{RST
  Tests the adaptive multi-dimensional (Genz-Malik) integrator.

  The degree of the cubature rule is tested directly, by limiting the number of integrand evaluations to exactly the number of
  points in one region so that the rule is applied once and never refined. Subdivision would otherwise mask a rule of the wrong
  degree, since an adaptive scheme built on a lower-degree rule still converges to the correct answer---it simply costs more.
  The remaining tests exercise the adaptive scheme itself against analytically-known integrals.
  !!}
  use :: Display                     , only : displayMessage      , displayVerbositySet  , verbosityLevelStandard
  use :: Error                       , only : errorStatusSuccess
  use :: ISO_Varying_String          , only : varying_string      , var_str              , operator(//)          , assignment(=)
  use :: Numerical_Constants_Math    , only : Pi
  use :: Numerical_Integration_nD    , only : integratorGenzMalikND
  use :: String_Handling             , only : operator(//)
  use :: Test_Integration_ND_Functions, only : integrandPolynomial5, integrandPolynomial7 , integrandProduct4D    , &
       &                                       integrandRootXYSquare,integrandSeparable3D , integrandSinCos       , &
       &                                       integrandSphere     , integrandXSquaredCos , integrandYOverRootX
  use :: Unit_Tests                  , only : Assert              , Unit_Tests_Begin_Group, Unit_Tests_End_Group  , &
       &                                      Unit_Tests_Finish
  implicit none
  type            (integratorGenzMalikND)               :: integrator_
  double precision                                      :: integral
  integer                                               :: status

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Numerical integration (nD)")

  ! Degree of the cubature rule. The evaluation budget is set to exactly the number of points in a single region, 2ᵈ+2d²+2d+1,
  ! so that the rule is applied once and the result is that of the rule alone.
  call Unit_Tests_Begin_Group("cubature rule")
  !! The embedded degree-5 rule integrates a polynomial of total degree five exactly, and so agrees with the degree-7 rule. The
  !! error estimate is their difference, so it must vanish, and the integration must therefore converge after a single region.
  integrator_=integratorGenzMalikND(3,integrandPolynomial5,toleranceRelative=1.0d-6)
  integral   =integrator_%evaluate([0.0d0,0.0d0,0.0d0],[1.0d0,1.0d0,1.0d0])
  call Assert("degree-5 polynomial over the unit cube is integrated exactly"      ,integral                     ,29.0d0/9.0d0     ,relTol=1.0d-12)
  call Assert("degree-5 polynomial converges after a single region"               ,integrator_%countEvaluated(),33                             )
  !! The degree-7 rule integrates a polynomial of total degree seven exactly.
  integrator_=integratorGenzMalikND(3,integrandPolynomial7,toleranceRelative=1.0d-6,countEvaluationsMaximum=33)
  integral   =integrator_%evaluate([0.0d0,0.0d0,0.0d0],[1.0d0,1.0d0,1.0d0],status=status)
  call Assert("degree-7 polynomial is integrated exactly by a single application of the rule",integral         ,829.0d0/2520.0d0 ,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! Two-dimensional integrals, using the integrands of the nested two-dimensional integrator's own tests.
  call Unit_Tests_Begin_Group("two dimensions")
  integrator_=integratorGenzMalikND(2,integrandSinCos      ,toleranceRelative=1.0d-8)
  integral   =integrator_%evaluate([-1.0d0,-1.0d0],[+1.0d0,+1.0d0])
  call report("∫sin(x²)cos(y)dxdy over -1<x<1, -1<y<1",integral,1.0443270936233882d0,integrator_%countEvaluated())
  integrator_=integratorGenzMalikND(2,integrandXSquaredCos ,toleranceRelative=1.0d-8)
  integral   =integrator_%evaluate([-1.0d0,-1.0d0],[+1.0d0,+1.0d0])
  call report("∫x²cos(y)dxdy over -1<x<1, -1<y<1"     ,integral,1.1219613130771953d0,integrator_%countEvaluated())
  integrator_=integratorGenzMalikND(2,integrandRootXYSquare,toleranceRelative=1.0d-6)
  integral   =integrator_%evaluate([+0.0d0,-1.0d0],[+1.0d0,+1.0d0])
  call report("∫√x y²dxdy over 0<x<1, -1<y<1"         ,integral,0.4444444444444444d0,integrator_%countEvaluated())
  !! This integrand is singular at x=0. The cubature never places a point on the boundary of a region, so the singularity is
  !! never evaluated, but convergence is slow and only a loose tolerance is asked for.
  integrator_=integratorGenzMalikND(2,integrandYOverRootX  ,toleranceRelative=1.0d-4)
  integral   =integrator_%evaluate([+0.0d0,+0.0d0],[+1.0d0,+2.0d0])
  call report("∫y/√x dxdy over 0<x<1, 0<y<2"          ,integral,4.0d0               ,integrator_%countEvaluated(),toleranceTest=1.0d-4)
  call Unit_Tests_End_Group()

  ! Three-dimensional integrals.
  call Unit_Tests_Begin_Group("three dimensions")
  integrator_=integratorGenzMalikND(3,integrandSeparable3D,toleranceRelative=1.0d-8)
  integral   =integrator_%evaluate([0.0d0,0.0d0,0.0d0],[1.0d0,1.0d0,1.0d0])
  call report("∫exp(x)sin(y)cos(z)dxdydz over the unit cube",integral,0.6646696797813771d0,integrator_%countEvaluated())
  !! A discontinuous integrand: a unit-radius sphere inscribed in the cube, whose integral is the volume of the sphere.
  !!
  !! Only a loose tolerance is asked for here, and deliberately so. A cubature refines d-dimensional boxes, but the
  !! discontinuity is a surface of dimension d-1, so the regions straddling it must be shrunk until their total volume is itself
  !! within tolerance: the cost grows as ε^-(d-1), which is 10⁹ evaluations for a relative tolerance of 10⁻³ in three
  !! dimensions. A nest of one-dimensional integrators does far better on such an integrand, because the innermost integral sees
  !! the discontinuity as a step function in one variable and resolves it to machine precision by bisection, leaving a
  !! continuous integrand for the levels outside it. This integrator should not be used for discontinuous integrands.
  integrator_=integratorGenzMalikND(3,integrandSphere,toleranceRelative=1.0d-1)
  integral   =integrator_%evaluate([-1.0d0,-1.0d0,-1.0d0],[+1.0d0,+1.0d0,+1.0d0])
  call report("∫dV over a unit-radius sphere inscribed in the cube",integral,4.0d0*Pi/3.0d0,integrator_%countEvaluated(),toleranceTest=1.0d-1)
  call Unit_Tests_End_Group()

  ! Four dimensions, to confirm the rule is correct in a dimensionality neither of the existing nested integrators reaches.
  call Unit_Tests_Begin_Group("four dimensions")
  integrator_=integratorGenzMalikND(4,integrandProduct4D,toleranceRelative=1.0d-8)
  integral   =integrator_%evaluate([0.0d0,0.0d0,0.0d0,0.0d0],[1.0d0,1.0d0,1.0d0,1.0d0])
  call report("∫x y² z³ w⁴ over the unit hypercube",integral,1.0d0/120.0d0,integrator_%countEvaluated())
  call Unit_Tests_End_Group()

  call Unit_Tests_End_Group()
  call Unit_Tests_Finish  ()

contains

  subroutine report(label,integral,integralExact,countEvaluations,toleranceTest)
    !!{RST
    Assert that  integral matches  integralExact, and report the number of integrand evaluations used.
    !!}
    implicit none
    character       (len=*         ), intent(in   )           :: label
    double precision                , intent(in   )           :: integral        , integralExact
    integer                         , intent(in   )           :: countEvaluations
    double precision                , intent(in   ), optional :: toleranceTest
    type            (varying_string)                          :: message
    !![
    <optionalArgument name="toleranceTest" defaultsTo="1.0d-6"/>
    !!]

    message=var_str('  ')//label//': '//countEvaluations//' evaluations'
    call displayMessage(message,verbosityLevelStandard)
    call Assert(label,integral,integralExact,relTol=toleranceTest_)
    return
  end subroutine report

end program Test_Integration_ND
