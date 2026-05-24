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
Module providing root-callback functions and shared state for the \mono{rootFinder}
benchmark program. Kept separate from the program unit so the callbacks can be
passed as \mono{procedure(...)} actual arguments to \refClass{rootFinder} without
running into Fortran restrictions on procedure pointers to internal procedures.
!!}

module Benchmark_Root_Finder_Functions
  !!{
  Callbacks and shared state for the \refClass{rootFinder} benchmark.

  The functions read a single module-level shift, \mono{bmShift}, which is
  rewritten by the benchmark driver before every call. This perturbs the
  problem (so the optimizer cannot constant-fold the solver chain) and keeps
  the cache state of \mono{rootFinder}'s wrapper realistic from one call to
  the next.
  !!}
  implicit none
  private
  public :: bmShift, fQuad, fQuadDerivative, fQuadBoth, &
       &    fFar

  ! Module-level shift consumed by every callback. Updated by the benchmark
  ! driver immediately before each find() invocation.
  double precision, save :: bmShift=0.0d0

contains

  double precision function fQuad(x)
    !!{
    Quadratic with a single root in the interval $[0,3]$.

    $f(x) = (x - (\mathrm{bmShift}+1))(x + 5)$ has roots at $x=\mathrm{bmShift}+1$
    and $x=-5$. With $\mathrm{bmShift} \in [-0.4,0.4]$ the in-bracket root sits
    in $[0.6,1.4]$.
    !!}
    double precision, intent(in   ) :: x

    fQuad=(x-(bmShift+1.0d0))*(x+5.0d0)
    return
  end function fQuad

  double precision function fQuadDerivative(x)
    !!{
    Derivative of \mono{fQuad}.
    !!}
    double precision, intent(in   ) :: x

    fQuadDerivative=2.0d0*x+5.0d0-(bmShift+1.0d0)
    return
  end function fQuadDerivative

  subroutine fQuadBoth(x,f,df)
    !!{
    Combined value and derivative of \mono{fQuad}.
    !!}
    double precision, intent(in   ) :: x
    double precision, intent(  out) :: f, df

    f =(x-(bmShift+1.0d0))*(x+5.0d0)
    df=2.0d0*x+5.0d0-(bmShift+1.0d0)
    return
  end subroutine fQuadBoth

  double precision function fFar(x)
    !!{
    Linear function with its root far above a typical initial guess.

    $f(x) = x - (\mathrm{bmShift}+10)$. With $\mathrm{bmShift} \in [-0.4,0.4]$
    the root lies in $[9.6,10.4]$. Used in expansion benchmarks where the
    initial bracket/guess is intentionally below the root so the bracket
    expansion loop is exercised.
    !!}
    double precision, intent(in   ) :: x

    fFar=x-(bmShift+10.0d0)
    return
  end function fFar

end module Benchmark_Root_Finder_Functions
