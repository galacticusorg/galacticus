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
Contains a module that implements a simple two-dimensional integrator.
!!}

module Numerical_Integration_2D
  !!{
  Implements a simple two-dimensional integrator.
  !!}
  use :: Numerical_Integration, only : integrator
  implicit none
  private
  public :: integrator2D

  type :: integrator2D
     !!{
     A simple 2D integrator. The integrand and tolerances are fixed at construction, and the (one-dimensional) integrator
     objects used for the inner and outer integrals are built once in the constructor and reused on each call to
     {\normalfont \ttfamily integrate}.
     !!}
     private
     procedure       (integrandTemplate), nopass        , pointer :: integrand   => null()
     double precision                   , dimension(2,2)          :: boundaries
     type            (integrator        )                         :: integratorX           , integratorY
   contains
     !![
     <methods>
       <method method="integrate" description="Perform integration over a 2D region."/>
     </methods>
     !!]
     procedure :: integrate => integrator2DIntegrate
  end type integrator2D

  interface integrator2D
     !!{
     Constructors for the \refClass{integrator2D} class.
     !!}
     module procedure :: integrator2DConstructor
  end interface integrator2D

  abstract interface
     !!{
     Template for 2D integrands.
     !!}
     double precision function integrandTemplate(x,y)
       double precision, intent(in) :: x, y
     end function integrandTemplate
  end interface

  ! Thread-private state used to communicate the active object and the current value of the outer integration variable to the
  ! module-level integrand functions. The underlying 1D integrator invokes its integrand through a plain procedure pointer which
  ! can carry no object state, so this state is passed here instead. As a consequence a given integrator2D object must not be
  ! shared between threads (each thread should construct its own).
  class           (integrator2D), pointer :: self_ => null()
  double precision                        :: x_
  !$omp threadprivate(self_,x_)

contains

  function integrator2DConstructor(integrand,toleranceRelative,toleranceAbsolute) result(self)
    !!{
    Constructor for the \refClass{integrator2D} class. The {\normalfont \ttfamily integrand} is the two-dimensional function to
    be integrated. The (optional) {\normalfont \ttfamily toleranceRelative} and {\normalfont \ttfamily toleranceAbsolute}
    arguments set the relative and absolute tolerances used for both the inner and outer integrals.
    !!}
    implicit none
    type            (integrator2D     )                          :: self
    procedure       (integrandTemplate)                          :: integrand
    double precision                   , intent(in   ), optional :: toleranceRelative, toleranceAbsolute
    !![
    <optionalArgument name="toleranceRelative" defaultsTo="1.0d-2"/>
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0" />
    !!]

    self%integrand   => integrand
    self%integratorX =  integrator(integrator2DIntegrandX,toleranceRelative=toleranceRelative_,toleranceAbsolute=toleranceAbsolute_)
    self%integratorY =  integrator(integrator2DIntegrandY,toleranceRelative=toleranceRelative_,toleranceAbsolute=toleranceAbsolute_)
    return
  end function integrator2DConstructor

  double precision function integrator2DIntegrate(self,boundaries)
    !!{
    Perform integration over a 2D region. {\normalfont \ttfamily boundaries} holds the integration limits, with
    {\normalfont \ttfamily boundaries(1,:)} the limits of the outer ($x$) integral and {\normalfont \ttfamily boundaries(2,:)}
    the limits of the inner ($y$) integral.
    !!}
    implicit none
    class           (integrator2D), intent(inout), target        :: self
    double precision              , intent(in   ), dimension(2,2) :: boundaries

    self%boundaries      =  boundaries
    self_                => self
    integrator2DIntegrate=  self%integratorX%integrate(boundaries(1,1),boundaries(1,2))
    return
  end function integrator2DIntegrate

  double precision function integrator2DIntegrandX(x)
    !!{
    Evaluate the inner ($y$) integral at fixed $x$.
    !!}
    implicit none
    double precision, intent(in) :: x

    x_                    =x
    integrator2DIntegrandX=self_%integratorY%integrate(self_%boundaries(2,1),self_%boundaries(2,2))
    return
  end function integrator2DIntegrandX

  double precision function integrator2DIntegrandY(y)
    !!{
    Evaluate the 2D integrand at the current outer-integration $x$ and the given $y$.
    !!}
    implicit none
    double precision, intent(in) :: y

    integrator2DIntegrandY=self_%integrand(x_,y)
    return
  end function integrator2DIntegrandY

end module Numerical_Integration_2D
