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
  implicit none
  private
  public :: integrator2D

  type :: integrator2D
     !!{
     A simple 2D integrator.
     !!}
     private
     procedure       (integrandTemplate), nopass        , pointer :: integrand
     double precision                   , dimension(2,2)          :: boundaries
   contains
     !![
     <methods>
       <method method="integrate"    description="Perform integration over a 2D region."/>
       <method method="setIntegrand" description="The the integrand for 2D integration."/>
     </methods>
     !!]
     procedure :: integrate    => integrator2DIntegrate
     procedure :: setIntegrand => integrator2DSetIntegrand
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
    double precision function integrandTemplate(x, y)
      double precision, intent(in) :: x, y
    end function integrandTemplate
  end interface

contains

  function integrator2DConstructor() result(self)
    !!{
    Constructor for the \refClass{integrator2D} class.
    !!}
    implicit none
    type(integrator2D) :: self

    self%integrand => null()
    return
  end function integrator2DConstructor
  
  subroutine integrator2DSetIntegrand(self,integrand)
    !!{
    Set the integrand for 2D integration.
    !!}
    class    (integrator2D     ), intent(inout) :: self
    procedure(integrandTemplate)                :: integrand
    
    self%integrand => integrand
    return
  end subroutine integrator2DSetIntegrand

  double precision function integrator2DIntegrate(self,boundaries)
    !!{
    Perform integration over a 2D region.
    !!}
    use :: Numerical_Integration, only: integrator
    use :: Coordinates          , only: coordinateCartesian
    implicit none
    class           (integrator2D), intent(inout), target         :: self
    double precision              , intent(in   ), dimension(2,2) :: boundaries
    type            (integrator  )                                :: integrator_
    double precision                                              :: x_         , y_

    self%boundaries      =boundaries
    integrator_          =integrator(integrator2DIntegrandX,toleranceRelative=1.0d-2)
    integrator2DIntegrate=integrator_%integrate(boundaries(1,1),boundaries(1,2))
    return

  contains

    double precision function integrator2DIntegrandX(x)
      !!{
      Perform the integrand integrated over $y$.
      !!}
      implicit none
      double precision            , intent(in) :: x
      type            (integrator)             :: integrator_

      x_                    =x
      integrator_           =integrator           (integrator2DIntegrandY,toleranceRelative=1.0d-2)
      integrator2DIntegrandX=integrator_%integrate(self%boundaries(2,1), self%boundaries(2,2))
      return
    end function integrator2DIntegrandX

    double precision function integrator2DIntegrandY(y)
      !!{
      Evaluate the 2D integrand.
      !!}
      implicit none
      double precision, intent(in) :: y
      
      y_                    =y
      integrator2DIntegrandY=self%integrand(x_,y_)
      return
    end function integrator2DIntegrandY
    
  end function integrator2DIntegrate

end module Numerical_Integration_2D
