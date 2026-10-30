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
Contains a module which provides an independent reference calculation for the tidally truncated and heated exponential disk
mass distribution unit tests.
!!}

module Test_Exponential_Disk_Tidal_Reference
  !!{RST
  Provides an independent calculation of the mid-plane rotation curve and potential of a truncated exponential disk, by direct numerical
  evaluation of the Hankel transforms (rather than by the :term:`FFTLog` method used by
  :galacticus-class:`massDistributionExponentialDiskTidal`).
  !!}
  use :: Mass_Distributions_Cylindrical_Utilities, only : kernelVerticalSechSquared
  use :: Numerical_Integration                   , only : integrator
  implicit none
  private
  public :: referenceRotationCurveSquared, referencePotential, referenceRadiusTruncate

  type            (kernelVerticalSechSquared) :: kernel_
  double precision                            :: radiusTruncate_, radiusDecayOverRadiusTruncate_, &
       &                                         wavenumber_    , radius_
  type            (integrator               ) :: integratorInner, integratorOuter               , &
       &                                         integratorPotential

contains

  double precision function surfaceDensity(radius)
    !!{RST
    The surface density of a truncated exponential disk of unit mass and scale length.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: u

    surfaceDensity=exp(-radius)/2.0d0/Pi
    if (radius > radiusTruncate_) then
       u             =radius/radiusTruncate_
       surfaceDensity=surfaceDensity*exp((log(u)-(u-1.0d0))/radiusDecayOverRadiusTruncate_)
    end if
    return
  end function surfaceDensity

  double precision function massRemovedIntegrand(radius)
    !!{RST
    Integrand for the mass removed by truncation.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    massRemovedIntegrand=2.0d0*Pi*radius*(exp(-radius)/2.0d0/Pi-surfaceDensity(radius))
    return
  end function massRemovedIntegrand

  double precision function referenceRadiusTruncate(truncation,radiusDecayOverRadiusTruncate) result(radiusTruncate)
    !!{RST
    Find, by bisection, the truncation radius at which the truncated disk retains a fraction :math:`\exp(-x)` of its mass.
    !!}
    implicit none
    double precision            , intent(in   ) :: truncation        , radiusDecayOverRadiusTruncate
    type            (integrator)                :: integrator_
    double precision                            :: radiusLower       , radiusUpper                  , &
         &                                         massRemovedTarget
    integer                                     :: i

    radiusDecayOverRadiusTruncate_=radiusDecayOverRadiusTruncate
    massRemovedTarget             =1.0d0-exp(-truncation)
    integrator_                   =integrator(massRemovedIntegrand,toleranceAbsolute=1.0d-14,toleranceRelative=1.0d-10)
    radiusLower                   =1.0d-3
    radiusUpper                   =1.0d+2
    do i=1,80
       radiusTruncate_=sqrt(radiusLower*radiusUpper)
       if (integrator_%integrate(radiusTruncate_,radiusTruncate_*(1.0d0+60.0d0*radiusDecayOverRadiusTruncate)+60.0d0) > massRemovedTarget) then
          radiusLower=radiusTruncate_
       else
          radiusUpper=radiusTruncate_
       end if
    end do
    radiusTruncate=radiusTruncate_
    return
  end function referenceRadiusTruncate

  double precision function innerIntegrand(radius)
    !!{RST
    Integrand for the Hankel transform of the surface density.
    !!}
    use :: Bessel_Functions, only : Bessel_Function_J0
    implicit none
    double precision, intent(in   ) :: radius

    innerIntegrand=surfaceDensity(radius)*Bessel_Function_J0(wavenumber_*radius)*radius
    return
  end function innerIntegrand

  double precision function outerIntegrand(wavenumber)
    !!{RST
    Integrand for the mid-plane rotation curve.
    !!}
    use :: Bessel_Functions                        , only : Bessel_Function_J1
    use :: Mass_Distributions_Cylindrical_Utilities, only : Hankel_Transform_Integrate
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: transform

    if (wavenumber <= 0.0d0) then
       outerIntegrand=0.0d0
       return
    end if
    wavenumber_   =wavenumber
    transform     =Hankel_Transform_Integrate(integratorInner,0,wavenumber)
    outerIntegrand=wavenumber*Bessel_Function_J1(wavenumber*radius_)*transform*kernel_%value(wavenumber,1,1)
    return
  end function outerIntegrand

  double precision function potentialIntegrand(wavenumber)
    !!{RST
    Integrand for the mid-plane potential.
    !!}
    use :: Bessel_Functions                        , only : Bessel_Function_J0
    use :: Mass_Distributions_Cylindrical_Utilities, only : Hankel_Transform_Integrate
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: transform

    if (wavenumber <= 0.0d0) then
       potentialIntegrand=0.0d0
       return
    end if
    wavenumber_       =wavenumber
    transform         =Hankel_Transform_Integrate(integratorInner,0,wavenumber)
    potentialIntegrand=Bessel_Function_J0(wavenumber*radius_)*transform*kernel_%value(wavenumber,1,1)
    return
  end function potentialIntegrand

  double precision function referenceRotationCurveSquared(radius,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
    !!{RST
    Return the square of the mid-plane rotation curve (for unit gravitational constant) of a truncated exponential disk of unit
    (untruncated) mass and scale length.
    !!}
    use :: Mass_Distributions_Cylindrical_Utilities, only : Hankel_Transform_Integrate, Kernel_Vertical_Integrators_Allocate, Kernel_Vertical_Integrators_Deallocate
    use :: Numerical_Constants_Math                , only : Pi
    implicit none
    double precision, intent(in   ) :: radius, radiusTruncate, radiusDecayOverRadiusTruncate, scaleHeight

    radius_                       =radius
    radiusTruncate_               =radiusTruncate
    radiusDecayOverRadiusTruncate_=radiusDecayOverRadiusTruncate
    kernel_                       =kernelVerticalSechSquared(scaleHeight,[0.0d0])
    integratorInner               =integrator(innerIntegrand,toleranceAbsolute=1.0d-14,toleranceRelative=1.0d-9)
    integratorOuter               =integrator(outerIntegrand,toleranceAbsolute=1.0d-10,toleranceRelative=1.0d-6)
    call Kernel_Vertical_Integrators_Allocate  ()
    referenceRotationCurveSquared =2.0d0*Pi*radius*Hankel_Transform_Integrate(integratorOuter,1,radius)
    call Kernel_Vertical_Integrators_Deallocate()
    return
  end function referenceRotationCurveSquared

  double precision function referencePotential(radius,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
    !!{RST
    Return the mid-plane potential (for unit gravitational constant) of a truncated exponential disk of unit (untruncated) mass
    and scale length, :math:`\Phi(R,0) = -2\pi \int_0^\infty J_0(kR) \tilde{\Sigma}(k) K(k,0) \mathrm{d}k`.
    !!}
    use :: Mass_Distributions_Cylindrical_Utilities, only : Hankel_Transform_Integrate, Kernel_Vertical_Integrators_Allocate, Kernel_Vertical_Integrators_Deallocate
    use :: Numerical_Constants_Math                , only : Pi
    implicit none
    double precision, intent(in   ) :: radius, radiusTruncate, radiusDecayOverRadiusTruncate, scaleHeight

    radius_                       =radius
    radiusTruncate_               =radiusTruncate
    radiusDecayOverRadiusTruncate_=radiusDecayOverRadiusTruncate
    kernel_                       =kernelVerticalSechSquared(scaleHeight,[0.0d0])
    integratorInner               =integrator(innerIntegrand    ,toleranceAbsolute=1.0d-14,toleranceRelative=1.0d-9)
    integratorPotential           =integrator(potentialIntegrand,toleranceAbsolute=1.0d-10,toleranceRelative=1.0d-6)
    call Kernel_Vertical_Integrators_Allocate  ()
    referencePotential            =-2.0d0*Pi*Hankel_Transform_Integrate(integratorPotential,0,radius)
    call Kernel_Vertical_Integrators_Deallocate()
    return
  end function referencePotential

end module Test_Exponential_Disk_Tidal_Reference
