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

!!{RST
Contains a module which implements calculations of error functions.
!!}

! Add dependency on GSL library.
!; gsl

module Error_Functions
  !!{RST
  Implements calculations of error functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Error_Function           , Error_Function_Complementary  , erfApproximate         , Faddeeva, &
       &    Error_Function_Difference, Error_Function_Difference_Fast, Error_Function_Tabulate

  ! Tabulation of the error function for fast (approximate) evaluation. erfTable(i)=erf(i·step) on a uniform
  ! grid spanning x∈[0,erfTableMaximum], with one ghost point beyond each end so the cubic interpolation
  ! stencil is always valid; erf(-x)=-erf(x) is used for negative arguments and erf saturates to ±1 beyond
  ! the table. The table is evaluated with a C¹-continuous Catmull-Rom cubic, giving ~10⁻¹² accuracy —
  ! comfortably below the absolute tolerances of the numerical integrands that use it, so those integrands
  ! remain smooth enough (and accurate enough) for the adaptive integrator to converge. Built once by
  ! Error_Function_Tabulate(); written once (serially, before any parallel use) and thereafter only read.
  integer         , parameter                         :: erfTableIntervals=3000
  double precision, parameter                         :: erfTableMaximum  =6.0d0
  double precision, parameter                         :: erfTableStep     =erfTableMaximum/dble(erfTableIntervals)
  double precision, dimension(-1:erfTableIntervals+1) :: erfTable
  logical                                             :: erfTabulated     =.false.

  interface Error_Function
     module procedure Error_Function_Real
     module procedure Error_Function_Quad
     module procedure Error_Function_Complex
  end interface Error_Function

  interface Error_Function_Complementary
     module procedure Error_Function_Complementary_Real
     module procedure Error_Function_Complementary_Quad
     module procedure Error_Function_Complementary_Complex
  end interface Error_Function_Complementary

  interface erfApproximate
     module procedure erfApproximateQuad
  end interface erfApproximate
  
contains

  elemental double precision function Error_Function_Real(argument)
    !!{RST
    Computes the error function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function_Real=erf(argument)
    return
  end function Error_Function_Real

  elemental function Error_Function_Quad(argument)
    !!{RST
    Computes the error function with quad precision.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real(kind=kind_quad)                :: Error_Function_Quad
    real(kind=kind_quad), intent(in   ) :: argument

    Error_Function_Quad=erf(argument)
    return
  end function Error_Function_Quad

  elemental double precision function Error_Function_Complementary_Real(argument)
    !!{RST
    Computes the complementary error function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function_Complementary_Real=erfc(argument)
    return
  end function Error_Function_Complementary_Real

  elemental function Error_Function_Complementary_Quad(argument)
    !!{RST
    Computes the complementary error function.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real(kind=kind_quad)                :: Error_Function_Complementary_Quad
    real(kind=kind_quad), intent(in   ) :: argument

    Error_Function_Complementary_Quad=erfc(argument)
    return
  end function Error_Function_Complementary_Quad

  elemental double complex function Error_Function_Complex(argument)
    !!{RST
    Computes the complex complementary error function.
    !!}
    implicit none
    double complex, intent(in   ) :: argument

    Error_Function_Complex=1.0d0-Error_Function_Complementary_Complex(argument)
    return
  end function Error_Function_Complex

  elemental double complex function Error_Function_Complementary_Complex(argument)
    !!{RST
    Computes the complex complementary error function, using the algorithm of :cite:t:`abrarov_rapid_2013`.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double complex  , intent(in   ) :: argument
    integer         , parameter     :: N        =23
    double precision, parameter     :: sigma    = 2.0d0
    double precision, parameter     :: tauM     =12.0d0
    integer                         :: i
    double complex                  :: summation       , z
    double precision                :: A               , B

    if (real(argument) < 0.0d0) then
       z=-argument
    else
       z=+argument
    end if
    summation=dcmplx(0.0d0,0.0d0)
    do i=1,N
       A=                      &
            & +2.0d0           &
            & *tauM            &
            & *exp(            &
            &      +sigma  **2 &
            &      -dble(i)**2 &
            &      *Pi     **2 &
            &      /tauM   **2 &
            &     )            &
            & *cos(            &
            &      +2.0d0      &
            &      *dble(i)    &
            &      *Pi         &
            &      *sigma      &
            &      /tauM       &
            &     )
       B=                      &
            & +2.0d0           &
            & *dble(i)         &
            & *Pi              &
            & *exp(            &
            &      +sigma  **2 &
            &      -dble(i)**2 &
            &      *Pi     **2 &
            &      /tauM   **2 &
            &     )            &
            & *sin(            &
            &      +2.0d0      &
            &      *dble(i)    &
            &      *Pi         &
            &      *sigma      &
            &      /tauM       &
            &     )
       summation=                            &
            & +summation                     &
            & +(                             &
            &   +A                           &
            &   *(                           &
            &     +sigma                     &
            &     +z                         &
            &    )                           &
            &   +B                           &
            &  )                             &
            & /(                             &
            &   +dcmplx(dble(i)*Pi,0.0d0)**2 &
            &   +(                           &
            &     +tauM                      &
            &     *(                         &
            &       +sigma                   &
            &       +z                       &
            &      )                         &
            &    )                       **2 &
            &  )
    end do
    Error_Function_Complementary_Complex= &
         & +exp(-z**2)                    &
         & *(                             &
         &   +exp(sigma**2)               &
         &   /tauM                        &
         &   /(                           &
         &     +sigma                     &
         &     +z                         &
         &    )                           &
         &   +summation                   &
         &  )
    if (real(argument) < 0.0d0) Error_Function_Complementary_Complex=2.0d0-Error_Function_Complementary_Complex
    return
  end function Error_Function_Complementary_Complex

  elemental double complex function Faddeeva(argument)
    !!{RST
    The `Fadeeva function <http://en.wikipedia.org/wiki/Faddeeva_function>`_.
    !!}
    implicit none
    double complex  , intent(in   ) :: argument

    Faddeeva=exp(-argument**2)*Error_Function_Complementary(-dcmplx(0.0d0,1.0d0)*argument)
    return
  end function Faddeeva

 function erfApproximateQuad(x)
    !!{RST
    An approximation to the error function due to :cite:t:`winitzki_uniform_2003`. that is designed to be very accurate in the vicinity of zero and infinity.
    !!}
    use :: Kind_Numbers            , only : kind_quad
    use :: Numerical_Constants_Math, only : PiQuadPrecision
    implicit none
    real(kind=kind_quad)                :: erfApproximateQuad
    real(kind=kind_quad), intent(in   ) :: x
    real(kind=kind_quad), parameter     :: a               =8.00_kind_quad*(PiQuadPrecision-3.0_kind_quad)/3.0_kind_quad/PiQuadPrecision/(4.0_kind_quad-PiQuadPrecision)
    ! Value above which erf(x)=1 to quad precision.
    real(kind=kind_quad), parameter     :: xMaximum        =8.75_kind_quad
    
    if (x > xMaximum) then
       erfApproximateQuad=1.0_kind_quad
    else
       erfApproximateQuad=sqrt(1.0_kind_quad-exp(-x**2*(4.0_kind_quad/PiQuadPrecision+a*x**2)/(1.0_kind_quad+a*x**2)))
    end if
    return
  end function erfApproximateQuad

  double precision function Error_Function_Difference(x1,x2) result(difference)
    !!{RST
    Evaluates the difference in the error function at the given arguments, :math:`\mathrm{erf}(x_2)-\mathrm{erf}(x_1)`. Utiliizes symmetries of the error function and the complementary error function to maintain accuracy.
    !!}
    use :: Error       , only : Error_Report
    use :: Kind_Numbers, only : kind_quad   , kind_dble
    implicit none
    double precision, intent(in   ) :: x1              , x2
    double precision, parameter     :: tolerance=1.0d-6
    double precision                :: erf1            , erf2

    ! Validate input.
    difference=0.0d0
    if (x2 == x1) return
    if (x2 <  x1) call Error_Report('x₂ ≥ x₁ is required'//{introspection:location})
    ! First try simply evaluating the difference.
    erf1      =+erf(x1)
    erf2      =+erf(x2)
    difference=+erf2    &
         &     -erf1
    ! Check if a non-zero result was found.
    if (difference > 0.0d0) return
    ! Check if we are close to +∞.
    if (erf2 > +1.0d0-tolerance) then
       erf2      =+erfc(x2)
       erf1      =+erfc(x1)
       difference=-erf2 &
            &     +erf1
       if (difference > 0.0d0) return
       ! Try using quad precision.
       difference=real(                                                &
            &               +erfApproximate(real(+x2,kind=kind_quad))  &
            &               -erfApproximate(real(+x1,kind=kind_quad)), &
            &          kind=kind_dble                                  &
            &         )
       if (difference > 0.0d0) return
    end if
    ! Check if we are close to -∞.
    if (erf1 < -1.0d0+tolerance) then
       erf1      =+erfc(-x1)
       erf2      =+erfc(-x2)
       difference=+erf2      &
            &     -erf1
       if (difference > 0.0d0) return

       ! Try using quad precision.
       difference=real(                                                &
            &               -erfApproximate(real(-x2,kind=kind_quad))  &
            &               +erfApproximate(real(-x1,kind=kind_quad)), &
            &          kind=kind_dble                                  &
            &         )
       if (difference > 0.0d0) return
    end if
    ! Nothing has worked - return 0.
    return
  end function Error_Function_Difference

  subroutine Error_Function_Tabulate()
    !!{RST
    Build the tabulation of the error function used by {\normalfont \ttfamily Error\_Function\_Difference\_Fast}
    (and {\normalfont \ttfamily erfFast}). This must be called once before any fast evaluation. It is intended
    to be called from a serial initialization context (e.g.~an object constructor), so no locking is performed;
    the resulting table is written once and thereafter only read, so it is safe to use from parallel regions
    entered after this routine has completed.
    !!}
    implicit none
    integer :: i

    if (erfTabulated) return
    do i=-1,erfTableIntervals+1
       erfTable(i)=erf(dble(i)*erfTableStep)
    end do
    erfTabulated=.true.
    return
  end subroutine Error_Function_Tabulate

  elemental double precision function erfFast(x) result(erfValue)
    !!{RST
    A fast, tabulated (approximate) evaluation of the error function, :math:`\mathrm{erf}(x)`, accurate to
    :math:`\sim 10^{-12}` for :math:`|x|<` {\normalfont \ttfamily erfTableMaximum} and saturating to
    :math:`\pm 1` beyond. Uses C¹-continuous Catmull-Rom cubic interpolation of the tabulation so that
    integrands built from it remain smooth. {\normalfont \ttfamily Error\_Function\_Tabulate()} must have
    been called first.
    !!}
    implicit none
    double precision, intent(in   ) :: x
    double precision                :: xAbs, t , &
         &                             p0  , p1, &
         &                             p2  , p3
    integer                         :: i

    xAbs=abs(x)
    if (xAbs >= erfTableMaximum) then
       erfValue=1.0d0
    else
       t =xAbs/erfTableStep
       i =int(t)                                  ! interval index, in [0,erfTableIntervals-1]
       t =t-dble(i)                               ! offset within the interval, in [0,1)
       p0=erfTable(i-1)
       p1=erfTable(i  )
       p2=erfTable(i+1)
       p3=erfTable(i+2)
       ! Catmull-Rom cubic through (p0,p1,p2,p3), evaluated at t between p1 and p2.
       erfValue=p1+0.5d0*t*(                                    &
            &               (p2-p0)                             &
            &               +t*(                                &
            &                   (2.0d0*p0-5.0d0*p1+4.0d0*p2-p3) &
            &                   +t*(3.0d0*(p1-p2)+p3-p0)        &
            &                  )                                &
            &              )
    end if
    if (x < 0.0d0) erfValue=-erfValue
    return
  end function erfFast

  double precision function Error_Function_Difference_Fast(x1,x2) result(difference)
    !!{RST
    A fast (approximate) evaluation of :math:`\mathrm{erf}(x_2)-\mathrm{erf}(x_1)` intended for use inside
    numerical integrands where full precision is unnecessary. When both arguments lie within the tabulated
    range the tabulated error function is used; if either argument lies outside that range---the deep tails,
    where the direct difference suffers catastrophic cancellation and the exact treatment is required---the
    exact {\normalfont \ttfamily Error\_Function\_Difference} is used instead. {\normalfont \ttfamily
    Error\_Function\_Tabulate()} must have been called first. As for {\normalfont \ttfamily
    Error\_Function\_Difference}, :math:`x_2 \ge x_1` is required.
    !!}
    implicit none
    double precision, intent(in   ) :: x1, x2

    if (max(abs(x1),abs(x2)) < erfTableMaximum) then
       difference=erfFast(x2)-erfFast(x1)
    else
       difference=Error_Function_Difference(x1,x2)
    end if
    return
  end function Error_Function_Difference_Fast

end module Error_Functions
