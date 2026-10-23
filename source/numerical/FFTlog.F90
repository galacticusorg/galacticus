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
Contains a module which wraps the :term:`FFTLog` functions.
!!}

module FFTLogs
  !!{RST
  Wraps the :term:`FFTLog` functions.
  !!}
  private
  public :: FFTLog, FFTLogSineTransform

  !: $(BUILDPATH)/external/FFTlog/cdgamma.o
  !: $(BUILDPATH)/external/FFTlog/drfftb.o
  !: $(BUILDPATH)/external/FFTlog/drfftf.o
  !: $(BUILDPATH)/external/FFTlog/drffti.o
  !: $(BUILDPATH)/external/FFTlog/fftlog.o

  ! Labels for forward/backward FFTs.
  integer         , parameter, public :: fftLogForward =+1
  integer         , parameter, public :: fftLogBackward=-1

  ! Value of mu for sine transforms.
  double precision, parameter, public :: fftLogSine    =+0.5d0

contains

  subroutine FFTLogSineTransform(r,k,f,ft,direction)
    !!{RST
    Wrapper function for ``FFTLog()`` which performs a Fourier sine transform,

    .. math::

       \tilde{f}(k) = \int_0^\infty f(r) \sin(kr) \mathrm{d}r,

    for ``direction=fftLogForward``, or its inverse, :math:`f(r) = (2/\pi) \int_0^\infty \tilde{f}(k) \sin(kr) \mathrm{d}k`,
    for ``direction=fftLogBackward``. Since :math:`(kr)^{1/2} J_{1/2}(kr) = (2/\pi)^{1/2} \sin(kr)`, the transform computed by
    ``FFTLog()`` with :math:`\mu=1/2` is :math:`(2/\pi)^{1/2}` times the sine transform.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ), dimension(     : ) :: r        , f
    double precision, intent(inout), dimension(     : ) :: k
    double precision, intent(  out), dimension(     : ) :: ft
    integer         , intent(in   )                     :: direction

    call FFTLog(r,k,f,ft,fftLogSine,direction)
    select case (direction)
    case (fftLogForward )
       ft=ft*sqrt(Pi   /2.0d0)
    case (fftLogBackward)
       ft=ft*sqrt(2.0d0/Pi   )
    case default
       call Error_Report('direction must be -1 or +1'//{introspection:location})
    end select
    return
  end subroutine FFTLogSineTransform
  
  subroutine FFTLog(r,k,f,ft,mu,direction)
    !!{RST
    Perform a discrete Hankel transform of the function ``f``, tabulated at the logarithmically spaced points ``r``, using
    :term:`FFTLog`. The transform computed is

    .. math::

       \tilde{f}(k) = \int_0^\infty f(r) (kr)^{1/2} J_\mu(kr) \mathrm{d}r,

    which is returned in ``ft`` at the logarithmically spaced points ``k``, which have the same spacing as ``r``. This transform is
    its own inverse, so the forward (``direction=fftLogForward``) and backward (``direction=fftLogBackward``) transforms are
    identical: in either case the input points (``r``) may be radii or wavenumbers, and a forward transform followed by a
    backward transform recovers the original function.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   ), dimension(:                          ) :: r                        , f
    double precision, intent(inout), dimension(:                          ) :: k
    double precision, intent(  out), dimension(:                          ) :: ft
    double precision, intent(in   )                                         :: mu
    integer         , intent(in   )                                         :: direction
    double precision, parameter                                             :: bias               =0.0d0
    integer         , parameter                                             :: krOption           =1 ! Ajdust krCentral to closest low-ringing value.
    double precision                , dimension(2*size(r)+3*(size(r)/2)+19) :: workSpace
    double precision                                                        :: deltaLogR                , krCentral    , &
         &                                                                     kCentralLogarithmic      , iCentral     , &
         &                                                                     rCentralLogarithmic      , normalization
    integer                                                                 :: i
    logical                                                                 :: errorCode

    ! Compute point separation.
    deltaLogR=log(r(size(r))/r(1))/dble(size(r)-1)
    ! Compute central points.
    iCentral =dble(size(r)+1)/2.0d0
    rCentralLogarithmic=log(r(size(r))*r(1))/2.0d0
    krCentral=1.0d0
    ! Call the FFTLog initialization function.
    call fhti(size(r),mu,bias,deltaLogR,krCentral,krOption,workSpace,errorCode)
    if (.not.errorCode) call Error_Report('FFTLog initialization failed'//{introspection:location})
    ! Compute central points.
    kCentralLogarithmic=log(krCentral)-rCentralLogarithmic
    ! Perform the FFT.
    if (direction /= -1 .and. direction /= 1) call Error_Report('direction must be -1 or +1'//{introspection:location})
    ft=f
    ! The normalization required by fftl() is the ratio of the central input point to the central output point, r_c/k_c. Note
    ! that k_c r_c is not necessarily unity, as fhti() adjusts it to a low-ringing value. Since the transform is its own inverse,
    ! fftl() is always called for a forward transform. (With the zero bias used here, a backward transform in fftl() does not
    ! invert its forward transform.)
    normalization=exp(rCentralLogarithmic-kCentralLogarithmic)
    call fftl(size(r),ft,normalization,fftLogForward,workSpace)
    ! Compute actual k values.
    forall(i=1:size(r))
       k(i)=exp(kCentralLogarithmic+(dble(i)-iCentral)*deltaLogR)
    end forall
  return
  end subroutine FFTLog

end module FFTLogs
