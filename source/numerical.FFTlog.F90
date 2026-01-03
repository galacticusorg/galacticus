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
Contains a module which wraps the \gls{FFTLog} functions.
!!}

module FFTLogs
  !!{
  Wraps the \gls{FFTLog} functions.
  !!}
  private
  public :: FFTLog, FFTLogSineTransform, FFTLogCosineTransform

  !: $(BUILDPATH)/FFTlog/cdgamma.o
  !: $(BUILDPATH)/FFTlog/drfftb.o
  !: $(BUILDPATH)/FFTlog/drfftf.o
  !: $(BUILDPATH)/FFTlog/drffti.o
  !: $(BUILDPATH)/FFTlog/fftlog.o

  ! Labels for forward/backward FFTs.
  integer         , parameter, public :: fftLogForward =+1
  integer         , parameter, public :: fftLogBackward=-1

  ! Values of mu for sine/cosine transforms.
  double precision, parameter, public :: fftLogSine    =+0.5d0
  double precision, parameter, public :: fftLogCosine  =-0.5d0

contains

  subroutine FFTLogSineTransform(r,k,f,ft,direction)
    !!{
    Wrapper function for \hyperlink{numerical.FFTlog.F90:fftlogs:fftlog}{{\normalfont \ttfamily FFTLog()}} which performs a
    Fourier sine transform. Since \hyperlink{numerical.FFTlog.F90:fftlogs:fftlog}{{\normalfont \ttfamily FFTLog()}} achieves
    this by using the $J_{1/2}(x)=(2/\pi x)^{1/2} \sin(x)$ Bessel function we apply the inverse of these factors to get a sine
    transform.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ), dimension(     : ) :: r        , f
    double precision, intent(inout), dimension(     : ) :: k
    double precision, intent(  out), dimension(     : ) :: ft
    integer         , intent(in   )                     :: direction
    double precision               , dimension(size(f)) :: fScaled

    fScaled=f*sqrt(r)
    call FFTLog(r,k,fScaled,ft,fftLogSine,direction)
    ft=ft*sqrt(k)*sqrt(Pi/2.0d0)
    return
  end subroutine FFTLogSineTransform
  
  subroutine FFTLogCosineTransform(r,k,f,ft,direction)
    !!{
    Wrapper function for \hyperlink{numerical.FFTlog.F90:fftlogs:fftlog}{{\normalfont \ttfamily FFTLog()}} which performs a
    Fourier cosine transform. Since \hyperlink{numerical.FFTlog.F90:fftlogs:fftlog}{{\normalfont \ttfamily FFTLog()}} achieves
    this by using the $J_{1/2}(x)=(2/\pi x)^{1/2} \cos(x)$ Bessel function we apply the inverse of these factors to get a
    cosine transform.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ), dimension(     : ) :: r        , f
    double precision, intent(inout), dimension(     : ) :: k
    double precision, intent(  out), dimension(     : ) :: ft
    integer         , intent(in   )                     :: direction
    double precision               , dimension(size(f)) :: fScaled

    fScaled=f*sqrt(r)
    call FFTLog(r,k,fScaled,ft,fftLogCosine,direction)
    ft=ft*sqrt(k)*sqrt(Pi/2.0d0)
    return
  end subroutine FFTLogCosineTransform
  
  subroutine FFTLog(r,k,f,ft,mu,direction)
    !!{
    Perform a discrete FFT on logarithmically spaced data.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   ), dimension(:                          ) :: r, f
    double precision, intent(inout), dimension(:                          ) :: k
    double precision, intent(  out), dimension(:                          ) :: ft
    double precision, intent(in   )                                         :: mu
    integer         , intent(in   )                                         :: direction
    double precision, parameter                                             :: bias    =0.0d0
    integer         , parameter                                             :: krOption=1     ! Ajdust krCentral to closest low-ringing value.
    double precision                , dimension(2*size(r)+3*(size(r)/2)+19) :: workSpace
    double precision                                                        :: deltaLogR , krCentral, kCentralLogarithmic, iCentral, rCentralLogarithmic, normalization
    integer                                                                 :: i
    logical                                                                 :: errorCode

    ! Compute point separation.
    deltaLogR=log(r(size(r))/r(1))/dble(size(r))
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
    normalization=exp(2.0d0*rCentralLogarithmic)
    call fftl(size(r),ft,normalization,direction,workSpace)
    ! Compute actual k values.
    forall(i=1:size(r))
       k(i)=exp(kCentralLogarithmic+(dble(i)-iCentral)*deltaLogR)
    end forall
  return
  end subroutine FFTLog

end module FFTLogs
