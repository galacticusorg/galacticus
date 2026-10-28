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
Contains a program to test the :term:`FFTLog` wrapper functions.
!!}

program Test_FFTLog
  !!{RST
  Tests the :term:`FFTLog` wrapper functions against Hankel and Fourier sine transforms of Gaussians, which are known
  analytically. For :math:`\mu > -1`,

  .. math::

     \int_0^\infty r^{\mu+1} \mathrm{e}^{-r^2/2} J_\mu(kr) \mathrm{d}r = k^\mu \mathrm{e}^{-k^2/2},

  so, since ``FFTLog()`` computes :math:`\tilde{f}(k) = \int_0^\infty f(r) (kr)^{1/2} J_\mu(kr) \mathrm{d}r`, the transform of
  :math:`f(r) = r^{\mu+1/2} \mathrm{e}^{-r^2/2}` is :math:`k^{\mu+1/2} \mathrm{e}^{-k^2/2}`. Similarly, the Fourier sine
  transform of :math:`r \mathrm{e}^{-r^2/2}` is :math:`\int_0^\infty r \mathrm{e}^{-r^2/2} \sin(kr) \mathrm{d}r =
  (\pi/2)^{1/2} k \mathrm{e}^{-k^2/2}`.
  !!}
  use :: Display                 , only : displayVerbositySet, verbosityLevelStandard
  use :: FFTLogs                 , only : FFTLog             , FFTLogSineTransform   , fftLogBackward      , fftLogForward
  use :: Numerical_Constants_Math, only : Pi
  use :: Numerical_Ranges        , only : Make_Range         , rangeTypeLogarithmic
  use :: Unit_Tests              , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                         , parameter              :: countPoints      =2048
  ! The grid is deliberately not centered on unity, so that the normalization of the transforms (which depends on the
  ! central radius and wavenumber) is tested. It extends to small radii so that the functions transformed are small at both ends
  ! of the grid (as FFTLog treats them as periodic).
  double precision                , parameter              :: radiusMinimum    =1.0d-7             , radiusMaximum    =1.0d+9, &
       &                                                      wavenumberMinimum=1.0d-2             , wavenumberMaximum=3.0d+0
  double precision                , dimension(3          ) :: orders           =[0.0d0,1.0d0,1.5d0]
  double precision                , dimension(countPoints) :: radii                                , wavenumbers             , &
       &                                                      function_                            , transform               , &
       &                                                      roundTrip                            , radiiRoundTrip
  logical                         , dimension(countPoints) :: inRange
  integer                                                  :: i
  character       (len=16        )                         :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("FFTLog")
  radii=Make_Range(radiusMinimum,radiusMaximum,countPoints,rangeTypeLogarithmic)
  ! Hankel transforms of Gaussians.
  do i=1,size(orders)
     write (label,'(f3.1)') orders(i)
     function_=radii**(orders(i)+0.5d0)*exp(-0.5d0*radii**2)
     call FFTLog(radii,wavenumbers,function_,transform,orders(i),fftLogForward)
     inRange  =wavenumbers >= wavenumberMinimum .and. wavenumbers <= wavenumberMaximum
     call Assert(                                                                         &
          &      'Hankel transform of a Gaussian, μ='//trim(label)                      , &
          &      pack(transform                                                ,inRange), &
          &      pack(wavenumbers**(orders(i)+0.5d0)*exp(-0.5d0*wavenumbers**2),inRange), &
          &      relTol=1.0d-6                                                            &
          &     )
     ! A backward transform should recover the original function.
     call FFTLog(wavenumbers,radiiRoundTrip,transform,roundTrip,orders(i),fftLogBackward)
     call Assert(                                                                                              &
          &      'Hankel transform of a Gaussian, μ='//trim(label)//': transformed radii are the input radii', &
          &      radiiRoundTrip                                                                              , &
          &      radii                                                                                       , &
          &      relTol=1.0d-9                                                                                 &
          &     )
     inRange  =radii >= wavenumberMinimum .and. radii <= wavenumberMaximum
     call Assert(                                                                                        &
          &      'Hankel transform of a Gaussian, μ='//trim(label)//': forward then backward transform', &
          &      pack(roundTrip,inRange)                                                               , &
          &      pack(function_,inRange)                                                               , &
          &      relTol=1.0d-6                                                                           &
          &     )
  end do
  ! Fourier sine transform of a Gaussian.
  function_=radii*exp(-0.5d0*radii**2)
  call FFTLogSineTransform(radii,wavenumbers,function_,transform,fftLogForward)
  inRange  =wavenumbers >= wavenumberMinimum .and. wavenumbers <= wavenumberMaximum
  call Assert(                                                                     &
       &      'Fourier sine transform of a Gaussian'                             , &
       &      pack(transform                                            ,inRange), &
       &      pack(sqrt(0.5d0*Pi)*wavenumbers*exp(-0.5d0*wavenumbers**2),inRange), &
       &      relTol=1.0d-6                                                        &
       &     )
  ! A backward sine transform should recover the original function.
  call FFTLogSineTransform(wavenumbers,radiiRoundTrip,transform,roundTrip,fftLogBackward)
  inRange  =radii >= wavenumberMinimum .and. radii <= wavenumberMaximum
  call Assert(                                                                         &
       &      'Fourier sine transform of a Gaussian: forward then backward transform', &
       &      pack(roundTrip,inRange)                                                , &
       &      pack(function_,inRange)                                                , &
       &      relTol=1.0d-6                                                            &
       &     )
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_FFTLog
