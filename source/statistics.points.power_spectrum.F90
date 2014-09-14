!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which computes power spectra from point distributions.

module Statistics_Points_Power_Spectra
  !% Compute power spectra from point distributions.
  private
  public :: Statistics_Points_Power_Spectrum

contains

  subroutine Statistics_Points_Power_Spectrum(dataPosition,boxLength,wavenumberMinimum,wavenumberMaximum,wavenumberCount,wavenumber,powerSpectrum)
    !% Compute the power spectrum from a set of points in a periodic cube.
    use, intrinsic :: ISO_C_Binding
    use               Memory_Management
    use               Numerical_Constants_Math
    use               Meshes
    use               Numerical_Ranges
    use               Galacticus_Display
    use               Galacticus_Error
    use               FFTW3
    use               ISO_Varying_String
    use               String_Handling

use vectors

    implicit none
    double precision                  , intent(in   ), dimension(:,:  )              :: dataPosition
    double precision                  , intent(in   )                                :: wavenumberMinimum          , wavenumberMaximum       , &
         &                                                                              boxLength
    double precision                  , intent(  out), dimension(:    ), allocatable :: powerSpectrum              , wavenumber
    integer                           , intent(in   )                                :: wavenumberCount
    complex         (c_double_complex)               , dimension(:,:,:), allocatable :: density                    , densityFourier
    integer                                          , dimension(:    ), allocatable :: powerSpectrumCount
    integer                                                                          :: iBin                       , i                       , &
         &                                                                              gridCount                  , u                       , &
         &                                                                              v                          , w
    complex         (c_double_complex)                                               :: normalization
    type            (c_ptr           )                                               :: plan
    double precision                                                                 :: waveNumberU                , waveNumberV             , &
         &                                                                              waveNumberW                , waveNumberSquared       , &
         &                                                                              logarithmicWavenumberOffset, inverseDeltaLogarithmicWavenumber
    type            (varying_string  )                                               :: message


double precision :: mag, wavemin, wavemax, a
integer :: j

    ! Allocate arrays.
    if (allocated(wavenumber   )) call Dealloc_Array(wavenumber   )
    if (allocated(powerSpectrum)) call Dealloc_Array(powerSpectrum)
    call Alloc_Array(wavenumber        ,[wavenumberCount])
    call Alloc_Array(powerSpectrum     ,[wavenumberCount])
    allocate        (powerSpectrumCount (wavenumberCount))
    ! Generate the array of wavenumbers.
    wavenumber                       =Make_Range(wavenumberMinimum,wavenumberMaximum,wavenumberCount,rangeType=rangeTypeLogarithmic)
    powerSpectrum                    =0.0d0
    powerSpectrumCount               =0
    inverseDeltaLogarithmicWavenumber=1.0d0/(log(wavenumber(2))-log(wavenumber(1)))
    logarithmicWavenumberOffset      =log(wavenumber(1))-0.5d0/inverseDeltaLogarithmicWavenumber
    ! Generate density grid.
    gridCount=2**int(log(wavenumberMaximum*boxLength/sqrt(3.0d0)/Pi)/log(2.0d0)+1.0d0)
    message='Constructing power spectrum using grid of '
    message=message//gridCount//"³ cells"
    call Galacticus_Display_Message(message)
    allocate(density       (gridCount,gridCount,gridCount))
    allocate(densityFourier(gridCount,gridCount,gridCount))
    density=0.0d0
    do i=1,size(dataPosition,dim=2)
       call Meshes_Apply_Point(density,boxLength,dataPosition(:,i),pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypePoint)
    end do
    ! Take the Fourier transform of the selection function.
    plan=fftw_plan_dft_3d(gridCount,gridCount,gridCount,density,densityFourier,FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft (plan,density,densityFourier)
    call fftw_destroy_plan(plan                       )
    ! Normalize the window function.
    normalization=densityFourier(1,1,1)
    if (real(normalization) > 0.0d0) densityFourier=densityFourier/normalization
    ! Accumulate power spectrum to bins. From FFTW, each cell of the Fourier density field has side of length 2π/L.
    !$omp parallel do private(u,v,w,waveNumberU,waveNumberV,waveNumberW,waveNumberSquared,iBin)
    do u=1,gridCount
       waveNumberU      =FFTW_Wavenumber(u,gridCount)*2.0d0*Pi/boxLength
       do v=1,gridCount
          waveNumberV   =FFTW_Wavenumber(v,gridCount)*2.0d0*Pi/boxLength
          do w=1,gridCount
             waveNumberW=FFTW_Wavenumber(w,gridCount)*2.0d0*Pi/boxLength
             ! Compute the net wavenumber.
             waveNumberSquared=waveNumberU**2+waveNumberV**2+waveNumberW**2
             if (waveNumberSquared <= 0.0d0) cycle
             ! Find the bin corresponding to this wavenumber.
             iBin=int((0.5d0*log(waveNumberSquared)-logarithmicWavenumberOffset)*inverseDeltaLogarithmicWavenumber)+1
             ! Accumulate power.
             if (iBin > 0 .and. iBin <= wavenumberCount) then
                !$omp atomic
                powerSpectrum     (iBin)=powerSpectrum     (iBin)+real(densityFourier(u,v,w)*conjg(densityFourier(u,v,w)))
                !$omp atomic
                powerSpectrumCount(iBin)=powerSpectrumCount(iBin)+1
             end if
          end do
       end do
    end do
    !$omp end parallel do
    ! Find average power in each bin, subtract off shot-noise term, and normalize.
    where (powerSpectrumCount > 0)
       powerSpectrum=+(                                &
            &          +powerSpectrum                  &
            &          /dble(powerSpectrumCount      ) &
            &          -1.0d0                          &
            &          /dble(size(dataPosition,dim=2)) &
            &         )                                &
            &        *(                                &
            &          +2.0d0                          &
            &          *Pi                             &
            &         )**3                             &
            &        /(                                &
            &          +2.0d0                          &
            &          *Pi                             &
            &          /boxLength                      &
            &         )**3
    end where

! wavemax=wavenumber(22)
! a=0.0d0
! !$omp parallel do private(i,j,mag) reduction(+:a)
! do i=1,size(dataPosition,dim=2)
! do j=i+1,size(dataPosition,dim=2)
! mag=wavemax*vector_magnitude(dataposition(:,i)-dataposition(:,j))
! a=a+2.0d0*sin(mag)/mag
! enddo
! enddo
! !$omp end parallel do
! a=a/dble(size(dataPosition,dim=2))**2
! write (0,*) wavemax,a*wavemax**3/2.0d0/Pi**2

    return
  end subroutine Statistics_Points_Power_Spectrum

end module Statistics_Points_Power_Spectra
