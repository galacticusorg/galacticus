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
Contains a module which imports the FFTW3 library Fortran interface.
!!}

module FFTW3
  !!{
  Imports the FFTW3 library Fortran interface.
  !!}
#ifdef FFTW3AVAIL
  use, intrinsic :: ISO_C_Binding, only : c_ptr           , c_int     , c_funptr, c_float        , &
       &                                  c_double        , c_int32_t , c_size_t, c_float_complex, &
       &                                  c_double_complex, c_intptr_t, c_char
#endif
  private
  public  :: FFTW_Wavenumber , FFTW_Forward     , FFTW_Estimate, FFTW_Plan_DFT_3D, &
       &     FFTW_Execute_DFT, FFTW_Destroy_Plan
#ifdef FFTW3AVAIL
  include 'fftw3.f03'
#endif

contains

  double precision function FFTW_Wavenumber(k,n)
    !!{
    Return the wavenumber (in units of $1/L$ where $L$ is the box length) corresponding to element {\normalfont \ttfamily k} out of {\normalfont \ttfamily n} of a
    1-D FFT using the FFTW convention.
    !!}
    implicit none
    integer, intent(in   ) :: k , n
    integer                :: kk

    kk=k-1
    if (kk < n/2) then
       FFTW_Wavenumber=dble(kk  )
    else
       FFTW_Wavenumber=dble(kk-n)
    end if
    return
  end function FFTW_Wavenumber

end module FFTW3
