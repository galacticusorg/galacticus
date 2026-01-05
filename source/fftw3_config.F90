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

! Exclude from "make all".
!/ exclude

! Contains a test code used to see if FFTW3 is available.

program fftw3_config
  ! Test code used to see if FFTW3 is available.
  use, intrinsic :: ISO_C_Binding, only : c_ptr           , c_int     , c_funptr, c_float        , &
       &                                  c_double        , c_int32_t , c_size_t, c_float_complex, &
       &                                  c_double_complex, c_intptr_t, c_char
  include 'fftw3.f03'
end program fftw3_config
