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

!!{
Contains a program which determines the C-interoperable Fortran types that correspond to the HDF5 \mono{hid\_t},
\mono{hsize\_t}, and \mono{size\_t} types on the current system.
!!}

program hdf5FCInterop
  !!{
  Determines C-interoperable Fortran types that correspond to HDF5 types. The results are written to standard output and are
  used during the build process to select appropriate \mono{ISO\_C\_Binding} type constants for each HDF5 type, avoiding
  potential compiler warnings about C interoperability that would arise from using the HDF5-provided type kinds directly.
  !!}
  use            :: HDF5         , only : hid_t, hsize_t, size_t
  use, intrinsic :: ISO_C_Binding, only : c_int, c_long , c_long_long, c_size_t
  implicit none
  integer(hid_t      ) :: type_hid_t
  integer(hsize_t    ) :: type_hsize_t
  integer(size_t     ) :: type_size_t
  integer(c_int      ) :: type_c_int
  integer(c_size_t   ) :: type_c_size_t
  integer(c_long     ) :: type_c_long
  integer(c_long_long) :: type_c_long_long

  ! hid_t type.
  if      (sizeof(type_hid_t) == sizeof(type_c_int      )) then
     write (*,*) "hid_t = c_int"
  else if (sizeof(type_hid_t) == sizeof(type_c_size_t   )) then
     write (*,*) "hid_t = c_size_t"
  else if (sizeof(type_hid_t) == sizeof(type_c_long     )) then
     write (*,*) "hid_t = c_long"
  else if (sizeof(type_hid_t) == sizeof(type_c_long_long)) then
     write (*,*) "hid_t = c_long_long"
  end if
  ! hsize_t type.
  if      (sizeof(type_hsize_t) == sizeof(type_c_int      )) then
     write (*,*) "hsize_t = c_int"
  else if (sizeof(type_hsize_t) == sizeof(type_c_size_t   )) then
     write (*,*) "hsize_t = c_size_t"
  else if (sizeof(type_hsize_t) == sizeof(type_c_long     )) then
     write (*,*) "hsize_t = c_long"
  else if (sizeof(type_hsize_t) == sizeof(type_c_long_long)) then
     write (*,*) "hsize_t = c_long_long"
  end if
  ! size_t type.
  if      (sizeof(type_size_t) == sizeof(type_c_int      )) then
     write (*,*) "size_t = c_int"
  else if (sizeof(type_size_t) == sizeof(type_c_size_t   )) then
     write (*,*) "size_t = c_size_t"
  else if (sizeof(type_size_t) == sizeof(type_c_long     )) then
     write (*,*) "size_t = c_long"
  else if (sizeof(type_size_t) == sizeof(type_c_long_long)) then
     write (*,*) "size_t = c_long_long"
  end if
end program hdf5FCInterop

