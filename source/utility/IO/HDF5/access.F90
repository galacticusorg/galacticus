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
Contains a module which provides an OpenMP lock to serialize access to the HDF5 library, which is not fully thread-safe in all configurations.
!!}

module HDF5_Access
  !!{RST
  Provides an OpenMP lock to serialize access to the HDF5 library, preventing race conditions when multiple threads attempt concurrent HDF5 operations.
  !!}
  use :: Locks, only : mutex
  implicit none
  public

  ! Lock object to coordinate access to HDF5.
  type   (mutex) :: hdf5Access
  logical        :: hdf5AccessInitialized=.false.
  
end module HDF5_Access
