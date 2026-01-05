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
Contains a module which manages HDF5 output from \glc.
!!}

module Output_HDF5
  !!{
  Manages HDF5 output from \glc.
  !!}
  use :: HDF5   , only : HSIZE_T
  use :: IO_HDF5, only : hdf5Object
  implicit none
  public

  ! Flag indicating if output file has been opened.
  logical                            :: outputFileIsOpen      =.false.

  ! Galacticus output file object.
  type   (hdf5Object  ), allocatable :: outputFile                   , outputGroup

  ! Chunk size.
  integer(kind=HSIZE_T)              :: hdf5ChunkSize         =+1

  ! Compression level (-1 means no compression, 0-9 means GNU gzip compression with higher numbers giving more compression).
  integer                            :: hdf5CompressionLevel  =-1

  ! Sieve buffer size.
  integer(kind=HSIZE_T)              :: hdf5SieveBufferSize

  ! File format to use.
  logical                            :: hdf5UseLatestFormat

  ! Cache size.
  integer(kind=size_t )              :: hdf5CacheElementsCount       , hdf5CacheSizeBytes

end module Output_HDF5
