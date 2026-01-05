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
Contains a module which provides buffer types for merger tree outputters.
!!}

module Merger_Tree_Outputter_Buffer_Types
  !!{
  Provides buffer types for merger tree outputters.
  !!}
  use :: Kind_Numbers      , only : kind_int8
  use :: Hashes            , only : doubleHash    , rank1DoubleHash
  use :: IO_HDF5           , only : hdf5VarDouble , hdf5VarDouble2D, hdf5VarInteger8
  use :: ISO_Varying_String, only : varying_string
  public

  ! Maximum length of names and comments.
  integer, parameter :: propertyNameLengthMax=256, propertyCommentLengthMax=256
  
  type :: outputPropertyInteger
     !!{
     A type used to store integer data for output.
     !!}
     character       (len=propertyNameLengthMax   )                              :: name
     character       (len=propertyCommentLengthMax)                              :: comment
     type            (doubleHash                  ), allocatable                 :: metaDataRank0
     type            (rank1DoubleHash             ), allocatable                 :: metaDataRank1
     double precision                                                            :: unitsInSI
     integer         (kind_int8                   ), allocatable, dimension(:  ) :: scalar
     integer         (kind_int8                   ), allocatable, dimension(:,:) :: rank1
     type            (hdf5VarInteger8             ), allocatable, dimension(:  ) :: rank1VarLen
     type            (varying_string              ), allocatable, dimension(:  ) :: rank1Descriptors
  end type outputPropertyInteger
  
  type :: outputPropertyDouble
     !!{
     A type used to store double precision data for output.
     !!}
     character       (len=propertyNameLengthMax   )                              :: name
     character       (len=propertyCommentLengthMax)                              :: comment
     type            (doubleHash                  ), allocatable                 :: metaDataRank0
     type            (rank1DoubleHash             ), allocatable                 :: metaDataRank1
     double precision                                                            :: unitsInSI
     double precision                              , allocatable, dimension(:  ) :: scalar
     double precision                              , allocatable, dimension(:,:) :: rank1
     type            (hdf5VarDouble               ), allocatable, dimension(:  ) :: rank1VarLen
     type            (hdf5VarDouble2D             ), allocatable, dimension(:  ) :: rank2VarLen
     type            (varying_string              ), allocatable, dimension(:  ) :: rank1Descriptors
     double precision                              , allocatable, dimension(:  ) :: rank1DescriptorValues
     type            (varying_string               )                             :: rank1DescriptorComment
     double precision                                                            :: rank1DescriptorUnitsInSI
  end type outputPropertyDouble

end module Merger_Tree_Outputter_Buffer_Types
