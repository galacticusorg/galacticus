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

!+    Contributions to this file made by:  Alex Merson.

!!{
Contains a module which holds data for the stellar luminosities class.
!!}

module Stellar_Luminosities_Data
  !!{
  Holds data for the stellar luminosities class.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  implicit none
  public

  ! Used to store a list of output redshifts - this avoids needing to access the outputTimesClass directly from within the stellar
  ! luminosities module which would lead to a circular dependency.
  integer         (c_size_t)                            :: outputCount
  double precision          , dimension(:), allocatable :: outputRedshifts
  
end module Stellar_Luminosities_Data
