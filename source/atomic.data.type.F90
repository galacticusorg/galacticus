!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which defines a data type for storing atomic data.

module Atomic_Data_Type
  !% Defines a data type for storing atomic data.
  private
  public :: atomicData
  
  type atomicData
     !% Data type for storing atomic data.
     integer                                     :: atomicNumber
     double precision                            :: atomicMass
     double precision, allocatable, dimension(:) :: abundanceByMass
     character(len=3)                            :: shortLabel
     character(len=20)                           :: name
  end type atomicData

end module Atomic_Data_Type
