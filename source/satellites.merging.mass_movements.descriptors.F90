!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines descriptors for satellite merger mass movements.

module Satellite_Merging_Mass_Movements_Descriptors
  !% Defines descriptors for satellite merger mass movements.
  implicit none
  public

  integer, parameter :: doesNotMove    =0
  integer, parameter :: movesToDisk    =1
  integer, parameter :: movesToSpheroid=2
  
  ! Stored mass movement descriptors for the current merging event.
  integer            :: thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo
  !$omp threadprivate(thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo)
  logical            :: thisMergerIsMajor
  !$omp threadprivate(thisMergerIsMajor)

end module Satellite_Merging_Mass_Movements_Descriptors
