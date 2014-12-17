!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the bars dynamics statistics component.

subroutine Node_Component_Dynamics_Statistics_Bars_Record(self,time,barInstabilityTimescale,adiabaticRatio)
  !% Record the dynamical state.
  use Memory_Management
  implicit none
  class           (nodeComponentDynamicsStatisticsBars), intent(inout)               :: self
  double precision                                     , intent(in   )               :: time          , barInstabilityTimescale, &
       &                                                                                adiabaticRatio
  double precision                                     , allocatable  , dimension(:) :: arrayTemporary
  integer                                                                            :: entryIndex

  ! Expand the arrays as necessary.
  if (allocated(self%barInstabilityTimescaleData%value)) then
     call Move_Alloc (self%barInstabilityTimescaleData%value,        arrayTemporary  )
     call Alloc_Array(self%barInstabilityTimescaleData%value,[1+size(arrayTemporary)])
     self%barInstabilityTimescaleData%value(1:size(arrayTemporary))=arrayTemporary
     call Dealloc_Array(arrayTemporary)
  else
     call Alloc_Array(self%barInstabilityTimescaleData%value,[1])
  end if
  if (allocated(self%adiabaticRatioData         %value)) then
     call Move_Alloc (self%adiabaticRatioData         %value,        arrayTemporary  )
     call Alloc_Array(self%adiabaticRatioData         %value,[1+size(arrayTemporary)])
     self%adiabaticRatioData         %value(1:size(arrayTemporary))=arrayTemporary
     call Dealloc_Array(arrayTemporary)
  else
     call Alloc_Array(self%barInstabilityTimescaleData%value,[1])
  end if
  if (allocated(self%timeData                   %value)) then
     call Move_Alloc (self%timeData                   %value,        arrayTemporary  )
     call Alloc_Array(self%timeData                   %value,[1+size(arrayTemporary)])
     self%timeData                   %value(1:size(arrayTemporary))=arrayTemporary
     call Dealloc_Array(arrayTemporary)
  else
     call Alloc_Array(self%timeData                   %value,[1])
  end if
  ! Find the current extent of the arrays.
  entryIndex=size(self%timeData%value)
  ! Record the current time.
  self%timeData                   %value(entryIndex)=time
  ! Record the bar instability timescale
  self%barInstabilityTimescaleData%value(entryIndex)=barInstabilityTimescale
  ! Record the adiabatic ratio
  self%adiabaticRatioData         %value(entryIndex)=adiabaticRatio
  return
end subroutine Node_Component_Dynamics_Statistics_Bars_Record
