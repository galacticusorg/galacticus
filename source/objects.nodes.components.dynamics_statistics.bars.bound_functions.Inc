!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains custom functions for the bars dynamics statistics component.
!!}

subroutine Node_Component_Dynamics_Statistics_Bars_Record(self,time,barInstabilityTimescale,adiabaticRatio)
  !!{
  Record the dynamical state.
  !!}
  use :: Memory_Management, only : allocateArray, deallocateArray
  implicit none
  class           (nodeComponentDynamicsStatisticsBars), intent(inout)               :: self
  double precision                                     , intent(in   )               :: time          , barInstabilityTimescale, &
       &                                                                                adiabaticRatio
  double precision                                     , allocatable  , dimension(:) :: arrayTemporary
  integer                                                                            :: entryIndex

  ! Expand the arrays as necessary.
  if (allocated(self%barInstabilityTimescaleData)) then
     call Move_Alloc (self%barInstabilityTimescaleData,        arrayTemporary  )
     call allocateArray(self%barInstabilityTimescaleData,[1+size(arrayTemporary)])
     self%barInstabilityTimescaleData(1:size(arrayTemporary))=arrayTemporary
     call deallocateArray(arrayTemporary)
  else
     call allocateArray(self%barInstabilityTimescaleData,[1])
  end if
  if (allocated(self%adiabaticRatioData         )) then
     call Move_Alloc (self%adiabaticRatioData         ,        arrayTemporary  )
     call allocateArray(self%adiabaticRatioData         ,[1+size(arrayTemporary)])
     self%adiabaticRatioData         (1:size(arrayTemporary))=arrayTemporary
     call deallocateArray(arrayTemporary)
  else
     call allocateArray(self%barInstabilityTimescaleData,[1])
  end if
  if (allocated(self%timeData                   )) then
     call Move_Alloc (self%timeData                   ,        arrayTemporary  )
     call allocateArray(self%timeData                   ,[1+size(arrayTemporary)])
     self%timeData                   (1:size(arrayTemporary))=arrayTemporary
     call deallocateArray(arrayTemporary)
  else
     call allocateArray(self%timeData                   ,[1])
  end if
  ! Find the current extent of the arrays.
  entryIndex=size(self%timeData)
  ! Record the current time.
  self%timeData                   (entryIndex)=time
  ! Record the bar instability timescale
  self%barInstabilityTimescaleData(entryIndex)=barInstabilityTimescale
  ! Record the adiabatic ratio
  self%adiabaticRatioData         (entryIndex)=adiabaticRatio
  return
end subroutine Node_Component_Dynamics_Statistics_Bars_Record
