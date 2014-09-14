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

!% Defines methods for a linear mapping.

function mappingLinearConstructor()
  !% Constructor for linear mapping class.
  implicit none
  type(mappingLinear) :: mappingLinearConstructor

  return
end function mappingLinearConstructor

subroutine mappingLinearMapPrior(self,thisPrior)
  !% Map a prior through a linear mapping (i.e. do nothing).
  implicit none
  class(mappingLinear), intent(in   ) :: self
  class(prior        ), intent(inout) :: thisPrior

  return
end subroutine mappingLinearMapPrior

double precision function mappingLinearMap(self,x)
  !% Map a value via a linear mapping (i.e. do nothing).
  implicit none
  class           (mappingLinear), intent(in   ) :: self
  double precision               , intent(in   ) :: x
  
  mappingLinearMap=x
  return
end function mappingLinearMap

double precision function mappingLinearUnmap(self,x)
  !% Unmap a value via a linear mapping (i.e. do nothing).
  implicit none
  class           (mappingLinear), intent(in   ) :: self
  double precision               , intent(in   ) :: x
  
  mappingLinearUnmap=x
  return
end function mappingLinearUnmap

