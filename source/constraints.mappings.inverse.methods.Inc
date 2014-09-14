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

!% Defines methods for a inverse mapping.

function mappingInverseConstructor()
  !% Constructor for inverse mapping class.
  implicit none
  type(mappingInverse) :: mappingInverseConstructor

  return
end function mappingInverseConstructor

subroutine mappingInverseMapPrior(self,thisPrior)
  !% Map a prior through a inverse mapping.
  use Galacticus_Error
  use Statistics_Distributions
  implicit none
  class(mappingInverse), intent(in   ) :: self
  class(prior         ), intent(inout) :: thisPrior

  select type (thisDistribution => thisPrior%priorDistribution)
  type is (distributionUniform)
     ! Simply map the upper and lower limits of this uniform distribution.
     if     (                                                                                                      &
          &   thisDistribution%limitLower <= 0.0d0                                                                 &
          &  .and.                                                                                                 &
          &   thisDistribution%limitUpper >= 0.0d0                                                                 &
          & ) call Galacticus_Error_Report(                                                                        &
          &                                'mappingInverseMapPrior'                                              , &
          &                                'inverse mapping cannot be applied to uniform distribution spans zero'  &
          &                               )
     thisDistribution%limitLower=1.0d0/thisDistribution%limitLower
     thisDistribution%limitUpper=1.0d0/thisDistribution%limitUpper
  class default
     ! Report an error.
     call Galacticus_Error_Report('mappingInverseMapPrior','mapping can not handle this prior distribution')
  end select
  return
end subroutine mappingInverseMapPrior

double precision function mappingInverseMap(self,x)
  !% Map a value via a inverse mapping.
  implicit none
  class           (mappingInverse), intent(in   ) :: self
  double precision                , intent(in   ) :: x
  
  mappingInverseMap=1.0d0/x
  return
end function mappingInverseMap

double precision function mappingInverseUnmap(self,x)
  !% Unmap a value via a inverse mapping.
  implicit none
  class           (mappingInverse), intent(in   ) :: self
  double precision                , intent(in   ) :: x
  
  mappingInverseUnmap=1.0d0/x
  return
end function mappingInverseUnmap

