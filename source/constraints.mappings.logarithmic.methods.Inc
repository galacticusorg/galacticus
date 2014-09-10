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

!% Defines methods for a logarithmic mapping.

function mappingLogarithmicConstructor()
  !% Constructor for logarithmic mapping class.
  implicit none
  type(mappingLogarithmic) :: mappingLogarithmicConstructor

  return
end function mappingLogarithmicConstructor

subroutine mappingLogarithmicMapPrior(self,thisPrior)
  !% Map a prior through a logarithmic mapping (i.e. do nothing).
  use Galacticus_Error
  use Statistics_Distributions
  implicit none
  class(mappingLogarithmic), intent(in   ) :: self
  class(prior             ), intent(inout) :: thisPrior

  select type (thisDistribution => thisPrior%priorDistribution)
  type is (distributionUniform)
     ! Simply map the upper and lower limits of this uniform distribution.
     if     (                                                                                                                   &
          &   thisDistribution%limitLower <= 0.0d0                                                                              &
          &  .or.                                                                                                               &
          &   thisDistribution%limitUpper <= 0.0d0                                                                              &
          & ) call Galacticus_Error_Report(                                                                                     &
          &                                'mappingLogarithmicMapPrior'                                                       , &
          &                                'logarithmic mapping cannot be applied to uniform distribution with negative range'  &
          &                               )
     thisDistribution%limitLower=log(thisDistribution%limitLower)
     thisDistribution%limitUpper=log(thisDistribution%limitUpper)
  class default
     ! Report an error.
     call Galacticus_Error_Report('mappingLogarithmicMapPrior','mapping can not handle this prior distribution')
  end select
  return
end subroutine mappingLogarithmicMapPrior

double precision function mappingLogarithmicMap(self,x)
  !% Map a value via a logarithmic mapping.
  implicit none
  class           (mappingLogarithmic), intent(in   ) :: self
  double precision                    , intent(in   ) :: x
  
  mappingLogarithmicMap=log(x)
  return
end function mappingLogarithmicMap

double precision function mappingLogarithmicUnmap(self,x)
  !% Unmap a value via a logarithmic mapping.
  implicit none
  class           (mappingLogarithmic), intent(in   ) :: self
  double precision                    , intent(in   ) :: x
  
  mappingLogarithmicUnmap=exp(x)
  return
end function mappingLogarithmicUnmap

