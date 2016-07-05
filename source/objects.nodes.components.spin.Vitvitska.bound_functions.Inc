!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+ Contributions to this file made by: Chrsitoph Behrens.

!% Contains custom functions for the Vitvitska spin component.

double precision function spinVitvitskaSpinGrowthRate(self)
  !% Return the {\normalfont \ttfamily spinGrowthRate} property of the {\normalfont \ttfamily
  !% SpinVitvitska} component. Spin does not undergo differential evolution along branches in
  !% this model.
  implicit none
  class(nodeComponentSpinVitvitska), intent(inout) :: self

  spinVitvitskaSpinGrowthRate=0.0d0
  return
end function spinVitvitskaSpinGrowthRate
