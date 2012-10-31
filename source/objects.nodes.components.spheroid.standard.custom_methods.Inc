!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the standard spheroid component.

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self)
  !% Return the half-mass radius of the standard spheroid.
  use Node_Component_Spheroid_Standard_Data
  use Mass_Distributions
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self
  
  select type (spheroidMassDistribution)
  class is (massDistributionSpherical)
     Node_Component_Spheroid_Standard_Half_Mass_Radius=self%radius()*spheroidMassDistribution%halfMassRadius()
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius
