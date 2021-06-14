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

!% Contains custom functions for the standard hot halo component.

subroutine Node_Component_Hot_Halo_Standard_Mass_Removal_Rate(self,massRemovalRate)
  !% Handle instances of mass removal from the standard hot halo component
  !% class. For the standard hot halo component type, we do nothing.
  implicit none
  class          (nodeComponentHotHaloStandard), intent(inout) :: self
  double precision                             , intent(in   ) :: massRemovalRate
  !$GLC attributes unused :: self, massRemovalRate

  return
end subroutine Node_Component_Hot_Halo_Standard_Mass_Removal_Rate

double precision function Node_Component_Hot_Halo_Standard_Mass_Total(self)
  !% Return the total active mass in the hot halo.
  implicit none
  class(nodeComponentHotHaloStandard), intent(inout) :: self

  Node_Component_Hot_Halo_Standard_Mass_Total=self%mass()
  return
end function Node_Component_Hot_Halo_Standard_Mass_Total

