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

!% Contains custom functions for the exponential disk component.

subroutine Node_Component_Disk_Exponential_Attach_Pipes(self)
  !% Attach cooling pipes to the exponential disk component.
  use Galacticus_Error
  implicit none
  class(nodeComponentDiskExponential), intent(in) :: self
  type (nodeComponentHotHalo        )             :: hotHalo
  
  if (hotHalo%coolingMassRateIsAttached           ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo mass cooling pipe'            )
  if (hotHalo%coolingAngularMomentumRateIsAttached()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo angular momentum cooling pipe')
  if (hotHalo%coolingAbundancesRateIsAttached     ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo abundances cooling pipe'      )
  call hotHalo%coolingMassRateFunction           (DiskExponentialMassGasRateGeneric        )
  call hotHalo%coolingAngularMomentumRateFunction(DiskExponentialAngularMomentumRateGeneric)
  call hotHalo%coolingAbundancesRateFunction     (DiskExponentialAbundancesGasRateGeneric  )
  return
end subroutine Node_Component_Disk_Exponential_Attach_Pipes

double precision function Node_Component_Disk_Exponential_Half_Mass_Radius(self)
  !% Return the half-mass radius of the exponential disk.
  implicit none
  class(nodeComponentDiskExponential), intent(inout) :: self
  double precision                   , parameter     :: halfMassRadiusToScaleRadius=1.678346990d0

  Node_Component_Disk_Exponential_Half_Mass_Radius=self%radius()*halfMassRadiusToScaleRadius
  return
end function Node_Component_Disk_Exponential_Half_Mass_Radius
