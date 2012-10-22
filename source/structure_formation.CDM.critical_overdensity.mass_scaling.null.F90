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

!% Contains a module which implements a null scaling of critical overdensities for collapse.

module Critical_Overdensity_Mass_Scalings_Null
  !% Implements a null scaling of critical overdensities for collapse.
  use Galacticus_Nodes
  implicit none
  private
  public :: Critical_Overdensity_Mass_Scaling_Null_Initialize
  
contains

  !# <criticalOverdensityMassScalingMethod>
  !#  <unitName>Critical_Overdensity_Mass_Scaling_Null_Initialize</unitName>
  !# </criticalOverdensityMassScalingMethod>
  subroutine Critical_Overdensity_Mass_Scaling_Null_Initialize(criticalOverdensityMassScalingMethod&
       &,Critical_Overdensity_Mass_Scaling_Get,Critical_Overdensity_Mass_Scaling_Gradient_Get)
    !% Initializes the ``null'' critical overdensity mass scaling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: criticalOverdensityMassScalingMethod
    procedure(double precision), pointer, intent(inout) :: Critical_Overdensity_Mass_Scaling_Get,Critical_Overdensity_Mass_Scaling_Gradient_Get
    
    if (criticalOverdensityMassScalingMethod == 'null') then
       Critical_Overdensity_Mass_Scaling_Get          => Critical_Overdensity_Mass_Scaling_Null
       Critical_Overdensity_Mass_Scaling_Gradient_Get => Critical_Overdensity_Mass_Scaling_Gradient_Null
    end if
    return
  end subroutine Critical_Overdensity_Mass_Scaling_Null_Initialize

  double precision function Critical_Overdensity_Mass_Scaling_Null(mass)
    !% Returns a mass scaling for critical overdensities that is always unity.
    implicit none
    double precision, intent(in) :: mass

    Critical_Overdensity_Mass_Scaling_Null=1.0d0
    return
  end function Critical_Overdensity_Mass_Scaling_Null
  
  double precision function Critical_Overdensity_Mass_Scaling_Gradient_Null(mass)
    !% Returns the gradient of a mass scaling for critical overdensities that is always unity.
    implicit none
    double precision, intent(in) :: mass

    Critical_Overdensity_Mass_Scaling_Gradient_Null=0.0d0
    return
  end function Critical_Overdensity_Mass_Scaling_Gradient_Null
  
end module Critical_Overdensity_Mass_Scalings_Null
