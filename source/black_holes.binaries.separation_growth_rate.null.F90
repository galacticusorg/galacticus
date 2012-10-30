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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements a black hole binary separation growth rate which is always zero.

module Black_Hole_Binary_Separations_Null
  !% Implements a black hole binary initial separation growth rate which is always zero.
  implicit none
  private
  public :: Black_Hole_Binary_Separation_Growth_Rate_Null_Initialize

contains

  !# <blackHoleBinarySeparationGrowthRateMethod>
  !#  <unitName>Black_Hole_Binary_Separation_Growth_Rate_Null_Initialize</unitName>
  !# </blackHoleBinarySeparationGrowthRateMethod>
  subroutine Black_Hole_Binary_Separation_Growth_Rate_Null_Initialize(blackHoleBinarySeparationGrowthRateMethod,Black_Hole_Binary_Separation_Growth_Rate_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinarySeparationGrowthRateMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Separation_Growth_Rate_Get
    
    if (blackHoleBinarySeparationGrowthRateMethod == 'null') Black_Hole_Binary_Separation_Growth_Rate_Get => Black_Hole_Binary_Separation_Growth_Rate_Null
    return
  end subroutine Black_Hole_Binary_Separation_Growth_Rate_Null_Initialize

  double precision function Black_Hole_Binary_Separation_Growth_Rate_Null(thisBlackHoleComponent)
    !% Returns a separation growth rate for a binary black hole that is always zero.
    use Galacticus_Nodes
    implicit none
    class(nodeComponentBlackHole), intent(inout), pointer :: thisBlackHoleComponent

    Black_Hole_Binary_Separation_Growth_Rate_Null=0.0d0 
    return
  end function Black_Hole_Binary_Separation_Growth_Rate_Null

end module Black_Hole_Binary_Separations_Null
