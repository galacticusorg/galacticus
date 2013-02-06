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

!% Contains a module which implements a zero cooling rate calculation.

module Cooling_Rates_Zero
  !% Implements zero cooling rate calculation.
  implicit none
  private
  public :: Cooling_Rate_Zero_Initialize

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_Zero_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_Zero_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``zero'' cooling rate module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in   ) :: coolingRateMethod
    procedure(double precision), pointer, intent(inout) :: Cooling_Rate_Get
    
    if (coolingRateMethod == 'zero') Cooling_Rate_Get => Cooling_Rate_Zero
    return
  end subroutine Cooling_Rate_Zero_Initialize

  double precision function Cooling_Rate_Zero(thisNode)
    !% Returns a zero mass cooling rate in a hot gas halos.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Cooling_Rate_Zero=0.0d0
    return
  end function Cooling_Rate_Zero

end module Cooling_Rates_Zero
