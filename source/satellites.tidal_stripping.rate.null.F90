!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module with a null implementation of calculations of satellite mass loss due to tidal stripping.

module Tidal_Stripping_Rate_Null
  !% Implements null value of calculations of satellite mass loss due to tidal stripping.
  implicit none
  private
  public :: Satellite_Tidal_Stripping_Rate_Null_Initialize

contains

  !# <satelliteTidalStrippingMethod>
  !#  <unitName>Satellite_Tidal_Stripping_Rate_Null_Initialize</unitName>
  !# </satelliteTidalStrippingMethod>
  subroutine Satellite_Tidal_Stripping_Rate_Null_Initialize(satelliteTidalStrippingMethod,Satellite_Tidal_Stripping_Rate)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string                     ),          intent(in   ) :: satelliteTidalStrippingMethod
    procedure(Satellite_Tidal_Stripping_Rate_Null), pointer, intent(inout) :: Satellite_Tidal_Stripping_Rate

    if (satelliteTidalStrippingMethod == 'null') Satellite_Tidal_Stripping_Rate => Satellite_Tidal_Stripping_Rate_Null
    return
  end subroutine Satellite_Tidal_Stripping_Rate_Null_Initialize

  double precision function Satellite_Tidal_Stripping_Rate_Null(thisNode)
    !% Return the null mass loss for satellites due to tidal stripping.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), target :: thisNode
    !GCC$ attributes unused :: thisNode
    
    Satellite_Tidal_Stripping_Rate_Null=0.0d0
    return
  end function Satellite_Tidal_Stripping_Rate_Null

end module Tidal_Stripping_Rate_Null
