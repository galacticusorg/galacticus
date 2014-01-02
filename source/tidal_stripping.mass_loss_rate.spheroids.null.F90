!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null mass loss rates from spheroids due to tidal stripping.

module Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null
  !% Implements a null mass loss rates from spheroids due to tidal stripping.
  use Galacticus_Nodes
  implicit none
  private
  public :: Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null_Init

contains

  !# <tidalStrippingMassLossRateSpheroidsMethod>
  !#  <unitName>Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null_Init</unitName>
  !# </tidalStrippingMassLossRateSpheroidsMethod>
  subroutine Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null_Init(tidalStrippingMassLossRateSpheroidsMethod,Tidal_Stripping_Mass_Loss_Rate_Spheroid_Get)
    !% Initializes the ``null'' tidal stripping mass loss rate from spheroids module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                              ), intent(in   )          :: tidalStrippingMassLossRateSpheroidsMethod
    procedure(Tidal_Stripping_Mass_Loss_Rate_Spheroid_Null), intent(inout), pointer :: Tidal_Stripping_Mass_Loss_Rate_Spheroid_Get

    if (tidalStrippingMassLossRateSpheroidsMethod == 'null') Tidal_Stripping_Mass_Loss_Rate_Spheroid_Get => Tidal_Stripping_Mass_Loss_Rate_Spheroid_Null
    return
  end subroutine Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null_Init

  double precision function Tidal_Stripping_Mass_Loss_Rate_Spheroid_Null(thisNode)
    !% Computes the mass loss rate from spheroids due to tidal stripping. Always returns zero.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Tidal_Stripping_Mass_Loss_Rate_Spheroid_Null=0.0d0
    return
  end function Tidal_Stripping_Mass_Loss_Rate_Spheroid_Null

end module Tidal_Stripping_Mass_Loss_Rate_Spheroids_Null
