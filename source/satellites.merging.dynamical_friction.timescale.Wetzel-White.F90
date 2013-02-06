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

!+    Contributions to this file made by:  Martin White.

!% Contains a module which implements calculations of satellite merging times using the \cite{wetzel_what_2010} method.

module Dynamical_Friction_Wetzel_White
  !% Implements calculations of satellite merging times using the \cite{wetzel_what_2010} method.
  implicit none
  private
  public :: Satellite_Time_Until_Merging_Wetzel_White_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Wetzel_White_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Wetzel_White_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: satelliteMergingMethod
    procedure(double precision), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Wetzel-White2010') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Wetzel_White
    return
  end subroutine Satellite_Time_Until_Merging_Wetzel_White_Initialize

  double precision function Satellite_Time_Until_Merging_Wetzel_White(thisNode)
    !% Return the timescale for merging satellites using the \cite{wetzel_what_2010} method.
    use Galacticus_Nodes
    use Dynamical_Friction_Timescale_Utilities
    use Cosmology_Functions
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    type (treeNode          ), pointer                :: hostNode
    class(nodeComponentBasic), pointer                :: thisBasicComponent,hostBasicComponent
    double precision,          parameter              :: timeScaleNormalization=0.2d0 ! C_dyn from Wetzel & White (2010).
    double precision                                  :: massRatio

    ! Find the host node.
    hostNode => thisNode%parent
    ! Compute mass ratio.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    massRatio=hostBasicComponent%mass()/thisBasicComponent%mass()
    ! Compute dynamical friction timescale using eqn. (2) from Wetzel & White (2010).
    Satellite_Time_Until_Merging_Wetzel_White= Dynamical_Friction_Timescale_Multiplier()                   &
         &                                    *timeScaleNormalization                                      &
         &                                    /Expansion_Rate(Expansion_Factor(thisBasicComponent%time())) &
         &                                    *           massRatio                                        &
         &                                    /dlog(1.0d0+massRatio)
    return
  end function Satellite_Time_Until_Merging_Wetzel_White

end module Dynamical_Friction_Wetzel_White
