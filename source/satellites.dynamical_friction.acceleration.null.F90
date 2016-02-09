!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module with a null implementation of calculations of satellite acceleration due to dynamical friction.

module Dynamical_Friction_Acceleration_Null
  !% Implements null value of calculations of satellite acceleration due to dynamical friction.
  implicit none
  private
  public :: Satellite_Dynamical_Friction_Acceleration_Null_Initialize

contains

  !# <satelliteDynamicalFrictionMethod>
  !#  <unitName>Satellite_Dynamical_Friction_Acceleration_Null_Initialize</unitName>
  !# </satelliteDynamicalFrictionMethod>
  subroutine Satellite_Dynamical_Friction_Acceleration_Null_Initialize(satelliteDynamicalFrictionMethod,Satellite_Dynamical_Friction_Acceleration)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                ),          intent(in   ) :: satelliteDynamicalFrictionMethod
    procedure(Satellite_Dynamical_Friction_Acceleration_Null), pointer, intent(inout) :: Satellite_Dynamical_Friction_Acceleration

    if (satelliteDynamicalFrictionMethod == 'null') Satellite_Dynamical_Friction_Acceleration => Satellite_Dynamical_Friction_Acceleration_Null
    return
  end subroutine Satellite_Dynamical_Friction_Acceleration_Null_Initialize

  function Satellite_Dynamical_Friction_Acceleration_Null(thisNode)
    !% Return a null acceleration for satellites due to dynamical friction.
    use Galacticus_Nodes
    implicit none
    double precision          , dimension(3)                :: Satellite_Dynamical_Friction_Acceleration_Null
    type            (treeNode), pointer     , intent(inout) :: thisNode

    Satellite_Dynamical_Friction_Acceleration_Null=[0.0d0,0.0d0,0.0d0]
    return
  end function Satellite_Dynamical_Friction_Acceleration_Null

end module Dynamical_Friction_Acceleration_Null
