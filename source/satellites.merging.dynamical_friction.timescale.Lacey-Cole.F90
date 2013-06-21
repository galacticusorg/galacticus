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

!% Contains a module which implements calculations of satellite merging times using the \cite{lacey_merger_1993} method.

module Dynamical_Friction_Lacey_Cole
  !% Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method.
  implicit none
  private
  public :: Satellite_Time_Until_Merging_Lacey_Cole_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string                         ), intent(in   )          :: satelliteMergingMethod
    procedure(Satellite_Time_Until_Merging_Lacey_Cole), intent(inout), pointer :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Lacey-Cole') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Lacey_Cole
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Initialize

  double precision function Satellite_Time_Until_Merging_Lacey_Cole(thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (keplerOrbit       ), intent(inout)          :: thisOrbit
    type            (treeNode          )               , pointer :: hostNode
    class           (nodeComponentBasic)               , pointer :: hostBasicComponent                         , thisBasicComponent
    double precision                    , parameter              :: inverseTwoB1                 =1.169335453d0                     !  1/2/B(1).
    double precision                                             :: equivalentCircularOrbitRadius              , massRatio                       , &
         &                                                          orbitalCircularity                         , radialScale                     , &
         &                                                          velocityScale

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute radius of orbit with same energy.
    equivalentCircularOrbitRadius=exp(thisOrbit%energy()/velocityScale**2+0.5d0)
    ! Compute orbital circularity.
    orbitalCircularity=thisOrbit%angularMomentum()/velocityScale/radialScale/equivalentCircularOrbitRadius
    ! Compute mass ratio.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    massRatio=hostBasicComponent%mass()/thisBasicComponent%mass()
    ! Check for a greater than unity mass ratio.
    if (massRatio <= 1.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       Satellite_Time_Until_Merging_Lacey_Cole=0.0d0
    else
       ! Compute dynamical friction timescale.
       Satellite_Time_Until_Merging_Lacey_Cole=Dynamical_Friction_Timescale_Multiplier()*(orbitalCircularity**0.78d0)&
            &*(equivalentCircularOrbitRadius**2) *Dark_Matter_Halo_Dynamical_Timescale(hostNode)*inverseTwoB1*massRatio&
            &/log(massRatio)
    end if
    return
  end function Satellite_Time_Until_Merging_Lacey_Cole

end module Dynamical_Friction_Lacey_Cole
