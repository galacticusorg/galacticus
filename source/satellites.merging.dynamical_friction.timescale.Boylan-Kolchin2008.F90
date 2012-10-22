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

!% Contains a module which implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.

module Dynamical_Friction_BoylanKolchin2008
  !% Implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.
  implicit none
  private
  public :: Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(in)    :: satelliteMergingMethod
    procedure(double precision), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'BoylanKolchin2008') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_BoylanKolchin2008
    return
  end subroutine Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize

  double precision function Satellite_Time_Until_Merging_BoylanKolchin2008(thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{boylan-kolchin_dynamical_2008} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits
    use Satellite_Orbits
    implicit none
    type (treeNode   ),        pointer, intent(inout) :: thisNode
    type (keplerOrbit),                 intent(inout) :: thisOrbit
    type (treeNode   ),        pointer                :: hostNode
    class(nodeComponentBasic), pointer                :: thisBasicComponent,hostBasicComponent
    logical,                   parameter              :: acceptUnboundOrbits=.false.
    double precision,          parameter              :: A=0.216d0, b=1.3d0, c=1.9d0, d=1.0d0 ! Fitting parameters from eqn. (6) of Boylan-Kolchin et al.
    double precision                                  :: equivalentCircularOrbitRadius,orbitalCircularity,velocityScale,radialScale&
         &,massRatio

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit)
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute orbital circularity.
    orbitalCircularity=thisOrbit%angularMomentum()/equivalentCircularOrbitRadius/Dark_Matter_Profile_Circular_Velocity(hostNode&
         &,equivalentCircularOrbitRadius)
    ! Compute mass ratio (mass in host [not including satellite] divided by mass in satellite).
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    massRatio=hostBasicComponent%mass()/thisBasicComponent%mass()-1.0d0
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       Satellite_Time_Until_Merging_BoylanKolchin2008=0.0d0
    else
       ! Compute dynamical friction timescale.
       Satellite_Time_Until_Merging_BoylanKolchin2008=Dynamical_Friction_Timescale_Multiplier() &
            &*Dark_Matter_Halo_Dynamical_Timescale(hostNode)*A*((massRatio**b)/dlog(1.0d0 +massRatio))*dexp(c*orbitalCircularity)&
            & *((equivalentCircularOrbitRadius/radialScale)**d)
    end if
    return
  end function Satellite_Time_Until_Merging_BoylanKolchin2008

end module Dynamical_Friction_BoylanKolchin2008
