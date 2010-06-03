!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.

module Dynamical_Friction_Jiang2008
  !% Implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.
  private
  public :: Satellite_Time_Until_Merging_Jiang2008_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Jiang2008_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Jiang2008_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in)    :: satelliteMergingMethod
    procedure(), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Jiang2008') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Jiang2008
    return
  end subroutine Satellite_Time_Until_Merging_Jiang2008_Initialize

  double precision function Satellite_Time_Until_Merging_Jiang2008(thisNode)
    !% Return the timescale for merging satellites using the \cite{jiang_fitting_2008} method.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Virial_Orbits
    use Numerical_Constants_Math
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    type(treeNode),   pointer                :: hostNode
    logical,          parameter              :: acceptUnboundOrbits=.false.
    double precision, parameter              :: a=0.94d0, b=0.60d0, C=0.43d0, d=0.60d0 ! Fitting parameters from Jiang's paper.
    double precision                         :: angularMomentum,orbitalEnergy,equivalentCircularOrbitRadius,orbitalCircularity &
         &,velocityScale,radialScale,massRatio

    ! Find the host node.
    hostNode => thisNode%parentNode
    ! Get orbital parameters for this satellite.
    call Virial_Orbital_Parameters(thisNode,acceptUnboundOrbits,angularMomentum=angularMomentum,orbitalEnergy=orbitalEnergy&
         &,equivalentCircularOrbitRadius=equivalentCircularOrbitRadius)
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute orbital circularity.
    orbitalCircularity=angularMomentum/equivalentCircularOrbitRadius/Dark_Matter_Profile_Circular_Velocity(hostNode&
         &,equivalentCircularOrbitRadius)
    ! Compute mass ratio (mass in host [not including satellite] divided by mass in satellite).
    massRatio=Tree_Node_Mass(hostNode)/Tree_Node_Mass(thisNode)-1.0d0
    ! Compute dynamical friction timescale.
    Satellite_Time_Until_Merging_Jiang2008=Dark_Matter_Halo_Dynamical_Timescale(hostNode)*dsqrt(equivalentCircularOrbitRadius&
         &/radialScale)*((a*(orbitalCircularity**b)+d)/2.0/C)*massRatio/dlog(1.0d0+massRatio)
    return
  end function Satellite_Time_Until_Merging_Jiang2008

end module Dynamical_Friction_Jiang2008
