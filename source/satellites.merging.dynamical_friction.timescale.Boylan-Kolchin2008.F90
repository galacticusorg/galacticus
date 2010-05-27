!% Contains a module which implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.

module Dynamical_Friction_BoylanKolchin2008
  !% Implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.
  private
  public :: Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in)    :: satelliteMergingMethod
    procedure(), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'BoylanKolchin2008') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_BoylanKolchin2008
    return
  end subroutine Satellite_Time_Until_Merging_BoylanKolchin2008_Initialize

  double precision function Satellite_Time_Until_Merging_BoylanKolchin2008(thisNode)
    !% Return the timescale for merging satellites using the \cite{boylan-kolchin_dynamical_2008} method.
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
    double precision, parameter              :: A=0.216d0, b=1.3d0, c=1.9d0, d=1.0d0 ! Fitting parameters from eqn. (6) of Boylan-Kolchin et al.
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
    Satellite_Time_Until_Merging_BoylanKolchin2008=Dark_Matter_Halo_Dynamical_Timescale(hostNode)*A*((massRatio**b)/dlog(1.0d0&
         &+massRatio))*dexp(c*orbitalCircularity)*((equivalentCircularOrbitRadius/radialScale)**d)
    return
  end function Satellite_Time_Until_Merging_BoylanKolchin2008

end module Dynamical_Friction_BoylanKolchin2008
