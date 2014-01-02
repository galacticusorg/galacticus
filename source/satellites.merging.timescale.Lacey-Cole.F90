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

  !% Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesLaceyCole1993" />

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesLaceyCole1993
     !% A class implementing the \cite{lacey_merger_1993} method for satellite merging timescales.
     private
   contains
     !@ <objectMethods>
     !@   <object>satelliteMergingTimescalesLaceyCole1993</object>
     !@   <objectMethod>
     !@     <method>timeUntilMergingMassDependence</method>
     !@     <type>double precision</type>
     !@     <arguments></arguments>
     !@     <description>Return the mass-dependent part of the time (in Gyr) until the satellite will merge with its host.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                   laceyCole1993Destructor
     procedure :: timeUntilMerging               => laceyCole1993TimeUntilMerging
     procedure :: timeUntilMergingMassDependence => laceyCole1993TimeUntilMergingMassDependence
  end type satelliteMergingTimescalesLaceyCole1993

  interface satelliteMergingTimescalesLaceyCole1993
     !% Constructors for the \cite{lacey_merger_1993} merging timescale class.
     module procedure laceyCole1993DefaultConstructor
  end interface satelliteMergingTimescalesLaceyCole1993

contains

  function laceyCole1993DefaultConstructor()
    !% Default constructor for the \cite{lacey_merger_1993} merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesLaceyCole1993) :: laceyCole1993DefaultConstructor

    return
  end function laceyCole1993DefaultConstructor

  elemental subroutine laceyCole1993Destructor(self)
    !% Default constructor for the \cite{lacey_merger_1993} merging timescale class.
    implicit none
    type(satelliteMergingTimescalesLaceyCole1993), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine laceyCole1993Destructor

  double precision function laceyCole1993TimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    use Galacticus_Nodes
    use Kepler_Orbits
    use Dark_Matter_Halo_Scales
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993), intent(inout)          :: self
    type            (treeNode                               ), intent(inout), pointer :: thisNode
    type            (keplerOrbit                            ), intent(inout)          :: thisOrbit
    type            (treeNode                               )               , pointer :: hostNode
    double precision                                                                  :: equivalentCircularOrbitRadius, orbitalCircularity, &
         &                                                                               radialScale                  , velocityScale

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute radius of orbit with same energy.
    equivalentCircularOrbitRadius=exp(thisOrbit%energy()/velocityScale**2+0.5d0)
    ! Compute orbital circularity.
    orbitalCircularity=thisOrbit%angularMomentum()/velocityScale/radialScale/equivalentCircularOrbitRadius
    ! Compute the merging timescale.
    laceyCole1993TimeUntilMerging                         &
         & =orbitalCircularity           **0.78d0         &
         & *equivalentCircularOrbitRadius**2              &
         & *self%timeUntilMergingMassDependence(thisNode)
    return
  end function laceyCole1993TimeUntilMerging

  double precision function laceyCole1993TimeUntilMergingMassDependence(self,thisNode)
    !% Return the mass-dependent part of the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Dynamical_Friction_Timescale_Utilities
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993), intent(inout)          :: self
    type            (treeNode                               ), intent(inout), pointer :: thisNode
    type            (treeNode                               )               , pointer :: hostNode
    class           (nodeComponentBasic                     )               , pointer :: hostBasic                 , thisBasic
    double precision                                         , parameter              :: inverseTwoB1=1.169335453d0            !  1/2/B(1).
    double precision                                                                  :: massRatio

    ! Find the host node.
    hostNode => thisNode%parent
    ! Compute mass ratio.
    thisBasic => thisNode%basic()
    hostBasic => hostNode%basic()
    massRatio=hostBasic%mass()/thisBasic%mass()
    ! Check for a greater than unity mass ratio.
    if (massRatio <= 1.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       laceyCole1993TimeUntilMergingMassDependence=0.0d0
    else
       ! Compute dynamical friction timescale.
       laceyCole1993TimeUntilMergingMassDependence            &
            & =Dynamical_Friction_Timescale_Multiplier()      &
            & *Dark_Matter_Halo_Dynamical_Timescale(hostNode) &
            & *inverseTwoB1                                   &
            & *    massRatio                                  &
            & /log(massRatio)
    end if
    return
  end function laceyCole1993TimeUntilMergingMassDependence
