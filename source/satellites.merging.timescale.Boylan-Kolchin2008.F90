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

  !% Implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.
 
  !# <satelliteMergingTimescales name="satelliteMergingTimescalesBoylanKolchin2008" />

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesBoylanKolchin2008
     !% A class implementing the \cite{boylan-kolchin_dynamical_2008} method for satellite merging timescales.
     private
   contains
     final     ::                     boylanKolchin2008Destructor
     procedure :: timeUntilMerging => boylanKolchin2008TimeUntilMerging
  end type satelliteMergingTimescalesBoylanKolchin2008

  interface satelliteMergingTimescalesBoylanKolchin2008
     !% Constructors for the \cite{boylan-kolchin_dynamical_2008} merging timescale class.
     module procedure boylanKolchin2008DefaultConstructor
  end interface satelliteMergingTimescalesBoylanKolchin2008

contains

  function boylanKolchin2008DefaultConstructor()
    !% Default constructor for the \cite{boylan-kolchin_dynamical_2008} merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesBoylanKolchin2008) :: boylanKolchin2008DefaultConstructor

    return
  end function boylanKolchin2008DefaultConstructor

  elemental subroutine boylanKolchin2008Destructor(self)
    !% Default constructor for the \cite{boylan-kolchin_dynamical_2008} merging timescale class.
    implicit none
    type(satelliteMergingTimescalesBoylanKolchin2008), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine boylanKolchin2008Destructor

  double precision function boylanKolchin2008TimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{boylan-kolchin_dynamical_2008} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits
    use Satellite_Orbits
    implicit none
    class           (satelliteMergingTimescalesBoylanKolchin2008)           , intent(inout)          :: self
    type            (treeNode                                   )           , intent(inout), pointer :: thisNode
    type            (keplerOrbit                                )           , intent(inout)          :: thisOrbit
    type            (treeNode                                   )                          , pointer :: hostNode
    class           (nodeComponentBasic                         )                          , pointer :: hostBasic                            , thisBasic
    logical                                                      , parameter                         :: acceptUnboundOrbits          =.false.
    double precision                                             , parameter                         :: A                            =0.216d0, b                 =1.3d0, &  !   Fitting parameters from eqn. (6) of Boylan-Kolchin et al.
         &                                                                                              c                            =1.9d0  , d                 =1.0d0
    double precision                                                                                 :: equivalentCircularOrbitRadius        , massRatio               , &
         &                                                                                              orbitalCircularity                   , radialScale             , &
         &                                                                                              velocityScale

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit)
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute orbital circularity.
    orbitalCircularity                                                                    &
         & =thisOrbit%angularMomentum()                                                   &
         & /equivalentCircularOrbitRadius                                                 &
         & /Dark_Matter_Profile_Circular_Velocity(hostNode,equivalentCircularOrbitRadius)
    ! Compute mass ratio (mass in host [not including satellite] divided by mass in satellite).
    thisBasic => thisNode%basic()
    hostBasic => hostNode%basic()
    massRatio=hostBasic%mass()/thisBasic%mass()-1.0d0
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       boylanKolchin2008TimeUntilMerging=0.0d0
    else
       ! Compute dynamical friction timescale.
       boylanKolchin2008TimeUntilMerging &
            & =Dynamical_Friction_Timescale_Multiplier()           &
            & *Dark_Matter_Halo_Dynamical_Timescale(hostNode)      &
            & *A                                                   &
            & *          massRatio**b                              &
            & /log(1.0d0+massRatio   )                             &
            & *exp(c*orbitalCircularity)                           &
            & *(equivalentCircularOrbitRadius/radialScale)**d
    end if
    return
  end function boylanKolchin2008TimeUntilMerging
