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

!% Contains a module which implements a model of the tidal field acting on a satellite assuming spherical symmetry in the host.

module Satellites_Tidal_Fields_Spherical_Symmetry
  !% Implements a module which implements a model of the tidal field acting on a satellite assuming spherical symmetry in the host.
  use Galacticus_Nodes
  implicit none
  private
  public :: Satellites_Tidal_Fields_Spherical_Symmetry_Initialize

  ! Boost factor for the tidal field strength.
  double precision :: satelliteTidalFieldBoostFactor

contains

  !# <satellitesTidalFieldMethod>
  !#  <unitName>Satellites_Tidal_Fields_Spherical_Symmetry_Initialize</unitName>
  !# </satellitesTidalFieldMethod>
  subroutine Satellites_Tidal_Fields_Spherical_Symmetry_Initialize(satellitesTidalFieldMethod,Satellites_Tidal_Field_Get)
    !% Initializes the ``spherical symmetry'' satellite tidal field module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                ), intent(in   )          :: satellitesTidalFieldMethod
    procedure(Satellites_Tidal_Fields_Spherical_Symmetry_Get), intent(inout), pointer :: Satellites_Tidal_Field_Get

    if (satellitesTidalFieldMethod == 'sphericalSymmetry') then
       Satellites_Tidal_Field_Get => Satellites_Tidal_Fields_Spherical_Symmetry_Get
       !@ <inputParameter>
       !@   <name>satelliteTidalFieldBoostFactor</name>
       !@   <defaultValue>1.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which to boost satellite tidal fields in the {\tt sphericalSymmetry} tidal field method.
       !@   </description>
       !@   <type>float</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteTidalFieldBoostFactor',satelliteTidalFieldBoostFactor,defaultValue=1.0d0)
    end if
    return
  end subroutine Satellites_Tidal_Fields_Spherical_Symmetry_Initialize

  double precision function Satellites_Tidal_Fields_Spherical_Symmetry_Get(thisNode)
    !% Computes the tidal field acting on a satellite assuming a spherically symmetric host.
    use Galacticus_Nodes
    use Kepler_Orbits
    use Satellite_Orbits
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    type            (treeNode              )               , pointer :: hostNode
    class           (nodeComponentSatellite)               , pointer :: thisSatellite
    type            (keplerOrbit           )                         :: thisOrbit
    double precision                                                 :: densityHost  , enclosedMassHost, orbitalRadius, orbitalVelocity

    ! For isolated halos, always return zero tidal field.
    if (thisNode%isSatellite()) then
       ! Find the host node.
       hostNode      => thisNode     %parent
       ! Get the satellite component.
       thisSatellite => thisNode     %satellite  ()
       ! Get the orbit for this node.
       thisOrbit     =  thisSatellite%virialOrbit()
       ! Get the orbital radius and velocity at pericenter.
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(hostNode,thisOrbit,extremumPericenter,orbitalRadius,orbitalVelocity)
       ! Find the mass and density of the host halo at pericenter.
       densityHost     =Galactic_Structure_Density      (                                              &
            &                                            hostNode                                    , &
            &                                            [orbitalRadius,0.0d0,0.0d0]                 , &
            &                                            coordinateSystem=coordinateSystemCylindrical  &
            &                                           )
       enclosedMassHost=Galactic_Structure_Enclosed_Mass(                                              &
            &                                            hostNode                                    , &
            &                                            orbitalRadius                                 &
            &                                           )
       ! Compute the tidal field.
       Satellites_Tidal_Fields_Spherical_Symmetry_Get=                                                        &
            &             gravitationalConstantGalacticus*enclosedMassHost/                 orbitalRadius **3 &
            &   -4.0d0*Pi*gravitationalConstantGalacticus*densityHost                                         &
            &   +                                                          (orbitalVelocity/orbitalRadius)**2
       ! Boost the tidal field.
       Satellites_Tidal_Fields_Spherical_Symmetry_Get=satelliteTidalFieldBoostFactor*Satellites_Tidal_Fields_Spherical_Symmetry_Get
    else
       Satellites_Tidal_Fields_Spherical_Symmetry_Get=0.0d0
    end if
    return
  end function Satellites_Tidal_Fields_Spherical_Symmetry_Get

end module Satellites_Tidal_Fields_Spherical_Symmetry
