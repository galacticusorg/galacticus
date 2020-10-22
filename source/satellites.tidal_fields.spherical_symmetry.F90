!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !# <satelliteTidalField name="satelliteTidalFieldSphericalSymmetry">
  !#  <description>
  !#   A satellite tidal field class which assumes a spherically-symmetric host halo, and computes the tidal field accordingly
  !#   using:
  !#   \begin{equation}
  !#    \mathcal{F} = {\mathrm{G} M_\mathrm{host}(<r_\mathrm{p}) \over r_\mathrm{p}^3} - 4 \pi \mathrm{G}
  !#    \rho_\mathrm{host}(r_\mathrm{p}) + \omega_\mathrm{p}^2,
  !#   \end{equation}
  !#   where $r_\mathrm{p}$ is the pericentric radius. $M_\mathrm{host}(<r)$ is the mass of the host halo enclosed within a sphere
  !#   of radius $r$, $\rho_\mathrm{host}(r)$ is the host density at radius $r$, and $\omega_\mathrm{p}$ is the orbital angular
  !#   velocity at pericenter.
  !#  </description>
  !# </satelliteTidalField>
  type, extends(satelliteTidalFieldClass) :: satelliteTidalFieldSphericalSymmetry
     !% Implementation of a satellite tidal friction class which assumes spherical symmetry.
     private
     double precision :: factorBoost
   contains
     procedure :: tidalTensorRadial => sphericalSymmetryTidalTensorRadial
  end type satelliteTidalFieldSphericalSymmetry

  interface satelliteTidalFieldSphericalSymmetry
     !% Constructors for the sphericalSymmetry satellite tidal field class.
     module procedure sphericalSymmetryConstructorParameters
     module procedure sphericalSymmetryConstructorInternal
  end interface satelliteTidalFieldSphericalSymmetry

contains

  function sphericalSymmetryConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalSymmetry} satellite tidal field class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalFieldSphericalSymmetry)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    double precision                                                      :: factorBoost

    !# <inputParameter>
    !#   <name>factorBoost</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The factor by which to boost satellite tidal fields in the {\normalfont \ttfamily sphericalSymmetry} tidal field class.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=satelliteTidalFieldSphericalSymmetry(factorBoost)
    !# <inputParametersValidate source="parameters"/>
    return
  end function sphericalSymmetryConstructorParameters

  function sphericalSymmetryConstructorInternal(factorBoost) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalSymmetry} satellite tidal field class.
    implicit none
    type            (satelliteTidalFieldSphericalSymmetry)            :: self
    double precision                                     , intent(in) :: factorBoost
    !# <constructorAssign variables="factorBoost"/>

    return
  end function sphericalSymmetryConstructorInternal

  double precision function sphericalSymmetryTidalTensorRadial(self,node)
    !% Return the radial part of the tidal tensor for satellite halos assuming spherical symmetry of the host.
    use :: Galactic_Structure_Densities      , only : Galactic_Structure_Density
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : coordinateSystemCylindrical
    use :: Galacticus_Nodes                  , only : nodeComponentSatellite                          , treeNode
    use :: Kepler_Orbits                     , only : keplerOrbit
    use :: Numerical_Constants_Math          , only : Pi
    use :: Numerical_Constants_Astronomical      , only : gravitationalConstantGalacticus
    use :: Satellite_Orbits                  , only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumPericenter
    implicit none
    class           (satelliteTidalFieldSphericalSymmetry), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    type            (treeNode                            ), pointer       :: nodeHost
    class           (nodeComponentSatellite              ), pointer       :: satellite
    type            (keplerOrbit                         )                :: orbit
    double precision                                                      :: densityHost  , enclosedMassHost, &
         &                                                                   radiusOrbital, velocityOrbital

    ! For isolated halos, always return zero tidal field.
    if (node%isSatellite()) then
       ! Find the host node.
       nodeHost  => node     %parent
       ! Get the satellite component.
       satellite => node     %satellite  ()
       ! Get the orbit for this node.
       orbit     =  satellite%virialOrbit()
       ! Get the orbital radius and velocity at pericenter.
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost,orbit,extremumPericenter,radiusOrbital,velocityOrbital)
       ! Find the mass and density of the host halo at pericenter.
       densityHost     =Galactic_Structure_Density      (                                              &
            &                                            nodeHost                                    , &
            &                                            [radiusOrbital,0.0d0,0.0d0]                 , &
            &                                            coordinateSystem=coordinateSystemCylindrical  &
            &                                           )
       enclosedMassHost=Galactic_Structure_Enclosed_Mass(                                              &
            &                                            nodeHost                                    , &
            &                                            radiusOrbital                                 &
            &                                           )
       ! Compute the tidal field.
       sphericalSymmetryTidalTensorRadial=+         gravitationalConstantGalacticus*enclosedMassHost/                 radiusOrbital **3 &
            &                             -4.0d0*Pi*gravitationalConstantGalacticus*densityHost                                         &
            &                             +                                                          (velocityOrbital/radiusOrbital)**2
       ! Boost the tidal field.
       sphericalSymmetryTidalTensorRadial=+self%factorBoost                   &
            &                             *sphericalSymmetryTidalTensorRadial
    else
       sphericalSymmetryTidalTensorRadial=+0.0d0
    end if
    return
  end function sphericalSymmetryTidalTensorRadial
