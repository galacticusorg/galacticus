!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a model of ram pressure stripping of hot halos based on the methods of
!% \cite{font_colours_2008}.

module Hot_Halo_Ram_Pressure_Force_Font2008
  !% Implements a module which implements the calculation of the ram pressure force on hot halos based on the methods of
  !% \cite{font_colours_2008}.
  use Galacticus_Nodes
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Force_Font2008_Initialize
  
contains

  !# <hotHaloRamPressureForceMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Force_Font2008_Initialize</unitName>
  !# </hotHaloRamPressureForceMethod>
  subroutine Hot_Halo_Ram_Pressure_Force_Font2008_Initialize(hotHaloRamPressureForceMethod,Hot_Halo_Ram_Pressure_Force_Get)
    !% Initializes the ``Font2008'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                          ), intent(in   )          :: hotHaloRamPressureForceMethod
    procedure(Hot_Halo_Ram_Pressure_Force_Font2008_Get), intent(inout), pointer :: Hot_Halo_Ram_Pressure_Force_Get

    if (hotHaloRamPressureForceMethod == 'Font2008') then
       Hot_Halo_Ram_Pressure_Force_Get => Hot_Halo_Ram_Pressure_Force_Font2008_Get
    end if
    return
  end subroutine Hot_Halo_Ram_Pressure_Force_Font2008_Initialize

  double precision function Hot_Halo_Ram_Pressure_Force_Font2008_Get(thisNode)
    !% Computes the hot halo ram pressure force
    use Kepler_Orbits
    use Satellite_Orbits
    use Hot_Halo_Mass_Distributions
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: thisNode
    class           (nodeComponentSatellite      )               , pointer :: thisSatellite

    type            (treeNode                    )               , pointer :: hostNode
    type            (keplerOrbit                 )                         :: thisOrbit
    class           (hotHaloMassDistributionClass)               , pointer :: hotHaloMassDistribution_
    double precision                                                       :: densityHotHaloHost      , orbitalRadius, orbitalVelocity
    
    ! Find the host node.
    hostNode      => thisNode%parent
    ! Get the satellite component.
    thisSatellite => thisNode%satellite()
    ! Get the orbit for this node.
    thisOrbit=thisSatellite%virialOrbit()
    ! Get the orbital radius and velocity at pericenter.
    call Satellite_Orbit_Extremum_Phase_Space_Coordinates(hostNode,thisOrbit,extremumPericenter,orbitalRadius,orbitalVelocity)
    ! Get the hot halo mass distribution.
    hotHaloMassDistribution_ => hotHaloMassDistribution()
    ! Find the density of the host node hot halo at the pericentric radius.
    densityHotHaloHost=hotHaloMassDistribution_%density(hostNode,orbitalRadius)
    ! Find the ram pressure force at pericenter.
    Hot_Halo_Ram_Pressure_Force_Font2008_Get=densityHotHaloHost*orbitalVelocity**2
    return
  end function Hot_Halo_Ram_Pressure_Force_Font2008_Get

end module Hot_Halo_Ram_Pressure_Force_Font2008
