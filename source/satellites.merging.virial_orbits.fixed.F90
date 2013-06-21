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

!% Contains a module which implements a fixed orbital parameter distribution for merging subhalos.

module Virial_Orbits_Fixed
  !% Implements a fixed orbital parameter distribution for merging subhalos.
  implicit none
  private
  public :: Virial_Orbital_Parameters_Fixed_Initialize

  ! Fixed radial and tangential velocities to use (in units of host node virial velocity).
  double precision :: virialOrbitsFixedRadialVelocity, virialOrbitsFixedTangentialVelocity

contains

  !# <virialOrbitsMethod>
  !#  <unitName>Virial_Orbital_Parameters_Fixed_Initialize</unitName>
  !# </virialOrbitsMethod>
  subroutine Virial_Orbital_Parameters_Fixed_Initialize(virialOrbitsMethod,Virial_Orbital_Parameters_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Kepler_Orbits
    implicit none
    type     (varying_string                 ), intent(in   )          :: virialOrbitsMethod
    procedure(Virial_Orbital_Parameters_Fixed), intent(inout), pointer :: Virial_Orbital_Parameters_Get

    if (virialOrbitsMethod == 'fixed') then
       ! Return a pointer to our implementation.
       Virial_Orbital_Parameters_Get => Virial_Orbital_Parameters_Fixed
       ! Read parameters of the fixed orbits.
       !@ <inputParameter>
       !@   <name>virialOrbitsFixedRadialVelocity</name>
       !@   <defaultValue>0.90</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('virialOrbitsFixedRadialVelocity'    ,virialOrbitsFixedRadialVelocity    ,defaultValue=-0.90d0)
       !@ <inputParameter>
       !@   <name>virialOrbitsFixedTangentialVelocity</name>
       !@   <defaultValue>0.75</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('virialOrbitsFixedTangentialVelocity',virialOrbitsFixedTangentialVelocity,defaultValue= 0.75d0)
    end if
    return
  end subroutine Virial_Orbital_Parameters_Fixed_Initialize

  function Virial_Orbital_Parameters_Fixed(thisNode,hostNode,acceptUnboundOrbits) result (thisOrbit)
    !% Return fixed orbital parameters for a satellite.
    use Kepler_Orbits
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type            (keplerOrbit       )                         :: thisOrbit
    type            (treeNode          ), intent(inout), pointer :: hostNode           , thisNode
    logical                             , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentBasic)               , pointer :: hostBasicComponent , thisBasicComponent
    double precision                                             :: velocityScale

    ! Reset the orbit.
    call thisOrbit%reset()
    ! Set masses and radius of the orbit.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    call thisOrbit%massesSet(thisBasicComponent%mass(),hostBasicComponent%mass())
    call thisOrbit%radiusSet(Dark_Matter_Halo_Virial_Radius(hostNode))
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    call thisOrbit%velocityRadialSet    (virialOrbitsFixedRadialVelocity    *velocityScale)
    call thisOrbit%velocityTangentialSet(virialOrbitsFixedTangentialVelocity*velocityScale)
    return
  end function Virial_Orbital_Parameters_Fixed

end module Virial_Orbits_Fixed
