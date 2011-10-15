!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a fixed orbital parameter distribution for merging subhalos.

module Virial_Orbits_Fixed
  !% Implements a fixed orbital parameter distribution for merging subhalos.
  implicit none
  private
  public :: Virial_Orbital_Parameters_Fixed_Initialize
  
  ! Fixed radial and tangential velocities to use (in units of host node virial velocity).
  double precision :: virialOrbitsFixedRadialVelocity,virialOrbitsFixedTangentialVelocity

contains

  !# <virialOrbitsMethod>
  !#  <unitName>Virial_Orbital_Parameters_Fixed_Initialize</unitName>
  !# </virialOrbitsMethod>
  subroutine Virial_Orbital_Parameters_Fixed_Initialize(virialOrbitsMethod,Virial_Orbital_Parameters_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Kepler_Orbits_Structure
    implicit none
    type(varying_string),                  intent(in)    :: virialOrbitsMethod
    procedure(type(keplerOrbit)), pointer, intent(inout) :: Virial_Orbital_Parameters_Get
    
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
    use Kepler_Orbits_Structure
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(keplerOrbit)                         :: thisOrbit
    type(treeNode),   intent(inout), pointer  :: thisNode,hostNode
    logical,          intent(in)              :: acceptUnboundOrbits
    double precision                          :: velocityScale

    ! Reset the orbit.
    call thisOrbit%reset()
    ! Set masses and radius of the orbit.
    call thisOrbit%massesSet(Tree_Node_Mass(thisNode),Tree_Node_Mass(hostNode))
    call thisOrbit%radiusSet(Dark_Matter_Halo_Virial_Radius(hostNode))
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    call thisOrbit%velocityRadialSet    (virialOrbitsFixedRadialVelocity    *velocityScale)
    call thisOrbit%velocityTangentialSet(virialOrbitsFixedTangentialVelocity*velocityScale)
    return
  end function Virial_Orbital_Parameters_Fixed
  
end module Virial_Orbits_Fixed
