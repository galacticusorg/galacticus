!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements satellite orbital parameters at virial radius crossing.

module Virial_Orbits
  !% Implements satellite orbital parameters at virial radius crossing.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Virial_Orbital_Parameters

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: virialOrbitsInitialized=.false.

  ! Name of virial overdensity method used.
  type(varying_string)                           :: virialOrbitsMethod

  ! Pointer to the function that returns virial orbital parameters.
  procedure(Virial_Orbital_Parameters), pointer :: Virial_Orbital_Parameters_Get => null()
 
contains

  function Virial_Orbital_Parameters(thisNode,hostNode,acceptUnboundOrbits) result (thisOrbit)
    !% Returns virial orbital parameters.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use Input_Parameters
    use Root_Finder
    use FGSL
    use Dark_Matter_Halo_Scales
    use Kepler_Orbits_Structure
    !# <include directive="virialOrbitsMethod" type="moduleUse">
    include 'satellites.merging.virial_orbits.modules.inc'
    !# </include>
    implicit none
    type(keplerOrbit)                         :: thisOrbit
    type(treeNode),   intent(inout), pointer  :: thisNode,hostNode
    logical,          intent(in)              :: acceptUnboundOrbits
    
    if (.not.virialOrbitsInitialized) then
       !$omp critical(virialOrbitsInitialized)
       if (.not.virialOrbitsInitialized) then
          ! Get the virial orbits method parameter.
          !@ <inputParameter>
          !@   <name>virialOrbitsMethod</name>
          !@   <defaultValue>Benson2005</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects the method to be used for finding orbital parameters of satellites at virial radius crossing.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialOrbitsMethod',virialOrbitsMethod,defaultValue='Benson2005')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="virialOrbitsMethod" type="code" action="subroutine">
          !#  <subroutineArgs>virialOrbitsMethod,Virial_Orbital_Parameters_Get</subroutineArgs>
          include 'satellites.merging.virial_orbits.inc'
          !# </include>
          if (.not.associated(Virial_Orbital_Parameters_Get)) call Galacticus_Error_Report('Virial_Orbital_Parameters','method ' &
               &//char(virialOrbitsMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          virialOrbitsInitialized=.true.
       end if
       !$omp end critical(virialOrbitsInitialized)
    end if

    ! Call the routine to get the orbital parameters.
    thisOrbit=Virial_Orbital_Parameters_Get(thisNode,hostNode,acceptUnboundOrbits)

    return
  end function Virial_Orbital_Parameters

end module Virial_Orbits
