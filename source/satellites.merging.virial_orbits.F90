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

!% Contains a module which implements satellite orbital parameters at virial radius crossing.

module Virial_Orbits
  !% Implements satellite orbital parameters at virial radius crossing.
  use ISO_Varying_String
  use Galacticus_Nodes
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
    use Galacticus_Error
    use Input_Parameters
    use Kepler_Orbits
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
          !# <include directive="virialOrbitsMethod" type="functionCall" functionType="void">
          !#  <functionArgs>virialOrbitsMethod,Virial_Orbital_Parameters_Get</functionArgs>
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
