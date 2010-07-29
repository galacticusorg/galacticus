!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: Virial_Orbital_Parameters

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: virialOrbitsInitialized=.false.

  ! Name of virial overdensity method used.
  type(varying_string)                           :: virialOrbitsMethod

  ! Pointer to the subroutine that returns virial orbital parameters.
  procedure(Virial_Orbital_Parameters), pointer :: Virial_Orbital_Parameters_Get => null()
 
  ! Orbital energy - used for finding radius of equivalent circular orbit.
  double precision                              :: orbitalEnergyInternal
  !$omp threadprivate(orbitalEnergyInternal)

  ! Node used in root finding calculations.
  type(treeNode),                       pointer :: hostNode
  !$omp threadprivate(hostNode)

contains

  subroutine Virial_Orbital_Parameters(thisNode,acceptUnboundOrbits,velocityRadial,velocityTangential,angularMomentum &
       &,orbitalEnergy,eccentricity,semimajorAxis,equivalentCircularOrbitRadius)
    !% Returns virial orbital parameters.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use Input_Parameters
    use Root_Finder
    use FGSL
    use Dark_Matter_Halo_Scales
    !# <include directive="virialOrbitsMethod" type="moduleUse">
    include 'satellites.merging.virial_orbits.modules.inc'
    !# </include>
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    logical,                 intent(in)              :: acceptUnboundOrbits
    double precision,        intent(out),   optional :: velocityRadial,velocityTangential,angularMomentum,orbitalEnergy&
         &,eccentricity ,semimajorAxis,equivalentCircularOrbitRadius
    type(fgsl_function),     save                    :: rootFunction
    type(fgsl_root_fsolver), save                    :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision,        parameter               :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6
    type(c_ptr)                                      :: parameterPointer
    double precision                                 :: velocityRadialInternal,velocityTangentialInternal,angularMomentumInternal&
         & ,eccentricityInternal,semimajorAxisInternal,radiusMinimum,radiusMaximum
    
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

    ! Call the routine to get the orbital parameters.
    call Virial_Orbital_Parameters_Get(thisNode,acceptUnboundOrbits,velocityRadialInternal,velocityTangentialInternal&
         &,angularMomentumInternal ,orbitalEnergyInternal,eccentricityInternal,semimajorAxisInternal)
    if (present(velocityRadial    )) velocityRadial    =velocityRadialInternal
    if (present(velocityTangential)) velocityTangential=velocityTangentialInternal
    if (present(angularMomentum   )) angularMomentum   =angularMomentumInternal
    if (present(orbitalEnergy     )) orbitalEnergy     =orbitalEnergyInternal
    if (present(eccentricity      )) eccentricity      =eccentricityInternal
    if (present(semimajorAxis     )) semimajorAxis     =semimajorAxisInternal

    ! If the radius of the equivalent circular orbit was requested then search for it.
    if (present(equivalentCircularOrbitRadius)) then
       if (orbitalEnergyInternal >= 0.0d0) then
          ! Orbit is unbound, return unphysical value.
          equivalentCircularOrbitRadius=-1.0d0
       else
          hostNode => thisNode%parentNode
          radiusMinimum=Dark_Matter_Halo_Virial_Radius(hostNode)
          radiusMaximum=radiusMinimum
          do while (Equivalent_Orbit_Solver(radiusMinimum,parameterPointer) >= 0.0d0)
             radiusMinimum=0.5d0*radiusMinimum
          end do
          do while (Equivalent_Orbit_Solver(radiusMaximum,parameterPointer) <= 0.0d0)
             radiusMaximum=2.0d0*radiusMaximum
          end do
          equivalentCircularOrbitRadius=Root_Find(radiusMinimum,radiusMaximum,Equivalent_Orbit_Solver,parameterPointer &
               &,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
       end if
    end if    
    return
  end subroutine Virial_Orbital_Parameters
 
  function Equivalent_Orbit_Solver(radius,parameterPointer) bind(c)
    !% Root function used in finding equivalent circular orbits.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    implicit none
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer
    real(c_double)        :: Equivalent_Orbit_Solver

    Equivalent_Orbit_Solver=Dark_Matter_Profile_Potential(hostNode,radius)+0.5d0&
         &*Dark_Matter_Profile_Circular_Velocity(hostNode,radius)**2-orbitalEnergyInternal
    return
  end function Equivalent_Orbit_Solver

end module Virial_Orbits
