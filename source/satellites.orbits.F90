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


!% Contains a module which implements calculations related to satellite orbits.

module Satellite_Orbits
  !% Implements calculations related to satellite orbits.
  use Tree_Nodes
  implicit none
  private
  public :: Satellite_Orbit_Equivalent_Circular_Orbit_Radius
  
  ! Orbital energy - used for finding radius of equivalent circular orbit.
  double precision          :: orbitalEnergyInternal
  !$omp threadprivate(orbitalEnergyInternal)
  
  ! Node used in root finding calculations.
  type(treeNode),   pointer :: activeNode
  !$omp threadprivate(activeNode)

contains
  
  double precision function Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit)
    !% Solves for the equivalent circular orbit radius for {\tt thisOrbit} in {\tt hostNode}.
    use, intrinsic :: ISO_C_Binding
    use Root_Finder
    use FGSL
    use Kepler_Orbits_Structure
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),          pointer, intent(inout) :: hostNode
    type(keplerOrbit),                intent(inout) :: thisOrbit
    double precision,        parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6
    type(fgsl_function),     save                   :: rootFunction
    type(fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    type(c_ptr)                                     :: parameterPointer
    double precision                                :: radiusMinimum,radiusMaximum

    orbitalEnergyInternal=thisOrbit%energy()
    if (orbitalEnergyInternal >= 0.0d0) then
       ! Orbit is unbound, return unphysical value.
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=-1.0d0
    else
       activeNode => hostNode
       radiusMinimum=Dark_Matter_Halo_Virial_Radius(hostNode)
       radiusMaximum=radiusMinimum
       do while (Equivalent_Circular_Orbit_Solver(radiusMinimum,parameterPointer) >= 0.0d0)
          radiusMinimum=0.5d0*radiusMinimum
       end do
       do while (Equivalent_Circular_Orbit_Solver(radiusMaximum,parameterPointer) <= 0.0d0)
          radiusMaximum=2.0d0*radiusMaximum
       end do
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=Root_Find(radiusMinimum,radiusMaximum,Equivalent_Circular_Orbit_Solver&
            &,parameterPointer ,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
    end if
    return
  end function Satellite_Orbit_Equivalent_Circular_Orbit_Radius

  function Equivalent_Circular_Orbit_Solver(radius,parameterPointer) bind(c)
    !% Root function used in finding equivalent circular orbits.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    implicit none
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer
    real(c_double)        :: Equivalent_Circular_Orbit_Solver

    Equivalent_Circular_Orbit_Solver=Dark_Matter_Profile_Potential(activeNode,radius)+0.5d0&
         &*Dark_Matter_Profile_Circular_Velocity(activeNode,radius)**2-orbitalEnergyInternal
    return
  end function Equivalent_Circular_Orbit_Solver

end module Satellite_Orbits
