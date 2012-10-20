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

!% Contains a module which implements calculations related to satellite orbits.

module Satellite_Orbits
  !% Implements calculations related to satellite orbits.
  use Tree_Nodes
  implicit none
  private
  public :: Satellite_Orbit_Equivalent_Circular_Orbit_Radius, Satellite_Orbit_Pericenter_Phase_Space_Coordinates
  
  ! Orbital energy and angular momentum - used for finding radius of equivalent circular orbit.
  double precision          :: orbitalEnergyInternal,orbitalAngularMomentumInternal
  !$omp threadprivate(orbitalEnergyInternal,orbitalAngularMomentumInternal)
  
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

  subroutine Satellite_Orbit_Pericenter_Phase_Space_Coordinates(hostNode,thisOrbit,radius,velocity)
    !% Solves for the pericentric radius and velocity of {\tt thisOrbit} in {\tt hostNode}.
    use, intrinsic :: ISO_C_Binding
    use Root_Finder
    use FGSL
    use Kepler_Orbits_Structure
    implicit none
    type(treeNode),          pointer, intent(inout) :: hostNode
    type(keplerOrbit),                intent(inout) :: thisOrbit
    double precision,                 intent(  out) :: radius,velocity
    double precision,        parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6
    type(fgsl_function),     save                   :: rootFunction
    type(fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    type(c_ptr)                                     :: parameterPointer
    double precision                                :: radiusMinimum,radiusMaximum

    ! Extract the orbital energy and angular momentum.
    orbitalEnergyInternal         =  thisOrbit%energy         ()
    orbitalAngularMomentumInternal=  thisOrbit%angularMomentum()
    ! Set a pointer to the host node.
    activeNode                    => hostNode
    ! Find the radial range within which the pericenter must lie. The pericenter must be smaller than (or equal to) the current radius of the orbit.
    radiusMinimum=thisOrbit%radius()
    radiusMaximum=radiusMinimum

    ! Catch orbits which are close to being circular.
    if      (Pericenter_Solver(radiusMinimum,parameterPointer) == 0.0d0) then
       ! Orbit is at pericenter.
       radius=radiusMinimum
    else if (Pericenter_Solver(radiusMinimum,parameterPointer) >  0.0d0) then
       ! No solution exists, assume a circular orbit.
       radius=radiusMinimum
    else
       ! Find the pericenter of the orbit.
       do while (Pericenter_Solver(radiusMinimum,parameterPointer) <= 0.0d0)
          radiusMinimum=0.5d0*radiusMinimum
       end do
       ! Solve for the pericentric radius.
       radius=Root_Find(radiusMinimum,radiusMaximum,Pericenter_Solver,parameterPointer,rootFunction,rootFunctionSolver&
            &,toleranceAbsolute,toleranceRelative)
    end if
    ! Get the orbital velocity at this radius.
    velocity=orbitalAngularMomentumInternal/radius
    return
  end subroutine Satellite_Orbit_Pericenter_Phase_Space_Coordinates

  function Pericenter_Solver(radius,parameterPointer) bind(c)
    !% Root function used in finding orbital pericentric radius.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    implicit none
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer
    real(c_double)        :: Pericenter_Solver

    Pericenter_Solver=Dark_Matter_Profile_Potential(activeNode,radius)+0.5d0*(orbitalAngularMomentumInternal/radius)**2-orbitalEnergyInternal
    return
  end function Pericenter_Solver

end module Satellite_Orbits
