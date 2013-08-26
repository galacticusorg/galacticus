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

!% Contains a module which implements calculations related to satellite orbits.

module Satellite_Orbits
  !% Implements calculations related to satellite orbits.
  use Galacticus_Nodes
  implicit none
  private
  public :: Satellite_Orbit_Equivalent_Circular_Orbit_Radius, Satellite_Orbit_Extremum_Phase_Space_Coordinates

  ! Orbital energy and angular momentum - used for finding radius of equivalent circular orbit.
  double precision                    :: orbitalAngularMomentumInternal, orbitalEnergyInternal
  !$omp threadprivate(orbitalEnergyInternal,orbitalAngularMomentumInternal)
  ! Node used in root finding calculations.
  type            (treeNode), pointer :: activeNode
  !$omp threadprivate(activeNode)

  ! Enumeratation used to indicate type of extremum.
  integer, public, parameter :: extremumPericenter=-1
  integer, public, parameter :: extremumApocenter =+1

contains

  double precision function Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit)
    !% Solves for the equivalent circular orbit radius for {\tt thisOrbit} in {\tt hostNode}.
    use Root_Finder
    use Kepler_Orbits
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode         ), intent(inout), pointer :: hostNode
    type            (keplerOrbit      ), intent(inout)          :: thisOrbit
    double precision                   , parameter              :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder       ), save                   :: finder
    !$omp threadprivate(finder)
    type            (keplerOrbit      )                         :: currentOrbit

    ! Convert the orbit to the potential of the current halo in which the satellite finds itself.
    currentOrbit=Satellite_Orbit_Convert_To_Current_Potential(thisOrbit,hostNode)

    orbitalEnergyInternal=thisOrbit%energy()
    if (orbitalEnergyInternal >= 0.0d0) then
       ! Orbit is unbound, return unphysical value.
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=-1.0d0
    else
       activeNode => hostNode
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(Equivalent_Circular_Orbit_Solver   )
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
       end if
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=finder%find(rootGuess=Dark_Matter_Halo_Virial_Radius(hostNode))
    end if
    return
  end function Satellite_Orbit_Equivalent_Circular_Orbit_Radius

  double precision function Equivalent_Circular_Orbit_Solver(radius)
    !% Root function used in finding equivalent circular orbits.
    use Dark_Matter_Profiles
    implicit none
    double precision, intent(in   ) :: radius

    Equivalent_Circular_Orbit_Solver=Dark_Matter_Profile_Potential(activeNode,radius)+0.5d0&
         &*Dark_Matter_Profile_Circular_Velocity(activeNode,radius)**2-orbitalEnergyInternal
    return
  end function Equivalent_Circular_Orbit_Solver

  subroutine Satellite_Orbit_Extremum_Phase_Space_Coordinates(hostNode,thisOrbit,extremumType,radius,velocity)
    !% Solves for the pericentric radius and velocity of {\tt thisOrbit} in {\tt hostNode}.
    use Root_Finder
    use Kepler_Orbits
    implicit none
    type            (treeNode         ), intent(inout), pointer :: hostNode
    type            (keplerOrbit      ), intent(inout)          :: thisOrbit
    integer                            , intent(in   )          :: extremumType
    double precision                   , intent(  out)          :: radius                 , velocity
    double precision                   , parameter              :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder       ), save                   :: finder
    !$omp threadprivate(finder)
    type            (keplerOrbit      )                         :: currentOrbit

    ! Convert the orbit to the potential of the current halo in which the satellite finds itself.
    currentOrbit=Satellite_Orbit_Convert_To_Current_Potential(thisOrbit,hostNode)
    ! Extract the orbital energy and angular momentum.
    orbitalEnergyInternal         =  currentOrbit%energy         ()
    orbitalAngularMomentumInternal=  currentOrbit%angularMomentum()
    ! Set a pointer to the host node.
    activeNode                    => hostNode
    ! Catch orbits which are close to being circular.
    if      (   Extremum_Solver(currentOrbit%radius()) == 0.0d0             ) then
       ! Orbit is at extremum.
       radius=currentOrbit%radius()
    else if (&
         &    (&
         &      extremumType                           == extremumPericenter       &
         &     .and.                                                               &
         &      Extremum_Solver(currentOrbit%radius()) >  0.0d0                    &
         &    )                                                                    &
         &   .or.                                                                  &
         &    (                                                                    &
         &      extremumType                           == extremumApocenter        &
         &     .and.                                                               &
         &      Extremum_Solver(currentOrbit%radius()) <  0.0d0                    &
         &    )                                                                    &
         &  ) then
       ! No solution exists, assume a circular orbit.
       radius=currentOrbit%radius()
    else
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(Extremum_Solver                    )
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          select case (extremumType)
          case (extremumPericenter)
             call finder%rangeExpand (                                                             &
                  &                   rangeExpandDownward          =0.5d0                        , &
                  &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                  &                   rangeExpandType              =rangeExpandMultiplicative      &
                  &                  )
          case (extremumApocenter )
             call finder%rangeExpand (                                                             &
                  &                   rangeExpandUpward            =2.0d0                        , &
                  &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
                  &                   rangeExpandType              =rangeExpandMultiplicative      &
                  &                  )
          end select
       end if
       radius=finder%find(rootGuess=currentOrbit%radius())
    end if
    ! Get the orbital velocity at this radius.
    velocity=orbitalAngularMomentumInternal/radius
    return
  end subroutine Satellite_Orbit_Extremum_Phase_Space_Coordinates
  
  double precision function Extremum_Solver(radius)
    !% Root function used in finding orbital extremum radius.
    use Dark_Matter_Profiles
    implicit none
    double precision, intent(in   ) :: radius

    Extremum_Solver=Dark_Matter_Profile_Potential(activeNode,radius)+0.5d0*(orbitalAngularMomentumInternal/radius)**2-orbitalEnergyInternal
    return
  end function Extremum_Solver

  function Satellite_Orbit_Convert_To_Current_Potential(thisOrbit,currentHost)
    !% Takes a virial orbit and adjusts the energy to account for the change in the definition of potential between the original
    !% halo in which the orbit was defined and the current halo. Since the potential at the virial radius of halos is always
    !% defined to be $\Phi(r_{\rm vir}) = - V_{\rm vir}^2$ then the specific energy transforms as:
    !% \begin{equation}
    !% e \rightarrow e + V^2_{\rm vir,0} + \Phi(r_{\rm vir,0}),
    !% \end{equation}
    !% where subscript $0$ refers to the original halo in which the orbit was defined and $\Phi(r)$ is the potential of the
    !% current halo.
    use Galactic_Structure_Potentials
    use Numerical_Constants_Physical
    use Kepler_Orbits
    implicit none
    type            (keplerOrbit)                         :: Satellite_Orbit_Convert_To_Current_Potential
    type            (keplerOrbit), intent(inout)          :: thisOrbit
    type            (treeNode   ), intent(inout), pointer :: currentHost
    double precision                                      :: potentialHost                               , radiusVirialOriginal, &
         &                                                   velocityVirialOriginal

    ! Compute the properties of the initial orbit, and the current potential.
    radiusVirialOriginal  =gravitationalConstantGalacticus*thisOrbit%hostMass()/thisOrbit%velocityScale()**2
    velocityVirialOriginal=                                                     thisOrbit%velocityScale()
    potentialHost         =Galactic_Structure_Potential(currentHost,radiusVirialOriginal)
    ! Create a new orbit with an adjusted energy.
    Satellite_Orbit_Convert_To_Current_Potential=thisOrbit
    call Satellite_Orbit_Convert_To_Current_Potential%energySet(thisOrbit%energy()+velocityVirialOriginal**2+potentialHost)
    return
  end function Satellite_Orbit_Convert_To_Current_Potential

end module Satellite_Orbits
