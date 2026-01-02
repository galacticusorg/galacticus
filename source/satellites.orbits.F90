!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which implements calculations related to satellite orbits.
!!}

module Satellite_Orbits
  !!{
  Implements calculations related to satellite orbits.
  !!}
  use :: Galacticus_Nodes  , only : treeNode
  use :: Kind_Numbers      , only : kind_int8
  use :: Mass_Distributions, only : massDistributionClass
  implicit none
  private
  public :: Satellite_Orbit_Equivalent_Circular_Orbit_Radius, Satellite_Orbit_Extremum_Phase_Space_Coordinates

  ! Orbital energy and angular momentum - used for finding radius of equivalent circular orbit.
  double precision                                               :: orbitalAngularMomentumInternal        , orbitalEnergyInternal
  !$omp threadprivate(orbitalEnergyInternal,orbitalAngularMomentumInternal)

  ! Objects used in root finding calculations.
  type            (treeNode                 ), pointer           :: activeNode
  class           (massDistributionClass    ), pointer           :: massDistribution__
  double precision                                               :: radiusVirial__                        , massVirial__
  !$omp threadprivate(activeNode,radiusVirial__,massVirial__,massDistribution__)

  ! Enumeration used to indicate type of extremum.
  integer                                    , parameter, public :: extremumPericenter            =-1
  integer                                    , parameter, public :: extremumApocenter             =+1

  ! Error codes.
  integer                                    , parameter, public :: errorCodeSuccess              =0
  integer                                    , parameter, public :: errorCodeOrbitUnbound         =1
  integer                                    , parameter, public :: errorCodeNoEquivalentOrbit    =2

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8           )                    :: lastUniqueID                  =-1
  logical                                                        :: pericenterCalculated          =.false.
  logical                                                        :: apocenterCalculated           =.false.
  double precision                                               :: timePrevious
  double precision                                               :: orbitalEnergyPrevious
  double precision                                               :: orbitalAngularMomentumPrevious
  double precision                                               :: pericenterRadius
  double precision                                               :: pericenterVelocity
  double precision                                               :: apocenterRadius
  double precision                                               :: apocenterVelocity
  !$omp threadprivate(lastUniqueID,pericenterCalculated,apocenterCalculated,timePrevious,orbitalEnergyPrevious,orbitalAngularMomentumPrevious,pericenterRadius,pericenterVelocity,apocenterRadius,apocenterVelocity)

contains

  double precision function Satellite_Orbit_Equivalent_Circular_Orbit_Radius(nodeHost,orbit,darkMatterHaloScale_,errorCode)
    !!{
    Solves for the equivalent circular orbit radius for {\normalfont \ttfamily orbit} in {\normalfont \ttfamily nodeHost}.
    !!}
    use :: Galacticus_Nodes       , only : nodeComponentBasic
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
    use :: Kepler_Orbits          , only : keplerOrbit
    use :: Root_Finder            , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    type            (treeNode                 ), intent(inout), target   :: nodeHost
    type            (keplerOrbit              ), intent(inout)           :: orbit
    integer                                    , intent(  out), optional :: errorCode
    class           (darkMatterHaloScaleClass ), intent(inout)           :: darkMatterHaloScale_
    class           (nodeComponentBasic       )               , pointer  :: basicHost
    double precision                           , parameter               :: toleranceAbsolute    =0.0d0, toleranceRelative=1.0d-6, &
         &                                                                  factorRadiusLarge    =1.0d6
    type            (rootFinder               )                          :: finder
    type            (keplerOrbit              )                          :: orbitCurrent
    double precision                                                     :: potential

    ! Assign the active node.
    activeNode         => nodeHost
    ! Get the mass distribution.
    massDistribution__ => nodeHost%massDistribution()
    ! Convert the orbit to the potential of the current halo in which the satellite finds itself.
    orbitCurrent=Satellite_Orbit_Convert_To_Current_Potential(orbit,nodeHost,darkMatterHaloScale_)
    ! Get virial properties.
    basicHost      => nodeHost             %basic       (        )
    radiusVirial__ =  darkMatterHaloScale_ %radiusVirial(nodeHost)
    massVirial__   =  basicHost            %mass        (        )
    ! Store the orbital energy.
    orbitalEnergyInternal=orbit%energy()
    ! Test for conditions that an equivalent circular orbit exists.
    potential=Satellite_Orbit_Potential(factorRadiusLarge*radiusVirial__,radiusVirial__,massVirial__)
    if (orbitalEnergyInternal >= potential) then
       ! Orbit is unbound, return unphysical value.
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=-1.0d0
       if (present(errorCode)) errorCode=errorCodeOrbitUnbound
    else if (Equivalent_Circular_Orbit_Solver(0.0d0) > 0.0d0) then
       ! No equivalent circular orbit exists (i.e. the orbital energy is less [i.e. more negative] than the gravitational
       ! potential at zero radius. Return an unphysical value.
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=-1.0d0
       if (present(errorCode)) errorCode=errorCodeNoEquivalentOrbit
    else
       finder        =rootFinder(                                                                &
            &                    rootFunction                 =Equivalent_Circular_Orbit_Solver, &
            &                    toleranceAbsolute            =toleranceAbsolute               , &
            &                    toleranceRelative            =toleranceRelative               , &
            &                    rangeExpandUpward            =2.0d0                           , &
            &                    rangeExpandDownward          =0.5d0                           , &
            &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative   , &
            &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive   , &
            &                    rangeExpandType              =rangeExpandMultiplicative         &
            &                   )
       Satellite_Orbit_Equivalent_Circular_Orbit_Radius=finder%find(rootGuess=radiusVirial__)
       if (present(errorCode)) errorCode=errorCodeSuccess
    end if
    !![
    <objectDestructor name="massDistribution__"/>
    !!]
    return
  end function Satellite_Orbit_Equivalent_Circular_Orbit_Radius

  double precision function Equivalent_Circular_Orbit_Solver(radius) result(radiusCircular)
    !!{
    Root function used in finding equivalent circular orbits.
    !!}
    use :: Galactic_Structure_Options, only : enumerationStructureErrorCodeType, structureErrorCodeInfinite, structureErrorCodeSuccess
    use :: Error                     , only : Error_Report
    implicit none
    double precision                                   , intent(in   ) :: radius
    double precision                                   , parameter     :: potentialInfinite=huge(1.0d0)
    double precision                                                   :: potential
    type            (enumerationStructureErrorCodeType)                :: status

    ! Get potential.
    potential=Satellite_Orbit_Potential(radius,radiusVirial__,massVirial__,status)
    select case (status%ID)
    case (structureErrorCodeSuccess %ID)
       radiusCircular=+potential+0.5d0*massDistribution__%rotationCurve(radius)**2-orbitalEnergyInternal
    case (structureErrorCodeInfinite%ID)
       ! The gravitational potential is negative infinity at this radius (most likely zero radius). Since all we care about in
       ! this root-finding function is the sign of the function, return a large negative value.
       radiusCircular=-potentialInfinite
    case default
       radiusCircular=+0.0d0
       call Error_Report('dark matter potential evaluation failed'//{introspection:location})
    end select
    return
  end function Equivalent_Circular_Orbit_Solver

  subroutine Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost,orbit,extremumType,radius,velocity,darkMatterHaloScale_)
    !!{
    Solves for the pericentric radius and velocity of {\normalfont \ttfamily orbit} in {\normalfont \ttfamily nodeHost}.
    !!}
    use :: Galactic_Structure_Options  , only : structureErrorCodeInfinite, structureErrorCodeSuccess    ,enumerationStructureErrorCodeType, radiusLarge
    use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleClass
    use :: Error                       , only : Error_Report
    use :: Galacticus_Nodes            , only : nodeComponentBasic        , treeNode
    use :: Kepler_Orbits               , only : keplerOrbit
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Root_Finder                 , only : rangeExpandMultiplicative , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive    , rootFinder
    implicit none
    type            (treeNode                         ), intent(inout), target :: nodeHost
    type            (keplerOrbit                      ), intent(inout)         :: orbit
    integer                                            , intent(in   )         :: extremumType
    double precision                                   , intent(  out)         :: radius                   , velocity
    class           (darkMatterHaloScaleClass         ), intent(inout)         :: darkMatterHaloScale_
    class           (nodeComponentBasic               ), pointer               :: basicHost
    double precision                                   , parameter             :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-6
    type            (rootFinder                       ), save                  :: finder
    logical                                            , save                  :: finderConstructed=.false.
    !$omp threadprivate(finder,finderConstructed)
    type            (keplerOrbit                      )                        :: orbitCurrent
    type            (enumerationStructureErrorCodeType)                        :: status
    double precision                                                           :: potential                , energyKinetic

    ! Get the mass distribution.
    massDistribution__ => nodeHost%massDistribution()
    ! Convert the orbit to the potential of the current halo in which the satellite finds itself.
    orbitCurrent=Satellite_Orbit_Convert_To_Current_Potential(orbit,nodeHost,darkMatterHaloScale_)
    ! Extract the orbital energy and angular momentum.
    orbitalEnergyInternal         =orbitCurrent%energy         ()
    orbitalAngularMomentumInternal=orbitCurrent%angularMomentum()
    ! Check if node or orbit differs from previous one for which we performed calculations.
    basicHost => nodeHost%basic()
    if     (                                                                  &
         &   nodeHost %uniqueID()           /= lastUniqueID                   &
         &  .or.                                                              &
         &   basicHost%time    ()           /= timePrevious                   &
         &  .or.                                                              &
         &   orbitalEnergyInternal          /= orbitalEnergyPrevious          &
         &  .or.                                                              &
         &   orbitalAngularMomentumInternal /= orbitalAngularMomentumPrevious &
         & ) call Satellite_Orbit_Reset(nodeHost)
    ! Determine if we need to compute the extremum properties.
    if     (                                                                      &
         &   (extremumType == extremumPericenter .and. .not.pericenterCalculated) &
         &  .or.                                                                  &
         &   (extremumType == extremumApocenter  .and. .not. apocenterCalculated) &
         & ) then
       ! Set a pointer to the host node.
       activeNode                    => nodeHost
       radiusVirial__                =  darkMatterHaloScale_%radiusVirial(nodeHost)
       massVirial__                  =  basicHost           %mass        (        )
       ! Record previous orbital properties.
       lastUniqueID                  =nodeHost %uniqueID()
       timePrevious                  =basicHost%time    ()
       orbitalEnergyPrevious         =orbitalEnergyInternal
       orbitalAngularMomentumPrevious=orbitalAngularMomentumInternal
       ! Catch orbits which are close to being circular.
       if      (   Extremum_Solver(orbitCurrent%radius()) >= 0.0d0             ) then
          ! Orbit is at extremum.
          radius=orbitCurrent%radius()
       else if (                                                                      &
            &      extremumType                           == extremumPericenter       &
            &   .and.                                                                 &
            &      orbitalAngularMomentumInternal         <= 0.0d0                    &
            &  ) then
          ! Orbit is radial, so pericenter is zero.
          radius=0.0d0
       else if (                                                                      &
            &      extremumType                           == extremumApocenter        &
            &   .and.                                                                 &
            &      Extremum_Solver(radiusLarge)           <  0.0d0                    &
            &  ) then
          ! Orbit is unbound, so apocenter is infinite.
          radius=huge(0.0d0)
       else
          if (.not.finderConstructed) then
             finder           =rootFinder(                                     &
                  &                       rootFunction     =Extremum_Solver  , &
                  &                       toleranceAbsolute=toleranceAbsolute, &
                  &                       toleranceRelative=toleranceRelative  &
                  &                      )
             finderConstructed=.true.
          end if
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
                  &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                  &                   rangeExpandType              =rangeExpandMultiplicative      &
                  &                  )
          end select
          radius=finder%find(rootGuess=orbitCurrent%radius())
       end if
       ! Get the orbital velocity at this radius.
       if (orbitalAngularMomentumInternal > 0.0d0) then
          ! Orbit is non-radial - use angular momentum to find velocity.
          velocity=orbitalAngularMomentumInternal/radius
       else
          ! Orbit is radial - use energy to find velocity.
          potential=Satellite_Orbit_Potential(radius,darkMatterHaloScale_%radiusVirial(activeNode),basicHost%mass(),status)
          select case (status%ID)
          case (structureErrorCodeSuccess %ID)
             energyKinetic=max(orbitalEnergyInternal-potential,0.0d0)
             velocity     =sqrt(2.0d0*energyKinetic)
          case (structureErrorCodeInfinite%ID)
             ! The gravitational potential is negative infinity at this radius (most likely zero
             ! radius). Velocity is formally infinite. Return speed of light as a suitably fast
             ! value.
             velocity=speedLight/kilo
          case default
             call Error_Report('dark matter potential evaluation failed'//{introspection:location})
          end select
       end if
       ! Store values and record that they are computed.
       select case (extremumType)
       case (extremumPericenter)
          pericenterRadius    =radius
          pericenterVelocity  =velocity
          pericenterCalculated=.true.
       case (extremumApocenter )
          apocenterRadius     =radius
          apocenterVelocity   =velocity
          apocenterCalculated =.true.
       end select
    else
       ! Use the previously computed values.
       select case (extremumType)
       case (extremumPericenter)
          radius  =pericenterRadius
          velocity=pericenterVelocity
       case (extremumApocenter )
          radius  = apocenterRadius
          velocity= apocenterVelocity
       end select
    end if
    !![
    <objectDestructor name="massDistribution__"/>
    !!]
    return
  end subroutine Satellite_Orbit_Extremum_Phase_Space_Coordinates

  double precision function Extremum_Solver(radius)
    !!{
    Root function used in finding orbital extremum radius.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: potential

    potential=Satellite_Orbit_Potential(radius,radiusVirial__,massvirial__)
    Extremum_Solver=potential+0.5d0*(orbitalAngularMomentumInternal/radius)**2-orbitalEnergyInternal
    return
  end function Extremum_Solver

  function Satellite_Orbit_Convert_To_Current_Potential(orbit,currentHost,darkMatterHaloScale_) result(orbitCurrent)
    !!{
    Takes a virial orbit and adjusts the energy to account for the change in the definition of potential between the original
    halo in which the orbit was defined and the current halo. Since the potential at the virial radius of halos is always
    defined to be $\Phi(r_\mathrm{vir}) = - V_\mathrm{vir}^2$ then the specific energy transforms as:
    \begin{equation}
    e \rightarrow e + V^2_\mathrm{vir,0} + \Phi(r_\mathrm{vir,0}),
    \end{equation}
    where subscript $0$ refers to the original halo in which the orbit was defined and $\Phi(r)$ is the potential of the
    current halo.
    !!}
    use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (keplerOrbit             )                :: orbitCurrent
    type            (keplerOrbit             ), intent(inout) :: orbit
    type            (treeNode                ), intent(inout) :: currentHost
    class           (darkMatterHaloScaleClass), intent(inout) :: darkMatterHaloScale_
    double precision                                          :: potentialHost         , radiusVirialOriginal, &
         &                                                       velocityVirialOriginal, radiusVirial

    ! Compute the properties of the initial orbit, and the current potential.
    radiusVirialOriginal   =gravitationalConstant_internal*orbit%massHost()/orbit%velocityScale()**2
    velocityVirialOriginal =                                                orbit%velocityScale()    
    radiusVirial           =darkMatterHaloScale_%radiusVirial    (currentHost)
    potentialHost          =Satellite_Orbit_Potential(radiusVirialOriginal,radiusVirial,orbit%massHost())
    ! Create a new orbit with an adjusted energy.
    orbitCurrent=orbit
    call orbitCurrent%energySet(orbit%energy()+velocityVirialOriginal**2+potentialHost)
    return
  end function Satellite_Orbit_Convert_To_Current_Potential

  double precision function Satellite_Orbit_Potential(radius,radiusVirial,massVirial,status) result(potential)
    !!{
    Evaluate the gravitational potential under the convention that the potential at the virial radius is always
    $\Phi(r_\mathrm{vir}) = - V_\mathrm{vir}^2$, as is assumed for {\normalfont \ttfamily keplerOrbit} objects.
    !!}
    use :: Galactic_Structure_Options      , only : enumerationStructureErrorCodeType
    use :: Coordinates                     , only : coordinateCartesian              , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision                                   , intent(in   )           :: radius     , radiusVirial     , &
         &                                                                          massVirial
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type            (coordinateCartesian              )                          :: coordinates, coordinatesVirial

    coordinates      =[radius      ,0.0d0,0.0d0]
    coordinatesVirial=[radiusVirial,0.0d0,0.0d0]
    potential        =+massDistribution__%potentialDifference(coordinates,coordinatesVirial,status) &
         &            -gravitationalConstant_internal                                               &
         &            *  massVirial                                                                 &
         &            /radiusVirial
    return
  end function Satellite_Orbit_Potential
  
  subroutine Satellite_Orbit_Reset(node)
    !!{
    Reset the satellite orbit calculations.
    !!}
    implicit none
    type(treeNode), intent(inout) :: node

    pericenterCalculated=.false.
    apocenterCalculated =.false.
    lastUniqueID        =node%uniqueID()
    return
  end subroutine Satellite_Orbit_Reset

end module Satellite_Orbits
