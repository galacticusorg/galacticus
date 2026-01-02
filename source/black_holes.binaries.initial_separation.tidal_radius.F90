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

  !+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !!{
  Implements a class for black hole binary initial separation based on tidal disruption of the satellite galaxy.
  !!}

  use :: Mass_Distributions, only : massDistributionClass
  use :: Root_Finder       , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder

  !![
  <blackHoleBinaryInitialSeparation name="blackHoleBinaryInitialSeparationTidalRadius">
   <description>
    A black hole binary initial separation class that assumes an initial separation that corresponds to the distance at which
    the satellite galaxy is tidally stripped to its half-mass radius, thus only leaving the central massive black
    hole. Specifically, the initial radius is given by:
    \begin{equation}
    {M_\mathrm{sat} \over 2 r_\mathrm{sat,1/2}^3 } = - {\mathrm{d} \over \mathrm{d} r} {M_\mathrm{host}(r_\mathrm{initial})
    \over r_\mathrm{initial}^2},
    \end{equation}
    where $M_\mathrm{sat}$ is the mass of the satellite galaxy, $r_\mathrm{sat,1/2}$ is its half mass radius,
    $M_\mathrm{host}(r)$ is the mass of the host galaxy within radius $r$ and $r_\mathrm{initial}$ is the initial radius.
   </description>
  </blackHoleBinaryInitialSeparation>
  !!]
  type, extends(blackHoleBinaryInitialSeparationClass) :: blackHoleBinaryInitialSeparationTidalRadius
     !!{
     A black hole binary initial separation class in which the radius is based on tidal disruption of the satellite galaxy.
     !!}
     private
     type (rootFinder) :: finder
  contains
    procedure :: separationInitial => tidalRadiusSeparationInitial
  end type blackHoleBinaryInitialSeparationTidalRadius

  interface blackHoleBinaryInitialSeparationTidalRadius
     !!{
     Constructors for the \refClass{blackHoleBinaryInitialSeparationTidalRadius} black hole binary initial radius class.
     !!}
     module procedure tidalRadiusConstructorParameters
     module procedure tidalRadiusConstructorInternal
  end interface blackHoleBinaryInitialSeparationTidalRadius

  ! Module-scope variables used in root finding.
  class           (massDistributionClass), pointer :: massDistribution_
  double precision                                 :: massHalf         , radiusMassHalf
  !$omp threadprivate(radiusMassHalf,massHalf,massDistribution_)

contains

  function tidalRadiusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryInitialSeparationTidalRadius} black hole binary recoil class which takes a parameter list as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (blackHoleBinaryInitialSeparationTidalRadius)                :: self
    type (inputParameters                            ), intent(inout) :: parameters

    self=blackHoleBinaryInitialSeparationTidalRadius()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalRadiusConstructorParameters

  function tidalRadiusConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{blackHoleBinaryInitialSeparationTidalRadius} black hole binary recoil class.
    !!}
    implicit none
    type (blackHoleBinaryInitialSeparationTidalRadius) :: self

    self%finder=rootFinder(                                                             &
            &              rootFunction=tidalRadiusRoot                               , &
            &              rangeExpandDownward          =0.5d+0                       , &
            &              rangeExpandUpward            =2.0d+0                       , &
            &              rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &              rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &              rangeExpandType              =rangeExpandMultiplicative    , &
            &              toleranceAbsolute            =0.0d+0                       , &
            &              toleranceRelative            =1.0d-6                         &
            &             )
     return
  end function tidalRadiusConstructorInternal

  double precision function tidalRadiusSeparationInitial(self,node,nodeHost)
    !!{
    Returns an initial separation for a binary black holes through tidal disruption.
    !!}
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Galacticus_Nodes          , only : nodeComponentBlackHole
    implicit none
    class(blackHoleBinaryInitialSeparationTidalRadius), intent(inout), target :: self
    type (treeNode                                   ), intent(inout), target :: nodeHost , node
    class(nodeComponentBlackHole                     ), pointer               :: blackHole
    !$GLC attributes unused :: self

    ! Assume zero separation by default.
    tidalRadiusSeparationInitial=0.0d0
    ! Get the black hole component.
    blackHole => node%blackHole(instance=1)
    ! If the primary black hole has zero mass (i.e. has been ejected), then return immediately.
    if (blackHole%mass() <= 0.0d0) return
    ! Get the mass distribution.
    massDistribution_ => node%massDistribution(massType=massTypeGalactic)
    ! Get the half-mass radius of the satellite galaxy.
    radiusMassHalf=massDistribution_%radiusEnclosingMass (massFractional=0.5d0         )
    ! Get the mass within the half-mass radius.
    massHalf      =massDistribution_%massEnclosedBySphere(               radiusMassHalf)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Solve for the radius around the host at which the satellite gets disrupted.
    if (radiusMassHalf > 0.0d0 .and. massHalf > 0.0d0) then
       massDistribution_            => nodeHost       %massDistribution(massType =massTypeGalactic                                           )
       tidalRadiusSeparationInitial =  self    %finder%find            (rootGuess=massDistribution_%radiusEnclosingMass(massFractional=0.5d0))
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end if
    return
  end function tidalRadiusSeparationInitial

  double precision function tidalRadiusRoot(radius)
    !!{
    Root function used in solving for the radius of tidal disruption of a satellite galaxy.
    !!}
    use :: Galactic_Structure_Options, only : massTypeGalactic
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    tidalRadiusRoot=+massDistribution_%massEnclosedBySphere(radius) &
         &          /massHalf                                       &
         &          -(                                              &
         &            +radius                                       &
         &            /radiusMassHalf                               &
         &           )**3
    return
  end function tidalRadiusRoot
