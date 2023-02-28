!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of an abstract mass distribution class for spherically symmetric distributions.
  !!}

  !![
  <massDistribution name="massDistributionSpherical" abstract="yes">
   <description>An abstract mass distribution class for spherically symmetric distributions.</description>
  </massDistribution>
  !!]
  type, extends(massDistributionClass), abstract :: massDistributionSpherical
     !!{
     Implementation of an abstract mass distribution class for spherically symmetric distributions.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Returns the radius enclosing half of the mass of the mass distribution." method="radiusHalfMass"     />
       <method description="Returns the radius enclosing the given mass."                            method="radiusEnclosingMass"/>
     </methods>
     !!]
     procedure :: symmetry             => sphericalSymmetry
     procedure :: massEnclosedBySphere => sphericalMassEnclosedBySphere
     procedure :: radiusHalfMass       => sphericalRadiusHalfMass
     procedure :: acceleration         => sphericalAcceleration
     procedure :: tidalTensor          => sphericalTidalTensor
     procedure :: radiusEnclosingMass  => sphericalRadiusEnclosingMass
     procedure :: positionSample       => sphericalPositionSample
  end type massDistributionSpherical

  ! Module scope variables used in integration and root finding.
  class           (massDistributionSpherical), pointer :: self_
  double precision                                     :: massTarget
  !$omp threadprivate(self_,massTarget)

contains

  function sphericalSymmetry(self)
    !!{
    Returns symmetry label for mass distributions with spherical symmetry.
    !!}
    implicit none
    type (enumerationMassDistributionSymmetryType)                :: sphericalSymmetry
    class(massDistributionSpherical              ), intent(inout) :: self
    !$GLC attributes unused :: self

    sphericalSymmetry=massDistributionSymmetrySpherical
    return
  end function sphericalSymmetry

  double precision function sphericalMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions using numerical integration.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (massDistributionSpherical   ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (integrator                  )                          :: integrator_

    if (.not.self%matches(componentType,massType)) then
       sphericalMassEnclosedBySphere=0.0d0
       return
    end if
    self_                         =>  self
    integrator_                   =   integrator(sphericalMassEnclosedBySphereIntegrand,toleranceRelative=1.0d-6)
    sphericalMassEnclosedBySphere =  +4.0d0                              &
         &                           *Pi                                 &
         &                           *integrator_%integrate(0.0d0,radius)
    return
  end function sphericalMassEnclosedBySphere

  double precision function sphericalMassEnclosedBySphereIntegrand(radius)
    !!{
    Enclosed mass integrand for spherical mass distributions.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: position

    position                              =[radius,0.0d0,0.0d0]
    sphericalMassEnclosedBySphereIntegrand=+radius                 **2 &
         &                                 *self_%density(position)
    return
  end function sphericalMassEnclosedBySphereIntegrand

  double precision function sphericalRadiusEnclosingMass(self,mass,componentType,massType)
    !!{
    Computes the radius enclosing a given mass in a spherically symmetric mass distribution using numerical root finding.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder, &
         &                     GSL_Root_fSolver_Brent
    implicit none
    class           (massDistributionSpherical   ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: mass
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (rootFinder                  ), save                    :: finder
    logical                                       , save                    :: finderConstructed=.false.
    !$omp threadprivate(finder,finderConstructed)
    double precision                              , parameter               :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-6

    if (mass <= 0.0d0 .or. .not.self%matches(componentType,massType)) then
       sphericalRadiusEnclosingMass=0.0d0
       return
    end if
    if (.not.finderConstructed) then
       finder           =rootFinder(                                                             &
            &                       rootFunction                 =sphericalMassRoot            , &
            &                       toleranceAbsolute            =toleranceAbsolute            , &
            &                       toleranceRelative            =toleranceRelative            , &
            &                       solverType                   =GSL_Root_fSolver_Brent       , &
            &                       rangeExpandUpward            =2.0d0                        , &
            &                       rangeExpandDownward          =0.5d0                        , &
            &                       rangeExpandType              =rangeExpandMultiplicative    , &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
            &                      )
       finderConstructed=.true.
    end if
    self_              => self
    massTarget          =  mass
    sphericalRadiusEnclosingMass =  finder%find(rootGuess=1.0d0)
    return
  end function sphericalRadiusEnclosingMass

  double precision function sphericalRadiusHalfMass(self,componentType,massType)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    implicit none
    class(massDistributionSpherical   ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       sphericalRadiusHalfMass=0.0d0
       return
    end if
    sphericalRadiusHalfMass=self%radiusEnclosingMass(0.5d0*self%massTotal())
    return
  end function sphericalRadiusHalfMass

  double precision function sphericalMassRoot(radius)
    !!{
    Root function used in finding half mass radii of spherically symmetric mass distributions.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    sphericalMassRoot=+self_%massEnclosedBySphere(radius) &
         &            -      massTarget
    return
  end function sphericalMassRoot
  
  function sphericalAcceleration(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=), coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gigaYear     , megaParsec         , gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                              , dimension(3)            :: sphericalAcceleration
    class           (massDistributionSpherical   ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (coordinateSpherical         )                          :: coordinatesSpherical
    type            (coordinateCartesian         )                          :: coordinatesCartesian
    double precision                                                        :: radius
    double precision                              , dimension(3)            :: positionCartesian

    if (.not.self%matches(componentType,massType)) then
       sphericalAcceleration=0.0d0
       return
    end if
    ! Get position in spherical and Cartesian coordinate systems.
    coordinatesSpherical=coordinates
    coordinatesCartesian=coordinates
    ! Compute the density at this position.
    positionCartesian    = coordinatesCartesian
    radius               =+coordinatesSpherical%r                   (                 )
    sphericalAcceleration=-self                %massEnclosedBySphere(radius           )    &
         &                *                                          positionCartesian     &
         &                /                                          radius            **3
    if (.not.self%isDimensionless())                              &
         & sphericalAcceleration=+sphericalAcceleration           &
         &                       *kilo                            &
         &                       *gigaYear                        &
         &                       /megaParsec                      &
         &                       *gravitationalConstantGalacticus
    return
  end function sphericalAcceleration

  function sphericalTidalTensor(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=)                  , coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    use :: Tensors                         , only : tensorIdentityR2D3Sym          , tensorNullR2D3Sym   ,assignment(=)      , operator(*)
    use :: Vectors                         , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                          :: sphericalTidalTensor
    class           (massDistributionSpherical     ), intent(inout)           :: self
    class           (coordinate                    ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType  ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType       ), intent(in   ), optional :: massType
    type            (coordinateSpherical           )                          :: coordinatesSpherical
    type            (coordinateCartesian           )                          :: coordinatesCartesian
    double precision                                                          :: radius              , massEnclosed, &
         &                                                                       density
    double precision                                , dimension(3)            :: positionCartesian
    type            (tensorRank2Dimension3Symmetric)                          :: positionTensor

    if (.not.self%matches(componentType,massType)) then
       sphericalTidalTensor=tensorNullR2D3Sym
       return
    end if
    ! Get position in spherical and Cartesian coordinate systems.
    coordinatesSpherical=coordinates
    coordinatesCartesian=coordinates
    ! Compute the enclosed mass and density at this position.
    positionCartesian= coordinatesCartesian
    radius           =+coordinatesSpherical%r                   (           )
    massEnclosed     =+self                %massEnclosedBySphere(radius     ) 
    density          =+self                %density             (coordinates) 
    positionTensor   =Vector_Outer_Product(positionCartesian,symmetrize=.true.)
    ! Find the gravitational tidal tensor.
    sphericalTidalTensor=-(massEnclosed         /radius**3)*tensorIdentityR2D3Sym &
         &               +(massEnclosed*3.0d0   /radius**5)*positionTensor        &
         &               -(density     *4.0d0*Pi/radius**2)*positionTensor
    ! For dimensionful profiles, add the appropriate normalization.
    if (.not.self%isDimensionless())                              &
         & sphericalTidalTensor=+sphericalTidalTensor             &
         &                       *gravitationalConstantGalacticus
    return
  end function sphericalTidalTensor
  
  function sphericalPositionSample(self,randomNumberGenerator_,componentType,massType)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                              , dimension(3)            :: sphericalPositionSample
    class           (massDistributionSpherical   ), intent(inout)           :: self
    class           (randomNumberGeneratorClass  ), intent(inout)           :: randomNumberGenerator_
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: mass                   , radius, &
         &                                                                     theta                  , phi

    if (.not.self%matches(componentType,massType)) then
       sphericalPositionSample=0.0d0
       return
    end if
    ! Choose an enclosed mass and find the radius enclosing that mass. Choose angular
    ! coordinates at random and finally convert to Cartesian.
    mass  =+              self                  %massTotal          (    )        &
         & *              randomNumberGenerator_%uniformSample      (    )
    radius=+              self                  %radiusEnclosingMass(mass)
    phi   =+     2.0d0*Pi*randomNumberGenerator_%uniformSample      (    )
    theta =+acos(2.0d0   *randomNumberGenerator_%uniformSample      (    )-1.0d0)
    sphericalPositionSample=+radius                &
         &                  *[                     &
         &                    sin(theta)*cos(phi), &
         &                    sin(theta)*sin(phi), &
         &                    cos(theta)           &       
         &                   ]
    return
  end function sphericalPositionSample
