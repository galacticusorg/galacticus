!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
       <method description="Returns the radius enclosing half of the mass of the mass distribution." method="radiusHalfMass"                     />
       <method description="Compute the potential energy of mass distribution."                      method="energyPotential"                    />
       <method description="Compute the kinetic energy of the mass distribution."                    method="energyKinetic"                      />
       <method description="Compute (numerically) the potential energy of mass distribution."        method="energyPotentialNumerical"           />
       <method description="Compute (numerically) the kinetic energy of the mass distribution."      method="energyKineticNumerical"             />
       <method description="Compute (numerically) a radial densitymoment of the mass distribution."  method="densityRadialMomentNumerical"       />
       <method description="Compute (numerically) the radial density gradient."                      method="densityGradientRadialNumerical"     />
       <method description="Compute (numerically) the mass enclosed by a sphere."                    method="massEnclosedBySphereNumerical"      />
       <method description="Compute (numerically) the rotation curve."                               method="rotationCurveNumerical"             />
       <method description="Compute (numerically) the gravitational potential."                      method="potentialNumerical"                 />
       <method description="Compute (numerically) the Fourier transform of the density profile."     method="fourierTransformNumerical"          />
       <method description="Compute (numerically) the freefall radius."                              method="radiusFreefallNumerical"            />
       <method description="Compute (numerically) the growth rate of the freefall radius."           method="radiusFreefallIncreaseRateNumerical"/>
       <method description="Compute (numerically) the energy."                                       method="energyNumerical"                    />
     </methods>
     !!]
     procedure :: symmetry                            => sphericalSymmetry
     procedure :: isSphericallySymmetric              => sphericalIsSphericallySymmetric
     procedure :: densityGradientRadial               => sphericalDensityGradientRadial
     procedure :: densityGradientRadialNumerical      => sphericalDensityGradientRadialNumerical
     procedure :: massEnclosedBySphere                => sphericalMassEnclosedBySphere
     procedure :: massEnclosedBySphereNumerical       => sphericalMassEnclosedBySphereNumerical
     procedure :: densityRadialMoment                 => sphericalDensityRadialMoment
     procedure :: densityRadialMomentNumerical        => sphericalDensityRadialMomentNumerical
     procedure :: potential                           => sphericalPotential
     procedure :: potentialNumerical                  => sphericalPotentialNumerical
     procedure :: fourierTransform                    => sphericalFourierTransform
     procedure :: fourierTransformNumerical           => sphericalFourierTransformNumerical
     procedure :: radiusFreefall                      => sphericalRadiusFreefall
     procedure :: radiusFreefallNumerical             => sphericalRadiusFreefallNumerical
     procedure :: radiusFreefallIncreaseRate          => sphericalRadiusFreefallIncreaseRate
     procedure :: radiusFreefallIncreaseRateNumerical => sphericalRadiusFreefallIncreaseRateNumerical
     procedure :: energy                              => sphericalEnergy
     procedure :: energyNumerical                     => sphericalEnergyNumerical
     procedure :: energyPotential                     => sphericalEnergyPotential
     procedure :: energyKinetic                       => sphericalEnergyKinetic
     procedure :: energyPotentialNumerical            => sphericalEnergyPotentialNumerical
     procedure :: energyKineticNumerical              => sphericalEnergyKineticNumerical
     procedure :: densitySphericalAverage             => sphericalDensitySphericalAverage
     procedure :: potentialDifferenceNumerical        => sphericalPotentialDifferenceNumerical
     procedure :: surfaceDensity                      => sphericalSurfaceDensity
     procedure :: chandrasekharIntegral               => sphericalChandrasekharIntegral
     procedure :: acceleration                        => sphericalAcceleration
     procedure :: tidalTensor                         => sphericalTidalTensor
     procedure :: positionSample                      => sphericalPositionSample
     procedure :: rotationCurve                       => sphericalRotationCurve
     procedure :: rotationCurveNumerical              => sphericalRotationCurve
     procedure :: rotationCurveGradient               => sphericalRotationCurveGradient
     procedure :: radiusHalfMass                      => sphericalRadiusHalfMass
  end type massDistributionSpherical

  ! Module scope variables used in integration and root finding.
  class           (massDistributionSpherical), pointer :: self_
  double precision                                     :: time_, radiusFreefall_
  !$omp threadprivate(self_,time_,radiusFreefall_)

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

  logical function sphericalIsSphericallySymmetric(self) result(isSphericallySymmetric)
    !!{
    Return true if the distribution is spherically symmetric.
    !!}
    implicit none
    class(massDistributionSpherical), intent(inout) :: self

    isSphericallySymmetric=.true.
    return
  end function sphericalIsSphericallySymmetric

  double precision function sphericalDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSpherical), intent(inout), target   :: self
    class  (coordinate               ), intent(in   )           :: coordinates
    logical                           , intent(in   ), optional :: logarithmic
    
    densityGradient=self%densityGradientRadialNUmerical(coordinates,logarithmic)
    return
  end function sphericalDensityGradientRadial

  double precision function sphericalDensityGradientRadialNumerical(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution using a numerical calculation.
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (massDistributionSpherical   ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    double precision                              , parameter               :: radiusLogarithmicStep=0.1d0
    type            (differentiator              )                          :: differentiator_
    double precision                                                        :: radius
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    self_           =>  self
    radius          =   coordinates    %rSpherical(                                     )
    differentiator_ =   differentiator            (densityEvaluate                      )
    densityGradient =  +differentiator_%derivative(log(radius)    ,radiusLogarithmicStep)
    if (.not.logarithmic_)                                    &
         & densityGradient=+     densityGradient              &
         &                 *self%density        (coordinates) &
         &                 /     radius
    return
  end function sphericalDensityGradientRadialNumerical

  double precision function densityEvaluate(radiusLogarithmic) result(density)
    !!{
    GSL-callable function to evaluate the density of the dark matter profile.
    !!}
      use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ), value :: radiusLogarithmic
    type            (coordinateSpherical)                       :: coordinates

    coordinates=[exp(radiusLogarithmic),0.0d0,0.0d0]
    density    =log(self_%density(coordinates))
    return
  end function densityEvaluate

  double precision function sphericalMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: radius

    mass=self%massEnclosedBySphereNumerical(radius)
    return
  end function sphericalMassEnclosedBySphere

  double precision function sphericalMassEnclosedBySphereNumerical(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions using numerical integration.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: radius
    type            (integrator               )                        :: integrator_

    self_       =>  self
    integrator_ =   integrator(sphericalMassEnclosedBySphereIntegrand,toleranceRelative=1.0d-6)
    mass        =  +4.0d0                              &
         &         *Pi                                 &
         &         *integrator_%integrate(0.0d0,radius)
    return
  end function sphericalMassEnclosedBySphereNumerical

  double precision function sphericalDensitySphericalAverage(self,radius)
    !!{
    Computes the density averaged over a spherical shell.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (coordinateSpherical      )                :: position

    ! For a spherical mass distribution, the density averaged over a spherical shell, is just the regular density at that radius.
    position                        =[radius,0.0d0,0.0d0]
    sphericalDensitySphericalAverage=self%density(position)
    return
  end function sphericalDensitySphericalAverage

  double precision function sphericalMassEnclosedBySphereIntegrand(radius) result(integrand)
    !!{
    Enclosed mass integrand for spherical mass distributions.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: position

    if (radius > 0.0d0) then
       position =[radius,0.0d0,0.0d0]
       integrand=+radius**2               &
            &    *self_%density(position)
    else
      integrand=+0.0d0 
    end if
    return
  end function sphericalMassEnclosedBySphereIntegrand

  double precision function sphericalRadiusHalfMass(self)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    implicit none
    class(massDistributionSpherical), intent(inout) :: self

    sphericalRadiusHalfMass=self%radiusEnclosingMass(0.5d0*self%massTotal())
    return
  end function sphericalRadiusHalfMass

  double precision function sphericalPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSpherical        ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status

    potential=self%potentialNumerical(coordinates,status)
    return
  end function sphericalPotential

  double precision function sphericalPotentialNumerical(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration           , only : integrator
    use :: Numerical_Comparison            , only : Values_Agree
    implicit none
    class           (massDistributionSpherical        ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                   , parameter               :: toleranceRelative  =1.0d-3
    double precision                                   , parameter               :: radiusMaximumFactor=1.0d+1
    type            (integrator                       ), save                    :: integrator_
    logical                                            , save                    :: initialized        =.false.
    !$omp threadprivate(integrator_,initialized)
    double precision                                                             :: radiusMaximum              , potentialPrevious, &
         &                                                                          radius

    if (present(status)) status=structureErrorCodeSuccess
    if (.not.initialized) then
       integrator_=integrator(integrandPotential,toleranceRelative=toleranceRelative)
       initialized=.true.
    end if
    self_             =>  self
    potential         =  +0.0d0
    potentialPrevious =  +1.0d0
    radius            =  +coordinates%rSpherical()
    radiusMaximum     =  +            radius
    do while (.not.Values_Agree(potential,potentialPrevious,relTol=toleranceRelative))
       potentialPrevious=+potential
       radiusMaximum    =+radiusMaximum                        &
            &            *radiusMaximumFactor
       potential        =+integrator_%integrate(               &
            &                                   radius       , &
            &                                   radiusMaximum  &
            &                                  )
    end do
    ! Convert to dimensionful units.    
    if (.not.self%isDimensionless()) potential=+gravitationalConstantGalacticus &
         &                                     *potential
    return
  end function sphericalPotentialNumerical

  double precision function sphericalPotentialDifferenceNumerical(self,coordinates1,coordinates2,status) result(potential)
    !!{
    Return the potential difference between the two specified
    {\normalfont \ttfamily coordinates} in a spherical mass
    distribution using a numerical calculation.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration           , only : integrator
    use :: Numerical_Comparison            , only : Values_Agree
    implicit none
    class           (massDistributionSpherical        ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates1               , coordinates2
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                   , parameter               :: toleranceRelative  =1.0d-3
    double precision                                   , parameter               :: radiusMaximumFactor=1.0d+1
    type            (integrator                       ), save                    :: integrator_
    logical                                            , save                    :: initialized        =.false.
    !$omp threadprivate(integrator_,initialized)

    if (present(status)) status=structureErrorCodeSuccess
    if (.not.initialized) then
       integrator_=integrator(integrandPotential,toleranceRelative=toleranceRelative)
       initialized=.true.
    end if
    self_     => self
    potential =  integrator_%integrate(                           &
         &                             coordinates1%rSpherical(), &
         &                             coordinates2%rSpherical()  &
         &                            )
    ! Convert to dimensionful units.    
    if (.not.self%isDimensionless()) potential=+gravitationalConstantGalacticus &
         &                                     *potential
    return
  end function sphericalPotentialDifferenceNumerical
  
  double precision function integrandPotential(radius)
    !!{
    Integrand for gravitational potential in a spherical mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
      
    if (radius > 0.0d0) then
       integrandPotential=-self_%massEnclosedBySphere(radius)    &
            &             /                           radius **2
    else
       integrandPotential=0.0d0
    end if
    return
  end function integrandPotential

  function sphericalAcceleration(self,coordinates) result(acceleration)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=), coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gigaYear     , megaParsec         , gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                              , dimension(3)  :: acceleration
    class           (massDistributionSpherical   ), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    type            (coordinateSpherical         )                :: coordinatesSpherical
    type            (coordinateCartesian         )                :: coordinatesCartesian
    double precision                                              :: radius
    double precision                              , dimension(3)  :: positionCartesian

    ! Get position in spherical and Cartesian coordinate systems.
    coordinatesSpherical=coordinates
    coordinatesCartesian=coordinates
    ! Compute the density at this position.
    positionCartesian=coordinatesCartesian
    radius           =coordinatesSpherical%r()
    if (radius > 0.0d0) then
       acceleration=-self%massEnclosedBySphere(radius           )    &
            &       *                          positionCartesian     &
            &       /                          radius            **3
       if (.not.self%isDimensionless())                     &
            & acceleration=+acceleration                    &
            &              *kilo                            &
            &              *gigaYear                        &
            &              /megaParsec                      &
            &              *gravitationalConstantGalacticus
    else
       acceleration=0.0d0
    end if
    return
  end function sphericalAcceleration

  function sphericalTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=)                  , coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    use :: Tensors                         , only : tensorIdentityR2D3Sym          , assignment(=)      , operator(*)
    use :: Vectors                         , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: sphericalTidalTensor
    class           (massDistributionSpherical     ), intent(inout) :: self
    class           (coordinate                    ), intent(in   ) :: coordinates
    type            (coordinateSpherical           )                :: coordinatesSpherical
    type            (coordinateCartesian           )                :: coordinatesCartesian
    double precision                                                :: radius              , massEnclosed, &
         &                                                             density
    double precision                                , dimension(3)  :: positionCartesian
    type            (tensorRank2Dimension3Symmetric)                :: positionTensor

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

  double precision function sphericalRotationCurve(self,radius)
    !!{
    Return the rotation curve for a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius

    if (radius <= 0.0d0) then
       sphericalRotationCurve=+0.0d0
    else
       sphericalRotationCurve=+sqrt(                                   &
            &                       +self%massEnclosedBySphere(radius) &
            &                       /                          radius  &
            &                      )
       ! Make dimensionful if necessary.
       if (.not.self%dimensionless) sphericalRotationCurve= &
            & +sqrt(gravitationalConstantGalacticus)        &
            & *sphericalRotationCurve
    end if
    return
  end function sphericalRotationCurve

  double precision function sphericalRotationCurveGradient(self,radius)
    !!{
    Return the rotation curve gradient for a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)                  , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (coordinateSpherical      )                :: position

    position                      =[radius,0.0d0,0.0d0]
    sphericalRotationCurveGradient=+4.0d0                                  &
         &                         *Pi                                     &
         &                         *                          radius       &
         &                         *self%density             (position)    &
         &                         -self%massEnclosedBySphere(radius  )    &
         &                         /                          radius   **2
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) sphericalRotationCurveGradient= &
         &  +gravitationalConstantGalacticus                     &
         &  *sphericalRotationCurveGradient
    return
  end function sphericalRotationCurveGradient

  function sphericalPositionSample(self,randomNumberGenerator_)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                              , dimension(3)  :: sphericalPositionSample
    class           (massDistributionSpherical   ), intent(inout) :: self
    class           (randomNumberGeneratorClass  ), intent(inout) :: randomNumberGenerator_
    double precision                                              :: mass                   , radius, &
         &                                                           theta                  , phi

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

  double precision function sphericalSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(massDistributionSpherical), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates

    sphericalSurfaceDensity=0.0d0
    call Error_Report('surface density is not defined for spherically-symmetric distributions'//{introspection:location})
    return
  end function sphericalSurfaceDensity

  function sphericalChandrasekharIntegral(self,massDistributionEmbedding,massDistributionPerturber,massPerturber,coordinates,velocity)
    !!{
    Compute the Chandrasekhar integral at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    use :: Coordinates               , only : coordinateCartesian  , assignment(=)
    use :: Numerical_Constants_Math  , only : Pi
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeAll
    use :: Ideal_Gases_Thermodynamics, only : Ideal_Gas_Sound_Speed
    use :: Error                     , only : Error_Report
    implicit none
    double precision                              , dimension(3)  :: sphericalChandrasekharIntegral
    class           (massDistributionSpherical   ), intent(inout) :: self
    class           (massDistributionClass       ), intent(inout) :: massDistributionEmbedding            , massDistributionPerturber
    double precision                              , intent(in   ) :: massPerturber
    class           (coordinate                  ), intent(in   ) :: coordinates                          , velocity
    double precision                              , dimension(3)  :: velocityCartesian_
    double precision                              , parameter     :: XvMaximum                     =10.0d0
    type            (coordinateCartesian         )                :: velocityCartesian
    double precision                                              :: radius                               , velocity_                , &
         &                                                           density                              , velocityDispersion       , &
         &                                                           xV
    !$GLC attributes unused :: massDistributionPerturber, massPerturber
    
    sphericalChandrasekharIntegral=0.0d0
    velocity_=velocity%rSpherical()
    if (velocity_ <= 0.0d0) return
    radius =coordinates%rSpherical(           )
    density=self       %density   (coordinates)
    if (density  <= 0.0d0) return
    if (.not.associated(self%kinematicsDistribution_)) call Error_Report('a kinematics distribution is needed to compute the Chandrasekhar integral'//{introspection:location})
    if (self%kinematicsDistribution_%isCollisional()) then
       velocityDispersion=Ideal_Gas_Sound_Speed(self%kinematicsDistribution_%temperature         (coordinates                          ))
    else
       velocityDispersion=                      self%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionEmbedding)
    end if
    if (velocityDispersion > 0.0d0) then    
       xV                         =+velocity_             &
            &                      /velocityDispersion    &
            &                      /sqrt(2.0d0)
    else
       
       xV                         =+huge(0.0d0)
    end if
    velocityCartesian             = velocity
    velocityCartesian_            = velocityCartesian
    sphericalChandrasekharIntegral=-density               &
         &                         *velocityCartesian_    &
         &                         /velocity_         **3
    if (Xv <= XvMaximum)                                                  &
         & sphericalChandrasekharIntegral=+sphericalChandrasekharIntegral &
         &                                *(                              &
         &                                  +erf ( xV   )                 &
         &                                  -2.0d0                        &
         &                                  *      xV                     &
         &                                  *exp (-xV**2)                 &
         &                                  /sqrt( Pi   )                 &
         &                                 )
    return
  end function sphericalChandrasekharIntegral

  double precision function sphericalFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter, wavenumber
 
    fourierTransform=self%fourierTransformNumerical(radiusOuter,wavenumber)
    return
  end function sphericalFourierTransform

  double precision function sphericalFourierTransformNumerical(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{   
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical mass
    distribution using a numerical calculation.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter, wavenumber
    type            (integrator               )                :: integrator_

    fourierTransform=0.0d0
    integrator_     = integrator                      (integrandFourierTransform,toleranceRelative=1.0d-3)
    fourierTransform=+integrator_%integrate           (0.0d0                    ,radiusOuter             ) &
         &           /self       %massEnclosedBySphere(                          radiusOuter             )
    return

  contains

    double precision function integrandFourierTransform(radius)
      !!{
      Integrand for Fourier transform of a spherical mass distribution.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      use :: Coordinates             , only : assignment(=), coordinateSpherical
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates
      
      if (radius > 0.0d0) then
         coordinates              =[radius,0.0d0,0.0d0]
         integrandFourierTransform=+4.0d0                     &
              &                    *Pi                        &
              &                    *               radius **2 &
              &                    *sin(wavenumber*radius)    &
              &                    /   (wavenumber*radius)    &
              &                    *self%density(coordinates)
      else
         integrandFourierTransform=0.0d0
      end if
      return
    end function integrandFourierTransform
    
  end function sphericalFourierTransformNumerical
  
  double precision function sphericalRadiusFreefall(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time
 
    radius=self%radiusFreefallNumerical(time)
    return
  end function sphericalRadiusFreefall
  
  double precision function sphericalRadiusFreefallNumerical(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily
    time} in a spherical mass distribution using a numerical
    calculation.
    !!}
    use :: Root_Finder                     , only : rangeExpandMultiplicative      , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus, Mpc_per_km_per_s_To_Gyr
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: time
    double precision                           , parameter             :: toleranceAbsolute  =0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder               )                        :: finder
    type            (coordinateSpherical      )                        :: coordinates
    double precision                                                   :: timeFreefallMinimum

    radius=0.0d0
    coordinates=[0.0d0,0.0d0,0.0d0]
    if (self%densityGradientRadial(coordinates,logarithmic=.true.) == 0.0d0) then
       ! For mass distributions with a constant density core, the potential in the center is harmonic. This means there is a
       ! minimum to the freefall time as a function of radius. Compute that minimum here so that we can return a zero radius for
       ! times less than this.
       timeFreefallMinimum=+sqrt(                                 &
            &                    + 3.0d0                          &
            &                    /16.0d0                          &
            &                    *Pi                              &
            &                    /gravitationalConstantGalacticus &
            &                    /self%density(coordinates)       &
            &                   )                                 &
            &              *Mpc_per_km_per_s_To_Gyr
    else
       timeFreefallMinimum=+0.0d0
    end if   
    if (time < timeFreefallMinimum) return
    self_  => self
    time_  =  time
    finder =  rootFinder(                                                             &
         &               rootFunction                 =rootRadiusFreefall           , &
         &               toleranceAbsolute            =toleranceAbsolute            , &
         &               toleranceRelative            =toleranceRelative            , &
         &               rangeExpandDownward          =0.5d0                        , &
         &               rangeExpandUpward            =2.0d0                        , &
         &               rangeExpandType              =rangeExpandMultiplicative    , &
         &               rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &               rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &              )
    radius=finder%find(rootGuess=1.0d0)
    return
  end function sphericalRadiusFreefallNumerical
  
  double precision function rootRadiusFreefall(radiusFreefall)
    !!{
    Root function used in finding the radius corresponding to a given freefall time.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radiusFreefall
    type            (integrator)                :: integrator_
    
    radiusFreefall_   =+radiusFreefall
    integrator_       = integrator              (integrandTimeFreefall,toleranceRelative=1.0d-3)
    rootRadiusFreefall=+integrator_   %integrate(0.0d0                ,radiusFreefall          ) &
         &             -time_
    return
  end function rootRadiusFreefall

  double precision function integrandTimeFreefall(radius)
    !!{
    Integrand for freefall time in a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)          , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    double precision                     , intent(in   ) :: radius
    double precision                                     :: potentialDifference
    type            (coordinateSpherical)                :: coordinates        , coordinatesFreefall

    coordinates        =[radius         ,0.0d0,0.0d0]
    coordinatesFreefall=[radiusFreefall_,0.0d0,0.0d0]
    potentialDifference=+self_%potentialDifference(coordinates,coordinatesFreefall)
    if (potentialDifference < 0.0d0) then
       integrandTimeFreefall=+Mpc_per_km_per_s_To_Gyr   &
            &                /sqrt(                     &
            &                      -2.0d0               &
            &                      *potentialDifference &
            &                     )
    else
       ! Avoid floating point errors arising from rounding errors.
       integrandTimeFreefall=0.0d0
    end if
    return
  end function integrandTimeFreefall

  double precision function sphericalRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time

    radiusIncreaseRate=self%radiusFreefallIncreaseRateNumerical(time)
    return
  end function sphericalRadiusFreefallIncreaseRate

  double precision function sphericalRadiusFreefallIncreaseRateNumerical(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an spherical mass
    distribution using a numerical calculation.
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time
    double precision                           , parameter     :: timeLogarithmicStep=1.0d-2
    type            (differentiator           )                :: differentiator_

    differentiator_   = differentiator            (radiusFreefallEvaluate                    )
    radiusIncreaseRate=+differentiator_%derivative(log(time)             ,timeLogarithmicStep) &
         &             /                               time
    return
  end function sphericalRadiusFreefallIncreaseRateNumerical

  double precision function radiusFreefallEvaluate(timeLogarithmic)
    !!{
    GSL-callable function to evaluate the freefall radius of the mass distribution.
    !!}
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    radiusFreefallEvaluate=self_%radiusFreefall(exp(timeLogarithmic))
    return
  end function radiusFreefallEvaluate

  double precision function sphericalEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: radiusOuter
    class           (massDistributionClass    ), intent(inout), target :: massDistributionEmbedding

    energy=self%energyNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergy

  double precision function sphericalEnergyNumerical(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution using a numerical calculation.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding

    energy=+self%energyPotential(radiusOuter                          ) &
         & +self%energyKinetic  (radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergyNumerical

  double precision function sphericalEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter

    energy=self%energyPotentialNumerical(radiusOuter)
    return
  end function sphericalEnergyPotential

  double precision function sphericalEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding

    energy=self%energyKineticNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergyKinetic
  
  double precision function sphericalEnergyPotentialNumerical(self,radiusOuter) result(energy)
    !!{
    Compute (numerically) the potential energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration           , only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    type            (integrator               )                :: integrator_

    integrator_= integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
    energy     =-0.5d0                                                    &
         &      *gravitationalConstantGalacticus                          &
         &      *(                                                        &
         &        +integrator_%integrate           (0.0d0,radiusOuter)    &
         &        +self       %massEnclosedBySphere(      radiusOuter)**2 &
         &        /            radiusOuter                                &
         &       )
    return

  contains

    double precision function integrandEnergyPotential(radius)
      !!{
      Integrand for potential energy of a spherical mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandEnergyPotential=(                                   &
              &                    +self%massEnclosedBySphere(radius) &
              &                    /                          radius  &
              &                   )**2
      else
         integrandEnergyPotential=0.0d0
      end if
      return
    end function integrandEnergyPotential

  end function sphericalEnergyPotentialNumerical

  double precision function sphericalEnergyKineticNumerical(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute (numerically) the kinetic energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding
    type            (integrator               )                :: integrator_

    integrator_= integrator(integrandEnergyKinetic,toleranceRelative=1.0d-3)
    energy     =+6.0d0                                    &
         &      *Pi                                       &
         &      *integrator_%integrate(0.0d0,radiusOuter)
    return

  contains

    double precision function integrandEnergyKinetic(radius)
      !!{
      Integrand for kinetic energy of the halo.
      !!}
      use :: Coordinates, only : coordinateSpherical, assignment(=)
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates

      if (radius > 0.0d0) then
         coordinates           =[radius,0.0d0,0.0d0]
         integrandEnergyKinetic=+self                        %density             (coordinates                          )    &
              &                 *self%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionEmbedding)**2 &
              &                 *                             radius                                                     **2
      else
         integrandEnergyKinetic=0.0d0
      end if
      return
    end function integrandEnergyKinetic

  end function sphericalEnergyKineticNumerical

  double precision function sphericalDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite

    densityRadialMoment=self%densityRadialMomentNumerical(moment,radiusMinimum,radiusMaximum,isInfinite)
    return
  end function sphericalDensityRadialMoment

  double precision function sphericalDensityRadialMomentNumerical(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for a spherical mass distribution using a numerical calculation.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite
    type            (integrator               )                          :: integrator_
    double precision                                                     :: radiusMinimum_
    
    densityRadialMoment=0.0d0
    if (.not.present(radiusMaximum)) call Error_Report('a maximum radius must be provided'//{introspection:location})
    if (present(isInfinite)) isInfinite=.false.
    radiusMinimum_=0.0d0
    if (present(radiusMinimum)) radiusMinimum_=radiusMinimum
    integrator_        = integrator           (integrandMoment,toleranceRelative=1.0d-3)
    densityRadialMoment=+integrator_%integrate(radiusMinimum_ ,radiusMaximum           )
    return

  contains

    double precision function integrandMoment(radius)
      !!{
      Integrand for radial density moment in a spherical mass distribution.
      !!}
      use :: Coordinates, only : assignment(=), coordinateSpherical
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates
      
      if (radius > 0.0d0) then
         coordinates     = [radius,0.0d0,0.0d0]
         integrandMoment=+self%density(coordinates)         &
              &          *             radius      **moment
      else
         integrandMoment=+0.0d0
      end if
      return
    end function integrandMoment
    
  end function sphericalDensityRadialMomentNumerical
