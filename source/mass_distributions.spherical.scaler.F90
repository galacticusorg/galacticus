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
  Implementation of a scaling spherical mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionSphericalScaler">
    <description>
      A mass distribution class for scaling spherical mass distributions. Specifically, the density at position $\mathbf{x}$ is
      given by
      \begin{equation}
      \rho(\mathbf{x}) = \frac{f_\mathrm{M}}{f_\mathrm{r}^3} \rho^\prime(\mathbf{x}/f_\mathrm{r}),
      \end{equation}      
      where $\rho^\prime(\mathbf{x})$ is the original mass distribution, and $f_\mathrm{r}=${\normalfont \ttfamily
      [factorScalingLength]}, and $f_\mathrm{M}=${\normalfont \ttfamily [factorScalingMass]}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionSphericalScaler
     !!{
     A mass distribution class for scaling spherical mass distributions.
     !!}
     class           (massDistributionSpherical), pointer :: massDistribution_   => null()
     double precision                                     :: factorScalingLength          , factorScalingMass
   contains
     final     ::                                      sphericalScalerDestructor
     procedure :: massTotal                         => sphericalScalerMassTotal
     procedure :: density                           => sphericalScalerDensity
     procedure :: densityGradientRadial             => sphericalScalerDensityGradientRadial
     procedure :: densityRadialMoment               => sphericalScalerDensityRadialMoment
     procedure :: massEnclosedBySphere              => sphericalScalerMassEnclosedBySphere
     procedure :: velocityRotationCurveMaximum      => sphericalScalerVelocityRotationCurveMaximum
     procedure :: radiusRotationCurveMaximum        => sphericalScalerRadiusRotationCurveMaximum
     procedure :: radiusEnclosingMass               => sphericalScalerRadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => sphericalScalerRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => sphericalScalerRadiusFromSpecificAngularMomentum
     procedure :: fourierTransform                  => sphericalScalerFourierTransform
     procedure :: radiusFreefall                    => sphericalScalerRadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => sphericalScalerRadiusFreefallIncreaseRate
     procedure :: energyPotential                   => sphericalScalerEnergyPotential
     procedure :: densitySphericalAverage           => sphericalScalerDensitySphericalAverage
     procedure :: rotationCurve                     => sphericalScalerRotationCurve
     procedure :: rotationCurveGradient             => sphericalScalerRotationCurveGradient
     procedure :: potential                         => sphericalScalerPotential
     procedure :: radiusHalfMass                    => sphericalScalerRadiusHalfMass
     procedure :: tidalTensor                       => sphericalScalerTidalTensor
     procedure :: acceleration                      => sphericalScalerAcceleration
     procedure :: positionSample                    => sphericalScalerPositionSample
  end type massDistributionSphericalScaler

  interface massDistributionSphericalScaler
     !!{
     Constructors for the {\normalfont \ttfamily sphericalScaler} mass distribution class.
     !!}
     module procedure sphericalScalerConstructorParameters
     module procedure sphericalScalerConstructorInternal
  end interface massDistributionSphericalScaler

contains

  function sphericalScalerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sphericalScaler} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionSphericalScaler)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (massDistributionClass          ), pointer       :: massDistribution_
    double precision                                                 :: factorScalingLength, factorScalingMass

    !![
    <inputParameter>
      <name>factorScalingLength</name>
      <description>The factor by which to scale lengths.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>factorScalingMass</name>
      <description>The factor by which to scale the mass.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalScaler(factorScalingLength,factorScalingMass,massDistribution_)
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalScalerConstructorParameters
  
  function sphericalScalerConstructorInternal(factorScalingLength,factorScalingMass,massDistribution_) result(self)
    !!{
    Constructor for ``sphericalScaler'' convergence class.
    !!}
    implicit none
    type            (massDistributionSphericalScaler)                        :: self
    class           (massDistributionSpherical      ), intent(in   ), target :: massDistribution_
    double precision                                 , intent(in   )         :: factorScalingLength, factorScalingMass
    !![
    <constructorAssign variables="factorScalingLength, factorScalingMass, *massDistribution_"/>
    !!]
 
    self%componentType=self%massDistribution_%componentType
    self%     massType=self%massDistribution_%     massType
    self%dimensionless=.false.
    return
  end function sphericalScalerConstructorInternal

  subroutine sphericalScalerDestructor(self)
    !!{
    Destructor for the ``sphericalScaler'' mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalScaler), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalScalerDestructor

  double precision function sphericalScalerMassTotal(self,componentType,massType)
    !!{
    Return the total mass in a scaled spherical distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout)           :: self
    type (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerMassTotal=+self%massDistribution_%massTotal        (               &
         &                                                             componentType, &
         &                                                             massType       &
         &                                                            )               &
         &                   *self                  %factorScalingMass
    return
  end function sphericalScalerMassTotal

  double precision function sphericalScalerDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout)              :: self
    class(coordinate                     ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType   ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional    :: massType
    class(coordinate                     )               , allocatable :: coordinatesScaled

    if (self%factorScalingMass > 0.0d0) then
       call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
       sphericalScalerDensity=+self%massDistribution_%density          (                           &
            &                                                           coordinatesScaled        , &
            &                                                           componentType            , &
            &                                                           massType                   &
            &                                                          )                           &
            &                 *self                  %factorScalingMass                            &
            &                 /self                  %factorScalingLength**3
    else
       sphericalScalerDensity=+0.0d0
    end if
    return
  end function sphericalScalerDensity

  double precision function sphericalScalerDensityGradientRadial(self,coordinates,logarithmic,componentType,massType)
    !!{
    Return the density gradient in the radial direction in a scaled spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalScaler), intent(inout), target      :: self
    class  (coordinate                     ), intent(in   )              :: coordinates
    logical                                 , intent(in   ), optional    :: logarithmic
    type   (enumerationComponentTypeType   ), intent(in   ), optional    :: componentType
    type   (enumerationMassTypeType        ), intent(in   ), optional    :: massType
    class  (coordinate                     )               , allocatable :: coordinatesScaled
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    sphericalScalerDensityGradientRadial=+self%massDistribution_%densityGradientRadial(                           &
         &                                                                             coordinatesScaled        , &
         &                                                                             logarithmic              , &
         &                                                                             componentType            , &
         &                                                                             massType                   &
         &                                                                            )
    if (.not.logarithmic)                                                             &
         & sphericalScalerDensityGradientRadial=+sphericalScalerDensityGradientRadial &
         &                                      *self%factorScalingMass               &
         &                                      /self%factorScalingLength**4
    return
  end function sphericalScalerDensityGradientRadial

  double precision function sphericalScalerMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: radius
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerMassEnclosedBySphere=+self%massDistribution_%massEnclosedBySphere(                           &
         &                                                                                 radius               &
         &                                                                           /self%factorScalingLength, &
         &                                                                           componentType            , &
         &                                                                           massType                   &
         &                                                                          )                           &
         &                              *self                  %factorScalingMass
    return
  end function sphericalScalerMassEnclosedBySphere

  double precision function sphericalScalerPotential(self,coordinates,componentType,massType,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class(massDistributionSphericalScaler  ), intent(inout), target      :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    type (enumerationStructureErrorCodeType), intent(  out), optional    :: status
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    sphericalScalerPotential=+self%massDistribution_%potential                     (                   &
         &                                                                          coordinatesScaled, &
         &                                                                          componentType    , &
         &                                                                          massType         , &
         &                                                                          status             &
         &                                                                         )                   &
         &                  *self                  %factorScalingMass                                  &
         &                  /self                  %factorScalingLength
    if (self%massDistribution_%isDimensionless())                    &
         & sphericalScalerPotential=+sphericalScalerPotential        &
         &                          *gravitationalConstantGalacticus
    return
  end function sphericalScalerPotential

  double precision function sphericalScalerDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Computes radial moments of the density in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                          , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerDensityRadialMoment=0.0d0
    !![
    <conditionalCall>
      <call>sphericalScalerDensityRadialMoment=self%massDistribution_%densityRadialMoment(moment=moment,isInfinite=isInfinite,componentType=componentType,massType=massType{conditions})</call>
      <argument name="radiusMinimum" value="radiusMinimum/self%factorScalingLength" condition="present(radiusMinimum)"/>
      <argument name="radiusMaximum" value="radiusMaximum/self%factorScalingLength" condition="present(radiusMaximum)"/>
    </conditionalCall>
    !!]
    sphericalScalerDensityRadialMoment=+     sphericalScalerDensityRadialMoment  &
         &                             *self%factorScalingMass                   &
         &                             /self%factorScalingLength**(3.0d0-moment)
    return    
  end function sphericalScalerDensityRadialMoment

  double precision function sphericalScalerRadiusHalfMass(self,componentType,massType)
    !!{
    Computes the half-mass radius in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout)           :: self
    type (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerRadiusHalfMass=+self%massDistribution_%radiusHalfMass     (componentType,massType) &
         &                        *self                  %factorScalingLength
    return
  end function sphericalScalerRadiusHalfMass

  function sphericalScalerAcceleration(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, gravitationalConstantGalacticus, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                 , dimension(3)               :: sphericalScalerAcceleration
    class           (massDistributionSphericalScaler), intent(inout)              :: self
    class           (coordinate                     ), intent(in   )              :: coordinates
    type            (enumerationComponentTypeType   ), intent(in   ), optional    :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional    :: massType
    class           (coordinate                     )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    sphericalScalerAcceleration=+self%massDistribution_%acceleration         (                     &
         &                                                                    coordinatesScaled  , &
         &                                                                          componentType, &
         &                                                                          massType       &
         &                                                                   )                     &
         &                      *self                  %factorScalingMass                          &
         &                      /self                  %factorScalingLength**2
    if (self%massDistribution_%isDimensionless())                       &
         & sphericalScalerAcceleration=+sphericalScalerAcceleration     &
         &                             *kilo                            &
         &                             *gigaYear                        &
         &                             /megaParsec                      &
         &                             *gravitationalConstantGalacticus
    return
  end function sphericalScalerAcceleration

  double precision function sphericalScalerDensitySphericalAverage(self,radius,componentType,massType)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: radius
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerDensitySphericalAverage=+self%massDistribution_%densitySphericalAverage(                           &
         &                                                                                 +     radius               &
         &                                                                                 /self%factorScalingLength, &
         &                                                                                 componentType            , &
         &                                                                                 massType                   &
         &                                                                                )                           &
         &                                 *self                  %factorScalingMass                                  &
         &                                 /self                  %factorScalingLength**3
    return
  end function sphericalScalerDensitySphericalAverage

  double precision function sphericalScalerRotationCurve(self,radius,componentType,massType)
    !!{
    Return the mid-plane rotation curve for a scaled spherical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: radius
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerRotationCurve=+      self%massDistribution_%rotationCurve      (                           &
         &                                                                         +     radius               &
         &                                                                         /self%factorScalingLength, &
         &                                                                         componentType            , &
         &                                                                         massType                   &
         &                                                                        )                           &
         &                       *sqrt(                                                                       &
         &                             +self                  %factorScalingMass                              &
         &                             /self                  %factorScalingLength                            &
         &                            )
    if (self%massDistribution_%isDimensionless())                              &
         & sphericalScalerRotationCurve=+sphericalScalerRotationCurve          &
         &                              *sqrt(gravitationalConstantGalacticus)
    return
  end function sphericalScalerRotationCurve

  double precision function sphericalScalerRotationCurveGradient(self,radius,componentType,massType)
    !!{
    Return the mid-plane rotation curve gradient (specifically, $\mathrm{d}V^2_\mathrm{c}/\mathrm{d}r$) for a scaled spherical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: radius
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerRotationCurveGradient=+self%massDistribution_%rotationCurveGradient(                           &
         &                                                                             +     radius               &
         &                                                                             /self%factorScalingLength, &
         &                                                                             componentType            , &
         &                                                                             massType                   &
         &                                                                            )                           &
         &                               *self%factorScalingMass                                                  &
         &                               /self%factorScalingLength**2
    if (self%massDistribution_%isDimensionless())                                     &
         & sphericalScalerRotationCurveGradient=+sphericalScalerRotationCurveGradient &
         &                                      *gravitationalConstantGalacticus
    return
  end function sphericalScalerRotationCurveGradient

  function sphericalScalerTidalTensor(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type (tensorRank2Dimension3Symmetric )                             :: sphericalScalerTidalTensor
    class(massDistributionSphericalScaler), intent(inout)              :: self
    class(coordinate                     ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType   ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional    :: massType
    class(coordinate                     )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    sphericalScalerTidalTensor=+self%massDistribution_%tidalTensor           (                   &
         &                                                                    coordinatesScaled, &
         &                                                                    componentType    , &
         &                                                                    massType           &
         &                                                                   )                   &
         &                     *self                  %factorScalingMass                         &
         &                     /self                  %factorScalingLength**3
    if (self%massDistribution_%isDimensionless())                      &
         & sphericalScalerTidalTensor=+sphericalScalerTidalTensor      &
         &                            *gravitationalConstantGalacticus
    return
  end function sphericalScalerTidalTensor
  
  function sphericalScalerPositionSample(self,randomNumberGenerator_,componentType,massType)
    !!{
    Sample a position from a scaled spherical mass distribution.
    !!}
    implicit none
    double precision                                 , dimension(3)            :: sphericalScalerPositionSample
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    class           (randomNumberGeneratorClass     ), intent(inout)           :: randomNumberGenerator_
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    sphericalScalerPositionSample=+self%massDistribution_%positionSample     (randomNumberGenerator_,componentType,massType) &
         &                        *self                  %factorScalingLength
    return
  end function sphericalScalerPositionSample

  double precision function sphericalScalerFourierTransform(self,radiusOuter,wavenumber,componentType,massType) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical, scaled mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: radiusOuter  , wavenumber
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    fourierTransform=self%massDistribution_%fourierTransform(radiusOuter/self%factorScalingLength,wavenumber*self%factorScalingLength,componentType,massType)
    return
  end function sphericalScalerFourierTransform

  double precision function sphericalScalerRadiusFreefall(self,time,componentType,massType) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: time
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
    
    radius=+self%massDistribution_%radiusFreefall(                                    &
         &                                        +time                               &
         &                                        *sqrt(                              &
         &                                              +self%factorScalingMass       &
         &                                              /self%factorScalingLength**3  &
         &                                             )                            , &
         &                                        componentType                     , &
         &                                        massType                            &
         &                                       )                                    &
         & *                                       self%factorScalingLength
    return
  end function sphericalScalerRadiusFreefall
  
  double precision function sphericalScalerRadiusFreefallIncreaseRate(self,time,componentType,massType) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: time
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    radiusIncreaseRate=+self%massDistribution_%radiusFreefallIncreaseRate(                                    &
         &                                                                +time                               &
         &                                                                *sqrt(                              &
         &                                                                      +self%factorScalingMass       &
         &                                                                      /self%factorScalingLength**3  &
         &                                                                     )                            , &
         &                                                                componentType                     , &
         &                                                                massType                            &
         &                                                               )                                    &
         &             *                                                   sqrt(                              &
         &                                                                      +self%factorScalingLength**5  &
         &                                                                      /self%factorScalingMass       &
         &                                                                     )
    return
  end function sphericalScalerRadiusFreefallIncreaseRate

  double precision function sphericalScalerEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout) :: self
    double precision                                 , intent(in   ) :: radiusOuter

    energy =+self%massDistribution_%energyPotential(radiusOuter/self%factorScalingLength   ) &
         &  *                                                   self%factorScalingMass  **2  &
         &  /                                                   self%factorScalingLength
    return
  end function sphericalScalerEnergyPotential

  double precision function sphericalScalerVelocityRotationCurveMaximum(self,componentType,massType) result(velocity)
    !!{
    Return the peak velocity in the rotation curve for an spherical scaled mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout)           :: self
    type (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional :: massType

    velocity=+self%massDistribution_%velocityRotationCurveMaximum(componentType,massType) &
         &   *sqrt(                                                                       &
         &         +self%factorScalingMass                                                &
         &         /self%factorScalingLength                                              &
         &        )
    return
  end function sphericalScalerVelocityRotationCurveMaximum

  double precision function sphericalScalerRadiusRotationCurveMaximum(self,componentType,massType) result(radius)
    !!{
    Return the peak velocity in the rotation curve for an spherical scaled mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout), target   :: self
    type (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType        ), intent(in   ), optional :: massType

    radius=+self%massDistribution_%radiusRotationCurveMaximum(componentType,massType) &
         & *self%factorScalingLength
    return
  end function sphericalScalerRadiusRotationCurveMaximum

  double precision function sphericalScalerRadiusEnclosingMass(self,mass,massFractional,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for spherical scaled mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target   :: self
    double precision                                 , intent(in   ), optional :: mass         , massFractional
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    if (present(massFractional)) then
       radius=+self%massDistribution_%radiusEnclosingMass(                                  &
            &                                             massFractional=massFractional   , &
            &                                             componentType=componentType     , &
            &                                             massType=massType                 &
            &                                            )
    else if (present(mass)) then
       radius=+self%massDistribution_%radiusEnclosingMass(                                  &
            &                                             mass=mass*self%factorScalingMass, &
            &                                             componentType=componentType     , &
            &                                             massType=massType                 &
            &                                            )
    else
       radius=+0.0d0
       call Error_Report('either "mass" or "massFractional" must be provided'//{introspection:location})
    end if
    radius=+radius                   &
         & *self%factorScalingLength
    return
  end function sphericalScalerRadiusEnclosingMass
  
  double precision function sphericalScalerRadiusEnclosingDensity(self,density,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mean density for spherical scaled mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: density
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType

    radius=+self%massDistribution_%radiusEnclosingDensity(                              &
         &                                                +density                      &
         &                                                *self%factorScalingLength**3  &
         &                                                /self%factorScalingMass     , &
         &                                                componentType               , &
         &                                                massType                      &
         &                                               )                              &
         & *                                               self%factorScalingLength
    return
  end function sphericalScalerRadiusEnclosingDensity

  double precision function sphericalScalerRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific,componentType,massType) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for sphericalScaler mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: angularMomentumSpecific
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
    
    radius=+self%massDistribution_%radiusFromSpecificAngularMomentum(                                 &
         &                                                           +angularMomentumSpecific         &
         &                                                           /sqrt(                           &
         &                                                                 +self%factorScalingMass    &
         &                                                                 *self%factorScalingLength  &
         &                                                           )                              , &
         &                                                           componentType                  , &
         &                                                           massType                         &
         &                                                          )                                 &
         & *                                                                self%factorScalingLength
    return
  end function sphericalScalerRadiusFromSpecificAngularMomentum
