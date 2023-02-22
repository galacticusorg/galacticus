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
  Implementation of an scaling cylindrical mass distribution class.
  !!}
  
  !![
  <massDistribution name="massDistributionCylindricalScaler">
    <description>
      A mass distribution class for scaling cylindrical mass distributions. Specifically, the density at position $\mathbf{x}$ is
      given by
      \begin{equation}
      \rho(\mathbf{x}) = \frac{f_\mathrm{M}}{f_\mathrm{r}^3} \rho^\prime(\mathbf{x}/f_\mathrm{r}),
      \end{equation}      
      where $\rho^\prime(\mathbf{x})$ is the original mass distribution, and $f_\mathrm{r}=${\normalfont \ttfamily
      [factorScalingLength]}, and $f_\mathrm{M}=${\normalfont \ttfamily [factorScalingMass]}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionCylindrical) :: massDistributionCylindricalScaler
     !!{
     A mass distribution class for scaling cylindrical mass distributions.
     !!}
     private
     class           (massDistributionCylindrical), pointer :: massDistribution_   => null()
     double precision                                       :: factorScalingLength          , factorScalingMass
   contains
     final     ::                                      cylindricalScalerDestructor
     procedure :: massTotal                         => cylindricalScalerMassTotal
     procedure :: density                           => cylindricalScalerDensity
     procedure :: densitySphericalAverage           => cylindricalScalerDensitySphericalAverage
     procedure :: surfaceDensity                    => cylindricalScalerSurfaceDensity
     procedure :: radiusHalfMass                    => cylindricalScalerRadiusHalfMass
     procedure :: massEnclosedBySphere              => cylindricalScalerMassEnclosedBySphere
     procedure :: potential                         => cylindricalScalerPotential
     procedure :: rotationCurve                     => cylindricalScalerRotationCurve
     procedure :: rotationCurveGradient             => cylindricalScalerRotationCurveGradient
     procedure :: surfaceDensityRadialMoment        => cylindricalScalerSurfaceDensityRadialMoment
     procedure :: acceleration                      => cylindricalScalerAcceleration
     procedure :: tidalTensor                       => cylindricalScalerTidalTensor
     procedure :: positionSample                    => cylindricalScalerPositionSample
     procedure :: isDimensionless                   => cylindricalScalerIsDimensionless
  end type massDistributionCylindricalScaler

  interface massDistributionCylindricalScaler
     !!{
     Constructors for the {\normalfont \ttfamily cylindricalScaler} mass distribution class.
     !!}
     module procedure cylindricalScalerConstructorParameters
     module procedure cylindricalScalerConstructorInternal
  end interface massDistributionCylindricalScaler

contains

  function cylindricalScalerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily cylindricalScaler} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionCylindricalScaler)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (massDistributionClass            ), pointer       :: massDistribution_
    double precision                                                   :: factorScalingLength , factorScalingMass

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
    class is (massDistributionCylindrical)
       self=massDistributionCylindricalScaler(factorScalingLength,factorScalingMass,massDistribution_)
    class default
       call Error_Report('a cylindrically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function cylindricalScalerConstructorParameters

  function cylindricalScalerConstructorInternal(factorScalingLength,factorScalingMass,massDistribution_) result(self)
    !!{
    Internal constructor for ``cylindricalScaler'' mass distribution class.
    !!}
    implicit none
    type            (massDistributionCylindricalScaler)                        :: self
    class           (massDistributionCylindrical      ), intent(in   ), target :: massDistribution_
    double precision                                   , intent(in   )         :: factorScalingLength, factorScalingMass
    !![
    <constructorAssign variables="factorScalingLength, factorScalingMass, *massDistribution_"/>
    !!]

    self%componentType=massDistribution_%componentType
    self%     massType=massDistribution_%     massType
    return
  end function cylindricalScalerConstructorInternal

  subroutine cylindricalScalerDestructor(self)
    !!{
    Destructor for the ``cylindricalScaler'' mass distribution class.
    !!}
    implicit none
    type(massDistributionCylindricalScaler), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine cylindricalScalerDestructor

  logical function cylindricalScalerIsDimensionless(self)
    !!{
    Return the dimensional status.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self

    cylindricalScalerIsDimensionless=.false.
    return
  end function cylindricalScalerIsDimensionless

  double precision function cylindricalScalerMassTotal(self,componentType,massType)
    !!{
    Return the total mass in a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout)           :: self
    type (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerMassTotal=+self%massDistribution_%massTotal        (               &
         &                                                               componentType, &
         &                                                               massType       &
         &                                                              )               &
         &                     *self                  %factorScalingMass
    return
  end function cylindricalScalerMassTotal

  double precision function cylindricalScalerDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout)              :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerDensity=+self%massDistribution_%density            (                           &
         &                                                               coordinatesScaled        , &
         &                                                               componentType            , &
         &                                                               massType                   &
         &                                                              )                           &
         &                   *self                  %factorScalingMass                              &
         &                   /self                  %factorScalingLength**3
    return
  end function cylindricalScalerDensity

  double precision function cylindricalScalerDensitySphericalAverage(self,radius,componentType,massType)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: radius
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerDensitySphericalAverage=+self%massDistribution_%densitySphericalAverage(                           &
         &                                                                                   +     radius               &
         &                                                                                   /self%factorScalingLength, &
         &                                                                                   componentType            , &
         &                                                                                   massType                   &
         &                                                                                  )                           &
         &                                   *self                  %factorScalingMass                                  &
         &                                   /self                  %factorScalingLength**3
    return
  end function cylindricalScalerDensitySphericalAverage

  double precision function cylindricalScalerRadiusHalfMass(self,componentType,massType)
    !!{
    Interface for cylindrically symmetric mass distribution half mass radii functions.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout)           :: self
    type (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerRadiusHalfMass=+self%massDistribution_%radiusHalfMass     (componentType,massType) &
         &                          *self                  %factorScalingLength
    return
  end function cylindricalScalerRadiusHalfMass

  double precision function cylindricalScalerMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a scaled cylindrical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout), target   :: self
    double precision                                   , intent(in   )           :: radius
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerMassEnclosedBySphere=+self%massDistribution_%massEnclosedBySphere(                           &
         &                                                                                   radius               &
         &                                                                             /self%factorScalingLength, &
         &                                                                             componentType            , &
         &                                                                             massType                   &
         &                                                                            )                           &
         &                                *self                  %factorScalingMass
    return
  end function cylindricalScalerMassEnclosedBySphere

  double precision function cylindricalScalerSurfaceDensity(self,coordinates,componentType,massType)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class(massDistributionCylindricalScaler), intent(inout)              :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerSurfaceDensity=+self%massDistribution_%surfaceDensity     (                           &
         &                                                                      coordinatesScaled        , &
         &                                                                      componentType            , &
         &                                                                      massType                   &
         &                                                                     )                           &
         &                          *self                  %factorScalingMass                              &
         &                          /self                  %factorScalingLength**2
    return
  end function cylindricalScalerSurfaceDensity

  double precision function cylindricalScalerRotationCurve(self,radius,componentType,massType)
    !!{
    Return the mid-plane rotation curve for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: radius
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerRotationCurve=+      self%massDistribution_%rotationCurve                  (                           &
         &                                                                                       +     radius               &
         &                                                                                       /self%factorScalingLength, &
         &                                                                                       componentType            , &
         &                                                                                       massType                   &
         &                                                                                      )                           &
         &                         *sqrt(                                                                                   &
         &                               +                       gravitationalConstantGalacticus                            &
         &                               *self                  %factorScalingMass                                          &
         &                               /self                  %factorScalingLength                                        &
         &                              )
    return
  end function cylindricalScalerRotationCurve

  double precision function cylindricalScalerRotationCurveGradient(self,radius,componentType,massType)
    !!{
    Return the mid-plane rotation curve gradient (specifically, $\mathrm{d}V^2_\mathrm{c}/\mathrm{d}r$) for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: radius
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerRotationCurveGradient=+self%massDistribution_%rotationCurveGradient(                           &
         &                                                                               +     radius               &
         &                                                                               /self%factorScalingLength, &
         &                                                                               componentType            , &
         &                                                                               massType                   &
         &                                                                              )                           &
         &                                 *     gravitationalConstantGalacticus                                    &
         &                                 *self%factorScalingMass                                                  &
         &                                 /self%factorScalingLength**2
    return
  end function cylindricalScalerRotationCurveGradient

  double precision function cylindricalScalerPotential(self,coordinates,componentType,massType)
    !!{
    Return the gravitational potential for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class(massDistributionCylindricalScaler), intent(inout)              :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerPotential=+self%massDistribution_%potential                      (                           &
         &                                                                             coordinatesScaled        , &
         &                                                                             componentType            , &
         &                                                                             massType                   &
         &                                                                            )                           &
         &                     *                       gravitationalConstantGalacticus                            &
         &                     *self                  %factorScalingMass                                          &
         &                     /self                  %factorScalingLength
    return
  end function cylindricalScalerPotential

  double precision function cylindricalScalerSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Compute radial moments of a scaled cylindrical distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: moment
    double precision                                   , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                            , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerSurfaceDensityRadialMoment=0.0d0
    !![
    <conditionalCall>
      <call>cylindricalScalerSurfaceDensityRadialMoment=self%massDistribution_%surfaceDensityRadialMoment(moment=moment,isInfinite=isInfinite,componentType=componentType,massType=massType{conditions})</call>
      <argument name="radiusMinimum" value="radiusMinimum/self%factorScalingLength" condition="present(radiusMinimum)"/>
      <argument name="radiusMaximum" value="radiusMaximum/self%factorScalingLength" condition="present(radiusMaximum)"/>
    </conditionalCall>
    !!]
    cylindricalScalerSurfaceDensityRadialMoment=+     cylindricalScalerSurfaceDensityRadialMoment &
         &                                      *self%factorScalingMass                           &
         &                                      /self%factorScalingLength**(2.0d0-moment)
    return
  end function cylindricalScalerSurfaceDensityRadialMoment

  function cylindricalScalerAcceleration(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, gravitationalConstantGalacticus, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                   , dimension(3  )             :: cylindricalScalerAcceleration
    class           (massDistributionCylindricalScaler), intent(inout)              :: self
    class           (coordinate                       ), intent(in   )              :: coordinates
    type            (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    class           (coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerAcceleration=+self%massDistribution_%acceleration       (                           &
         &                                                                    coordinatesScaled        , &
         &                                                                    componentType            , &
         &                                                                    massType                   &
         &                                                                   )                           &
         &                        *self                  %factorScalingMass                              &
         &                        /self                  %factorScalingLength**2                         &
         &                        *                       kilo                                           &
         &                        *                       gigaYear                                       &
         &                        /                       megaParsec                                     &
         &                        *                       gravitationalConstantGalacticus
     return
  end function cylindricalScalerAcceleration

  function cylindricalScalerTidalTensor(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type (tensorRank2Dimension3Symmetric   )                             :: cylindricalScalerTidalTensor
    class(massDistributionCylindricalScaler), intent(inout)              :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional    :: massType
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerTidalTensor=+self%massDistribution_%tidalTensor                    (                   &
         &                                                                               coordinatesScaled, &
         &                                                                               componentType    , &
         &                                                                               massType           &
         &                                                                              )                   &
         &                       *                       gravitationalConstantGalacticus                    &
         &                       *self                  %factorScalingMass                                  &
         &                       /self                  %factorScalingLength**3
    return
  end function cylindricalScalerTidalTensor
  
  function cylindricalScalerPositionSample(self,randomNumberGenerator_,componentType,massType)
    !!{
    Sample a position from a scaled cylindrical distribution.
    !!}
    implicit none
    double precision                                   , dimension(3)            :: cylindricalScalerPositionSample
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    class           (randomNumberGeneratorClass       ), intent(inout)           :: randomNumberGenerator_
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType

    cylindricalScalerPositionSample=+self%massDistribution_%positionSample     (randomNumberGenerator_,componentType,massType) &
         &                          *self                  %factorScalingLength
    return
  end function cylindricalScalerPositionSample
