!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
     !![
     <methods>
       <method method="unscaled" description="Return a pointer to the unscaled mass distribution."/>
     </methods>
     !!]
     final     ::                                            cylindricalScalerDestructor
     procedure :: unscaled                                => cylindricalScalerUnscaled
     procedure :: assumeMonotonicDecreasingSurfaceDensity => cylindricalScalerAssumeMonotonicDecreasingSurfaceDensity
     procedure :: massTotal                               => cylindricalScalerMassTotal
     procedure :: density                                 => cylindricalScalerDensity
     procedure :: densitySphericalAverage                 => cylindricalScalerDensitySphericalAverage
     procedure :: densityGradientRadial                   => cylindricalScalerDensityGradientRadial
     procedure :: surfaceDensity                          => cylindricalScalerSurfaceDensity
     procedure :: radiusHalfMass                          => cylindricalScalerRadiusHalfMass
     procedure :: massEnclosedBySphere                    => cylindricalScalerMassEnclosedBySphere
     procedure :: radiusEnclosingMass                     => cylindricalScalerRadiusEnclosingMass
     procedure :: radiusEnclosingDensity                  => cylindricalScalerRadiusEnclosingDensity
     procedure :: radiusEnclosingSurfaceDensity           => cylindricalScalerRadiusEnclosingSurfaceDensity
     procedure :: potentialIsAnalytic                     => cylindricalScalerPotentialIsAnalytic
     procedure :: potential                               => cylindricalScalerPotential
     procedure :: rotationCurve                           => cylindricalScalerRotationCurve
     procedure :: rotationCurveGradient                   => cylindricalScalerRotationCurveGradient
     procedure :: surfaceDensityRadialMoment              => cylindricalScalerSurfaceDensityRadialMoment
     procedure :: acceleration                            => cylindricalScalerAcceleration
     procedure :: tidalTensor                             => cylindricalScalerTidalTensor
     procedure :: positionSample                          => cylindricalScalerPositionSample
     procedure :: isDimensionless                         => cylindricalScalerIsDimensionless
  end type massDistributionCylindricalScaler

  interface massDistributionCylindricalScaler
     !!{
     Constructors for the \refClass{massDistributionCylindricalScaler} mass distribution class.
     !!}
     module procedure cylindricalScalerConstructorParameters
     module procedure cylindricalScalerConstructorInternal
  end interface massDistributionCylindricalScaler

contains

  function cylindricalScalerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionCylindricalScaler} mass distribution class which builds the object from a parameter
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
    Internal constructor for \refClass{massDistributionCylindricalScaler} mass distribution class.
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
    Destructor for the \refClass{massDistributionCylindricalScaler} mass distribution class.
    !!}
    implicit none
    type(massDistributionCylindricalScaler), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine cylindricalScalerDestructor

  function cylindricalScalerUnscaled(self) result(massDistribution_)
    !!{
    Return a pointer to the unscaled mass distribution.
    !!}
    implicit none
    class(massDistributionClass            ), pointer       :: massDistribution_
    class(massDistributionCylindricalScaler), intent(inout) :: self

    massDistribution_ => self%massDistribution_
    return
  end function cylindricalScalerUnscaled
  
  logical function cylindricalScalerAssumeMonotonicDecreasingSurfaceDensity(self) result(assumeMonotonicDecreasingSurfaceDensity)
    !!{
    Return true indicating that this distribution has a monotonically-decreasing surface density.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self

    assumeMonotonicDecreasingSurfaceDensity=self%massDistribution_%assumeMonotonicDecreasingSurfaceDensity()
    return
  end function cylindricalScalerAssumeMonotonicDecreasingSurfaceDensity
  
  logical function cylindricalScalerIsDimensionless(self)
    !!{
    Return the dimensional status.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self

    cylindricalScalerIsDimensionless=.false.
    return
  end function cylindricalScalerIsDimensionless

  double precision function cylindricalScalerMassTotal(self)
    !!{
    Return the total mass in a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout):: self

    cylindricalScalerMassTotal=+self%massDistribution_%massTotal        () &
         &                     *self                  %factorScalingMass
    return
  end function cylindricalScalerMassTotal

  double precision function cylindricalScalerDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates
    class(coordinate                       ), allocatable   :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerDensity=+self%massDistribution_%density            (coordinatesScaled) &
         &                   *self                  %factorScalingMass                      &
         &                   /self                  %factorScalingLength**3
    return
  end function cylindricalScalerDensity

  double precision function cylindricalScalerDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density gradient in the radial direction in a scaled cylindrical mass distribution.
    !!}
    implicit none
    class  (massDistributionCylindricalScaler), intent(inout), target      :: self
    class  (coordinate                       ), intent(in   )              :: coordinates
    logical                                   , intent(in   ), optional    :: logarithmic
    class  (coordinate                       )               , allocatable :: coordinatesScaled
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerDensityGradientRadial=+self%massDistribution_%densityGradientRadial(                   &
         &                                                                               coordinatesScaled, &
         &                                                                               logarithmic        &
         &                                                                              )
    if (.not.logarithmic)                                                                 &
         & cylindricalScalerDensityGradientRadial=+cylindricalScalerDensityGradientRadial &
         &                                        *self%factorScalingMass                 &
         &                                        /self%factorScalingLength**4
    return
  end function cylindricalScalerDensityGradientRadial

  double precision function cylindricalScalerDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius

    cylindricalScalerDensitySphericalAverage=+self%massDistribution_%densitySphericalAverage(                           &
         &                                                                                   +     radius               &
         &                                                                                   /self%factorScalingLength  &
         &                                                                                  )                           &
         &                                   *self                  %factorScalingMass                                  &
         &                                   /self                  %factorScalingLength**3
    return
  end function cylindricalScalerDensitySphericalAverage

  double precision function cylindricalScalerRadiusHalfMass(self)
    !!{
    Interface for cylindrically symmetric mass distribution half mass radii functions.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self

    cylindricalScalerRadiusHalfMass=+self%massDistribution_%radiusHalfMass     () &
         &                          *self                  %factorScalingLength
    return
  end function cylindricalScalerRadiusHalfMass

  double precision function cylindricalScalerMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a scaled cylindrical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout), target   :: self
    double precision                                   , intent(in   )           :: radius

    cylindricalScalerMassEnclosedBySphere=+self%massDistribution_%massEnclosedBySphere(                           &
         &                                                                                   radius               &
         &                                                                             /self%factorScalingLength  &
         &                                                                            )                           &
         &                                *self                  %factorScalingMass
    return
  end function cylindricalScalerMassEnclosedBySphere

  double precision function cylindricalScalerRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for cylindrically-scaled mass distributions.
    !!}    
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout), target   :: self
    double precision                                   , intent(in   ), optional :: mass, massFractional

    if (present(mass)) then
       radius=+self%massDistribution_%radiusEnclosingMass(mass/self%factorScalingMass,massFractional) &
            & *self                  %factorScalingLength
    else
       radius=+self%massDistribution_%radiusEnclosingMass(mass                       ,massFractional) &
            & *self                  %factorScalingLength
    end if
    return
  end function cylindricalScalerRadiusEnclosingMass

  double precision function cylindricalScalerRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for cylindrically-scaled mass distributions.
    !!}    
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout), target   :: self
    double precision                                   , intent(in   )           :: density
    double precision                                   , intent(in   ), optional :: radiusGuess

    if (present(radiusGuess)) then
       radius=+self%massDistribution_%radiusEnclosingDensity(                              &
            &                                                +     density                 &
            &                                                /self%factorScalingMass       &
            &                                                *self%factorScalingLength**3, &
            &                                                +     radiusGuess             &
            &                                                /self%factorScalingLength     &
            &                                               )                              &
            & *self                  %factorScalingLength
    else
       radius=+self%massDistribution_%radiusEnclosingDensity(                              &
            &                                                +     density                 &
            &                                                /self%factorScalingMass       &
            &                                                *self%factorScalingLength**3  &
            &                                               )                              &
            & *self                  %factorScalingLength
    end if
    return
  end function cylindricalScalerRadiusEnclosingDensity
  
  double precision function cylindricalScalerRadiusEnclosingSurfaceDensity(self,densitySurface,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given surface density for cylindrically-scaled mass distributions.
    !!}    
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout), target   :: self
    double precision                                   , intent(in   )           :: densitySurface
    double precision                                   , intent(in   ), optional :: radiusGuess
    
    if (present(radiusGuess)) then
       radius=+self%massDistribution_%radiusEnclosingSurfaceDensity(                              &
            &                                                       +     densitySurface          &
            &                                                       /self%factorScalingMass       &
            &                                                       *self%factorScalingLength**2, &
            &                                                       +     radiusGuess             &
            &                                                       /self%factorScalingLength     &
            &                                                      )                              &
            & *self                  %factorScalingLength
    else
       radius=+self%massDistribution_%radiusEnclosingSurfaceDensity(                              &
            &                                                       +     densitySurface          &
            &                                                       /self%factorScalingMass       &
            &                                                       *self%factorScalingLength**2  &
            &                                                      )                              &
            & *self                  %factorScalingLength
    end if
    return
  end function cylindricalScalerRadiusEnclosingSurfaceDensity
  
  double precision function cylindricalScalerSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates
    class(coordinate                       ), allocatable   :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerSurfaceDensity=+self%massDistribution_%surfaceDensity     (coordinatesScaled) &
         &                          *self                  %factorScalingMass                      &
         &                          /self                  %factorScalingLength**2
    return
  end function cylindricalScalerSurfaceDensity

  double precision function cylindricalScalerRotationCurve(self,radius)
    !!{
    Return the mid-plane rotation curve for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: radius

    cylindricalScalerRotationCurve=+      self%massDistribution_%rotationCurve      (                          &
         &                                                                           +     radius              &
         &                                                                           /self%factorScalingLength &
         &                                                                          )                          &
         &                         *sqrt(                                                                      &
         &                               +self                  %factorScalingMass                             &
         &                               /self                  %factorScalingLength                           &
         &                              )
    if (self%massDistribution_%isDimensionless())                               &
         & cylindricalScalerRotationCurve=+cylindricalScalerRotationCurve       &
         &                                *sqrt(gravitationalConstant_internal)
    return
  end function cylindricalScalerRotationCurve

  double precision function cylindricalScalerRotationCurveGradient(self,radius)
    !!{
    Return the mid-plane rotation curve gradient (specifically, $\mathrm{d}V^2_\mathrm{c}/\mathrm{d}r$) for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius

    cylindricalScalerRotationCurveGradient=+self%massDistribution_%rotationCurveGradient(                          &
         &                                                                               +     radius              &
         &                                                                               /self%factorScalingLength &
         &                                                                              )                          &
         &                                 *self%factorScalingMass                                                 &
         &                                 /self%factorScalingLength**2
    if (self%massDistribution_%isDimensionless())                                         &
         & cylindricalScalerRotationCurveGradient=+cylindricalScalerRotationCurveGradient &
         &                                        *gravitationalConstant_internal
    return
  end function cylindricalScalerRotationCurveGradient

  logical function cylindricalScalerPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return if the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self

    isAnalytic=self%massDistribution_%potentialIsAnalytic()
    return
  end function cylindricalScalerPotentialIsAnalytic

  double precision function cylindricalScalerPotential(self,coordinates,status)
    !!{
    Return the gravitational potential for a scaled cylindrical distribution.
    !!}
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(massDistributionCylindricalScaler), intent(inout), target      :: self
    class(coordinate                       ), intent(in   )              :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional    :: status
    class(coordinate                       )               , allocatable :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerPotential=+self%massDistribution_%potential          (                   &
         &                                                                 coordinatesScaled, &
         &                                                                 status             &
         &                                                                )                   &
         &                     *self                  %factorScalingMass                      &
         &                     /self                  %factorScalingLength
    if (self%massDistribution_%isDimensionless())                     &
         & cylindricalScalerPotential=+cylindricalScalerPotential     &
         &                            *gravitationalConstant_internal
    return
  end function cylindricalScalerPotential

  double precision function cylindricalScalerSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Compute radial moments of a scaled cylindrical distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout)           :: self
    double precision                                   , intent(in   )           :: moment
    double precision                                   , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                            , intent(  out), optional :: isInfinite
 
    cylindricalScalerSurfaceDensityRadialMoment=0.0d0
    !![
    <conditionalCall>
      <call>cylindricalScalerSurfaceDensityRadialMoment=self%massDistribution_%surfaceDensityRadialMoment(moment=moment,isInfinite=isInfinite{conditions})</call>
      <argument name="radiusMinimum" value="radiusMinimum/self%factorScalingLength" condition="present(radiusMinimum)"/>
      <argument name="radiusMaximum" value="radiusMaximum/self%factorScalingLength" condition="present(radiusMaximum)"/>
    </conditionalCall>
    !!]
    cylindricalScalerSurfaceDensityRadialMoment=+     cylindricalScalerSurfaceDensityRadialMoment &
         &                                      *self%factorScalingMass                           &
         &                                      /self%factorScalingLength**(2.0d0-moment)
    return
  end function cylindricalScalerSurfaceDensityRadialMoment

  function cylindricalScalerAcceleration(self,coordinates)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                   , dimension(3)  :: cylindricalScalerAcceleration
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    class           (coordinate                       ), intent(in   ) :: coordinates
    class           (coordinate                       ), allocatable   :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerAcceleration=+self%massDistribution_%acceleration       (coordinatesScaled) &
         &                        *self                  %factorScalingMass                      &
         &                        /self                  %factorScalingLength**2
    if (self%massDistribution_%isDimensionless())                        &
         & cylindricalScalerAcceleration=+cylindricalScalerAcceleration  &
         &                               *kilo                           &
         &                               *gigaYear                       &
         &                               /megaParsec                     &
         &                               *gravitationalConstant_internal
       return
  end function cylindricalScalerAcceleration

  function cylindricalScalerTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for a scaled cylindrical distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type (tensorRank2Dimension3Symmetric   )                :: cylindricalScalerTidalTensor
    class(massDistributionCylindricalScaler), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates
    class(coordinate                       ), allocatable   :: coordinatesScaled

    call coordinates%scale(1.0d0/self%factorScalingLength,coordinatesScaled)
    cylindricalScalerTidalTensor=+self%massDistribution_%tidalTensor           (coordinatesScaled) &
         &                       *self                  %factorScalingMass                         &
         &                       /self                  %factorScalingLength**3
    if (self%massDistribution_%isDimensionless())                       &
         & cylindricalScalerTidalTensor=+cylindricalScalerTidalTensor   &
         &                              *gravitationalConstant_internal
     return
  end function cylindricalScalerTidalTensor
  
  function cylindricalScalerPositionSample(self,randomNumberGenerator_)
    !!{
    Sample a position from a scaled cylindrical distribution.
    !!}
    implicit none
    double precision                                   , dimension(3)  :: cylindricalScalerPositionSample
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    class           (randomNumberGeneratorClass       ), intent(inout) :: randomNumberGenerator_

    cylindricalScalerPositionSample=+self%massDistribution_%positionSample     (randomNumberGenerator_) &
         &                          *self                  %factorScalingLength
    return
  end function cylindricalScalerPositionSample
