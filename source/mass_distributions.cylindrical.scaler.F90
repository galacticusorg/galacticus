!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
      \rho(\mathbf{x}) = f_\mathrm{M} \rho^\prime(\mathbf{x}/f_\mathrm{r}),
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

  double precision function cylindricalScalerDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates
 
    cylindricalScalerDensity=+self%massDistribution_%density            (                          &
         &                                                                     coordinates         &
         &                                                               /self%factorScalingLength &
         &                                                              )                          &
         &                   *self                  %factorScalingMass                             &
         &                   /self                  %factorScalingLength**3
    return
  end function cylindricalScalerDensity

  double precision function cylindricalScalerDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius

    cylindricalScalerDensitySphericalAverage=+self%massDistribution_%densitySphericalAverage(                          &
         &                                                                                   +     radius              &
         &                                                                                   /self%factorScalingLength &
         &                                                                                  )                          &
         &                                   *self                  %factorScalingMass                                 &
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
    class           (massDistributionCylindricalScaler), intent(inout), target :: self
    double precision                                   , intent(in   )         :: radius

    cylindricalScalerMassEnclosedBySphere=+self%massDistribution_%massEnclosedBySphere(                          &
         &                                                                                   radius              &
         &                                                                             /self%factorScalingLength &
         &                                                                            )                          &
         &                                *self                  %factorScalingMass
    return
  end function cylindricalScalerMassEnclosedBySphere

  double precision function cylindricalScalerSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a scaled cylindrical distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    class           (coordinate                       ), intent(in   ) :: coordinates

    cylindricalScalerSurfaceDensity=+self%massDistribution_%surfaceDensity     (                          &
         &                                                                            coordinates         &
         &                                                                      /self%factorScalingLength &
         &                                                                     )                          &
         &                          *self                  %factorScalingMass                             &
         &                          /self                  %factorScalingLength**2
    return
  end function cylindricalScalerSurfaceDensity

  double precision function cylindricalScalerRotationCurve(self,radius)
    !!{
    Return the mid-plane rotation curve for a scaled cylindrical distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius

    cylindricalScalerRotationCurve=+      self%massDistribution_%rotationCurve      (                          &
         &                                                                           +     radius              &
         &                                                                           /self%factorScalingLength &
         &                                                                          )                          &
         &                         *sqrt(                                                                      &
         &                               +self                  %factorScalingMass                             &
         &                               /self                  %factorScalingLength                           &
         &                              )
    return
  end function cylindricalScalerRotationCurve

  double precision function cylindricalScalerRotationCurveGradient(self,radius)
    !!{
    Return the mid-plane rotation curve gradient for a scaled cylindrical distribution.
    !!}
    implicit none
    class           (massDistributionCylindricalScaler), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius

    cylindricalScalerRotationCurveGradient=+      self%massDistribution_%rotationCurveGradient(                          &
         &                                                                                     +     radius              &
         &                                                                                     /self%factorScalingLength &
         &                                                                                    )                          &
         &                                 *sqrt(                                                                        &
         &                                       +self                  %factorScalingMass                               &
         &                                       /self                  %factorScalingLength                             &
         &                                      )                                                                        &
         &                                 /      self                  %factorScalingLength   
    return
  end function cylindricalScalerRotationCurveGradient

  double precision function cylindricalScalerPotential(self,coordinates)
    !!{
    Return the gravitational potential for a scaled cylindrical distribution.
    !!}
    implicit none
    class(massDistributionCylindricalScaler), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates

    cylindricalScalerPotential=+self%massDistribution_%potential          (                          &
         &                                                                       coordinates         &
         &                                                                 /self%factorScalingLength &
         &                                                                )                          &
         &                     *self                  %factorScalingMass                             &
         &                     /self                  %factorScalingLength
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
    implicit none
    double precision                                   , dimension(3  ) :: cylindricalScalerAcceleration
    class           (massDistributionCylindricalScaler), intent(inout)  :: self
    class           (coordinate                       ), intent(in   )  :: coordinates

    cylindricalScalerAcceleration=+self%massDistribution_%acceleration       (                          &
         &                                                                          coordinates         &
         &                                                                    /self%factorScalingLength &
         &                                                                   )                          &
         &                        *self                  %factorScalingMass                             &
         &                        /self                  %factorScalingLength**2
    return
  end function cylindricalScalerAcceleration

  function cylindricalScalerTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for a scaled cylindrical distribution.
    !!}
    implicit none
    type (tensorRank2Dimension3Symmetric   )                 :: cylindricalScalerTidalTensor
    class(massDistributionCylindricalScaler), intent(inout)  :: self
    class(coordinate                       ), intent(in   )  :: coordinates

    cylindricalScalerTidalTensor=+self%massDistribution_%tidalTensor        (                          &
         &                                                                         coordinates         &
         &                                                                   /self%factorScalingLength &
         &                                                                  )                          &
         &                       *self                  %factorScalingMass                             &
         &                       /self                  %factorScalingLength**3
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
