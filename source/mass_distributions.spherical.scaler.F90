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
  Implementation of a scaling spherical mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionSphericalScaler">
    <description>
      A mass distribution class for scaling spherical mass distributions. Specifically, the density at position $\mathbf{x}$ is
      given by
      \begin{equation}
      \rho(\mathbf{x}) = f_\mathrm{M} \rho^\prime(\mathbf{x}/f_\mathrm{r}),
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
     final     ::                          sphericalScalerDestructor
     procedure :: density               => sphericalScalerDensity
     procedure :: densityGradientRadial => sphericalScalerDensityGradientRadial
     procedure :: densityRadialMoment   => sphericalScalerDensityRadialMoment
     procedure :: massEnclosedBySphere  => sphericalScalerMassEnclosedBySphere
     procedure :: potential             => sphericalScalerPotential
     procedure :: radiusHalfMass        => sphericalScalerRadiusHalfMass
     procedure :: tidalTensor           => sphericalScalerTidalTensor
     procedure :: radiusEnclosingMass   => sphericalScalerRadiusEnclosingMass
     procedure :: positionSample        => sphericalScalerPositionSample
     procedure :: isDimensionless       => sphericalScalerIsDimensionless
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

  logical function sphericalScalerIsDimensionless(self)
    !!{
    Return the dimensional status.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout) :: self

    sphericalScalerIsDimensionless=.false.
    return
  end function sphericalScalerIsDimensionless

  double precision function sphericalScalerDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout) :: self
    class(coordinate                     ), intent(in   ) :: coordinates

    sphericalScalerDensity=+self%massDistribution_%density          (                          &
         &                                                                 coordinates         &
         &                                                           /self%factorScalingLength &
         &                                                          )                          &
         &                 *self                  %factorScalingMass                           &
         &                 /self                  %factorScalingLength**3
    return
  end function sphericalScalerDensity

  double precision function sphericalScalerDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density gradient in the radial direction in a scaled spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalScaler), intent(inout)           :: self
    class  (coordinate                     ), intent(in   )           :: coordinates
    logical                                 , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    sphericalScalerDensityGradientRadial=+self%massDistribution_%densityGradientRadial(                           &
         &                                                                                   coordinates          &
         &                                                                             /self%factorScalingLength, &
         &                                                                             logarithmic                &
         &                                                                            )
    if (.not.logarithmic)                                                             &
         & sphericalScalerDensityGradientRadial=+sphericalScalerDensityGradientRadial &
         &                                      *self%factorScalingMass               &
         &                                      /self%factorScalingLength**4
    return
  end function sphericalScalerDensityGradientRadial

  double precision function sphericalScalerMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target :: self
    double precision                                 , intent(in   )         :: radius

    sphericalScalerMassEnclosedBySphere=+self%massDistribution_%massEnclosedBySphere(                          &
         &                                                                                 radius              &
         &                                                                           /self%factorScalingLength &
         &                                                                          )                          &
         &                              *self                  %factorScalingMass
    return
  end function sphericalScalerMassEnclosedBySphere

  double precision function sphericalScalerPotential(self,coordinates)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout) :: self
    class(coordinate                     ), intent(in   ) :: coordinates

    sphericalScalerPotential=+self%massDistribution_%potential         (                          &
         &                                                                    coordinates         &
         &                                                              /self%factorScalingLength &
         &                                                             )                          &
         &                  *self                  %factorScalingMass                             &
         &                  /self                  %factorScalingLength
    return
  end function sphericalScalerPotential

  double precision function sphericalScalerDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes radial moments of the density in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout)           :: self
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                          , intent(  out), optional :: isInfinite
    
    sphericalScalerDensityRadialMoment=0.0d0
    !![
    <conditionalCall>
      <call>sphericalScalerDensityRadialMoment=self%massDistribution_%densityRadialMoment(moment=moment,isInfinite=isInfinite{conditions})</call>
      <argument name="radiusMinimum" value="radiusMinimum/self%factorScalingLength" condition="present(radiusMinimum)"/>
      <argument name="radiusMaximum" value="radiusMaximum/self%factorScalingLength" condition="present(radiusMaximum)"/>
    </conditionalCall>
    !!]
    sphericalScalerDensityRadialMoment=+     sphericalScalerDensityRadialMoment  &
         &                             *self%factorScalingMass                   &
         &                             /self%factorScalingLength**(3.0d0-moment)
    return    
  end function sphericalScalerDensityRadialMoment

  double precision function sphericalScalerRadiusEnclosingMass(self,mass)
    !!{
    Computes the radius enclosing a given mass in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalScaler), intent(inout), target :: self
    double precision                                 , intent(in   )         :: mass

   sphericalScalerRadiusEnclosingMass=+self%massDistribution_%radiusEnclosingMass(                        &
         &                                                                        +     mass              &
         &                                                                        /self%factorScalingMass &
         &                                                                       )                        &
         &                            *self                  %factorScalingLength
    return
  end function sphericalScalerRadiusEnclosingMass

  double precision function sphericalScalerRadiusHalfMass(self)
    !!{
    Computes the half-mass radius in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalScaler), intent(inout) :: self

    sphericalScalerRadiusHalfMass=+self%massDistribution_%radiusHalfMass     () &
         &                        *self                  %factorScalingLength
    return
  end function sphericalScalerRadiusHalfMass

  function sphericalScalerAcceleration(self,coordinates)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    implicit none
    double precision                                 , dimension(3)  :: sphericalScalerAcceleration
    class           (massDistributionSphericalScaler), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates

    sphericalScalerAcceleration=+self%massDistribution_%acceleration   (                          &
         &                                                                    coordinates         &
         &                                                              /self%factorScalingLength &
         &                                                             )                          &
         &                  *self                  %factorScalingMass                             &
         &                  /self                  %factorScalingLength**2
    return
  end function sphericalScalerAcceleration

  function sphericalScalerTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    type (tensorRank2Dimension3Symmetric )                :: sphericalScalerTidalTensor
    class(massDistributionSphericalScaler), intent(inout) :: self
    class(coordinate                     ), intent(in   ) :: coordinates

    sphericalScalerTidalTensor=+self%massDistribution_%tidalTensor     (                          &
         &                                                                    coordinates         &
         &                                                              /self%factorScalingLength &
         &                                                             )                          &
         &                  *self                  %factorScalingMass                             &
         &                  /self                  %factorScalingLength**3
    return
  end function sphericalScalerTidalTensor
  
  function sphericalScalerPositionSample(self,randomNumberGenerator_)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution in a scaled spherical mass distribution.
    !!}
    implicit none
    double precision                                 , dimension(3)  :: sphericalScalerPositionSample
    class           (massDistributionSphericalScaler), intent(inout) :: self
    class           (randomNumberGeneratorClass     ), intent(inout) :: randomNumberGenerator_

    sphericalScalerPositionSample=+self%massDistribution_%positionSample     (randomNumberGenerator_) &
         &                        *self                  %factorScalingLength
    return
  end function sphericalScalerPositionSample
