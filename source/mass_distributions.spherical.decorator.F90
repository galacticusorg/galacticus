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
  Implements an abstract spherical mass distribution decorator class.
  !!}

  !![
  <massDistribution name="massDistributionSphericalDecorator" abstract="yes">
   <description>
     An abstract mass distribution class for decorators of other mass distributions. ``Fallthrough'' functions are provided that
     all the decorated class or numerical solutions to be used.
   </description>
  </massDistribution>
  !!]
  type, abstract, extends(massDistributionSpherical) :: massDistributionSphericalDecorator
     !!{
     Implementation of a decorator spherical mass distribution.
     !!}
     private
     type (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     class(massDistributionSpherical        ), pointer :: massDistribution_ => null()
  contains
     !![
     <methods>
       <method method="massEnclosedBySphereNonAnalytic"              description="Compute mass enclosed by a sphere for non-analytic cases."                                />
       <method method="radiusEnclosingMassNonAnalytic"               description="Compute radius enclosing a mass for non-analytic cases."                                  />
       <method method="densityGradientRadialNonAnalytic"             description="Compute radial density gradient for non-analytic cases."                                  />
       <method method="densityRadialMomentNonAnalytic"               description="Compute radial density moment for non-analytic cases."                                    />
       <method method="radiusEnclosingDensityNonAnalytic"            description="Compute radius enclosing a mean density for non-analytic cases."                          />
       <method method="radiusFromSpecificAngularMomentumNonAnalytic" description="Compute radius from specific angular momentum for non-analytic cases."                    />
       <method method="fourierTransformNonAnalytic"                  description="Compute Fourier transform for non-analytic cases."                                        />
       <method method="radiusFreefallNonAnalytic"                    description="Compute freefall radius for non-analytic cases."                                          />
       <method method="radiusFreefallIncreaseRateNonAnalytic"        description="Compute freefall radius growth rate for non-analytic cases."                              />
       <method method="potentialNonAnalytic"                         description="Compute gravitational potential for non-analytic cases."                                  />
       <method method="energyPotentialNonAnalytic"                   description="Compute gravitational potential energy for non-analytic cases."                           />
       <method method="energyKinetciNonAnalytic"                     description="Compute kinetic energy for non-analytic cases."                                           />
       <method method="useUndecorated"                               description="Return true if the undecorated solution (instead of a numerical solution) should be used."/>
     </methods>
     !!]
     procedure :: massEnclosedBySphere                         => sphericalDecoratorMassEnclosedBySphere
     procedure :: radiusEnclosingMass                          => sphericalDecoratorRadiusEnclosingMass
     procedure :: densityGradientRadial                        => sphericalDecoratorDensityGradientRadial
     procedure :: densityRadialMoment                          => sphericalDecoratorDensityRadialMoment
     procedure :: radiusEnclosingDensity                       => sphericalDecoratorRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum            => sphericalDecoratorRadiusFromSpecificAngularMomentum
     procedure :: fourierTransform                             => sphericalDecoratorFourierTransform
     procedure :: radiusFreefall                               => sphericalDecoratorRadiusFreefall
     procedure :: radiusFreefallIncreaseRate                   => sphericalDecoratorRadiusFreefallIncreaseRate 
     procedure :: potential                                    => sphericalDecoratorPotential
     procedure :: energyPotential                              => sphericalDecoratorEnergyPotential
     procedure :: energyKinetic                                => sphericalDecoratorEnergyKinetic
     procedure :: massEnclosedBySphereNonAnalytic              => sphericalDecoratorMassEnclosedBySphereNonAnalytic
     procedure :: radiusEnclosingMassNonAnalytic               => sphericalDecoratorRadiusEnclosingMassNonAnalytic
     procedure :: densityGradientRadialNonAnalytic             => sphericalDecoratorDensityGradientRadialNonAnalytic
     procedure :: densityRadialMomentNonAnalytic               => sphericalDecoratorDensityRadialMomentNonAnalytic
     procedure :: radiusEnclosingDensityNonAnalytic            => sphericalDecoratorRadiusEnclosingDensityNonAnalytic
     procedure :: radiusFromSpecificAngularMomentumNonAnalytic => sphericalDecoratorRadiusFromSpecificAngularMomentumNonAnalytic
     procedure :: fourierTransformNonAnalytic                  => sphericalDecoratorFourierTransformNonAnalytic
     procedure :: radiusFreefallNonAnalytic                    => sphericalDecoratorRadiusFreefallNonAnalytic
     procedure :: radiusFreefallIncreaseRateNonAnalytic        => sphericalDecoratorRadiusFreefallIncreaseRateNonAnalytic 
     procedure :: potentialNonAnalytic                         => sphericalDecoratorPotentialNonAnalytic
     procedure :: energyPotentialNonAnalytic                   => sphericalDecoratorEnergyPotentialNonAnalytic
     procedure :: energyKineticNonAnalytic                     => sphericalDecoratorEnergyKineticNonAnalytic
     procedure :: useUndecorated                               => sphericalDecoratorUseUndecorated
  end type massDistributionSphericalDecorator

contains

  logical function sphericalDecoratorUseUndecorated(self) result(useUndecorated)
    !!{
    Determines whether to use the undecorated solution.
    !!}
    implicit none
    class(massDistributionSphericalDecorator), intent(inout) :: self

    useUndecorated=self%nonAnalyticSolver == nonAnalyticSolversFallThrough
    return
  end function sphericalDecoratorUseUndecorated

  double precision function sphericalDecoratorMassEnclosedBySphere(self,radius,componentType,massType) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for decorator mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: radius
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    mass=self%massEnclosedBySphereNonAnalytic(radius,componentType,massType)
    return
  end function sphericalDecoratorMassEnclosedBySphere
  
  double precision function sphericalDecoratorMassEnclosedBySphereNonAnalytic(self,radius,componentType,massType) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for decorator mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: radius
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       mass=self%massDistribution_%massEnclosedBySphere         (radius,componentType,massType)
    else
       mass=self                  %massEnclosedBySphereNumerical(radius,componentType,massType)
    end if
    return
  end function sphericalDecoratorMassEnclosedBySphereNonAnalytic
  
  double precision function sphericalDecoratorDensityGradientRadial(self,coordinates,logarithmic,componentType,massType) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalDecorator), intent(inout), target   :: self
    class  (coordinate                        ), intent(in   )           :: coordinates
    logical                                    , intent(in   ), optional :: logarithmic
    type   (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType           ), intent(in   ), optional :: massType

    densityGradient=self%densityGradientRadialNonAnalytic(coordinates,logarithmic,componentType,massType)
    return
  end function sphericalDecoratorDensityGradientRadial
  
  double precision function sphericalDecoratorDensityGradientRadialNonAnalytic(self,coordinates,logarithmic,componentType,massType) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalDecorator), intent(inout), target   :: self
    class  (coordinate                        ), intent(in   )           :: coordinates
    logical                                    , intent(in   ), optional :: logarithmic
    type   (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       densityGradient=self%massDistribution_%densityGradientRadial         (coordinates,logarithmic,componentType,massType)
    else
       densityGradient=self                  %densityGradientRadialNumerical(coordinates,logarithmic,componentType,massType)
    end if 
    return
  end function sphericalDecoratorDensityGradientRadialNonAnalytic
  
  double precision function sphericalDecoratorRadiusEnclosingDensity(self,density,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mean density for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: density
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    radius=self%radiusEnclosingDensityNonAnalytic(density,componentType,massType)
    return
  end function sphericalDecoratorRadiusEnclosingDensity
  
  double precision function sphericalDecoratorRadiusEnclosingDensityNonAnalytic(self,density,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mean density for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: density
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusEnclosingDensity         (density,componentType,massType)
    else
       radius=self                  %radiusEnclosingDensityNumerical(density,componentType,massType)
    end if
    return
  end function sphericalDecoratorRadiusEnclosingDensityNonAnalytic
  
  double precision function sphericalDecoratorRadiusEnclosingMass(self,mass,massFractional,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for heated spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass         , massFractional
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    radius=self%radiusEnclosingMassNonAnalytic(mass,massFractional,componentType,massType)
    return
  end function sphericalDecoratorRadiusEnclosingMass

  double precision function sphericalDecoratorRadiusEnclosingMassNonAnalytic(self,mass,massFractional,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for heated spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass         , massFractional
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusEnclosingMass         (mass,massFractional,componentType,massType)
    else
       radius=self                  %radiusEnclosingMassNumerical(mass,massFractional,componentType,massType)
    end if
    return
  end function sphericalDecoratorRadiusEnclosingMassNonAnalytic

  double precision function sphericalDecoratorRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific,componentType,massType) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: angularMomentumSpecific
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    radius=self%radiusFromSpecificAngularMomentumNonAnalytic(angularMomentumSpecific,componentType,massType)
    return
  end function sphericalDecoratorRadiusFromSpecificAngularMomentum
  
  double precision function sphericalDecoratorRadiusFromSpecificAngularMomentumNonAnalytic(self,angularMomentumSpecific,componentType,massType) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: angularMomentumSpecific
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusFromSpecificAngularMomentum         (angularMomentumSpecific,componentType,massType)
    else
       radius=self                  %radiusFromSpecificAngularMomentumNumerical(angularMomentumSpecific,componentType,massType)
    end if
    return
  end function sphericalDecoratorRadiusFromSpecificAngularMomentumNonAnalytic
  
  double precision function sphericalDecoratorDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType) result(densityRadialMoment)
    !!{
    Returns a radial density moment for the decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: moment
    double precision                                    , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                             , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    densityRadialMoment=self%densityRadialMomentNonAnalytic(moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    return
  end function sphericalDecoratorDensityRadialMoment

  double precision function sphericalDecoratorDensityRadialMomentNonAnalytic(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType) result(densityRadialMoment)
    !!{
    Returns a radial density moment for the decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: moment
    double precision                                    , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                             , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       densityRadialMoment=self%massDistribution_%densityRadialMoment         (moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    else
       densityRadialMoment=self                  %densityRadialMomentNumerical(moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    end if
    return
  end function sphericalDecoratorDensityRadialMomentNonAnalytic

  double precision function sphericalDecoratorPotential(self,coordinates,componentType,massType,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalDecorator  ), intent(inout), target   :: self
    class(coordinate                          ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType        ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType             ), intent(in   ), optional :: massType
    type (enumerationStructureErrorCodeType   ), intent(  out), optional :: status

    potential=self%potentialNonAnalytic(coordinates,componentType,massType,status)
    return
  end function sphericalDecoratorPotential

  double precision function sphericalDecoratorPotentialNonAnalytic(self,coordinates,componentType,massType,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalDecorator  ), intent(inout), target   :: self
    class(coordinate                          ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType        ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType             ), intent(in   ), optional :: massType
    type (enumerationStructureErrorCodeType   ), intent(  out), optional :: status

    if (self%useUndecorated()) then
       potential=self%massDistribution_%potential         (coordinates,componentType,massType,status)
    else
       potential=self                  %potentialNumerical(coordinates,componentType,massType,status)
    end if
    return
  end function sphericalDecoratorPotentialNonAnalytic

  double precision function sphericalDecoratorFourierTransform(self,radiusOuter,wavenumber,componentType,massType) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: radiusOuter  , wavenumber
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType
    
    fourierTransform=self%fourierTransformNonAnalytic(radiusOuter,wavenumber,componentType,massType)
    return
  end function sphericalDecoratorFourierTransform
  
  double precision function sphericalDecoratorFourierTransformNonAnalytic(self,radiusOuter,wavenumber,componentType,massType) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: radiusOuter  , wavenumber
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType
    
    if (self%useUndecorated()) then
       fourierTransform=self%massDistribution_%fourierTransform         (radiusOuter,wavenumber,componentType,massType)
    else
       fourierTransform=self                  %fourierTransformNumerical(radiusOuter,wavenumber,componentType,massType)
    end if
    return
  end function sphericalDecoratorFourierTransformNonAnalytic
  
  double precision function sphericalDecoratorRadiusFreefall(self,time,componentType,massType) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: time
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    radius=self%radiusFreefallNonAnalytic(time,componentType,massType)
    return
  end function sphericalDecoratorRadiusFreefall
  
  double precision function sphericalDecoratorRadiusFreefallNonAnalytic(self,time,componentType,massType) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: time
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusFreefall         (time,componentType,massType)
    else
       radius=self                  %radiusFreefallNumerical(time,componentType,massType)
    end if
    return
  end function sphericalDecoratorRadiusFreefallNonAnalytic
  
  double precision function sphericalDecoratorRadiusFreefallIncreaseRate(self,time,componentType,massType) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: time
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    radiusIncreaseRate=self%radiusFreefallIncreaseRateNonAnalytic(time,componentType,massType)
    return
  end function sphericalDecoratorRadiusFreefallIncreaseRate

  double precision function sphericalDecoratorRadiusFreefallIncreaseRateNonAnalytic(self,time,componentType,massType) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: time
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%useUndecorated()) then
       radiusIncreaseRate=self%massDistribution_%radiusFreefallIncreaseRate         (time,componentType,massType)
    else
       radiusIncreaseRate=self                  %radiusFreefallIncreaseRateNumerical(time,componentType,massType)
    end if
    return
  end function sphericalDecoratorRadiusFreefallIncreaseRateNonAnalytic

  double precision function sphericalDecoratorEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter

    energy=self%energyPotentialNonAnalytic(radiusOuter)
    return
  end function sphericalDecoratorEnergyPotential

  double precision function sphericalDecoratorEnergyPotentialNonAnalytic(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter

    if (self%useUndecorated()) then
       energy=self%massDistribution_%energyPotential         (radiusOuter)
    else
       energy=self                  %energyPotentialNumerical(radiusOuter)
    end if
    return
  end function sphericalDecoratorEnergyPotentialNonAnalytic

  double precision function sphericalDecoratorEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter
    class           (massDistributionClass             ), intent(inout) :: massDistributionEmbedding

    energy=self%energyKineticNonAnalytic(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalDecoratorEnergyKinetic

  double precision function sphericalDecoratorEnergyKineticNonAnalytic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter
    class           (massDistributionClass             ), intent(inout) :: massDistributionEmbedding

    if (self%useUndecorated()) then
       energy=self%massDistribution_%energyKinetic         (radiusOuter,massDistributionEmbedding)
    else
       energy=self                  %energyKineticNumerical(radiusOuter,massDistributionEmbedding)
    end if
    return
  end function sphericalDecoratorEnergyKineticNonAnalytic
