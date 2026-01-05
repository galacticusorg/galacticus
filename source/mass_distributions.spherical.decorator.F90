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
       <method method="energyNonAnalytic"                            description="Compute total energy for non-analytic cases."                                             />
       <method method="energyPotentialNonAnalytic"                   description="Compute gravitational potential energy for non-analytic cases."                           />
       <method method="energyKineticNonAnalytic"                     description="Compute kinetic energy for non-analytic cases."                                           />
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
     procedure :: potentialIsAnalytic                          => sphericalDecoratorPotentialIsAnalytic
     procedure :: potential                                    => sphericalDecoratorPotential
     procedure :: energy                                       => sphericalDecoratorEnergy
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
     procedure :: energyNonAnalytic                            => sphericalDecoratorEnergyNonAnalytic
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

  double precision function sphericalDecoratorMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for decorator mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target :: self
    double precision                                    , intent(in   )         :: radius

    mass=self%massEnclosedBySphereNonAnalytic(radius)
    return
  end function sphericalDecoratorMassEnclosedBySphere
  
  double precision function sphericalDecoratorMassEnclosedBySphereNonAnalytic(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for decorator mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target :: self
    double precision                                    , intent(in   )         :: radius
 
    if (self%useUndecorated()) then
       mass=self%massDistribution_%massEnclosedBySphere         (radius)
    else
       mass=self                  %massEnclosedBySphereNumerical(radius)
    end if
    return
  end function sphericalDecoratorMassEnclosedBySphereNonAnalytic
  
  double precision function sphericalDecoratorDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalDecorator), intent(inout), target   :: self
    class  (coordinate                        ), intent(in   )           :: coordinates
    logical                                    , intent(in   ), optional :: logarithmic

    densityGradient=self%densityGradientRadialNonAnalytic(coordinates,logarithmic)
    return
  end function sphericalDecoratorDensityGradientRadial
  
  double precision function sphericalDecoratorDensityGradientRadialNonAnalytic(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalDecorator), intent(inout), target   :: self
    class  (coordinate                        ), intent(in   )           :: coordinates
    logical                                    , intent(in   ), optional :: logarithmic

    if (self%useUndecorated()) then
       densityGradient=self%massDistribution_%densityGradientRadial         (coordinates,logarithmic)
    else
       densityGradient=self                  %densityGradientRadialNumerical(coordinates,logarithmic)
    end if 
    return
  end function sphericalDecoratorDensityGradientRadialNonAnalytic
  
  double precision function sphericalDecoratorRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: density
    double precision                                    , intent(in   ), optional :: radiusGuess

    radius=self%radiusEnclosingDensityNonAnalytic(density,radiusGuess)
    return
  end function sphericalDecoratorRadiusEnclosingDensity
  
  double precision function sphericalDecoratorRadiusEnclosingDensityNonAnalytic(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: density
    double precision                                    , intent(in   ), optional :: radiusGuess
 
    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusEnclosingDensity         (density,radiusGuess)
    else
       radius=self                  %radiusEnclosingDensityNumerical(density,radiusGuess)
    end if
    return
  end function sphericalDecoratorRadiusEnclosingDensityNonAnalytic
  
  double precision function sphericalDecoratorRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for heated spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass, massFractional

    radius=self%radiusEnclosingMassNonAnalytic(mass,massFractional)
    return
  end function sphericalDecoratorRadiusEnclosingMass

  double precision function sphericalDecoratorRadiusEnclosingMassNonAnalytic(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for heated spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass, massFractional

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusEnclosingMass         (mass,massFractional)
    else
       radius=self                  %radiusEnclosingMassNumerical(mass,massFractional)
    end if
    return
  end function sphericalDecoratorRadiusEnclosingMassNonAnalytic

  double precision function sphericalDecoratorRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target :: self
    double precision                                    , intent(in   )         :: angularMomentumSpecific
    
    radius=self%radiusFromSpecificAngularMomentumNonAnalytic(angularMomentumSpecific)
    return
  end function sphericalDecoratorRadiusFromSpecificAngularMomentum
  
  double precision function sphericalDecoratorRadiusFromSpecificAngularMomentumNonAnalytic(self,angularMomentumSpecific) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for decorator spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target :: self
    double precision                                    , intent(in   )         :: angularMomentumSpecific

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusFromSpecificAngularMomentum         (angularMomentumSpecific)
    else
       radius=self                  %radiusFromSpecificAngularMomentumNumerical(angularMomentumSpecific)
    end if
    return
  end function sphericalDecoratorRadiusFromSpecificAngularMomentumNonAnalytic
  
  double precision function sphericalDecoratorDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for the decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: moment
    double precision                                    , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                             , intent(  out), optional :: isInfinite
 
    densityRadialMoment=self%densityRadialMomentNonAnalytic(moment,radiusMinimum,radiusMaximum,isInfinite)
    return
  end function sphericalDecoratorDensityRadialMoment

  double precision function sphericalDecoratorDensityRadialMomentNonAnalytic(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for the decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout)           :: self
    double precision                                    , intent(in   )           :: moment
    double precision                                    , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                             , intent(  out), optional :: isInfinite

    if (self%useUndecorated()) then
       densityRadialMoment=self%massDistribution_%densityRadialMoment         (moment,radiusMinimum,radiusMaximum,isInfinite)
    else
       densityRadialMoment=self                  %densityRadialMomentNumerical(moment,radiusMinimum,radiusMaximum,isInfinite)
    end if
    return
  end function sphericalDecoratorDensityRadialMomentNonAnalytic

  logical function sphericalDecoratorPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return if the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionSphericalDecorator), intent(inout) :: self

    isAnalytic=self%useUndecorated() .and. self%massDistribution_%potentialIsAnalytic()
    return
  end function sphericalDecoratorPotentialIsAnalytic

  double precision function sphericalDecoratorPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalDecorator), intent(inout), target   :: self
    class(coordinate                        ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType ), intent(  out), optional :: status

    potential=self%potentialNonAnalytic(coordinates,status)
    return
  end function sphericalDecoratorPotential

  double precision function sphericalDecoratorPotentialNonAnalytic(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a decorator spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalDecorator), intent(inout), target   :: self
    class(coordinate                        ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType ), intent(  out), optional :: status

    if (self%useUndecorated()) then
       potential=self%massDistribution_%potential         (coordinates,status)
    else
       potential=self                  %potentialNumerical(coordinates,status)
    end if
    return
  end function sphericalDecoratorPotentialNonAnalytic

  double precision function sphericalDecoratorFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter, wavenumber
    
    fourierTransform=self%fourierTransformNonAnalytic(radiusOuter,wavenumber)
    return
  end function sphericalDecoratorFourierTransform
  
  double precision function sphericalDecoratorFourierTransformNonAnalytic(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter, wavenumber
    
    if (self%useUndecorated()) then
       fourierTransform=self%massDistribution_%fourierTransform         (radiusOuter,wavenumber)
    else
       fourierTransform=self                  %fourierTransformNumerical(radiusOuter,wavenumber)
    end if
    return
  end function sphericalDecoratorFourierTransformNonAnalytic
  
  double precision function sphericalDecoratorRadiusFreefall(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
 
    radius=self%radiusFreefallNonAnalytic(time)
    return
  end function sphericalDecoratorRadiusFreefall
  
  double precision function sphericalDecoratorRadiusFreefallNonAnalytic(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: time

    if (self%useUndecorated()) then
       radius=self%massDistribution_%radiusFreefall         (time)
    else
       radius=self                  %radiusFreefallNumerical(time)
    end if
    return
  end function sphericalDecoratorRadiusFreefallNonAnalytic
  
  double precision function sphericalDecoratorRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
 
    radiusIncreaseRate=self%radiusFreefallIncreaseRateNonAnalytic(time)
    return
  end function sphericalDecoratorRadiusFreefallIncreaseRate

  double precision function sphericalDecoratorRadiusFreefallIncreaseRateNonAnalytic(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in a decorator spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: time

    if (self%useUndecorated()) then
       radiusIncreaseRate=self%massDistribution_%radiusFreefallIncreaseRate         (time)
    else
       radiusIncreaseRate=self                  %radiusFreefallIncreaseRateNumerical(time)
    end if
    return
  end function sphericalDecoratorRadiusFreefallIncreaseRateNonAnalytic

  double precision function sphericalDecoratorEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the total energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout), target :: self
    double precision                                    , intent(in   )         :: radiusOuter
    class           (massDistributionClass             ), intent(inout), target :: massDistributionEmbedding

    energy=self%energyNonAnalytic(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalDecoratorEnergy

  double precision function sphericalDecoratorEnergyNonAnalytic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the total energy within a given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalDecorator), intent(inout) :: self
    double precision                                    , intent(in   ) :: radiusOuter
    class           (massDistributionClass             ), intent(inout) :: massDistributionEmbedding

    if (self%useUndecorated()) then
       energy=self%massDistribution_%energy         (radiusOuter,self%massDistribution_        )
    else
       energy=self                  %energyNumerical(radiusOuter,     massDistributionEmbedding)
    end if
    return
  end function sphericalDecoratorEnergyNonAnalytic
  
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
       energy=self%massDistribution_%energyKinetic         (radiusOuter,self%massDistribution_        )
    else
       energy=self                  %energyKineticNumerical(radiusOuter,     massDistributionEmbedding)
    end if
    return
  end function sphericalDecoratorEnergyKineticNonAnalytic
