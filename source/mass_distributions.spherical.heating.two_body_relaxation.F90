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
  Implements a mass distribution heating class that computes heating due to two-body relaxation.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingTwoBodyRelaxation">
    <description>A mass distribution heating class that computes heating due to two-body relaxation.</description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingTwoBodyRelaxation
     !!{
     Implementation of a mass distribution heating class that computes heating due to two-body relaxation.
     !!}
     private
     double precision :: massParticle, lengthSoftening, &
          &              timeRelaxing, efficiency
   contains
     procedure :: specificEnergy                 => twoBodyRelaxationSpecificEnergy
     procedure :: specificEnergyGradient         => twoBodyRelaxationSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => twoBodyRelaxationSpecificEnergyIsEverywhereZero
  end type massDistributionHeatingTwoBodyRelaxation

  interface massDistributionHeatingTwoBodyRelaxation
     !!{
     Constructors for the \refClass{massDistributionHeatingTwoBodyRelaxation} mass distribution class.
     !!}
     module procedure twoBodyRelaxationConstructorParameters
     module procedure twoBodyRelaxationConstructorInternal
  end interface massDistributionHeatingTwoBodyRelaxation

contains

  function twoBodyRelaxationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingTwoBodyRelaxation} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionHeatingTwoBodyRelaxation)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    double precision                                                          :: massParticle, lengthSoftening, &
           &                                                                     timeRelaxing, efficiency

    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>The particle mass to use for two-body relaxation calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthSoftening</name>
      <source>parameters</source>
      <description>The softening length to use for two-body relaxation calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>timeRelaxing</name>
      <source>parameters</source>
      <description>The time for which the system has been relaxing.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiency</name>
      <source>parameters</source>
      <description>The fractional efficiency of two-body relaxation heating.</description>
    </inputParameter>
    !!]
    self=massDistributionHeatingTwoBodyRelaxation(massParticle,lengthSoftening,timeRelaxing,efficiency)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function twoBodyRelaxationConstructorParameters
  
  function twoBodyRelaxationConstructorInternal(massParticle,lengthSoftening,timeRelaxing,efficiency) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingTwoBodyRelaxation} dark matter profile heating class.
    !!}
    implicit none
    type             (massDistributionHeatingTwoBodyRelaxation)                :: self
    double precision                                           , intent(in   ) :: massParticle, lengthSoftening, &
         &                                                                        timeRelaxing, efficiency
    !![
    <constructorAssign variables="massParticle, lengthSoftening, timeRelaxing, efficiency"/>
    !!]
 
    return
  end function twoBodyRelaxationConstructorInternal

  double precision function twoBodyRelaxationSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}. The assumption here is that the mean
    fractional change in energy for a particle per crossing time is $8 \log \Lambda / N$ where $N$ is the number of particles
    within radius $r=${\normalfont \ttfamily radius}. The crossing time is approximated by $r/V(r)$ where $V(r)$ is the
    circular velocity at $r$. The Coulomb logarithm is given by $\log\Lambda=\hbox{max}(\epsilon,b_{90})$ where $\epsilon$ is
    the softening length, $b_{90}=2\mathrm{G}m_\mathrm{p}/V^2(r)$, and $m_\mathrm{p}$ is the particle mass. Finally, the
    specific energy is assumed to be $\sigma^2(r)/2\approx V^2(r)/4$.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, megaParsec, gravitationalConstant_internal
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (massDistributionHeatingTwoBodyRelaxation), intent(inout) :: self
    double precision                                          , intent(in   ) :: radius
    class           (massDistributionClass                   ), intent(inout) :: massDistribution_
    double precision                                                          :: particleCount    , velocity               , &
         &                                                                       logarithmCoulomb , impactParameterCritical
    
    if (self%timeRelaxing > 0.0d0) then
       velocity               =+massDistribution_                 %rotationCurve     (radius)
       impactParameterCritical=+2.0d0                                                                &
            &                  *gravitationalConstant_internal                                       &
            &                  *self                              %massParticle                      &
            &                  /velocity                                                        **2
       logarithmCoulomb       =+0.5d0                                                                &
            &                  *log(                                                                 &
            &                       +1.0d0                                                           &
            &                       +(                                                               &
            &                         +radius                                                        &
            &                         /max(                                                          &
            &                              self                   %lengthSoftening                 , &
            &                              impactParameterCritical                                   &
            &                             )                                                          &
            &                        )                                                          **2  &
            &                      )
       particleCount          =+massDistribution_                 %massEnclosedBySphere(radius)      &
            &                  /self                              %massParticle
       energySpecific         =+2.0d0                                                                &
            &                  *self                              %efficiency                        &
            &                  *logarithmCoulomb                                                     &
            &                  *self                              %timeRelaxing                      &
            &                  *velocity                                                        **3  &
            &                  /radius                                                               &
            &                  /particleCount                                                        &
            &                  *kilo                                                                 &
            &                  *gigaYear                                                             &
            &                  /megaParsec
    else
       energySpecific         =+0.0d0
    end if
    return
  end function twoBodyRelaxationSpecificEnergy

  double precision function twoBodyRelaxationSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingTwoBodyRelaxation), intent(inout) :: self
    double precision                                          , intent(in   ) :: radius
    class           (massDistributionClass                   ), intent(inout) :: massDistribution_
    double precision                                                          :: particleCount    , velocity               , &
         &                                                                       logarithmCoulomb , impactParameterCritical, &
         &                                                                       gradientCoulomb
    type            (coordinateSpherical                     )                :: coordinates
    
    if (self%timeRelaxing > 0.0d0) then
       coordinates             =[radius,0.0d0,0.0d0]
       velocity                =+massDistribution_                 %rotationCurve       (radius)
       impactParameterCritical =+2.0d0                                                               &
            &                   *gravitationalConstant_internal                                      &
            &                   *self                              %massParticle                     &
            &                   /velocity                                                      **2
       logarithmCoulomb        =+0.5d0                                                               &
            &                   *log(                                                                &
            &                        +1.0d0                                                          &
            &                        +(                                                              &
            &                          +radius                                                       &
            &                          /max(                                                         &
            &                               self                   %lengthSoftening                , &
            &                               impactParameterCritical                                  &
            &                              )                                                         &
            &                         )                                                         **2  &
            &                       )
       particleCount           =+massDistribution_                 %massEnclosedBySphere(radius)     &
            &                   /self                              %massParticle
       if (self%lengthSoftening > impactParameterCritical) then
          gradientCoulomb=+radius                                                            &
               &          /self                             %lengthSoftening
       else
          gradientCoulomb=+2.0d0                                                             &
               &          *radius                                                            &
               &          /impactParameterCritical                                           &
               &          *(                                                                 &
               &            -1.0d0                                                           &
               &            +8.0d0                                                           &
               &            *Pi                                                              &
               &            *gravitationalConstant_internal                                  &
               &            *radius                                                      **2 &
               &            *massDistribution_              %density        (coordinates)    &
               &            /velocity                                                    **2 &
               &           )
       end if
       energySpecificGradient=+self                             %specificEnergy(radius     ,massDistribution_)    &
            &                 /                                                 radius                            &
            &                 *(                                                                                  &
            &                   -2.5d0                                                                            &
            &                   +6.0d0                                                                            &
            &                   *Pi                                                                               &
            &                   *gravitationalConstant_internal                                                   &
            &                   *massDistribution_              %density       (coordinates                  )    &
            &                   *radius                                                                       **2 &
            &                   /velocity                                                                     **2 &
            &                   -                                               radius                            &
            &                   *                gradientCoulomb                                                  &
            &                   /               logarithmCoulomb                                                  &
            &                   *sqrt(exp(2.0d0*logarithmCoulomb)-1.0d0)                                          &
            &                   /     exp(2.0d0*logarithmCoulomb)                                                 &
            &                  )
    else
       energySpecificGradient=+0.0d0
    end if
    return
  end function twoBodyRelaxationSpecificEnergyGradient

  logical function twoBodyRelaxationSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingTwoBodyRelaxation), intent(inout) :: self

    energySpecificIsEverywhereZero=self%timeRelaxing <= 0.0d0
    return
  end function twoBodyRelaxationSpecificEnergyIsEverywhereZero
