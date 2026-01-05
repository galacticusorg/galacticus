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
  Implementation of a kinematic distribution class for NFW mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionNFW">
   <description>A kinematic distribution class for NFW mass distributions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionNFW
     !!{
     A kinematics distribution for NFW distributions.
     !!}
     logical :: useSeriesApproximation
   contains
     procedure :: isCollisional        => nfwIsCollisional
     procedure :: velocityDispersion1D => nfwVelocityDispersion1D
  end type kinematicsDistributionNFW

  interface kinematicsDistributionNFW
     !!{
     Constructors for the \refClass{kinematicsDistributionNFW} kinematic distribution class.
     !!}
     module procedure nfwConstructorParameters
     module procedure nfwConstructorInternal
  end interface kinematicsDistributionNFW

contains

  function nfwConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionNFW} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (kinematicsDistributionNFW)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    logical                                           :: useSeriesApproximation

    !![
    <inputParameter>
    <name>useSeriesApproximation</name>
    <defaultValue>.false.</defaultValue>
    <description>If true, use a fast series approximation to the velocity dispersion profile in an NFW mass distribution.</description>
    <source>parameters</source>
    </inputParameter>
    !!]
    self=kinematicsDistributionNFW(useSeriesApproximation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nfwConstructorParameters

  function nfwConstructorInternal(useSeriesApproximation) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionNFW} kinematic distribution class.
    !!}
    implicit none
    type   (kinematicsDistributionNFW)                :: self
    logical                           , intent(in   ) :: useSeriesApproximation
    !![
    <constructorAssign variables="useSeriesApproximation"/>
    !!]
    
    return
  end function nfwConstructorInternal
  
  logical function nfwIsCollisional(self)
    !!{
    Return true indicating that the nfw kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionNFW), intent(inout) :: self
    
    nfwIsCollisional=.false.
    return
  end function nfwIsCollisional

  double precision function nfwVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an NFW kinematic distribution.
    !!}
    use :: Dilogarithms                    , only : Dilogarithm
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (kinematicsDistributionNFW), intent(inout)                      :: self
    class           (coordinate               ), intent(in   )                      :: coordinates
    class           (massDistributionClass    ), intent(inout), target              :: massDistribution_                      , massDistributionEmbedding
    class           (massDistributionClass    )               , pointer             :: massDistribution__
    double precision                           , parameter                          :: minimumRadiusForExactSolution   =1.0d-2
    double precision                           , parameter                          :: maximumRadiusForExactSolution   =1.0d+2    
    double precision                           , parameter                          :: nfwNormalizationFactorUnitRadius=-8.5d0+Pi**2-6.0d0*log(2.0d0)+6.0d0*log(2.0d0)**2 ! Precomputed NFW normalization factor for unit radius.
    integer                                    , parameter                          :: maximumExpansionOrder           =7
    double precision                           , dimension(maximumExpansionOrder+1) :: coefficient                           , radiusPower
    double precision                                                                :: logRadius                             , onePlusRadius           , &
         &                                                                             logOnePlusRadius                      , velocityDispersionSquare, &
         &                                                                             radius
    integer                                                                         :: i

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating NFW distribution we have an analytic solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionNFW)
          radius =+coordinates              %rSpherical () &
               &  /massDistributionEmbedding%scaleLength
          if  (self%useSeriesApproximation) then
             if (radius == 0.0d0) then
                velocityDispersionSquare=0.0d0
             else
                if      (radius < 0.33d0) then
                   ! Expand around 0.
                   radiusPower(1)= 1.0d0
                   radiusPower(2)= radius
                   logRadius     = log(radius)
                   coefficient(1)=  0.0d0
                   coefficient(2)=  1.0d0/   4.0d0*(-23.0d0       + 2.0d0*Pi**2- 2.0d0*logRadius)
                   coefficient(3)=                 (-59.0d0/6.0d0 +       Pi**2-       logRadius)
                   coefficient(4)=  1.0d0/  24.0d0*(-101.0d0      +12.0d0*Pi**2-12.0d0*logRadius)
                   coefficient(5)= 11.0d0/  60.0d0
                   coefficient(6)=-13.0d0/ 240.0d0
                   coefficient(7)= 37.0d0/1400.0d0
                   coefficient(8)=-17.0d0/1050.0d0
                else if (radius <  0.68d0) then
                   ! Expand around 1/2.
                   radiusPower(1)= 1.0d0
                   radiusPower(2)= radius-0.5d0
                   coefficient(1)= 9.2256912491493508d-2
                   coefficient(2)= 1.8995942538987498d-2
                   coefficient(3)=-6.1247239215578800d-2
                   coefficient(4)= 9.7544538830827322d-2
                   coefficient(5)=-1.4457663797045428d-1
                   coefficient(6)= 2.1545129876370470d-1
                   coefficient(7)=-3.2824371986452579d-1
                   coefficient(8)= 5.1242111712986012d-1
                else if (radius < 1.35d0) then
                   ! Expand around 1.
                   radiusPower(1)= 1.0d0
                   radiusPower(2)= radius-1.0d0
                   coefficient(1)= 9.3439401238895310d-2
                   coefficient(2)=-6.2683780821546887d-3
                   coefficient(3)=-8.2007484513808621d-3
                   coefficient(4)= 1.0119593363084506d-2
                   coefficient(5)=-9.2481085050239271d-3
                   coefficient(6)= 7.8754354146912774d-3
                   coefficient(7)=-6.5855139302751235d-3
                   coefficient(8)= 5.5035102596088475d-3
                else if (radius < 2.66d0) then
                   ! Expand around 2.
                   radiusPower(1)= 1.0d0
                   radiusPower(2)= radius-2.0d0
                   coefficient(1)= 8.4126434467263518d-2
                   coefficient(2)=-9.8388986218866523d-3
                   coefficient(3)= 6.1288152708705594d-4
                   coefficient(4)= 4.3464937545102683d-4
                   coefficient(5)=-3.4479664620159904d-4
                   coefficient(6)= 1.8815165134120623d-4
                   coefficient(7)=-9.2066324234421410d-5
                   coefficient(8)= 4.3068151103206337d-5
                else
                   ! Expand around infinity.
                   radiusPower(1)= 1.0d0
                   radiusPower(2)= 1.0d0/radius
                   logRadius     = log(radius)
                   coefficient(1)=     0.0d0
                   coefficient(2)=(-   3.0d0+   4.0d0*logRadius)/    16.0d0
                   coefficient(3)=(   69.0d0+  20.0d0*logRadius)/   200.0d0
                   coefficient(4)=(-  97.0d0-  60.0d0*logRadius)/  1200.0d0
                   coefficient(5)=(   71.0d0+ 105.0d0*logRadius)/  3675.0d0
                   coefficient(6)=(-   1.0d0-  56.0d0*logRadius)/  3136.0d0
                   coefficient(7)=(-1271.0d0+2520.0d0*logRadius)/211680.0d0
                   coefficient(8)=(  341.0d0- 360.0d0*logRadius)/ 43200.0d0
                end if
                do i=3,maximumExpansionOrder+1
                   radiusPower(i)=radiusPower(i-1)*radiusPower(2)
                end do
                velocityDispersionSquare=sum(coefficient*radiusPower)
             end if
          else
             if (radius == 1.0d0) then
                velocityDispersionSquare=nfwNormalizationFactorUnitRadius
             else if (radius >= maximumRadiusForExactSolution) then
                logRadius                      = log(radius)
                velocityDispersionSquare=+(-   3.0d0+   4.0d0*logRadius)/(    16.0d0*radius   ) &
                     &                   +(   69.0d0+  20.0d0*logRadius)/(   200.0d0*radius**2) &
                     &                   +(-  97.0d0-  60.0d0*logRadius)/(  1200.0d0*radius**3) &
                     &                   +(   71.0d0+ 105.0d0*logRadius)/(  3675.0d0*radius**4) &
                     &                   +(-   1.0d0-  56.0d0*logRadius)/(  3136.0d0*radius**5) &
                     &                   +(-1271.0d0+2520.0d0*logRadius)/(211680.0d0*radius**6)
             else if (radius >= minimumRadiusForExactSolution) then
                onePlusRadius                 =      1.0d0+radius
                logRadius                     = log(       radius)
                logOnePlusRadius              = log(onePlusRadius)
                velocityDispersionSquare=+0.5d0                      &
                     &                   *       radius              &
                     &                   *onePlusRadius**2           &
                     &                   *(                          &
                     &                     +Pi**2                    &
                     &                     -logRadius                &
                     &                     -1.0d0/       radius      &
                     &                     -1.0d0/onePlusRadius**2   &
                     &                     -6.0d0/onePlusRadius      &
                     &                     +(                        &
                     &                       +1.0d0+ 1.0d0/radius**2 &
                     &                             - 4.0d0/radius    &
                     &                       -2.0d0/onePlusRadius    &
                     &                      )                        &
                     &                     *logOnePlusRadius         &
                     &                     +3.0d0                    &
                     &                     *logOnePlusRadius**2      &
                     &                     +6.0d0                    &
                     &                     *Dilogarithm(-radius)     &
                     &                    )
             else if (radius > 0.0d0) then
                logRadius                     = log(radius)
                velocityDispersionSquare=+ 1.0d0/   4.0d0*(-23.0d0       + 2.0d0*Pi**2- 2.0d0*logRadius)*radius    &
                     &                   +                (-59.0d0/6.0d0 +       Pi**2-       logRadius)*radius**2 &
                     &                   + 1.0d0/  24.0d0*(-101.0d0      +12.0d0*Pi**2-12.0d0*logRadius)*radius**3 &
                     &                   +11.0d0/  60.0d0                                               *radius**4 &
                     &                   -13.0d0/ 240.0d0                                               *radius**5 &
                     &                   +37.0d0/1400.0d0                                               *radius**6
             else
                velocityDispersionSquare=0.0d0
             end if
          end if
          velocityDispersion=+sqrt(                                                &
               &                   +4.0d0                                          &
               &                   *Pi                                             &
               &                   *velocityDispersionSquare                       &
               &                   *gravitationalConstant_internal                 &
               &                   *massDistributionEmbedding%densityNormalization &
               &                  )                                                &
               &             *      massDistributionEmbedding%scaleLength
       class default
          velocityDispersion=0.0d0
          call Error_Report('expecting an NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our NFW distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function nfwVelocityDispersion1D
