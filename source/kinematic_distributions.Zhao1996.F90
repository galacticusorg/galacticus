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

!+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

  !!{
  Implementation of a kinematic distribution class for \cite{zhao_analytical_1996} mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionZhao1996">
   <description>A kinematic distribution class for \cite{zhao_analytical_1996} mass distributions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionZhao1996
     !!{
     A kinematics distribution for \cite{zhao_analytical_1996} distributions.
     !!}
   contains
     procedure :: isCollisional        => zhao1996IsCollisional
     procedure :: velocityDispersion1D => zhao1996VelocityDispersion1D
  end type kinematicsDistributionZhao1996

  interface kinematicsDistributionZhao1996
     !!{
     Constructors for the \refClass{kinematicsDistributionZhao1996} kinematic distribution class.
     !!}
     module procedure zhao1996ConstructorParameters
  end interface kinematicsDistributionZhao1996

contains

  function zhao1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionZhao1996} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(kinematicsDistributionZhao1996)                :: self
    type(inputParameters               ), intent(inout) :: parameters
    
    self=kinematicsDistributionZhao1996()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zhao1996ConstructorParameters
  
  logical function zhao1996IsCollisional(self)
    !!{
    Return false indicating that the zhao1996 kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionZhao1996), intent(inout) :: self
    
    zhao1996IsCollisional=.false.
    return
  end function zhao1996IsCollisional

  double precision function zhao1996VelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an Zhao1996 kinematic distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Dilogarithms                    , only : Dilogarithm
    implicit none
    class           (kinematicsDistributionZhao1996), intent(inout)          :: self
    class           (coordinate                    ), intent(in   )          :: coordinates
    class           (massDistributionClass         ), intent(inout), target  :: massDistribution_        , massDistributionEmbedding
    class           (massDistributionClass         )               , pointer :: massDistribution__
    double precision                                , parameter              :: radiusTiny        =1.0d-3, radiusLarge              =1.0d2
    double precision                                                         :: radius

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating Zhao1996 distribution we have an analytic solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionZhao1996)
          radius            =+coordinates              %rSpherical () &
               &             /massDistributionEmbedding%scaleLength
          select case (massDistributionEmbedding%specialCase%ID)
          case (specialCaseGeneral    %ID)
             ! No analytic solution is available - fall back to the numerical solution.
             velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
             return
          case (specialCaseNFW        %ID)
             if      (radius < radiusTiny ) then
                ! Use series solution for small radii.
                velocityDispersion=+11.0d0*Pi/15.0d0*radius**4                                            &
                     &             +       Pi/ 6.0d0*radius**3*(-101.0d0+12.0d0*Pi**2-12.0d0*log(radius)) &
                     &             + 2.0d0*Pi/ 3.0d0*radius**2*(- 59.0d0+ 6.0d0*Pi**2- 6.0d0*log(radius)) &
                     &             +       Pi       *radius   *(- 23.0d0+ 2.0d0*Pi**2- 2.0d0*log(radius))
             else if (radius > radiusLarge) then
                ! Use series solution for large radii.
                velocityDispersion=+Pi*(-  3.0d0+  4.0d0*log(radius))/(   4.0d0*radius   ) &
                     &             +Pi*(+ 69.0d0+ 20.0d0*log(radius))/(  50.0d0*radius**2) &
                     &             +Pi*(- 97.0d0- 60.0d0*log(radius))/( 300.0d0*radius**3) &
                     &             +Pi*(+284.0d0+420.0d0*log(radius))/(3675.0d0*radius**4)
             else
                velocityDispersion=+2.0d0                                                                                                         &
                     &             *Pi                                                                                                            &
                     &             *(                                                                                                             &
                     &               +radius                                                                                                      &
                     &               *(                                                                                                           &
                     &                 -1.0d0                                                                                                     &
                     &                 +radius                                                                                                    &
                     &                 *(                                                                                                         &
                     &                   -9.0d0                                                                                                   &
                     &                   -7.0d0   *        radius                                                                                 &
                     &                   +Pi   **2*(+1.0d0+radius)**2                                                                             &
                     &                  )                                                                                                         &
                     &                )                                                                                                           &
                     &               +radius**4*log(+1.0d0+1.0d0/radius)                                                                          &
                     &               +log(+1.0d0+radius)+radius*(-(radius*(+1.0d0+2.0d0*radius)*log(radius))                                      &
                     &               +log(+1.0d0+radius)*(-2.0d0-4.0d0*radius*(2.0d0+radius)+3.0d0*radius*(+1.0d0+radius)**2*log(+1.0d0+radius))) &
                     &               +6.0d0*radius**2*(+1.0d0+radius)**2*Dilogarithm(-radius)                                                     &
                     &              )                                                                                                             &
                     &             /radius
             end if
          case (specialCaseCoredNFW   %ID)
             if      (radius < radiusTiny ) then
                ! Use series solution for small radii.
                velocityDispersion=+(119.0d0*Pi-12.0d0*Pi**3)/ 6.0d0           &
                     &             +(119.0d0*Pi-12.0d0*Pi**3)/ 2.0d0*radius    &
                     &             +(353.0d0*Pi-36.0d0*Pi**3)/ 6.0d0*radius**2 &
                     &             +(121.0d0*Pi-12.0d0*Pi**3)/ 6.0d0*radius**3 &
                     &                         - 9.0d0*Pi**4 /20.0d0*radius**4
             else if (radius > radiusLarge) then
                ! Use series solution for large radii.
                velocityDispersion=+Pi*(-   5.0d0+   4.0d0*log(radius))/(    4.0d0*radius   ) &
                     &             +Pi*(+ 177.0d0+  60.0d0*log(radius))/(  100.0d0*radius**2) &
                     &             +Pi*(- 157.0d0-  60.0d0*log(radius))/(  300.0d0*radius**3) &
                     &             +Pi*(+5857.0d0+1260.0d0*log(radius))/(14700.0d0*radius**4)
             else
                ! Use full solution.
                velocityDispersion=+(+1.0d0+radius)**3                                                                               &
                     &             *(                                                                                                &
                     &               +Pi                                                                                             &
                     &               *(                                                                                              &
                     &                 +95.0d0                                                                                       &
                     &                 -12.0d0*Pi**2*(+1.0d0+radius)**4                                                              &
                     &                 + 2.0d0      *        radius    *(130.0d0+9.0d0*radius*(13.0d0+4.0d0*radius))                 &
                     &                )                                                                                              &
                     &               /6.0d0                                                                                          &
                     &               /(+1.0d0+radius)**4                                                                             &
                     &               +(2.0d0*Pi*(2.0d0+9.0d0*radius+6.0d0*radius**2)*log(+1.0d0+radius))/(radius*(+1.0d0+radius)**2) &
                     &               - 6.0d0*Pi*log        (+1.0d0+radius)**2                                                        &
                     &               -12.0d0*Pi*Dilogarithm(      -radius)                                                           &
                     &              )
             end if
          case (specialCaseGamma0_5NFW%ID)
             if      (radius < radiusTiny ) then
                ! Use series solution for small radii.
                velocityDispersion=+16.0d0*Pi/  3.0d0*sqrt(radius      )*(- 11.0d0+  16.0d0*log(2.0d0)) &
                     &             +16.0d0*Pi/ 15.0d0*     radius**1.5d0*(-139.0d0+ 200.0d0*log(2.0d0)) &
                     &             + 2.0d0*Pi/  7.0d0*     radius**2.5d0*(-387.0d0+ 560.0d0*log(2.0d0)) &
                     &             + 4.0d0*Pi/189.0d0*     radius**3.5d0*(-887.0d0+1260.0d0*log(2.0d0))
             else if (radius > radiusLarge) then
                ! Use series solution for large radii.
                velocityDispersion=+Pi*(- 29.0d0+ 24.0d0*log(2.0d0)+ 12.0d0*log(radius))/(  12.0d0*radius   ) &
                     &             +Pi*(+107.0d0+120.0d0*log(2.0d0)+ 60.0d0*log(radius))/( 120.0d0*radius**2) &
                     &             +Pi*(- 11.0d0- 40.0d0*log(2.0d0)- 20.0d0*log(radius))/(  96.0d0*radius**3) &
                     &             +Pi*(+ 57.0d0+280.0d0*log(2.0d0)+140.0d0*log(radius))/(1344.0d0*radius**4)
             else
                ! Use full solution.
                velocityDispersion=+8.0d0                                                                                                                  &
                     &             /9.0d0                                                                                                                  &
                     &             *Pi                                                                                                                     &
                     &             *(                                                                                                                      &
                     &               - 6.0d0*        sqrt(radius   *(+1.0d0+radius))                                                                       &
                     &               -38.0d0*        sqrt(radius**3*(+1.0d0+radius))                                                                       &
                     &               -57.0d0*        sqrt(radius**5*(+1.0d0+radius))                                                                       &
                     &               -24.0d0*        sqrt(radius**7*(+1.0d0+radius))                                                                       &
                     &               +24.0d0*(                                                                                                             &
                     &                        +      sqrt(radius**3*(+1.0d0+radius))                                                                       &
                     &                        +3.0d0*sqrt(radius**5*(+1.0d0+radius))                                                                       &
                     &                        +3.0d0*sqrt(radius**7*(+1.0d0+radius))                                                                       &
                     &                        +      sqrt(radius**9*(+1.0d0+radius))                                                                       &
                     &               )                                                                                                                     &
                     &               * (log(radius) + log(+1.0d0+radius))                                                                                  &
                     &               - 6.0d0*(+1.0d0+radius)**2    *(+1.0d0+2.0d0*radius)       *(-1.0d0+8.0d0*radius*(+1.0d0+radius))*asinh(sqrt(radius)) &
                     &               -24.0d0*        radius **1.5d0*(+1.0d0      +radius)**3.5d0*(log(+radius)-log(+16.0d0*(+1.0d0+radius)))               &
                     &              )                                                                                                                      &
                     &             /        radius                                                                                                         &
                     &             /(+1.0d0+radius)
             end if
          case (specialCaseGamma1_5NFW%ID)
             if      (radius < radiusTiny ) then
                ! Use series solution for small radii.
                velocityDispersion=+8.0d0*Pi/   3.0d0*sqrt(radius       )                                                        &
                     &             -      Pi/3150.0d0*     radius**3.5d0 *(-19861.0d0+60480.0d0*log(2.0d0)-7560.0d0*log(radius)) &
                     &             -      Pi/ 175.0d0*     radius**2.5d0 *(- 8683.0d0+13440.0d0*log(2.0d0)-1680.0d0*log(radius)) &
                     &             -4.0d0*Pi/  75.0d0*     radius**1.5d0 *(-  817.0d0+  960.0d0*log(2.0d0)- 120.0d0*log(radius))
             else if (radius > radiusLarge) then
                ! Use series solution for large radii.
                velocityDispersion=+Pi*(-   7.0d0+   8.0d0*log(2.0d0)+   4.0d0*log(radius))/(    4.0d0*radius   ) &
                     &             +Pi*(+ 147.0d0+ 120.0d0*log(2.0d0)+  60.0d0*log(radius))/(  200.0d0*radius**2) &
                     &             +Pi*(-  79.0d0- 840.0d0*log(2.0d0)- 420.0d0*log(radius))/( 2400.0d0*radius**3) &
                     &             +Pi*(-3589.0d0+7560.0d0*log(2.0d0)+3780.0d0*log(radius))/(33600.0d0*radius**4)
             else
                ! Use full solution.
                velocityDispersion=+8.0d0                                                                                                               &
                     &             /5.0d0                                                                                                               &
                     &             *Pi                                                                                                                  &
                     &             *        radius **1.5d0                                                                                              &
                     &             *(+1.0d0+radius)**1.5d0                                                                                              &
                     &             *(                                                                                                                   &
                     &               +2.0d0*(+1.0d0+2*radius*(-1.0d0+4.0d0*radius+8.0d0*radius**2))*asinh(sqrt(radius))/sqrt(radius**5*(+1.0d0+radius)) &
                     &               +(                                                                                                                 &
                     &                 -2.0d0+radius   *(5.0d0+12.0d0*radius-32.0d0*radius*(+1.0d0+radius)* log(2.0d0 )                          )      &
                     &                 +4.0d0*radius**2*                                   (+1.0d0+radius)*(log(radius)-5.0d0*log(+1.0d0+radius))       &
                     &                )                                                                                                                 &
                     &               /radius**2                                                                                                         &
                     &               /(+1.0d0+radius)                                                                                                   &
                     &              )
             end if
          case default
             velocityDispersion=+0.0d0
             call Error_Report('unknown special case'//{introspection:location})
          end select
          velocityDispersion=+sqrt(                                                &
               &                   +velocityDispersion                             &
               &                   *gravitationalConstant_internal                 &
               &                   *massDistributionEmbedding%densityNormalization &
               &                  )                                                &
               &                   *      massDistributionEmbedding%scaleLength
       class default
          velocityDispersion=0.0d0
          call Error_Report('expecting a Zhao1996 mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our Zhao1996 distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function zhao1996VelocityDispersion1D
