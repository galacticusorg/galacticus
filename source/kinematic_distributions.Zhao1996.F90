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
     Constructors for the {\normalfont \ttfamily zhao1996} kinematic distribution class.
     !!}
     module procedure zhao1996ConstructorParameters
     module procedure zhao1996ConstructorInternal
  end interface kinematicsDistributionZhao1996

contains

  function zhao1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily zhao1996} kinematic distribution class which builds the object from a parameter
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

  function zhao1996ConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily zhao1996} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionZhao1996) :: self
    
    return
  end function zhao1996ConstructorInternal
  
  logical function zhao1996IsCollisional(self)
    !!{
    Return true indicating that the zhao1996 kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionZhao1996), intent(inout) :: self
    
    zhao1996IsCollisional=.false.
    return
  end function zhao1996IsCollisional

  double precision function zhao1996VelocityDispersion1D(self,coordinates,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an Zhao1996 kinematic distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Dilogarithms                    , only : Dilogarithm
    implicit none
    class           (kinematicsDistributionZhao1996), intent(inout) :: self
    class           (coordinate                    ), intent(in   ) :: coordinates
    class           (massDistributionClass         ), intent(inout) :: massDistributionEmbedding
    double precision                                                :: radius

    select type (massDistributionEmbedding)
    class is (massDistributionZhao1996)
       radius            =+coordinates              %rSpherical () &
            &             /massDistributionEmbedding%scaleLength
       select case (massDistributionEmbedding%specialCase%ID)
       case (specialCaseGeneral    %ID)
          ! No analytic solution is available - fall back to the numerical solution.
          velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistributionEmbedding)
          return
       case (specialCaseNFW        %ID)
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
               &               +radius**4*log(+1.0d0+1/radius)                                                                              &
               &               +log(+1.0d0+radius)+radius*(-(radius*(+1.0d0+2*radius)*log(radius))                                          &
               &               +log(+1.0d0+radius)*(-2.0d0-4.0d0*radius*(2.0d0+radius)+3.0d0*radius*(+1.0d0+radius)**2*log(+1.0d0+radius))) &
               &               +6.0d0*radius**2*(+1.0d0+radius)**2*Dilogarithm(-radius)                                                     &
               &              )                                                                                                             &
               &             /radius
       case (specialCaseCoredNFW   %ID)
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
       case (specialCaseGamma0_5NFW%ID)
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
       case (specialCaseGamma1_5NFW%ID)
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
       case default
          velocityDispersion=+0.0d0
          call Error_Report('unknown special case'//{introspection:location})
       end select
       velocityDispersion=+sqrt(                                                &
            &                   +velocityDispersion                             &
            &                   *gravitationalConstantGalacticus                &
            &                   *massDistributionEmbedding%densityNormalization &
            &                  )                                                &
            &                   *      massDistributionEmbedding%scaleLength
    class default
       velocityDispersion=0.0d0
       call Error_Report('expecting a Zhao1996 mass distribution'//{introspection:location})
    end select
    return
  end function zhao1996VelocityDispersion1D
