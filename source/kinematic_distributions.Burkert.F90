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
  Implementation of a kinematic distribution class for Burkert mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionBurkert">
   <description>A kinematic distribution class for Burkert mass distributions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionBurkert
     !!{
     A kinematics distribution for Burkert distributions.
     !!}
   contains
     procedure :: isCollisional        => burkertIsCollisional
     procedure :: velocityDispersion1D => burkertVelocityDispersion1D
  end type kinematicsDistributionBurkert

  interface kinematicsDistributionBurkert
     !!{
     Constructors for the \refClass{kinematicsDistributionBurkert} kinematic distribution class.
     !!}
     module procedure burkertConstructorParameters
  end interface kinematicsDistributionBurkert

contains

  function burkertConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionBurkert} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (kinematicsDistributionBurkert)                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    
    self=kinematicsDistributionBurkert()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function burkertConstructorParameters
  
  logical function burkertIsCollisional(self)
    !!{
    Return true indicating that the burkert kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionBurkert), intent(inout) :: self
    
    burkertIsCollisional=.false.
    return
  end function burkertIsCollisional

  double precision function burkertVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an Burkert kinematic distribution.
    !!}
    use :: Dilogarithms                    , only : Dilogarithm
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (kinematicsDistributionBurkert), intent(inout)          :: self
    class           (coordinate                   ), intent(in   )          :: coordinates
    class           (massDistributionClass        ), intent(inout), target  :: massDistribution_ , massDistributionEmbedding
    class           (massDistributionClass        )               , pointer :: massDistribution__
    double precision                                                        :: radius

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating Burkert distribution we have an analytic solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionBurkert)
          radius            =+coordinates              %rSpherical ()                                                                                                                                               &
               &             /massDistributionEmbedding%scaleLength
          velocityDispersion=real(                                                                                                                                                                                  &
               &                  +(                                                                                                                                                                                &
               &                    +sqrt(                                                                                                                                                                          &
               &                          +Pi                                                                                                                                                                       &
               &                          /3.0d0                                                                                                                                                                    &
               &                         )                                                                                                                                                                          &
               &                    *sqrt(                                                                                                                                                                          &
               &                          +(                                                                                                                                                                        &
               &                            +(+1.0d0+radius   )                                                                                                                                                     &
               &                            *(+1.0d0+radius**2)                                                                                                                                                     &
               &                            *(                                                                                                                                                                      &
               &                                      +48.0d0         *Pi   *radius                                                                                                                                 &
               &                                      -11.0d0         *Pi**2*radius                                                                                                                                 &
               &                                      -96.0d0                      *atan(radius)                                                                                                                    &
               &                                      -96.0d0               *radius*atan(radius)                                                                                                                    &
               &                                      +12.0d0               *radius*log(             +2.0d0         )**2                                                                                            &
               &                                      - 4.0d0         *Pi   *radius*log(             +8.0d0         )                                                                                               &
               &                              +dcmplx(+ 6.0d0,- 6.0d0)      *radius*log(dcmplx(0.0d0,-1.0d0)- radius)**2                                                                                            &
               &                              +dcmplx(+ 6.0d0,+ 6.0d0)      *radius*log(dcmplx(0.0d0,+1.0d0)- radius)**2                                                                                            &
               &                              +dcmplx(+12.0d0,-12.0d0)      *radius*log(dcmplx(0.0d0,-1.0d0)- radius)   *log((+1.0d0+dcmplx(0.0d0,1.0d0)*radius)/2.0d0)                                             &
               &                              -dcmplx(+24.0d0,+24.0d0)      *radius*atan(radius)*log(dcmplx(0.0d0,-2.0d0)/(dcmplx(0.0d0,-1.0d0)+radius))                                                            &
               &                              -dcmplx(+24.0d0,-24.0d0)      *radius*atan(radius)*log(dcmplx(0.0d0,+2.0d0)/(dcmplx(0.0d0,+1.0d0)+radius))                                                            &
               &                              +dcmplx(+12.0d0,+12.0d0)      *radius*log(                       dcmplx(+0.0d0,+1.0d0)-radius)           *log(dcmplx(+0.0d0,-0.5d0)*(dcmplx(+0.0d0,+1.0d0)+radius  )) &
               &                              -dcmplx(+ 0.0d0,+24.0d0)      *radius*log(              +1.0d0 + dcmplx(+0.0d0,+1.0d0)*radius)           *log(dcmplx(+0.5d0,-0.5d0)*(              +1.0d0 +radius  )) &
               &                              +dcmplx(+ 0.0d0,+24.0d0)      *radius*log(              +1.0d0 - dcmplx(+0.0d0,+1.0d0)*radius)           *log(dcmplx(+0.5d0,+0.5d0)*(              +1.0d0 +radius  )) &
               &                                      +96.0d0                      *log(                                     +1.0d0 +radius )                                                                       &
               &                                      +96.0d0               *radius*log(                                     +1.0d0 +radius )                                                                       &
               &                              -dcmplx(+ 0.0d0,+24.0d0)      *radius*log(dcmplx(-0.5d0,+0.5d0)*(dcmplx(+0.0d0,-1.0d0)+radius))          *log(                                     +1.0d0+radius    ) &
               &                              +dcmplx(+ 0.0d0,+24.0d0)      *radius*log(dcmplx(-0.5d0,-0.5d0)*(dcmplx(+0.0d0,+1.0d0)+radius))          *log(                                     +1.0d0+radius    ) &
               &                                      -24.0d0               *radius*log(                                     +1.0d0 +radius )**2+48.0d0*log(                                     +1.0d0+radius**2 ) &
               &                                              -48.0d0       *radius                                                                    *log(                                     +1.0d0+radius**2 ) &
               &                              -dcmplx(+12.0d0,-12.0d0)      *radius*log(                       dcmplx(+0.0d0,-1.0d0)-radius )   *       log(                                     +1.0d0+radius**2 ) &
               &                              -dcmplx(+12.0d0,+12.0d0)      *radius*log(                       dcmplx(+0.0d0,+1.0d0)-radius )   *       log(                                     +1.0d0+radius**2 ) &
               &                                      -24.0d0               *radius*log(                                     +1.0d0 +radius )   *       log(                                     +1.0d0+radius**2 ) &
               &                              +dcmplx(+12.0d0,+12.0d0)*radius*Dilogarithm(        +0.5d0                + dcmplx(0.0d0,+0.5d0)*        radius    )                                                  &
               &                                      -96.0d0         *radius*Dilogarithm(                                                            -radius    )                                                  &
               &                              -dcmplx(+ 0.0d0,+48.0d0)*radius*Dilogarithm( dcmplx(+0.0d0,-1.0d0)                              *        radius    )                                                  &
               &                              +dcmplx(+ 0.0d0,+48.0d0)*radius*Dilogarithm( dcmplx(+0.0d0,+1.0d0)                              *        radius    )                                                  &
               &                                     -24.0d0          *radius*Dilogarithm(                                                            -radius**2 )                                                  &
               &                              -dcmplx(+ 0.0d0,+24.0d0)*radius*Dilogarithm( dcmplx(-0.5d0,+0.5d0)        *(dcmplx(0.0d0,-1.0d0)        +radius   ))                                                  &
               &                              +dcmplx(+12.0d0,+12.0d0)*radius*Dilogarithm((dcmplx(+0.0d0,-1.0d0)+radius)/(dcmplx(0.0d0,+1.0d0)        +radius   ))                                                  &
               &                              +dcmplx(+ 0.0d0,+24.0d0)*radius*Dilogarithm( dcmplx(-0.5d0,-0.5d0)        *(dcmplx(0.0d0,+1.0d0)        +radius   ))                                                  &
               &                              +dcmplx(+12.0d0,-12.0d0)*radius*Dilogarithm( dcmplx(+0.0d0,-0.5d0)        *(dcmplx(0.0d0,+1.0d0)        +radius   ))                                                  &
               &                              +dcmplx(+12.0d0,-12.0d0)*radius*Dilogarithm((dcmplx(+0.0d0,+1.0d0)+radius)/(dcmplx(0.0d0,-1.0d0)        +radius   ))                                                  &
               &                              -dcmplx(+ 0.0d0,+24.0d0)*radius*Dilogarithm( dcmplx(+0.5d0,-0.5d0)                              *(+1.0d0+radius   ))                                                  &
               &                              +dcmplx(+ 0.0d0,+24.0d0)*radius*Dilogarithm( dcmplx(+0.5d0,+0.5d0)                              *(+1.0d0+radius   ))                                                  &
               &                             )                                                                                                                                                                      &
               &                           )                                                                                                                                                                        &
               &                          /radius                                                                                                                                                                   &
               &                         )                                                                                                                                                                          &
               &                       )                                                                                                                                                                            &
               &                      /4.0d0                                                                                                                                                                        &
               &                     )                                                                                                                                                                              &
               &                    *sqrt(                                                                                                                                                                          &
               &                          +gravitationalConstant_internal                                                                                                                                           &
               &                          *massDistributionEmbedding%densityNormalization                                                                                                                           &
               &                         )                                                                                                                                                                          &
               &                    *      massDistributionEmbedding%scaleLength
          class default
             velocityDispersion=0.0d0
             call Error_Report('expecting a Burkert mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
          end select
    else
       ! Our Burkert distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function burkertVelocityDispersion1D
