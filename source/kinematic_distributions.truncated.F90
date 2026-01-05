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
  Implementation of a kinematic distribution class for truncated mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionTruncated">
   <description>An truncated kinematic distribution class masses.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionTruncated
     !!{
     A kinematic distribution for truncated mass distributions.
     !!}
     double precision :: velocityDispersionDecoratedTruncateMinimum, velocityDispersionTruncateMinimum, &
          &              densityTruncateMinimum
     logical          :: velocityDispersionTruncateMinimumComputed
   contains
     procedure :: isCollisional        => truncatedIsCollisional
     procedure :: velocityDispersion1D => truncatedVelocityDispersion1D
  end type kinematicsDistributionTruncated

  interface kinematicsDistributionTruncated
     !!{
     Constructors for the \refClass{kinematicsDistributionTruncated} kinematic distribution class.
     !!}
     module procedure truncatedConstructorParameters
     module procedure truncatedConstructorInternal
  end interface kinematicsDistributionTruncated

contains

  function truncatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionTruncated} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(kinematicsDistributionTruncated)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=kinematicsDistributionTruncated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function truncatedConstructorParameters

  function truncatedConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionTruncated} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionTruncated) :: self

    self%velocityDispersionTruncateMinimumComputed=.false.
    return
  end function truncatedConstructorInternal
  
  logical function truncatedIsCollisional(self)
    !!{
    Return false indicating that the truncated kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionTruncated), intent(inout) :: self
    
    truncatedIsCollisional=.false.
    return
  end function truncatedIsCollisional

  double precision function truncatedVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an truncated kinematic distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (kinematicsDistributionTruncated), intent(inout)          :: self
    class           (coordinate                     ), intent(in   )          :: coordinates
    class           (massDistributionClass          ), intent(inout), target  :: massDistribution_         , massDistributionEmbedding
    class           (massDistributionClass          )               , pointer :: massDistribution__
    class           (kinematicsDistributionClass    )               , pointer :: kinematicsDistribution_
    type            (coordinateSpherical            )                         :: coordinatesTruncateMinimum
    logical                                                                   :: analytic
    double precision                                                          :: density                   , velocityDispersionDecorated

    analytic=.false.
    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating truncated distribution we have a piecewise solution for the velocity dispersion in some cases.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalTruncated)
          if (coordinates%rSpherical() < massDistributionEmbedding%radiusTruncateMinimum) then
             ! Use the decorated mass distribution solution, adjusted for the outer truncation shell.
             analytic                                           =   .true.
             kinematicsDistribution_                            =>  massDistributionEmbedding%massDistribution_%kinematicsDistribution       (                                                                                                                  )
             density                                            =   massDistributionEmbedding                  %density                      (coordinates                                                                                                       )
             velocityDispersionDecorated                        =   kinematicsDistribution_                    %velocityDispersion1D         (coordinates               ,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
             if (.not.self%velocityDispersionTruncateMinimumComputed) then
                coordinatesTruncateMinimum                      =  [massDistributionEmbedding                  %radiusTruncateMinimum,0.0d0,0.0d0]
                self%densityTruncateMinimum                     =   massDistributionEmbedding                  %density                      (coordinatesTruncateMinimum                                                                                        )
                self%velocityDispersionDecoratedTruncateMinimum =   kinematicsDistribution_                    %velocityDispersion1D         (coordinatesTruncateMinimum,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
                self%velocityDispersionTruncateMinimum          =   self                                       %velocityDispersion1DNumerical(coordinatesTruncateMinimum,massDistributionEmbedding                  ,massDistributionEmbedding                  )
                self%velocityDispersionTruncateMinimumComputed  =   .true.
             end if
             !![
             <objectDestructor name="kinematicsDistribution_"/>
             !!]
          end if
       class is (massDistributionSphericalTruncatedExponential)
          if (coordinates%rSpherical() < massDistributionEmbedding%radiusTruncateMinimum) then
             ! Use the decorated mass distribution solution, adjusted for the outer truncation shell.
             analytic                                           =   .true.
             kinematicsDistribution_                            =>  massDistributionEmbedding%massDistribution_%kinematicsDistribution       (                                                                                                                  )
             density                                            =   massDistributionEmbedding                  %density                      (coordinates                                                                                                       )
             velocityDispersionDecorated                        =   kinematicsDistribution_                    %velocityDispersion1D         (coordinates               ,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
             if (.not.self%velocityDispersionTruncateMinimumComputed) then                                                                                                                                                                                      
                coordinatesTruncateMinimum                      =  [massDistributionEmbedding                  %radiusTruncateMinimum,0.0d0,0.0d0]                                                                                                              
                self%densityTruncateMinimum                     =   massDistributionEmbedding                  %density                      (coordinatesTruncateMinimum                                                                                        )
                self%velocityDispersionDecoratedTruncateMinimum =   kinematicsDistribution_                    %velocityDispersion1D         (coordinatesTruncateMinimum,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
                self%velocityDispersionTruncateMinimum          =   self                                       %velocityDispersion1DNumerical(coordinatesTruncateMinimum,massDistributionEmbedding                  ,massDistributionEmbedding                  )
                self%velocityDispersionTruncateMinimumComputed  =   .true.
             end if
             !![
             <objectDestructor name="kinematicsDistribution_"/>
             !!]
          end if
       class default
          velocityDispersion=+0.0d0
          call Error_Report('mass distribution must be of the `massDistributionSphericalTruncated` or `massDistributionSphericalTruncatedExponential` class'//{introspection:location})
       end select
    end if
    ! Use a numerical solution if no analytic solution was available.
    if (analytic) then
       velocityDispersion=+sqrt(                                                      &
            &                   +       velocityDispersionDecorated               **2 &
            &                   +self%densityTruncateMinimum                          &
            &                   /       density                                       &
            &                   *(                                                    &
            &                     +self%velocityDispersionTruncateMinimum         **2 &
            &                     -self%velocityDispersionDecoratedTruncateMinimum**2 &
            &                    )                                                    &
            &                  )
    else
       velocityDispersion=+self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function truncatedVelocityDispersion1D
  
