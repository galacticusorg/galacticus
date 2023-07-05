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
     Constructors for the {\normalfont \ttfamily burkert} kinematic distribution class.
     !!}
     module procedure burkertConstructorParameters
     module procedure burkertConstructorInternal
  end interface kinematicsDistributionBurkert

contains

  function burkertConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily burkert} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (kinematicsDistributionBurkert)                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    logical                                               :: useSeriesApproximation

    self=kinematicsDistributionBurkert()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function burkertConstructorParameters

  function burkertConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily burkert} kinematic distribution class.
    !!}
    implicit none
    type   (kinematicsDistributionBurkert)                :: self
    
    return
  end function burkertConstructorInternal
  
  logical function burkertIsCollisional(self)
    !!{
    Return true indicating that the burkert kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionBurkert), intent(inout) :: self
    
    burkertIsCollisional=.false.
    return
  end function burkertIsCollisional

  double precision function burkertVelocityDispersion1D(self,coordinates,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an Burkert kinematic distribution.
    !!}
    implicit none
    class           (kinematicsDistributionBurkert), intent(inout)                      :: self
    class           (coordinate               ), intent(in   )                      :: coordinates
    class           (massDistributionClass    ), intent(inout)                      :: massDistributionEmbedding

!! AJB TODO
    return
  end function burkertVelocityDispersion1D
