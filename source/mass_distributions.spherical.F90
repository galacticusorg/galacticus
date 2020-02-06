!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of an abstract mass distribution class for spherically symmetric distributions.

  !# <massDistribution name="massDistributionSpherical" abstract="yes">
  !#  <description>An abstract mass distribution class for spherically symmetric distributions.</description>
  !# </massDistribution>
  type, extends(massDistributionClass), abstract :: massDistributionSpherical
     !% Implementation of an abstract mass distribution class for spherically symmetric distributions.
     private
   contains
     !@ <objectMethods>
     !@   <object>massDistributionSpherical</object>
     !@   <objectMethod>
     !@     <method>radiusHalfMass</method>
     !@     <description>Returns the radius enclosing half of the mass of the mass distribution.</description>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusEnclosingMass</method>
     !@     <description>Returns the radius enclosing the given mass.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ mass\argin</arguments>
     !@   </objectMethod>
      !@ </objectMethods>
     procedure :: symmetry             => sphericalSymmetry
     procedure :: massEnclosedBySphere => sphericalMassEnclosedBySphere
     procedure :: radiusHalfMass       => sphericalRadiusHalfMass
     procedure :: radiusEnclosingMass  => sphericalRadiusEnclosingMass
     procedure :: positionSample       => sphericalPositionSample
  end type massDistributionSpherical

  ! Module scope variables used in integration and root finding.
  class           (massDistributionSpherical), pointer :: sphericalActive
  double precision                                     :: sphericalMassTarget
  !$omp threadprivate(sphericalActive,sphericalMassTarget)

contains

  integer function sphericalSymmetry(self)
    !% Returns symmetry label for mass dsitributions with spherical symmetry.
    implicit none
    class(massDistributionSpherical), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalSymmetry=massDistributionSymmetrySpherical
    return
  end function sphericalSymmetry

  double precision function sphericalMassEnclosedBySphere(self,radius)
    !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    !% distributions using numerical integration.
    use :: FGSL                    , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : Integrate    , Integrate_Done
    implicit none
    class           (massDistributionSpherical ), intent(inout), target :: self
    double precision                            , intent(in   )         :: radius
    type            (fgsl_function             )                        :: integrandFunction
    type            (fgsl_integration_workspace)                        :: integrationWorkspace
    logical                                                             :: integrationReset

    sphericalActive  => self
    integrationReset =  .true.
    sphericalMassEnclosedBySphere=+4.0d0                                                               &
         &                        *Pi                                                                  &
         &                        *Integrate(                                                          &
         &                                                     0.0d0                                 , &
         &                                                     radius                                , &
         &                                                     sphericalMassEnclosedBySphereIntegrand, &
         &                                                     integrandFunction                     , &
         &                                                     integrationWorkspace                  , &
         &                                   reset            =integrationReset                      , &
         &                                   toleranceAbsolute=0.0d+0                                , &
         &                                   toleranceRelative=1.0d-6                                  &
         &                       )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function sphericalMassEnclosedBySphere

  double precision function sphericalMassEnclosedBySphereIntegrand(radius)
    !% Enclosed mass integrand for spherical mass distributions.
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: position

    position                              =[radius,0.0d0,0.0d0]
    sphericalMassEnclosedBySphereIntegrand=+radius                           **2 &
         &                                 *sphericalActive%density(position)
    return
  end function sphericalMassEnclosedBySphereIntegrand

  double precision function sphericalRadiusEnclosingMass(self,mass)
    !% Computes the radius enclosing a given mass in a spherically symmetric mass distribution using numerical root finding.
    use :: FGSL       , only : FGSL_Root_fSolver_Brent
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: mass
    type            (rootFinder               ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                           , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6

    if (.not.finder%isInitialized()) then
       call finder%rootFunction(                                                             &
            &                                                 sphericalMassRoot              &
            &                  )
       call finder%tolerance   (                                                             &
            &                                                 toleranceAbsolute            , &
            &                                                 toleranceRelative              &
            &                  )
       call finder%type        (                                                             &
            &                                                 FGSL_Root_fSolver_Brent        &
            &                  )
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandType              =rangeExpandMultiplicative    , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
            &                  )
    end if
    sphericalActive              => self
    sphericalMassTarget          =  mass
    sphericalRadiusEnclosingMass =  finder%find(rootGuess=1.0d0)
    return
  end function sphericalRadiusEnclosingMass

  double precision function sphericalRadiusHalfMass(self)
    !% Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    implicit none
    class(massDistributionSpherical), intent(inout) :: self

    sphericalRadiusHalfMass=self%radiusEnclosingMass(0.5d0*self%massTotal())
    return
  end function sphericalRadiusHalfMass

  double precision function sphericalMassRoot(radius)
    !% Root function used in finding half mass radii of spherically symmetric mass distributions.
    implicit none
    double precision, intent(in   ) :: radius

    sphericalMassRoot=+sphericalActive%massEnclosedBySphere(radius) &
         &            -sphericalMassTarget
    return
  end function sphericalMassRoot

  function sphericalPositionSample(self,randomNumberGenerator_)
    !% Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                            , dimension(3)  :: sphericalPositionSample
    class           (massDistributionSpherical ), intent(inout) :: self
    class           (randomNumberGeneratorClass), intent(inout) :: randomNumberGenerator_
    double precision                                            :: mass                   , radius, &
         &                                                         theta                  , phi

    ! Choose an enclosed mass and find the radius enclosing that mass. Choose angular
    ! coordinates at random and finally convert to Cartesian.
    mass  =+              self                  %massTotal          (    )        &
         & *              randomNumberGenerator_%uniformSample      (    )
    radius=+              self                  %radiusEnclosingMass(mass)
    phi   =+     2.0d0*Pi*randomNumberGenerator_%uniformSample      (    )
    theta =+acos(2.0d0   *randomNumberGenerator_%uniformSample      (    )-1.0d0)
    sphericalPositionSample=+radius                &
         &                  *[                     &
         &                    sin(theta)*cos(phi), &
         &                    sin(theta)*sin(phi), &
         &                    cos(theta)           &       
         &                   ]
    return
  end function sphericalPositionSample
