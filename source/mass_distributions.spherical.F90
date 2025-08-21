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
  Implementation of an abstract mass distribution class for spherically symmetric distributions.
  !!}

  !![
  <massDistribution name="massDistributionSpherical" abstract="yes">
   <description>An abstract mass distribution class for spherically symmetric distributions.</description>
  </massDistribution>
  !!]
  type, extends(massDistributionClass), abstract :: massDistributionSpherical
     !!{
     Implementation of an abstract mass distribution class for spherically symmetric distributions.
     !!}
     private
     ! Memoized solutions for the enclosed mass.
     double precision              , allocatable, dimension(:) :: massProfileMass__                                          , massProfileRadius__
     double precision                                          :: massProfileRadiusMinimum__                    =+huge(0.0d0), massProfileRadiusMaximum__     =-huge(0.0d0)
     type            (interpolator), allocatable               :: massProfile__
     logical                                                   :: tolerateEnclosedMassIntegrationFailure        =.false.
     ! Memoized solutions for the potential.
     double precision              , allocatable, dimension(:) :: potentialProfilePotential__                                , potentialProfileRadius__
     double precision                                          :: potentialProfileRadiusMinimum__               =+huge(0.0d0), potentialProfileRadiusMaximum__=-huge(0.0d0), &
          &                                                       potentialProfileRadiusMinimumActual__         =+huge(0.0d0), potentialRadiusZeroPoint__     =-huge(0.0d0)
     type            (interpolator), allocatable               :: potentialProfile__
     logical                                                   :: toleratePotentialIntegrationFailure           =.false.
     double precision                                          :: toleranceRelativePotential                    =1.0d-6
     ! Options controlling implementation.
     logical                                                   :: chandrasekharIntegralComputeVelocityDispersion=.true.
   contains
     !![
     <methods>
       <method description="Returns the radius enclosing half of the mass of the mass distribution."                               method="radiusHalfMass"                     />
       <method description="Compute the potential energy of mass distribution."                                                    method="energyPotential"                    />
       <method description="Compute the kinetic energy of the mass distribution."                                                  method="energyKinetic"                      />
       <method description="Compute (numerically) the potential energy of mass distribution."                                      method="energyPotentialNumerical"           />
       <method description="Compute (numerically) the kinetic energy of the mass distribution."                                    method="energyKineticNumerical"             />
       <method description="Compute (numerically) a radial densitymoment of the mass distribution."                                method="densityRadialMomentNumerical"       />
       <method description="Compute (numerically) the radial density gradient."                                                    method="densityGradientRadialNumerical"     />
       <method description="Compute (numerically) the mass enclosed by a sphere."                                                  method="massEnclosedBySphereNumerical"      />
       <method description="Compute (numerically) the rotation curve."                                                             method="rotationCurveNumerical"             />
       <method description="Compute (numerically) the gravitational potential."                                                    method="potentialNumerical"                 />
       <method description="Compute (numerically) the Fourier transform of the density profile."                                   method="fourierTransformNumerical"          />
       <method description="Compute (numerically) the freefall radius."                                                            method="radiusFreefallNumerical"            />
       <method description="Compute (numerically) the growth rate of the freefall radius."                                         method="radiusFreefallIncreaseRateNumerical"/>
       <method description="Compute (numerically) the energy."                                                                     method="energyNumerical"                    />
       <method description="Integrand for dark matter profile potential."                                                          method="potentialSolverIntegrand"           />
       <method description="Return the radius variable used in solving the potential that corresponds to a given physical radius." method="potentialSolverRadius"              />
       <method description="Set sub-module scope pointers on a stack to allow recursive calls to functions."                       method="solverSphericalSet"                 />
       <method description="Unset sub-module scope pointers on a stack."                                                           method="solverSphericalUnset"               />
     </methods>
     !!]
     procedure :: symmetry                            => sphericalSymmetry
     procedure :: isSphericallySymmetric              => sphericalIsSphericallySymmetric
     procedure :: densityGradientRadial               => sphericalDensityGradientRadial
     procedure :: densityGradientRadialNumerical      => sphericalDensityGradientRadialNumerical
     procedure :: massEnclosedBySphere                => sphericalMassEnclosedBySphere
     procedure :: massEnclosedBySphereNumerical       => sphericalMassEnclosedBySphereNumerical
     procedure :: radiusEnclosingMassNumerical        => sphericalRadiusEnclosingMassNumerical
     procedure :: densityRadialMoment                 => sphericalDensityRadialMoment
     procedure :: densityRadialMomentNumerical        => sphericalDensityRadialMomentNumerical
     procedure :: potential                           => sphericalPotential
     procedure :: potentialNumerical                  => sphericalPotentialNumerical
     procedure :: potentialSolverIntegrand            => sphericalPotentialSolverIntegrand
     procedure :: potentialSolverRadius               => sphericalPotentialSolverRadius
     procedure :: fourierTransform                    => sphericalFourierTransform
     procedure :: fourierTransformNumerical           => sphericalFourierTransformNumerical
     procedure :: radiusFreefall                      => sphericalRadiusFreefall
     procedure :: radiusFreefallNumerical             => sphericalRadiusFreefallNumerical
     procedure :: radiusFreefallIncreaseRate          => sphericalRadiusFreefallIncreaseRate
     procedure :: radiusFreefallIncreaseRateNumerical => sphericalRadiusFreefallIncreaseRateNumerical
     procedure :: energy                              => sphericalEnergy
     procedure :: energyNumerical                     => sphericalEnergyNumerical
     procedure :: energyPotential                     => sphericalEnergyPotential
     procedure :: energyKinetic                       => sphericalEnergyKinetic
     procedure :: energyPotentialNumerical            => sphericalEnergyPotentialNumerical
     procedure :: energyKineticNumerical              => sphericalEnergyKineticNumerical
     procedure :: densitySphericalAverage             => sphericalDensitySphericalAverage
     procedure :: potentialDifferenceNumerical        => sphericalPotentialDifferenceNumerical
     procedure :: surfaceDensity                      => sphericalSurfaceDensity
     procedure :: chandrasekharIntegral               => sphericalChandrasekharIntegral
     procedure :: acceleration                        => sphericalAcceleration
     procedure :: tidalTensor                         => sphericalTidalTensor
     procedure :: positionSample                      => sphericalPositionSample
     procedure :: rotationCurve                       => sphericalRotationCurve
     procedure :: rotationCurveNumerical              => sphericalRotationCurve
     procedure :: rotationCurveGradient               => sphericalRotationCurveGradient
     procedure :: radiusHalfMass                      => sphericalRadiusHalfMass
     procedure :: solverSphericalSet                  => sphericalSolverSphericalSet
     procedure :: solverSphericalUnset                => sphericalSolverSphericalUnset
  end type massDistributionSpherical

  ! Submodule-scope pointers used in integrand functions and root finding.
  type :: massSolverSpherical
     class           (massDistributionSpherical), pointer :: self => null()
     double precision                                     :: time          , radiusFreefall
  end type massSolverSpherical
  type   (massSolverSpherical), allocatable, dimension(:) :: massSphericalSolvers
  integer                     , parameter                 :: massSphericalSolversIncrement=10
  integer                                                 :: massSphericalSolversCount    = 0
  !$omp threadprivate(massSphericalSolvers,massSphericalSolversCount)
  
contains

  function sphericalSymmetry(self)
    !!{
    Returns symmetry label for mass distributions with spherical symmetry.
    !!}
    implicit none
    type (enumerationMassDistributionSymmetryType)                :: sphericalSymmetry
    class(massDistributionSpherical              ), intent(inout) :: self
    !$GLC attributes unused :: self

    sphericalSymmetry=massDistributionSymmetrySpherical
    return
  end function sphericalSymmetry

  logical function sphericalIsSphericallySymmetric(self) result(isSphericallySymmetric)
    !!{
    Return true if the distribution is spherically symmetric.
    !!}
    implicit none
    class(massDistributionSpherical), intent(inout) :: self

    isSphericallySymmetric=.true.
    return
  end function sphericalIsSphericallySymmetric

  double precision function sphericalDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSpherical), intent(inout), target   :: self
    class  (coordinate               ), intent(in   )           :: coordinates
    logical                           , intent(in   ), optional :: logarithmic
    
    densityGradient=self%densityGradientRadialNUmerical(coordinates,logarithmic)
    return
  end function sphericalDensityGradientRadial

  double precision function sphericalDensityGradientRadialNumerical(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution using a numerical calculation.
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (massDistributionSpherical   ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    double precision                              , parameter               :: radiusLogarithmicStep=0.1d0
    type            (differentiator              )                          :: differentiator_
    double precision                                                        :: radius
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    call self%solverSet  ()
    radius          =   coordinates    %rSpherical(                                     )
    differentiator_ =   differentiator            (densityEvaluate                      )
    densityGradient =  +differentiator_%derivative(log(radius)    ,radiusLogarithmicStep)
    call self%solverUnset()
    if (.not.logarithmic_)                                    &
         & densityGradient=+     densityGradient              &
         &                 *self%density        (coordinates) &
         &                 /     radius
    return
  end function sphericalDensityGradientRadialNumerical

  double precision function densityEvaluate(radiusLogarithmic) result(density)
    !!{
    GSL-callable function to evaluate the density of the dark matter profile.
    !!}
      use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ), value :: radiusLogarithmic
    type            (coordinateSpherical)                       :: coordinates

    coordinates=[exp(radiusLogarithmic),0.0d0,0.0d0]
    density    =log(massSolvers(massSolversCount)%self%density(coordinates))
    return
  end function densityEvaluate

  double precision function sphericalMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: radius

    mass=self%massEnclosedBySphereNumerical(radius)
    return
  end function sphericalMassEnclosedBySphere

  double precision function sphericalMassEnclosedBySphereNumerical(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions using numerical integration.
    !!}
    use, intrinsic :: ISO_C_Binding           , only : c_size_t
    use            :: Numerical_Constants_Math, only : Pi
    use            :: Numerical_Integration   , only : integrator
    use            :: Numerical_Ranges        , only : Make_Range                  , rangeTypeLogarithmic
    use            :: Table_Labels            , only : extrapolationTypeExtrapolate
    use            :: Numerical_Interpolation , only : gsl_interp_linear
    use            :: Error                   , only : Error_Report                , errorStatusSuccess
    implicit none
    class           (massDistributionSpherical), intent(inout), target     :: self
    double precision                           , intent(in   )             :: radius
    double precision                                         , parameter   :: countPointsPerOctave     =4.0d+00
    double precision                                         , parameter   :: radiusVirialFractionSmall=1.0d-12
    double precision                          , dimension(:) , allocatable :: masses                           , radii              , &
         &                                                                    masses_                          , radii_
    type            (integrator              )                             :: integrator_
    integer         (c_size_t                )                             :: countRadii                       , iMinimum           , &
         &                                                                    iMaximum                         , i                  , &
         &                                                                    iPrevious
    logical                                                                :: remakeTable
    double precision                                                       :: radiusIntegralLower              , radiusIntegralUpper, &
         &                                                                    radiusMinimum                    , radiusMaximum
    integer                                                                :: status

    if (radius <= 0.0d0) then
       mass=0.0d0
       return
    end if
    ! Determine if the table must be rebuilt.
    remakeTable=.false.
    if (.not.allocated(self%massProfileMass__)) then
       remakeTable=.true.
    else
       remakeTable= radius < self%massProfileRadiusMinimum__ &
            &      .or.                                      &
            &       radius > self%massProfileRadiusMaximum__
    end if
    if (remakeTable) then
       ! Set sub-module-scope pointers.
       call self%solverSet()
       ! Initialize integrator.
       integrator_=integrator(sphericalMassEnclosedBySphereIntegrand,toleranceRelative=1.0d-2)
       ! Find the range of radii at which to compute the enclosed mass, and construct the arrays.
       !! Set an initial range of radii that brackets the requested radii.
       radiusMinimum=0.5d0*radius
       radiusMaximum=2.0d0*radius
       !! Round to the nearest factor of 2.
       radiusMinimum=2.0d0**floor  (log(radiusMinimum)/log(2.0d0))
       radiusMaximum=2.0d0**ceiling(log(radiusMaximum)/log(2.0d0))
       !! Expand to encompass any pre-existing range.
       if (allocated(self%massProfileRadius__)) then
          radiusMinimum=min(radiusMinimum,self%massProfileRadiusMinimum__)
          radiusMaximum=max(radiusMaximum,self%massProfileRadiusMaximum__)
       end if
       !! Construct arrays.
       countRadii=nint(log(radiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave+1.0d0)
       allocate(radii (countRadii))
       allocate(masses(countRadii))
       radii=Make_Range(radiusMinimum,radiusMaximum,int(countRadii),rangeTypeLogarithmic)
       ! Copy in any usable results from any previous solution.
       !! Assume by default that no previous solutions are usable.
       iMinimum=+huge(0_c_size_t)
       iMaximum=-huge(0_c_size_t)
       !! Check that a pre-existing solution exists.
       if (allocated(self%massProfileRadius__)) then
          iMinimum=nint(log(self%massProfileRadiusMinimum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
          iMaximum=nint(log(self%massProfileRadiusMaximum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
          masses(iMinimum:iMaximum)=self%massProfileMass__
       end if
       !! If ignoring errors in mass profile tabulation we must always compute the profile completely if we have a new minimum
       !! radius. Failure to do so can lead to non-monotonically increasing mass profiles (which lead to negative densities) due
       !! to prior failures in computing the mass profile.
       if (self%tolerateEnclosedMassIntegrationFailure .and. iMinimum > 1_c_size_t) then
          iMinimum=+huge(0_c_size_t)
          iMaximum=-huge(0_c_size_t)
       end if
       ! Solve for the enclosed mass where old results were unavailable.
       do i=1,countRadii
          ! Skip cases for which we have a pre-existing solution.
          if (i >= iMinimum .and. i <= iMaximum) cycle
          ! Find the limits for the integral.
          if (i == 1) then
             radiusIntegralLower=0.0d0
          else
             radiusIntegralLower=radii(i-1)
          end if
          radiusIntegralUpper   =radii(i  )
          ! Evaluate the integral.
          masses           (i)=+4.0d0*Pi*integrator_%integrate(radiusIntegralLower,radiusIntegralUpper,status=status)
          if (status /= errorStatusSuccess .and. .not.self%tolerateEnclosedMassIntegrationFailure) &
               & call Error_Report('enclosed mass profile integration failed'//{introspection:location})
          if (i > 1) masses(i)=+masses(i  ) &
               &               +masses(i-1)
       end do
       ! Build the interpolator.
       if (allocated(self%massProfile__)) deallocate(self%massProfile__)
       if (all(masses > 0.0d0)) then
          allocate(self%massProfile__)
          self%massProfile__=interpolator(log(radii),log(masses),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
       else if (all(masses == 0.0d0)) then
          ! This is a fully-destroyed profile - we leave the mass profile unallocated to indicate this.
       else
          ! We have a partial profile. Attempt to work with this by selecting out the physically-reasonable values.
          countRadii=+0_c_size_t
          iPrevious =-1_c_size_t
          do i=size(masses),1,-1
             if     (                                         &
                  &     masses   (i) >  0.0d0                 &
                  &  .and.                                    &
                  &   (                                       &
                  &     iPrevious    <= 0_c_size_t            &
                  &    .or.                                   &
                  &     masses   (i) <  masses    (iPrevious) &
                  &   )                                       &
                  & ) then
                countRadii=countRadii+1_c_size_t
                iPrevious =           i
             end if
          end do
          if (countRadii > 1) then
             allocate(radii_ (countRadii))
             allocate(masses_(countRadii))
             iPrevious=-1_c_size_t
             do i=size(masses),1,-1
                if     (                                         &
                     &     masses   (i) >  0.0d0                 &
                     &  .and.                                    &
                     &   (                                       &
                     &     iPrevious    <= 0_c_size_t            &
                     &    .or.                                   &
                     &     masses   (i) <  masses    (iPrevious) &
                     &   )                                       &
                     & ) then
                   radii_    (countRadii)=radii     (i)
                   masses_   (countRadii)=masses    (i)
                   iPrevious             =           i
                   countRadii            =countRadii-1_c_size_t
                end if
             end do
             allocate(self%massProfile__)
             self%massProfile__=interpolator(log(radii_),log(masses_),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
          end if
       end if
       call self%solverUnset()
       ! Store the current results for future re-use.
       if (allocated(self%massProfileRadius__)) deallocate(self%massProfileRadius__)
       if (allocated(self%massProfileMass__  )) deallocate(self%massProfileMass__  )
       allocate(self%massProfileRadius__(countRadii))
       allocate(self%massProfileMass__  (countRadii))
       self%massProfileRadius__       =radii
       self%massProfileMass__         =masses
       self%massProfileRadiusMinimum__=radiusMinimum
       self%massProfileRadiusMaximum__=radiusMaximum
    end if
    ! Interpolate in the table to find the mass difference.
    if (allocated(self%massProfile__)) then
       mass=+exp(self%massProfile__%interpolate(log(radius)))
    else
       ! Fully-destroyed profile.
       mass=+0.0d0
    end if    
    return
  end function sphericalMassEnclosedBySphereNumerical

  double precision function sphericalMassEnclosedBySphereIntegrand(radius) result(integrand)
    !!{
    Enclosed mass integrand for spherical mass distributions.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: position

    if (radius > 0.0d0) then
       position =[radius,0.0d0,0.0d0]
       integrand=+radius**2               &
            &    *massSolvers(massSolversCount)%self%density(position)
    else
      integrand=+0.0d0 
    end if
    return
  end function sphericalMassEnclosedBySphereIntegrand

  double precision function sphericalRadiusEnclosingMassNumerical(self,mass,massFractional) result(radius)
    !!{
    Return the radius enclosing a specified mass using a numerical calculation.
    !!}
    use :: Root_Finder, only : rootFinder                   , GSL_Root_fSolver_Brent, rangeExpandMultiplicative, rangeExpandSignExpectNegative, &
         &                     rangeExpandSignExpectPositive
    use :: Error      , only : Error_Report
    implicit none
    class           (massDistributionSpherical), intent(inout), target   :: self
    double precision                           , intent(in   ), optional :: mass                   , massFractional
    type            (rootFinder               )                          :: finder
    double precision                           , parameter               :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    double precision                                                     :: massTarget
    
    if      (present(mass          )) then
       massTarget=     mass
    else if (present(massFractional)) then
       massTarget=self%massTotal()*massFractional
    else
       massTarget=0.0d0
       call Error_Report('either "mass" or "massFractional" must be provided'//{introspection:location})
    end if
    if (massTarget <= 0.0d0 .or. self%massEnclosedBySphere(0.0d0) >= massTarget) then
       radius=0.0d0
       return
    end if
    finder =rootFinder(                                                             &
         &             rootFunction                 =sphericalMassEnclosedRoot    , &
         &             toleranceAbsolute            =toleranceAbsolute            , &
         &             toleranceRelative            =toleranceRelative            , &
         &             solverType                   =GSL_Root_fSolver_Brent       , &
         &             rangeExpandUpward            =2.0d0                        , &
         &             rangeExpandDownward          =0.5d0                        , &
         &             rangeExpandType              =rangeExpandMultiplicative    , &
         &             rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
         &             rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
         &            )
    call self%solverSet  (massTarget=massTarget)
    if (allocated(self%massProfile__)) then
       ! A mass profile exists - use its extent as the initial search range. If the required radius lies within this range we then
       ! avoid triggering any retabulation.
       radius=finder%find(rootRange=[self%massProfileRadiusMinimum__,self%massProfileRadiusMaximum__])
    else
       ! No mass profile exists - make a reasonable guess at a radius to begin.
       radius=finder%find(rootGuess=1.0d0                                                            )
    end if
    call self%solverUnset(                     )
    return
  end function sphericalRadiusEnclosingMassNumerical
  
  double precision function sphericalMassEnclosedRoot(radius) result(root)
    !!{
    Root function used in finding radii enclosing a target mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    root   =+massSolvers(massSolversCount)%self%massEnclosedBySphere(radius) &
         &  -massSolvers(massSolversCount)     %massTarget
    return
  end function sphericalMassEnclosedRoot
  
  double precision function sphericalDensitySphericalAverage(self,radius)
    !!{
    Computes the density averaged over a spherical shell.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (coordinateSpherical      )                :: position

    ! For a spherical mass distribution, the density averaged over a spherical shell, is just the regular density at that radius.
    position                        =[radius,0.0d0,0.0d0]
    sphericalDensitySphericalAverage=self%density(position)
    return
  end function sphericalDensitySphericalAverage

  double precision function sphericalRadiusHalfMass(self)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    implicit none
    class(massDistributionSpherical), intent(inout) :: self

    sphericalRadiusHalfMass=self%radiusEnclosingMass(0.5d0*self%massTotal())
    return
  end function sphericalRadiusHalfMass

  double precision function sphericalPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSpherical        ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status

    potential=self%potentialNumerical(coordinates,status)
    return
  end function sphericalPotential

  double precision function sphericalPotentialNumerical(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Coordinates                     , only : assignment(=)
    use            :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use            :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use            :: Numerical_Integration           , only : integrator
    use            :: Numerical_Comparison            , only : Values_Agree
    use            :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLogarithmic
    use            :: Table_Labels                    , only : extrapolationTypeExtrapolate
    use            :: Numerical_Interpolation         , only : gsl_interp_linear
    use            :: Error                           , only : Error_Report                  , errorStatusSuccess
    implicit none
    class           (massDistributionSpherical        ), intent(inout), target      :: self
    class           (coordinate                       ), intent(in   )              :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional    :: status
    double precision                                   , dimension(:) , allocatable :: potentials                  , radii               , &
         &                                                                             potentials_                 , radii_
    double precision                                   , parameter                  :: countPointsPerOctave=4.0d+0
    type            (integrator                       )                             :: integrator_
    integer         (c_size_t                )                                      :: countRadii                  , iMinimum            , &
         &                                                                             iMaximum                    , iPrevious           , &
         &                                                                             i
    logical                                                                         :: remakeTable
    double precision                                                                :: radiusMaximum               , radiusMinimum       , &
         &                                                                             radius                      , massRadiusMaximum   , &
         &                                                                             radiusLowerPotential        , radiusUpperPotential
    integer                                                                         :: status_

    if (present(status)) status=structureErrorCodeSuccess
    potential         =+0.0d0
    radius            =+coordinates%rSpherical()
    radiusMaximum     =+            radius
    ! Round to the nearest factor of 2.
    radiusMaximum=2.0d0**ceiling(log(radiusMaximum)/log(2.0d0))
    ! Determine if the table must be rebuilt.
    remakeTable=.false.
    if (.not.allocated(self%potentialProfilePotential__)) then
       remakeTable=.true.
    else
       remakeTable= radius < self%potentialProfileRadiusMinimum__ &
            &      .or.                                           &
            &       radius > self%potentialProfileRadiusMaximum__
    end if
    if (remakeTable) then
       integrator_=integrator(integrandPotential,toleranceRelative=self%toleranceRelativePotential)
       ! Find the range of radii at which to compute the enclosed mass, and construct the arrays.
       !! Check the mass at the maximum radius is non-zero.
       massRadiusMaximum=self%massEnclosedBySphere(radiusMaximum)
       if (massRadiusMaximum > 0.0d0) then
          !! Set an initial range of radii that brackets the requested radii.
          radiusMinimum=0.5d0*radius
          !! Round to the nearest factor of 2.
          radiusMinimum=2.0d0**floor(log(radiusMinimum)/log(2.0d0))
          !! Expand to encompass any pre-existing range.
          if (allocated(self%potentialProfileRadius__)) then
             radiusMinimum=min(radiusMinimum,self%potentialProfileRadiusMinimum__)
             radiusMaximum=max(radiusMaximum,self%potentialProfileRadiusMaximum__)
          end if
          !! Construct arrays.
          countRadii=nint(log(radiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave+1.0d0)
          allocate(radii     (countRadii))
          allocate(potentials(countRadii))
          radii=Make_Range(radiusMinimum,radiusMaximum,int(countRadii),rangeTypeLogarithmic)
          ! Copy in any usable results from any previous solution.
          !! Assume by default that no previous solutions are usable.
          iMinimum=+huge(0_c_size_t)
          iMaximum=-huge(0_c_size_t)
          !! Check that a pre-existing solution exists.
          if (allocated(self%potentialProfileRadius__)) then
             iMinimum=nint(log(self%potentialProfileRadiusMinimum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             iMaximum=nint(log(self%potentialProfileRadiusMaximum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             potentials(iMinimum:iMaximum)=self%potentialProfilePotential__
          end if
          ! Set the radius for the zero point of the potential.
          if (self%potentialRadiusZeroPoint__ < 0.0d0) self%potentialRadiusZeroPoint__=radiusMaximum
          ! Solve for the enclosed mass where old results were unavailable.
          call self%solverSet  ()
          do i=countRadii,1,-1
             ! Skip cases for which we have a pre-existing solution.
             if (i >= iMinimum .and. i <= iMaximum) cycle  
             ! Evaluate the integral.
             radiusLowerPotential   =self%potentialSolverRadius(     radii                     (i  ))
             ! If the prior point in the grid was successfully evaluated, we need only compute the potential difference relative
             ! to that point.
             if (i == countRadii .or. potentials(i+1) == -huge(0.0d0)) then
                radiusUpperPotential=self%potentialSolverRadius(self%potentialRadiusZeroPoint__     )
             else
                radiusUpperPotential=self%potentialSolverRadius(     radii                     (i+1))
             end if
             potentials(i)=integrator_%integrate(radiusLowerPotential,radiusUpperPotential,status_)
             ! Convert a potential difference to the actual potential if necessary.
             if (i < countRadii .and. potentials(i+1) /= -huge(0.0d0)) &
                  & potentials(i)=+potentials(i  ) &
                  &               +potentials(i+1)
             if (status_ /= errorStatusSuccess) then
                if (self%toleratePotentialIntegrationFailure) then
                   potentials(i)=-huge(0.0d0)
                else
                   call Error_Report('potential integration failed'//{introspection:location})
                end if
             end if
          end do
          call self%solverUnset()
          ! Build the interpolator.
          if (allocated(self%potentialProfile__)) deallocate(self%potentialProfile__)
          if (all(potentials > -huge(0.0d0))) then
             radii_     =radii
             potentials_=potentials
          else if (all(potentials == -huge(0.0d0))) then
             ! A fully-destroyed profile.
             allocate(radii_     (0))
             allocate(potentials_(0))
          else
             !  We have a partial profile. Attempt to work with this by selecting out the physically-reasonable values.
             countRadii=+0_c_size_t
             iPrevious =-1_c_size_t
             do i=size(potentials),1,-1
                if     (                                          &
                     &     potentials(i) >  -huge(0.0d0)          &
                     &  .and.                                     &
                     &   (                                        &
                     &     iPrevious     <= 0_c_size_t            &
                     &    .or.                                    &
                     &     potentials(i) <  potentials(iPrevious) &
                     &   )                                        &
                     & ) then
                   countRadii=countRadii+1_c_size_t
                   iPrevious =           i
                end if
             end do
             if (countRadii > 1) then
                allocate(radii_     (countRadii))
                allocate(potentials_(countRadii))
                iPrevious=-1_c_size_t
                do i=size(potentials),1,-1
                   if     (                                          &
                        &     potentials(i) >  -huge(0.0d0)          &
                        &  .and.                                     &
                        &   (                                        &
                        &     iPrevious     <= 0_c_size_t            &
                        &    .or.                                    &
                        &     potentials(i) <  potentials(iPrevious) &
                        &   )                                        &
                        & ) then
                      radii_     (countRadii)=radii     (i)
                      potentials_(countRadii)=potentials(i)
                      iPrevious              =           i
                      countRadii             =countRadii-1_c_size_t
                   end if
                end do
             else
                ! Only one point is available - can not make an interpolator.
                allocate(radii_     (0))
                allocate(potentials_(0))
             end if
          end if
          if (size(radii_) > 0) then
             allocate(self%potentialProfile__)
             self%potentialProfile__=interpolator(radii_,potentials_,interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)  
             ! Store the current results for future re-use.
             if (allocated(self%potentialProfileRadius__   )) deallocate(self%potentialProfileRadius__   )
             if (allocated(self%potentialProfilePotential__)) deallocate(self%potentialProfilePotential__)
             allocate(self%potentialProfileRadius__   (countRadii))
             allocate(self%potentialProfilePotential__(countRadii))
             self%potentialProfileRadius__             =radii
             self%potentialProfilePotential__          =potentials
             self%potentialProfileRadiusMinimum__      =radiusMinimum
             self%potentialProfileRadiusMinimumActual__=radii_       (1)
             self%potentialProfileRadiusMaximum__      =radiusMaximum
          else
             ! The profile is undefined. Leave the table unallocated.
             if (allocated(self%potentialProfile__)) deallocate(self%potentialProfile__)
          end if
       else
          ! The profile is undefined. Leave the table unallocated.
          if (allocated(self%potentialProfile__)) deallocate(self%potentialProfile__)
       end if
    end if
    if (allocated(self%potentialProfile__)) then
       potential=self%potentialProfile__%interpolate(max(radius,self%potentialProfileRadiusMinimumActual__))
    else
       potential=+0.0d0
    end if
    ! Convert to dimensionful units.    
    if (.not.self%isDimensionless()) potential=+gravitationalConstant_internal &
         &                                     *potential
    return
  end function sphericalPotentialNumerical

  double precision function sphericalPotentialDifferenceNumerical(self,coordinates1,coordinates2,status) result(potential)
    !!{
    Return the potential difference between the two specified {\normalfont \ttfamily coordinates} in a spherical mass distribution
    using a numerical calculation.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Integration           , only : integrator
    use :: Numerical_Comparison            , only : Values_Agree
    implicit none
    class           (massDistributionSpherical        ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates1            , coordinates2
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                   , parameter               :: toleranceRelative=1.0d-3
    type            (integrator                       )                          :: integrator_

    if (present(status)) status=structureErrorCodeSuccess
    integrator_=integrator(integrandPotential,toleranceRelative=toleranceRelative)
    call self%solverSet  ()
    potential =  integrator_%integrate(                           &
         &                             coordinates1%rSpherical(), &
         &                             coordinates2%rSpherical()  &
         &                            )
    call self%solverUnset()
    ! Convert to dimensionful units.    
    if (.not.self%isDimensionless()) potential=+gravitationalConstant_internal &
         &                                     *potential
    return
  end function sphericalPotentialDifferenceNumerical
  
  double precision function integrandPotential(radius)
    !!{
    Integrand for gravitational potential in a generic dark matter profile.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    select type (self => massSolvers(massSolversCount)%self)
    class is (massDistributionSpherical)
       integrandPotential=self%potentialSolverIntegrand(radius)
    class default
       integrandPotential=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function integrandPotential

  double precision function sphericalPotentialSolverIntegrand(self,radius) result(integrandPotential)
    !!{
    Integrand for gravitational potential in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
      
    if (radius > 0.0d0) then
       integrandPotential=-self%massEnclosedBySphere(radius)    &
            &             /                          radius **2
    else
       integrandPotential=0.0d0
    end if
    return
  end function sphericalPotentialSolverIntegrand

  double precision function sphericalPotentialSolverRadius(self,radius) result(radiusSolver)
    !!{
    Return the radius variable used in computing the potential that corresponds to a given physical radius.
    In some cases, it is easier to do the integration with respect to another variable which is a function of
    the physical radius.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    !$GLC attributes unused :: self

    radiusSolver=radius
    return
  end function sphericalPotentialSolverRadius

  function sphericalAcceleration(self,coordinates) result(acceleration)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=), coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gigaYear     , megaParsec         , gravitationalConstant_internal
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                              , dimension(3)  :: acceleration
    class           (massDistributionSpherical   ), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    type            (coordinateSpherical         )                :: coordinatesSpherical
    type            (coordinateCartesian         )                :: coordinatesCartesian
    double precision                                              :: radius
    double precision                              , dimension(3)  :: positionCartesian

    ! Get position in spherical and Cartesian coordinate systems.
    coordinatesSpherical=coordinates
    coordinatesCartesian=coordinates
    ! Compute the density at this position.
    positionCartesian=coordinatesCartesian
    radius           =coordinatesSpherical%r()
    if (radius > 0.0d0) then
       acceleration=-self%massEnclosedBySphere(radius           )    &
            &       *                          positionCartesian     &
            &       /                          radius            **3
       if (.not.self%isDimensionless())                    &
            & acceleration=+acceleration                   &
            &              *kilo                           &
            &              *gigaYear                       &
            &              /megaParsec                     &
            &              *gravitationalConstant_internal
    else
       acceleration=0.0d0
    end if
    return
  end function sphericalAcceleration

  function sphericalTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for spherically-symmetric mass
    distributions.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateSpherical, coordinateCartesian
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Tensors                         , only : tensorIdentityR2D3Sym         , assignment(=)      , operator(*)
    use :: Vectors                         , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: sphericalTidalTensor
    class           (massDistributionSpherical     ), intent(inout) :: self
    class           (coordinate                    ), intent(in   ) :: coordinates
    type            (coordinateSpherical           )                :: coordinatesSpherical
    type            (coordinateCartesian           )                :: coordinatesCartesian
    double precision                                                :: radius              , massEnclosed, &
         &                                                             density
    double precision                                , dimension(3)  :: positionCartesian
    type            (tensorRank2Dimension3Symmetric)                :: positionTensor

    ! Get position in spherical and Cartesian coordinate systems.
    coordinatesSpherical=coordinates
    coordinatesCartesian=coordinates
    ! Compute the enclosed mass and density at this position.
    positionCartesian= coordinatesCartesian
    radius           =+coordinatesSpherical%r                   (           )
    massEnclosed     =+self                %massEnclosedBySphere(radius     ) 
    density          =+self                %density             (coordinates) 
    positionTensor   =Vector_Outer_Product(positionCartesian,symmetrize=.true.)
    ! Find the gravitational tidal tensor.
    sphericalTidalTensor=-(massEnclosed         /radius**3)*tensorIdentityR2D3Sym &
         &               +(massEnclosed*3.0d0   /radius**5)*positionTensor        &
         &               -(density     *4.0d0*Pi/radius**2)*positionTensor
    ! For dimensionful profiles, add the appropriate normalization.
    if (.not.self%isDimensionless())                             &
         & sphericalTidalTensor=+sphericalTidalTensor            &
         &                       *gravitationalConstant_internal
    return
  end function sphericalTidalTensor

  double precision function sphericalRotationCurve(self,radius)
    !!{
    Return the rotation curve for a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius

    if (radius <= 0.0d0) then
       sphericalRotationCurve=+0.0d0
    else
       sphericalRotationCurve=+sqrt(                                   &
            &                       +self%massEnclosedBySphere(radius) &
            &                       /                          radius  &
            &                      )
       ! Make dimensionful if necessary.
       if (.not.self%dimensionless) sphericalRotationCurve= &
            & +sqrt(gravitationalConstant_internal)         &
            & *sphericalRotationCurve
    end if
    return
  end function sphericalRotationCurve

  double precision function sphericalRotationCurveGradient(self,radius)
    !!{
    Return the rotation curve gradient for a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (coordinateSpherical      )                :: position

    position                      =[radius,0.0d0,0.0d0]
    sphericalRotationCurveGradient=+4.0d0                                  &
         &                         *Pi                                     &
         &                         *                          radius       &
         &                         *self%density             (position)    &
         &                         -self%massEnclosedBySphere(radius  )    &
         &                         /                          radius   **2
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) sphericalRotationCurveGradient= &
         &  +gravitationalConstant_internal                      &
         &  *sphericalRotationCurveGradient
    return
  end function sphericalRotationCurveGradient

  function sphericalPositionSample(self,randomNumberGenerator_)
    !!{
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                              , dimension(3)  :: sphericalPositionSample
    class           (massDistributionSpherical   ), intent(inout) :: self
    class           (randomNumberGeneratorClass  ), intent(inout) :: randomNumberGenerator_
    double precision                                              :: mass                   , radius, &
         &                                                           theta                  , phi

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

  double precision function sphericalSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(massDistributionSpherical), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates

    sphericalSurfaceDensity=0.0d0
    call Error_Report('surface density is not defined for spherically-symmetric distributions'//{introspection:location})
    return
  end function sphericalSurfaceDensity

  function sphericalChandrasekharIntegral(self,massDistributionEmbedding,massDistributionPerturber,massPerturber,coordinates,velocity) result(integral)
    !!{
    Compute the Chandrasekhar integral at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    use :: Coordinates               , only : coordinateCartesian  , assignment(=)
    use :: Numerical_Constants_Math  , only : Pi
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeAll
    use :: Ideal_Gases_Thermodynamics, only : Ideal_Gas_Sound_Speed
    use :: Error                     , only : Error_Report
    implicit none
    double precision                           , dimension(3)  :: integral
    class           (massDistributionSpherical), intent(inout) :: self
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding           , massDistributionPerturber
    double precision                           , intent(in   ) :: massPerturber
    class           (coordinate               ), intent(in   ) :: coordinates                         , velocity
    double precision                           , dimension(3)  :: velocityCartesian_
    double precision                           , parameter     :: XvMaximum                    =10.0d0
    type            (coordinateCartesian      )                :: velocityCartesian
    double precision                                           :: radius                              , velocity_                , &
         &                                                        density                             , velocityDispersion       , &
         &                                                        factorSuppressionExtendedMass       , xV
    
    integral =0.0d0
    velocity_=velocity%rSpherical()
    if (velocity_ <= 0.0d0) return
    radius =coordinates%rSpherical(           )
    density=self       %density   (coordinates)
    if (density  <= 0.0d0) return
    if (.not.associated(self%kinematicsDistribution_)) call Error_Report('a kinematics distribution is needed to compute the Chandrasekhar integral'//{introspection:location})
    if (self%kinematicsDistribution_%isCollisional()) then
       velocityDispersion=Ideal_Gas_Sound_Speed(self%kinematicsDistribution_%temperature         (coordinates                               ))
    else
       velocityDispersion=                      self%kinematicsDistribution_%velocityDispersion1D(coordinates,self,massDistributionEmbedding)
    end if
    if (velocityDispersion > 0.0d0) then    
       xV             =+velocity_             &
            &          /velocityDispersion    &
            &          /sqrt(2.0d0)
    else
       xV             =+huge(0.0d0)
    end if
    velocityCartesian = velocity
    velocityCartesian_= velocityCartesian
    integral          =-density               &
         &             *velocityCartesian_    &
         &             /velocity_         **3
    if (Xv <= XvMaximum)                      &
         & integral   =+integral              &
         &             *(                     &
         &               +erf ( xV   )        &
         &               -2.0d0               &
         &               *      xV            &
         &               *exp (-xV**2)        &
         &               /sqrt( Pi   )        &
         &              )
    ! Compute suppression factor due to satellite being an extended mass distribution. This is largely untested - it is meant to
    ! simply avoid extremely large accelerations for subhalo close to the center of its host when that subhalo is much more
    ! extended than the host.
    factorSuppressionExtendedMass=min(1.0d0,massDistributionPerturber%massEnclosedBySphere(radius)/massPerturber)
    ! Evaluate the integral.
    integral=+integral                      &
         &   *factorSuppressionExtendedMass
    return
  end function sphericalChandrasekharIntegral

  double precision function sphericalFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter, wavenumber
 
    fourierTransform=self%fourierTransformNumerical(radiusOuter,wavenumber)
    return
  end function sphericalFourierTransform

  double precision function sphericalFourierTransformNumerical(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{   
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical mass
    distribution using a numerical calculation.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter, wavenumber
    type            (integrator               )                :: integrator_

    fourierTransform=0.0d0
    integrator_     = integrator                      (integrandFourierTransform,toleranceRelative=1.0d-3)
    fourierTransform=+integrator_%integrate           (0.0d0                    ,radiusOuter             ) &
         &           /self       %massEnclosedBySphere(                          radiusOuter             )
    return

  contains

    double precision function integrandFourierTransform(radius)
      !!{
      Integrand for Fourier transform of a spherical mass distribution.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      use :: Coordinates             , only : assignment(=), coordinateSpherical
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates
      
      if (radius > 0.0d0) then
         coordinates              =[radius,0.0d0,0.0d0]
         integrandFourierTransform=+4.0d0                     &
              &                    *Pi                        &
              &                    *               radius **2 &
              &                    *sin(wavenumber*radius)    &
              &                    /   (wavenumber*radius)    &
              &                    *self%density(coordinates)
      else
         integrandFourierTransform=0.0d0
      end if
      return
    end function integrandFourierTransform
    
  end function sphericalFourierTransformNumerical
  
  double precision function sphericalRadiusFreefall(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time
 
    radius=self%radiusFreefallNumerical(time)
    return
  end function sphericalRadiusFreefall
  
  double precision function sphericalRadiusFreefallNumerical(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily
    time} in a spherical mass distribution using a numerical
    calculation.
    !!}
    use :: Root_Finder                     , only : rangeExpandMultiplicative     , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: time
    double precision                           , parameter             :: toleranceAbsolute  =0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder               )                        :: finder
    type            (coordinateSpherical      )                        :: coordinates
    double precision                                                   :: timeFreefallMinimum

    radius=0.0d0
    coordinates=[0.0d0,0.0d0,0.0d0]
    if (self%densityGradientRadial(coordinates,logarithmic=.true.) == 0.0d0) then
       ! For mass distributions with a constant density core, the potential in the center is harmonic. This means there is a
       ! minimum to the freefall time as a function of radius. Compute that minimum here so that we can return a zero radius for
       ! times less than this.
       timeFreefallMinimum=+sqrt(                                &
            &                    + 3.0d0                         &
            &                    /16.0d0                         &
            &                    *Pi                             &
            &                    /gravitationalConstant_internal &
            &                    /self%density(coordinates)      &
            &                   )                                &
            &              *MpcPerKmPerSToGyr
    else
       timeFreefallMinimum=+0.0d0
    end if   
    if (time < timeFreefallMinimum) return
    call self%solverSphericalSet(time)
    finder =  rootFinder(                                                             &
         &               rootFunction                 =rootRadiusFreefall           , &
         &               toleranceAbsolute            =toleranceAbsolute            , &
         &               toleranceRelative            =toleranceRelative            , &
         &               rangeExpandDownward          =0.5d0                        , &
         &               rangeExpandUpward            =2.0d0                        , &
         &               rangeExpandType              =rangeExpandMultiplicative    , &
         &               rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &               rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &              )
    radius=finder%find(rootGuess=1.0d0)
    call self%solverSphericalUnset()
    return
  end function sphericalRadiusFreefallNumerical
  
  double precision function rootRadiusFreefall(radiusFreefall)
    !!{
    Root function used in finding the radius corresponding to a given freefall time.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radiusFreefall
    type            (integrator)                :: integrator_
      
    massSphericalSolvers(massSphericalSolversCount)%radiusFreefall=+                                                                               radiusFreefall
    integrator_                                                   = integrator                                              (integrandTimeFreefall,toleranceRelative=1.0d-3)
    rootRadiusFreefall                                            =+integrator_                                   %integrate(0.0d0                ,radiusFreefall          ) &
         &                                                         -massSphericalSolvers(massSphericalSolversCount)%time
    return
  end function rootRadiusFreefall

  double precision function integrandTimeFreefall(radius)
    !!{
    Integrand for freefall time in a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)    , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    double precision                     , intent(in   ) :: radius
    double precision                                     :: potentialDifference
    type            (coordinateSpherical)                :: coordinates        , coordinatesFreefall

    coordinates        =[                                                radius        ,0.0d0,0.0d0]
    coordinatesFreefall=[massSphericalSolvers(massSphericalSolversCount)%radiusFreefall,0.0d0,0.0d0]
    potentialDifference=+massSphericalSolvers(massSphericalSolversCount)%self%potentialDifference(coordinates,coordinatesFreefall)
    if (potentialDifference < 0.0d0) then
       integrandTimeFreefall=+MpcPerKmPerSToGyr         &
            &                /sqrt(                     &
            &                      -2.0d0               &
            &                      *potentialDifference &
            &                     )
    else
       ! Avoid floating point errors arising from rounding errors.
       integrandTimeFreefall=0.0d0
    end if
    return
  end function integrandTimeFreefall

  double precision function sphericalRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an spherical mass
    distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time

    radiusIncreaseRate=self%radiusFreefallIncreaseRateNumerical(time)
    return
  end function sphericalRadiusFreefallIncreaseRate

  double precision function sphericalRadiusFreefallIncreaseRateNumerical(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an spherical mass
    distribution using a numerical calculation.
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: time
    double precision                           , parameter     :: timeLogarithmicStep=1.0d-2
    type            (differentiator           )                :: differentiator_

    call self%solverSphericalSet(time)
    differentiator_   = differentiator            (radiusFreefallEvaluate                    )
    radiusIncreaseRate=+differentiator_%derivative(log(time)             ,timeLogarithmicStep) &
         &             /                               time
    call self%solverSphericalUnset()
    return
  end function sphericalRadiusFreefallIncreaseRateNumerical

  double precision function radiusFreefallEvaluate(timeLogarithmic)
    !!{
    GSL-callable function to evaluate the freefall radius of the mass distribution.
    !!}
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    radiusFreefallEvaluate=massSphericalSolvers(massSphericalSolversCount)%self%radiusFreefall(exp(timeLogarithmic))
    return
  end function radiusFreefallEvaluate

  double precision function sphericalEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout), target :: self
    double precision                           , intent(in   )         :: radiusOuter
    class           (massDistributionClass    ), intent(inout), target :: massDistributionEmbedding

    energy=self%energyNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergy

  double precision function sphericalEnergyNumerical(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution using a numerical calculation.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding

    energy=+self%energyPotential(radiusOuter                          ) &
         & +self%energyKinetic  (radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergyNumerical

  double precision function sphericalEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter

    energy=self%energyPotentialNumerical(radiusOuter)
    return
  end function sphericalEnergyPotential

  double precision function sphericalEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding

    energy=self%energyKineticNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalEnergyKinetic
  
  double precision function sphericalEnergyPotentialNumerical(self,radiusOuter) result(energy)
    !!{
    Compute (numerically) the potential energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Integration           , only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    type            (integrator               )                :: integrator_

    integrator_= integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
    energy     =-0.5d0                                                    &
         &      *gravitationalConstant_internal                           &
         &      *(                                                        &
         &        +integrator_%integrate           (0.0d0,radiusOuter)    &
         &        +self       %massEnclosedBySphere(      radiusOuter)**2 &
         &        /            radiusOuter                                &
         &       )
    return

  contains

    double precision function integrandEnergyPotential(radius)
      !!{
      Integrand for potential energy of a spherical mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandEnergyPotential=(                                   &
              &                    +self%massEnclosedBySphere(radius) &
              &                    /                          radius  &
              &                   )**2
      else
         integrandEnergyPotential=0.0d0
      end if
      return
    end function integrandEnergyPotential

  end function sphericalEnergyPotentialNumerical

  double precision function sphericalEnergyKineticNumerical(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute (numerically) the kinetic energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout) :: self
    double precision                           , intent(in   ) :: radiusOuter
    class           (massDistributionClass    ), intent(inout) :: massDistributionEmbedding
    type            (integrator               )                :: integrator_

    integrator_= integrator(integrandEnergyKinetic,toleranceRelative=1.0d-3)
    energy     =+6.0d0                                    &
         &      *Pi                                       &
         &      *integrator_%integrate(0.0d0,radiusOuter)
    return

  contains

    double precision function integrandEnergyKinetic(radius)
      !!{
      Integrand for kinetic energy of the halo.
      !!}
      use :: Coordinates, only : coordinateSpherical, assignment(=)
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates

      if (radius > 0.0d0) then
         coordinates           =[radius,0.0d0,0.0d0]
         integrandEnergyKinetic=+self                        %density             (coordinates                               )    &
              &                 *self%kinematicsDistribution_%velocityDispersion1D(coordinates,self,massDistributionEmbedding)**2 &
              &                 *                             radius                                                          **2
      else
         integrandEnergyKinetic=0.0d0
      end if
      return
    end function integrandEnergyKinetic

  end function sphericalEnergyKineticNumerical

  double precision function sphericalDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite

    densityRadialMoment=self%densityRadialMomentNumerical(moment,radiusMinimum,radiusMaximum,isInfinite)
    return
  end function sphericalDensityRadialMoment

  double precision function sphericalDensityRadialMomentNumerical(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Returns a radial density moment for a spherical mass distribution using a numerical calculation.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionSpherical), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite
    type            (integrator               )                          :: integrator_
    double precision                                                     :: radiusMinimum_
    
    densityRadialMoment=0.0d0
    if (.not.present(radiusMaximum)) call Error_Report('a maximum radius must be provided'//{introspection:location})
    if (present(isInfinite)) isInfinite=.false.
    radiusMinimum_=0.0d0
    if (present(radiusMinimum)) radiusMinimum_=radiusMinimum
    integrator_        = integrator           (integrandMoment,toleranceRelative=1.0d-3)
    densityRadialMoment=+integrator_%integrate(radiusMinimum_ ,radiusMaximum           )
    return

  contains

    double precision function integrandMoment(radius)
      !!{
      Integrand for radial density moment in a spherical mass distribution.
      !!}
      use :: Coordinates, only : assignment(=), coordinateSpherical
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates
      
      if (radius > 0.0d0) then
         coordinates     = [radius,0.0d0,0.0d0]
         integrandMoment=+self%density(coordinates)         &
              &          *             radius      **moment
      else
         integrandMoment=+0.0d0
      end if
      return
    end function integrandMoment
    
  end function sphericalDensityRadialMomentNumerical

  subroutine sphericalSolverSphericalSet(self,time)
    !!{
    Unset a sub-module scope pointers on the stack.
    !!}
    implicit none
    class           (massDistributionSpherical), intent(inout), target       :: self
    double precision                           , intent(in   )               :: time
    integer                                                                  :: i
    type            (massSolverSpherical      ), allocatable  , dimension(:) :: solvers_

    if (allocated(massSphericalSolvers)) then
       if (massSphericalSolversCount == size(massSphericalSolvers)) then
          call move_alloc(massSphericalSolvers,solvers_)
          allocate(massSphericalSolvers(size(solvers_)+massSphericalSolversIncrement))
          massSphericalSolvers(1:size(solvers_))=solvers_
          do i=1,size(solvers_)
             nullify(solvers_(i)%self)
          end do
          deallocate(solvers_)
       end if
    else
       allocate(massSphericalSolvers(massSphericalSolversIncrement))
    end if
    massSphericalSolversCount=massSphericalSolversCount+1
    massSphericalSolvers(massSphericalSolversCount)%self => self
    massSphericalSolvers(massSphericalSolversCount)%time =  time
    return
  end subroutine sphericalSolverSphericalSet

  subroutine sphericalSolverSphericalUnset(self)
    !!{
    Unset a sub-module scope pointers on the stack.
    !!}
    implicit none
    class(massDistributionSpherical), intent(inout) :: self
    !$GLC attributes unused :: self
    
    massSphericalSolvers(massSphericalSolversCount)%self => null()
    massSphericalSolversCount=massSphericalSolversCount-1
    return
  end subroutine sphericalSolverSphericalUnset
