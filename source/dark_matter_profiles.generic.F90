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
Contains a module which implements a base class for dark matter profiles from which both dark-matter-only and
non-dark-matter-only profiles inherit.
!!}

module Dark_Matter_Profiles_Generic
  !!{
  A base class for dark matter profiles from which both dark-matter-only and non-dark-matter-only profiles inherit. Implements
  numerical calculations of certain halo properties which are to be used as a fall-back option when no analytical solution
  exists.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Function_Classes       , only : functionClass
  use :: Galacticus_Nodes       , only : nodeComponentBasic      , nodeComponentDarkMatterProfile, treeNode
  use :: Kind_Numbers           , only : kind_int8
  use :: Numerical_Interpolation, only : interpolator
  private
  public :: darkMatterProfileGeneric

  !![
  <functionClassType name="darkMatterProfileGeneric"/>
  !!]
  type, extends(functionClass), abstract :: darkMatterProfileGeneric
     !!{
     A dark matter halo profile class implementing numerical calculations for generic dark matter halos.
     !!}
     ! Note that the following components can not be "private", as private components of parent types which are accessed through a
     ! "USE" association are inaccessible to the child type
     ! (e.g. https://www.ibm.com/support/knowledgecenter/SSGH4D_15.1.3/com.ibm.xlf1513.aix.doc/language_ref/extensible.html).
     class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_                         => null()
     ! Tolerances used in numerical solutions.
     double precision                                                      :: toleranceRelativeVelocityDispersion          =  1.0d-6
     double precision                                                      :: toleranceRelativeVelocityDispersionMaximum   =  1.0d-3
     double precision                                                      :: toleranceRelativePotential                   =  1.0d-6
     ! Unique ID for memoization
     integer         (kind_int8               )                            :: genericLastUniqueID
     ! Memoized solutions for the radial velocity dispersion.
     double precision                          , allocatable, dimension(:) :: genericVelocityDispersionRadialVelocity               , genericVelocityDispersionRadialRadius
     double precision                                                      :: genericVelocityDispersionRadialRadiusMinimum          , genericVelocityDispersionRadialRadiusMaximum, &
          &                                                                   genericVelocityDispersionRadialRadiusOuter
     type            (interpolator            ), allocatable               :: genericVelocityDispersionRadial
     ! Memoized solutions for the enclosed mass.
     double precision                          , allocatable, dimension(:) :: genericEnclosedMassMass                               , genericEnclosedMassRadius
     double precision                                                      :: genericEnclosedMassRadiusMinimum                      , genericEnclosedMassRadiusMaximum
     type            (interpolator            ), allocatable               :: genericEnclosedMass
   contains 
     !![
     <methods>
       <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="enclosedMass" />
       <method description="Returns the density (in $M_\odot/$Mpc$^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="density" />
       <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation." method="enclosedMassNumerical" />
       <method description="Returns the enclosed mass difference (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation." method="enclosedMassDifferenceNumerical" />
       <method description="Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation." method="potentialNumerical" />
       <method description="Returns the gravitational potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation." method="potentialDifferenceNumerical" />
       <method description="Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="circularVelocityNumerical" />
       <method description="Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="radialVelocityDispersionNumerical" />
       <method description="Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc)." method="radialMomentNumerical" />
       <method description="Return the normalization of the rotation velocity vs. specific angular momentum relation." method="rotationNormalizationNumerical" />
       <method description="Returns the Fourier transform of the density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$)." method="kSpaceNumerical" />
       <method description="Return the energy of the dark matter density profile." method="energyNumerical" />
       <method description="Returns the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr)." method="freefallRadiusNumerical" />
       <method description="Returns the rate of increase of the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr)." method="freefallRadiusIncreaseRateNumerical" />
       <method description="Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$)." method="radiusEnclosingDensityNumerical" />
       <method description="Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily mass} (given in units of $M_\odot$)." method="radiusEnclosingMassNumerical" />
       <method description="Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}." method="circularVelocityMaximumNumerical" />
       <method description="Returns the radius (in Mpc) at which the maximum circular velocity is reached in the dark matter profile of {\normalfont \ttfamily node}." method="radiusCircularVelocityMaximumNumerical" />
       <method description="Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given in units of km s$^{-1}$ Mpc)." method="radiusFromSpecificAngularMomentumNumerical" />
       <method description="Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="densityLogSlopeNumerical" />
       <method description="Set a sub-module scope pointers on a stack to allow recursive calls to functions." method="solverSet"/>
       <method description="Unset the sub-module scope pointer." method="solverUnset"/>
       <method description="Reset generic profile memoized calculations." method="calculationResetGeneric"/>
       <method description="Integrand for dark matter profile Jeans equation." method="jeansEquationIntegrand"/>
       <method description="Return the radius variable used in solving the Jeans equation that corresponds to a given physical radius." method="jeansEquationRadius"/>
     </methods>
     !!]
     procedure(genericDensityInterface     ), deferred :: density
     procedure(genericEnclosedMassNumerical), deferred :: enclosedMass
     procedure                                         :: enclosedMassNumerical                      => genericEnclosedMassNumerical
     procedure                                         :: enclosedMassDifferenceNumerical            => genericEnclosedMassDifferenceNumerical
     procedure                                         :: potentialNumerical                         => genericPotentialNumerical
     procedure                                         :: potentialDifferenceNumerical               => genericPotentialDifferenceNumerical
     procedure                                         :: circularVelocityNumerical                  => genericCircularVelocityNumerical
     procedure                                         :: radialVelocityDispersionNumerical          => genericRadialVelocityDispersionNumerical
     procedure                                         :: jeansEquationIntegrand                     => genericJeansEquationIntegrand
     procedure                                         :: jeansEquationRadius                        => genericJeansEquationRadius
     procedure                                         :: radialMomentNumerical                      => genericRadialMomentNumerical
     procedure                                         :: rotationNormalizationNumerical             => genericRotationNormalizationNumerical
     procedure                                         :: kSpaceNumerical                            => genericKSpaceNumerical
     procedure                                         :: energyNumerical                            => genericEnergyNumerical
     procedure                                         :: freefallRadiusNumerical                    => genericFreefallRadiusNumerical
     procedure                                         :: freefallRadiusIncreaseRateNumerical        => genericFreefallRadiusIncreaseRateNumerical
     procedure                                         :: radiusEnclosingDensityNumerical            => genericRadiusEnclosingDensityNumerical
     procedure                                         :: radiusEnclosingMassNumerical               => genericRadiusEnclosingMassNumerical
     procedure                                         :: circularVelocityMaximumNumerical           => genericCircularVelocityMaximumNumerical
     procedure                                         :: radiusCircularVelocityMaximumNumerical     => genericRadiusCircularVelocityMaximumNumerical
     procedure                                         :: radiusFromSpecificAngularMomentumNumerical => genericRadiusFromSpecificAngularMomentumNumerical
     procedure                                         :: densityLogSlopeNumerical                   => genericDensityLogSlopeNumerical
     procedure                                         :: solverSet                                  => genericSolverSet
     procedure                              , nopass   :: solverUnset                                => genericSolverUnset
     procedure                                         :: calculationResetGeneric                    => genericCalculationResetGeneric
  end type darkMatterProfileGeneric

  abstract interface
     double precision function genericDensityInterface(self,node,radius)
       !!{
       Returns the density (in $M_\odot/$Mpc$^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
       units of Mpc).
       !!}
       import darkMatterProfileGeneric, treeNode
       class           (darkMatterProfileGeneric), intent(inout) :: self
       type            (treeNode                ), intent(inout) :: node
       double precision                          , intent(in   ) :: radius
     end function genericDensityInterface
  end interface

  ! Module-scope pointers used in integrand functions and root finding.
  type :: genericSolver
     class(darkMatterProfileGeneric), pointer :: self => null()
     type (treeNode                ), pointer :: node => null()
  end type genericSolver
  type            (genericSolver                 ), allocatable, dimension(:) :: solvers
  integer                                         , parameter                 :: solversIncrement              =10
  integer                                                                     :: solversCount                  = 0
  class           (nodeComponentBasic            ), pointer                   :: genericBasic
  class           (nodeComponentDarkMatterProfile), pointer                   :: genericDarkMatterProfile
  double precision                                                            :: genericTime                      , genericRadiusFreefall , genericDensity        , genericMass , &
       &                                                                         genericSpecificAngularMomentum   , genericMassGrowthRate , genericScaleGrowthRate, genericScale, &
       &                                                                         genericShape                     , genericShapeGrowthRate
  !$omp threadprivate(solvers,solversCount,genericBasic,genericTime,genericRadiusFreefall,genericDensity,genericMass,genericSpecificAngularMomentum,genericMassGrowthRate,genericDarkMatterProfile,genericScaleGrowthRate,genericScale,genericShape,genericShapeGrowthRate)

contains

  double precision function genericEnclosedMassNumerical(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    genericEnclosedMassNumerical=self%enclosedMassDifferenceNumerical(node,0.0d0,radius)
    return
  end function genericEnclosedMassNumerical

  double precision function genericEnclosedMassDifferenceNumerical(self,node,radiusLower,radiusUpper)
    !!{
    Returns the enclosed mass difference (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} between the
    given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical
    calculation.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Numerical_Integration  , only : integrator
    use            :: Numerical_Ranges       , only : Make_Range                  , rangeTypeLogarithmic
    use            :: Table_Labels           , only : extrapolationTypeExtrapolate
    use            :: Numerical_Interpolation, only : gsl_interp_linear
    use            :: Error                  , only : Error_Report
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target      :: self
    type            (treeNode                ), intent(inout), target      :: node
    double precision                          , intent(in   )              :: radiusLower                      , radiusUpper
    double precision                                         , parameter   :: countPointsPerOctave     =4.0d+00
    double precision                                         , parameter   :: radiusVirialFractionSmall=1.0d-12
    double precision                          , dimension(:) , allocatable :: masses                           , radii
    type            (integrator              ), save                       :: integrator_
    logical                                   , save                       :: initialized              =.false.
    !$omp threadprivate(integrator_,initialized)
    integer         (c_size_t                )                             :: countRadii                       , iMinimum           , &
         &                                                                    iMaximum                         , i
    logical                                                                :: remakeTable
    double precision                                                       :: radiusIntegralLower              , radiusIntegralUpper, &
         &                                                                    radiusMinimum                    , radiusMaximum      , &
         &                                                                    radiusVirial

    ! Validate input.
    if (radiusUpper < radiusLower) call Error_Report('radiusUpper ≥ radiusLower is required'//{introspection:location})
    if (radiusUpper <= 0.0d0) then
       genericEnclosedMassDifferenceNumerical=0.0d0
       return
    end if
    ! Reset calculations if necessary.
    if (node%uniqueID() /= self%genericLastUniqueID) call self%calculationResetGeneric(node)
    ! Determine if the table must be rebuilt.
    remakeTable=.false.
    if (.not.allocated(self%genericEnclosedMassMass)) then
       remakeTable=.true.
    else
       remakeTable= radiusLower < self%genericEnclosedMassRadiusMinimum &
            &      .or.                                                 &
            &       radiusUpper > self%genericEnclosedMassRadiusMaximum
    end if
    if (remakeTable) then
       ! Initialize integrator if necessary.
    if (.not.initialized) then
       integrator_=integrator(genericMassIntegrand,toleranceRelative=1.0d-2)
       initialized=.true.
    end if
       ! Find the range of radii at which to compute the enclosed mass, and construct the arrays.
       call self%solverSet  (node)
       radiusVirial =self%darkMatterHaloScale_%radiusVirial(node)
       !! Set an initial range of radii that brackets the requested radii.
       if (radiusLower <= 0.0d0) then
          radiusMinimum=max(0.5d0*radiusUpper,radiusVirial*radiusVirialFractionSmall)
       else
          radiusMinimum=max(0.5d0*radiusLower,radiusVirial*radiusVirialFractionSmall)
       end if
       radiusMaximum=2.0d0*radiusUpper
       !! Round to the nearest factor of 2.
       radiusMinimum=2.0d0**floor  (log(radiusMinimum)/log(2.0d0))
       radiusMaximum=2.0d0**ceiling(log(radiusMaximum)/log(2.0d0))
       !! Expand to encompass any pre-existing range.
       if (allocated(self%genericEnclosedMassRadius)) then
          radiusMinimum=min(radiusMinimum,self%genericEnclosedMassRadiusMinimum)
          radiusMaximum=max(radiusMaximum,self%genericEnclosedMassRadiusMaximum)
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
       if (allocated(self%genericEnclosedMassRadius)) then
          iMinimum=nint(log(self%genericEnclosedMassRadiusMinimum/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
          iMaximum=nint(log(self%genericEnclosedMassRadiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
          masses(iMinimum:iMaximum)=self%genericEnclosedMassMass
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
          masses           (i)= integrator_%integrate(radiusIntegralLower,radiusIntegralUpper)
          if (i > 1) masses(i)=+masses(i  ) &
               &               +masses(i-1)
       end do
       call self%solverUnset(   )
       ! Build the interpolator.
       if (allocated(self%genericEnclosedMass)) deallocate(self%genericEnclosedMass)
       allocate(self%genericEnclosedMass)
       self%genericEnclosedMass=interpolator(log(radii),log(masses),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
       ! Store the current results for future re-use.
       if (allocated(self%genericEnclosedMassRadius)) deallocate(self%genericEnclosedMassRadius)
       if (allocated(self%genericEnclosedMassMass  )) deallocate(self%genericEnclosedMassMass  )
       allocate(self%genericEnclosedMassRadius(countRadii))
       allocate(self%genericEnclosedMassMass  (countRadii))
       self%genericEnclosedMassRadius       =radii
       self%genericEnclosedMassMass         =masses
       self%genericEnclosedMassRadiusMinimum=radiusMinimum
       self%genericEnclosedMassRadiusMaximum=radiusMaximum
    end if
    ! Interpolate in the table to find the mass difference.
    genericEnclosedMassDifferenceNumerical       =+exp(self%genericEnclosedMass                   %interpolate(log(radiusUpper)))
    if (radiusLower > 0.0d0)                                                                                                      &
         & genericEnclosedMassDifferenceNumerical=+         genericEnclosedMassDifferenceNumerical                                &
         &                                        +exp(self%genericEnclosedMass                   %interpolate(log(radiusLower)))
    return
  end function genericEnclosedMassDifferenceNumerical
  
  double precision function genericMassIntegrand(radius)
    !!{
    Integrand for mass in generic dark matter profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       genericMassIntegrand=4.0d0*Pi*radius**2*solvers(solversCount)%self%density(solvers(solversCount)%node,radius)
    else
       genericMassIntegrand=0.0d0
    end if
    return
  end function genericMassIntegrand

  double precision function genericPotentialNumerical(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc) using a numerical calculation.
    !!}
    use :: Galactic_Structure_Options      , only : enumerationStructureErrorCodeType, structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration           , only : integrator
    implicit none
    class           (darkMatterProfileGeneric         ), intent(inout), target   :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                   , parameter               :: radiusMaximumFactor=1.0d2
    type            (integrator                       ), save                    :: integrator_
    logical                                            , save                    :: initialized        =.false.
    !$omp threadprivate(integrator_,initialized)
    double precision                                                             :: radiusMaximum

    if (present(status)) status=structureErrorCodeSuccess
    if (.not.initialized) then
       integrator_=integrator(integrandPotential,toleranceRelative=self%toleranceRelativePotential)
       initialized=.true.
    end if
    call self%solverSet  (node)
    radiusMaximum             =  +radiusMaximumFactor                          &
         &                       *self%darkMatterHaloScale_%radiusVirial(node)
    if (radius < radiusMaximum) then
       genericPotentialNumerical =   integrator_%integrate(               &
            &                                              radius       , &
            &                                              radiusMaximum  &
            &                                             )
    else
       ! Beyond some large radius approximate as a point mass.
       genericPotentialNumerical=+gravitationalConstantGalacticus                                                   &
            &                    *solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radiusMaximum) &
            &                    *(                                                                                 &
            &                      +1.0d0/radiusMaximum                                                             &
            &                      -1.0d0/radius                                                                    &
            &                     )
    end if
    call self%solverUnset(   )
    return
  end function genericPotentialNumerical

  double precision function genericPotentialDifferenceNumerical(self,node,radiusLower,radiusUpper)
    !!{
    Returns the potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the
    given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical
    calculation.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target   :: self
    type            (treeNode                ), intent(inout), pointer  :: node
    double precision                          , intent(in   )           :: radiusLower        , radiusUpper
    type            (integrator              ), save                    :: integrator_
    logical                                   , save                    :: initialized=.false.
    !$omp threadprivate(integrator_,initialized)
 
    if (.not.initialized) then
       integrator_=  integrator(integrandPotential,toleranceRelative=1.0d-6)
       initialized=.true.
    end if
    call self%solverSet  (node)
    genericPotentialDifferenceNumerical=integrator_%integrate(             &
         &                                                    radiusLower, &
         &                                                    radiusUpper  &
         &                                                   )
    call self%solverUnset(   )
    return
  end function genericPotentialDifferenceNumerical

  double precision function integrandPotential(radius)
    !!{
    Integrand for gravitational potential in a generic dark matter profile.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius

    if (radius > 0.0d0) then
       integrandPotential=-gravitationalConstantGalacticus                                               &
            &             *solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius)    &
            &             /                                                                   radius **2
    else
       integrandPotential=0.0d0
    end if
    return
  end function integrandPotential

  double precision function genericCircularVelocityNumerical(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: radius

    if (radius > 0.0d0) then
       genericCircularVelocityNumerical=sqrt(                                 &
            &                                +gravitationalConstantGalacticus &
            &                                *self%enclosedMass(node,radius)  &
            &                                /                       radius   &
            &                               )
    else
       genericCircularVelocityNumerical=0.0d0
    end if
    return
  end function genericCircularVelocityNumerical
 
  double precision function genericRadialVelocityDispersionNumerical(self,node,radius,radiusOuter)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Error                  , only : Error_Report        , errorStatusSuccess
    use            :: Numerical_Integration  , only : integrator
    use            :: Numerical_Ranges       , only : Make_Range          , rangeTypeLogarithmic
    use            :: Table_Labels           , only : extrapolationTypeFix
    use            :: Numerical_Interpolation, only : gsl_interp_linear
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target      :: self
    type            (treeNode                ), intent(inout), target      :: node
    double precision                          , intent(in   )              :: radius
    double precision                          , intent(in   ), optional    :: radiusOuter
    double precision                                         , parameter   :: radiusTinyFactor     =1.0d-9 , radiusLargeFactor=5.0d2
    double precision                                         , parameter   :: countPointsPerOctave =2.0d0
    double precision                                         , parameter   :: toleranceFactor      =2.0d0
    double precision                          , dimension(:) , allocatable :: velocityDispersions          , radii
    double precision                                                       :: radiusMinimum                , radiusMaximum           , &
         &                                                                    radiusVirial                 , density                 , &
         &                                                                    jeansIntegral                , radiusOuter_            , &
         &                                                                    radiusLower                  , radiusUpper             , &
         &                                                                    radiusLowerJeansEquation     , radiusUpperJeansEquation, &
         &                                                                    jeansIntegralPrevious        , toleranceRelative
    integer         (c_size_t                )                             :: countRadii                   , iMinimum                , &
         &                                                                    iMaximum                     , i
    integer                                                                :: status
    type            (integrator              ), save                       :: integrator_
    logical                                   , save                       :: initialized          =.false.
    logical                                                                :: remakeTable
    !$omp threadprivate(integrator_,initialized)

    ! Reset calculations if necessary.
    if (node%uniqueID() /= self%genericLastUniqueID) call self%calculationResetGeneric(node)
    ! Determine if the table must be rebuilt.
    remakeTable=.false.
    if (.not.allocated(self%genericVelocityDispersionRadialVelocity)) then
       remakeTable=.true.
    else
       remakeTable= radius < self%genericVelocityDispersionRadialRadiusMinimum &
            &      .or.                                                        &
            &       radius > self%genericVelocityDispersionRadialRadiusMaximum
    end if
    if (remakeTable) then
       ! Initialize integrator if necessary.
       if (.not.initialized) then
          integrator_=integrator(genericJeansEquationIntegrand_,toleranceRelative=self%toleranceRelativeVelocityDispersion)
          initialized=.true.
       end if
       ! Find the range of radii at which to compute the velocity dispersion, and construct the arrays.
       call self%solverSet  (node)
       radiusVirial =self%darkMatterHaloScale_%radiusVirial(node)
       !! Set an initial range of radii that brackets the requested radius, but avoids tiny radii.
       radiusMinimum=max(0.5d0*radius,radiusTinyFactor*radiusVirial)
       radiusMaximum=max(2.0d0*radius,           2.0d0*radiusVirial)
       !! Round to the nearest factor of 2.
       radiusMinimum=2.0d0**floor  (log(radiusMinimum)/log(2.0d0))
       radiusMaximum=2.0d0**ceiling(log(radiusMaximum)/log(2.0d0))
       !! Expand to encompass any pre-existing range.
       if (allocated(self%genericVelocityDispersionRadialRadius)) then
          radiusMinimum=min(radiusMinimum,self%genericVelocityDispersionRadialRadiusMinimum)
          radiusMaximum=max(radiusMaximum,self%genericVelocityDispersionRadialRadiusMaximum)
       end if
       !! Set a suitable outer radius for integration.
       if (present(radiusOuter)) then
          radiusOuter_=radiusOuter
       else
          radiusOuter_=max(10.0d0*radiusMaximum,radiusLargeFactor*radiusVirial)
       end if
       !! Construct arrays.
       countRadii=nint(log(radiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave+1.0d0)
       allocate(radii              (countRadii))
       allocate(velocityDispersions(countRadii))
       radii=Make_Range(radiusMinimum,radiusMaximum,int(countRadii),rangeTypeLogarithmic)
       ! Copy in any usable results from any previous solution.
       !! Assume by default that no previous solutions are usable.
       iMinimum=+huge(0_c_size_t)
       iMaximum=-huge(0_c_size_t)
       !! Check that a pre-existing solution exists.
       if (allocated(self%genericVelocityDispersionRadialRadius)) then
          !! Check that the outer radius for integration has not changed - if it has we need to recompute the full solution for
          !! consistency.
          if (radiusOuter_ == self%genericVelocityDispersionRadialRadiusOuter) then
             iMinimum=nint(log(self%genericVelocityDispersionRadialRadiusMinimum/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             iMaximum=nint(log(self%genericVelocityDispersionRadialRadiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             velocityDispersions(iMinimum:iMaximum)=self%genericVelocityDispersionRadialVelocity
          end if
       end if
       ! Solve for the velocity dispersion where old results were unavailable.
       jeansIntegralPrevious=0.0d0
       do i=countRadii,1,-1
          ! Skip cases for which we have a pre-existing solution.
          if (i >= iMinimum .and. i <= iMaximum) cycle
          ! Find the limits for the integral.
          if (i == countRadii) then
             radiusUpper=radiusOuter_
          else
             radiusUpper=radii(i+1)
          end if
          radiusLower   =radii(i  )
          ! Reset the accumulated Jeans integral if necessary.
          if (i == iMinimum-1) jeansIntegralPrevious=+     velocityDispersions(           iMinimum )**2 &
               &                                     *self%density            (node,radii(iMinimum))
          ! If the interval is wholly outside of the outer radius, the integral is zero.
          if (radiusLower > radiusOuter_) then
             jeansIntegral         =0.0d0
             velocityDispersions(i)=0.0d0
         else
             ! Evaluate the integral.
             density                 =self       %density            (node,radiusLower                                             )
             radiusLowerJeansEquation=self       %jeansEquationRadius(node,radiusLower                                             )
             radiusUpperJeansEquation=self       %jeansEquationRadius(node,radiusUpper                                             )
             jeansIntegral           =integrator_%integrate          (     radiusLowerJeansEquation,radiusUpperJeansEquation,status)
             if (status /= errorStatusSuccess) then
                ! Integration failed.
                toleranceRelative=+     toleranceFactor                     &
                     &            *self%toleranceRelativeVelocityDispersion
                do while (toleranceRelative < self%toleranceRelativeVelocityDispersionMaximum)
                   call integrator_%toleranceSet(toleranceRelative=toleranceRelative)
                  jeansIntegral=integrator_%integrate(radiusLowerJeansEquation,radiusUpperJeansEquation,status)
                  if (status == errorStatusSuccess) then
                      exit
                   else
                      toleranceRelative=+toleranceFactor   &
                           &            *toleranceRelative
                   end if
                end do
                if (status /= errorStatusSuccess) call Error_Report('integration of Jeans equation failed'//{introspection:location})
                call integrator_%toleranceSet(toleranceRelative=self%toleranceRelativeVelocityDispersion)
             end if
             if (density <= 0.0d0) then
                ! Density is zero - the velocity dispersion is undefined. If the Jeans integral is also zero this is acceptable - we've
                ! been asked for the velocity dispersion in a region of zero density, so we simply return zero dispersion as it should have
                ! no consequence. If the Jeans integral is non-zero however, then something has gone wrong.
                velocityDispersions(i)=0.0d0
                if (jeansIntegral+jeansIntegralPrevious > 0.0d0) call Error_Report('undefined velocity dispersion'//{introspection:location})
             else
                velocityDispersions(i)=sqrt(                         &
                     &                      +(                       &
                     &                        +jeansIntegral         &
                     &                        +jeansIntegralPrevious &
                     &                       )                       &
                     &                      /density                 &
                     &                     )
             end if
          end if
          jeansIntegralPrevious=+jeansIntegralPrevious &
               &                +jeansIntegral
       end do
       call self%solverUnset(   )
       ! Build the interpolator.
       if (allocated(self%genericVelocityDispersionRadial)) deallocate(self%genericVelocityDispersionRadial)
       allocate(self%genericVelocityDispersionRadial)
       self%genericVelocityDispersionRadial=interpolator(log(radii),velocityDispersions,interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeFix)
       ! Store the current results for future re-use.
       if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
       if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
       allocate(self%genericVelocityDispersionRadialRadius  (countRadii))
       allocate(self%genericVelocityDispersionRadialVelocity(countRadii))
       self%genericVelocityDispersionRadialRadius       =radii
       self%genericVelocityDispersionRadialVelocity     =velocityDispersions
       self%genericVelocityDispersionRadialRadiusMinimum=radiusMinimum
       self%genericVelocityDispersionRadialRadiusMaximum=radiusMaximum
       self%genericVelocityDispersionRadialRadiusOuter  =radiusOuter_
    end if
    ! Interpolate in the table to find the velocity dispersion.
    genericRadialVelocityDispersionNumerical=self%genericVelocityDispersionRadial%interpolate(log(radius))
    return
  end function genericRadialVelocityDispersionNumerical
  
  double precision function genericJeansEquationIntegrand_(radius)
    !!{
    Integrand for generic dark matter profile Jeans equation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    genericJeansEquationIntegrand_=solvers(solversCount)%self%jeansEquationIntegrand(solvers(solversCount)%node,radius)
    return
  end function genericJeansEquationIntegrand_

  double precision function genericJeansEquationIntegrand(self,node,radius)
    !!{
    Integrand for generic dark matter profile Jeans equation.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: radius

    if (radius > 0.0d0) then
       genericJeansEquationIntegrand=+gravitationalConstantGalacticus   &
            &                        *self%enclosedMass(node,radius)    &
            &                        *self%density     (node,radius)    &
            &                        /                       radius **2
    else
       genericJeansEquationIntegrand=0.0d0
    end if
    return
  end function genericJeansEquationIntegrand

  double precision function genericJeansEquationRadius(self,node,radius)
    !!{
    Return the radius variable used in solving the Jeans equation that corresponds to a given physical radius.
    In some cases, it is easier to do the integration with respect to another variable which is a function of
    the physical radius.
    !!}
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: radius
    !$GLC attributes unused :: self, node

    genericJeansEquationRadius=radius
    return
  end function genericJeansEquationRadius
  
  double precision function genericRadialMomentNumerical(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given
    {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout)           :: self
    type            (treeNode                ), intent(inout)           :: node
    double precision                          , intent(in   )           :: moment
    double precision                          , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    type            (integrator              )                          :: integrator_
    double precision                                                    :: radiusMinimumActual, radiusMaximumActual

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%radiusVirial(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    integrator_=integrator(integrandRadialMoment,toleranceRelative=1.0d-3)
    genericRadialMomentNumerical=integrator_%integrate(radiusMinimumActual,radiusMaximumActual)
    return

  contains

    double precision function integrandRadialMoment(radius)
      !!{
      Integrand for radial moment in a generic dark matter profile.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandRadialMoment=+                  radius **moment &
              &                *self%density(node,radius)
      else
         integrandRadialMoment=0.0d0
      end if
      return
    end function integrandRadialMoment

  end function genericRadialMomentNumerical

  double precision function genericRotationNormalizationNumerical(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                                          :: radiusVirial

    radiusVirial                         =+self%darkMatterHaloScale_%radiusVirial         (                            &
         &                                                                                               node          &
         &                                                                                )
    genericRotationNormalizationNumerical=+self                     %enclosedMass         (                            &
         &                                                                                               node        , &
         &                                                                                               radiusVirial  &
         &                                                                                 )                           &
         &                                /4.0d0                                                                       &
         &                                /Pi                                                                          &
         &                                /self                     %radialMomentNumerical(                            &
         &                                                                                               node        , &
         &                                                                                 moment       =3.0d0       , &
         &                                                                                 radiusMinimum=0.0d0       , &
         &                                                                                 radiusMaximum=radiusVirial  &
         &                                                                                )
    return
  end function genericRotationNormalizationNumerical

  double precision function genericKSpaceNumerical(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the dark matter density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout)         :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: waveNumber
    type            (integrator              )                        :: integrator_
    double precision                                                  :: radiusVirial

    radiusVirial          =+self       %darkMatterHaloScale_%radiusVirial(node                                              )
    integrator_           = integrator                                   (integrandFourierTransform,toleranceRelative=1.0d-3)
    genericKSpaceNumerical=+integrator_%integrate                        (0.0d0                    ,radiusVirial            ) &
         &                 /self                            %enclosedMass(node                     ,radiusVirial            )
    return

  contains

    double precision function integrandFourierTransform(radius)
      !!{
      Integrand for Fourier transform of the generic dark matter profile.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandFourierTransform=+4.0d0                     &
              &                    *Pi                        &
              &                    *               radius **2 &
              &                    *sin(wavenumber*radius)    &
              &                    /   (waveNumber*radius)    &
              &                    *self%density(node,radius)
      else
         integrandFourierTransform=0.0d0
      end if
      return
    end function integrandFourierTransform

  end function genericKSpaceNumerical

  double precision function genericEnergyNumerical(self,node)
    !!{
    Return the energy of a generic dark matter density profile.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration           , only : integrator
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , parameter     :: multiplierRadius   =100.0d0
    type            (integrator                )                :: integratorPotential        , integratorKinetic, &
         &                                                         integratorPressure
    double precision                                            :: radiusVirial               , radiusLarge      , &
         &                                                         energyPotential            , energyKinetic    , &
         &                                                         pseudoPressure
    
    integratorPotential=integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
    integratorKinetic  =integrator(integrandEnergyKinetic  ,toleranceRelative=1.0d-3)
    integratorPressure =integrator(integrandPseudoPressure ,toleranceRelative=1.0d-3)
    radiusVirial          =+self%darkMatterHaloScale_%radiusVirial(node)
    radiusLarge           =+multiplierRadius                                          &
         &                 *radiusVirial
    energyPotential       =+integratorPotential%integrate(0.0d0       ,radiusVirial)
    energyKinetic         =+integratorKinetic  %integrate(0.0d0       ,radiusVirial)
    pseudoPressure        =+integratorPressure %integrate(radiusVirial,radiusLarge )
    genericEnergyNumerical=-0.5d0                                                     &
         &                 *gravitationalConstantGalacticus                           &
         &                 *(                                                         &
         &                   +energyPotential                                         &
         &                   +self%enclosedMass(node,radiusVirial)**2                 &
         &                   /                       radiusVirial                     &
         &                  )                                                         &
         &                 +2.0d0                                                     &
         &                 *Pi                                                        &
         &                 *gravitationalConstantGalacticus                           &
         &                 *(                                                         &
         &                   +radiusVirial**3                                         &
         &                   *pseudoPressure                                          &
         &                   +energyKinetic                                           &
         &                  )
    return

  contains

    double precision function integrandEnergyPotential(radius)
      !!{
      Integrand for potential energy of the halo.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandEnergyPotential=(                                &
              &                    +self%enclosedMass(node,radius) &
              &                    /                       radius  &
              &                   )**2
      else
         integrandEnergyPotential=0.0d0
      end if
      return
    end function integrandEnergyPotential

    double precision function integrandEnergyKinetic(radius)
      !!{
      Integrand for kinetic energy of the halo.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandEnergyKinetic=+self%enclosedMass(node,radius) &
              &                 *self%density     (node,radius) &
              &                 *                       radius
      else
         integrandEnergyKinetic=0.0d0
      end if
      return
    end function integrandEnergyKinetic

    double precision function integrandPseudoPressure(radius)
      !!{
      Integrand for pseudo-pressure ($\rho(r) \sigma^2(r)$) of the halo.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         integrandPseudoPressure=+self%enclosedMass(node,radius)    &
              &                  *self%density     (node,radius)    &
              &                  /                       radius **2
      else
         integrandPseudoPressure=0.0d0
      end if
      return
    end function integrandPseudoPressure

  end function genericEnergyNumerical

  double precision function genericEnergyEvaluate(timeLogarithmic)
    !!{
    GSL-callable function to evaluate the energy of the dark matter profile.
    !!}
    use :: Functions_Global, only : Calculations_Reset_
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic
    double precision                       :: time

    time=exp(timeLogarithmic)
    call genericBasic            %timeSet            (                                     time             )
    call genericBasic            %timeLastIsolatedSet(                                     time             )
    call genericBasic            %massSet            (genericMass +genericMassGrowthRate *(time-genericTime))
    call genericDarkMatterProfile%scaleSet           (genericScale+genericScaleGrowthRate*(time-genericTime))
    call genericDarkMatterProfile%shapeSet           (genericShape+genericShapeGrowthRate*(time-genericTime))
    call Calculations_Reset_(solvers(solversCount)%node)
    genericEnergyEvaluate=solvers(solversCount)%self%energyNumerical(solvers(solversCount)%node)
    return
  end function genericEnergyEvaluate

  double precision function genericFreefallRadiusNumerical(self,node,time)
    !!{
    Returns the freefall radius in the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: time
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    call self%solverSet  (node)
    genericTime =  time
    finder      =  rootFinder(                                                             &
         &                    rootFunction                 =rootRadiusFreefall           , &
         &                    toleranceAbsolute            =toleranceAbsolute            , &
         &                    toleranceRelative            =toleranceRelative            , &
         &                    rangeExpandDownward          =0.5d0                        , &
         &                    rangeExpandUpward            =2.0d0                        , &
         &                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                   )
    genericFreefallRadiusNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    call self%solverUnset(   )
    return
  end function genericFreefallRadiusNumerical

  double precision function rootRadiusFreefall(radiusFreefall)
    !!{
    Root function used in finding the radius corresponding to a given freefall time.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radiusFreefall
    type            (integrator)                :: integrator_

    genericRadiusFreefall=+radiusFreefall
    integrator_          = integrator              (integrandTimeFreefall,toleranceRelative=1.0d-3)
    rootRadiusFreefall   =+integrator_   %integrate(0.0d0                ,radiusFreefall          ) &
         &                -genericTime
    return
  end function rootRadiusFreefall

  double precision function integrandTimeFreefall(radius)
    !!{
    Integrand for freefall time in the halo.
    !!}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: potentialDifference

    potentialDifference=+solvers(solversCount)%self%potentialDifferenceNumerical(solvers(solversCount)%node,radius,genericRadiusFreefall)
    if (potentialDifference < 0.0d0) then
       integrandTimeFreefall=+Mpc_per_km_per_s_To_Gyr   &
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

  double precision function genericFreefallRadiusIncreaseRateNumerical(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the dark matter density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: time
    double precision                          , parameter             :: timeLogarithmicStep=0.1d0
    type            (differentiator          )                        :: differentiator_

    call self%solverSet  (node)
    differentiator_                            =   differentiator            (genericFreefallRadiusEvaluate                    )
    genericFreefallRadiusIncreaseRateNumerical =  +differentiator_%derivative(log(time)                    ,timeLogarithmicStep) &
         &                                        /                               time
    call self%solverUnset(   )
    return
  end function genericFreefallRadiusIncreaseRateNumerical

  double precision function genericFreefallRadiusEvaluate(timeLogarithmic)
    !!{
    GSL-callable function to evaluate the freefall radius of the dark matter profile.
    !!}
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    genericFreefallRadiusEvaluate=solvers(solversCount)%self%freefallRadiusNumerical(solvers(solversCount)%node,exp(timeLogarithmic))
    return
  end function genericFreefallRadiusEvaluate

  double precision function genericRadiusEnclosingDensityNumerical(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: density
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    call self%solverSet  (node)
    genericDensity=density
    finder        =rootFinder(                                                             &
         &                    rootFunction                 =rootDensity                  , &
         &                    toleranceAbsolute            =toleranceAbsolute            , &
         &                    toleranceRelative            =toleranceRelative            , &
         &                    rangeExpandDownward          =0.5d0                        , &
         &                    rangeExpandUpward            =2.0d0                        , &
         &                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                   )
    genericRadiusEnclosingDensityNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    call self%solverUnset(   )
    return
  end function genericRadiusEnclosingDensityNumerical

  double precision function rootDensity(radius)
    !!{
    Root function used in finding the radius enclosing a given mean density.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                                                         &
         &      *solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius)    &
         &      /4.0d0                                                                         &
         &      /Pi                                                                            &
         &      /                                                                   radius **3 &
         &      -genericDensity
    return
  end function rootDensity

  double precision function genericRadiusEnclosingMassNumerical(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: mass
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    call self%solverSet  (node)
    genericMass =  mass
    finder      =  rootFinder(                                                             &
         &                    rootFunction                 =rootMass                     , &
         &                    toleranceAbsolute            =toleranceAbsolute            , &
         &                    toleranceRelative            =toleranceRelative            , &
         &                    rangeExpandUpward            =2.0d0                        , &
         &                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                   )
    genericRadiusEnclosingMassNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%radiusVirial(node)])
    call self%solverUnset(   )
    return
  end function genericRadiusEnclosingMassNumerical

  double precision function rootMass(radius)
    !!{
    Root function used in finding the radius enclosing a given mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius) &
         &   -             genericMass
    return
  end function rootMass

  double precision function genericRadiusCircularVelocityMaximumNumerical(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    use :: Root_Finder         , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    class           (nodeComponentBasic      ), pointer       :: basic
    type            (rootFinder              )                :: finder

    call self%solverSet  (node)
    finder      =  rootFinder(                                                             &
         &                    rootFunction                 =rootCircularVelocityMaximum  , &
         &                    toleranceAbsolute            =toleranceAbsolute            , &
         &                    toleranceRelative            =toleranceRelative            , &
         &                    rangeExpandDownward          =0.5d0                        , &
         &                    rangeExpandUpward            =2.0d0                        , &
         &                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                   )
    ! Isothermal profiles have dVc²/dr=0 everywhere. To handle these profiles, first test if the root function is sufficiently
    ! close to zero at a few points throughout the halo (which it will be for an isothermal profile), and return the circular
    ! velocity at the virial radius if so. Otherwise solve for the radius corresponding to the maximum circular velocity.
    basic => node%basic()
    if     (                                                                                                        &
         &   Values_Agree(                                                                                          &
         &                       +rootCircularVelocityMaximum(1.0d+0*self%darkMatterHaloScale_%radiusVirial(node)), &
         &                       +0.0d0                                                                           , &
         &                absTol=+toleranceRelative                                                                 &
         &                       *basic%mass                 (                                                   )  &
         &                )                                                                                         &
         &  .and.                                                                                                   &
         &   Values_Agree(                                                                                          &
         &                       +rootCircularVelocityMaximum(3.0d-1*self%darkMatterHaloScale_%radiusVirial(node)), &
         &                       +0.0d0                                                                           , &
         &                absTol=+toleranceRelative                                                                 &
         &                       *basic%mass                 (                                                   )  &
         &                )                                                                                         &
         &  .and.                                                                                                   &
         &   Values_Agree(                                                                                          &
         &                       +rootCircularVelocityMaximum(1.0d-1*self%darkMatterHaloScale_%radiusVirial(node)), &
         &                       +0.0d0                                                                           , &
         &                absTol=+toleranceRelative                                                                 &
         &                       *basic%mass                 (                                                   )  &
         &                )                                                                                         &
         &  .and.                                                                                                   &
         &   Values_Agree(                                                                                          &
         &                       +rootCircularVelocityMaximum(3.0d-2*self%darkMatterHaloScale_%radiusVirial(node)), &
         &                       +0.0d0                                                                           , &
         &                absTol=+toleranceRelative                                                                 &
         &                       *basic%mass                 (                                                   )  &
         &                )                                                                                         &
         & ) then
       genericRadiusCircularVelocityMaximumNumerical=                      self%darkMatterHaloScale_%radiusVirial(node)
    else
       genericRadiusCircularVelocityMaximumNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    end if
    call self%solverUnset(   )
    return
  end function genericRadiusCircularVelocityMaximumNumerical

  double precision function genericCircularVelocityMaximumNumerical(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileGeneric), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node

    genericCircularVelocityMaximumNumerical=self%circularVelocityNumerical(node,self%radiusCircularVelocityMaximumNumerical(node))
    return
  end function genericCircularVelocityMaximumNumerical

  double precision function rootCircularVelocityMaximum(radius)
    !!{
    Root function used in finding the radius at which the maximum circular velocity occurs. Since for a spherical profile $V_\mathrm{c}^2(r)=\mathrm{G}M(r)/r$, then
    \begin{equation}
    {\mathrm{d} V_\mathrm{c}^2 \over \mathrm{d} r} = - {\mathrm{G} M(r) \over r^2} + 4 \pi \mathrm{G} \rho(r) r.
    \end{equation}
    Therefore, the peak of the rotation curve satisfies $4 \pi r^3 \rho(r) - M(r)=0$.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootCircularVelocityMaximum=+4.0d0                                                                         &
         &                      *Pi                                                                            &
         &                      *                                                                   radius **3 &
         &                      *solvers(solversCount)%self%density     (solvers(solversCount)%node,radius)    &
         &                      -solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius)
    return
  end function rootCircularVelocityMaximum

  double precision function genericRadiusFromSpecificAngularMomentumNumerical(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: specificAngularMomentum
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    call self%solverSet  (node)
    genericSpecificAngularMomentum =  specificAngularMomentum
    finder                         =  rootFinder(                                                             &
         &                                       rootFunction                 =rootSpecificAngularMomentum  , &
         &                                       toleranceAbsolute            =toleranceAbsolute            , &
         &                                       toleranceRelative            =toleranceRelative            , &
         &                                       rangeExpandUpward            =2.0d0                        , &
         &                                       rangeExpandType              =rangeExpandMultiplicative    , &
         &                                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                                      )
    genericRadiusFromSpecificAngularMomentumNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%radiusVirial(node)])
    call self%solverUnset(   )
    return
  end function genericRadiusFromSpecificAngularMomentumNumerical

  double precision function rootSpecificAngularMomentum(radius)
    !!{
    Root function used in finding the radius enclosing a given specific angular momentum.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rootSpecificAngularMomentum=+solvers(solversCount)%self%circularVelocityNumerical(solvers(solversCount)%node,radius) &
         &                      *                                                  radius &
         &                      -genericSpecificAngularMomentum
    return
  end function rootSpecificAngularMomentum

  double precision function genericDensityLogSlopeNumerical(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: radius
    double precision                          , parameter             :: radiusLogarithmicStep=0.1d0
    type            (differentiator          )                        :: differentiator_

    call self%solverSet  (node)
    differentiator_                =    differentiator                   (genericDensityEvaluate                      )
    genericDensityLogSlopeNumerical=   +differentiator_       %derivative(log(radius)           ,radiusLogarithmicStep) &
         &                             /genericDensityEvaluate           (log(radius)                                 )
    call self%solverUnset(   )
    return
  end function genericDensityLogSlopeNumerical

  double precision function genericDensityEvaluate(radiusLogarithmic)
    !!{
    GSL-callable function to evaluate the density of the dark matter profile.
    !!}
    implicit none
    double precision, intent(in   ), value :: radiusLogarithmic

    genericDensityEvaluate=solvers(solversCount)%self%density(solvers(solversCount)%node,exp(radiusLogarithmic))
    return
  end function genericDensityEvaluate

  subroutine genericSolverSet(self,node)
    !!{
    Set a sub-module scope pointers on a stack to allow recursive calls to functions.
    !!}
    implicit none
    class  (darkMatterProfileGeneric), intent(inout), target       :: self
    type   (treeNode                ), intent(inout), target       :: node
    type   (genericSolver           ), allocatable  , dimension(:) :: solvers_
    integer                                                        :: i

    ! Increment the state counter. This is necessary to ensure that this function can be called recursively.
    if (allocated(solvers)) then
       if (solversCount == size(solvers)) then
          call move_alloc(solvers,solvers_)
          allocate(solvers(size(solvers_)+solversIncrement))
          solvers(1:size(solvers_))=solvers_
          do i=1,size(solvers_)
             nullify(solvers_(i)%self)
             nullify(solvers_(i)%node)
          end do
          deallocate(solvers_)
       end if
    else
       allocate(solvers(solversIncrement))
    end if
    solversCount=solversCount+1
    solvers(solversCount)%self => self
    solvers(solversCount)%node => node
    return
  end subroutine genericSolverSet
  
  subroutine genericSolverUnset()
    !!{
    Unset a sub-module scope pointers on the stack.
    !!}
    implicit none

    solvers(solversCount)%self => null()
    solvers(solversCount)%node => null()
    solversCount=solversCount-1
    return
  end subroutine genericSolverUnset

  subroutine genericCalculationResetGeneric(self,node)
    !!{
    Reset generic profile memoized data.
    !!}
    implicit none
    class(darkMatterProfileGeneric), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node

    self%genericLastUniqueID=node%uniqueID()
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine genericCalculationResetGeneric
  
end module Dark_Matter_Profiles_Generic
