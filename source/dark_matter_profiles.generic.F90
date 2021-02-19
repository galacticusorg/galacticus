!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a base class for dark matter profiles from which both dark-matter-only and
!% non-dark-matter-only profiles inherit.

module Dark_Matter_Profiles_Generic
  !% A base class for dark matter profiles from which both dark-matter-only and non-dark-matter-only profiles inherit. Implements
  !% numerical calculations of certain halo properties which are to be used as a fall-back option when no analytical solution
  !% exists.
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Function_Classes       , only : functionClass
  use :: Galacticus_Nodes       , only : nodeComponentBasic      , nodeComponentDarkMatterProfile, treeNode
  private
  public :: darkMatterProfileGeneric

  !# <functionClassType name="darkMatterProfileGeneric"/>
  type, extends(functionClass), abstract :: darkMatterProfileGeneric
     !% A dark matter halo profile class implementing numerical calculations for generic dark matter halos.
     ! Note that the following components can not be "private", as private components of parent types which are accessed through a
     ! "USE" association are inaccessible to the child type
     ! (e.g. https://www.ibm.com/support/knowledgecenter/SSGH4D_15.1.3/com.ibm.xlf1513.aix.doc/language_ref/extensible.html).
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_                => null()
     ! Tolerances used in numerical solutions.
     double precision                                    :: toleranceRelativeVelocityDispersion =  1.0d-6
   contains 
     !# <methods>
     !#   <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="enclosedMass" />
     !#   <method description="Returns the density (in $M_\odot/$Mpc$^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="density" />
     !#   <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation." method="enclosedMassNumerical" />
     !#   <method description="Returns the enclosed mass difference (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation." method="enclosedMassDifferenceNumerical" />
     !#   <method description="Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation." method="potentialNumerical" />
     !#   <method description="Returns the gravitational potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation." method="potentialDifferenceNumerical" />
     !#   <method description="Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="circularVelocityNumerical" />
     !#   <method description="Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="radialVelocityDispersionNumerical" />
     !#   <method description="Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc)." method="radialMomentNumerical" />
     !#   <method description="Return the normalization of the rotation velocity vs. specific angular momentum relation." method="rotationNormalizationNumerical" />
     !#   <method description="Returns the Fourier transform of the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$)." method="kSpaceNumerical" />
     !#   <method description="Return the energy of the dark matter density profile." method="energyNumerical" />
     !#   <method description="Return the rate of growth of the energy of the dark matter density profile." method="energyGrowthRateNumerical" />
     !#   <method description="Returns the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr)." method="freefallRadiusNumerical" />
     !#   <method description="Returns the rate of increase of the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr)." method="freefallRadiusIncreaseRateNumerical" />
     !#   <method description="Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$)." method="radiusEnclosingDensityNumerical" />
     !#   <method description="Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily mass} (given in units of $M_\odot$)." method="radiusEnclosingMassNumerical" />
     !#   <method description="Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}." method="circularVelocityMaximumNumerical" />
     !#   <method description="Returns the radius (in Mpc) at which the maximum circular velocity is reached in the dark matter profile of {\normalfont \ttfamily node}." method="radiusCircularVelocityMaximumNumerical" />
     !#   <method description="Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given in units of km s$^{-1}$ Mpc)." method="radiusFromSpecificAngularMomentumNumerical" />
     !#   <method description="Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="densityLogSlopeNumerical" />
     !# </methods>
     procedure(genericDensityInterface     ), deferred :: density
     procedure(genericEnclosedMassNumerical), deferred :: enclosedMass
     procedure                                         :: enclosedMassNumerical                      => genericEnclosedMassNumerical
     procedure                                         :: enclosedMassDifferenceNumerical            => genericEnclosedMassDifferenceNumerical
     procedure                                         :: potentialNumerical                         => genericPotentialNumerical
     procedure                                         :: potentialDifferenceNumerical               => genericPotentialDifferenceNumerical
     procedure                                         :: circularVelocityNumerical                  => genericCircularVelocityNumerical
     procedure                                         :: radialVelocityDispersionNumerical          => genericRadialVelocityDispersionNumerical
     procedure                                         :: radialMomentNumerical                      => genericRadialMomentNumerical
     procedure                                         :: rotationNormalizationNumerical             => genericRotationNormalizationNumerical
     procedure                                         :: kSpaceNumerical                            => genericKSpaceNumerical
     procedure                                         :: energyNumerical                            => genericEnergyNumerical
     procedure                                         :: energyGrowthRateNumerical                  => genericEnergyGrowthRateNumerical
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
  end type darkMatterProfileGeneric

  abstract interface
     double precision function genericDensityInterface(self,node,radius)
       !% Returns the density (in $M_\odot/$Mpc$^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
       !% units of Mpc).
       import darkMatterProfileGeneric, treeNode
       class           (darkMatterProfileGeneric), intent(inout) :: self
       type            (treeNode                ), intent(inout) :: node
       double precision                          , intent(in   ) :: radius
     end function genericDensityInterface
  end interface

  ! Module-scope pointers used in integrand functions and root finding.
  type :: genericSolver
     class(darkMatterProfileGeneric), pointer :: self
     type (treeNode                ), pointer :: node
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

  !# <enumeration>
  !#  <name>nonAnalyticSolvers</name>
  !#  <description>Used to specify the type of solution to use when no analytic solution is available.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <entry label="fallThrough"/>
  !#  <entry label="numerical"  />
  !# </enumeration>

contains

  double precision function genericEnclosedMassNumerical(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    genericEnclosedMassNumerical=self%enclosedMassDifferenceNumerical(node,0.0d0,radius)
    return
  end function genericEnclosedMassNumerical

  double precision function genericEnclosedMassDifferenceNumerical(self,node,radiusLower,radiusUpper)
    !% Returns the enclosed mass difference (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} between the
    !% given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical
    !% calculation.
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: radiusLower        , radiusUpper
    type            (integrator              ), save                  :: integrator_
    logical                                   , save                  :: initialized=.false.
    !$omp threadprivate(integrator_,initialized)
    
    if (.not.initialized) then
       integrator_=integrator(genericMassIntegrand,toleranceRelative=1.0d-2)
       initialized=.true.
    end if
    call self%solverSet  (node)
    genericEnclosedMassDifferenceNumerical =  integrator_%integrate(radiusLower,radiusUpper)
    call self%solverUnset(   )
    return
  end function genericEnclosedMassDifferenceNumerical
  
  double precision function genericMassIntegrand(radius)
    !% Integrand for mass in generic dark matter profiles.
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
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc) using a numerical calculation.
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target   :: self
    type            (treeNode                ), intent(inout), target   :: node
    double precision                          , intent(in   )           :: radius
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: radiusMaximumFactor=1.0d2
    type            (integrator              ), save                    :: integrator_
    logical                                   , save                    :: initialized        =.false.
    !$omp threadprivate(integrator_,initialized)
    double precision                                                    :: radiusMaximum

    if (present(status)) status=structureErrorCodeSuccess
    if (.not.initialized) then
       integrator_=integrator(integrandPotential,toleranceRelative=1.0d-6)
       initialized=.true.
    end if
    call self%solverSet  (node)
    radiusMaximum             =  +radiusMaximumFactor                          &
         &                       *self%darkMatterHaloScale_%virialRadius(node)
    genericPotentialNumerical =   integrator_%integrate(               &
         &                                              radius       , &
         &                                              radiusMaximum  &
         &                                             )
    call self%solverUnset(   )
    return
  end function genericPotentialNumerical

  double precision function genericPotentialDifferenceNumerical(self,node,radiusLower,radiusUpper)
    !% Returns the potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the
    !% given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical
    !% calculation.
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
    !% Integrand for gravitational potential in a generic dark matter profile.
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
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
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

  double precision function genericRadialVelocityDispersionNumerical(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target    :: self
    type            (treeNode                ), intent(inout), target    :: node
    double precision                          , intent(in   )            :: radius
    double precision                                         , parameter :: radiusTinyFraction=1.0d-9 , radiusLargeFactor=5.0d2
    double precision                                                     :: radiusMinimum             , radiusMaximum          , &
         &                                                                  radiusVirial              , density                , &
         &                                                                  jeansIntegral
    type            (integrator              ), save                     :: integrator_
    logical                                   , save                     :: initialized       =.false.
    !$omp threadprivate(integrator_,initialized)
    
    if (.not.initialized) then
       integrator_=integrator(genericJeansEquationIntegrand,toleranceRelative=self%toleranceRelativeVelocityDispersion)
       initialized=.true.
    end if
    call self%solverSet  (node)
    radiusVirial  =  self       %darkMatterHaloScale_%virialRadius(node                            )
    radiusMinimum =  max(       radius,radiusTinyFraction*radiusVirial)
    radiusMaximum =  max(10.0d0*radius,radiusLargeFactor *radiusVirial)    
    jeansIntegral =  integrator_                     %integrate   (     radiusMinimum,radiusMaximum)
    density       =  self                            %density     (node,radiusMinimum              )
    if (density <= 0.0d0) then
       ! Density is zero - the velocity dispersion is undefined. If the Jeans integral is also zero this is acceptable - we've
       ! been asked for the velocity dispersion in a region of zero density, so we simply return zero dispersion as it should have
       ! no consequence. If the Jeans integral in non-zero however, then something has gone wrong.
       genericRadialVelocityDispersionNumerical=0.0d0
       if (jeansIntegral > 0.0d0) call Galacticus_Error_Report('undefined velocity dispersion'//{introspection:location})
    else
       genericRadialVelocityDispersionNumerical=sqrt(               &
            &                                        +jeansIntegral &
            &                                        /density       &
            &                                       )
    end if
    call self%solverUnset(   )
    return
  end function genericRadialVelocityDispersionNumerical
  
  double precision function genericJeansEquationIntegrand(radius)
    !% Integrand for generic drak matter profile Jeans equation.
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       genericJeansEquationIntegrand=+gravitationalConstantGalacticus                                            &
            &                        *solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius) &
            &                        *solvers(solversCount)%self%density     (solvers(solversCount)%node,radius) &
            &                        /radius**2
    else
       genericJeansEquationIntegrand=0.0d0
    end if
    return
  end function genericJeansEquationIntegrand
  
  double precision function genericRadialMomentNumerical(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given
    !% {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout)           :: self
    type            (treeNode                ), intent(inout)           :: node
    double precision                          , intent(in   )           :: moment
    double precision                          , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    type            (integrator              )                          :: integrator_
    double precision                                                    :: radiusMinimumActual, radiusMaximumActual

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%virialRadius(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    integrator_=integrator(integrandRadialMoment,toleranceRelative=1.0d-3)
    genericRadialMomentNumerical=integrator_%integrate(radiusMinimumActual,radiusMaximumActual)
    return

  contains

    double precision function integrandRadialMoment(radius)
      !% Integrand for radial moment in a generic dark matter profile.
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
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                                          :: radiusVirial

    radiusVirial                         =+self%darkMatterHaloScale_%virialRadius         (                            &
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
    !% Returns the Fourier transform of the dark matter density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout)         :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: waveNumber
    type            (integrator              )                        :: integrator_
    double precision                                                  :: radiusVirial

    radiusVirial          =+self       %darkMatterHaloScale_%virialRadius(node                                              )
    integrator_           = integrator                                   (integrandFourierTransform,toleranceRelative=1.0d-3)
    genericKSpaceNumerical=+integrator_%integrate                        (0.0d0                    ,radiusVirial            ) &
         &                 /self                            %enclosedMass(node                     ,radiusVirial            )
    return

  contains

    double precision function integrandFourierTransform(radius)
      !% Integrand for Fourier transform of the generic dark matter profile.
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
    !% Return the energy of a generic dark matter density profile.
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
    radiusVirial          =+self%darkMatterHaloScale_%virialRadius(node)
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
      !% Integrand for potential energy of the halo.
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
      !% Integrand for kinetic energy of the halo.
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
      !% Integrand for pseudo-pressure ($\rho(r) \sigma^2(r)$) of the halo.
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

  double precision function genericEnergyGrowthRateNumerical(self,node)
    !% Returns the rate of growth of the energy if the dark matter density profile.
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , parameter             :: timeLogarithmicStep=0.1d0
    type            (differentiator          )                        :: differentiator_
    double precision                                                  :: timeLastIsolated

    call self%solverSet  (node)
    differentiator_                  =   differentiator                            (genericEnergyEvaluate                    )
    genericBasic                     =>  node                    %basic            (                                         )
    genericDarkMatterProfile         =>  node                    %darkMatterProfile(                                         )
    genericTime                      =   genericBasic            %time             (                                         )
    timeLastIsolated                 =   genericBasic            %timeLastIsolated (                                         )
    genericMass                      =   genericBasic            %mass             (                                         )
    genericMassGrowthRate            =   genericBasic            %accretionRate    (                                         )
    genericScale                     =   genericDarkMatterProfile%scale            (                                         )
    genericScaleGrowthRate           =   genericDarkMatterProfile%scaleGrowthRate  (                                         )
    genericShape                     =   genericDarkMatterProfile%shape            (                                         )
    genericShapeGrowthRate           =   genericDarkMatterProfile%shapeGrowthRate  (                                         )
    genericEnergyGrowthRateNumerical =  +differentiator_         %derivative       (log(genericTime)     ,timeLogarithmicStep) &
         &                              /                                               genericTime
    call genericBasic%timeSet                       (genericTime           )
    call genericBasic%timeLastIsolatedSet           (timeLastIsolated      )
    call genericBasic%massSet                       (genericMass           )
    call genericDarkMatterProfile%scaleSet          (genericScale          )
    call genericDarkMatterProfile%scaleGrowthRateSet(genericScaleGrowthRate)
    call genericDarkMatterProfile%shapeSet          (genericShape          )
    call genericDarkMatterProfile%shapeGrowthRateSet(genericShapeGrowthRate)
    call self%solverUnset(   )
    return
  end function genericEnergyGrowthRateNumerical

  double precision function genericEnergyEvaluate(timeLogarithmic)
    !% GSL-callable function to evaluate the energy of the dark matter profile.
    use :: Functions_Global, only : Galacticus_Calculations_Reset_
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic
    double precision                       :: time

    time=exp(timeLogarithmic)
    call genericBasic            %timeSet            (                                     time             )
    call genericBasic            %timeLastIsolatedSet(                                     time             )
    call genericBasic            %massSet            (genericMass +genericMassGrowthRate *(time-genericTime))
    call genericDarkMatterProfile%scaleSet           (genericScale+genericScaleGrowthRate*(time-genericTime))
    call genericDarkMatterProfile%shapeSet           (genericShape+genericShapeGrowthRate*(time-genericTime))
    call Galacticus_Calculations_Reset_(solvers(solversCount)%node)
    genericEnergyEvaluate=solvers(solversCount)%self%energyNumerical(solvers(solversCount)%node)
    return
  end function genericEnergyEvaluate

  double precision function genericFreefallRadiusNumerical(self,node,time)
    !% Returns the freefall radius in the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
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
    genericFreefallRadiusNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node))
    call self%solverUnset(   )
    return
  end function genericFreefallRadiusNumerical

  double precision function rootRadiusFreefall(radiusFreefall)
    !% Root function used in finding the radius corresponding to a given freefall time.
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
    !% Integrand for freefall time in the halo.
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
    !% Returns the rate of increase of the freefall radius in the dark matter density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
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
    !% GSL-callable function to evaluate the freefall radius of the dark matter profile.
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    genericFreefallRadiusEvaluate=solvers(solversCount)%self%freefallRadiusNumerical(solvers(solversCount)%node,exp(timeLogarithmic))
    return
  end function genericFreefallRadiusEvaluate

  double precision function genericRadiusEnclosingDensityNumerical(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , intent(in   )         :: density
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    call self%solverSet  (node)
    genericDensity =  density
    finder      =  rootFinder(                                                             &
         &                    rootFunction                 =rootDensity                  , &
         &                    toleranceAbsolute            =toleranceAbsolute            , &
         &                    toleranceRelative            =toleranceRelative            , &
         &                    rangeExpandDownward          =0.5d0                        , &
         &                    rangeExpandUpward            =2.0d0                        , &
         &                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                   )
    genericRadiusEnclosingDensityNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node))
    call self%solverUnset(   )
    return
  end function genericRadiusEnclosingDensityNumerical

  double precision function rootDensity(radius)
    !% Root function used in finding the radius enclosing a given mean density.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                           &
         &      *solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius)    &
         &      /4.0d0                                           &
         &      /Pi                                              &
         &      /                                     radius **3 &
         &      -genericDensity
    return
  end function rootDensity

  double precision function genericRadiusEnclosingMassNumerical(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
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
    genericRadiusEnclosingMassNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%virialRadius(node)])
    call self%solverUnset(   )
    return
  end function genericRadiusEnclosingMassNumerical

  double precision function rootMass(radius)
    !% Root function used in finding the radius enclosing a given mass.
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius) &
         &   -             genericMass
    return
  end function rootMass

  double precision function genericRadiusCircularVelocityMaximumNumerical(self,node)
    !% Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.
    use :: Numerical_Comparison, only : Values_Agree
    use :: Root_Finder         , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
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
    ! close to zero at the virial radius (which it will be for an isothermal profile), and return the circular velocity at that
    ! radius if so. Otherwise solve for the radius corresponding to the maximum circular velocity.
    if     (                                                                                                           &
         &  Values_Agree(                                                                                              &
         &                      +rootCircularVelocityMaximum   (     self%darkMatterHaloScale_%virialRadius(node))   , &
         &                      +0.0d0                                                                               , &
         &               absTol=+toleranceRelative                                                                     &
         &                      *self%circularVelocityNumerical(node,self%darkMatterHaloScale_%virialRadius(node))**2  &
         &                      /                                    self%darkMatterHaloScale_%virialRadius(node)      &
         &                                                                                                             &
         &               )                                                                                             &
         & ) then
       genericRadiusCircularVelocityMaximumNumerical=                      self%darkMatterHaloScale_%virialRadius(node)
    else
       genericRadiusCircularVelocityMaximumNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node))
    end if
    call self%solverUnset(   )
    return
  end function genericRadiusCircularVelocityMaximumNumerical

  double precision function genericCircularVelocityMaximumNumerical(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileGeneric), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node

    genericCircularVelocityMaximumNumerical=self%circularVelocityNumerical(node,self%radiusCircularVelocityMaximumNumerical(node))
    return
  end function genericCircularVelocityMaximumNumerical

  double precision function rootCircularVelocityMaximum(radius)
    !% Root function used in finding the radius at which the maximum circular velocity occurs. Since for a spherical profile $V_\mathrm{c}^2(r)=\mathrm{G}M(r)/r$, then
    !% \begin{equation}
    !% {\mathrm{d} V_\mathrm{c}^2 \over \mathrm{d} r} = - {\mathrm{G} M(r) \over r^2} + 4 \pi \mathrm{G} \rho(r) r.
    !% \end{equation}
    !% Therefore, the peak of the rotation curve satisfies $4 \pi r^3 \rho(r) - M(r)=0$.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootCircularVelocityMaximum=+4.0d0                                           &
         &                      *Pi                                              &
         &                      *                                     radius **3 &
         &                      *solvers(solversCount)%self%density     (solvers(solversCount)%node,radius)    &
         &                      -solvers(solversCount)%self%enclosedMass(solvers(solversCount)%node,radius)
    return
  end function rootCircularVelocityMaximum

  double precision function genericRadiusFromSpecificAngularMomentumNumerical(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
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
    genericRadiusFromSpecificAngularMomentumNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%virialRadius(node)])
    call self%solverUnset(   )
    return
  end function genericRadiusFromSpecificAngularMomentumNumerical

  double precision function rootSpecificAngularMomentum(radius)
    !% Root function used in finding the radius enclosing a given specific angular momentum.
    implicit none
    double precision, intent(in   ) :: radius

    rootSpecificAngularMomentum=+solvers(solversCount)%self%circularVelocityNumerical(solvers(solversCount)%node,radius) &
         &                      *                                                  radius &
         &                      -genericSpecificAngularMomentum
    return
  end function rootSpecificAngularMomentum

  double precision function genericDensityLogSlopeNumerical(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
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
    !% GSL-callable function to evaluate the density of the dark matter profile.
    implicit none
    double precision, intent(in   ), value :: radiusLogarithmic

    genericDensityEvaluate=solvers(solversCount)%self%density(solvers(solversCount)%node,exp(radiusLogarithmic))
    return
  end function genericDensityEvaluate

  subroutine genericSolverSet(self,node)
    !% Set a sub-module scope pointers on a stack to allow recursive calls to functions.
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
    !% Unset a sub-module scope pointers on the stack.
    implicit none

    solvers(solversCount)%self => null()
    solvers(solversCount)%node => null()
    solversCount=solversCount-1
    return
  end subroutine genericSolverUnset

end module Dark_Matter_Profiles_Generic
