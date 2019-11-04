!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     !@ <objectMethods>
     !@   <object>darkMatterProfileGeneric</object>
     !@   <objectMethod>
     !@     <method>enclosedMass</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\argin</arguments>
     !@     <description>Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>density</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\argin</arguments>
     !@     <description>Returns the density (in $M_\odot/$Mpc$^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enclosedMassNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\argin</arguments>
     !@     <description>Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enclosedMassDifferenceNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout, \textcolor{red}{\textless double precision\textgreater} radiusLower\argin, \textcolor{red}{\textless double precision\textgreater} radiusUpper\argin</arguments>
     !@     <description>Returns the enclosed mass difference (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>potentialNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\argin,\textcolor{red}{\textless integer\textgreater} status\argout</arguments>
     !@     <description>Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc) using a numerical calculation.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>potentialDifferenceNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout, \textcolor{red}{\textless double precision\textgreater} radiusLower\argin, \textcolor{red}{\textless double precision\textgreater} radiusUpper\argin</arguments>
     !@     <description>Returns the gravitational potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical calculation.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>circularVelocityNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\\argin</arguments>
     !@     <description>Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radialVelocityDispersionNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\\argin</arguments>
     !@     <description>Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radialMomentNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout, \doublezero\ moment\argin, \doublezero\ radiusMinimum\argin, \doublezero\ radiusMaximum\argin</arguments>
     !@     <description>Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rotationNormalizationNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout</arguments>
     !@     <description>Return the normalization of the rotation velocity vs. specific angular momentum relation.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>kSpaceNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout, \doublezero\ waveNumber\argin</arguments>
     !@     <description>Returns the Fourier transform of the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>energyNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout</arguments>
     !@     <description>Return the energy of the dark matter density profile.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>energyGrowthRateNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout</arguments>
     !@     <description>Return the rate of growth of the energy of the dark matter density profile.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallRadiusNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} time\argin</arguments>
     !@     <description>Returns the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallRadiusIncreaseRateNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} time\argin</arguments>
     !@     <description>Returns the rate of increase of the freefall radius in the dark matter density profile at the specified {\normalfont \ttfamily time} (given in Gyr).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusEnclosingDensityNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} density\argin</arguments>
     !@     <description>Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusEnclosingMassNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} mass\argin</arguments>
     !@     <description>Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given {\normalfont \ttfamily mass} (given in units of $M_\odot$).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>circularVelocityMaximumNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout</arguments>
     !@     <description>Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusFromSpecificAngularMomentumNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} specificAngularMomentum\argin</arguments>
     !@     <description>Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given in units of km s$^{-1}$ Mpc).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityLogSlopeNumerical</method>
     !@     <type>double precision</type>
     !@     <arguments>\textcolor{red}{\textless type((treeNode))\textgreater} node\arginout,\textcolor{red}{\textless double precision\textgreater} radius\argin</arguments>
     !@     <description>Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).</description>
     !@   </objectMethod>
     !@ </objectMethods>
   contains
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
     procedure                                         :: radiusFromSpecificAngularMomentumNumerical => genericRadiusFromSpecificAngularMomentumNumerical
     procedure                                         :: densityLogSlopeNumerical                   => genericDensityLogSlopeNumerical
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
  class           (darkMatterProfileGeneric      ), pointer :: genericSelf
  type            (treeNode                      ), pointer :: genericNode
  class           (nodeComponentBasic            ), pointer :: genericBasic
  class           (nodeComponentDarkMatterProfile), pointer :: genericDarkMatterProfile
  double precision                                          :: genericTime                   , genericRadiusFreefall , genericDensity        , genericMass , &
       &                                                       genericSpecificAngularMomentum, genericMassGrowthRate , genericScaleGrowthRate, genericScale, &
       &                                                       genericShape                  , genericShapeGrowthRate
  !$omp threadprivate(genericSelf,genericNode,genericBasic,genericTime,genericRadiusFreefall,genericDensity,genericMass,genericSpecificAngularMomentum,genericMassGrowthRate,genericDarkMatterProfile,genericScaleGrowthRate,genericScale,genericShape,genericShapeGrowthRate)

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
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radiusLower         , radiusUpper
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    genericEnclosedMassDifferenceNumerical=+Integrate(                                        &
         &                                                              radiusLower         , &
         &                                                              radiusUpper         , &
         &                                                              genericMassIntegrand, &
         &                                                              integrandFunction   , &
         &                                                              integrationWorkspace, &
         &                                            toleranceAbsolute=0.0d+0              , &
         &                                            toleranceRelative=1.0d-6                &
         &                                           )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return

  contains

    double precision function genericMassIntegrand(radius)
      !% Integrand for mass in generic dark matter profiles.
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         genericMassIntegrand=4.0d0*Pi*radius**2*self%density(node,radius)
      else
         genericMassIntegrand=0.0d0
      end if
      return
    end function genericMassIntegrand

  end function genericEnclosedMassDifferenceNumerical

  double precision function genericPotentialNumerical(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc) using a numerical calculation.
    use :: FGSL                      , only : fgsl_function            , fgsl_integration_workspace
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    use :: Numerical_Integration     , only : Integrate                , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout), target   :: self
    type            (treeNode                  ), intent(inout), target   :: node
    double precision                            , intent(in   )           :: radius
    integer                                     , intent(  out), optional :: status
    double precision                            , parameter               :: radiusMaximumFactor =1.0d2
    type            (fgsl_function             )                          :: integrandFunction
    type            (fgsl_integration_workspace)                          :: integrationWorkspace
    double precision                                                      :: radiusMaximum

    if (present(status)) status=structureErrorCodeSuccess
    genericSelf               =>  self
    genericNode               =>  node
    radiusMaximum             =  +radiusMaximumFactor                          &
         &                       *self%darkMatterHaloScale_%virialRadius(node)
    genericPotentialNumerical =   Integrate(                             &
         &                                  radius                     , &
         &                                  radiusMaximum              , &
         &                                  integrandPotential         , &
         &                                  integrandFunction          , &
         &                                  integrationWorkspace       , &
         &                                  toleranceAbsolute   =0.0d+0, &
         &                                  toleranceRelative   =1.0d-6  &
         &                                 )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function genericPotentialNumerical

  double precision function genericPotentialDifferenceNumerical(self,node,radiusLower,radiusUpper)
    !% Returns the potential difference (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} between the
    !% given {\normalfont \ttfamily radiusLower} and {\normalfont \ttfamily radiusUpper} (given in units of Mpc) using a numerical
    !% calculation.
    use :: FGSL                      , only : fgsl_function            , fgsl_integration_workspace
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    use :: Numerical_Integration     , only : Integrate                , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout), target   :: self
    type            (treeNode                  ), intent(inout), pointer  :: node
    double precision                            , intent(in   )           :: radiusLower         , radiusUpper
    type            (fgsl_function             )                          :: integrandFunction
    type            (fgsl_integration_workspace)                          :: integrationWorkspace

    genericSelf                         => self
    genericNode                         => node
    genericPotentialDifferenceNumerical =  Integrate(                             &
         &                                           radiusLower                , &
         &                                           radiusUpper                , &
         &                                           integrandPotential         , &
         &                                           integrandFunction          , &
         &                                           integrationWorkspace       , &
         &                                           toleranceAbsolute   =0.0d+0, &
         &                                           toleranceRelative   =1.0d-6  &
         &                                          )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function genericPotentialDifferenceNumerical

  double precision function integrandPotential(radius)
    !% Integrand for gravitational potential in a generic dark matter profile.
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius

    if (radius > 0.0d0) then
       integrandPotential=-gravitationalConstantGalacticus                 &
            &             *genericSelf%enclosedMass(genericNode,radius)    &
            &             /                                     radius **2
    else
       integrandPotential=0.0d0
    end if
    return
  end function integrandPotential

  double precision function genericCircularVelocityNumerical(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
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
    use :: FGSL                        , only : fgsl_function                  , fgsl_integration_workspace
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration       , only : Integrate                      , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout)            :: self
    type            (treeNode                  ), intent(inout)            :: node
    double precision                            , intent(in   )            :: radius
    double precision                                           , parameter :: radiusTinyFraction =1.0d-9, radiusLargeFactor=5.0d2
    double precision                                                       :: radiusMinimum             , radiusMaximum          , &
         &                                                                    radiusVirial
    type            (fgsl_function             )                           :: integrandFunction
    type            (fgsl_integration_workspace)                           :: integrationWorkspace

    radiusVirial =self%darkMatterHaloScale_%virialRadius(node)
    radiusMinimum=max(       radius,radiusTinyFraction*radiusVirial)
    radiusMaximum=max(10.0d0*radius,radiusLargeFactor *radiusVirial)
    genericRadialVelocityDispersionNumerical=sqrt(                                              &
         &                                        +Integrate(                                   &
         &                                                   radiusMinimum                    , &
         &                                                   radiusMaximum                    , &
         &                                                   genericJeansEquationIntegrand    , &
         &                                                   integrandFunction                , &
         &                                                   integrationWorkspace             , &
         &                                                   toleranceAbsolute    =0.0d+0     , &
         &                                                   toleranceRelative    =1.0d-6       &
         &                                                  )                                   &
         &                                         /self%density(node,radius)                   &
         &                                        )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return

  contains

    double precision function genericJeansEquationIntegrand(radius)
      !% Integrand for generic drak matter profile Jeans equation.
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         genericJeansEquationIntegrand=+gravitationalConstantGalacticus &
              &                        *self%enclosedMass(node,radius)  &
              &                        *self%density     (node,radius)  &
              &                        /radius**2
      else
         genericJeansEquationIntegrand=0.0d0
      end if
      return
    end function genericJeansEquationIntegrand

  end function genericRadialVelocityDispersionNumerical

  double precision function genericRadialMomentNumerical(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given
    !% {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout)           :: self
    type            (treeNode                  ), intent(inout)           :: node
    double precision                            , intent(in   )           :: moment
    double precision                            , intent(in   ), optional :: radiusMinimum       , radiusMaximum
    type            (fgsl_function             )                          :: integrandFunction
    type            (fgsl_integration_workspace)                          :: integrationWorkspace
    double precision                                                      :: radiusMinimumActual , radiusMaximumActual

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%virialRadius(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    genericRadialMomentNumerical=Integrate(                              &
         &                                 radiusMinimumActual         , &
         &                                 radiusMaximumActual         , &
         &                                 integrandRadialMoment       , &
         &                                 integrandFunction           , &
         &                                 integrationWorkspace        , &
         &                                 toleranceAbsolute    =0.0d+0, &
         &                                 toleranceRelative    =1.0d-3  &
         &                                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
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
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout)         :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                            , intent(in   )         :: waveNumber
    type            (fgsl_function             )                        :: integrandFunction
    type            (fgsl_integration_workspace)                        :: integrationWorkspace
    double precision                                                    :: radiusVirial

    radiusVirial          =+self%darkMatterHaloScale_%virialRadius(node             )
    genericKSpaceNumerical=+Integrate(                                                &
         &                            lowerLimit          =0.0d0                    , &
         &                            upperLimit          =radiusVirial             , &
         &                            integrand           =integrandFourierTransform, &
         &                            integrandFunction   =integrandFunction        , &
         &                            integrationWorkspace=integrationWorkspace     , &
         &                            toleranceAbsolute   =0.0d+0                   , &
         &                            toleranceRelative   =1.0d-3                     &
         &                           )                                                &
         &                  /self                    %enclosedMass(node,radiusVirial)
    call Integrate_Done(integrandFunction,integrationWorkspace)
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
    use :: FGSL                        , only : fgsl_function                  , fgsl_integration_workspace
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    use :: Numerical_Integration       , only : Integrate                      , Integrate_Done
    implicit none
    class           (darkMatterProfileGeneric  ), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , parameter     :: multiplierRadius    =100.0d0
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: radiusVirial                , radiusLarge, energyPotential, energyKinetic, pseudoPressure

    radiusVirial          =+self%darkMatterHaloScale_%virialRadius(node)
    radiusLarge           =+multiplierRadius                                         &
         &                 *radiusVirial
    energyPotential       =+Integrate(                                               &
         &                            lowerLimit          =0.0d0                   , &
         &                            upperLimit          =radiusVirial            , &
         &                            integrand           =integrandEnergyPotential, &
         &                            integrandFunction   =integrandFunction       , &
         &                            integrationWorkspace=integrationWorkspace    , &
         &                            toleranceAbsolute   =0.0d+0                  , &
         &                            toleranceRelative   =1.0d-3                    &
         &                           )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    energyKinetic         =+Integrate(                                               &
         &                            lowerLimit          =0.0d0                   , &
         &                            upperLimit          =radiusVirial            , &
         &                            integrand           =integrandEnergyKinetic  , &
         &                            integrandFunction   =integrandFunction       , &
         &                            integrationWorkspace=integrationWorkspace    , &
         &                            toleranceAbsolute   =0.0d+0                  , &
         &                            toleranceRelative   =1.0d-3                    &
         &                           )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    pseudoPressure        =+Integrate(                                               &
         &                            lowerLimit          =radiusVirial            , &
         &                            upperLimit          =radiusLarge             , &
         &                            integrand           =integrandPseudoPressure , &
         &                            integrandFunction   =integrandFunction       , &
         &                            integrationWorkspace=integrationWorkspace    , &
         &                            toleranceAbsolute   =0.0d+0                  , &
         &                            toleranceRelative   =1.0d-3                    &
         &                           )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    genericEnergyNumerical=-0.5d0                                                    &
         &                 *gravitationalConstantGalacticus                          &
         &                 *(                                                        &
         &                   +energyPotential                                        &
         &                   +self%enclosedMass(node,radiusVirial)**2                &
         &                   /                       radiusVirial                    &
         &                  )                                                        &
         &                 +2.0d0                                                    &
         &                 *Pi                                                       &
         &                 *gravitationalConstantGalacticus                          &
         &                 *(                                                        &
         &                   +radiusVirial**3                                        &
         &                   *pseudoPressure                                         &
         &                   +energyKinetic                                          &
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

    differentiator_                  =   differentiator                            (genericEnergyEvaluate                    )
    genericSelf                      =>  self
    genericNode                      =>  node
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
    call Galacticus_Calculations_Reset_(genericNode)
    genericEnergyEvaluate=genericSelf%energyNumerical(genericNode)
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

    genericSelf => self
    genericNode => node
    genericTime =  time
    call finder%rootFunction(rootRadiusFreefall                  )
    call finder%tolerance   (toleranceAbsolute ,toleranceRelative)
    call finder%rangeExpand (                                                             &
         &                   rangeExpandDownward          =0.5d0                        , &
         &                   rangeExpandUpward            =2.0d0                        , &
         &                   rangeExpandType              =rangeExpandMultiplicative    , &
         &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                  )
    genericFreefallRadiusNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node))
    return
  end function genericFreefallRadiusNumerical

  double precision function rootRadiusFreefall(radiusFreefall)
    !% Root function used in finding the radius corresponding to a given freefall time.
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    double precision                            , intent(in   ) :: radiusFreefall
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    genericRadiusFreefall=radiusFreefall
    rootRadiusFreefall=+Integrate(                                            &
         &                        lowerLimit          =0.0d0                , &
         &                        upperLimit          =radiusFreefall       , &
         &                        integrand           =integrandTimeFreefall, &
         &                        integrandFunction   =integrandFunction    , &
         &                        integrationWorkspace=integrationWorkspace , &
         &                        toleranceAbsolute   =0.0d+0               , &
         &                        toleranceRelative   =1.0d-3                 &
         &                       )                                            &
         &                      -genericTime
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function rootRadiusFreefall

  double precision function integrandTimeFreefall(radius)
    !% Integrand for freefall time in the halo.
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: potentialDifference

    potentialDifference=+genericSelf%potentialDifferenceNumerical(genericNode,radius,genericRadiusFreefall)
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

    genericSelf                                =>  self
    genericNode                                =>  node
    differentiator_                            =   differentiator            (genericFreefallRadiusEvaluate                    )
    genericFreefallRadiusIncreaseRateNumerical =  +differentiator_%derivative(log(time)                    ,timeLogarithmicStep) &
         &                                        /                               time
    return
  end function genericFreefallRadiusIncreaseRateNumerical

  double precision function genericFreefallRadiusEvaluate(timeLogarithmic)
    !% GSL-callable function to evaluate the freefall radius of the dark matter profile.
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    genericFreefallRadiusEvaluate=genericSelf%freefallRadiusNumerical(genericNode,exp(timeLogarithmic))
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

    genericSelf    => self
    genericNode    => node
    genericDensity =  density
    call finder%rootFunction(rootDensity                         )
    call finder%tolerance   (toleranceAbsolute ,toleranceRelative)
    call finder%rangeExpand (                                                             &
         &                   rangeExpandDownward          =0.5d0                        , &
         &                   rangeExpandUpward            =2.0d0                        , &
         &                   rangeExpandType              =rangeExpandMultiplicative    , &
         &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                  )
    genericRadiusEnclosingDensityNumerical=finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node))
    return
  end function genericRadiusEnclosingDensityNumerical

  double precision function rootDensity(radius)
    !% Root function used in finding the radius enclosing a given mean density.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                           &
         &      *genericSelf%enclosedMass(genericNode,radius)    &
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

    genericSelf => self
    genericNode => node
    genericMass =  mass
    call finder%rootFunction(rootMass                            )
    call finder%tolerance   (toleranceAbsolute ,toleranceRelative)
    call finder%rangeExpand (                                                             &
         &                   rangeExpandUpward            =2.0d0                        , &
         &                   rangeExpandType              =rangeExpandMultiplicative    , &
         &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                  )
    genericRadiusEnclosingMassNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%virialRadius(node)])
    return
  end function genericRadiusEnclosingMassNumerical

  double precision function rootMass(radius)
    !% Root function used in finding the radius enclosing a given mass.
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+genericSelf%enclosedMass(genericNode,radius) &
         &   -             genericMass
    return
  end function rootMass

  double precision function genericCircularVelocityMaximumNumerical(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use :: Numerical_Comparison, only : Values_Agree
    use :: Root_Finder         , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target :: self
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                        :: finder

    genericSelf => self
    genericNode => node
    call finder%rootFunction(rootCircularVelocityMaximum                  )
    call finder%tolerance   (toleranceAbsolute          ,toleranceRelative)
    call finder%rangeExpand (                                                             &
         &                   rangeExpandDownward          =0.5d0                        , &
         &                   rangeExpandUpward            =2.0d0                        , &
         &                   rangeExpandType              =rangeExpandMultiplicative    , &
         &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                  )
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
       genericCircularVelocityMaximumNumerical=self%circularVelocityNumerical(node,                      self%darkMatterHaloScale_%virialRadius(node) )
    else
       genericCircularVelocityMaximumNumerical=self%circularVelocityNumerical(node,finder%find(rootGuess=self%darkMatterHaloScale_%virialRadius(node)))
    end if
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
         &                      *genericSelf%density     (genericNode,radius)    &
         &                      -genericSelf%enclosedMass(genericNode,radius)
    return
  end function rootCircularVelocityMaximum

  double precision function genericRadiusFromSpecificAngularMomentumNumerical(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileGeneric), intent(inout), target  :: self
    type            (treeNode                ), intent(inout), pointer :: node
    double precision                          , intent(in   )          :: specificAngularMomentum
    double precision                          , parameter              :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder              )                         :: finder

    genericSelf                    => self
    genericNode                    => node
    genericSpecificAngularMomentum =  specificAngularMomentum
    call finder%rootFunction(rootSpecificAngularMomentum                  )
    call finder%tolerance   (toleranceAbsolute          ,toleranceRelative)
    call finder%rangeExpand (                                                             &
         &                   rangeExpandUpward            =2.0d0                        , &
         &                   rangeExpandType              =rangeExpandMultiplicative    , &
         &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                  )
    genericRadiusFromSpecificAngularMomentumNumerical=finder%find(rootRange=[0.0d0,self%darkMatterHaloScale_%virialRadius(node)])
    return
  end function genericRadiusFromSpecificAngularMomentumNumerical

  double precision function rootSpecificAngularMomentum(radius)
    !% Root function used in finding the radius enclosing a given specific angular momentum.
    implicit none
    double precision, intent(in   ) :: radius

    rootSpecificAngularMomentum=+genericSelf%circularVelocityNumerical(genericNode,radius) &
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

    genericSelf                    =>   self
    genericNode                    =>   node
    differentiator_                =    differentiator                   (genericDensityEvaluate                      )
    genericDensityLogSlopeNumerical=   +differentiator_       %derivative(log(radius)           ,radiusLogarithmicStep) &
         &                             /genericDensityEvaluate           (log(radius)                                 )
    return
  end function genericDensityLogSlopeNumerical

  double precision function genericDensityEvaluate(radiusLogarithmic)
    !% GSL-callable function to evaluate the density of the dark matter profile.
    implicit none
    double precision, intent(in   ), value :: radiusLogarithmic

    genericDensityEvaluate=genericSelf%density(genericNode,exp(radiusLogarithmic))
    return
  end function genericDensityEvaluate

end module Dark_Matter_Profiles_Generic

