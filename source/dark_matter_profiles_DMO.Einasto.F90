!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of ``Einasto'' dark matter halo profiles.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  use :: Root_Finder            , only : rootFinder

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOEinasto">
   <description>
    A dark matter profile DMO class which implements the Einasto density profile (e.g. \citealt{cardone_spherical_2005}):
    \begin{equation}
      \rho_\mathrm{dark matter}(r) = \rho_{-2} \exp \left( - {2 \over \alpha} \left[ \left( {r \over r_{-2}} \right)^\alpha - 1
      \right] \right),
    \end{equation}
    normalized such that the total mass of the \gls{node} is enclosed with the virial radius and with the characteristic length
    $r_{-2} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see \refPhysics{darkMatterProfileConcentration}). The
    shape parameter, $\alpha$, is set using the density profile shape method (see \refPhysics{darkMatterProfileShape}).
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOEinasto
     !!{
     A dark matter halo profile class implementing ``Einasto'' dark matter halos.
     !!}
     private
     ! Tables for specific angular momentum vs. radius table
     double precision                                              :: angularMomentumTableRadiusMinimum
     double precision                                              :: angularMomentumTableRadiusMaximum
     double precision                                              :: angularMomentumTableAlphaMinimum
     double precision                                              :: angularMomentumTableAlphaMaximum
     logical                                                       :: angularMomentumTableInitialized
     integer                                                       :: angularMomentumTableAlphaCount                , angularMomentumTableRadiusCount
     double precision              , allocatable, dimension(:  )   :: angularMomentumTableAlpha                     , angularMomentumTableRadius
     double precision              , allocatable, dimension(:,:)   :: angularMomentumTable
     type            (interpolator), allocatable                   :: angularMomentumTableAlphaInterpolator
     type            (interpolator), allocatable, dimension(:  )   :: angularMomentumTableRadiusInterpolator
     ! Tables for freefall time vs. radius table
     double precision                                              :: freefallRadiusTableRadiusMinimum
     double precision                                              :: freefallRadiusTableRadiusMaximum
     double precision                                              :: freefallRadiusTableAlphaMinimum
     double precision                                              :: freefallRadiusTableAlphaMaximum
     logical                                                       :: freefallRadiusTableInitialized
     integer                                                       :: freefallRadiusTableAlphaCount                 , freefallRadiusTableRadiusCount
     double precision              , allocatable, dimension(:  )   :: freefallRadiusTableAlpha                      , freefallRadiusTableRadius
     double precision              , allocatable, dimension(:,:)   :: freefallRadiusTable
     double precision                                              :: freefallTimeMaximum                           , freefallTimeMinimum
     type            (interpolator), allocatable                   :: freefallRadiusTableAlphaInterpolator
     type            (interpolator), allocatable, dimension(:  )   :: freefallRadiusTableRadiusInterpolator
     ! Tables for radial velocity dispersion vs. radius table
     double precision                                              :: radialVelocityDispersionRadiusMinimum
     double precision                                              :: radialVelocityDispersionRadiusMaximum
     double precision                                              :: radialVelocityDispersionAlphaMinimum
     double precision                                              :: radialVelocityDispersionAlphaMaximum
     logical                                                       :: radialVelocityDispersionTableInitialized
     integer                                                       :: radialVelocityDispersionTableAlphaCount       , radialVelocityDispersionTableRadiusCount
     double precision              , allocatable, dimension(:  )   :: radialVelocityDispersionTableAlpha            , radialVelocityDispersionTableRadius
     double precision              , allocatable, dimension(:,:)   :: radialVelocityDispersionTable
     type            (interpolator), allocatable                   :: radialVelocityDispersionTableAlphaInterpolator, radialVelocityDispersionTableRadiusInterpolator
     ! Tables for energy as a function of concentration and alpha.
     double precision                                              :: energyTableConcentrationMinimum
     double precision                                              :: energyTableConcentrationMaximum
     double precision                                              :: energyTableAlphaMinimum
     double precision                                              :: energyTableAlphaMaximum
     logical                                                       :: energyTableInitialized
     integer                                                       :: energyTableAlphaCount                         , energyTableConcentrationCount
     double precision              , allocatable, dimension(:  )   :: energyTableAlpha                              , energyTableConcentration
     double precision              , allocatable, dimension(:,:)   :: energyTable
     type            (interpolator), allocatable                   :: energyTableAlphaInterpolator                  , energyTableConcentrationInterpolator
     ! Tables for specific Fourier transform of density profile as a function of alpha and radius.
     double precision                                              :: fourierProfileTableConcentrationMinimum

     double precision                                              :: fourierProfileTableConcentrationMaximum
     double precision                                              :: fourierProfileTableWavenumberMinimum
     double precision                                              :: fourierProfileTableWavenumberMaximum
     double precision                                              :: fourierProfileTableAlphaMinimum
     double precision                                              :: fourierProfileTableAlphaMaximum
     logical                                                       :: fourierProfileTableInitialized
     integer                                                       :: fourierProfileTableAlphaCount                 , fourierProfileTableConcentrationCount          , &
          &                                                           fourierProfileTableWavenumberCount
     double precision              , allocatable, dimension(:    ) :: fourierProfileTableAlpha                      , fourierProfileTableConcentration               , &
       &                                                              fourierProfileTableWavenumber
     double precision              , allocatable, dimension(:,:,:) :: fourierProfileTable
     type            (interpolator), allocatable                   :: fourierProfileTableAlphaInterpolator          , fourierProfileTableConcentrationInterpolator   , &
          &                                                           fourierProfileTableWavenumberInterpolator
     ! Root finders.
     type            (rootFinder  )                                :: finderEnclosedDensity                         , finderVelocityPeak
   contains
     !![
     <methods>
       <method description="Returns the density (in units such that the virial mass and scale length are unity) in an Einasto dark matter profile with given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="densityScaleFree" />
       <method description="Returns the enclosed mass (in units of the virial mass) in an Einasto dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="enclosedMassScaleFree" />
       <method description="Returns the scale-free radial velocity dispersion in an Einasto dark matter profile." method="radialVelocityDispersionScaleFree" />
       <method description="Tabulates the radial velocity dispersion vs. radius for Einasto halos." method="radialVelocityDispersionTabulate" />
       <method description="Create a tabulation of the energy of Einasto profiles as a function of their concentration of $\alpha$ parameter." method="energyTableMake" />
       <method description="Create a tabulation of the Fourier transform of Einasto profiles as a function of their $\alpha$ parameter and dimensionless wavenumber." method="fourierProfileTableMake" />
       <method description="Tabulates the freefall time vs. freefall radius for Einasto halos." method="freefallTabulate" />
       <method description="Compute the freefall time in a scale-free Einasto halo." method="freefallTimeScaleFree" />
       <method description="Returns the gravitational potential (in units where the virial mass and scale radius are unity) in an Einasto dark matter profile with given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="potentialScaleFree" />
       <method description=" Comptue the radius at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentumScaleFree} in a scale free Einasto profile." method="radiusFromSpecificAngularMomentumScaleFree" />
       <method description="Create a tabulation of the relation between specific angular momentum and radius in an Einasto profile." method="radiusFromSpecificAngularMomentumTableMake" />
     </methods>
     !!]
     final     ::                                               einastoDestructor
     procedure :: density                                    => einastoDensity
     procedure :: densityLogSlope                            => einastoDensityLogSlope
     procedure :: radialMoment                               => einastoRadialMoment
     procedure :: enclosedMass                               => einastoEnclosedMass
     procedure :: radiusEnclosingDensity                     => einastoRadiusEnclosingDensity
     procedure :: potential                                  => einastoPotential
     procedure :: circularVelocity                           => einastoCircularVelocity
     procedure :: circularVelocityMaximum                    => einastoCircularVelocityMaximum
     procedure :: radialVelocityDispersion                   => einastoRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum          => einastoRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization                      => einastoRotationNormalization
     procedure :: energy                                     => einastoEnergy
     procedure :: kSpace                                     => einastoKSpace
     procedure :: freefallRadius                             => einastoFreefallRadius
     procedure :: freefallRadiusIncreaseRate                 => einastoFreefallRadiusIncreaseRate
     procedure :: densityScaleFree                           => einastoDensityScaleFree
     procedure :: enclosedMassScaleFree                      => einastoEnclosedMassScaleFree
     procedure :: radialVelocityDispersionScaleFree          => einastoRadialVelocityDispersionScaleFree
     procedure :: radialVelocityDispersionTabulate           => einastoRadialVelocityDispersionTabulate
     procedure :: freefallTimeScaleFree                      => einastoFreefallTimeScaleFree
     procedure :: potentialScaleFree                         => einastoPotentialScaleFree
     procedure :: radiusFromSpecificAngularMomentumScaleFree => einastoRadiusFromSpecificAngularMomentumScaleFree
     procedure :: radiusFromSpecificAngularMomentumTableMake => einastoRadiusFromSpecificAngularMomentumTableMake
     procedure :: freefallTabulate                           => einastoFreefallTabulate
     procedure :: energyTableMake                            => einastoEnergyTableMake
     procedure :: fourierProfileTableMake                    => einastoFourierProfileTableMake
  end type darkMatterProfileDMOEinasto

  interface darkMatterProfileDMOEinasto
     !!{
     Constructors for the {\normalfont \ttfamily einasto} dark matter halo profile class.
     !!}
     module procedure einastoConstructorParameters
     module procedure einastoConstructorInternal
  end interface darkMatterProfileDMOEinasto

  ! Granularity parameters for tabulations.
  integer, parameter :: einastoAngularMomentumTableRadiusPointsPerDecade         =100
  integer, parameter :: einastoAngularMomentumTableAlphaPointsPerUnit            =100
  integer, parameter :: einastoFreefallRadiusTableRadiusPointsPerDecade          = 30
  integer, parameter :: einastoFreefallRadiusTableAlphaPointsPerUnit             = 30
  integer, parameter :: einastoEnergyTableConcentrationPointsPerDecade           =100
  integer, parameter :: einastoEnergyTableAlphaPointsPerUnit                     =100
  integer, parameter :: einastoFourierProfileTableConcentrationPointsPerDecade   =100
  integer, parameter :: einastoFourierProfileTableWavenumberPointsPerDecade      =100
  integer, parameter :: einastoFourierProfileTableAlphaPointsPerUnit             =100
  integer, parameter :: einastoRadialVelocityDispersionTableRadiusPointsPerDecade=100
  integer, parameter :: einastoRadialVelocityDispersionTableAlphaPointsPerUnit   =100

  ! Module-scope variables used in root finding.
  class           (darkMatterProfileDMOEinasto), pointer :: einastoSelf
  type            (treeNode                   ), pointer :: einastoNode
  double precision                                       :: einastoDensityEnclosed, einastoAlpha
  !$omp threadprivate(einastoSelf,einastoNode,einastoDensityEnclosed,einastoAlpha)

contains

  function einastoConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily einasto} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOEinasto )                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOEinasto(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function einastoConstructorParameters

  function einastoConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily einasto} dark matter halo profile class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Galacticus_Error, only : Galacticus_Component_List        , Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Root_Finder     , only : rangeExpandMultiplicative        , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
   implicit none
    type            (darkMatterProfileDMOEinasto)                        :: self
    class           (darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    double precision                             , parameter             :: toleranceAbsolute   =0.0d0, toleranceRelative=1.0d-3
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    ! Initialize table states.
    self%angularMomentumTableRadiusMinimum       = 1.0d-3
    self%angularMomentumTableRadiusMaximum       =20.0d+0
    self%angularMomentumTableAlphaMinimum        = 0.1d+0
    self%angularMomentumTableAlphaMaximum        = 0.3d+0
    self%angularMomentumTableInitialized         =.false.
    self%freefallRadiusTableRadiusMinimum        = 1.0d-3
    self%freefallRadiusTableRadiusMaximum        =20.0d+0
    self%freefallRadiusTableAlphaMinimum         = 0.1d+0
    self%freefallRadiusTableAlphaMaximum         = 0.3d+0
    self%freefallRadiusTableInitialized          =.false.
    self%energyTableConcentrationMinimum         = 2.0d0
    self%energyTableConcentrationMaximum         =20.0d0
    self%energyTableAlphaMinimum                 = 0.1d0
    self%energyTableAlphaMaximum                 = 0.3d0
    self%energyTableInitialized                  =.false.
    self%fourierProfileTableConcentrationMinimum = 2.0d0
    self%fourierProfileTableConcentrationMaximum =20.0d0
    self%fourierProfileTableWavenumberMinimum    = 1.0d-3
    self%fourierProfileTableWavenumberMaximum    = 1.0d+3
    self%fourierProfileTableAlphaMinimum         = 0.1d+0
    self%fourierProfileTableAlphaMaximum         = 0.3d+0
    self%fourierProfileTableInitialized          =.false.
    self%radialVelocityDispersionRadiusMinimum   = 1.0d-3
    self%radialVelocityDispersionRadiusMaximum   =20.0d+0
    self%radialVelocityDispersionAlphaMinimum    = 0.1d+0
    self%radialVelocityDispersionAlphaMaximum    = 0.3d+0
    self%radialVelocityDispersionTableInitialized=.false.
    ! Initialize root finders.
    self%finderEnclosedDensity=rootFinder(                                                                 &
         &                                rootFunction                 =einastoRadiusEnclosingDensityRoot, &
         &                                toleranceAbsolute            =toleranceAbsolute                , &
         &                                toleranceRelative            =toleranceRelative                , &
         &                                rangeExpandUpward            =2.0d0                            , &
         &                                rangeExpandDownward          =0.5d0                            , &
         &                                rangeExpandType              =rangeExpandMultiplicative        , &
         &                                rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative    , &
         &                                rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive      &
         &                               )
    self%finderVelocityPeak   =rootFinder(&
         &                                rootFunction                 =einastoCircularVelocityPeakRadius, &
         &                                toleranceRelative            =toleranceRelative                , &
         &                                rangeExpandUpward            =2.0d0                            , &
         &                                rangeExpandDownward          =0.5d0                            , &
         &                                rangeExpandType              =rangeExpandMultiplicative          &
         &                               )
    ! Ensure that the dark matter profile component supports both "scale" and "shape" properties. Since we've been called with
    ! a treeNode to process, it should have been initialized by now.
    if     (                                                                                                                 &
         &  .not.(                                                                                                           &
         &         defaultDarkMatterProfileComponent%scaleIsGettable()                                                       &
         &        .and.                                                                                                      &
         &         defaultDarkMatterProfileComponent%shapeIsGettable()                                                       &
         &       )                                                                                                           &
         & ) then
       call Galacticus_Error_Report                                                                                               &
            &        (                                                                                                            &
            &         'Einasto dark matter profile requires a dark matter profile component that supports gettable '          //  &
            &         '"scale" and "shape" properties.'                                                                       //  &
            &         Galacticus_Component_List(                                                                                  &
            &                                   'darkMatterProfile'                                                             , &
            &                                    defaultDarkMatterProfileComponent%shapeAttributeMatch(requireGettable=.true.)    &
            &                                   .intersection.                                                                    &
            &                                    defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)    &
            &                                  )                                                                              //  &
            &         {introspection:location}                                                                                    &
            &        )
    end if
    ! Initialize the tabulations.
    call self%radialVelocityDispersionTabulate()
    return
  end function einastoConstructorInternal

  subroutine einastoDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily einasto} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOEinasto), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine einastoDestructor

  double precision function einastoDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha            , radiusOverScaleRadius      , &
         &                                                             scaleRadius      , virialRadiusOverScaleRadius

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    radiusOverScaleRadius      =radius                                      /scaleRadius
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    einastoDensity             =self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha) &
         &                      *basic%mass()/scaleRadius**3
    return
  end function einastoDensity

  double precision function einastoDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha            , radiusOverScaleRadius      , &
         &                                                             scaleRadius
    !$GLC attributes unused :: self

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =  darkMatterProfile%scale()
    alpha                          =  darkMatterProfile%shape()
    radiusOverScaleRadius          =  radius/scaleRadius
    einastoDensityLogSlope         = -2.0d0                        &
         &                           *radiusOverScaleRadius**alpha
    return
  end function einastoDensityLogSlope

  double precision function einastoRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    in units of Mpc).
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile         , treeNode
    use :: Gamma_Functions         , only : Gamma_Function    , Gamma_Function_Incomplete_Complementary
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout)           :: node
    double precision                                , intent(in   )           :: moment
    double precision                                , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    class           (nodeComponentBasic            )               , pointer  :: basic
    class           (nodeComponentDarkMatterProfile)               , pointer  :: darkMatterProfile
    double precision                                                          :: scaleRadius        , virialRadiusOverScaleRadius, &
         &                                                                       radiusMinimumActual, radiusMaximumActual        , &
         &                                                                       alpha              , densityNormalization

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%radiusVirial(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    densityNormalization= (alpha/4.0d0/Pi)                                                                                    &
         &               *   ((2.0d0/alpha)                    **(3.0d0/alpha)                                              ) &
         &               *exp(-2.0d0/alpha                                                                                  ) &
         &               /Gamma_Function                         (3.0d0/alpha                                               ) &
         &               /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha)
    einastoRadialMoment=+densityNormalization                                            &
         &              *basic%mass()                                                    &
         &              *scaleRadius**(moment-2.0d0)                                     &
         &              *(                                                               &
         &                +einastoRadialMomentScaleFree(radiusMaximumActual/scaleRadius) &
         &                -einastoRadialMomentScaleFree(radiusMinimumActual/scaleRadius) &
         &               )
    return

  contains

    double precision function einastoRadialMomentScaleFree(radius)
      !!{
      Provides the scale-free part of the radial moment of the Einasto density profile.
      !!}
      use :: Gamma_Functions, only : Gamma_Function_Incomplete
      implicit none
      double precision, intent(in   ) :: radius

      einastoRadialMomentScaleFree=-exp(2.0d0/alpha)                                                          &
           &                       *0.5d0   **((1.0d0+moment)/alpha)                                          &
           &                       *alpha   **((1.0d0+moment)/alpha)                                          &
           &                       *Gamma_Function           ((1.0d0+moment)/alpha                          ) &
           &                       *Gamma_Function_Incomplete((1.0d0+moment)/alpha,2.0d0*radius**alpha/alpha) &
           &                       /alpha
      return
    end function einastoRadialMomentScaleFree

  end function einastoRadialMoment

  double precision function einastoEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha            , radiusOverScaleRadius      , &
         &                                                             scaleRadius      , virialRadiusOverScaleRadius

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    radiusOverScaleRadius      =radius                                      /scaleRadius
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    einastoEnclosedMass        =self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha) &
         &                      *basic%mass()
    return
  end function einastoEnclosedMass

  double precision function einastoCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius

    if (radius > 0.0d0) then
       einastoCircularVelocity=sqrt(gravitationalConstantGalacticus*self%enclosedMass(node,radius)/radius)
    else
       einastoCircularVelocity=0.0d0
    end if
    return
  end function einastoCircularVelocity

  double precision function einastoCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha            , radiusScale, &
         &                                                             radiusPeak

    ! Get the shape parameter for this halo.
    darkMatterProfile => node                 %darkMatterProfile(autoCreate=.true.)
    alpha             =  darkMatterProfile%shape            (                 )
    radiusScale       =  darkMatterProfile%scale            (                 )
    ! Solve for the radius (in units of the scale radius) at which the rotation curve peaks.
    einastoAlpha=alpha
    radiusPeak  =self%finderVelocityPeak%find(rootGuess=radiusScale)
    ! Find the peak velocity.
    einastoCircularVelocityMaximum=self%circularVelocity(node,radiusPeak*radiusScale)
    return
  end function einastoCircularVelocityMaximum
  
  double precision function einastoCircularVelocityPeakRadius(radius)
    !!{
    Computes the derivative of the square of circular velocity for an Einasto density profile.
    !!}
    use :: Gamma_Functions, only : Gamma_Function, Gamma_Function_Incomplete_Complementary
    implicit none
    double precision, intent(in   ) :: radius
    
    einastoCircularVelocityPeakRadius=+        2.0d0                             **(      +3.0d0/einastoAlpha)          &
         &                            *        radius                            **(-2.0d0+      einastoAlpha)          &
         &                            *(       radius**einastoAlpha/einastoAlpha)**(-1.0d0+3.0d0/einastoAlpha)          &
         &                            * exp(                                                                            &
         &                                  -(                                                                          &
         &                                    +2.0d0                                                                    &
         &                                    *radius**einastoAlpha                                                     &
         &                                   )                                                                          &
         &                                  /einastoAlpha                                                               &
         &                                 )                                                                            &
         &                            /Gamma_Function                         (                                         &
         &                                                                     3.0d0                     /einastoAlpha  &
         &                                                                    )                                         &
         &                            -Gamma_Function_Incomplete_Complementary(                                         &
         &                                                                     3.0d0                     /einastoAlpha, &
         &                                                                     2.0d0*radius**einastoAlpha/einastoAlpha  &
         &                                                                    )                                         &
         &                            /       radius              **  2
    return
  end function einastoCircularVelocityPeakRadius
  
  double precision function einastoRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use            :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile         , treeNode
    use            :: Gamma_Functions , only : Gamma_Function_Incomplete_Complementary
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout)  :: self
    type            (treeNode                      ), intent(inout)  :: node
    double precision                                , intent(in   )  :: radius
    class           (nodeComponentDarkMatterProfile), pointer        :: darkMatterProfile
    integer         (c_size_t                      ), dimension(0:1) :: jAlpha
    double precision                                , dimension(0:1) :: hAlpha
    integer                                                          :: iAlpha
    double precision                                                 :: alpha                , scaleRadius                , &
         &                                                              radiusOverScaleRadius, virialRadiusOverScaleRadius

    darkMatterProfile           => node%darkMatterProfile(autoCreate=.true.)
    ! Get the scale radius.
    scaleRadius                 =  darkMatterProfile%scale()
    radiusOverScaleRadius       =  radius                                      /scaleRadius
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    ! Get the shape parameter.
    alpha                       =  darkMatterProfile%shape()
    if (radius > 0.0d0) then
       ! Ensure table is sufficiently extensive.
       call self%radialVelocityDispersionTabulate(radiusOverScaleRadius,alpha)
       ! Interpolate to get the radial velocity dispersion.
       ! Get interpolating factors in alpha.
       call self%radialVelocityDispersionTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
       jAlpha(1)=jAlpha(0)+1
       einastoRadialVelocityDispersion=0.0d0
       do iAlpha=0,1
          einastoRadialVelocityDispersion=+einastoRadialVelocityDispersion                                                                                                              &
               &                          +self%radialVelocityDispersionTableRadiusInterpolator%interpolate(radiusOverScaleRadius,self%radialVelocityDispersionTable(:,jAlpha(iAlpha))) &
               &                          *                                                                                                                            hAlpha(iAlpha)
       end do
       einastoRadialVelocityDispersion=+einastoRadialVelocityDispersion                                                                           &
            &                          *self%darkMatterHaloScale_%velocityVirial(node)                                                            &
            &                          *sqrt(                                                                                                     &
            &                                +virialRadiusOverScaleRadius                                                                         &
            &                                /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha) &
            &                               )
    else
       einastoRadialVelocityDispersion=0.0d0
    end if
    return
  end function einastoRadialVelocityDispersion

  double precision function einastoPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentBasic            )               , pointer  :: basic
    class           (nodeComponentDarkMatterProfile)               , pointer  :: darkMatterProfile
    double precision                                                          :: alpha            , radiusOverScaleRadius      , &
         &                                                                       scaleRadius      , virialRadiusOverScaleRadius

    ! Assume success.
    if (present(status)) status=structureErrorCodeSuccess
    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    radiusOverScaleRadius      =radius                       /scaleRadius
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    einastoPotential=+self%potentialScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha) &
         &           *gravitationalConstantGalacticus                                                  &
         &           *basic%mass()                                                                     &
         &           /scaleRadius
    return
  end function einastoPotential

  double precision function einastoRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha                           , scaleRadius, &
         &                                                             specificAngularMomentumScaleFree

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! Get the scale radius of the halo.
    scaleRadius=darkMatterProfile%scale()
    ! Get the shape parameter of the halo.
    alpha      =darkMatterProfile%shape()
    ! Compute the specific angular momentum in scale free units.
    specificAngularMomentumScaleFree=specificAngularMomentum                    &
         &                           /sqrt(                                     &
         &                                 +gravitationalConstantGalacticus     &
         &                                 *scaleRadius                         &
         &                                 *self%enclosedMass(node,scaleRadius) &
         &                                )
    ! Compute the corresponding radius.
    einastoRadiusFromSpecificAngularMomentum=scaleRadius                                                                              &
         &                                   *self%radiusFromSpecificAngularMomentumScaleFree(alpha,specificAngularMomentumScaleFree)
    return
  end function einastoRadiusFromSpecificAngularMomentum

  double precision function einastoRadiusFromSpecificAngularMomentumScaleFree(self,alpha,specificAngularMomentumScaleFree)
    !!{
    Comptue the radius at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentumScaleFree} in a scale free Einasto
    profile.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout)  :: self
    double precision                             , intent(in   )  :: alpha , specificAngularMomentumScaleFree
    integer         (c_size_t                   ), dimension(0:1) :: jAlpha
    double precision                             , dimension(0:1) :: hAlpha
    integer                                                       :: iAlpha

    ! Return immediately for zero angular momentum.
    if (specificAngularMomentumScaleFree <= 0.0d0) then
       einastoRadiusFromSpecificAngularMomentumScaleFree=0.0d0
       return
    end if

    ! Ensure the table exists and is sufficiently tabulated.
    call self%radiusFromSpecificAngularMomentumTableMake(alpha,specificAngularMomentumScaleFree)

    ! Get interpolating factors in alpha.
    call self%angularMomentumTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
    jAlpha(1)=jAlpha(0)+1

    ! Interpolate in specific angular momentum to get radius.
    einastoRadiusFromSpecificAngularMomentumScaleFree=0.0d0
    do iAlpha=0,1
       einastoRadiusFromSpecificAngularMomentumScaleFree=                                                                &
            &  einastoRadiusFromSpecificAngularMomentumScaleFree                                                         &
            & +self%angularMomentumTableRadiusInterpolator(jAlpha(iAlpha))%interpolate(specificAngularMomentumScaleFree) &
            & *                                            hAlpha(iAlpha)
    end do
    return
  end function einastoRadiusFromSpecificAngularMomentumScaleFree

  subroutine einastoRadiusFromSpecificAngularMomentumTableMake(self,alphaRequired,specificAngularMomentumRequired)
    !!{
    Create a tabulation of the relation between specific angular momentum and radius in an Einasto profile.
    !!}
    use :: Gamma_Functions  , only : Gamma_Function_Incomplete_Complementary
    use :: Memory_Management, only : allocateArray                          , deallocateArray
    use :: Numerical_Ranges , only : Make_Range                             , rangeTypeLinear, rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alphaRequired, specificAngularMomentumRequired
    integer                                                      :: iAlpha       , iRadius
    logical                                                      :: makeTable
    double precision                                             :: alpha        , enclosedMass                   , &
         &                                                          radius

    ! Always check if we need to make the table.
    makeTable=.true.
    do while (makeTable)
       ! Assume table does not need remaking.
       makeTable=.false.
       ! Check for uninitialized table.
       if (.not.self%angularMomentumTableInitialized) then
          makeTable=.true.
          ! Check for alpha out of range.
       else if (alphaRequired < self%angularMomentumTableAlpha(1) .or. alphaRequired >&
            & self%angularMomentumTableAlpha(self%angularMomentumTableAlphaCount)) then
          makeTable=.true.
          ! Compute the range of tabulation and number of points to use.
          self%angularMomentumTableAlphaMinimum=min(self%angularMomentumTableAlphaMinimum,0.9d0*alphaRequired)
          self%angularMomentumTableAlphaMaximum=max(self%angularMomentumTableAlphaMaximum,1.1d0*alphaRequired)
          ! Check for angular momentum below minimum tabulated value.
       else if (any(specificAngularMomentumRequired < self%angularMomentumTable(1,:))) then
          makeTable=.true.
          self%angularMomentumTableRadiusMinimum=0.5d0*self%angularMomentumTableRadiusMinimum
          ! Check for angular momentum above maximum tabulated value.
       else if (any(specificAngularMomentumRequired > self%angularMomentumTable(self%angularMomentumTableRadiusCount,:))) then
          makeTable=.true.
          self%angularMomentumTableRadiusMaximum=2.0d0*self%angularMomentumTableRadiusMaximum
       end if
       ! Remake the table if necessary.
       if (makeTable) then
          ! Allocate arrays to the appropriate sizes.
          self%angularMomentumTableAlphaCount =int(     (self%angularMomentumTableAlphaMaximum -self%angularMomentumTableAlphaMinimum ) &
               &                                   *dble(einastoAngularMomentumTableAlphaPointsPerUnit   )                              &
               &                                  )+1
          self%angularMomentumTableRadiusCount=int(log10(self%angularMomentumTableRadiusMaximum/self%angularMomentumTableRadiusMinimum) &
               &                                   *dble(einastoAngularMomentumTableRadiusPointsPerDecade)                              &
               &                                  )+1
          if (allocated(self%angularMomentumTableAlpha )) call deallocateArray(self%angularMomentumTableAlpha )
          if (allocated(self%angularMomentumTableRadius)) call deallocateArray(self%angularMomentumTableRadius)
          if (allocated(self%angularMomentumTable      )) call deallocateArray(self%angularMomentumTable      )
          call allocateArray(self%angularMomentumTableAlpha ,[                                     self%angularMomentumTableAlphaCount])
          call allocateArray(self%angularMomentumTableRadius,[self%angularMomentumTableRadiusCount                                    ])
          call allocateArray(self%angularMomentumTable      ,[self%angularMomentumTableRadiusCount,self%angularMomentumTableAlphaCount])
          ! Create ranges of alpha and radius.
          self%angularMomentumTableAlpha =Make_Range(self%angularMomentumTableAlphaMinimum ,self%angularMomentumTableAlphaMaximum,  &
               &                                     self%angularMomentumTableAlphaCount   ,rangeType=rangeTypeLinear     )
          self%angularMomentumTableRadius=Make_Range(self%angularMomentumTableRadiusMinimum,self%angularMomentumTableRadiusMaximum, &
               &                                     self%angularMomentumTableRadiusCount  ,rangeType=rangeTypeLogarithmic)
          ! Tabulate the radius vs. specific angular momentum relation.
          do iAlpha=1,self%angularMomentumTableAlphaCount
             alpha=self%angularMomentumTableAlpha(iAlpha)
             do iRadius=1,self%angularMomentumTableRadiusCount
                radius=self%angularMomentumTableRadius(iRadius)
                enclosedMass= Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius**alpha/alpha) &
                     &       /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0              /alpha)
                self%angularMomentumTable(iRadius,iAlpha)=sqrt(enclosedMass*radius)
             end do
          end do
          ! Build interpolators.
          if (allocated(self%angularMomentumTableRadiusInterpolator)) deallocate(self%angularMomentumTableRadiusInterpolator)
          if (allocated(self%angularMomentumTableAlphaInterpolator )) deallocate(self%angularMomentumTableAlphaInterpolator )
          allocate(self%angularMomentumTableRadiusInterpolator(self%angularMomentumTableAlphaCount))
          allocate(self%angularMomentumTableAlphaInterpolator                                      )
          do iAlpha=1,self%angularMomentumTableAlphaCount
             self%angularMomentumTableRadiusInterpolator(iAlpha)=interpolator(self%angularMomentumTable     (:,iAlpha),self%angularMomentumTableRadius)
          end do
          self%angularMomentumTableAlphaInterpolator            =interpolator(self%angularMomentumTableAlpha                                          )
          ! Flag that the table is now initialized.
          self%angularMomentumTableInitialized=.true.
       end if
    end do
    return
  end subroutine einastoRadiusFromSpecificAngularMomentumTableMake

  double precision function einastoRotationNormalization(self,node)
    !!{
    Return the rotation normalization of an Einasto halo density profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    use :: Gamma_Functions , only : Gamma_Function                , Gamma_Function_Incomplete_Complementary
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: alpha                      , scaleRadius, &
         &                                                             virialRadiusOverScaleRadius

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! Get scale radius, shape and concentration.
    scaleRadius                 =+darkMatterProfile%scale()
    alpha                       =+darkMatterProfile%shape()
    virialRadiusOverScaleRadius =+self             %darkMatterHaloScale_%radiusVirial(node) &
         &                       /                                       scaleRadius
    einastoRotationNormalization=+(2.0d0/alpha)**(1.0d0/alpha)                                                                        &
         &                       *Gamma_Function                         (3.0d0/alpha                                               ) &
         &                       /Gamma_Function                         (4.0d0/alpha                                               ) &
         &                       *Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha) &
         &                       /Gamma_Function_Incomplete_Complementary(4.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha) &
         &                       /scaleRadius
    return
  end function einastoRotationNormalization

  double precision function einastoEnergy(self,node)
    !!{
    Return the energy of an Einasto halo density profile.
    !!}
    use            :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout)  :: self
    type            (treeNode                      ), intent(inout)  :: node
    class           (nodeComponentBasic            ), pointer        :: basic
    class           (nodeComponentDarkMatterProfile), pointer        :: darkMatterProfile
    integer         (c_size_t                      ), dimension(0:1) :: jAlpha
    double precision                                , dimension(0:1) :: hAlpha
    integer                                                          :: iAlpha
    double precision                                                 :: alpha                      , scaleRadius, &
         &                                                              virialRadiusOverScaleRadius

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call self%energyTableMake(virialRadiusOverScaleRadius,alpha)

    ! Get interpolating factors in alpha.
    call self%energyTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
    jAlpha(1)=jAlpha(0)+1

    ! Find the energy by interpolation.
    einastoEnergy=0.0d0
    do iAlpha=0,1
       einastoEnergy=+einastoEnergy                                                                                                         &
            &        +self%energyTableConcentrationInterpolator%interpolate(virialRadiusOverScaleRadius,self%energyTable(:,jAlpha(iAlpha))) &
            &        *                                                                                                     hAlpha(iAlpha)
    end do

    ! Scale to dimensionful units.
    einastoEnergy=einastoEnergy*basic%mass()                         &
         &        *self%darkMatterHaloScale_%velocityVirial(node)**2
    return
  end function einastoEnergy

  subroutine einastoEnergyTableMake(self,concentrationRequired,alphaRequired)
    !!{
    Create a tabulation of the energy of Einasto profiles as a function of their concentration of $\alpha$ parameter.
    !!}
    use :: Memory_Management       , only : allocateArray, deallocateArray
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range   , rangeTypeLinear, rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alphaRequired          , concentrationRequired
    integer                                                      :: iAlpha                 , iConcentration
    logical                                                      :: makeTable
    double precision                                             :: alpha                  , concentration         , &
         &                                                          jeansEquationIntegral  , kineticEnergy         , &
         &                                                          kineticEnergyIntegral  , potentialEnergy       , &
         &                                                          potentialEnergyIntegral, radiusMaximum         , &
         &                                                          radiusMinimum          , concentrationParameter, &
         &                                                          alphaParameter
    type            (integrator                 )                :: integratorPotential    , integratorKinetic     , &
         &                                                          integratorJeans

    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.self%energyTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired < self%energyTableAlphaMinimum .or. alphaRequired > self%energyTableAlphaMaximum) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%energyTableAlphaMinimum =min(self%energyTableAlphaMinimum,0.9d0*alphaRequired)
       self%energyTableAlphaMaximum =max(self%energyTableAlphaMaximum,1.1d0*alphaRequired)
    end if
    ! Check for concentration below minimum tabulated value.
    if (concentrationRequired < self%energyTableConcentrationMinimum .or. concentrationRequired > self%energyTableConcentrationMaximum) then
       makeTable=.true.
       self%energyTableConcentrationMinimum=min(self%energyTableConcentrationMinimum,0.5d0*self%energyTableConcentrationMinimum)
       self%energyTableConcentrationMaximum=max(self%energyTableConcentrationMaximum,2.0d0*self%energyTableConcentrationMaximum)
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Allocate arrays to the appropriate sizes.
       self%energyTableAlphaCount        =int(     (self%energyTableAlphaMaximum        -self%energyTableAlphaMinimum        ) &
            &                                 *dble(einastoEnergyTableAlphaPointsPerUnit                                     ) &
            &                                )+1
       self%energyTableConcentrationCount=int(log10(self%energyTableConcentrationMaximum/self%energyTableConcentrationMinimum) &
            &                                 *dble(einastoEnergyTableConcentrationPointsPerDecade                           ) &
            &                                )+1
       if (allocated(self%energyTableAlpha        )) call deallocateArray(self%energyTableAlpha        )
       if (allocated(self%energyTableConcentration)) call deallocateArray(self%energyTableConcentration)
       if (allocated(self%energyTable             )) call deallocateArray(self%energyTable             )
       call allocateArray(self%energyTableAlpha        ,[                                   self%energyTableAlphaCount])
       call allocateArray(self%energyTableConcentration,[self%energyTableConcentrationCount                           ])
       call allocateArray(self%energyTable             ,[self%energyTableConcentrationCount,self%energyTableAlphaCount])
       ! Create ranges of alpha and concentration.
       self%energyTableAlpha        =Make_Range(self%energyTableAlphaMinimum        ,self%energyTableAlphaMaximum        , &
            &                                   self%energyTableAlphaCount          ,rangeType=rangeTypeLinear     )
       self%energyTableConcentration=Make_Range(self%energyTableConcentrationMinimum,self%energyTableConcentrationMaximum, &
            &                                   self%energyTableConcentrationCount  ,rangeType=rangeTypeLogarithmic)
       ! Tabulate the radius vs. specific angular momentum relation.
       integratorPotential=integrator(einastoPotentialEnergyIntegrand,toleranceRelative=1.0d-3)
       integratorKinetic  =integrator(einastoKineticEnergyIntegrand  ,toleranceRelative=1.0d-3)
       integratorJeans    =integrator(einastoJeansEquationIntegrand  ,toleranceRelative=1.0d-3)
       do iAlpha=1,self%energyTableAlphaCount
          alpha=self%energyTableAlpha(iAlpha)
          do iConcentration=1,self%energyTableConcentrationCount
             concentration=self%energyTableConcentration(iConcentration)

             ! Compute the potential energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             potentialEnergyIntegral=integratorPotential%integrate(radiusMinimum,radiusMaximum)
             potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)

             ! Compute the velocity dispersion at the virial radius.
             radiusMinimum         =        concentration
             radiusMaximum         =100.0d0*concentration
             concentrationParameter=        concentration
             alphaParameter        =alpha
             jeansEquationIntegral=integratorJeans%integrate(radiusMinimum,radiusMaximum)

             ! Compute the kinetic energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             kineticEnergyIntegral=integratorKinetic%integrate(radiusMinimum,radiusMaximum)
             kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)

             ! Compute the total energy.
             self%energyTable(iConcentration,iAlpha)=(potentialEnergy+kineticEnergy)*concentration

          end do
       end do
       ! Build interpolators.
       if (allocated(self%energyTableConcentrationInterpolator)) deallocate(self%energyTableConcentrationInterpolator)
       if (allocated(self%energyTableAlphaInterpolator        )) deallocate(self%energyTableAlphaInterpolator        )
       allocate(self%energyTableConcentrationInterpolator)
       allocate(self%energyTableAlphaInterpolator        )
       self%energyTableConcentrationInterpolator=interpolator(self%energyTableConcentration)
       self%energyTableAlphaInterpolator        =interpolator(self%energyTableAlpha        )
       ! Flag that the table is now initialized.
       self%energyTableInitialized=.true.
    end if
    return

  contains

    double precision function einastoPotentialEnergyIntegrand(radius)
      !!{
      Integrand for Einasto profile potential energy.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      einastoPotentialEnergyIntegrand=(self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter)/radius)**2
      return
    end function einastoPotentialEnergyIntegrand

    double precision function einastoKineticEnergyIntegrand(radius)
      !!{
      Integrand for Einasto profile kinetic energy.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      einastoKineticEnergyIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter) &
           &                        *self%densityScaleFree    (radius,concentrationParameter,alphaParameter) &
           &                        *radius
      return
    end function einastoKineticEnergyIntegrand

    double precision function einastoJeansEquationIntegrand(radius)
      !!{
      Integrand for Einasto profile Jeans equation.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      einastoJeansEquationIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter) &
           &                        *self%densityScaleFree    (radius,concentrationParameter,alphaParameter) &
           &                        /radius**2
      return
    end function einastoJeansEquationIntegrand

  end subroutine einastoEnergyTableMake

  double precision function einastoEnclosedMassScaleFree(self,radius,concentration,alpha)
    !!{
    Returns the enclosed mass (in units of the virial mass) in an Einasto dark matter profile with given {\normalfont \ttfamily concentration} at the
    given {\normalfont \ttfamily radius} (given in units of the scale radius).
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Complementary
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alpha, concentration, radius
    !$GLC attributes unused :: self

    einastoEnclosedMassScaleFree=+Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius       **alpha/alpha) &
         &                       /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    return
  end function einastoEnclosedMassScaleFree

  double precision function einastoDensityScaleFree(self,radius,concentration,alpha)
    !!{
    Returns the density (in units such that the virial mass and scale length are unity) in an Einasto dark matter profile with
    given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    !!}
    use :: Gamma_Functions         , only : Gamma_Function, Gamma_Function_Incomplete_Complementary
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alpha               , concentration, radius
    double precision                                             :: densityNormalization
    !$GLC attributes unused :: self

    densityNormalization= (alpha/4.0d0/Pi)                                                                      &
         &               *    ((2.0d0/alpha)                   **(3.0d0/alpha)                                ) &
         &               *exp(-2.0d0/alpha                                                                    ) &
         &               /Gamma_Function                         (3.0d0/alpha                                 ) &
         &               /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    einastoDensityScaleFree=densityNormalization*exp(-(2.0d0/alpha)*(radius**alpha-1.0d0))
    return
  end function einastoDensityScaleFree

  double precision function einastoPotentialScaleFree(self,radius,concentration,alpha)
    !!{
    Returns the gravitational potential (in units where the virial mass and scale radius are unity) in an Einasto dark matter
    profile with given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont
    \ttfamily radius} (given in units of the scale radius). Uses the results from \cite{retana-montenegro_analytical_2012},
    their equations (19) and (20).
    !!}
    use :: Gamma_Functions, only : Gamma_Function, Gamma_Function_Incomplete, Gamma_Function_Incomplete_Complementary
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alpha , concentration, &
         &                                                          radius
    !$GLC attributes unused :: self

    if (radius <= 0.0d0) then
       einastoPotentialScaleFree=-((2.0d0/alpha)**(1.0d0/alpha))                                                          &
            &                    *  Gamma_Function                         (2.0d0/alpha                                 ) &
            &                    /  Gamma_Function                         (3.0d0/alpha                                 ) &
            &                    /  Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)

    else
       einastoPotentialScaleFree=-1.0d0                                                                                   &
            &                    /radius                                                                                  &
            &                    *(                                                                                       &
            &                      +Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius       **alpha/alpha) &
            &                      +(                                                                                     &
            &                        +radius                                                                              &
            &                        *(2.0d0/alpha)**(1.0d0/alpha)                                                        &
            &                       )                                                                                     &
            &                      *Gamma_Function_Incomplete              (2.0d0/alpha,2.0d0*radius       **alpha/alpha) &
            &                      *Gamma_Function                         (2.0d0/alpha                                 ) &
            &                      /Gamma_Function                         (3.0d0/alpha                                 ) &
            &                     )                                                                                       &
            &                    /  Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    end if
    return
  end function einastoPotentialScaleFree

  double precision function einastoKSpace(self,node,wavenumber)
    !!{
    Returns the Fourier transform of the Einasto density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$)).
    !!}
    use            :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    class           (darkMatterProfileDMOEinasto   )                , intent(inout)          :: self
    type            (treeNode                      )                , intent(inout), target  :: node
    double precision                                                , intent(in   )          :: wavenumber
    class           (nodeComponentDarkMatterProfile)                               , pointer :: darkMatterProfile
    integer         (c_size_t                      ), dimension(0:1)                         :: jAlpha                     , jConcentration
    double precision                                , dimension(0:1)                         :: hAlpha                     , hConcentration
    integer                                                                                  :: iAlpha                     , iConcentration
    double precision                                                                         :: alpha                      , scaleRadius        , &
         &                                                                                      virialRadiusOverScaleRadius, wavenumberScaleFree

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =darkMatterProfile%scale()
    alpha                      =darkMatterProfile%shape()
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    wavenumberScaleFree        =wavenumber*scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call self%fourierProfileTableMake(wavenumberScaleFree,virialRadiusOverScaleRadius,alpha)

    ! Get interpolating factors in alpha.
    call self%fourierProfileTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
    jAlpha(1)=jAlpha(0)+1

    ! Get interpolating factors in concentration.
    call self%fourierProfileTableConcentrationInterpolator%linearFactors(virialRadiusOverScaleRadius,jConcentration(0),hConcentration)
    jConcentration(1)=jConcentration(0)+1

    ! Find the Fourier profile by interpolation.
    einastoKSpace=0.0d0
    do iAlpha=0,1
       do iConcentration=0,1
          einastoKSpace=+einastoKSpace                                                                                                                                             &
               &        +self%fourierProfileTableWavenumberInterpolator%interpolate(wavenumberScaleFree,self%fourierProfileTable(:,jConcentration(iConcentration),jAlpha(iAlpha))) &
               &        *                                                                                                                                         hAlpha(iAlpha)   &
               &        *                                                                                                          hConcentration(iConcentration)
       end do
    end do
    return
  end function einastoKSpace

  subroutine einastoFourierProfileTableMake(self,wavenumberRequired,concentrationRequired,alphaRequired)
    !!{
    Create a tabulation of the Fourier transform of Einasto profiles as a function of their $\alpha$ parameter and
    dimensionless wavenumber.
    !!}
    use            :: Display              , only : displayCounter         , displayCounterClear  , displayIndent       , displayUnindent, &
          &                                         verbosityLevelInfo     , verbosityLevelWorking
    use            :: Galacticus_Error     , only : Galacticus_Error_Report, errorStatusSuccess
    use, intrinsic :: ISO_C_Binding        , only : c_size_t
    use            :: Memory_Management    , only : allocateArray          , deallocateArray
    use            :: Numerical_Integration, only : integrator
    use            :: Numerical_Ranges     , only : Make_Range             , rangeTypeLinear      , rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alphaRequired                , concentrationRequired, wavenumberRequired
    double precision                             , parameter     :: profileTruncateLevel  =6.0d-6
    integer                                                      :: iAlpha                       , iConcentration       , iWavenumber        , &
         &                                                          percentage                   , errorStatus
    logical                                                      :: makeTable
    double precision                                             :: alpha                        , concentration        , radiusMaximum      , &
         &                                                          radiusMinimum                , wavenumber           , wavenumberParameter, &
         &                                                          concentrationParameter       , alphaParameter
    type            (integrator                 )                :: integrator_
    character       (len=12                     )                :: label
    type            (varying_string             )                :: message

    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.self%fourierProfileTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired         < self%fourierProfileTableAlphaMinimum         .or. alphaRequired         > self%fourierProfileTableAlphaMaximum        ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%fourierProfileTableAlphaMinimum        =min(self%fourierProfileTableAlphaMinimum     ,0.9d0*alphaRequired           )
       self%fourierProfileTableAlphaMaximum        =max(self%fourierProfileTableAlphaMaximum     ,1.1d0*alphaRequired           )
    end if
    ! Check for concentration out of range.
    if (concentrationRequired < self%fourierProfileTableConcentrationMinimum .or. concentrationRequired > self%fourierProfileTableConcentrationMaximum ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%fourierProfileTableConcentrationMinimum=min(self%fourierProfileTableConcentrationMinimum,0.5d0*concentrationRequired)
       self%fourierProfileTableConcentrationMaximum=max(self%fourierProfileTableConcentrationMaximum,2.0d0*concentrationRequired)
    end if
    ! Check for wavenumber below minimum tabulated value.
    if (wavenumberRequired    < self%fourierProfileTableWavenumberMinimum    .or. wavenumberRequired    > self%fourierProfileTableWavenumberMaximum    ) then
       makeTable=.true.
       self%fourierProfileTableWavenumberMinimum   =min(self%fourierProfileTableWavenumberMinimum   ,0.5d0*wavenumberRequired   )
       self%fourierProfileTableWavenumberMaximum   =max(self%fourierProfileTableWavenumberMaximum   ,2.0d0*wavenumberRequired   )
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Display a message.
       call displayIndent('Constructing Einasto profile Fourier transform lookup table...',verbosityLevelInfo)
       ! Allocate arrays to the appropriate sizes.
       self%fourierProfileTableAlphaCount        =int(      (self%fourierProfileTableAlphaMaximum        -self%fourierProfileTableAlphaMinimum        ) &
            &*dble(einastoFourierProfileTableAlphaPointsPerUnit          ))+1
       self%fourierProfileTableConcentrationCount=int(log10(self%fourierProfileTableConcentrationMaximum/self%fourierProfileTableConcentrationMinimum) &
            &*dble(einastoFourierProfileTableConcentrationPointsPerDecade))+1
       self%fourierProfileTableWavenumberCount   =int(log10(self%fourierProfileTableWavenumberMaximum   /self%fourierProfileTableWavenumberMinimum   ) &
            &*dble(einastoFourierProfileTableWavenumberPointsPerDecade   ))+1
       if (allocated(self%fourierProfileTableAlpha        )) call deallocateArray(self%fourierProfileTableAlpha        )
       if (allocated(self%fourierProfileTableConcentration)) call deallocateArray(self%fourierProfileTableConcentration)
       if (allocated(self%fourierProfileTableWavenumber   )) call deallocateArray(self%fourierProfileTableWavenumber   )
       if (allocated(self%fourierProfileTable             )) call deallocateArray(self%fourierProfileTable             )
       call allocateArray(self%fourierProfileTableAlpha        ,[                                                                                   self%fourierProfileTableAlphaCount])
       call allocateArray(self%fourierProfileTableConcentration,[                                        self%fourierProfileTableConcentrationCount                                   ])
       call allocateArray(self%fourierProfileTableWavenumber   ,[self%fourierProfileTableWavenumberCount                                                                              ])
       call allocateArray(self%fourierProfileTable             ,[self%fourierProfileTableWavenumberCount,self%fourierProfileTableConcentrationCount,self%fourierProfileTableAlphaCount])
       ! Create ranges of alpha and wavenumber.
       self%fourierProfileTableAlpha        =Make_Range(self%fourierProfileTableAlphaMinimum        ,self%fourierProfileTableAlphaMaximum        , &
            &                                           self%fourierProfileTableAlphaCount          ,rangeType=rangeTypeLinear     )
       self%fourierProfileTableConcentration=Make_Range(self%fourierProfileTableConcentrationMinimum,self%fourierProfileTableConcentrationMaximum, &
            &                                           self%fourierProfileTableConcentrationCount  ,rangeType=rangeTypeLogarithmic)
       self%fourierProfileTableWavenumber   =Make_Range(self%fourierProfileTableWavenumberMinimum   ,self%fourierProfileTableWavenumberMaximum   , &
            &                                           self%fourierProfileTableWavenumberCount     ,rangeType=rangeTypeLogarithmic)
       ! Tabulate the Fourier profile.
       integrator_=integrator(einastoFourierProfileIntegrand,toleranceRelative=1.0d-3,intervalsMaximum=1000_c_size_t)
       do iAlpha=1,self%fourierProfileTableAlphaCount
          alpha=self%fourierProfileTableAlpha(iAlpha)
          do iConcentration=1,self%fourierProfileTableConcentrationCount
             concentration=self%fourierProfileTableConcentration(iConcentration)

             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*self%fourierProfileTableConcentrationCount+iConcentration-1       ) &
                  &                /dble(self%fourierProfileTableAlphaCount*self%fourierProfileTableConcentrationCount) &
                  &        )
             call displayCounter(percentage,iAlpha == 1 .and. iConcentration == 1,verbosityLevelWorking)

             do iWavenumber=1,self%fourierProfileTableWavenumberCount
                ! If the Fourier profile has fallen below some minimal level, simply truncate to zero to avoid numerical
                ! integration problems.
                if (iWavenumber > 1 .and. self%fourierProfileTable(iWavenumber-1,iConcentration,iAlpha) <= profileTruncateLevel) then
                   self%fourierProfileTable(iWavenumber,iConcentration,iAlpha)=0.0d0
                else
                   wavenumber=self%fourierProfileTableWavenumber(iWavenumber)
                   ! Compute the Fourier profile.
                   radiusMinimum         =0.0d0
                   radiusMaximum         =concentration
                   wavenumberParameter   =wavenumber
                   alphaParameter        =alpha
                   concentrationParameter=concentration
                   self%fourierProfileTable(iWavenumber,iConcentration,iAlpha)=integrator_%integrate(radiusMinimum,radiusMaximum,status=errorStatus)
                   if (errorStatus /= errorStatusSuccess) then
                      message="Integration of Einasto profile Fourier transform failed at:"//char(10)
                      write (label,'(e12.6)') wavenumber
                      message=message//"   wavenumber: k="//trim(adjustl(label))//"Mpc"//char(10)
                      if (iWavenumber == 1) then
                         message=message//"   no previous tabulated point"
                      else
                         write (label,'(e12.6)') self%fourierProfileTable(iWavenumber-1,iConcentration,iAlpha)
                         message=message//"   value at previous tabulated point was "//trim(adjustl(label))
                      end if
                      call Galacticus_Error_Report(message//{introspection:location})
                   end if
                end if
             end do
          end do
       end do
       call displayCounterClear(verbosityLevelWorking)
       ! Build interpolators.
       if (allocated(self%fourierProfileTableWavenumberInterpolator   )) deallocate(self%fourierProfileTableWavenumberInterpolator   )
       if (allocated(self%fourierProfileTableAlphaInterpolator        )) deallocate(self%fourierProfileTableAlphaInterpolator        )
       if (allocated(self%fourierProfileTableConcentrationInterpolator)) deallocate(self%fourierProfileTableConcentrationInterpolator)
       allocate(self%fourierProfileTableWavenumberInterpolator   )
       allocate(self%fourierProfileTableAlphaInterpolator        )
       allocate(self%fourierProfileTableConcentrationInterpolator)
       self%fourierProfileTableWavenumberInterpolator   =interpolator(self%fourierProfileTableWavenumber)
       self%fourierProfileTableAlphaInterpolator        =interpolator(self%fourierProfileTableAlpha)
       self%fourierProfileTableConcentrationInterpolator=interpolator(self%fourierProfileTableConcentration)
       ! Flag that the table is now initialized.
       self%fourierProfileTableInitialized=.true.
       ! Display a message.
       call displayUnindent('done',verbosityLevelInfo)
    end if
    return

  contains

    double precision function einastoFourierProfileIntegrand(radius)
      !!{
      Integrand for Einasto Fourier profile.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius

      einastoFourierProfileIntegrand=4.0d0*Pi*radius*sin(wavenumberParameter*radius)*self%densityScaleFree(radius&
           &,concentrationParameter,alphaParameter)/wavenumberParameter
      return
    end function einastoFourierProfileIntegrand

  end subroutine einastoFourierProfileTableMake

  double precision function einastoFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the Einasto density profile at the specified {\normalfont \ttfamily time} (given in Gyr).
    !!}
    use            :: Galacticus_Nodes                , only : nodeComponentBasic                     , nodeComponentDarkMatterProfile , treeNode
    use            :: Gamma_Functions                 , only : Gamma_Function_Incomplete_Complementary
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr                , gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) , target :: self
    type            (treeNode                      ), intent(inout) , target :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentBasic            ), pointer                :: basic
    class           (nodeComponentDarkMatterProfile), pointer                :: darkMatterProfile
    integer         (c_size_t                      ), dimension(0:1)         :: jAlpha
    double precision                                , dimension(0:1)         :: hAlpha
    integer                                                                  :: iAlpha
    double precision                                                         :: alpha            , freefallTimeScaleFree      , &
         &                                                                      radiusScale      , timeScale                  , &
         &                                                                      velocityScale    , virialRadiusOverScaleRadius

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       einastoFreefallRadius=0.0d0
       return
    end if

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get the shape parameter.
    alpha        =darkMatterProfile%shape()

    ! Get the scale radius.
    radiusScale                =darkMatterProfile%scale()
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/radiusScale

    ! Get the velocity scale.
    velocityScale=sqrt(gravitationalConstantGalacticus*basic%mass()/radiusScale)

    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                                                                                   &
         &    *radiusScale                                                                                               &
         &    /velocityScale                                                                                             &
         &    *sqrt(Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha)) &
         &    /sqrt(Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0                                   /alpha))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    ! Get interpolating factors in alpha.
    call self%freefallRadiusTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
    jAlpha(1)=jAlpha(0)+1

    einastoFreefallRadius=0.0d0
    do iAlpha=0,1
       einastoFreefallRadius=+einastoFreefallRadius                                                                                                        &
            &                +self%freefallRadiusTableRadiusInterpolator(jAlpha(iAlpha))%interpolate(freefallTimeScaleFree,self%freefallRadiusTableRadius) &
            &                *                                           hAlpha(iAlpha)
    end do
    einastoFreefallRadius=einastoFreefallRadius*radiusScale
    return
  end function einastoFreefallRadius

  double precision function einastoFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the Einasto density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    use            :: Galacticus_Nodes                , only : nodeComponentBasic                     , nodeComponentDarkMatterProfile , treeNode
    use            :: Gamma_Functions                 , only : Gamma_Function_Incomplete_Complementary
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr                , gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOEinasto   ), intent(inout) , target :: self
    type            (treeNode                      ), intent(inout) , target :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentBasic            ), pointer                :: basic
    class           (nodeComponentDarkMatterProfile), pointer                :: darkMatterProfile
    integer         (c_size_t                      ), dimension(0:1)         :: jAlpha
    double precision                                , dimension(0:1)         :: hAlpha
    integer                                                                  :: iAlpha
    double precision                                                         :: alpha            , freefallTimeScaleFree      , &
         &                                                                      radiusScale      , timeScale                  , &
         &                                                                      velocityScale    , virialRadiusOverScaleRadius

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       einastoFreefallRadiusIncreaseRate=0.0d0
       return
    end if

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get the shape parameter.
    alpha        =darkMatterProfile%shape()

    ! Get the scale radius.
    radiusScale                =darkMatterProfile%scale()
    virialRadiusOverScaleRadius=self%darkMatterHaloScale_%radiusVirial(node)/radiusScale

    ! Get the velocity scale.
    velocityScale=sqrt(gravitationalConstantGalacticus*basic%mass()/radiusScale)

    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                                                                                   &
         &    *radiusScale                                                                                               &
         &    /velocityScale                                                                                             &
         &    *sqrt(Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha)) &
         &    /sqrt(Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0                                   /alpha))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    ! Get interpolating factors in alpha.
    call self%freefallRadiusTableAlphaInterpolator%linearFactors(alpha,jAlpha(0),hAlpha)
    jAlpha(1)=jAlpha(0)+1

    einastoFreefallRadiusIncreaseRate=0.0d0
    do iAlpha=0,1
       einastoFreefallRadiusIncreaseRate=+einastoFreefallRadiusIncreaseRate                                                                                           &
            &                            +self%freefallRadiusTableRadiusInterpolator(jAlpha(iAlpha))%derivative(freefallTimeScaleFree,self%freefallRadiusTableRadius) &
            &                            *                                           hAlpha(iAlpha)
    end do
    einastoFreefallRadiusIncreaseRate=einastoFreefallRadiusIncreaseRate*radiusScale/timeScale
    return
  end function einastoFreefallRadiusIncreaseRate

  subroutine einastoFreefallTabulate(self,freefallTimeScaleFree,alphaRequired)
    !!{
    Tabulates the freefall time vs. freefall radius for Einasto halos.
    !!}
    use :: Display          , only : displayCounter, displayIndent  , displayUnindent     , verbosityLevelWorking
    use :: Memory_Management, only : allocateArray , deallocateArray
    use :: Numerical_Ranges , only : Make_Range    , rangeTypeLinear, rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alphaRequired, freefallTimeScaleFree
    logical                                                      :: retabulate
    integer                                                      :: iAlpha       , iRadius              , percentage
    double precision                                             :: alpha

    retabulate=.not.self%freefallRadiusTableInitialized
    ! If the table has not yet been made, compute and store the freefall times corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMinimum,alphaRequired)
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMaximum,alphaRequired)
    end if
    do while (freefallTimeScaleFree < self%freefallTimeMinimum)
       self%freefallRadiusTableRadiusMinimum=0.5d0*self%freefallRadiusTableRadiusMinimum
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMinimum,alphaRequired)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > self%freefallTimeMaximum)
       self%freefallRadiusTableRadiusMaximum=2.0d0*self%freefallRadiusTableRadiusMaximum
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMaximum,alphaRequired)
       retabulate=.true.
    end do
    ! Check for alpha out of range.
    if (alphaRequired < self%freefallRadiusTableAlphaMinimum .or. alphaRequired > self%freefallRadiusTableAlphaMaximum) then
       retabulate=.true.
       ! Compute the range of tabulation.
       self%freefallRadiusTableAlphaMinimum=min(self%freefallRadiusTableAlphaMinimum,0.9d0*alphaRequired)
       self%freefallRadiusTableAlphaMaximum=max(self%freefallRadiusTableAlphaMaximum,1.1d0*alphaRequired)
    end if

    if (retabulate) then
       ! Display a message.
       call displayIndent('Constructing Einasto profile freefall radius lookup table...',verbosityLevelWorking)
       ! Decide how many points to tabulate and allocate table arrays.
       self%freefallRadiusTableRadiusCount=int(log10(self%freefallRadiusTableRadiusMaximum/self%freefallRadiusTableRadiusMinimum)*dble(einastoFreefallRadiusTableRadiusPointsPerDecade))+1
       self%freefallRadiusTableAlphaCount =int(     (self%freefallRadiusTableAlphaMaximum -self%freefallRadiusTableAlphaMinimum )*dble(einastoFreefallRadiusTableAlphaPointsPerUnit   ))+1
       if (allocated(self%freefallRadiusTableRadius)) then
          call deallocateArray(self%freefallRadiusTableAlpha )
          call deallocateArray(self%freefallRadiusTableRadius)
          call deallocateArray(self%freefallRadiusTable      )
       end if
       call allocateArray(self%freefallRadiusTableAlpha ,[                                    self%freefallRadiusTableAlphaCount])
       call allocateArray(self%freefallRadiusTableRadius,[self%freefallRadiusTableRadiusCount                                   ])
       call allocateArray(self%freefallRadiusTable      ,[self%freefallRadiusTableRadiusCount,self%freefallRadiusTableAlphaCount])
       ! Create a range of radii and alpha.
       self%freefallRadiusTableAlpha =Make_Range(self%freefallRadiusTableAlphaMinimum ,self%freefallRadiusTableAlphaMaximum ,self%freefallRadiusTableAlphaCount ,rangeType=rangeTypeLinear     )
       self%freefallRadiusTableRadius=Make_Range(self%freefallRadiusTableRadiusMinimum,self%freefallRadiusTableRadiusMaximum,self%freefallRadiusTableRadiusCount,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and alpha and populate tables.
       do iAlpha=1,self%freefallRadiusTableAlphaCount
          alpha=self%freefallRadiusTableAlpha(iAlpha)
          do iRadius=1,self%freefallRadiusTableRadiusCount
             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*self%freefallRadiusTableRadiusCount+iRadius-1              ) &
                  &                /dble(self%freefallRadiusTableAlphaCount*self%freefallRadiusTableRadiusCount) &
                  &        )
             call displayCounter(percentage,iAlpha == 1 .and. iRadius == 1,verbosityLevelWorking)
             ! Compute the freefall radius.
             self%freefallRadiusTable(iRadius,iAlpha)=self%freefallTimeScaleFree(self%freefallRadiusTableRadius(iRadius),alpha)
          end do
       end do
       ! Build interpolators.
       if (allocated(self%freefallRadiusTableAlphaInterpolator )) deallocate(self%freefallRadiusTableAlphaInterpolator )
       if (allocated(self%freefallRadiusTableRadiusInterpolator)) deallocate(self%freefallRadiusTableRadiusInterpolator)
       allocate(self%freefallRadiusTableAlphaInterpolator                                     )
       allocate(self%freefallRadiusTableRadiusInterpolator(self%freefallRadiusTableAlphaCount))
       self   %freefallRadiusTableAlphaInterpolator         =interpolator(self%freefallRadiusTableAlpha          )
       do iAlpha=1,self%freefallRadiusTableAlphaCount
          self%freefallRadiusTableRadiusInterpolator(iAlpha)=interpolator(self%freefallRadiusTable     (:,iAlpha))
       end do
       ! Store the minimum and maximum tabulated freefall times across all alpha values.
       self%freefallTimeMinimum=maxval(self%freefallRadiusTable(                                  1,:))
       self%freefallTimeMaximum=minval(self%freefallRadiusTable(self%freefallRadiusTableRadiusCount,:))
       ! Display a message.
       call displayUnindent('...done',verbosityLevelWorking)
       ! Specify that tabulation has been made.
       self%freefallRadiusTableInitialized=.true.
    end if
    return
  end subroutine einastoFreefallTabulate

  double precision function einastoFreefallTimeScaleFree(self,radius,alpha)
    !!{
    Compute the freefall time in a scale-free Einasto halo.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout) :: self
    double precision                             , intent(in   ) :: alpha         , radius
    type            (integrator                 )                :: integrator_
    double precision                                             :: radiusStart   , radiusEnd, &
         &                                                          alphaParameter

    radiusStart                 =radius
    radiusEnd                   =0.0d0
    alphaParameter              =alpha
    integrator_                 =integrator           (einastoFreefallTimeScaleFreeIntegrand,toleranceRelative=1.0d-3)
    einastoFreefallTimeScaleFree=integrator_%integrate(radiusEnd                            ,radiusStart             )
    return

  contains

    double precision function einastoFreefallTimeScaleFreeIntegrand(radius)
      !!{
      Integrand function used for finding the free-fall time in Einasto halos.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      einastoFreefallTimeScaleFreeIntegrand=+1.0d0                                                             &
           &                                /sqrt(+2.0d0                                                       &
           &                                      *(                                                           &
           &                                        +self%potentialScaleFree(radiusStart,1.0d0,alphaParameter) &
           &                                        -self%potentialScaleFree(radius     ,1.0d0,alphaParameter) &
           &                                       )                                                           &
           &                                     )
      return
    end function einastoFreefallTimeScaleFreeIntegrand

  end function einastoFreefallTimeScaleFree

  double precision function einastoRadiusEnclosingDensity(self,node,density)
    !!{
    Implementation of function to compute the radius enclosing a given density for Einasto dark matter halo profiles. This
    function uses a numerical root finder to find the enclosing radius---this is likely not the most efficient solution\ldots
    !!}
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout), target :: self
    type            (treeNode                   ), intent(inout), target :: node
    double precision                             , intent(in   )         :: density

    einastoDensityEnclosed        =  density
    einastoNode                   => node
    einastoSelf                   => self
    einastoRadiusEnclosingDensity =  self%finderEnclosedDensity%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    return
  end function einastoRadiusEnclosingDensity

  double precision function einastoRadiusEnclosingDensityRoot(radius)
    !!{
    Root function used in finding the radius enclosing a given density in Einasto profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    einastoRadiusEnclosingDensityRoot=+einastoDensityEnclosed                       &
         &                            -einastoSelf%enclosedMass(einastoNode,radius) &
         &                            *3.0d0                                        &
         &                            /4.0d0                                        &
         &                            /Pi                                           &
         &                            /radius**3
    return
  end function einastoRadiusEnclosingDensityRoot

  double precision function einastoRadialVelocityDispersionScaleFree(self,radius,alpha)
    !!{
    Compute the radial velocity dispersion in a scale-free Einasto halo.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout)            :: self
    double precision                             , intent(in   )            :: radius              , alpha
    double precision                                            , parameter :: radiusTiny   =1.0d-9, radiusLarge  =5.0d3
    double precision                                                        :: radiusMinimum       , radiusMaximum
    type            (integrator                 )                           :: integrator_
    !$GLC attributes unused :: self

    radiusMinimum=max(       radius,radiusTiny )
    radiusMaximum=max(10.0d0*radius,radiusLarge)
    integrator_=integrator(einastoJeansEquationIntegrand,toleranceRelative=1.0d-6)
    einastoRadialVelocityDispersionScaleFree=sqrt(                                                    &
         &                                        +integrator_%integrate(radiusMinimum,radiusMaximum) &
         &                                        /exp(-(2.0d0/alpha)*(radius**alpha-1.0d0))          &
         &                                       )
    return

  contains

    double precision function einastoJeansEquationIntegrand(radius)
      !!{
      Integrand for Einasto drak matter profile Jeans equation.
      !!}
      use :: Gamma_Functions, only : Gamma_Function_Incomplete_Complementary
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         einastoJeansEquationIntegrand=+Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius**alpha/alpha) &
              &                        *exp(-(2.0d0/alpha)*(radius**alpha-1.0d0))                                      &
              &                        /radius**2
      else
         einastoJeansEquationIntegrand=0.0d0
      end if
      return
    end function einastoJeansEquationIntegrand

  end function einastoRadialVelocityDispersionScaleFree

  subroutine einastoRadialVelocityDispersionTabulate(self,radius,alphaRequired)
    !!{
    Tabulates the radial velocity dispersion vs. radius for Einasto halos.
    !!}
    use :: Display          , only : displayCounter, displayIndent  , displayUnindent     , verbosityLevelWorking
    use :: Memory_Management, only : allocateArray , deallocateArray
    use :: Numerical_Ranges , only : Make_Range    , rangeTypeLinear, rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOEinasto), intent(inout)           :: self
    double precision                             , intent(in   ), optional :: radius    , alphaRequired
    logical                                                                :: retabulate
    integer                                                                :: iAlpha    , iRadius      , percentage
    double precision                                                       :: alpha

    retabulate=.not.self%radialVelocityDispersionTableInitialized
    if (present(radius)) then
       do while (radius < self%radialVelocityDispersionRadiusMinimum)
          self%radialVelocityDispersionRadiusMinimum=0.5d0*self%radialVelocityDispersionRadiusMinimum
          retabulate=.true.
       end do
       do while (radius > self%radialVelocityDispersionRadiusMaximum)
          self%radialVelocityDispersionRadiusMaximum=2.0d0*self%radialVelocityDispersionRadiusMaximum
          retabulate=.true.
       end do
    end if
    if (present(alphaRequired)) then
       ! Check for alpha out of range.
       if (alphaRequired < self%radialVelocityDispersionAlphaMinimum .or. alphaRequired > self%radialVelocityDispersionAlphaMaximum) then
          retabulate=.true.
          ! Compute the range of tabulation.
          self%radialVelocityDispersionAlphaMinimum=min(self%radialVelocityDispersionAlphaMinimum,0.9d0*alphaRequired)
          self%radialVelocityDispersionAlphaMaximum=max(self%radialVelocityDispersionAlphaMaximum,1.1d0*alphaRequired)
       end if
    end if
    if (retabulate) then
       ! Display a message.
       call displayIndent('Constructing Einasto profile radial velocity dispersion lookup table...',verbosityLevelWorking)
       ! Decide how many points to tabulate and allocate table arrays.
       self%radialVelocityDispersionTableRadiusCount=int(log10(self%radialVelocityDispersionRadiusMaximum/self%radialVelocityDispersionRadiusMinimum) &
            &                                            *dble(einastoRadialVelocityDispersionTableRadiusPointsPerDecade                            ) &
            &                                           )+1
       self%radialVelocityDispersionTableAlphaCount =int(     (self%radialVelocityDispersionAlphaMaximum -self%radialVelocityDispersionAlphaMinimum ) &
            &                                            *dble(einastoRadialVelocityDispersionTableAlphaPointsPerUnit                               ) &
            &                                           )+1
       if (allocated(self%radialVelocityDispersionTableRadius)) then
          call deallocateArray(self%radialVelocityDispersionTableAlpha )
          call deallocateArray(self%radialVelocityDispersionTableRadius)
          call deallocateArray(self%radialVelocityDispersionTable      )
       end if
       call allocateArray(self%radialVelocityDispersionTableAlpha ,[                                              self%radialVelocityDispersionTableAlphaCount])
       call allocateArray(self%radialVelocityDispersionTableRadius,[self%radialVelocityDispersionTableRadiusCount                                             ])
       call allocateArray(self%radialVelocityDispersionTable      ,[self%radialVelocityDispersionTableRadiusCount,self%radialVelocityDispersionTableAlphaCount])
       ! Create a range of radii and alpha.
       self%radialVelocityDispersionTableAlpha =Make_Range(self%radialVelocityDispersionAlphaMinimum ,self%radialVelocityDispersionAlphaMaximum ,self%radialVelocityDispersionTableAlphaCount ,rangeType=rangeTypeLinear     )
       self%radialVelocityDispersionTableRadius=Make_Range(self%radialVelocityDispersionRadiusMinimum,self%radialVelocityDispersionRadiusMaximum,self%radialVelocityDispersionTableRadiusCount,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and alpha and populate tables.
       do iAlpha=1,self%radialVelocityDispersionTableAlphaCount
          alpha=self%radialVelocityDispersionTableAlpha(iAlpha)
          do iRadius=1,self%radialVelocityDispersionTableRadiusCount
             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*self%radialVelocityDispersionTableRadiusCount+iRadius-1                        ) &
                  &                /dble(self%radialVelocityDispersionTableAlphaCount*self%radialVelocityDispersionTableRadiusCount) &
                  &        )
             call displayCounter(percentage,iAlpha == 1 .and. iRadius == 1,verbosityLevelWorking)
             ! Compute the radial velocity dispersion.
             self%radialVelocityDispersionTable(iRadius,iAlpha)=self%radialVelocityDispersionScaleFree(self%radialVelocityDispersionTableRadius(iRadius),alpha)
          end do
       end do
       ! Build interpolators.
       if (allocated(self%radialVelocityDispersionTableAlphaInterpolator )) deallocate(self%radialVelocityDispersionTableAlphaInterpolator )
       if (allocated(self%radialVelocityDispersionTableRadiusInterpolator)) deallocate(self%radialVelocityDispersionTableRadiusInterpolator)
       allocate(self%radialVelocityDispersionTableAlphaInterpolator )
       allocate(self%radialVelocityDispersionTableRadiusInterpolator)
       self%radialVelocityDispersionTableAlphaInterpolator =interpolator(self%radialVelocityDispersionTableAlpha )
       self%radialVelocityDispersionTableRadiusInterpolator=interpolator(self%radialVelocityDispersionTableRadius)
       ! Display a message.
       call displayUnindent('...done',verbosityLevelWorking)
       ! Specify that tabulation has been made.
       self%radialVelocityDispersionTableInitialized=.true.
    end if
    return
  end subroutine einastoRadialVelocityDispersionTabulate
