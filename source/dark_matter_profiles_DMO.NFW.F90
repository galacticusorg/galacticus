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

  !% An implementation of \cite{navarro_universal_1997} dark matter halo profiles.

  use :: Kind_Numbers, only : kind_int8
  use :: Tables      , only : table1D  , table1DLogarithmicLinear

  !# <darkMatterProfileDMO name="darkMatterProfileDMONFW">
  !#  <description>
  !#   A dark matter profile DMO class which implements the \gls{nfw} density profile \citep{navarro_universal_1997}:
  !#   \begin{equation}
  !#     \rho_\mathrm{dark matter}(r) \propto \left({r\over r_\mathrm{s}}\right)^{-1} \left[1 + \left({r\over r_\mathrm{s}}\right)
  !#     \right]^{-2},
  !#   \end{equation}
  !#   normalized such that the total mass of the \gls{node} is enclosed with the virial radius and with the scale length
  !#   $r_\mathrm{s} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see \refPhysics{darkMatterProfileConcentration}).
  !#  </description>
  !# </darkMatterProfileDMO>
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMONFW
     !% A dark matter halo profile class implementing \cite{navarro_universal_1997} dark matter halos.
     private
     ! Minimum and maximum concentrations to tabulate.
     double precision                                        :: concentrationMinimum                   , concentrationMaximum
     ! Minimum and maximum radii to tabulate.
     double precision                                        :: freefallRadiusMinimum                  , radiusMinimum
     double precision                                        :: freefallRadiusMaximum                  , radiusMaximum
     double precision                                        :: freefallTimeMinimum                    , specificAngularMomentumMinimum
     double precision                                        :: freefallTimeMaximum                    , specificAngularMomentumMaximum
     double precision                                        :: enclosedDensityRadiusMinimum           , enclosedDensityRadiusMaximum
     double precision                                        :: enclosedDensityMinimum                 , enclosedDensityMaximum
     ! Tables of NFW properties.
     logical                                                 :: nfwFreefallTableInitialized            , nfwInverseTableInitialized         , &
          &                                                     nfwTableInitialized                    , nfwEnclosedDensityTableInitialized
     integer                                                 :: nfwFreefallTableNumberPoints           , nfwInverseTableNumberPoints        , &
          &                                                     nfwTableNumberPoints                   , nfwEnclosedDensityTableNumberPoints
     type            (table1DLogarithmicLinear)              :: nfwConcentrationTable
     ! Tables.
     type            (table1DLogarithmicLinear)              :: nfwFreeFall                            , nfwSpecificAngularMomentum         , &
          &                                                     nfwEnclosedDensity
     class           (table1D                 ), allocatable :: nfwFreefallInverse                     , nfwSpecificAngularMomentumInverse  , &
          &                                                     nfwEnclosedDensityInverse
     ! Module variables used in integrations.
     double precision                                        :: concentrationParameter                 , radiusStart
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )              :: lastUniqueID
     ! Record of whether or not quantities have been computed.
     logical                                                 :: specificAngularMomentumScalingsComputed, maximumVelocityComputed
     ! Stored values of computed quantities.
     double precision                                        :: specificAngularMomentumLengthScale     , specificAngularMomentumScale       , &
          &                                                     concentrationPrevious                  , nfwNormalizationFactorPrevious     , &
          &                                                     maximumVelocityPrevious                , enclosedDensityPrevious            , &
          &                                                     enclosingDensityRadiusPrevious         , densityScalePrevious               , &
          &                                                     enclosedMassPrevious                   , enclosingMassRadiusPrevious        , &
          &                                                     massScalePrevious                      , circularVelocityPrevious           , &
          &                                                     circularVelocityRadiusPrevious
   contains
     !# <methods>
     !#   <method description="Reset memoized calculations." method="calculationReset" />
     !#   <method description="Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="densityScaleFree" />
     !#   <method description="Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="enclosedMassScaleFree" />
     !#   <method description="Returns the density (in units of the virial mass per cubic scale radius) in an NFW dark matter profile with given {\normalfont \ttfamily concentration} which is enclosed a given radius (in units of the scale radius)." method="densityEnclosedByRadiusScaleFree" />
     !#   <method description="Returns the radial velocity dispersion (in units of the virial velocity) in an NFW dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="radialVelocityDispersionScaleFree" />
     !#   <method description="Tabulates the freefall time vs. freefall radius for NFW halos." method="freefallTabulate" />
     !#   <method description="Compute the freefall time in a scale-free NFW halo." method="freefallTimeScaleFree" />
     !#   <method description="Returns the total angular momentum in an NFW dark matter profile with given {\normalfont \ttfamily concentration}." method="angularMomentumScaleFree" />
     !#   <method description="Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion." method="inverseAngularMomentum" />
     !#   <method description="Computes the total energy of an NFW profile halo of given {\normalfont \ttfamily concentration}." method="profileEnergy" />
     !#   <method description="Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position {\normalfont \ttfamily radius} (in units of the scale radius) in an NFW profile." method="specificAngularMomentumScaleFree" />
     !#   <method description="Tabulate properties of the NFW halo profile which must be computed numerically." method="tabulate" />
     !#   <method description="Tabulate the density enclosed within a given radius for the NFW profile." method="enclosedDensityTabulate" />
     !# </methods>
     final     ::                                      nfwDestructor
     procedure :: autoHook                          => nfwAutoHook
     procedure :: calculationReset                  => nfwCalculationReset
     procedure :: density                           => nfwDensity
     procedure :: densityLogSlope                   => nfwDensityLogSlope
     procedure :: radialMoment                      => nfwRadialMoment
     procedure :: enclosedMass                      => nfwEnclosedMass
     procedure :: radiusEnclosingDensity            => nfwRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => nfwRadiusEnclosingMass
     procedure :: potential                         => nfwPotential
     procedure :: circularVelocity                  => nfwCircularVelocity
     procedure :: circularVelocityMaximum           => nfwCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => nfwRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => nfwRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => nfwRotationNormalization
     procedure :: energy                            => nfwEnergy
     procedure :: energyGrowthRate                  => nfwEnergyGrowthRate
     procedure :: kSpace                            => nfwKSpace
     procedure :: freefallRadius                    => nfwFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => nfwFreefallRadiusIncreaseRate
     procedure :: profileEnergy                     => nfwProfileEnergy
     procedure :: specificAngularMomentumScaleFree  => nfwSpecificAngularMomentumScaleFree
     procedure :: angularMomentumScaleFree          => nfwAngularMomentumScaleFree
     procedure :: enclosedMassScaleFree             => nfwEnclosedMassScaleFree
     procedure :: densityEnclosedByRadiusScaleFree  => nfwDensityEnclosedByRadiusScaleFree
     procedure :: densityScaleFree                  => nfwDensityScaleFree
     procedure :: radialVelocityDispersionScaleFree => nfwRadialVelocityDispersionScaleFree
     procedure :: tabulate                          => nfwTabulate
     procedure :: inverseAngularMomentum            => nfwInverseAngularMomentum
     procedure :: freefallTabulate                  => nfwFreefallTabulate
     procedure :: freefallTimeScaleFree             => nfwFreefallTimeScaleFree
     procedure :: enclosedDensityTabulate           => nfwEnclosedDensityTabulate
  end type darkMatterProfileDMONFW

  interface darkMatterProfileDMONFW
     !% Constructors for the {\normalfont \ttfamily nfw} dark matter halo profile class.
     module procedure nfwConstructorParameters
     module procedure nfwConstructorInternal
  end interface darkMatterProfileDMONFW

  ! Number of points per decade of concentration in NFW tabulations.
  integer, parameter   :: nfwTablePointsPerDecade               =100
  integer, parameter   :: nfwInverseTablePointsPerDecade        =100
  integer, parameter   :: nfwFreefallTablePointsPerDecade       =300
  integer, parameter   :: nfwEnclosedDensityTablePointsPerDecade=100
  ! Indices for tabulated quantities.
  integer, parameter   :: nfwConcentrationEnergyIndex           =  1, nfwConcentrationRotationNormalizationIndex=2

contains

  function nfwConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily nfw} dark matter halo profile class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMONFW )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_

    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=darkMatterProfileDMONFW(darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function nfwConstructorParameters

  function nfwConstructorInternal(darkMatterHaloScale_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily nfw} dark matter halo profile class.
    use :: Galacticus_Error, only : Galacticus_Component_List        , Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type (darkMatterProfileDMONFW )                        :: self
    class(darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    !# <constructorAssign variables="*darkMatterHaloScale_"/>

    self%concentrationPrevious             =-1.0d+0
    self%concentrationMinimum              = 1.0d+0
    self%concentrationMaximum              =20.0d+0
    self%freefallRadiusMinimum             = 1.0d-3
    self%freefallRadiusMaximum             = 1.0d+2
    self%radiusMinimum                     = 1.0d-3
    self%radiusMaximum                     = 1.0d+2
    self%enclosedDensityRadiusMinimum      = 1.0d-3
    self%enclosedDensityRadiusMaximum      = 1.0d+2
    self%nfwEnclosedDensityTableInitialized=.false.
    self%nfwFreefallTableInitialized       =.false.
    self%nfwInverseTableInitialized        =.false.
    self%nfwTableInitialized               =.false.
    self%lastUniqueID                      =-1
    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                          &
         & call Galacticus_Error_Report                                                                                    &
         &      (                                                                                                          &
         &       'NFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'  //  &
         &       Galacticus_Component_List(                                                                                &
         &                                 'darkMatterProfile'                                                          ,  &
         &                                 defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)   &
         &                                )                                                                             // &
         &      {introspection:location}                                                                                   &
         &      )
    ! Initialize the tabulations.
    call self%tabulate              ()
    call self%inverseAngularMomentum()
    return
  end function nfwConstructorInternal

  subroutine nfwAutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMONFW), intent(inout) :: self

    call calculationResetEvent%attach(self,nfwCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine nfwAutoHook

  subroutine nfwDestructor(self)
    !% Destructor for the {\normalfont \ttfamily nfw} dark matter halo profile class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMONFW), intent(inout) :: self

    if (self%nfwFreefallTableInitialized) then
       call self%nfwFreeFall                      %destroy()
       call self%nfwFreeFallInverse               %destroy()
       deallocate(self%nfwFreefallInverse               )
    end if
    if (self%nfwEnclosedDensityTableInitialized) then
       call self%nfwEnclosedDensity               %destroy()
       call self%nfwEnclosedDensityInverse        %destroy()
       deallocate(self%nfwEnclosedDensityInverse        )
    end if
   if (self%nfwInverseTableInitialized ) then
       call self%nfwSpecificAngularMomentum       %destroy()
       call self%nfwSpecificAngularMomentumInverse%destroy()
       deallocate(self%nfwSpecificAngularMomentumInverse)
    end if
    if (self%nfwTableInitialized        ) then
       call self%nfwConcentrationTable            %destroy()
    end if
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    call calculationResetEvent%detach(self,nfwCalculationReset)
    return
  end subroutine nfwDestructor

  subroutine nfwCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileDMONFW), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    self%specificAngularMomentumScalingsComputed=.false.
    self%maximumVelocityComputed                =.false.
    self%enclosedDensityPrevious                =-1.0d0
    self%densityScalePrevious                   =-1.0d0
    self%enclosedMassPrevious                   =-1.0d0
    self%massScalePrevious                      =-1.0d0
    self%circularVelocityRadiusPrevious         =-1.0d0
    self%lastUniqueID                           =node%uniqueID()
    return
  end subroutine nfwCalculationReset

  subroutine nfwTabulate(self,concentration)
    !% Tabulate properties of the NFW halo profile which must be computed numerically.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout)           :: self
    double precision                         , intent(in   ), optional :: concentration
    integer                                                            :: iConcentration
    logical                                                            :: retabulate
    double precision                                                   :: tableConcentration

    retabulate=.not.self%nfwTableInitialized
    if (present(concentration)) then
       if (concentration < self%concentrationMinimum) then
          self%concentrationMinimum=0.5d0*concentration
          retabulate=.true.
       end if
       if (concentration > self%concentrationMaximum) then
          self%concentrationMaximum=2.0d0*concentration
          retabulate=.true.
       end if
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwTableNumberPoints=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(nfwTablePointsPerDecade))+1
       call self%nfwConcentrationTable%destroy()
       call self%nfwConcentrationTable%create(self%concentrationMinimum,self%concentrationMaximum,self%nfwTableNumberPoints,2)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,self%nfwTableNumberPoints
          tableConcentration=self%nfwConcentrationTable%x(iConcentration)
          call self%nfwConcentrationTable%populate(                   self%profileEnergy           (tableConcentration),iConcentration,table=nfwConcentrationEnergyIndex              )
          call self%nfwConcentrationTable%populate(tableConcentration/self%angularMomentumScaleFree(tableConcentration),iConcentration,table=nfwConcentrationRotationNormalizationIndex)
       end do
       ! Specify that tabulation has been made.
       self%nfwTableInitialized=.true.
    end if
    return
  end subroutine nfwTabulate

  subroutine nfwInverseAngularMomentum(self,specificAngularMomentum)
    !% Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout)           :: self
    double precision                         , intent(in   ), optional :: specificAngularMomentum
    integer                                                            :: iRadius
    logical                                                            :: retabulate

    retabulate=.not.self%nfwInverseTableInitialized
    ! If the table has not yet been made, compute and store the specific angular momenta corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%specificAngularMomentumMinimum=self%specificAngularMomentumScaleFree(self%radiusMinimum)
       self%specificAngularMomentumMaximum=self%specificAngularMomentumScaleFree(self%radiusMaximum)
    end if
    if (present(specificAngularMomentum)) then
       do while (specificAngularMomentum < self%specificAngularMomentumMinimum)
          self%radiusMinimum                 =0.5d0*self%radiusMinimum
          self%specificAngularMomentumMinimum=self%specificAngularMomentumScaleFree(self%radiusMinimum)
          retabulate=.true.
       end do
       do while (specificAngularMomentum > self%specificAngularMomentumMaximum)
          self%radiusMaximum                 =2.0d0*self%radiusMaximum
          self%specificAngularMomentumMaximum=self%specificAngularMomentumScaleFree(self%radiusMaximum)
          retabulate=.true.
       end do
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwInverseTableNumberPoints=int(log10(self%radiusMaximum/self%radiusMinimum)*dble(nfwInverseTablePointsPerDecade))+1
       ! Create a range of radii.
       call self%nfwSpecificAngularMomentum%destroy(                                                                      )
       call self%nfwSpecificAngularMomentum%create (self%radiusMinimum,self%radiusMaximum,self%nfwInverseTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%nfwInverseTableNumberPoints
          call self%nfwSpecificAngularMomentum%populate(                                                                                   &
               &                                        self%specificAngularMomentumScaleFree(self%nfwSpecificAngularMomentum%x(iRadius)), &
               &                                        iRadius                                                                            &
               &                                       )
       end do
       call self%nfwSpecificAngularMomentum%reverse(self%nfwSpecificAngularMomentumInverse)
       ! Specify that tabulation has been made.
       self%nfwInverseTableInitialized=.true.
    end if
    return
  end subroutine nfwInverseAngularMomentum

  double precision function nfwDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    !% in units of Mpc).
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: radiusOverScaleRadius         , scaleRadius, &
         &                                                             virialRadiusOverScaleRadius

    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =darkMatterProfile%scale()
    radiusOverScaleRadius          =radius                       /scaleRadius
    virialRadiusOverScaleRadius    =self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
    nfwDensity=self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)&
         &*basic%mass()/scaleRadius**3
    return
  end function nfwDensity

  double precision function nfwDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: radiusOverScaleRadius, scaleRadius
    !$GLC attributes unused :: self

    basic                 => node%basic            (                 )
    darkMatterProfile     => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius           =  darkMatterProfile%scale()
    radiusOverScaleRadius =  radius/scaleRadius
    nfwDensityLogSlope    = -(1.0d0+3.0d0*radiusOverScaleRadius) &
         &                  /(1.0d0+      radiusOverScaleRadius)
    return
  end function nfwDensityLogSlope

  double precision function nfwRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    !% in units of Mpc).
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMONFW      ), intent(inout)           :: self
    type            (treeNode                     ), intent(inout)           :: node
    double precision                               , intent(in   )           :: moment
    double precision                               , intent(in   ), optional :: radiusMinimum                 , radiusMaximum
    class           (nodeComponentBasic            )              , pointer  :: basic
    class           (nodeComponentDarkMatterProfile)              , pointer  :: darkMatterProfile
    double precision                                                         :: scaleRadius                   , virialRadiusOverScaleRadius, &
         &                                                                      radiusMinimumActual           , radiusMaximumActual

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%virialRadius(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =darkMatterProfile%scale()
    virialRadiusOverScaleRadius    =self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
    nfwRadialMoment                =+basic%mass()                                                &
         &                          *scaleRadius**(moment-2.0d0)                                 &
         &                          /(                                                           &
         &                            +log(1.0d0+virialRadiusOverScaleRadius)                    &
         &                            -          virialRadiusOverScaleRadius                     &
         &                            /   (1.0d0+virialRadiusOverScaleRadius)                    &
         &                          )                                                            &
         &                          /4.0d0                                                       &
         &                          /Pi                                                          &
         &                          *(                                                           &
         &                            +nfwRadialMomentScaleFree(radiusMaximumActual/scaleRadius) &
         &                            -nfwRadialMomentScaleFree(radiusMinimumActual/scaleRadius) &
         &                           )
    return

  contains

    double precision function nfwRadialMomentScaleFree(radius)
      !% Provides the scale-free part of the radial moment of the NFW density profile.
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      use :: Numerical_Comparison    , only : Values_Agree
      implicit none
      double precision, intent(in   ) :: radius

      if (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
         nfwRadialMomentScaleFree=+1.0d0/     (1.0d0+      radius) &
              &                   -2.0d0*atanh(1.0d0+2.0d0*radius)
      else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
         nfwRadialMomentScaleFree=-1.0d0/     (1.0d0      +radius)
      else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
         nfwRadialMomentScaleFree=+1.0d0/     (1.0d0      +radius) &
              &                   +      log  (1.0d0      +radius)
      else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
         nfwRadialMomentScaleFree=+                        radius  &
              &                   -1.0d0/     (1.0d0      +radius) &
              &                   -2.0d0*log  (1.0d0      +radius)
      else
         nfwRadialMomentScaleFree=+(1.0d0+radius)**(moment-1.0d0)                                                     &
              &                   /moment                                                                             &
              &                   /                (moment-1.0d0)                                                     &
              &                   *(                                                                                  &
              &                     - moment                                                                          &
              &                     *  Hypergeometric_2F1([1.0d0-moment,-moment],[2.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                     +(1.0d0+radius)                                                                   &
              &                     *(moment-1.0d0)                                                                   &
              &                     *(                                                                                &
              &                       +(radius/(1.0d0+radius))**moment                                                &
              &                       -Hypergeometric_2F1([     -moment,-moment],[1.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                     )                                                                                 &
              &                    )
      end if
      return
    end function nfwRadialMomentScaleFree

  end function nfwRadialMoment

  double precision function nfwEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: radiusOverScaleRadius      , scaleRadius, &
         &                                                             virialRadiusOverScaleRadius

    basic                       => node%basic            (                 )
    darkMatterProfile           => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                 =  darkMatterProfile%scale()
    radiusOverScaleRadius       =  radius                                      /scaleRadius
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
    nfwEnclosedMass             =  self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius) &
         &                         *basic%mass()
    return
  end function nfwEnclosedMass

  double precision function nfwPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentDarkMatterProfile)               , pointer  :: darkMatterProfile
    double precision                                , parameter               :: radiusSmall                =1.0d-10
    double precision                                                          :: radiusOverScaleRadius              , radiusTerm, &
         &                                                                       virialRadiusOverScaleRadius

    if (present(status)) status=structureErrorCodeSuccess
    darkMatterProfile           => node%darkMatterProfile(autoCreate=.true.)
    radiusOverScaleRadius       =  radius                                      /darkMatterProfile%scale()
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%virialRadius(node)/darkMatterProfile%scale()
    if (radiusOverScaleRadius < radiusSmall) then
       ! Use a series solution for very small radii.
       radiusTerm=1.0d0-0.5d0*radiusOverScaleRadius
    else
       ! Use the full expression for larger radii.
       radiusTerm=log(1.0d0+radiusOverScaleRadius)/radiusOverScaleRadius
    end if
    nfwPotential=-virialRadiusOverScaleRadius              &
         &       *radiusTerm                               &
         &       /(                                        &
         &         +log(1.0d0+virialRadiusOverScaleRadius) &
         &         -          virialRadiusOverScaleRadius  &
         &         /   (1.0d0+virialRadiusOverScaleRadius) &
         &        )                                        &
         &       *self%darkMatterHaloScale_%virialVelocity(node)**2
    return
  end function nfwPotential

  double precision function nfwCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: radius

    if (radius > 0.0d0) then
       ! Check if node differs from previous one for which we performed calculations.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       ! Compute the circular velocity if the radius has changed.
       if (radius /= self%circularVelocityRadiusPrevious) then
          self%circularVelocityPrevious      =sqrt(gravitationalConstantGalacticus*self%enclosedMass(node,radius)/radius)
          self%circularVelocityRadiusPrevious=radius
       end if
       nfwCircularVelocity=self%circularVelocityPrevious
    else
       nfwCircularVelocity=0.0d0
    end if
    return
  end function nfwCircularVelocity

  double precision function nfwCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes            , only : nodeComponentBasic             , nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (nodeComponentBasic            ), pointer       :: basic
    ! The circular velocity (in scale-free units) at the peak of the NFW rotation curve. Numerical value found using Mathematica.
    double precision                                , parameter     :: circularVelocityMaximumScaleFree=0.4649909628174221d0
    double precision                                                :: scaleRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if maximum velocity is already computed. Compute and store if not.
    if (.not.self%maximumVelocityComputed) then
       basic             => node             %basic            (                 )
       darkMatterProfile => node             %darkMatterProfile(autoCreate=.true.)
       scaleRadius       =  darkMatterProfile%scale            (                 )
       ! Ensure mass profile normalization factor has been computed.
       call nfwMassNormalizationFactor(self,self%darkMatterHaloScale_%virialRadius(node)/scaleRadius)
       ! Evaluate the circular velocity at the peak of the rotation curve.
       self%maximumVelocityPrevious=+circularVelocityMaximumScaleFree                                       &
            &                       *sqrt(                                                                  &
            &                             +gravitationalConstantGalacticus                                  &
            &                             *basic                          %mass                          () &
            &                             *self                           %nfwNormalizationFactorPrevious   &
            &                             /scaleRadius                                                      &
            &                            )
       self%maximumVelocityComputed= .true.
    end if
    nfwCircularVelocityMaximum=self%maximumVelocityPrevious
    return
  end function nfwCircularVelocityMaximum

  double precision function nfwRadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout)          :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: radiusOverScaleRadius      , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    darkMatterProfile           => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                 =  darkMatterProfile%scale()
    radiusOverScaleRadius       =  radius                                      /scaleRadius
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
    nfwRadialVelocityDispersion =  +self%radialVelocityDispersionScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius) &
         &                         *self%darkMatterHaloScale_%virialVelocity(node)
    return
  end function nfwRadialVelocityDispersion

  double precision function nfwRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: specificAngularMomentumScaleFree

    ! Return immediately with zero radius for non-positive specific angular momenta.
    if (specificAngularMomentum <= 0.0d0) then
       nfwRadiusFromSpecificAngularMomentum=0.0d0
       return
    end if
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if scalings are already computed. Compute and store if not.
    if (.not.self%specificAngularMomentumScalingsComputed) then
       ! Flag that scale quantities are now computed.
       self%specificAngularMomentumScalingsComputed=.true.

       ! Get the dark matter profile.
       darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

       ! Get the scale radius.
       self%specificAngularMomentumLengthScale=darkMatterProfile%scale()

       ! Get the specific angular momentum scale.
       self%specificAngularMomentumScale=self%specificAngularMomentumLengthScale*self%circularVelocity(node&
            &,self%specificAngularMomentumLengthScale)
    end if

    ! Compute the specific angular momentum in scale free units (using the scale length for distances the sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/self%specificAngularMomentumScale
    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%inverseAngularMomentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    nfwRadiusFromSpecificAngularMomentum=self%nfwSpecificAngularMomentumInverse%interpolate(specificAngularMomentumScaleFree)

    ! Convert to a physical radius.
    nfwRadiusFromSpecificAngularMomentum=nfwRadiusFromSpecificAngularMomentum*self%specificAngularMomentumLengthScale
    return
  end function nfwRadiusFromSpecificAngularMomentum

  double precision function nfwRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: concentration

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/darkMatterProfile%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the rotation normalization by interpolation.
    nfwRotationNormalization=self%nfwConcentrationTable%interpolate(concentration,table&
         &=nfwConcentrationRotationNormalizationIndex)/self%darkMatterHaloScale_%virialRadius(node)
    return
  end function nfwRotationNormalization

  double precision function nfwEnergy(self,node)
    !% Return the energy of an NFW halo density profile.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (nodeComponentBasic            ), pointer       :: basic
    double precision                                                :: concentration

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/darkMatterProfile%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy by interpolation.
    nfwEnergy=+self %nfwConcentrationTable%interpolate   (concentration,table=nfwConcentrationEnergyIndex)    &
         &    *self %darkMatterHaloScale_ %virialVelocity(node                                           )**2 &
         &    *basic                      %mass          (                                               )
    return
  end function nfwEnergy

  double precision function nfwEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of an NFW halo density profile.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), target  :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    class           (nodeComponentBasic            )               , pointer :: basic
    double precision                                                         :: concentration    , energy, &
         &                                                                      energyGradient

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/darkMatterProfile%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy gradient by interpolation.
    energy             =+self%nfwConcentrationTable%interpolate        (concentration,table=nfwConcentrationEnergyIndex)
    energyGradient     =+self%nfwConcentrationTable%interpolateGradient(concentration,table=nfwConcentrationEnergyIndex)

    nfwEnergyGrowthRate=+    self                     %energy                   (node)&
         &              *(                                                            &
         &                +       basic               %accretionRate           (    ) &
         &                /       basic               %mass                    (    ) &
         &                +2.0d0                                                      &
         &                *  self%darkMatterHaloScale_%virialVelocityGrowthRate(node) &
         &                /  self%darkMatterHaloScale_%virialVelocity          (node) &
         &                +(                                                          &
         &                  +energyGradient                                           &
         &                  *concentration                                            &
         &                  /energy                                                   &
         &                 )                                                          &
         &                *(                                                          &
         &                  +self%darkMatterHaloScale_%virialRadiusGrowthRate  (node) &
         &                  /self%darkMatterHaloScale_%virialRadius            (node) &
         &                  -     darkMatterProfile   %scaleGrowthRate         (    ) &
         &                  /     darkMatterProfile   %scale                   (    ) &
         &                 )                                                          &
         &               )
    return
  end function nfwEnergyGrowthRate

  double precision function nfwAngularMomentumScaleFree(self,concentration)
    !% Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    !% in an NFW dark matter profile with given {\normalfont \ttfamily concentration}. This is given by:
    !% \begin{equation}
    !% J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    !% \end{equation}
    !% where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    !% \begin{equation}
    !% J = \left. \left[ 1 + c - 2 \ln (1+c) - {1 \over 1+c} \right] \right/ \left[ \ln(1+c)-{c\over 1+c} \right].
    !% \end{equation}
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: concentration
    !$GLC attributes unused :: self

    nfwAngularMomentumScaleFree=(1.0d0+concentration-2.0d0*log(1.0d0+concentration)-1.0d0/(1.0d0+concentration)) &
         &/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function nfwAngularMomentumScaleFree

  double precision function nfwSpecificAngularMomentumScaleFree(self,radius)
    !% Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    !% {\normalfont \ttfamily radius} (in units of the scale radius) in an NFW profile.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: radius
    double precision                         , parameter     :: radiusSmall=1.0d-9
    !$GLC attributes unused :: self

    if (radius < radiusSmall) then
       ! Use a series expenasion solution for accuracy.
       nfwSpecificAngularMomentumScaleFree=+radius**1.5d0                &
            &                              /sqrt(                        &
            &                                    +    2.0d0              &
            &                                    *log(2.0d0)             &
            &                                    -    1.0d0              &
            &                                   )                        &
            &                              *(                            &
            &                                        +   1.0d0           &
            &                                -      radius               &
            &                                *(                          &
            &                                        -   2.0d0/    3.0d0 &
            &                                  +    radius               &
            &                                  *(                        &
            &                                        +  19.0d0/   36.0d0 &
            &                                    +  radius               &
            &                                    *(                      &
            &                                        - 121.0d0/  270.0d0 &
            &                                      +radius               &
            &                                      *(                    &
            &                                        +5123.0d0/12960.0d0 &
            &                                       )                    &
            &                                     )                      &
            &                                   )                        &
            &                                 )                          &
            &                               )
    else
       ! Use the full solution.
       nfwSpecificAngularMomentumScaleFree=sqrt(radius*self%enclosedMassScaleFree(radius,1.0d0))
    end if
    return
  end function nfwSpecificAngularMomentumScaleFree

  double precision function nfwEnclosedMassScaleFree(self,radius,concentration)
    !% Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\normalfont \ttfamily concentration} at the
    !% given {\normalfont \ttfamily radius} (given in units of the scale radius).
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: concentration                                                   , radius
    double precision                         , parameter     :: minimumRadiusForExactSolution          =1.0d-6
    ! Precomputed NFW normalization factor for unit radius.
    double precision                         , parameter     :: nfwNormalizationFactorUnitRadius       =log(2.0d0)-0.5d0

    if (radius == 1.0d0) then
       nfwEnclosedMassScaleFree=nfwNormalizationFactorUnitRadius
    else if (radius >= minimumRadiusForExactSolution) then
       nfwEnclosedMassScaleFree=(log(1.0d0+radius)-radius/(1.0d0+radius))
    else
       nfwEnclosedMassScaleFree=(radius**2)*(0.5d0+radius*(-2.0d0/3.0d0+radius*(0.75d0+radius*(-0.8d0))))
    end if
    ! Compute the mass profile normalization factor.
    call nfwMassNormalizationFactor(self,concentration)
    ! Evaluate the scale-free enclosed mass.
    nfwEnclosedMassScaleFree=nfwEnclosedMassScaleFree*self%nfwNormalizationFactorPrevious
    return
  end function nfwEnclosedMassScaleFree

  subroutine nfwMassNormalizationFactor(self,concentration)
    !% Compute the normalization factor for the NFW mass profile.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: concentration
    ! Precomputed NFW normalization factor for unit concentration.
    double precision                         , parameter     :: nfwNormalizationFactorUnitConcentration=1.0d0/(log(2.0d0)-0.5d0)

    ! Check if we were called with a different concentration compared to the previous call.
    if (concentration /= self%concentrationPrevious) then
       ! We were, so recompute the normalization factor.
       if (concentration == 1.0d0) then
          self%nfwNormalizationFactorPrevious=nfwNormalizationFactorUnitConcentration
       else
          self%nfwNormalizationFactorPrevious=1.0d0/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))
       end if
       self%concentrationPrevious=concentration
    end if
    return
  end subroutine nfwMassNormalizationFactor

  double precision function nfwRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in units of the scale radius) in an NFW dark matter profile with given {\normalfont \ttfamily
    !% concentration} which encloses a given density (in units of the virial mass per cubic scale radius).
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: density
    class           (nodeComponentBasic            ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile
    double precision                                                        :: scaleRadius                , densityScaleFree, &
         &                                                                     virialRadiusOverScaleRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Get scale radius if required.
    if (self%densityScalePrevious < 0.0d0 .or. density /= self%enclosedDensityPrevious) then
       darkMatterProfile => node             %darkMatterProfile(autoCreate=.true.)
       scaleRadius       =  darkMatterProfile%scale            (                 )
       ! Compute the density scale if necessary.
       if (self%densityScalePrevious < 0.0d0) then
          ! Extract profile parameters.
          basic                       => node                     %basic       (    )
          virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
          ! Compute density normalization scale.
          self%densityScalePrevious=+scaleRadius                                                         **3 &
               &                    /basic      %mass                 (                                 )    &
               &                    /self       %enclosedMassScaleFree(1.0d0,virialRadiusOverScaleRadius)
       end if
       ! Compute radius enclosing density if necessary.
       if (density /= self%enclosedDensityPrevious) then
          self%enclosedDensityPrevious=density
          ! Compute scaled density.
          densityScaleFree=+     density              &
               &           *self%densityScalePrevious
          ! Ensure density table spans required range.
          call self%enclosedDensityTabulate(densityScaleFree)
          ! Interpolate in density table to find the required radius.
          self%enclosingDensityRadiusPrevious=self%nfwEnclosedDensityInverse%interpolate(-densityScaleFree)*scaleRadius
       end if
    end if
    nfwRadiusEnclosingDensity=self%enclosingDensityRadiusPrevious
    return
  end function nfwRadiusEnclosingDensity

  double precision function nfwRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in an NFW dark matter profile with given {\normalfont \ttfamily
    !% concentration} which encloses a given mass (in $M_\odot$).
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    use :: Lambert_Ws      , only : Lambert_W0
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: mass
    class           (nodeComponentBasic            ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile
    double precision                                                        :: scaleRadius                , massScaleFree, &
         &                                                                     virialRadiusOverScaleRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Get scale radius if required.
    if (self%massScalePrevious < 0.0d0 .or. mass /= self%enclosedMassPrevious) then
       darkMatterProfile => node             %darkMatterProfile(autoCreate=.true.)
       scaleRadius       =  darkMatterProfile%scale            (                 )
       ! Compute the mass scale if necessary.
       if (self%massScalePrevious < 0.0d0) then
          ! Extract profile parameters.
          basic                       => node%basic()
          virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%virialRadius(node)/scaleRadius
          ! Compute the mass profile normalization factor.
          call nfwMassNormalizationFactor(self,virialRadiusOverScaleRadius)
          ! Compute mass normalization scale.
          self%massScalePrevious      =  1.0d0/basic%mass()/self%nfwNormalizationFactorPrevious
       end if
       ! Compute radius enclosing mass if necessary.
       if (mass /= self%enclosedMassPrevious) then
          self%enclosedMassPrevious       = mass
          ! Compute scaled mass.
          massScaleFree                   =+mass                                           &
               &                           *self%massScalePrevious
          ! Compute radius.
          self%enclosingMassRadiusPrevious=-(                                              &
               &                             +1.0d0/Lambert_W0(-exp(-1.0d0-massScaleFree)) &
               &                             +1.0d0                                        &
               &                            )                                              &
               &                           *scaleRadius
       end if
    end if
    nfwRadiusEnclosingMass=self%enclosingMassRadiusPrevious
    return
  end function nfwRadiusEnclosingMass

  double precision function nfwDensityEnclosedByRadiusScaleFree(self,radius)
    !% Returns the density (in units of the virial mass per cubic scale radius) in an NFW dark matter profile with given {\normalfont \ttfamily
    !% concentration} which is enclosed a given radius (in units of the scale radius).
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: radius

    nfwDensityEnclosedByRadiusScaleFree=+3.0d0                                       &
         &                              *self%enclosedMassScaleFree(radius,1.0d0)    &
         &                              /4.0d0                                       &
         &                              /Pi                                          &
         &                              /                           radius       **3
    return
  end function nfwDensityEnclosedByRadiusScaleFree

  subroutine nfwEnclosedDensityTabulate(self,enclosedDensityScaleFree)
    !% Tabulates the enclosed density vs. radius for NFW halos.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: enclosedDensityScaleFree
    logical                                                  :: retabulate
    integer                                                  :: iRadius

    retabulate=.not.self%nfwEnclosedDensityTableInitialized
    ! If the table has not yet been made, compute and store the enclosed density corresponding to the minimum and maximum radii
    ! that will be tabulated by default.
    if (retabulate) then
       self%enclosedDensityMinimum=self%densityEnclosedByRadiusScaleFree(self%enclosedDensityRadiusMaximum)
       self%enclosedDensityMaximum=self%densityEnclosedByRadiusScaleFree(self%enclosedDensityRadiusMinimum)
    end if
    do while (enclosedDensityScaleFree < self%enclosedDensityMinimum)
       self%enclosedDensityRadiusMaximum=2.0d0*self%enclosedDensityRadiusMaximum
       self%enclosedDensityMinimum=self%densityEnclosedByRadiusScaleFree(self%enclosedDensityRadiusMaximum)
       retabulate=.true.
    end do
    do while (enclosedDensityScaleFree > self%enclosedDensityMaximum)
       self%enclosedDensityRadiusMinimum=0.5d0*self%enclosedDensityRadiusMinimum
       self%enclosedDensityMaximum=self%densityEnclosedByRadiusScaleFree(self%enclosedDensityRadiusMinimum)
       retabulate=.true.
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwEnclosedDensityTableNumberPoints=int(log10(self%enclosedDensityRadiusMaximum/self%enclosedDensityRadiusMinimum)*dble(nfwEnclosedDensityTablePointsPerDecade))+1
       ! Create the table.
       call self%nfwEnclosedDensity%destroy(                                                                                                            )
       call self%nfwEnclosedDensity%create (self%enclosedDensityRadiusMinimum,self%enclosedDensityRadiusMaximum,self%nfwEnclosedDensityTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%nfwEnclosedDensityTableNumberPoints
          call self%nfwEnclosedDensity%populate(                                                                            &
               &                                -self%densityEnclosedByRadiusScaleFree(self%nfwEnclosedDensity%x(iRadius)), &
               &                                                                                                 iRadius    &
               &                               )
       end do
       call self%nfwEnclosedDensity%reverse(self%nfwEnclosedDensityInverse)
       ! Specify that tabulation has been made.
       self%nfwEnclosedDensityTableInitialized=.true.
    end if
    return
  end subroutine nfwEnclosedDensityTabulate

  double precision function nfwDensityScaleFree(self,radius,concentration)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with
    !% given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: concentration, radius
    !$GLC attributes unused :: self

    nfwDensityScaleFree=1.0d0/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))/radius/(1.0d0+radius)**2/4.0d0/Pi
    return
  end function nfwDensityScaleFree

  double precision function nfwRadialVelocityDispersionScaleFree(self,radius,concentration)
    !% Returns the radial velocity dispersion (in units of the virial velocity) in an NFW dark matter profile with given
    !% {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)
    !% using the result derived by \citeauthor{lokas_properties_2001}~(\citeyear{lokas_properties_2001}; eqn.~14). Note that
    !% approximate solutions are used at small and large radii.
    use :: Dilogarithms            , only : Dilogarithm
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout)            :: self
    double precision                         , intent(in   )            :: concentration, radius
    double precision                                        , parameter :: minimumRadiusForExactSolution   =1.0d-6
    double precision                                        , parameter :: maximumRadiusForExactSolution   =5.0d2
    ! Precomputed NFW normalization factor for unit radius.
    double precision                                        , parameter :: nfwNormalizationFactorUnitRadius=-8.5d0+Pi**2-6.0d0*log(2.0d0)+6.0d0*log(2.0d0)**2
    double precision                                                    :: radialVelocityDispersionSquare

    if (radius == 1.0d0) then
       radialVelocityDispersionSquare=nfwNormalizationFactorUnitRadius
    else if (radius >= maximumRadiusForExactSolution) then
       radialVelocityDispersionSquare=+(2.0d0+radius)                                      &
            &                         *(                                                   &
            &                           +         188.0d0-75.0d0*radius                    &
            &                           +20.0d0*(-  8.0d0+ 5.0d0*radius)*log(1.0d0+radius) &
            &                          )                                                   &
            &                         /(400.0d0*radius**3)
    else if (radius >= minimumRadiusForExactSolution) then
       radialVelocityDispersionSquare=+0.5d0*radius*(1.0d0+radius)**2 &
            &                         *(                              &
            &                           +Pi**2                        &
            &                           -log(radius)                  &
            &                           -1.0d0/       radius          &
            &                           -1.0d0/(1.0d0+radius)**2      &
            &                           -6.0d0/(1.0d0+radius)         &
            &                           +(                            &
            &                             +1.0d0+ 1.0d0/radius**2     &
            &                                   - 4.0d0/radius        &
            &                             -2.0d0/(1.0d0+radius)       &
            &                            )                            &
            &                           *log(1.0d0+radius)            &
            &                           +3.0d0*log(1.0d0+radius)**2   &
            &                           +6.0d0*Dilogarithm(-radius)   &
            &                          )
    else if (radius > 0.0d0) then
       radialVelocityDispersionSquare=+0.25d0      *(-23.0d0       + 2.0d0*Pi**2- 2.0d0*log(radius))*radius    &
            &                         +             (-59.0d0/6.0d0 +       Pi**2-       log(radius))*radius**2 &
            &                         +1.0d0/24.0d0*(-101.0d0      +12.0d0*Pi**2-12.0d0*log(radius))*radius**3
    else
       radialVelocityDispersionSquare=0.0d0
    end if
    nfwRadialVelocityDispersionScaleFree=sqrt(radialVelocityDispersionSquare)
    ! Compute the normalization factor.
    call nfwMassNormalizationFactor(self,concentration)
    ! Evaluate the scale-free radial velocity dispersion.
    nfwRadialVelocityDispersionScaleFree=+nfwRadialVelocityDispersionScaleFree      &
         &                               *sqrt(                                     &
         &                                     +self%nfwNormalizationFactorPrevious &
         &                                     *concentration                       &
         &                                    )
    return
  end function nfwRadialVelocityDispersionScaleFree

  double precision function nfwProfileEnergy(self,concentration)
    !% Computes the total energy of an NFW profile halo of given {\normalfont \ttfamily concentration} using the methods of
    !% \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A), except for potential energy
    !% which is computed using the result derived by \citeauthor{mo_formation_1998}~(\citeyear{mo_formation_1998}; eqn.~23).
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: concentration
    type            (integrator             )                :: integratorJeans       , integratorKinetic
    double precision                                         :: jeansEquationIntegral , kineticEnergy    , &
         &                                                      kineticEnergyIntegral , potentialEnergy  , &
         &                                                      radiusMinimum         , radiusMaximum    , &
         &                                                      concentrationParameter

    ! Compute the potential energy.
    potentialEnergy=-0.5d0                               &
         &          *(                                   &
         &            +1.0d0                             &
         &            -1.0d0                             &
         &            /         (1.0d0+concentration)**2 &
         &            -2.0d0*log(1.0d0+concentration)    &
         &            /         (1.0d0+concentration)    &
         &           )                                   &
         &          /(                                   &
         &            +                concentration     &
         &            /         (1.0d0+concentration)    &
         &                  -log(1.0d0+concentration)    &
         &           )                               **2
    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum         =        concentration
    radiusMaximum         =100.0d0*concentration
    concentrationParameter=        concentration
    integratorJeans       =integrator               (nfwJeansEquationIntegrand,toleranceRelative=1.0d-3)
    jeansEquationIntegral =integratorJeans%integrate(radiusMinimum            ,radiusMaximum           )
    ! Compute the kinetic energy.
    radiusMinimum         =0.0d0
    radiusMaximum         =concentration
    concentrationParameter=concentration
    integratorKinetic     =integrator                 (nfwKineticEnergyIntegrand,toleranceRelative=1.0d-3)
    kineticEnergyIntegral =integratorKinetic%integrate(radiusMinimum            ,radiusMaximum           )
    kineticEnergy         =2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    ! Compute the total energy.
    nfwProfileEnergy=(potentialEnergy+kineticEnergy)*concentration
    return

  contains

    double precision function nfwKineticEnergyIntegrand(radius)
      !% Integrand for NFW profile kinetic energy.
      implicit none
      double precision, intent(in   ) :: radius

      nfwKineticEnergyIntegrand=self%EnclosedMassScaleFree(radius,concentrationParameter) &
           &                    *self%densityScaleFree    (radius,concentrationParameter) &
           &                    *radius
      return
    end function nfwKineticEnergyIntegrand

    double precision function nfwJeansEquationIntegrand(radius)
      !% Integrand for NFW profile Jeans equation.
      implicit none
      double precision, intent(in   ) :: radius

      nfwJeansEquationIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter) &
           &                    *self%densityScaleFree    (radius,concentrationParameter) &
           &                    /radius**2
      return
    end function nfwJeansEquationIntegrand

  end function nfwProfileEnergy

  double precision function nfwKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the NFW density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use :: Exponential_Integrals, only : Cosine_Integral               , Sine_Integral
    use :: Galacticus_Nodes     , only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), target  :: node
    double precision                                , intent(in   )          :: waveNumber
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: concentration      , radiusScale, &
         &                                                                      waveNumberScaleFree

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=darkMatterProfile%scale()

    ! Compute the concentration parameter.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/radiusScale

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    nfwKSpace= (                                                                                                                                          &
         &      +sin(              waveNumberScaleFree)*(Sine_Integral  ((1.0d0+concentration)*waveNumberScaleFree)-Sine_Integral  (waveNumberScaleFree)) &
         &      -sin(concentration*waveNumberScaleFree)/(1.0d0+concentration)/waveNumberScaleFree                                                         &
         &      +cos(              waveNumberScaleFree)*(Cosine_Integral((1.0d0+concentration)*waveNumberScaleFree)-Cosine_Integral(waveNumberScaleFree)) &
         &     )                                                                                                                                          &
         &    /(log(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function nfwKSpace

  double precision function nfwFreefallRadius(self,node,time)
    !% Returns the freefall radius in the NFW density profile at the specified {\normalfont \ttfamily time} (given in Gyr).
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: time
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: concentration    , freefallTimeScaleFree, &
         &                                                             radiusScale      , timeScale            , &
         &                                                             velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       nfwFreefallRadius=0.0d0
       return
    end if

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=darkMatterProfile%scale()

    ! Get the concentration.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/radiusScale

    ! Get the virial velocity.
    velocityScale=self%darkMatterHaloScale_%virialVelocity(node)

    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                                                            &
         &    *radiusScale                                                                        &
         &    /velocityScale                                                                      &
         &    /sqrt(concentration/(log(1.0d0+concentration)-concentration/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius.
    nfwFreefallRadius=self%nfwFreefallInverse%interpolate(freefallTimeScaleFree)*radiusScale
    return
  end function nfwFreefallRadius

  double precision function nfwFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the NFW density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMONFW       ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: time
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: concentration    , freefallTimeScaleFree, &
         &                                                             radiusScale      , timeScale            , &
         &                                                             velocityScale

    ! For non-positive freefall times, return the limiting value for small radii.
    if (time <= 0.0d0) then
       nfwFreefallRadiusIncreaseRate=0.0d0
       return
    end if

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=darkMatterProfile%scale()

    ! Get the concentration.
    concentration=self%darkMatterHaloScale_%virialRadius(node)/radiusScale

    ! Get the virial velocity.
    velocityScale=self%darkMatterHaloScale_%virialVelocity(node)

    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                                                            &
         &    *radiusScale                                                                        &
         &    /velocityScale                                                                      &
         &    /sqrt(concentration/(log(1.0d0+concentration)-concentration/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius growth rate.
    nfwFreefallRadiusIncreaseRate=self%nfwFreefallInverse%interpolateGradient(freefallTimeScaleFree)*radiusScale/timeScale
    return
  end function nfwFreefallRadiusIncreaseRate

  subroutine nfwFreefallTabulate(self,freefallTimeScaleFree)
    !% Tabulates the freefall time vs. freefall radius for NFW halos.
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: freefallTimeScaleFree
    logical                                                  :: retabulate
    integer                                                  :: iRadius

    retabulate=.not.self%nfwFreefallTableInitialized
    ! If the table has not yet been made, compute and store the freefall corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusMinimum)
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
    end if
    do while (freefallTimeScaleFree < self%freefallTimeMinimum)
       self%freefallRadiusMinimum=0.5d0*self%freefallRadiusMinimum
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusMinimum)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > self%freefallTimeMaximum)
       self%freefallRadiusMaximum=2.0d0*self%freefallRadiusMaximum
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
       retabulate=.true.
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwFreefallTableNumberPoints=int(log10(self%freefallRadiusMaximum/self%freefallRadiusMinimum)*dble(nfwFreefallTablePointsPerDecade))+1
       ! Create the table.
       call self%nfwFreefall%destroy(                                                                        )
       call self%nfwFreefall%create (self%freefallRadiusMinimum,self%freefallRadiusMaximum,self%nfwFreefallTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%nfwFreefallTableNumberPoints
          call self%nfwFreefall%populate(                                                         &
               &                         self%freefallTimeScaleFree(self%nfwFreefall%x(iRadius)), &
               &                                                                       iRadius    &
               &                        )
       end do
       call self%nfwFreefall%reverse(self%nfwFreefallInverse)
       ! Specify that tabulation has been made.
       self%nfwFreefallTableInitialized=.true.
    end if
    return
  end subroutine nfwFreefallTabulate

  double precision function nfwFreefallTimeScaleFree(self,radius)
    !% Compute the freefall time in a scale-free NFW halo.
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMONFW), intent(inout) :: self
    double precision                         , intent(in   ) :: radius
    double precision                         , parameter     :: radiusSmall=4.0d-6
    type            (integrator             )                :: integrator_
    double precision                                         :: radiusEnd         , radiusStart
    !$GLC attributes unused :: self

    if (radius > radiusSmall) then
       ! Use the full solution.
       radiusStart             =radius
       radiusEnd               =0.0d0
       integrator_             =integrator           (nfwFreefallTimeScaleFreeIntegrand,toleranceRelative=1.0d-3)
       nfwFreefallTimeScaleFree=integrator_%integrate(radiusEnd                        ,radiusStart             )
    else
       ! Use an approximation here, found by taking series expansions of the logarithms in the integrand and keeping only the
       ! first order terms.
       nfwFreefallTimeScaleFree=2.0d0*sqrt(radius)
    end if
    return

  contains

    double precision function nfwFreefallTimeScaleFreeIntegrand(radius)
      !% Integrand function used for finding the free-fall time in NFW halos.
      implicit none
      double precision, intent(in   ) :: radius
      double precision, parameter     :: radiusSmall        =1.0d-6
      double precision, parameter     :: radiusSmallFraction=1.0d-3
      double precision                :: x

      if (radius < radiusSmall) then
         ! Use a series approximation for small radii.
         nfwFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-1.0d0+radius*(0.5d0-radius/3.0d0)
      else if (radius > radiusStart*(1.0d0-radiusSmallFraction)) then
         ! Use a series approximation for radii close to the initial radius.
         x=1.0d0-radius/radiusStart
         nfwFreefallTimeScaleFreeIntegrand=+(1.0d0/(1.0d0+radiusStart)-log(1.0d0+radiusStart)/radiusStart)*x &
              &                            +(                                                                &
              &                              +0.5d0*radiusStart/(1.0d0+radiusStart)**2                       &
              &                              +(radiusStart-(1.0d0+radiusStart)*log(1.0d0+radiusStart))       &
              &                              /radiusStart                                                    &
              &                              /(1.0d0+radiusStart)                                            &
              &                             )*x**2
      else
         ! Use full expression for larger radii.
         nfwFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-log(1.0d0+radius)/radius
      end if
      nfwFreefallTimeScaleFreeIntegrand=1.0d0/sqrt(-2.0d0*nfwFreefallTimeScaleFreeIntegrand)
      return
    end function nfwFreefallTimeScaleFreeIntegrand

  end function nfwFreefallTimeScaleFree
