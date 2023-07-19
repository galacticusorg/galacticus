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
  An implementation of \cite{burkert_structure_1995} dark matter halo profiles.
  !!}

  use :: Kind_Numbers            , only : kind_int8
  use :: Numerical_Constants_Math, only : Pi
  use :: Tables                  , only : table1D  , table1DLogarithmicLinear

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOBurkert">
   <description>
    A dark matter profile DMO class which implements the \citep{burkert_structure_1995} density profile is used
    \begin{equation}
      \rho_\mathrm{dark matter}(r) \propto \left(1+{r\over r_\mathrm{s}}\right)^{-1} \left(1+[{r\over
      r_\mathrm{s}}]^2\right)^{-1},
    \end{equation}
    normalized such that the total mass of the \gls{node} is enclosed with the virial radius and with the scale length
    $r_\mathrm{s} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see
    \refPhysics{darkMatterProfileConcentration}). The mass enclosed within radius $r$ is given by
    \begin{equation}
    M(&lt;r) = M_\mathrm{virial} {2 \log(1 + R) + \log(1 + R^2) -2 \tan^{-1}(R) \over 2 \log(1 + c) + \log(1 + c^2) -2 \tan^{-1}(c)},
    \end{equation}
    where $R=r/r_\mathrm{s}$. The associated gravitational potential is
    \begin{equation}
    \Phi(r) = -\mathrm{G} \left(1+{1 \over R}\right) { 2 \tan^{-1}(R) - 2 \log(1 + R) + \log(1 + R^2) \over -2 \tan^{-1}(c) + 2
    \log(1 + c) + \log(1 + c^2) }.
    \end{equation}
    The peak of the rotation curve occurs at $R=3.2446257246042642$ (found by numerical solution), and the Fourier transform of
    the profile, $F(k) = \int_0^c 4 \pi r^2 \exp(-i k r) \rho(r) \mathrm{d} r / k r$ (needed in calculations of clustering
    using the halo model) is given by
    \begin{eqnarray}
      F(k) &amp;=&amp; \left\{2 \exp(-i k) \mathrm{C}_\mathrm{i}(k) - 2 \exp(-i k) \mathrm{C}_\mathrm{i}(k[1 + c]) + (1 + i)
      \left[-i \exp(k) \pi - \exp(k) \mathrm{E}_\mathrm{i}(-k) \right. \right. \nonumber \\
        &amp; &amp; +i \exp(-k) \mathrm{E}_\mathrm{i}(k) + \exp(k) \mathrm{E}_\mathrm{i}(i [i + c] k) - i \exp(-k)
        \mathrm{E}_\mathrm{i}(k [1 + i c ]) + (1 + i) \exp(-i k) \mathrm{S}_\mathrm{i}(k) \nonumber \\
        &amp; &amp; \left. \left. - (1 + i) \exp(-i k) \mathrm{S}_\mathrm{i}(k[1 + c])\right]\right\}/\left[k \left\{-2
        \tan^{-1}(c) + 2 \log(1 + c) + \log(1 + c^2)\right\}\right].
    \end{eqnarray}
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOBurkert
     !!{
     A dark matter halo profile class implementing \cite{burkert_structure_1995} dark matter halos.
     !!}
     private
     ! Minimum and maximum concentrations to tabulate.
     double precision                                        :: concentrationMinimum                                    , concentrationMaximum
     ! Minimum and maximum radii to tabulate.
     double precision                                        :: freefallRadiusMinimum                                   , radiusMinimum                                , &
          &                                                     densityRadiusMinimum                                    , radialVelocityDispersionRadiusMinimum
     double precision                                        :: freefallRadiusMaximum                                   , radiusMaximum                                , &
          &                                                     densityRadiusMaximum                                    , radialVelocityDispersionRadiusMaximum
     double precision                                        :: freefallTimeMinimum                                     , specificAngularMomentumMinimum               , &
          &                                                     densityMinimum
     double precision                                        :: freefallTimeMaximum                                     , specificAngularMomentumMaximum               , &
          &                                                     densityMaximum
     ! Tables of Burkert properties.
     logical                                                 :: burkertFreefallTableInitialized                 =.false., burkertInverseTableInitialized       =.false., &
          &                                                     burkertTableInitialized                         =.false., burkertDensityTableInitialized       =.false., &
          &                                                     burkertRadialVelocityDispersionTableInitialized =.false.
     integer                                                 :: burkertFreefallTableNumberPoints                        , burkertInverseTableNumberPoints              , &
          &                                                     burkertTableNumberPoints                                , burkertDensityTableNumberPoints              , &
          &                                                     burkertRadialVelocityDispersionTableNumberPoints
     type            (table1DLogarithmicLinear)              :: burkertConcentrationTable
     ! Tables.
     type            (table1DLogarithmicLinear)              :: burkertFreeFall                                         , burkertSpecificAngularMomentum               , &
          &                                                     burkertDensityTable                                     , burkertRadialVelocityDispersionTable
     class           (table1D                 ), allocatable :: burkertFreefallInverse                                  , burkertSpecificAngularMomentumInverse        , &
          &                                                     burkertDensityTableInverse
     ! Module variables used in integrations.
     double precision                                        :: concentrationParameter                                  , radiusStart
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )              :: lastUniqueID
     ! Record of whether or not quantities have been computed.
     logical                                                 :: specificAngularMomentumScalingsComputed                 , maximumVelocityComputed
     ! Stored values of computed quantities.
     double precision                                        :: specificAngularMomentumLengthScale                      , specificAngularMomentumScale                 , &
          &                                                     concentrationPrevious                                   , burkertNormalizationFactorPrevious           , &
          &                                                     maximumVelocityPrevious
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
       <method description="Returns the density (in units such that the virial mass and scale length are unity) in a Burkert dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="densityScaleFree" />
       <method description="Returns the enclosed mass (in units of the virial mass) in a Burkert dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius)." method="enclosedMassScaleFree" />
       <method description="Returns the scale-free radial velocity dispersion in a Burkert dark matter profile." method="radialVelocityDispersionScaleFree" />
       <method description="Tabulates the radial velocity dispersion vs. radius for Burkert halos." method="radialVelocityDispersionTabulate" />
       <method description="Tabulates the freefall time vs. freefall radius for Burkert halos." method="freefallTabulate" />
       <method description="Compute the freefall time in a scale-free Burkert halo." method="freefallTimeScaleFree" />
       <method description="Tabulates the radius vs. enclosed density for Burkert halos." method="radiusEnclosingDensityTabulate" />
       <method description="Returns the total angular momentum in an Burkert dark matter profile with given {\normalfont \ttfamily concentration}." method="angularMomentumScaleFree" />
       <method description="Tabulates the specific angular momentum vs. radius in an Burkert profile for rapid inversion." method="inverseAngularMomentum" />
       <method description="Computes the total energy of an Burkert profile halo of given {\normalfont \ttfamily concentration}." method="profileEnergy" />
       <method description="Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position {\normalfont \ttfamily radius} (in units of the scale radius) in an Burkert profile." method="specificAngularMomentumScaleFree" />
       <method description="Tabulate properties of the Burkert halo profile which must be computed numerically." method="tabulate" />
     </methods>
     !!]
     final     ::                                      burkertDestructor
     procedure :: autoHook                          => burkertAutoHook
     procedure :: calculationReset                  => burkertCalculationReset
     procedure :: get                               => burkertGet
     procedure :: density                           => burkertDensity
     procedure :: densityLogSlope                   => burkertDensityLogSlope
     procedure :: enclosedMass                      => burkertEnclosedMass
     procedure :: radiusEnclosingDensity            => burkertRadiusEnclosingDensity
     procedure :: radiusEnclosingDensityTabulate    => burkertRadiusEnclosingDensityTabulate
     procedure :: potential                         => burkertPotential
     procedure :: circularVelocity                  => burkertCircularVelocity
     procedure :: circularVelocityMaximum           => burkertCircularVelocityMaximum
     procedure :: radiusCircularVelocityMaximum     => burkertRadiusCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => burkertRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => burkertRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => burkertRotationNormalization
     procedure :: energy                            => burkertEnergy
     procedure :: kSpace                            => burkertKSpace
     procedure :: freefallRadius                    => burkertFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => burkertFreefallRadiusIncreaseRate
     procedure :: profileEnergy                     => burkertProfileEnergy
     procedure :: specificAngularMomentumScaleFree  => burkertSpecificAngularMomentumScaleFree
     procedure :: angularMomentumScaleFree          => burkertAngularMomentumScaleFree
     procedure :: enclosedMassScaleFree             => burkertEnclosedMassScaleFree
     procedure :: densityScaleFree                  => burkertDensityScaleFree
     procedure :: radialVelocityDispersionScaleFree => burkertRadialVelocityDispersionScaleFree
     procedure :: radialVelocityDispersionTabulate  => burkertRadialVelocityDispersionTabulate
     procedure :: tabulate                          => burkertTabulate
     procedure :: inverseAngularMomentum            => burkertInverseAngularMomentum
     procedure :: freefallTabulate                  => burkertFreefallTabulate
     procedure :: freefallTimeScaleFree             => burkertFreefallTimeScaleFree
     procedure :: radialMoment                      => burkertRadialMoment
  end type darkMatterProfileDMOBurkert

  interface darkMatterProfileDMOBurkert
     !!{
     Constructors for the {\normalfont \ttfamily burkert} dark matter halo profile class.
     !!}
     module procedure burkertConstructorParameters
     module procedure burkertConstructorInternal
  end interface darkMatterProfileDMOBurkert

  ! Number of points per decade of concentration in Burkert tabulations.
  integer         , parameter :: tablePointsPerDecade                        =100
  integer         , parameter :: densityTablePointsPerDecade                 =100
  integer         , parameter :: inverseTablePointsPerDecade                 =100
  integer         , parameter :: freefallTablePointsPerDecade                =100
  integer         , parameter :: radialVelocityDispersionTablePointsPerDecade=100

  ! Indices for tabulated quantities.
  integer         , parameter :: concentrationEnergyIndex                    =  1                 , concetrationRotationNormalizationIndex=2

  ! Minimum (scale-free) freefall time in the Burkert profile.
  double precision, parameter :: freefallTimeScaleFreeMinimum                =sqrt(3.0d0)*Pi/4.0d0

contains

  function burkertConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOBurkert)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass   ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOBurkert(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function burkertConstructorParameters

  function burkertConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type (darkMatterProfileDMOBurkert)                        :: self
    class(darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    self%concentrationPrevious                          =  -1.0d+0
    self%concentrationMinimum                           =   1.0d+0
    self%concentrationMaximum                           =  20.0d+0
    self%densityRadiusMinimum                           =   1.0d-3
    self%freefallRadiusMinimum                          =   1.0d-3
    self%radiusMinimum                                  =   1.0d-3
    self%radialVelocityDispersionRadiusMinimum          =   1.0d-3
    self%densityRadiusMaximum                           =   1.0d+2
    self%freefallRadiusMaximum                          =   1.0d+2
    self%radiusMaximum                                  =   1.0d+2
    self%radialVelocityDispersionRadiusMaximum          =   1.0d+2
    self%burkertDensityTableInitialized                 =  .false.
    self%burkertFreefallTableInitialized                =  .false.
    self%burkertInverseTableInitialized                 =  .false.
    self%burkertTableInitialized                        =  .false.
    self%burkertRadialVelocityDispersionTableInitialized=  .false.
    self%lastUniqueID                                   =  -1
    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                            &
         & call Error_Report                                                                                                 &
         &      (                                                                                                            &
         &       'Burkert dark matter profile requires a dark matter profile component with a gettable "scale" property.'//  &
         &       Component_List(                                                                                             &
         &                      'darkMatterProfile'                                                                        , &
         &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)                &
         &                     )                                                                                         //  &
         &       {introspection:location}                                                                                    &
         &      )
    ! Initialize the tabulations.
    call self%tabulate                        ()
    call self%inverseAngularMomentum          ()
    call self%radialVelocityDispersionTabulate()
    return
  end function burkertConstructorInternal

  subroutine burkertAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOBurkert), intent(inout) :: self

    call calculationResetEvent%attach(self,burkertCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOBurkert')
    return
  end subroutine burkertAutoHook

  subroutine burkertDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOBurkert), intent(inout) :: self

    if (self%burkertFreefallTableInitialized                ) then
       call self%burkertFreeFall                      %destroy()
       call self%burkertFreeFallInverse               %destroy()
       deallocate(self%burkertFreefallInverse               )
    end if
    if (self%burkertDensityTableInitialized                 ) then
       call self%burkertDensityTable                  %destroy()
       call self%burkertDensityTableInverse           %destroy()
       deallocate(self%burkertDensityTableInverse           )
    end if
    if (self%burkertInverseTableInitialized                 ) then
       call self%burkertSpecificAngularMomentum       %destroy()
       call self%burkertSpecificAngularMomentumInverse%destroy()
       deallocate(self%burkertSpecificAngularMomentumInverse)
    end if
    if (self%burkertTableInitialized                        ) then
       call self%burkertConcentrationTable            %destroy()
    end if
    if (self%burkertRadialVelocityDispersionTableInitialized) then
       call self%burkertRadialVelocityDispersionTable %destroy()
    end if
    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    if (calculationResetEvent%isAttached(self,burkertCalculationReset)) call calculationResetEvent%detach(self,burkertCalculationReset)
    return
  end subroutine burkertDestructor

  subroutine burkertCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOBurkert), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node

    self%specificAngularMomentumScalingsComputed=.false.
    self%maximumVelocityComputed                =.false.
    self%lastUniqueID                           =node%uniqueID()
    return
  end subroutine burkertCalculationReset

  function burkertGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic     , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo  , massTypeDark                  , weightByMass
    use :: Mass_Distributions        , only : massDistributionBurkert, kinematicsDistributionBurkert
    implicit none
    class  (massDistributionClass         ), pointer                 :: massDistribution_
    type   (kinematicsDistributionBurkert ), pointer                 :: kinematicsDistribution_
    class  (darkMatterProfileDMOBurkert   ), intent(inout)           :: self
    type   (treeNode                      ), intent(inout)           :: node
    type   (enumerationWeightByType       ), intent(in   ), optional :: weightBy
    integer                                , intent(in   ), optional :: weightIndex
    class  (nodeComponentBasic            ), pointer                 :: basic
    class  (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionBurkert :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionBurkert)
       basic             => node%basic            ()
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionBurkert(                                                                                  &amp;
           &amp;                   mass         =basic            %mass                                      (    ), &amp;
           &amp;                   radiusOuter  =self             %darkMatterHaloScale_%radiusVirial         (node), &amp;
           &amp;                   scaleLength  =darkMatterProfile%scale                                     (    ), &amp;
           &amp;                   componentType=                                       componentTypeDarkHalo      , &amp;
           &amp;                   massType     =                                       massTypeDark                 &amp;
           &amp;                  )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionBurkert( &amp;
	 &amp;                       )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function burkertGet

  subroutine burkertTabulate(self,concentration)
    !!{
    Tabulate properties of the Burkert halo profile which must be computed numerically.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout)           :: self
    double precision                             , intent(in   ), optional :: concentration
    integer                                                                :: iConcentration
    logical                                                                :: retabulate
    double precision                                                       :: tableConcentration

    retabulate=.not.self%burkertTableInitialized
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
       self%burkertTableNumberPoints=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(tablePointsPerDecade))+1
       call self%burkertConcentrationTable%destroy()
       call self%burkertConcentrationTable%create(self%concentrationMinimum,self%concentrationMaximum,self%burkertTableNumberPoints,2)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,self%burkertTableNumberPoints
          tableConcentration=self%burkertConcentrationTable%x(iConcentration)
          call self%burkertConcentrationTable%populate(                   self%profileEnergy           (tableConcentration),iConcentration,table=concentrationEnergyIndex              )
          call self%burkertConcentrationTable%populate(tableConcentration/self%angularMomentumScaleFree(tableConcentration),iConcentration,table=concetrationRotationNormalizationIndex)
       end do
       ! Specify that tabulation has been made.
       self%burkertTableInitialized=.true.
    end if
    return
  end subroutine burkertTabulate

  subroutine burkertInverseAngularMomentum(self,specificAngularMomentum)
    !!{
    Tabulates the specific angular momentum vs. radius in an Burkert profile for rapid inversion.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout)           :: self
    double precision                             , intent(in   ), optional :: specificAngularMomentum
    integer                                                                :: iRadius
    logical                                                                :: retabulate

    retabulate=.not.self%burkertInverseTableInitialized
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
       self%burkertInverseTableNumberPoints=int(log10(self%radiusMaximum/self%radiusMinimum)*dble(inverseTablePointsPerDecade))+1
       ! Create a range of radii.
       call self%burkertSpecificAngularMomentum%destroy(                                                       )
       call self%burkertSpecificAngularMomentum%create (self%radiusMinimum,self%radiusMaximum,self%burkertInverseTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertInverseTableNumberPoints
          call self%burkertSpecificAngularMomentum%populate(                                                                                       &
               &                                            self%specificAngularMomentumScaleFree(self%burkertSpecificAngularMomentum%x(iRadius)), &
               &                                            iRadius                                                                                &
               &                                           )
       end do
       call self%burkertSpecificAngularMomentum%reverse(self%burkertSpecificAngularMomentumInverse)
       ! Specify that tabulation has been made.
       self%burkertInverseTableInitialized=.true.
    end if
    return
  end subroutine burkertInverseAngularMomentum

  double precision function burkertDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: radiusOverScaleRadius      , scaleRadius, &
         &                                                             virialRadiusOverScaleRadius

    basic                       => node             %basic            (                 )
    darkMatterProfile           => node             %darkMatterProfile(autoCreate=.true.)
    scaleRadius                 =  darkMatterProfile%scale            (                 )
    radiusOverScaleRadius       =  radius                                      /scaleRadius
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    burkertDensity              =  self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)*basic%mass()/scaleRadius**3
    return
  end function burkertDensity

  double precision function burkertDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily
    radius} (given in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: radiusOverScaleRadius, scaleRadius
    !$GLC attributes unused :: self

    darkMatterProfile      =>  node             %darkMatterProfile(autoCreate=.true.)
    scaleRadius            =   darkMatterProfile%scale            (                 )
    radiusOverScaleRadius  =  +     radius &
         &                    /scaleRadius
    burkertDensityLogSlope =  -      radiusOverScaleRadius   /(1.0d0+radiusOverScaleRadius   ) &
         &                    -2.0d0*radiusOverScaleRadius**2/(1.0d0+radiusOverScaleRadius**2)
    return
  end function burkertDensityLogSlope

  double precision function burkertEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
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
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    burkertEnclosedMass         =  self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)*basic%mass()
    return
  end function burkertEnclosedMass

  double precision function burkertPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galactic_Structure_Options, only : enumerationStructureErrorCodeType, structureErrorCodeSuccess
    use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    class           (darkMatterProfileDMOBurkert      ), intent(inout)           :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    class           (nodeComponentDarkMatterProfile   )               , pointer  :: darkMatterProfile
    double precision                                   , parameter               :: radiusSmall          =1.0d-10
    double precision                                                             :: radiusOverScaleRadius        , virialRadiusOverScaleRadius

    if (present(status)) status=structureErrorCodeSuccess
    darkMatterProfile   => node%darkMatterProfile(autoCreate=.true.)
    radiusOverScaleRadius            =radius                       /darkMatterProfile%scale()
    virialRadiusOverScaleRadius      =self%darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
    if (radiusOverScaleRadius < radiusSmall) then
       burkertPotential=                                                                          &
            & +(                                                                                  &
            &   -Pi                                                                               &
            &   +2.0d0                                                                            &
            &   /3.0d0                                                                            &
            &   *radiusOverScaleRadius**2                                                         &
            &  )                                                                                  &
            & /(                                                                                  &
            &     -2.0d0                              *atan(      virialRadiusOverScaleRadius   ) &
            &     +2.0d0                              *log (1.0d0+virialRadiusOverScaleRadius   ) &
            &     +                                    log (1.0d0+virialRadiusOverScaleRadius**2) &
            &  )                                                                                  &
            & *virialRadiusOverScaleRadius                                                        &
            & *self%darkMatterHaloScale_%velocityVirial(node)**2
    else
       burkertPotential=                                                                          &
            & +(                                                                                  &
            &   -Pi                                                                               &
            &   *radiusOverScaleRadius                                                            &
            &   +(                                                                                &
            &     +2.0d0*(1.0d0+radiusOverScaleRadius)*atan(      radiusOverScaleRadius         ) &
            &     -2.0d0*(1.0d0+radiusOverScaleRadius)*log (1.0d0+radiusOverScaleRadius         ) &
            &     -      (1.0d0-radiusOverScaleRadius)*log (1.0d0+radiusOverScaleRadius      **2) &
            &    )                                                                                &
            &  )                                                                                  &
            & /radiusOverScaleRadius                                                              &
            & /(                                                                                  &
            &     -2.0d0                              *atan(      virialRadiusOverScaleRadius   ) &
            &     +2.0d0                              *log (1.0d0+virialRadiusOverScaleRadius   ) &
            &     +                                    log (1.0d0+virialRadiusOverScaleRadius**2) &
            &  )                                                                                  &
            & *virialRadiusOverScaleRadius                                                        &
            & *self%darkMatterHaloScale_%velocityVirial(node)**2
    end if
    return
  end function burkertPotential

  double precision function burkertCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius

    if (radius > 0.0d0) then
       burkertCircularVelocity=sqrt(gravitationalConstantGalacticus*self%enclosedMass(node,radius)/radius)
    else
       burkertCircularVelocity=0.0d0
    end if
    return
  end function burkertCircularVelocity

  double precision function burkertCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    ! The radius (in units of the scale radius) at which the rotation speed peaks in a Burkert halo.
    double precision                                , parameter     :: radiusMaximum    =3.2446257246042642d0
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: scaleRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if maximum velocity is already computed. Compute and store if not.
    if (.not.self%maximumVelocityComputed) then
       darkMatterProfile            => node             %darkMatterProfile(autoCreate=.true.                          )
       scaleRadius                  =  darkMatterProfile%scale            (                                           )
       self%maximumVelocityPrevious =  self             %circularVelocity (node             ,radiusMaximum*scaleRadius)
       self%maximumVelocityComputed =  .true.
    end if
    burkertCircularVelocityMaximum=self%maximumVelocityPrevious
    return
  end function burkertCircularVelocityMaximum

  double precision function burkertRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    ! The radius (in units of the scale radius) at which the rotation speed peaks in a Burkert halo.
    double precision                                , parameter     :: radiusMaximum    =3.2446257246042642d0
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    
    darkMatterProfile                    =>  node             %darkMatterProfile(autoCreate=.true.)
    burkertRadiusCircularVelocityMaximum =  +                  radiusMaximum                        &
         &                                  *darkMatterProfile%scale            (                 )
    return
  end function burkertRadiusCircularVelocityMaximum

  double precision function burkertRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout)          :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: radiusOverScaleRadius      , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    darkMatterProfile           => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                 =  darkMatterProfile%scale()
    radiusOverScaleRadius       =  radius                                      /scaleRadius
    virialRadiusOverScaleRadius =  self%darkMatterHaloScale_%radiusVirial(node)/scaleRadius
    if (radius > 0.0d0) then
       call self%radialVelocityDispersionTabulate(radiusOverScaleRadius)
       burkertRadialVelocityDispersion=self%burkertRadialVelocityDispersionTable%interpolate(radiusOverScaleRadius)
    else
       burkertRadialVelocityDispersion=0.0d0
    end if
    ! Compute the normalization factor.
    call burkertMassNormalizationFactor(self,virialRadiusOverScaleRadius)
    ! Evaluate the radial velocity dispersion.
    burkertRadialVelocityDispersion=+burkertRadialVelocityDispersion                &
         &                          *sqrt(                                          &
         &                                +self%burkertNormalizationFactorPrevious  &
         &                                *virialRadiusOverScaleRadius              &
         &                               )                                          &
         &                          *self%darkMatterHaloScale_%velocityVirial(node)
    return
  end function burkertRadialVelocityDispersion

  double precision function burkertRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily
    specificAngularMomentum} (given in units of km s$^{-1}$ Mpc)
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: specificAngularMomentumScaleFree

    ! Return immediately with zero radius for non-positive specific angular momenta.
    if (specificAngularMomentum <= 0.0d0) then
       burkertRadiusFromSpecificAngularMomentum=0.0d0
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
       self%specificAngularMomentumScale=self%specificAngularMomentumLengthScale                              &
            &                            *self%circularVelocity(node,self%specificAngularMomentumLengthScale)
    end if

    ! Compute the specific angular momentum in scale free units (using the scale length for distances and sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/self%specificAngularMomentumScale

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%inverseAngularMomentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    burkertRadiusFromSpecificAngularMomentum=self%burkertSpecificAngularMomentumInverse%interpolate(specificAngularMomentumScaleFree)

    ! Convert to a physical radius.
    burkertRadiusFromSpecificAngularMomentum=burkertRadiusFromSpecificAngularMomentum*self%specificAngularMomentumLengthScale
    return
  end function burkertRadiusFromSpecificAngularMomentum

  double precision function burkertRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: concentration

    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the rotation normalization by interpolation.
    burkertRotationNormalization=+self%burkertConcentrationTable%interpolate(concentration,table=concetrationRotationNormalizationIndex) &
         &                       /self%darkMatterHaloScale_%radiusVirial(node)
    return
  end function burkertRotationNormalization

  double precision function burkertEnergy(self,node)
    !!{
    Return the energy of an Burkert halo density profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (nodeComponentBasic            ), pointer       :: basic
    double precision                                                :: concentration

    ! Get components.
    basic             => node%basic            (                 )
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy by interpolation.
    burkertEnergy=self%burkertConcentrationTable%interpolate(concentration,table=concentrationEnergyIndex) &
         &        *basic%mass()*self%darkMatterHaloScale_%velocityVirial(node)**2
    return
  end function burkertEnergy

  double precision function burkertAngularMomentumScaleFree(self,concentration)
    !!{
    Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    in an Burkert dark matter profile with given {\normalfont \ttfamily concentration}. This is given by:
    \begin{equation}
    J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    \end{equation}
    where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    \begin{equation}
    J = \left. \left[ 2 \tan^{-1} c + 2 \log(1+c) + \log(1+c^2) - 4c \right] \right/ \left[ 2 \tan^{-1} c - 2 \log(1+c) - \log(1+c^2) \right].
    \end{equation}
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: concentration
    !$GLC attributes unused :: self

    burkertAngularMomentumScaleFree=+(                                    &
         &                            +2.0d0*atan(      concentration   ) &
         &                            +2.0d0*log (1.0d0+concentration   ) &
         &                            +      log (1.0d0+concentration**2) &
         &                            -4.0d0*           concentration     &
         &                           )                                    &
         &                          /(                                    &
         &                            +2.0d0*atan(      concentration   ) &
         &                            -2.0d0*log (1.0d0+concentration   ) &
         &                            -      log (1.0d0+concentration**2) &
         &                           )
    return
  end function burkertAngularMomentumScaleFree

  double precision function burkertSpecificAngularMomentumScaleFree(self,radius)
    !!{
    Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    {\normalfont \ttfamily radius} (in units of the scale radius) in an Burkert profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: radius

    burkertSpecificAngularMomentumScaleFree=sqrt(radius*self%enclosedMassScaleFree(radius,1.0d0))
    return
  end function burkertSpecificAngularMomentumScaleFree

  double precision function burkertEnclosedMassScaleFree(self,radius,concentration)
    !!{
    Returns the enclosed mass (in units of the virial mass) in an Burkert dark matter profile with given {\normalfont \ttfamily
    concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: concentration                       , radius
    double precision                             , parameter     :: minimumRadiusForExactSolution=1.0d-4

    if (radius < minimumRadiusForExactSolution) then
       ! Use a series solution for small radii.
       burkertEnclosedMassScaleFree=(4.0d0/3.0d0)*radius**3*(1.0d0-radius)
    else
       ! Use the exact solution.
       burkertEnclosedMassScaleFree=(                             &
            &                        +2.0d0*log (1.0d0+radius   ) &
            &                        +      log (1.0d0+radius**2) &
            &                        -2.0d0*atan(      radius   ) &
            &                       )
    end if
    ! Compute the mass profile normalization factor.
    call burkertMassNormalizationFactor(self,concentration)
    ! Evaluate the scale-free enclosed mass.
    burkertEnclosedMassScaleFree=burkertEnclosedMassScaleFree*self%burkertNormalizationFactorPrevious
    return
  end function burkertEnclosedMassScaleFree

  subroutine burkertMassNormalizationFactor(self,concentration)
    !!{
    Compute the normalization factor for the burkert mass profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: concentration
    ! Precomputed Burkert normalization factor for unit concentration.
    double precision                             , parameter     :: burkertNormalizationFactorUnitConcentration=1.0d0/(3.0d0*log(2.0d0)-2.0d0*atan(1.0d0))

    ! Check if we were called with a different concentration compared to the previous call.
    if (concentration /= self%concentrationPrevious) then
       ! We were, so recompute the normalization factor.
       if (concentration == 1.0d0) then
          self%burkertNormalizationFactorPrevious=burkertNormalizationFactorUnitConcentration
       else
          self%burkertNormalizationFactorPrevious=+1.0d0                                &
               &                                  /(                                    &
               &                                    +2.0d0*log (1.0d0+concentration   ) &
               &                                    +      log (1.0d0+concentration**2) &
               &                                    -2.0d0*atan(      concentration   ) &
               &                                   )
       end if
       self%concentrationPrevious=concentration
    end if
    return
  end subroutine burkertMassNormalizationFactor

  double precision function burkertDensityScaleFree(self,radius,concentration)
    !!{
    Returns the density (in units such that the virial mass and scale length are unity) in an Burkert dark matter profile with
    given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: concentration, radius
    !$GLC attributes unused :: self

    burkertDensityScaleFree=+1.0d0                                &
         &                  /(1.0d0+radius   )                    &
         &                  /(1.0d0+radius**2)                    &
         &                  /Pi                                   &
         &                  /(                                    &
         &                    +2.0d0*log (1.0d0+concentration   ) &
         &                    +      log (1.0d0+concentration**2) &
         &                    -2.0d0*atan(      concentration   ) &
         &                   )
    return
  end function burkertDensityScaleFree

  double precision function burkertProfileEnergy(self,concentration)
    !!{
    Computes the total energy of an Burkert profile halo of given {\normalfont \ttfamily concentration} using the methods of
    \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: concentration
    type            (integrator                 )                :: integratorPotential    , integratorJeans       , &
         &                                                          integratorKinetic
    double precision                                             :: jeansEquationIntegral  , kineticEnergy         , &
         &                                                          kineticEnergyIntegral  , potentialEnergy       , &
         &                                                          potentialEnergyIntegral, radiusMaximum         , &
         &                                                          radiusMinimum          , concentrationParameter

    ! Compute the potential energy.
    radiusMinimum          =0.0d0
    radiusMaximum          =concentration
    concentrationParameter =concentration
    integratorPotential    =integrator(burkertPotentialEnergyIntegrand,toleranceRelative=1.0d-3)
    potentialEnergyIntegral=integratorPotential%integrate(radiusMinimum,radiusMaximum)
    potentialEnergy        =-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)
    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum         =concentration
    radiusMaximum         =100.0d0*concentration
    concentrationParameter=concentration
    integratorJeans       =integrator(burkertJeansEquationIntegrand,toleranceRelative=1.0d-3)
    jeansEquationIntegral =integratorJeans%integrate(radiusMinimum,radiusMaximum)
    ! Compute the kinetic energy.
    radiusMinimum         =0.0d0
    radiusMaximum         =concentration
    concentrationParameter=concentration
    integratorKinetic     =integrator(burkertKineticEnergyIntegrand,toleranceRelative=1.0d-3)
    kineticEnergyIntegral =integratorKinetic%integrate(radiusMinimum,radiusMaximum)
    kineticEnergy         =2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    ! Compute the total energy.
    burkertProfileEnergy  =(potentialEnergy+kineticEnergy)*concentration
    return

  contains

    double precision function burkertPotentialEnergyIntegrand(radius)
      !!{
      Integrand for Burkert profile potential energy.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      burkertPotentialEnergyIntegrand=(self%enclosedMassScaleFree(radius,concentrationParameter)/radius)**2
      return
    end function burkertPotentialEnergyIntegrand

    double precision function burkertKineticEnergyIntegrand(radius)
      !!{
      Integrand for Burkert profile kinetic energy.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      burkertKineticEnergyIntegrand=+self%EnclosedMassScaleFree(radius,concentrationParameter) &
           &                        *self%densityScaleFree     (radius,concentrationParameter) &
           &                        *radius
      return
    end function burkertKineticEnergyIntegrand

    double precision function burkertJeansEquationIntegrand(radius)
      !!{
      Integrand for Burkert profile Jeans equation.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      burkertJeansEquationIntegrand=+self%enclosedMassScaleFree(radius,concentrationParameter) &
           &                        *self%densityScaleFree     (radius,concentrationParameter) &
           &                        /radius**2
      return
    end function burkertJeansEquationIntegrand

  end function burkertProfileEnergy

  double precision function burkertKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the Burkert density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$), using the
    expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    use :: Exponential_Integrals   , only : Cosine_Integral               , Exponential_Integral, Sine_Integral
    use :: Galacticus_Nodes        , only : nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout)          :: self
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
    concentration=self%darkMatterHaloScale_%radiusVirial(node)/radiusScale

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    burkertKSpace=dimag(                                                                                                                                              &
         &              +(                                                                                                                                            &
         &                +       2.0d0         *exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Cosine_Integral     (                             waveNumberScaleFree) &
         &                -       2.0d0         *exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Cosine_Integral     (      (+1.0d0+concentration)*waveNumberScaleFree) &
         &                +dcmplx(1.0d0,1.0d0)                                                                                                                        &
         &                *(                                                                                                                                          &
         &                  -dcmplx(0.0d0,1.0d0)*exp(+                    waveNumberScaleFree)*Pi                                                                     &
         &                  -                    exp(+                    waveNumberScaleFree)*Exponential_Integral(       -1.0d0               *waveNumberScaleFree) &
         &                  +dcmplx(0.0d0,1.0d0)*exp(-                    waveNumberScaleFree)*Exponential_Integral(       +1.0d0               *waveNumberScaleFree) &
         &                  +                    exp(+                    waveNumberScaleFree)*Exponential_Integral(dcmplx(-1.0d0,concentration)*waveNumberScaleFree) &
         &                  -dcmplx(0.0d0,1.0d0)*exp(-                    waveNumberScaleFree)*Exponential_Integral(dcmplx(+1.0d0,concentration)*waveNumberScaleFree) &
         &                  +dcmplx(1.0d0,1.0d0)*exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Sine_Integral       (                             waveNumberScaleFree) &
         &                  -dcmplx(1.0d0,1.0d0)*exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Sine_Integral       (      (+1.0d0+concentration)*waveNumberScaleFree) &
         &                 )                                                                                                                                          &
         &               )                                                                                                                                            &
         &              /(                                                                                                                                            &
         &                +waveNumberScaleFree                                                                                                                        &
         &                *(                                                                                                                                          &
         &                  -2.0d0*atan(      concentration   )                                                                                                       &
         &                  +2.0d0*log (1.0d0+concentration   )                                                                                                       &
         &                  +      log (1.0d0+concentration**2)                                                                                                       &
         &                 )                                                                                                                                          &
         &               )                                                                                                                                            &
         &             )
    return
  end function burkertKSpace

  double precision function burkertFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the Burkert density profile at the specified {\normalfont \ttfamily time} (given in Gyr).
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: time
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile
    double precision                                                        :: concentration    , freefallTimeScaleFree, &
         &                                                                     radiusScale      , timeScale            , &
         &                                                                     velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       burkertFreefallRadius=0.0d0
       return
    end if
    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! Get the scale radius.
    radiusScale=darkMatterProfile%scale()
    ! Get the concentration.
    concentration=+self%darkMatterHaloScale_%radiusVirial  (node) &
         &        /radiusScale
    ! Get the virial velocity.
    velocityScale=+self%darkMatterHaloScale_%velocityVirial(node)
    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                  &
         &    *radiusScale                              &
         &    /velocityScale                            &
         &    /sqrt(                                    &
         &          +                 concentration     &
         &         )                                    &
         &    *sqrt(                                    &
         &          -2.0d0*atan(      concentration   ) &
         &          +2.0d0*log (1.0d0+concentration   ) &
         &          +      log (1.0d0+concentration**2) &
         &         )
    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale
    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)
    ! The freefall time is finite at zero radius in this profile. If the requested time is less than this, return zero radius.
    if (freefallTimeScaleFree < freefallTimeScaleFreeMinimum) then
       burkertFreefallRadius=0.0d0
    else
       burkertFreefallRadius=self%burkertFreefallInverse%interpolate(freefallTimeScaleFree)*radiusScale
    end if
    return
  end function burkertFreefallRadius

  double precision function burkertFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the Burkert density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, treeNode
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: time
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile
    double precision                                                        :: concentration    , freefallTimeScaleFree, &
         &                                                                     radiusScale      , timeScale            , &
         &                                                                     velocityScale

    ! For non-positive freefall times, return the limiting value for small radii.
    if (time <= 0.0d0) then
       burkertFreefallRadiusIncreaseRate=0.0d0
       return
    end if
    ! Get components.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! Get the scale radius.
    radiusScale=darkMatterProfile%scale()
    ! Get the concentration.
    concentration=+self%darkMatterHaloScale_%radiusVirial  (node) &
         &        /radiusScale
    ! Get the virial velocity.
    velocityScale=+self%darkMatterHaloScale_%velocityVirial(node)
    ! Compute time scale.
    timeScale=+Mpc_per_km_per_s_To_Gyr                  &
         &    *radiusScale                              &
         &    /velocityScale                            &
         &    /sqrt(                                    &
         &          +                 concentration     &
         &         )                                    &
         &    *sqrt(                                    &
         &          -2.0d0*atan(      concentration   ) &
         &          +2.0d0*log (1.0d0+concentration   ) &
         &          +      log (1.0d0+concentration**2) &
         &         )
    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale
    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)
    ! The freefall time is finite at zero radius in this profile. If the requested time is less than this, return zero radius.
    if (freefallTimeScaleFree < freefallTimeScaleFreeMinimum) then
       burkertFreefallRadiusIncreaseRate=0.0d0
    else
       burkertFreefallRadiusIncreaseRate=self%burkertFreefallInverse%interpolateGradient(freefallTimeScaleFree)*radiusScale/timeScale
    end if
    return
  end function burkertFreefallRadiusIncreaseRate

  subroutine burkertFreefallTabulate(self,freefallTimeScaleFree)
    !!{
    Tabulates the freefall time vs. freefall radius for Burkert halos.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: freefallTimeScaleFree
    logical                                                      :: retabulate
    integer                                                      :: iRadius
    double precision                                             :: freefallTime

    retabulate=.not.self%burkertFreefallTableInitialized
    ! If the table has not yet been made, compute and store the freefall corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusMinimum)
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
    end if
    do while (freefallTimeScaleFree > self%freefallTimeMaximum)
       self%freefallRadiusMaximum=2.0d0*self%freefallRadiusMaximum
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
       retabulate=.true.
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%burkertFreefallTableNumberPoints=int(log10(self%freefallRadiusMaximum/self%freefallRadiusMinimum)*dble(freefallTablePointsPerDecade))+1
       ! Create the table.
       call self%burkertFreefall%destroy(                                                                                           )
       call self%burkertFreefall%create (self%freefallRadiusMinimum,self%freefallRadiusMaximum,self%burkertFreefallTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertFreefallTableNumberPoints
          freefallTime=self%freefallTimeScaleFree(self%burkertFreefall%x(iRadius))
          if (iRadius > 1) freefallTime=max(freefallTime,self%burkertFreefall%y(iRadius-1))
          call self%burkertFreefall%populate(freefallTime,iRadius)
       end do
       call self%burkertFreefall%reverse(self%burkertFreefallInverse)
       ! Specify that tabulation has been made.
       self%burkertFreefallTableInitialized=.true.
    end if
    return
  end subroutine burkertFreefallTabulate

  double precision function burkertFreefallTimeScaleFree(self,radius)
    !!{
    Compute the freefall time in a scale-free Burkert halo.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: radius
    type            (integrator                 )                :: integrator_
    double precision                                             :: radiusEnd  , radiusStart
    !$GLC attributes unused :: self

    radiusStart                 =radius
    radiusEnd                   =0.0d0
    integrator_                 =integrator           (burkertFreefallTimeScaleFreeIntegrand,toleranceRelative=1.0d-5)
    burkertFreefallTimeScaleFree=integrator_%integrate(radiusEnd                            ,radiusStart             )
    return

  contains

    double precision function burkertFreefallTimeScaleFreeIntegrand(radius)
      !!{
      Integrand function used for finding the free-fall time in Burkert halos.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius
      double precision, parameter     :: radiusSmall        =1.0d-2
      double precision                :: potential                 , potentialStart, &
           &                             potentialDifference

      if (radius < radiusSmall) then
         potential     =+Pi                         &
              &         -2.0d0*radius     **2/3.0d0 &
              &         +      radius     **3/3.0d0
      else
         potential     =-(                                                      &
              &           -Pi   *                               radius          &
              &           +2.0d0*(1.0d0+radius     )*atan(      radius        ) &
              &           -2.0d0*(1.0d0+radius     )*log (1.0d0+radius        ) &
              &           -      (1.0d0-radius     )*log (1.0d0+radius     **2) &
              &          )                                                      &
              &         /              radius
      end if
      if (radiusStart < radiusSmall) then
         potentialStart=+Pi                         &
              &         -2.0d0*radiusStart**2/3.0d0 &
              &         +      radiusStart**3/3.0d0
      else
         potentialStart=-(                                                      &
              &           -Pi   *                               radiusStart     &
              &           +2.0d0*(1.0d0+radiusStart)*atan(      radiusStart   ) &
              &           -2.0d0*(1.0d0+radiusStart)*log (1.0d0+radiusStart   ) &
              &           -      (1.0d0-radiusStart)*log (1.0d0+radiusStart**2) &
              &          )                                                      &
              &         /              radiusStart
      end if
      potentialDifference=+potential-potentialStart
      if (potentialDifference > 0.0d0) then
         burkertFreefallTimeScaleFreeIntegrand=1.0d0/sqrt(2.0d0*potentialDifference)
      else
         burkertFreefallTimeScaleFreeIntegrand=0.0d0
      end if
      return
    end function burkertFreefallTimeScaleFreeIntegrand

  end function burkertFreefallTimeScaleFree

  double precision function burkertRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout)           :: node
    double precision                                , intent(in   )           :: moment
    double precision                                , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    class           (nodeComponentBasic            )               , pointer  :: basic
    class           (nodeComponentDarkMatterProfile)               , pointer  :: darkMatterProfile
    double precision                                                          :: radiusMinimumActual, radiusMaximumActual, &
         &                                                                       radiusScale        , concentration

    basic             =>  node                                  %basic            (                 )
    darkMatterProfile =>  node                                  %darkMatterProfile(autoCreate=.true.)
    radiusScale       =   darkMatterProfile                     %scale            (                 )
    concentration     =  +self             %darkMatterHaloScale_%radiusVirial     (           node  ) &
         &               /radiusScale
    radiusMinimumActual=0.0d0
    radiusMaximumActual=concentration
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum/radiusScale
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum/radiusScale
    burkertRadialMoment=+basic%mass()                        &
         &              *radiusScale**(moment-2.0d0)         &
         &              *(                                   &
         &                +radialMoment(radiusMaximumActual) &
         &                -radialMoment(radiusMinimumActual) &
         &               )
    return

  contains

    double precision function radialMoment(radius)
      !!{
      Evaluate the radial moment in the Burkert profile.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      use :: Numerical_Comparison    , only : Values_Agree
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius

      if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
         radialMoment=+(                                                                                                   &
              &         -2.0d0* log(1.0d0+radius          )                                                                &
              &         +       log(1.0d0+radius       **2)                                                                &
              &         +2.0d0*atan(      radius          )                                                                &
              &        )                                                                                                   &
              &       /4.0d0                                                                                               &
              &       /Pi                                                                                                  &
              &       /(                                                                                                   &
              &         +2.0d0* log(1.0d0+concentration   )                                                                &
              &         +       log(1.0d0+concentration**2)                                                                &
              &         -2.0d0*atan(      concentration   )                                                                &
              &        )
      else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
         radialMoment=+(                                                                                                   &
              &         +2.0d0* log(1.0d0+radius          )                                                                &
              &         +       log(1.0d0+radius       **2)                                                                &
              &         -2.0d0*atan(      radius          )                                                                &
              &        )                                                                                                   &
              &       /4.0d0                                                                                               &
              &       /Pi                                                                                                  &
              &       /(                                                                                                   &
              &         +2.0d0* log(1.0d0+concentration   )                                                                &
              &         +       log(1.0d0+concentration**2)                                                                &
              &         -2.0d0*atan(      concentration   )                                                                &
              &        )
      else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
         radialMoment=+(                                                                                                   &
              &         +4.0d0*           radius                                                                           &
              &         -2.0d0* log(1.0d0+radius          )                                                                &
              &         -       log(1.0d0+radius       **2)                                                                &
              &         -2.0d0*atan(      radius          )                                                                &
              &        )                                                                                                   &
              &       /4.0d0                                                                                               &
              &       /Pi                                                                                                  &
              &       /(                                                                                                   &
              &         +2.0d0* log(1.0d0+concentration   )                                                                &
              &         +       log(1.0d0+concentration**2)                                                                &
              &         -2.0d0*atan(      concentration   )                                                                &
              &        )
      else
         radialMoment=+                                                                          radius**   (1.0d0+moment) &
              &       *(                                                                                                   &
              &         -                                                                        radius                    &
              &         *Hypergeometric_2F1([1.0d0,0.5d0*(2.0d0+moment)],[0.5d0*(4.0d0+moment)],-radius**2)/(2.0d0+moment) &
              &         +Hypergeometric_2F1([1.0d0,0.5d0*(1.0d0+moment)],[0.5d0*(3.0d0+moment)],-radius**2)/(1.0d0+moment) &
              &         +Hypergeometric_2F1([1.0d0,      (1.0d0+moment)],[      (2.0d0+moment)],-radius   )/(1.0d0+moment) &
              &        )                                                                                                   &
              &       /2.0d0                                                                                               &
              &       /Pi                                                                                                  &
              &       /(                                                                                                   &
              &         +2.0d0* log(1.0d0+concentration   )                                                                &
              &         +       log(1.0d0+concentration**2)                                                                &
              &         -2.0d0*atan(      concentration   )                                                                &
              &        )
      end if
      return
    end function radialMoment

  end function burkertRadialMoment

  double precision function burkertRadiusEnclosingDensity(self,node,density)
    !!{
    Null implementation of function to compute the radius enclosing a given density for Burkert dark matter halo profiles.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (darkMatterProfileDMOBurkert   ), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: density
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile
    class           (nodeComponentBasic            ), pointer               :: basic
    double precision                                                        :: densityScaleFree , radiusScale, &
         &                                                                     concentration

    basic             =>  node                                  %basic                (                 )
    darkMatterProfile =>  node                                  %darkMatterProfile    (autoCreate=.true.)
    radiusScale       =   darkMatterProfile                     %scale                (                 )
    concentration     =  +self             %darkMatterHaloScale_%radiusVirial         (            node )      &
         &               /radiusScale
    densityScaleFree  =  +density                                                                              &
         &               *radiusScale                                                                      **3 &
         &               /basic                                 %mass                 (                   )    &
         &               /self                                  %enclosedMassScaleFree(1.0d0,concentration)
    call self%radiusEnclosingDensityTabulate(densityScaleFree)
    burkertRadiusEnclosingDensity=+self%burkertDensityTableInverse%interpolate(-densityScaleFree) &
         &                        *radiusScale
    return
  end function burkertRadiusEnclosingDensity

  subroutine burkertRadiusEnclosingDensityTabulate(self,densityScaleFree)
    !!{
    Tabulates the radius vs. enclosed density for Burkert halos.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout) :: self
    double precision                             , intent(in   ) :: densityScaleFree
    logical                                                      :: retabulate
    integer                                                      :: iRadius

    retabulate=.not.self%burkertDensityTableInitialized
    ! If the table has not yet been made, compute and store the enclosed density corresponding to the minimum and maximum radii
    ! that will be tabulated by default.
    if (retabulate) then
       self%densityMaximum=3.0d0*self%enclosedMassScaleFree(self%densityRadiusMinimum,1.0d0)/self%densityRadiusMinimum**3/4.0d0/Pi
       self%densityMinimum=3.0d0*self%enclosedMassScaleFree(self%densityRadiusMaximum,1.0d0)/self%densityRadiusMaximum**3/4.0d0/Pi
    end if
    do while (densityScaleFree < self%densityMinimum)
       self%densityRadiusMaximum=2.0d0*self%densityRadiusMaximum
       self%densityMinimum      =3.0d0*self%enclosedMassScaleFree(self%densityRadiusMaximum,1.0d0)/self%densityRadiusMaximum**3/4.0d0/Pi
       retabulate               =.true.
    end do
    do while (densityScaleFree > self%densityMaximum)
       self%densityRadiusMinimum=0.5d0*self%densityRadiusMinimum
       self%densityMaximum      =3.0d0*self%enclosedMassScaleFree(self%densityRadiusMinimum,1.0d0)/self%densityRadiusMinimum**3/4.0d0/Pi
       retabulate               =.true.
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%burkertDensityTableNumberPoints=int(log10(self%densityRadiusMaximum/self%densityRadiusMinimum)*dble(densityTablePointsPerDecade))+1
       ! Create the table.
       call self%burkertDensityTable%destroy(                                                                                        )
       call self%burkertDensityTable%create (self%densityRadiusMinimum,self%densityRadiusMaximum,self%burkertDensityTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertDensityTableNumberPoints
          call self%burkertDensityTable%populate(                                                                           &
               &                                 -3.0d0                                                                     &
               &                                 /4.0d0                                                                     &
               &                                 /Pi                                                                        &
               &                                 *self%enclosedMassScaleFree(self%burkertDensityTable%x(iRadius),1.0d0)     &
               &                                 /                           self%burkertDensityTable%x(iRadius)       **3, &
               &                                                                                        iRadius             &
               &                                )
       end do
       call self%burkertDensityTable%reverse(self%burkertDensityTableInverse)
       ! Specify that tabulation has been made.
       self%burkertDensityTableInitialized=.true.
    end if
    return
  end subroutine burkertRadiusEnclosingDensityTabulate

  double precision function burkertRadialVelocityDispersionScaleFree(self,radius)
    !!{
    Compute the radial velocity dispersion in a scale-free Burkert halo.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout)            :: self
    double precision                             , intent(in   )            :: radius
    double precision                                            , parameter :: radiusTiny   =1.0d-9, radiusLarge  =5.0d3
    double precision                                                        :: radiusMinimum       , radiusMaximum
    type            (integrator                 )                           :: integrator_
    !$GLC attributes unused :: self

    radiusMinimum                           =max(       radius,radiusTiny )
    radiusMaximum                           =max(10.0d0*radius,radiusLarge)
    integrator_                             =integrator(burkertJeansEquationIntegrand,toleranceRelative=1.0d-6)
    burkertRadialVelocityDispersionScaleFree=sqrt(                                                    &
         &                                        +integrator_%integrate(radiusMinimum,radiusMaximum) &
         &                                        *(1.0d0+radius   )                                  &
         &                                        *(1.0d0+radius**2)                                  &
         &                                       )
     return

  contains

    double precision function burkertJeansEquationIntegrand(radius)
      !!{
      Integrand for Burkert dark matter profile Jeans equation.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius > 0.0d0) then
         burkertJeansEquationIntegrand=+(                             &
              &                          +2.0d0*log (1.0d0+radius   ) &
              &                          +      log (1.0d0+radius**2) &
              &                          -2.0d0*atan(      radius   ) &
              &                         )                             &
              &                        /(1.0d0+radius   )             &
              &                        /(1.0d0+radius**2)             &
              &                        /       radius**2
      else
         burkertJeansEquationIntegrand=0.0d0
      end if
      return
    end function burkertJeansEquationIntegrand

  end function burkertRadialVelocityDispersionScaleFree

  subroutine burkertRadialVelocityDispersionTabulate(self,radius)
    !!{
    Tabulates the radial velocity dispersion vs. radius for Burkert halos.
    !!}
    implicit none
    class           (darkMatterProfileDMOBurkert), intent(inout)           :: self
    double precision                             , intent(in   ), optional :: radius
    logical                                                                :: retabulate
    integer                                                                :: iRadius

    retabulate=.not.self%burkertRadialVelocityDispersionTableInitialized
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
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%burkertRadialVelocityDispersionTableNumberPoints=int(log10(self%radialVelocityDispersionRadiusMaximum/self%radialVelocityDispersionRadiusMinimum) &
            &                                                    *dble(radialVelocityDispersionTablePointsPerDecade                                         ) &
            &                                                   )+1
       ! Create the table.
       call self%burkertRadialVelocityDispersionTable%destroy(                                                      )
       call self%burkertRadialVelocityDispersionTable%create (self%radialVelocityDispersionRadiusMinimum           ,self%radialVelocityDispersionRadiusMaximum, &
            &                                                 self%burkertRadialVelocityDispersionTableNumberPoints )
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertRadialVelocityDispersionTableNumberPoints
          call self%burkertRadialVelocityDispersionTable%populate(self%radialVelocityDispersionScaleFree(self%burkertRadialVelocityDispersionTable%x(iRadius)),iRadius)
       end do
       ! Specify that tabulation has been made.
       self%burkertRadialVelocityDispersionTableInitialized=.true.
    end if
    return
  end subroutine burkertRadialVelocityDispersionTabulate
