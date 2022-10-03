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
  An implementation of dark matter halo profiles with finite resolution (to mimic the effects of resolution in N-body
  simulations for example).
  !!}

  use :: Cosmology_Functions         , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversType, enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough
  use :: Numerical_Interpolation     , only : interpolator

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOFiniteResolutionNFW">
   <description>
    A dark matter profile DMO class which applies a finite resolution to an NFW density profile, typically to mimic the effects
    of finite resolution in an N-body simulation. Specifically, the density profile is given by
    \begin{equation}
     \rho(r) = \rho_\mathrm{NFW}(r) \left( 1 + \left[ \frac{\Delta x}{r} \right]^2 \right)^{-1/2},
    \end{equation}
    where $\Delta x$ is the larger of the resolution length, {\normalfont \ttfamily [lengthResolution]}, and the radius in the
    original profile enclosing the mass resolution, {\normalfont \ttfamily [massResolution]}.
   </description>

  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOFiniteResolution) :: darkMatterProfileDMOFiniteResolutionNFW
     !!{
     A dark matter halo profile class implementing finiteResolutionNFW dark matter halos.
     !!}
     private
     double precision                                            :: potentialRadiusPrevious                            , potentialPrevious                              , &
          &                                                         velocityDispersionRadialRadiusPrevious             , velocityDispersionRadialPrevious               , &
          &                                                         massNormalizationPrevious                          , lengthResolutionScaleFreePrevious              , &
          &                                                         lengthResolutionScaleFreePreviousSquared           , lengthResolutionScaleFreePreviousCubed         , &
          &                                                         lengthResolutionScaleFreePreviousSqrtTerm          , lengthResolutionScaleFreePreviousSqrt2Term     , &
          &                                                         lengthResolutionScaleFreePreviousSqrtCubedTerm     , lengthResolutionScaleFreePreviousLowerTerm     , &
          &                                                         lengthResolutionScaleFreePreviousOnePlusTerm       , lengthResolutionScaleFreePreviousOnePlus2Term  , &
          &                                                         densityRadiusPrevious                              , densityPrevious                                , &
          &                                                         densityNormalizationPrevious                       , radiusEnclosingDensityDensityPrevious          , &
          &                                                         radiusEnclosingDensityPrevious                     , radiusEnclosingMassMassPrevious                , &
          &                                                         radiusEnclosingMassPrevious                        , energyPrevious
     ! Velocity dispersion tabulation.
     logical                                                     :: velocityDispersionRadialTableInitialized
     integer                                                     :: velocityDispersionRadialTableRadiusCoreCount       , velocityDispersionRadialTableRadiusCount
     double precision              , allocatable, dimension(:  ) :: velocityDispersionRadialTableRadiusCore            , velocityDispersionRadialTableRadius
     double precision              , allocatable, dimension(:,:) :: velocityDispersionRadialTable
     type            (interpolator), allocatable                 :: velocityDispersionRadialTableRadiusCoreInterpolator, velocityDispersionRadialTableRadiusInterpolator
     double precision                                            :: velocityDispersionRadialRadiusMinimum              , velocityDispersionRadialRadiusMaximum          , &
          &                                                         velocityDispersionRadialRadiusCoreMinimum          , velocityDispersionRadialRadiusCoreMaximum
     ! Radius-enclosing-density tabulation.
     logical                                                     :: radiusEnclosingDensityTableInitialized
     integer                                                     :: radiusEnclosingDensityTableRadiusCoreCount         , radiusEnclosingDensityTableDensityCount
     double precision              , allocatable, dimension(:  ) :: radiusEnclosingDensityTableRadiusCore              , radiusEnclosingDensityTableDensity
     double precision              , allocatable, dimension(:,:) :: radiusEnclosingDensityTable
     type            (interpolator), allocatable                 :: radiusEnclosingDensityTableRadiusCoreInterpolator  , radiusEnclosingDensityTableDensityInterpolator
     double precision                                            :: radiusEnclosingDensityDensityMinimum               , radiusEnclosingDensityDensityMaximum           , &
          &                                                         radiusEnclosingDensityRadiusCoreMinimum            , radiusEnclosingDensityRadiusCoreMaximum
     ! Radius-enclosing-mass tabulation.
     logical                                                     :: radiusEnclosingMassTableInitialized
     integer                                                     :: radiusEnclosingMassTableRadiusCoreCount            , radiusEnclosingMassTableMassCount
     double precision              , allocatable, dimension(:  ) :: radiusEnclosingMassTableRadiusCore                 , radiusEnclosingMassTableMass
     double precision              , allocatable, dimension(:,:) :: radiusEnclosingMassTable
     type            (interpolator), allocatable                 :: radiusEnclosingMassTableRadiusCoreInterpolator     , radiusEnclosingMassTableMassInterpolator
     double precision                                            :: radiusEnclosingMassMassMinimum                     , radiusEnclosingMassMassMaximum                 , &
          &                                                         radiusEnclosingMassRadiusCoreMinimum               , radiusEnclosingMassRadiusCoreMaximum
     ! Energy tabulation.
     logical                                                     :: energyTableInitialized
     integer                                                     :: energyTableRadiusCoreCount                         , energyTableConcentrationCount
     double precision              , allocatable, dimension(:  ) :: energyTableRadiusCore                              , energyTableConcentration
     double precision              , allocatable, dimension(:,:) :: energyTable
     type            (interpolator), allocatable                 :: energyTableRadiusCoreInterpolator                  , energyTableConcentrationInterpolator
     double precision                                            :: energyConcentrationMinimum                         , energyConcentrationMaximum                     , &
          &                                                         energyRadiusCoreMinimum                            , energyRadiusCoreMaximum
  contains
     !![
     <methods>
       <method method="velocityDispersionRadialTabulate" description="Tabulate the enclosed mass as a function of radius and core radius."                               />
       <method method="radiusEnclosingDensityTabulate"   description="Tabulate the radius enclosing a given density as a function of density and core radius."           />
       <method method="radiusEnclosingMassTabulate"      description="Tabulate the radius enclosing a given mass as a function of density and core radius."              />
       <method method="energyTabulate"                   description="Tabulate the energy as a function of concentration and core radius."                               />
       <method method="densityScaleFree"                 description="The density of the profile in units where the mass and scale length are both 1."                   />
       <method method="massEnclosedScaleFree"            description="The mass enclosed of the profile in units where the mass and scale length are both 1."             />
       <method method="storeVelocityDispersionTable"     description="Store the tabulated velocity dispersion to file."                                                  />
       <method method="restoreVelocityDispersionTable"   description="Attempt to restore the tabulated velocity dispersion from file, returning true if successful."     />
       <method method="storeDensityTable"                description="Store the tabulated radius-enclosing-density to file."                                             />
       <method method="restoreDensityTable"              description="Attempt to restore the tabulated radius-enclosing-density from file, returning true if successful."/>
       <method method="storeMassTable"                   description="Store the tabulated radius-enclosing-mass to file."                                                />
       <method method="restoreMassTable"                 description="Attempt to restore the tabulated radius-enclosing-mass from file, returning true if successful."   />
       <method method="storeEnergyTable"                 description="Store the tabulated energy to file."                                                               />
       <method method="restoreEnergyTable"               description="Attempt to restore the tabulated energy from file, returning true if successful."                  />
     </methods>
     !!]
     procedure :: autoHook                         => finiteResolutionNFWAutoHook
     procedure :: calculationReset                 => finiteResolutionNFWCalculationReset
     procedure :: density                          => finiteResolutionNFWDensity
     procedure :: enclosedMass                     => finiteResolutionNFWEnclosedMass
     procedure :: potential                        => finiteResolutionNFWPotential
     procedure :: radiusEnclosingDensity           => finiteResolutionNFWRadiusEnclosingDensity
     procedure :: radiusEnclosingMass              => finiteResolutionNFWRadiusEnclosingMass
     procedure :: energy                           => finiteResolutionNFWEnergy
     procedure :: radialVelocityDispersion         => finiteResolutionNFWRadialVelocityDispersion
     procedure :: velocityDispersionRadialTabulate => finiteResolutionNFWVelocityDispersionRadialTabulate
     procedure :: radiusEnclosingDensityTabulate   => finiteResolutionNFWRadiusEnclosingDensityTabulate
     procedure :: radiusEnclosingMassTabulate      => finiteResolutionNFWRadiusEnclosingMassTabulate
     procedure :: energyTabulate                   => finiteResolutionNFWEnergyTabulate
     procedure :: densityScaleFree                 => finiteResolutionNFWDensityScaleFree
     procedure :: massEnclosedScaleFree            => finiteResolutionNFWMassEnclosedScaleFree
     procedure :: storeVelocityDispersionTable     => finiteResolutionNFWStoreVelocityDispersionTable
     procedure :: restoreVelocityDispersionTable   => finiteResolutionNFWRestoreVelocityDispersionTable
     procedure :: storeDensityTable                => finiteResolutionNFWStoreDensityTable
     procedure :: restoreDensityTable              => finiteResolutionNFWRestoreDensityTable
     procedure :: storeMassTable                   => finiteResolutionNFWStoreMassTable
     procedure :: restoreMassTable                 => finiteResolutionNFWRestoreMassTable
     procedure :: storeEnergyTable                 => finiteResolutionNFWStoreEnergyTable
     procedure :: restoreEnergyTable               => finiteResolutionNFWRestoreEnergyTable
  end type darkMatterProfileDMOFiniteResolutionNFW

  interface darkMatterProfileDMOFiniteResolutionNFW
     !!{
     Constructors for the {\normalfont \ttfamily finiteResolutionNFW} dark matter halo profile class.
     !!}
     module procedure finiteResolutionNFWConstructorParameters
     module procedure finiteResolutionNFWConstructorInternal
  end interface darkMatterProfileDMOFiniteResolutionNFW

  ! Tabulation resolution parameters.
  integer, parameter :: velocityDispersionRadialTableRadiusPointsPerDecade    =100
  integer, parameter :: velocityDispersionRadialTableRadiusCorePointsPerDecade=100
  integer, parameter :: radiusEnclosingDensityTableDensityPointsPerDecade     =100
  integer, parameter :: radiusEnclosingDensityTableRadiusCorePointsPerDecade  =100
  integer, parameter :: radiusEnclosingMassTableMassPointsPerDecade           =100
  integer, parameter :: radiusEnclosingMassTableRadiusCorePointsPerDecade     =100
  integer, parameter :: energyTableConcentrationPointsPerDecade               =100
  integer, parameter :: energyTableRadiusCorePointsPerDecade                  =100

  ! Sub-module-scope variables used in integrations.
  class  (darkMatterProfileDMOFiniteResolutionNFW), pointer :: self_
  integer                                                   :: iRadiusCore_, iDensity_, iMass_
  !$omp threadprivate(self_,iRadiusCore_,iDensity_,iMass_)
  
contains

  function finiteResolutionNFWConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily finiteResolutionNFW} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOFiniteResolutionNFW)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    double precision                                                         :: lengthResolution    , massResolution
    type            (varying_string                         )                :: nonAnalyticSolver
    logical                                                                  :: resolutionIsComoving

    !![
    <inputParameter>
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>The resolution length, $\Delta x$.</description>
    </inputParameter>
    <inputParameter>
      <name>massResolution</name>
      <source>parameters</source>
      <description>The resolution mass, $\Delta M$.</description>
    </inputParameter>
    <inputParameter>
      <name>resolutionIsComoving</name>
      <source>parameters</source>
      <description>If true, the resolution length is assumed to be fixed in comoving coordinates, otherwise in physical coordinates.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !!]
    self=darkMatterProfileDMOFiniteResolutionNFW(lengthResolution,massResolution,resolutionIsComoving,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterHaloScale_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="cosmologyFunctions_"  />
    !!]
    return
  end function finiteResolutionNFWConstructorParameters

  function finiteResolutionNFWConstructorInternal(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterHaloScale_,cosmologyFunctions_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily finiteResolutionNFW} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMOFiniteResolutionNFW)                        :: self
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    double precision                                         , intent(in   )         :: lengthResolution    , massResolution
    type            (enumerationNonAnalyticSolversType      ), intent(in   )         :: nonAnalyticSolver
    logical                                                  , intent(in   )         :: resolutionIsComoving
    !![
    <constructorAssign variables="lengthResolution, massResolution, resolutionIsComoving, nonAnalyticSolver, *darkMatterHaloScale_, *cosmologyFunctions_"/>
    !!]
    
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    self%lastUniqueID                                  =-1_kind_int8
    self%genericLastUniqueID                           =-1_kind_int8
    self%lengthResolutionPrevious                      =-huge(0.0d0)
    self%massNormalizationPrevious                     =-huge(0.0d0)
    self%enclosedMassPrevious                          =-huge(0.0d0)
    self%enclosedMassRadiusPrevious                    =-huge(0.0d0)
    self%potentialPrevious                             =-huge(0.0d0)
    self%potentialRadiusPrevious                       =-huge(0.0d0)
    self%velocityDispersionRadialPrevious              =-huge(0.0d0)
    self%velocityDispersionRadialRadiusPrevious        =-huge(0.0d0)
    self%lengthResolutionScaleFreePrevious             =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousSquared      =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousCubed        =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousOnePlusTerm  =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousOnePlus2Term =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousSqrtTerm     =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousSqrt2Term    =-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousSqrtCubedTerm=-huge(0.0d0)
    self%lengthResolutionScaleFreePreviousLowerTerm    =-huge(0.0d0)
    self%densityRadiusPrevious                         =-huge(0.0d0)
    self%densityPrevious                               =-huge(0.0d0)
    self%densityNormalizationPrevious                  =-huge(0.0d0)
    self%radiusEnclosingDensityDensityPrevious         =-huge(0.0d0)
    self%radiusEnclosingDensityPrevious                =-huge(0.0d0)
    self%radiusEnclosingMassMassPrevious               =-huge(0.0d0)
    self%radiusEnclosingMassPrevious                   =-huge(0.0d0)
    self%energyPrevious                                =+huge(0.0d0)
    ! Velocity dispersion table initialization.
    self%velocityDispersionRadialRadiusMinimum         =+huge(0.0d0)
    self%velocityDispersionRadialRadiusMaximum         =-huge(0.0d0)
    self%velocityDispersionRadialRadiusCoreMinimum     =+huge(0.0d0)
    self%velocityDispersionRadialRadiusCoreMaximum     =-huge(0.0d0)
    self%velocityDispersionRadialTableInitialized      =.false.
    ! Radius enclosing density table initialization.
    self%radiusEnclosingDensityDensityMinimum          =+huge(0.0d0)
    self%radiusEnclosingDensityDensityMaximum          =-huge(0.0d0)
    self%radiusEnclosingDensityRadiusCoreMinimum       =+huge(0.0d0)
    self%radiusEnclosingDensityRadiusCoreMaximum       =-huge(0.0d0)
    self%radiusEnclosingDensityTableInitialized        =.false.
    ! Radius enclosing mass table initialization.
    self%radiusEnclosingMassMassMinimum                =+huge(0.0d0)
    self%radiusEnclosingMassMassMaximum                =-huge(0.0d0)
    self%radiusEnclosingMassRadiusCoreMinimum          =+huge(0.0d0)
    self%radiusEnclosingMassRadiusCoreMaximum          =-huge(0.0d0)
    self%radiusEnclosingMassTableInitialized           =.false.
    ! Energy table initialization.
    self%energyConcentrationMinimum                    =+huge(0.0d0)
    self%energyConcentrationMaximum                    =-huge(0.0d0)
    self%energyRadiusCoreMinimum                       =+huge(0.0d0)
    self%energyRadiusCoreMaximum                       =-huge(0.0d0)
    self%energyTableInitialized                        =.false.
    allocate(darkMatterProfileDMONFW :: self%darkMatterProfileDMO_)
    select type (darkMatterProfileDMO_ => self%darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       !![
       <referenceConstruct owner="self" isResult="yes" object="darkMatterProfileDMO_" nameAssociated="darkMatterProfileDMO_">
        <constructor>
         darkMatterProfileDMONFW(                                                                &amp;
            &amp;                velocityDispersionUseSeriesExpansion=.false.                  , &amp;
            &amp;                darkMatterHaloScale_                =self%darkMatterHaloScale_  &amp;
            &amp;               )
        </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function finiteResolutionNFWConstructorInternal

  subroutine finiteResolutionNFWAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self

    call calculationResetEvent%attach(self,finiteResolutionNFWCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOFiniteResolutionNFW')
    return
  end subroutine finiteResolutionNFWAutoHook

  subroutine finiteResolutionNFWCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (treeNode                               ), intent(inout) :: node

    call self%darkMatterProfileDMOFiniteResolution%calculationReset(node)
    self%potentialPrevious                     =-huge(0.0d0)
    self%potentialRadiusPrevious               =-huge(0.0d0)
    self%velocityDispersionRadialPrevious      =-huge(0.0d0)
    self%velocityDispersionRadialRadiusPrevious=-huge(0.0d0)
    self%massNormalizationPrevious             =-huge(0.0d0)
    self%densityRadiusPrevious                 =-huge(0.0d0)
    self%densityPrevious                       =-huge(0.0d0)
    self%densityNormalizationPrevious          =-huge(0.0d0)
    self%radiusEnclosingDensityDensityPrevious =-huge(0.0d0)
    self%radiusEnclosingDensityPrevious        =-huge(0.0d0)
    self%radiusEnclosingMassMassPrevious       =-huge(0.0d0)
    self%radiusEnclosingMassPrevious           =-huge(0.0d0)
    self%energyPrevious                        =+huge(0.0d0)
    return
  end subroutine finiteResolutionNFWCalculationReset  

  double precision function finiteResolutionNFWDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot/\mathrm{MPc}^3$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    double precision                                         , intent(in   ) :: radius
    class           (nodeComponentBasic                     ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer       :: darkMatterProfile
    double precision                                                         :: concentration            , radiusScaleFree, &
         &                                                                      lengthResolutionScaleFree
    
    if (node%uniqueID() /= self%lastUniqueID         ) call self%calculationReset(node)
    if (     radius     /= self%densityRadiusPrevious) then
       darkMatterProfile          => node%darkMatterProfile       (    )
       radiusScaleFree            =       radius                        /darkMatterProfile%scale()
       lengthResolutionScaleFree  =  self%lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%densityRadiusPrevious =       radius
       if (self%densityNormalizationPrevious < 0.0d0) then
          basic                             =>  node %basic                            (    )
          concentration                     =   self %darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
          self%densityNormalizationPrevious =  +basic                     %mass        (    )    &
               &                               /darkMatterProfile         %scale       (    )**3 &
               &                               /(                                                &
               &                                 -          concentration                        &
               &                                 /   (1.0d0+concentration)                       &
               &                                 +log(1.0d0+concentration)                       &
               &                                )                                                &
               &                               /4.0d0                                            &
               &                               /Pi
       end if
       self%densityPrevious=+self%densityNormalizationPrevious                                            &
            &               *self%densityScaleFree            (radiusScaleFree,lengthResolutionScaleFree)
    end if
    finiteResolutionNFWDensity=self%densityPrevious
    return
  end function finiteResolutionNFWDensity
       
  double precision function finiteResolutionNFWDensityScaleFree(self,radius,radiusCore)
    !!{
    Returns the sclae-free density in the dark matter profile at the given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius, radiusCore
    
    finiteResolutionNFWDensityScaleFree=1.0d0/(1.0d0+radius)**2/sqrt(radius**2+radiusCore**2)
    return
  end function finiteResolutionNFWDensityScaleFree
  
  double precision function finiteResolutionNFWEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc). The analytic solution (computed using Mathematica) is
    \begin{equation}
    M(x) = M \frac{-\frac{\sqrt{x^2+X^2}}{(1+x) \left(1+X^2\right)}+\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)+\frac{\left(1+2X^2\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{1+X^2} \sqrt{x^2+X^2}}\right)}{\left(1+X^2\right)^{3/2}} -\frac{\left(1 + 2 X^2\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{1 + X^2}}\right)}{\left(1+ X^2\right)^{3/2}}+\frac{\sqrt{X^2}}{1 + X^2}}{\log (1+c)-\frac{c}{1+c}},
    \end{equation}
    where $x=r/r_\mathrm{s}$, $X = \Delta x/r_\mathrm{s}$, and $r_\mathrm{s}$ is the NFW scale length.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    double precision                                         , intent(in   ) :: radius
    class           (nodeComponentBasic                     ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer       :: darkMatterProfile
    double precision                                                         :: concentration            , radiusScaleFree, &
         &                                                                      lengthResolutionScaleFree
    
    if (node%uniqueID() /= self%lastUniqueID              ) call self%calculationReset(node)
    if (     radius     /= self%enclosedMassRadiusPrevious) then
       darkMatterProfile               => node%darkMatterProfile       (    )
       radiusScaleFree                 =       radius                        /darkMatterProfile%scale()
       lengthResolutionScaleFree       =  self%lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%enclosedMassRadiusPrevious =       radius
       if (self%massNormalizationPrevious < 0.0d0) then
          basic                          =>  node %basic                            (    )
          concentration                  =   self %darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
          self%massNormalizationPrevious =  +basic                     %mass        (    ) &
               &                            /(                                             &
               &                              -          concentration                     &
               &                              /   (1.0d0+concentration)                    &
               &                              +log(1.0d0+concentration)                    &
               &                              )
       end if
       self%enclosedMassPrevious=+self%massNormalizationPrevious                                            &
            &                    *self%massEnclosedScaleFree    (radiusScaleFree,lengthResolutionScaleFree)
    end if
    finiteResolutionNFWEnclosedMass=self%enclosedMassPrevious
    return
  end function finiteResolutionNFWEnclosedMass

  double precision function finiteResolutionNFWRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    type            (treeNode                               ), intent(inout), target :: node
    double precision                                         , intent(in   )         :: density
    class           (nodeComponentBasic                     ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer               :: darkMatterProfile
    double precision                                         , parameter             :: epsilonDensity           =1.0d-3
    double precision                                                                 :: concentration                   , densityScaleFree       , &
         &                                                                              lengthResolutionScaleFree       , densityScaleFreeMaximum
    integer         (c_size_t                               ), dimension(0:1)        :: jRadiusCore
    double precision                                         , dimension(0:1)        :: hRadiusCore
    integer                                                                          :: iRadiusCore
    
    if (node%uniqueID() /= self%lastUniqueID              ) call self%calculationReset(node)
    if (     density    /= self%radiusEnclosingDensityDensityPrevious) then
       basic                                      =>  node%basic                                       (    )
       darkMatterProfile                          =>  node%darkMatterProfile                           (    )
       concentration                              =  self%darkMatterHaloScale_%radiusVirial            (node)/darkMatterProfile%scale()
       lengthResolutionScaleFree                  =  self                     %lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%radiusEnclosingDensityDensityPrevious =                            density
       if (self%massNormalizationPrevious < 0.0d0) then
          basic                          =>  node %basic                            (    )
          concentration                  =   self %darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
          self%massNormalizationPrevious =  +basic                     %mass        (    ) &
               &                            /(                                             &
               &                              -          concentration                     &
               &                              /   (1.0d0+concentration)                    &
               &                              +log(1.0d0+concentration)                    &
               &                              )
       end if
       ! Find scale free density, and the maximum such density reached in the profile.
       densityScaleFree       =+                  density                        &
            &                  /self             %massNormalizationPrevious      &
            &                  *darkMatterProfile%scale                    ()**3
       densityScaleFreeMaximum=+1.0d0                                            &
            &                  /4.0d0                                            &
            &                  /Pi                                               &
            &                  /lengthResolutionScaleFree
       if      (densityScaleFree >= densityScaleFreeMaximum) then
          ! Maximum density is exceeded - return zero radius.
          self%radiusEnclosingDensityPrevious=0.0d0
       else if (densityScaleFree >= densityScaleFreeMaximum*(1.0d0-epsilonDensity)) then
          ! For densities close to the maximum density, use a series solution.
          self%radiusEnclosingDensityPrevious=+0.5d0                     &
          &                                   *(                         &
          &                                     +1.0d0                   &
          &                                     -densityScaleFree        &
          &                                     /densityScaleFreeMaximum &
          &                                    )                         &
          &                                   *darkMatterProfile%scale()
       else
          ! Use a tabulated solution in other regimes.   
          ! Ensure table is sufficiently extensive.
          call self%radiusEnclosingDensityTabulate(densityScaleFree,lengthResolutionScaleFree)
          ! Interpolate to get the scale free radius enclosing the scale free density.
          call self%radiusEnclosingDensityTableRadiusCoreInterpolator%linearFactors(lengthResolutionScaleFree,jRadiusCore(0),hRadiusCore)
          jRadiusCore(1)=jRadiusCore(0)+1
          self%radiusEnclosingDensityPrevious=0.0d0
          do iRadiusCore=0,1
             self%radiusEnclosingDensityPrevious=+self%radiusEnclosingDensityPrevious                                                                                                            &
                  &                              +self%radiusEnclosingDensityTableDensityInterpolator%interpolate(densityScaleFree,self%radiusEnclosingDensityTable(:,jRadiusCore(iRadiusCore))) &
                  &                              *                                                                                                                    hRadiusCore(iRadiusCore)
         end do
          self%radiusEnclosingDensityPrevious=+self             %radiusEnclosingDensityPrevious   &
               &                              *darkMatterProfile%scale                         ()
       end if
    end if
    finiteResolutionNFWRadiusEnclosingDensity=self%radiusEnclosingDensityPrevious
    return
  end function finiteResolutionNFWRadiusEnclosingDensity
  
  subroutine finiteResolutionNFWRadiusEnclosingDensityTabulate(self,density,radiusCore)
    !!{
    Tabulates the radius enclosing a given density for finite resolution NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    double precision                                         , intent(in   )         :: density                , radiusCore
    double precision                                         , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                          :: retabulate
    integer                                                                          :: iRadiusCore            , iDensity                , &
         &                                                                              i
    type            (rootFinder                             )                        :: finder

    do i=1,2
       retabulate=.false.
       if (.not.self%radiusEnclosingDensityTableInitialized) then
          retabulate=.true.
       else if (                                                           &
            &    density    < self%radiusEnclosingDensityDensityMinimum    &
            &   .or.                                                       &
            &    density    > self%radiusEnclosingDensityDensityMaximum    &
            &   .or.                                                       &
            &    radiusCore < self%radiusEnclosingDensityRadiusCoreMinimum &
            &   .or.                                                       &
            &    radiusCore > self%radiusEnclosingDensityRadiusCoreMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreDensityTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%radiusEnclosingDensityDensityMinimum      =min(0.5d0*density   ,self%radiusEnclosingDensityDensityMinimum   )
       self%radiusEnclosingDensityDensityMaximum      =max(2.0d0*density   ,self%radiusEnclosingDensityDensityMaximum   )
       self%radiusEnclosingDensityRadiusCoreMinimum   =min(0.5d0*radiusCore,self%radiusEnclosingDensityRadiusCoreMinimum)
       self%radiusEnclosingDensityRadiusCoreMaximum   =max(2.0d0*radiusCore,self%radiusEnclosingDensityRadiusCoreMaximum)
       self%radiusEnclosingDensityTableDensityCount   =int(log10(self%radiusEnclosingDensityDensityMaximum   /self%radiusEnclosingDensityDensityMinimum   )*dble(radiusEnclosingDensityTableDensityPointsPerDecade   ))+1
       self%radiusEnclosingDensityTableRadiusCoreCount=int(log10(self%radiusEnclosingDensityRadiusCoreMaximum/self%radiusEnclosingDensityRadiusCoreMinimum)*dble(radiusEnclosingDensityTableRadiusCorePointsPerDecade))+1
       if (allocated(self%radiusEnclosingDensityTableDensity)) then
          deallocate(self%radiusEnclosingDensityTableRadiusCore)
          deallocate(self%radiusEnclosingDensityTableDensity   )
          deallocate(self%radiusEnclosingDensityTable          )
       end if
       allocate(self%radiusEnclosingDensityTableRadiusCore(                                             self%radiusEnclosingDensityTableRadiusCoreCount))
       allocate(self%radiusEnclosingDensityTableDensity   (self%radiusEnclosingDensityTableDensityCount                                                ))
       allocate(self%radiusEnclosingDensityTable          (self%radiusEnclosingDensityTabledensityCount,self%radiusEnclosingDensityTableRadiusCoreCount))
       ! Create a range of radii and core radii.
       self%radiusEnclosingDensityTableDensity   =Make_Range(self%radiusEnclosingDensityDensityMinimum   ,self%radiusEnclosingDensityDensityMaximum   ,self%radiusEnclosingDensityTableDensityCount   ,rangeType=rangeTypeLogarithmic)
       self%radiusEnclosingDensityTableRadiusCore=Make_Range(self%radiusEnclosingDensityRadiusCoreMinimum,self%radiusEnclosingDensityRadiusCoreMaximum,self%radiusEnclosingDensityTableRadiusCoreCount,rangeType=rangeTypeLogarithmic)
       ! Initialize our root finder.
       finder=rootFinder(                                                             &
            &            rootFunction                 =rootDensity                  , &
            &            toleranceAbsolute            =toleranceAbsolute            , &
            &            toleranceRelative            =toleranceRelative            , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &           )
       ! Loop over density and core radius and populate tables.
       self_ => self
       do iRadiusCore=1,self%radiusEnclosingDensityTableRadiusCoreCount
          iRadiusCore_=iRadiusCore
          do iDensity=1,self%radiusEnclosingDensityTableDensityCount
             iDensity_=iDensity
             if (self%radiusEnclosingDensityTableDensity(iDensity) > 1.0d0/self%radiusEnclosingDensityTableRadiusCore(iRadiusCore)/4.0d0/Pi) then
                ! Density exceeds the maximum density in the profile - so set zero radius.
                self%radiusEnclosingDensityTable(iDensity,iRadiusCore)=0.0d0
             else
                self%radiusEnclosingDensityTable(iDensity,iRadiusCore)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%radiusEnclosingDensityTableRadiusCoreInterpolator)) deallocate(self%radiusEnclosingDensityTableRadiusCoreInterpolator)
       if (allocated(self%radiusEnclosingDensityTableDensityInterpolator   )) deallocate(self%radiusEnclosingDensityTableDensityInterpolator   )
       allocate(self%radiusEnclosingDensityTableRadiusCoreInterpolator)
       allocate(self%radiusEnclosingDensityTableDensityInterpolator   )
       self%radiusEnclosingDensityTableRadiusCoreInterpolator=interpolator(self%radiusEnclosingDensityTableRadiusCore)
       self%radiusEnclosingDensityTableDensityInterpolator   =interpolator(self%radiusEnclosingDensityTableDensity   )
       ! Specify that tabulation has been made.
       self%radiusEnclosingDensityTableInitialized=.true.
       call self%storeDensityTable()
    end if
    return
  end subroutine finiteResolutionNFWRadiusEnclosingDensityTabulate

  double precision function rootDensity(radius)
    !!{
    Root function used in finding the radius enclosing a given mean density.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                                                                                         &
         &      *self_%massEnclosedScaleFree             (radius,self_%radiusEnclosingDensityTableRadiusCore(iRadiusCore_))    &
         &      /4.0d0                                                                                                         &
         &      /Pi                                                                                                            &
         &      /                                         radius                                                           **3 &
         &      -self_%radiusEnclosingDensityTableDensity(                                                   iDensity_    )
    return
  end function rootDensity

  double precision function finiteResolutionNFWRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    type            (treeNode                               ), intent(inout), target :: node
    double precision                                         , intent(in   )         :: mass
    class           (nodeComponentBasic                     ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer               :: darkMatterProfile
    double precision                                                                 :: concentration            , massScaleFree, &
         &                                                                              lengthResolutionScaleFree
    integer         (c_size_t                               ), dimension(0:1)        :: jRadiusCore
    double precision                                         , dimension(0:1)        :: hRadiusCore
    integer                                                                          :: iRadiusCore
    
    if (node%uniqueID() /= self%lastUniqueID                   ) call self%calculationReset(node)
    if (     mass       /= self%radiusEnclosingMassMassPrevious) then
       basic                                =>  node%basic                                       (    )
       darkMatterProfile                    =>  node%darkMatterProfile                           (    )
       concentration                        =  self%darkMatterHaloScale_%radiusVirial            (node)/darkMatterProfile%scale()
       lengthResolutionScaleFree            =  self                     %lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%radiusEnclosingMassMassPrevious =                            mass
       if (self%massNormalizationPrevious < 0.0d0) then
          basic                          =>  node %basic                            (    )
          concentration                  =   self %darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
          self%massNormalizationPrevious =  +basic                     %mass        (    ) &
               &                            /(                                             &
               &                              -          concentration                     &
               &                              /   (1.0d0+concentration)                    &
               &                              +log(1.0d0+concentration)                    &
               &                              )
       end if
       ! Find scale free mass, and the maximum such mass reached in the profile.
       massScaleFree=+     mass                      &
            &        /self%massNormalizationPrevious
       ! Ensure table is sufficiently extensive.
       call self%radiusEnclosingMassTabulate(massScaleFree,lengthResolutionScaleFree)
       ! Interpolate to get the scale free radius enclosing the scale free mass.
       call self%radiusEnclosingMassTableRadiusCoreInterpolator%linearFactors(lengthResolutionScaleFree,jRadiusCore(0),hRadiusCore)
       jRadiusCore(1)=jRadiusCore(0)+1
       self%radiusEnclosingMassPrevious=0.0d0
       do iRadiusCore=0,1
          self%radiusEnclosingMassPrevious=+self%radiusEnclosingMassPrevious                                                                                                   &
               &                           +self%radiusEnclosingMassTableMassInterpolator%interpolate(massScaleFree,self%radiusEnclosingMassTable(:,jRadiusCore(iRadiusCore))) &
               &                           *                                                                                                        hRadiusCore(iRadiusCore)
       end do
       self%radiusEnclosingMassPrevious=+self             %radiusEnclosingMassPrevious   &
            &                           *darkMatterProfile%scale                      ()
    end if
    finiteResolutionNFWRadiusEnclosingMass=self%radiusEnclosingMassPrevious
    return
  end function finiteResolutionNFWRadiusEnclosingMass
  
  subroutine finiteResolutionNFWRadiusEnclosingMassTabulate(self,mass,radiusCore)
    !!{
    Tabulates the radius enclosing a given mass for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    double precision                                         , intent(in   )         :: mass                   , radiusCore
    double precision                                         , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                          :: retabulate
    integer                                                                          :: iRadiusCore            , iMass                   , &
         &                                                                              i
    type            (rootFinder                             )                        :: finder
    
    do i=1,2
       retabulate=.false.
       if (.not.self%radiusEnclosingMassTableInitialized) then
          retabulate=.true.
       else if (                                                        &
            &    mass       < self%radiusEnclosingMassMassMinimum       &
            &   .or.                                                    &
            &    mass       > self%radiusEnclosingMassMassMaximum       &
            &   .or.                                                    &
            &    radiusCore < self%radiusEnclosingMassRadiusCoreMinimum &
            &   .or.                                                    &
            &    radiusCore > self%radiusEnclosingMassRadiusCoreMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreMassTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%radiusEnclosingMassMassMinimum         =min(0.5d0*mass      ,self%radiusEnclosingMassMassMinimum      )
       self%radiusEnclosingMassMassMaximum         =max(2.0d0*mass      ,self%radiusEnclosingMassMassMaximum      )
       self%radiusEnclosingMassRadiusCoreMinimum   =min(0.5d0*radiusCore,self%radiusEnclosingMassRadiusCoreMinimum)
       self%radiusEnclosingMassRadiusCoreMaximum   =max(2.0d0*radiusCore,self%radiusEnclosingMassRadiusCoreMaximum)
       self%radiusEnclosingMassTableMassCount      =int(log10(self%radiusEnclosingMassMassMaximum      /self%radiusEnclosingMassMassMinimum      )*dble(radiusEnclosingMassTableMassPointsPerDecade      ))+1
       self%radiusEnclosingMassTableRadiusCoreCount=int(log10(self%radiusEnclosingMassRadiusCoreMaximum/self%radiusEnclosingMassRadiusCoreMinimum)*dble(radiusEnclosingMassTableRadiusCorePointsPerDecade))+1
       if (allocated(self%radiusEnclosingMassTableMass)) then
          deallocate(self%radiusEnclosingMassTableRadiusCore)
          deallocate(self%radiusEnclosingMassTableMass      )
          deallocate(self%radiusEnclosingMassTable          )
       end if
       allocate(self%radiusEnclosingMassTableRadiusCore(                                       self%radiusEnclosingMassTableRadiusCoreCount))
       allocate(self%radiusEnclosingMassTableMass      (self%radiusEnclosingMassTableMassCount                                             ))
       allocate(self%radiusEnclosingMassTable          (self%radiusEnclosingMassTablemassCount,self%radiusEnclosingMassTableRadiusCoreCount))
       ! Create a range of radii and core radii.
       self%radiusEnclosingMassTableMass      =Make_Range(self%radiusEnclosingMassMassMinimum      ,self%radiusEnclosingMassMassMaximum      ,self%radiusEnclosingMassTableMassCount      ,rangeType=rangeTypeLogarithmic)
       self%radiusEnclosingMassTableRadiusCore=Make_Range(self%radiusEnclosingMassRadiusCoreMinimum,self%radiusEnclosingMassRadiusCoreMaximum,self%radiusEnclosingMassTableRadiusCoreCount,rangeType=rangeTypeLogarithmic)
       ! Initialize our root finder.
       finder=rootFinder(                                                             &
            &            rootFunction                 =rootMass                     , &
            &            toleranceAbsolute            =toleranceAbsolute            , &
            &            toleranceRelative            =toleranceRelative            , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &           )
       ! Loop over mass and core radius and populate tables.
       self_ => self
       do iRadiusCore=1,self%radiusEnclosingMassTableRadiusCoreCount
          iRadiusCore_=iRadiusCore
          do iMass=1,self%radiusEnclosingMassTableMassCount
             iMass_=iMass
             ! Check that the root condition is satisfied at infinitely large radius. If it is not, then no radius encloses the
             ! required mass. Simply set the radius to an infinitely large value in such case.
             if (rootMass(radius=huge(0.0d0)) < 0.0d0) then
                self%radiusEnclosingMassTable(iMass,iRadiusCore)=huge(0.0d0)
             else
                self%radiusEnclosingMassTable(iMass,iRadiusCore)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%radiusEnclosingMassTableRadiusCoreInterpolator)) deallocate(self%radiusEnclosingMassTableRadiusCoreInterpolator)
       if (allocated(self%radiusEnclosingMassTableMassInterpolator      )) deallocate(self%radiusEnclosingMassTableMassInterpolator   )
       allocate(self%radiusEnclosingMassTableRadiusCoreInterpolator)
       allocate(self%radiusEnclosingMassTableMassInterpolator      )
       self%radiusEnclosingMassTableRadiusCoreInterpolator=interpolator(self%radiusEnclosingMassTableRadiusCore)
       self%radiusEnclosingMassTableMassInterpolator      =interpolator(self%radiusEnclosingMassTableMass      )
       ! Specify that tabulation has been made.
       self%radiusEnclosingMassTableInitialized=.true.
       call self%storeMassTable()
    end if
    return
  end subroutine finiteResolutionNFWRadiusEnclosingMassTabulate

  double precision function rootMass(radius)
    !!{
    Root function used in finding the radius enclosing a given mean mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+self_%massEnclosedScaleFree       (radius,self_%radiusEnclosingMassTableRadiusCore(iRadiusCore_)) &
         &   -self_%radiusEnclosingMassTableMass(                                                   iMass_    )
    return
  end function rootMass

  double precision function finiteResolutionNFWEnergy(self,node)
    !!{
    Returns the energy (in $M_\odot$ km$^2$/s$^2$) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout)  :: self
    type            (treeNode                               ), intent(inout)  :: node
    class           (nodeComponentBasic                     ), pointer        :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer        :: darkMatterProfile
    double precision                                                          :: concentration    , lengthResolutionScaleFree
    integer         (c_size_t                               ), dimension(0:1) :: jRadiusCore
    double precision                                         , dimension(0:1) :: hRadiusCore
    integer                                                                   :: iRadiusCore
    
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    if (self%energyPrevious > 0.0d0) then
       basic                     =>  node%basic                                       (    )
       darkMatterProfile         =>  node%darkMatterProfile                           (    )
       concentration             =  self%darkMatterHaloScale_%radiusVirial            (node)/darkMatterProfile%scale()
       lengthResolutionScaleFree =  self                     %lengthResolutionPhysical(node)/darkMatterProfile%scale()
       if (self%massNormalizationPrevious < 0.0d0) then
          basic                          =>  node %basic                            (    )
          concentration                  =   self %darkMatterHaloScale_%radiusVirial(node)/darkMatterProfile%scale()
          self%massNormalizationPrevious =  +basic                     %mass        (    ) &
               &                            /(                                             &
               &                              -          concentration                     &
               &                              /   (1.0d0+concentration)                    &
               &                              +log(1.0d0+concentration)                    &
               &                              )
       end if
       ! Ensure table is sufficiently extensive.
       call self%energyTabulate(concentration,lengthResolutionScaleFree)
       ! Interpolate to get the scale free energy.
       call self%energyTableRadiusCoreInterpolator%linearFactors(lengthResolutionScaleFree,jRadiusCore(0),hRadiusCore)
       jRadiusCore(1)=jRadiusCore(0)+1
       self%energyPrevious=0.0d0
       do iRadiusCore=0,1
          self%energyPrevious=+self%energyPrevious                                                                                               &
               &              +self%energyTableConcentrationInterpolator%interpolate(concentration,self%energyTable(:,jRadiusCore(iRadiusCore))) &
               &              *                                                                                       hRadiusCore(iRadiusCore)
       end do
       self%energyPrevious=+self             %energyPrevious                 &
            &              *gravitationalConstantGalacticus                  &
            &              *self             %massNormalizationPrevious  **2 &
            &              /darkMatterProfile%scale                    ()
    end if
    finiteResolutionNFWEnergy=self%energyPrevious
    return
  end function finiteResolutionNFWEnergy
  
  subroutine finiteResolutionNFWEnergyTabulate(self,concentration,radiusCore)
    !!{
    Tabulates the energy for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    double precision                                         , intent(in   )         :: concentration              , radiusCore
    double precision                                         , parameter             :: multiplierRadius   =100.0d0
    type            (integrator                             )                        :: integratorPotential        , integratorKinetic  , &
         &                                                                              integratorPressure
    double precision                                                                 :: pseudoPressure             , energyKinetic      , &
         &                                                                              energyPotential            , concentration_
    logical                                                                          :: retabulate
    integer                                                                          :: iRadiusCore                , iConcentration     , &
         &                                                                              i

    do i=1,2
       retabulate=.false.
       if (.not.self%energyTableInitialized) then
          retabulate=.true.
       else if (                                                 &
            &    concentration < self%energyConcentrationMinimum &
            &   .or.                                             &
            &    concentration > self%energyConcentrationMaximum &
            &   .or.                                             &
            &    radiusCore   < self%energyRadiusCoreMinimum     &
            &   .or.                                             &
            &    radiusCore   > self%energyRadiusCoreMaximum     &
            &  ) then
          retabulate=.true.
       end if
       if (     retabulate.and.i==1) call self%restoreEnergyTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%energyConcentrationMinimum   =min(0.5d0*concentration,self%energyConcentrationMinimum)
       self%energyConcentrationMaximum   =max(2.0d0*concentration,self%energyConcentrationMaximum)
       self%energyRadiusCoreMinimum      =min(0.5d0*radiusCore   ,self%energyRadiusCoreMinimum   )
       self%energyRadiusCoreMaximum      =max(2.0d0*radiusCore   ,self%energyRadiusCoreMaximum   )
       self%energyTableConcentrationCount=int(log10(self%energyConcentrationMaximum/self%energyConcentrationMinimum)*dble(energyTableConcentrationPointsPerDecade))+1
       self%energyTableRadiusCoreCount   =int(log10(self%energyRadiusCoreMaximum   /self%energyRadiusCoreMinimum   )*dble(energyTableRadiusCorePointsPerDecade   ))+1
       if (allocated(self%energyTableConcentration)) then
          deallocate(self%energyTableRadiusCore   )
          deallocate(self%energyTableConcentration)
          deallocate(self%energyTable             )
       end if
       allocate(self%energyTableRadiusCore   (                                   self%energyTableRadiusCoreCount))
       allocate(self%energyTableConcentration(self%energyTableConcentrationCount                                ))
       allocate(self%energyTable             (self%energyTableconcentrationCount,self%energyTableRadiusCoreCount))
       ! Create a range of radii and core radii.
       self%energyTableConcentration=Make_Range(self%energyConcentrationMinimum,self%energyConcentrationMaximum,self%energyTableConcentrationCount,rangeType=rangeTypeLogarithmic)
       self%energyTableRadiusCore   =Make_Range(self%energyRadiusCoreMinimum   ,self%energyRadiusCoreMaximum   ,self%energyTableRadiusCoreCount   ,rangeType=rangeTypeLogarithmic)
       ! Initialize integrators.
       integratorPotential=integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
       integratorKinetic  =integrator(integrandEnergyKinetic  ,toleranceRelative=1.0d-3)
       integratorPressure =integrator(integrandPseudoPressure ,toleranceRelative=1.0d-3)
       ! Loop over concentration and core radius and populate tables.
       self_ => self
       do iRadiusCore=1,self%energyTableRadiusCoreCount
          iRadiusCore_=iRadiusCore
          do iConcentration=1,self%energyTableConcentrationCount
             concentration_=self%energyTableConcentration(iConcentration)
             energyPotential                             =+integratorPotential%integrate(        0.0d0,                  concentration_)
             energyKinetic                               =+integratorKinetic  %integrate(        0.0d0,                  concentration_)/4.0d0/Pi
             pseudoPressure                              =+integratorPressure %integrate(concentration_,multiplierRadius*concentration_)/4.0d0/Pi
             self%energyTable(iConcentration,iRadiusCore)=-0.5d0                                                                                   &
                  &                                       *(                                                                                       &
                  &                                         +energyPotential                                                                       &
                  &                                         +self%massEnclosedScaleFree(concentration_,self%energyTableRadiusCore(iRadiusCore))**2 &
                  &                                         /concentration_                                                                        &
                  &                                        )                                                                                       &
                  &                                       +2.0d0                                                                                   &
                  &                                       *Pi                                                                                      &
                  &                                       *(                                                                                       &
                  &                                         +concentration_                                                                    **3 &
                  &                                         *pseudoPressure                                                                        &
                  &                                         +energyKinetic                                                                         &
                  &                                        )
            end do
       end do
       ! Build interpolators.
       if (allocated(self%energyTableRadiusCoreInterpolator   )) deallocate(self%energyTableRadiusCoreInterpolator   )
       if (allocated(self%energyTableConcentrationInterpolator)) deallocate(self%energyTableConcentrationInterpolator)
       allocate(self%energyTableRadiusCoreInterpolator   )
       allocate(self%energyTableConcentrationInterpolator)
       self%energyTableRadiusCoreInterpolator   =interpolator(self%energyTableRadiusCore   )
       self%energyTableConcentrationInterpolator=interpolator(self%energyTableConcentration)
       ! Specify that tabulation has been made.
       self%energyTableInitialized=.true.
       call self%storeEnergyTable()
    end if
    return
  end subroutine finiteResolutionNFWEnergyTabulate

  double precision function integrandEnergyPotential(radius)
    !!{
    Integrand for potential energy of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandEnergyPotential=(                                                                               &
            &                    +self_%massEnclosedScaleFree(radius,self_%energyTableRadiusCore(iRadiusCore_)) &
            &                    /                            radius                                            &
            &                   )**2
    else
       integrandEnergyPotential=+0.0d0
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
       integrandEnergyKinetic=+self_%massEnclosedScaleFree(radius,self_%energyTableRadiusCore(iRadiusCore_)) &
            &                 *self_%densityScaleFree     (radius,self_%energyTableRadiusCore(iRadiusCore_)) &
            &                 *                            radius
    else
       integrandEnergyKinetic=+0.0d0
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
       integrandPseudoPressure=+self_%massEnclosedScaleFree(radius,self_%energyTableRadiusCore(iRadiusCore_))    &
            &                  *self_%densityScaleFree     (radius,self_%energyTableRadiusCore(iRadiusCore_))    &
            &                  /                            radius                                           **2
    else
       integrandPseudoPressure=+0.0d0
    end if
    return
  end function integrandPseudoPressure

  double precision function finiteResolutionNFWMassEnclosedScaleFree(self,radiusScaleFree,lengthResolutionScaleFree)
    !!{
    Returns the scale-free enclosed mass in the dark matter profile at the given {\normalfont \ttfamily radius}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    double precision                                         , intent(in   ) :: radiusScaleFree                 , lengthResolutionScaleFree
    double precision                                         , parameter     :: radiusScaleFreeSmall     =1.0d-3, radiusScaleFreeLarge     =1.0d4, &
         &                                                                      radiusScaleFreeLargeATanh=1.0d+6
    double precision                                                         :: radiusScaleFreeEffective        , arctanhTerm                    , &
         &                                                                      arctanhTerm1
    
    if (lengthResolutionScaleFree /= self%lengthResolutionScaleFreePrevious) then
       self%lengthResolutionScaleFreePrevious             =lengthResolutionScaleFree
       self%lengthResolutionScaleFreePreviousSquared      =lengthResolutionScaleFree**2
       self%lengthResolutionScaleFreePreviousCubed        =lengthResolutionScaleFree**3
       self%lengthResolutionScaleFreePreviousOnePlusTerm  =+1.0d0+      self%lengthResolutionScaleFreePreviousSquared
       self%lengthResolutionScaleFreePreviousOnePlus2Term =+1.0d0+2.0d0*self%lengthResolutionScaleFreePreviousSquared
       self%lengthResolutionScaleFreePreviousSqrtTerm     =sqrt(self%lengthResolutionScaleFreePreviousOnePlusTerm )
       self%lengthResolutionScaleFreePreviousSqrt2Term    =sqrt(self%lengthResolutionScaleFreePreviousOnePlus2Term)
       self%lengthResolutionScaleFreePreviousSqrtCubedTerm=self%lengthResolutionScaleFreePreviousSqrtTerm**3
       ! For large values of the argument to arctanh(), use a series solution to avoiding floating point errors.
       if (lengthResolutionScaleFree > radiusScaleFreeLargeATanh) then
          arctanhTerm=-log(                           &
               &           +2.0d0                     &
               &           *lengthResolutionScaleFree &
               &          )                           &
               &      /2.0d0                          &
               &      +1.0d0                          &
               &      /2.0d0                          &
               &      /lengthResolutionScaleFree      &
               &      +1.0d0                          &
               &      /8.0d0                          &
               &      /lengthResolutionScaleFree**2
       else
          arctanhTerm=+atanh(                                                &
               &             +(+1.0d0-lengthResolutionScaleFree)             &
               &             /self%lengthResolutionScaleFreePreviousSqrtTerm &
               &            ) 
       end if
       self%lengthResolutionScaleFreePreviousLowerTerm    =+            lengthResolutionScaleFree                      &
            &                                              /       self%lengthResolutionScaleFreePreviousOnePlusTerm   &
            &                                              +       2.0d0                                               &
            &                                              *       self%lengthResolutionScaleFreePreviousOnePlus2Term  &
            &                                              *       arctanhTerm                                         &
            &                                              /       self%lengthResolutionScaleFreePreviousSqrtCubedTerm
    end if
    if (radiusScaleFree < radiusScaleFreeSmall) then
       ! Series expansion for small radii.
       finiteResolutionNFWMassEnclosedScaleFree=+  radiusScaleFree**3                                                                                                                          &
            &                                   *(                                                                                                                                             &
            &                                                                                                                       +1.0d0 /self%lengthResolutionScaleFreePrevious     / 3.0d0 &
            &                                     +radiusScaleFree   *(                                                             +1.0d0 /self%lengthResolutionScaleFreePrevious     / 2.0d0 &
            &                                     +radiusScaleFree   * ( 1.0d0+(+6.0d0*self%lengthResolutionScaleFreePreviousSquared-1.0d0)/self%lengthResolutionScaleFreePreviousCubed/10.0d0 &
            &                                     +radiusScaleFree   *  (1.0d0-(+4.0d0*self%lengthResolutionScaleFreePreviousSquared-1.0d0)/self%lengthResolutionScaleFreePreviousCubed/ 6.0d0 &
            &                                                           )                                                                                                                      &
            &                                                          )                                                                                                                       &
            &                                                         )                                                                                                                        &
            &                                    )
    else
       ! Full analytic solution.
       !! Limit the evaluation to some large radius.
       radiusScaleFreeEffective=min(radiusScaleFree,radiusScaleFreeLarge)
       if (radiusScaleFreeEffective > radiusScaleFreeLargeATanh*self%lengthResolutionScaleFreePrevious) then
          arctanhTerm1=+log  (                                               &
               &              +4.0d0                                         &
               &              *      radiusScaleFreeEffective**2             &
               &              /self%lengthResolutionScaleFreePreviousSquared &
               &             )                                               &
               &       /2.0d0                                                &
               &       -self%lengthResolutionScaleFreePreviousSquared        &
               &       /8.0d0                                                &
               &       /     radiusScaleFreeEffective**2
       else
          arctanhTerm1=+atanh(                                                                                  &
            &                       +radiusScaleFreeEffective                                                   &
            &                 /sqrt(+radiusScaleFreeEffective**2+self%lengthResolutionScaleFreePreviousSquared) &
            &                )
       end if
       finiteResolutionNFWMassEnclosedScaleFree=-         sqrt(+radiusScaleFreeEffective**2+self%lengthResolutionScaleFreePrevious**2) &
            &                                   /(+1.0d0+radiusScaleFreeEffective)                                                     &
            &                                   /self%lengthResolutionScaleFreePreviousOnePlusTerm                                     &
            &                                   -2.0d0                                                                                 &
            &                                   *self%lengthResolutionScaleFreePreviousOnePlus2Term                                    &
            &                                   *atanh(                                                                                &
            &                                          +(                                                                              &
            &                                            +1.0d0                                                                        &
            &                                            +radiusScaleFreeEffective                                                     &
            &                                            -sqrt(+radiusScaleFreeEffective**2+self%lengthResolutionScaleFreePrevious**2) &
            &                                           )                                                                              &
            &                                          /self%lengthResolutionScaleFreePreviousSqrtTerm                                 &
            &                                         )                                                                                &
            &                                   /self%lengthResolutionScaleFreePreviousSqrtCubedTerm                                   &
            &                                   +arctanhTerm1                                                                          &
            &                                   +self%lengthResolutionScaleFreePreviousLowerTerm
       !! Beyond the limiting radius assume logarithmic growth in mass as appropriate for an r⁻³ profile.
       if (radiusScaleFree > radiusScaleFreeEffective)                                           &
            & finiteResolutionNFWMassEnclosedScaleFree=+finiteResolutionNFWMassEnclosedScaleFree &
            &                                          *log(                                     &
            &                                               +radiusScaleFree                     &
            &                                               /radiusScaleFreeEffective            &
            &                                              )
    end if
    return
  end function finiteResolutionNFWMassEnclosedScaleFree
  
  double precision function finiteResolutionNFWPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc). The analytic solution (computed using Mathematica) is
    \begin{eqnarray}
    \Phi(x) &=& -\frac{\mathrm{G} M}{r_\mathrm{s}}  \nonumber \\
            & & \left\{ +\frac{\sqrt{x^2+X^2}}{x \left(X^2+1\right)} \right. \nonumber \\
            & & -\frac{X^2 \log \left(\sqrt{X^2+1} \sqrt{x^2+X^2}-x+X^2\right)}{\left(X^2+1\right)^{3/2}} \nonumber \\
            & & -\frac{\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)}{x} \nonumber \\
            & & -\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{X^2+1} \sqrt{x^2+X^2}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
            & & -\frac{\sqrt{X^2}}{x \left(X^2+1\right)}+\frac{X^2 \log (x+1)}{\left(X^2+1\right)^{3/2}} \nonumber \\
            & & +\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{X^2+1}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
            & & \left. +\frac{ \left(\sqrt{X^2+1}-X^2 \log \left(\sqrt{X^2+1}-1\right)\right)}{\left(X^2+1\right)^{3/2}} \right\} \nonumber \\
            & & /\left[\log (1+c)-\frac{c}{1+c}\right]
    \end{eqnarray}
    !!}
    use :: Galactic_Structure_Options      , only : enumerationStructureErrorCodeType, structureErrorCodeSuccess
    use :: Galacticus_Nodes                , only : nodeComponentBasic              , nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout)           :: self
    type            (treeNode                               ), intent(inout), target   :: node
    double precision                                         , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType      ), intent(  out), optional :: status
    class           (nodeComponentBasic                     ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer                 :: darkMatterProfile
    double precision                                                                   :: concentration                   , radiusScaleFree, &
         &                                                                                lengthResolutionScaleFree
    double precision                                         , parameter               :: radiusScaleFreeSmall     =1.0d-3

    if (present(status)) status=structureErrorCodeSuccess
    if (node%uniqueID() /= self%lastUniqueID              ) call self%calculationReset(node)
    if (     radius     /= self%potentialRadiusPrevious) then
       basic                        =>  node%basic                                        (    )
       darkMatterProfile            =>  node%darkMatterProfile                            (    )
       radiusScaleFree              =                             radius                        /darkMatterProfile%scale()
       concentration                =   self%darkMatterHaloScale_%radiusVirial            (node)/darkMatterProfile%scale()
       lengthResolutionScaleFree    =   self                     %lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%potentialRadiusPrevious =                             radius
       if (radiusScaleFree < radiusScaleFreeSmall) then
          ! Series expansion for small radii.
          self%potentialPrevious       =  -gravitationalConstantGalacticus                              &
               &                          *basic            %mass ()                                    &
               &                          /darkMatterProfile%scale()                                    &
               &                          *(                                                            &
               &                            +(+1.0d0-lengthResolutionScaleFree   )                      &
               &                            /(+1.0d0+lengthResolutionScaleFree**2)                      &
               &                            +        lengthResolutionScaleFree**2                       &
               &                            *(                                                          &
               &                              +asinh(lengthResolutionScaleFree   )                      &
               &                              +log  (                                                   &
               &                                     +(1.0d0+sqrt(+1.0d0+lengthResolutionScaleFree**2)) &
               &                                     /                   lengthResolutionScaleFree      &
               &                                    )                                                   &
               &                             )                                                          &
               &                            /(+1.0d0+lengthResolutionScaleFree**2)**1.5d0               &
               &                            -        radiusScaleFree**2                                 &
               &                            *(+1.0d0-radiusScaleFree             )                      &
               &                            /        lengthResolutionScaleFree                          &
               &                            /6.0d0                                                      &
               &                           )                                                            &
               &                          /(                                                            &
               &                            -          concentration                                    &
               &                            /   (1.0d0+concentration)                                   &
               &                            +log(1.0d0+concentration)                                   &
               &                           )
       else
          self%potentialPrevious       =  -gravitationalConstantGalacticus                                                                                          &
               &                          *basic            %mass ()                                                                                                &
               &                          /darkMatterProfile%scale()                                                                                                &
               &                          *(                                                                                                                        &
               &                            +       sqrt(                          lengthResolutionScaleFree**2                                           )         &
               &                            /                   radiusScaleFree                                /(+1.0d0+      lengthResolutionScaleFree**2)         &
               &                            -       sqrt(      +radiusScaleFree**2+lengthResolutionScaleFree**2                                           )         &
               &                            /                   radiusScaleFree                                /(+1.0d0+      lengthResolutionScaleFree**2)         &
               &                            -                                                                   (+1.0d0+2.0d0*lengthResolutionScaleFree**2)         &
               &                            *atanh(                                                                                                                 &
               &                                   +sqrt(                         +lengthResolutionScaleFree**2/(+1.0d0+      lengthResolutionScaleFree**2)       ) &
               &                                  )                                                                                                                 &
               &                            /                   radiusScaleFree                                /(+1.0d0+      lengthResolutionScaleFree**2)**1.5d0  &
               &                            +atanh(                                                                                                                 &
               &                                   +            radiusScaleFree                                                                                     &
               &                                   /sqrt(      +radiusScaleFree**2+lengthResolutionScaleFree**2                                           )         &
               &                                  )                                                                                                                 &
               &                            /                   radiusScaleFree                                                                                     &
               &                            +                                                                   (+1.0d0+2.0d0*lengthResolutionScaleFree**2)         &
               &                            *atanh(                                                                                                                 &
               &                                        (      -radiusScaleFree   +lengthResolutionScaleFree**2                                           )         &
               &                                   /sqrt                                                        (+1.0d0+      lengthResolutionScaleFree**2)         &
               &                                   /sqrt(      +radiusScaleFree**2+lengthResolutionScaleFree**2                                           )         &
               &                                  )                                                                                                                 &
               &                            /                   radiusScaleFree                                /(+1.0d0+      lengthResolutionScaleFree**2)**1.5d0  &
               &                            -                                      lengthResolutionScaleFree**2                                                     &
               &                            *log  (                                                                                                                 &
               &                                         +1.0d0+radiusScaleFree                                                                                     &
               &                                  )                                                                                                                 &
               &                            /                                                                   (+1.0d0+      lengthResolutionScaleFree**2)**1.5d0  &
               &                            +                                      lengthResolutionScaleFree**2/(+1.0d0+      lengthResolutionScaleFree**2)**1.5d0  &
               &                            *log  (                                                                                                                 &
               &                                               -radiusScaleFree   +lengthResolutionScaleFree**2                                                     &
               &                                   +sqrt(+1.0d0                   +lengthResolutionScaleFree**2                                           )         &
               &                                   *sqrt(      +radiusScaleFree**2+lengthResolutionScaleFree**2                                           )         &
               &                                  )                                                                                                                 &
               &                            +(                                                                                                                      &
               &                              +     sqrt(+1.0d0                   +lengthResolutionScaleFree**2                                           )         &
               &                              -                                    lengthResolutionScaleFree**2                                                     &
               &                              *log(                                                                                                                 &
               &                                         -1.0d0                                                                                                     &
               &                                   +sqrt(+1.0d0                   +lengthResolutionScaleFree**2                                           )         &
               &                                  )                                                                                                                 &
               &                             )                                                                                                                      &
               &                            /                                                                   (+1.0d0+      lengthResolutionScaleFree**2)**1.5d0  &
               &                           )                                                                                                                        &
               &                          /(                                                                                                                        &
               &                            -          concentration                                                                                                &
               &                            /   (1.0d0+concentration)                                                                                               &
               &                            +log(1.0d0+concentration)                                                                                               &
               &                           )
       end if
    end if
    finiteResolutionNFWPotential=self%potentialPrevious
    return
  end function finiteResolutionNFWPotential

  double precision function finiteResolutionNFWRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout)  :: self
    type            (treeNode                               ), intent(inout)  :: node
    double precision                                         , intent(in   )  :: radius
    class           (nodeComponentBasic                     ), pointer        :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer        :: darkMatterProfile
    double precision                                                          :: concentration                        , radiusScaleFree         , &
         &                                                                       lengthResolutionScaleFree            , radiusScaleFreeEffective
    double precision                                         , parameter      :: lengthResolutionScaleFreeSmall=1.0d-3
    integer         (c_size_t                               ), dimension(0:1) :: jRadiusCore
    double precision                                         , dimension(0:1) :: hRadiusCore
    integer                                                                   :: iRadiusCore
    
    if (node%uniqueID() /= self%lastUniqueID              ) call self%calculationReset(node)
    if (     radius     /= self%velocityDispersionRadialRadiusPrevious) then
       basic                                       =>  node%basic                                       (    )
       darkMatterProfile                           =>  node%darkMatterProfile                           (    )
       radiusScaleFree                             =                            radius                        /darkMatterProfile%scale()
       concentration                               =  self%darkMatterHaloScale_%radiusVirial            (node)/darkMatterProfile%scale()
       lengthResolutionScaleFree                   =  self                     %lengthResolutionPhysical(node)/darkMatterProfile%scale()
       self%velocityDispersionRadialRadiusPrevious =                            radius
       ! Compute the effective radius. In the core of the profile the velocity dispersion must become constant. Therefore, we
       ! limit the smallest radius we consider to a small fraction of the core radius. Below this radius a constant velocity
       ! dispersion is assumed.
       radiusScaleFreeEffective=max(radiusScaleFree,lengthResolutionScaleFreeSmall*lengthResolutionScaleFree)
       ! Ensure table is sufficiently extensive.
       call self%velocityDispersionRadialTabulate(radiusScaleFreeEffective,lengthResolutionScaleFree)
       ! Interpolate to get the velocity dispersion.
       call self%velocityDispersionRadialTableRadiusCoreInterpolator%linearFactors(lengthResolutionScaleFree,jRadiusCore(0),hRadiusCore)
       jRadiusCore(1)=jRadiusCore(0)+1
       self%velocityDispersionRadialPrevious=0.0d0
       do iRadiusCore=0,1
          self%velocityDispersionRadialPrevious=+self%velocityDispersionRadialPrevious                                                                                                                     &
               &                                +self%velocityDispersionRadialTableRadiusInterpolator%interpolate(radiusScaleFreeEffective,self%velocityDispersionRadialTable(:,jRadiusCore(iRadiusCore))) &
               &                                *                                                                                                                               hRadiusCore(iRadiusCore)
       end do
       self%velocityDispersionRadialPrevious=+self%velocityDispersionRadialPrevious &
            &                                *sqrt(                                 &
            &                                      +gravitationalConstantGalacticus &
            &                                      *basic            %mass ()       &
            &                                      /darkMatterProfile%scale()       &
            &                                      /(                               &
            &                                        -          concentration       &
            &                                        /   (1.0d0+concentration)      &
            &                                        +log(1.0d0+concentration)      &
            &                                       )                               &
            &                                     )
    end if
    finiteResolutionNFWRadialVelocityDispersion=self%velocityDispersionRadialPrevious
    return
  end function finiteResolutionNFWRadialVelocityDispersion
  
  subroutine finiteResolutionNFWVelocityDispersionRadialTabulate(self,radius,radiusCore)
    !!{
    Tabulates the mass enclosed within a given radius for finite resolution NFW density profiles.
    !!}
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (darkMatterProfileDMOFiniteResolutionNFW), intent(inout), target :: self
    double precision                                         , intent(in   )         :: radius               , radiusCore
    double precision                                         , parameter             :: radiusTiny   =1.0d-1
    type            (integrator                             ), save                  :: integrator_
    logical                                                  , save                  :: initialized  =.false.
    !$omp threadprivate(integrator_,initialized)
    logical                                                                          :: retabulate
    integer                                                                          :: iRadiusCore          , iRadius              , &
         &                                                                              i
    double precision                                                                 :: jeansIntegral        , jeansIntegralPrevious, &
         &                                                                              radiusLower          , radiusUpper          , &
         &                                                                              radiusOuter          , density

    do i=1,2
       retabulate=.false.
       if (.not.self%velocityDispersionRadialTableInitialized) then
          retabulate=.true.
       else if (                                                             &
            &    radius     < self%velocityDispersionRadialRadiusMinimum     &
            &   .or.                                                         &
            &    radius     > self%velocityDispersionRadialRadiusMaximum     &
            &   .or.                                                         &
            &    radiusCore < self%velocityDispersionRadialRadiusCoreMinimum &
            &   .or.                                                         &
            &    radiusCore > self%velocityDispersionRadialRadiusCoreMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreVelocityDispersionTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%velocityDispersionRadialRadiusMinimum       =min(0.5d0*radius    ,self%velocityDispersionRadialRadiusMinimum    )
       self%velocityDispersionRadialRadiusMaximum       =max(2.0d0*radius    ,self%velocityDispersionRadialRadiusMaximum    )
       self%velocityDispersionRadialRadiusCoreMinimum   =min(0.5d0*radiusCore,self%velocityDispersionRadialRadiusCoreMinimum)
       self%velocityDispersionRadialRadiusCoreMaximum   =max(2.0d0*radiusCore,self%velocityDispersionRadialRadiusCoreMaximum)
       self%velocityDispersionRadialTableRadiusCount    =int(log10(self%velocityDispersionRadialRadiusMaximum    /self%velocityDispersionRadialRadiusMinimum    )*dble(velocityDispersionRadialTableRadiusPointsPerDecade    ))+1
       self%velocityDispersionRadialTableRadiusCoreCount=int(log10(self%velocityDispersionRadialRadiusCoreMaximum/self%velocityDispersionRadialRadiusCoreMinimum)*dble(velocityDispersionRadialTableRadiusCorePointsPerDecade))+1
       if (allocated(self%velocityDispersionRadialTableRadius)) then
          deallocate(self%velocityDispersionRadialTableRadiusCore)
          deallocate(self%velocityDispersionRadialTableRadius    )
          deallocate(self%velocityDispersionRadialTable          )
       end if
       allocate(self%velocityDispersionRadialTableRadiusCore(                                              self%velocityDispersionRadialTableRadiusCoreCount))
       allocate(self%velocityDispersionRadialTableRadius    (self%velocityDispersionRadialTableRadiusCount                                                  ))
       allocate(self%velocityDispersionRadialTable          (self%velocityDispersionRadialTableRadiusCount,self%velocityDispersionRadialTableRadiusCoreCount))
       ! Create a range of radii and core radii.
       self%velocityDispersionRadialTableRadius    =Make_Range(self%velocityDispersionRadialRadiusMinimum    ,self%velocityDispersionRadialRadiusMaximum    ,self%velocityDispersionRadialTableRadiusCount    ,rangeType=rangeTypeLogarithmic)
       self%velocityDispersionRadialTableRadiusCore=Make_Range(self%velocityDispersionRadialRadiusCoreMinimum,self%velocityDispersionRadialRadiusCoreMaximum,self%velocityDispersionRadialTableRadiusCoreCount,rangeType=rangeTypeLogarithmic)
       ! Initialize integrator if necessary.
       if (.not.initialized) then
          integrator_=integrator(jeansEquationIntegrand,toleranceRelative=1.0d-2)
          initialized=.true.
       end if
       ! Loop over radii and alpha and populate tables.
       self_       => self
       radiusOuter =  max(10.0d0*self%velocityDispersionRadialRadiusMaximum,1000.0d0)
       do iRadiusCore=1,self%velocityDispersionRadialTableRadiusCoreCount
          iRadiusCore_         =iRadiusCore
          jeansIntegralPrevious=0.0d0
          do iRadius=self%velocityDispersionRadialTableRadiusCount,1,-1
             ! For radii that are tiny compared to the core radius the velocity dispersion become almost constant. Simply assume this to avoid floating point errors.
             if     (                                                                                                                          &
                  &   self%velocityDispersionRadialTableRadius(iRadius) < radiusTiny                                                           &
                  &  .and.                                                                                                                     &
                  &   self%velocityDispersionRadialTableRadius(iRadius) < radiusTiny*self%velocityDispersionRadialTableRadiusCore(iRadiusCore) &
                  &  .and.                                                                                                                     &
                  &   iRadius                                           < self%velocityDispersionRadialTableRadiusCount                        &
                  & ) then
                self%velocityDispersionRadialTable(iRadius,iRadiusCore)=self%velocityDispersionRadialTable(iRadius+1,iRadiusCore)
             else
                ! Find the limits for the integral.
                if (iRadius == self%velocityDispersionRadialTableRadiusCount) then
                   radiusUpper=radiusOuter
                else
                   radiusUpper=self%velocityDispersionRadialTableRadius(iRadius+1)
                end if
                radiusLower                                            =self       %velocityDispersionRadialTableRadius(                                                         iRadius     )
                density                                                =self       %densityScaleFree                   (radiusLower,self%velocityDispersionRadialTableRadiusCore(iRadiusCore))
                jeansIntegral                                          =integrator_%integrate                          (radiusLower,     radiusUpper                                         )
                self%velocityDispersionRadialTable(iRadius,iRadiusCore)=+sqrt(                         &
                     &                                                        +(                       &
                     &                                                          +jeansIntegral         &
                     &                                                          +jeansIntegralPrevious &
                     &                                                         )                       &
                     &                                                        /density                 &
                     &                                                       )
                jeansIntegralPrevious                                  =+jeansIntegralPrevious &
                     &                                                  +jeansIntegral
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%velocityDispersionRadialTableRadiusCoreInterpolator)) deallocate(self%velocityDispersionRadialTableRadiusCoreInterpolator)
       if (allocated(self%velocityDispersionRadialTableRadiusInterpolator    )) deallocate(self%velocityDispersionRadialTableRadiusInterpolator    )
       allocate(self%velocityDispersionRadialTableRadiusCoreInterpolator)
       allocate(self%velocityDispersionRadialTableRadiusInterpolator    )
       self%velocityDispersionRadialTableRadiusCoreInterpolator=interpolator(self%velocityDispersionRadialTableRadiusCore)
       self%velocityDispersionRadialTableRadiusInterpolator    =interpolator(self%velocityDispersionRadialTableRadius    )
       ! Specify that tabulation has been made.
       self%velocityDispersionRadialTableInitialized=.true.
       call self%storeVelocityDispersionTable()
    end if
    return
  end subroutine finiteResolutionNFWVelocityDispersionRadialTabulate
  
  double precision function jeansEquationIntegrand(radius)
    !!{
    Integrand for dark matter profile Jeans equation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       jeansEquationIntegrand=+self_%massEnclosedScaleFree(radius,self_%velocityDispersionRadialTableRadiusCore(iRadiusCore_))    &
            &                 *self_%densityScaleFree     (radius,self_%velocityDispersionRadialTableRadiusCore(iRadiusCore_))    &
            &                 /                            radius                                                             **2
    else
       jeansEquationIntegrand=0.0d0
    end if
    return
  end function jeansEquationIntegrand

  subroutine finiteResolutionNFWStoreVelocityDispersionTable(self)
    !!{
    Store the tabulated velocity dispersion data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'VelocityDispersion_'                            // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%velocityDispersionRadialTableRadiusCore,'radiusCore'        )
    call file%writeDataset(self%velocityDispersionRadialTableRadius    ,'radius'            )
    call file%writeDataset(self%velocityDispersionRadialTable          ,'velocityDispersion')
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine finiteResolutionNFWStoreVelocityDispersionTable

  subroutine finiteResolutionNFWRestoreVelocityDispersionTable(self)
    !!{
    Restore the tabulated velocity dispersion data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock               , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'VelocityDispersion_'                            // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%velocityDispersionRadialTableRadius)) then
          deallocate(self%velocityDispersionRadialTableRadiusCore)
          deallocate(self%velocityDispersionRadialTableRadius    )
          deallocate(self%velocityDispersionRadialTable          )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('radiusCore'        ,self%velocityDispersionRadialTableRadiusCore)
       call file%readDataset('radius'            ,self%velocityDispersionRadialTableRadius    )
       call file%readDataset('velocityDispersion',self%velocityDispersionRadialTable          )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%velocityDispersionRadialTableRadiusCount    =size(self%velocityDispersionRadialTableRadius    )
       self%velocityDispersionRadialTableRadiusCoreCount=size(self%velocityDispersionRadialTableRadiusCore)
       self%velocityDispersionRadialRadiusMinimum       =self%velocityDispersionRadialTableRadius    (                                                1)
       self%velocityDispersionRadialRadiusMaximum       =self%velocityDispersionRadialTableRadius    (self%velocityDispersionRadialTableRadiusCount    )
       self%velocityDispersionRadialRadiusCoreMinimum   =self%velocityDispersionRadialTableRadiusCore(                                                1)
       self%velocityDispersionRadialRadiusCoreMaximum   =self%velocityDispersionRadialTableRadiusCore(self%velocityDispersionRadialTableRadiusCoreCount)
       if (allocated(self%velocityDispersionRadialTableRadiusCoreInterpolator)) deallocate(self%velocityDispersionRadialTableRadiusCoreInterpolator)
       if (allocated(self%velocityDispersionRadialTableRadiusInterpolator    )) deallocate(self%velocityDispersionRadialTableRadiusInterpolator    )
       allocate(self%velocityDispersionRadialTableRadiusCoreInterpolator)
       allocate(self%velocityDispersionRadialTableRadiusInterpolator    )
       self%velocityDispersionRadialTableRadiusCoreInterpolator=interpolator(self%velocityDispersionRadialTableRadiusCore)
       self%velocityDispersionRadialTableRadiusInterpolator    =interpolator(self%velocityDispersionRadialTableRadius    )
       self%velocityDispersionRadialTableInitialized           =.true.
    end if    
    return
  end subroutine finiteResolutionNFWRestoreVelocityDispersionTable

  subroutine finiteResolutionNFWStoreDensityTable(self)
    !!{
    Store the tabulated radius-enclosing-density data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Density_'                                       // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%radiusEnclosingDensityTableRadiusCore,'radiusCore')
    call file%writeDataset(self%radiusEnclosingDensityTableDensity   ,'density'   )
    call file%writeDataset(self%radiusEnclosingDensityTable          ,'radius'    )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine finiteResolutionNFWStoreDensityTable

  subroutine finiteResolutionNFWRestoreDensityTable(self)
    !!{
    Restore the tabulated radius-enclosing-density data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Density_'                                       // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%radiusEnclosingDensityTableDensity)) then
          deallocate(self%radiusEnclosingDensityTableRadiusCore)
          deallocate(self%radiusEnclosingDensityTableDensity   )
          deallocate(self%radiusEnclosingDensityTable          )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('radiusCore',self%radiusEnclosingDensityTableRadiusCore)
       call file%readDataset('density'   ,self%radiusEnclosingDensityTableDensity   )
       call file%readDataset('radius'    ,self%radiusEnclosingDensityTable          )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%radiusEnclosingDensityTableDensityCount   =size(self%radiusEnclosingDensityTableDensity   )
       self%radiusEnclosingDensityTableRadiusCoreCount=size(self%radiusEnclosingDensityTableRadiusCore)
       self%radiusEnclosingDensityDensityMinimum      =self%radiusEnclosingDensityTableDensity   (                                              1)
       self%radiusEnclosingDensityDensityMaximum      =self%radiusEnclosingDensityTableDensity   (self%radiusEnclosingDensityTableDensityCount   )
       self%radiusEnclosingDensityRadiusCoreMinimum   =self%radiusEnclosingDensityTableRadiusCore(                                              1)
       self%radiusEnclosingDensityRadiusCoreMaximum   =self%radiusEnclosingDensityTableRadiusCore(self%radiusEnclosingDensityTableRadiusCoreCount)
       if (allocated(self%radiusEnclosingDensityTableRadiusCoreInterpolator)) deallocate(self%radiusEnclosingDensityTableRadiusCoreInterpolator)
       if (allocated(self%radiusEnclosingDensityTableDensityInterpolator   )) deallocate(self%radiusEnclosingDensityTableDensityInterpolator   )
       allocate(self%radiusEnclosingDensityTableRadiusCoreInterpolator)
       allocate(self%radiusEnclosingDensityTableDensityInterpolator   )
       self%radiusEnclosingDensityTableRadiusCoreInterpolator=interpolator(self%radiusEnclosingDensityTableRadiusCore)
       self%radiusEnclosingDensityTableDensityInterpolator   =interpolator(self%radiusEnclosingDensityTableDensity   )
       self%radiusEnclosingDensityTableInitialized           =.true.
    end if    
    return
  end subroutine finiteResolutionNFWRestoreDensityTable

  subroutine finiteResolutionNFWStoreMassTable(self)
    !!{
    Store the tabulated radius-enclosing-mass data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Mass_'                                          // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%radiusEnclosingMassTableRadiusCore,'radiusCore')
    call file%writeDataset(self%radiusEnclosingMassTableMass      ,'mass'      )
    call file%writeDataset(self%radiusEnclosingMassTable          ,'radius'    )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine finiteResolutionNFWStoreMassTable

  subroutine finiteResolutionNFWRestoreMassTable(self)
    !!{
    Restore the tabulated rdius-enclosing-mass data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Mass_'                                          // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%radiusEnclosingMassTableMass)) then
          deallocate(self%radiusEnclosingMassTableRadiusCore)
          deallocate(self%radiusEnclosingMassTableMass      )
          deallocate(self%radiusEnclosingMassTable          )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('radiusCore',self%radiusEnclosingMassTableRadiusCore)
       call file%readDataset('mass'      ,self%radiusEnclosingMassTableMass      )
       call file%readDataset('radius'    ,self%radiusEnclosingMassTable          )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%radiusEnclosingMassTableMassCount      =size(self%radiusEnclosingMassTableMass      )
       self%radiusEnclosingMassTableRadiusCoreCount=size(self%radiusEnclosingMassTableRadiusCore)
       self%radiusEnclosingMassMassMinimum         =self%radiusEnclosingMassTableMass      (                                           1)
       self%radiusEnclosingMassMassMaximum         =self%radiusEnclosingMassTableMass      (self%radiusEnclosingMassTableMassCount      )
       self%radiusEnclosingMassRadiusCoreMinimum   =self%radiusEnclosingMassTableRadiusCore(                                           1)
       self%radiusEnclosingMassRadiusCoreMaximum   =self%radiusEnclosingMassTableRadiusCore(self%radiusEnclosingMassTableRadiusCoreCount)
       if (allocated(self%radiusEnclosingMassTableRadiusCoreInterpolator)) deallocate(self%radiusEnclosingMassTableRadiusCoreInterpolator)
       if (allocated(self%radiusEnclosingMassTableMassInterpolator      )) deallocate(self%radiusEnclosingMassTableMassInterpolator      )
       allocate(self%radiusEnclosingMassTableRadiusCoreInterpolator)
       allocate(self%radiusEnclosingMassTableMassInterpolator      )
       self%radiusEnclosingMassTableRadiusCoreInterpolator=interpolator(self%radiusEnclosingMassTableRadiusCore)
       self%radiusEnclosingMassTableMassInterpolator      =interpolator(self%radiusEnclosingMassTableMass      )
       self%radiusEnclosingMassTableInitialized           =.true.
    end if    
    return
  end subroutine finiteResolutionNFWRestoreMassTable
  
  subroutine finiteResolutionNFWStoreEnergyTable(self)
    !!{
    Store the tabulated energy data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Energy_'                                        // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%energyTableRadiusCore   ,'radiusCore'   )
    call file%writeDataset(self%energyTableConcentration,'concentration')
    call file%writeDataset(self%energyTable             ,'energy'       )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine finiteResolutionNFWStoreEnergyTable

  subroutine finiteResolutionNFWRestoreEnergyTable(self)
    !!{
    Restore the tabulated rdius-enclosing-mass data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                         )                :: fileLock
    type (hdf5Object                             )                :: file
    type (varying_string                         )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Energy_'                                        // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%energyTableConcentration)) then
          deallocate(self%energyTableRadiusCore   )
          deallocate(self%energyTableConcentration)
          deallocate(self%energyTable             )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('radiusCore'   ,self%energyTableRadiusCore   )
       call file%readDataset('concentration',self%energyTableConcentration)
       call file%readDataset('energy'       ,self%energyTable             )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%energyTableConcentrationCount=size(self%energyTableConcentration      )
       self%energyTableRadiusCoreCount  =size(self%energyTableRadiusCore)
       self%energyConcentrationMinimum  =self%energyTableConcentration(                                 1)
       self%energyConcentrationMaximum  =self%energyTableConcentration(self%energyTableConcentrationCount)
       self%energyRadiusCoreMinimum     =self%energyTableRadiusCore   (                                 1)
       self%energyRadiusCoreMaximum     =self%energyTableRadiusCore   (self%energyTableRadiusCoreCount   )
       if (allocated(self%energyTableRadiusCoreInterpolator   )) deallocate(self%energyTableRadiusCoreInterpolator   )
       if (allocated(self%energyTableConcentrationInterpolator)) deallocate(self%energyTableConcentrationInterpolator)
       allocate(self%energyTableRadiusCoreInterpolator   )
       allocate(self%energyTableConcentrationInterpolator)
       self%energyTableRadiusCoreInterpolator   =interpolator(self%energyTableRadiusCore   )
       self%energyTableConcentrationInterpolator=interpolator(self%energyTableConcentration)
       self%energyTableInitialized              =.true.
    end if    
    return
  end subroutine finiteResolutionNFWRestoreEnergyTable
