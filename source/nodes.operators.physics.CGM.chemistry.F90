!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a node operator class that solves for chemical evolution in the \gls{cgm}.
  !!}

  use :: Atomic_Cross_Sections_Ionization_Photo, only : atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Rates_Ionization_Collisional   , only : atomicIonizationRateCollisionalClass
  use :: Atomic_Rates_Recombination_Radiative  , only : atomicRecombinationRateRadiativeClass
  use :: Chemical_Reaction_Rates               , only : chemicalReactionRateClass
  use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales               , only : darkMatterHaloScaleClass
  use :: Radiation_Fields                      , only : radiationFieldClass                   , crossSectionFunctionTemplate
  use :: Numerical_Constants_Physical          , only : plancksConstant                       , speedLight
  use :: Numerical_Constants_Units             , only : metersToAngstroms                     , electronVolt

  !![
  <nodeOperator name="nodeOperatorCGMChemistry">
    <description>
      A node operator class solves for chemical evolution in the \gls{cgm}. Chemical abundances are evolved according to a
      \refClass{chemicalReactionRateClass} object, with the option of the ionization state of atomic hydrogen being set to
      equilibrium values. This can be advantageous as the timescales for the reactions controlling the ionization state of atomic
      hydrogen can become extremely small, resulting in extremely slow evolution of the ODE system.

      The parameter {\normalfont \ttfamily fractionTimescaleEquilibrium} controls when the equilibrium assumption should be made. Specifically, equilibrium is assume if:
      \begin{equation}
      \tau_\mathrm{H} &lt; f_\mathrm{dyn} \tau_\mathrm{dyn},
      \end{equation}
      where $f_\mathrm{dyn}=${\normalfont \ttfamily [fractionTimescaleEquilibrium]}, $\tau_\mathrm{dyn}$ is the dynamical time in
      the halo, and
      \begin{equation}
      \tau_\mathrm{H} = \mathrm{min}\left(\tau_\alpha,\tau_\beta,\tau_\Gamma\right),
      \end{equation}
      where $\tau_\alpha=1/\alpha n$, $\tau_\beta=1/\beta n$,$\tau_\Gamma=1/\Gamma$, $n$ is the number density of hydrogen, and
      $\alpha$, $\beta$, and $\Gamma$ are the collisional ionization, radiative recombination, and photoionization rate
      coefficients for hydrogen respectively.

      If the system is judged to be in equilibrium then the neutral fraction of hydrogen is computed as:
      \begin{equation}
      x_\mathrm{H} = \frac{ \tau_\Gamma^{-1} + \tau_\alpha^{-1} +2 \tau_\beta^{-1} - \sqrt{ \tau_\Gamma^{-2} +2 \tau_\alpha^{-1} \tau_\Gamma^{-1}  + \tau_\alpha^{-2}  + 4 \tau_\beta^{-1} \tau_\Gamma^{-1} } }{ 2 \tau_\alpha^{-1} + \tau_\beta^{-1} }.
      \end{equation}
      The abundances of H, H$^+$, and e$^-$ are then fixed according to this fraction, and reaction rates for them are set to
      zero.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMChemistry
     !!{
     A node operator class that solves for chemical evolution in the \gls{cgm}.
     !!}
     private
     logical                                                                              :: chemicalsPresent                            , assumeEquilibrium
     class           (atomicIonizationRateCollisionalClass   ), pointer                   :: atomicIonizationRateCollisional_   => null()
     class           (atomicRecombinationRateRadiativeClass  ), pointer                   :: atomicRecombinationRateRadiative_  => null()
     class           (atomicCrossSectionIonizationPhotoClass ), pointer                   :: atomicCrossSectionIonizationPhoto_ => null()
     class           (chemicalReactionRateClass              ), pointer                   :: chemicalReactionRate_              => null()
     class           (darkMatterHaloScaleClass               ), pointer                   :: darkMatterHaloScale_               => null()
     class           (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_                => null()
     class           (radiationFieldClass                    ), pointer                   :: radiation_                         => null()
     logical                                                  , allocatable, dimension(:) :: maskAnalytic
     integer                                                                              :: atomicHydrogenIndex                         , atomicHydrogenCationIndex, &
          &                                                                                  electronIndex
     double precision                                                                     :: fractionTimescaleEquilibrium
   contains
     !![
     <methods>
       <method method="atomicEquilibrium" description="Determine if equilibrium should be assumed for atomic abundances."/>
       <method method="computeState"      description="Compute the state of the chemical system."                        />
     </methods>
     !!]
     final     ::                                        cgmChemistryDestructor
     procedure :: differentialEvolutionAnalytics      => cgmChemistryDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => cgmChemistryDifferentialEvolutionSolveAnalytics
     procedure :: differentialEvolution               => cgmChemistryDifferentialEvolution
     procedure :: atomicEquilibrium                   => cgmChemistryAtomicEquilibrium
     procedure :: computeState                        => cgmChemistryComputeState
  end type nodeOperatorCGMChemistry

  interface nodeOperatorCGMChemistry
     !!{
     Constructors for the \refClass{nodeOperatorCGMChemistry} node operator class.
     !!}
     module procedure cgmChemistryConstructorParameters
     module procedure cgmChemistryConstructorInternal
  end interface nodeOperatorCGMChemistry

  ! Ionization edge wavelength for HI.
  double precision                              , parameter :: energyIonizationAtomicHydrogen     =  +13.6d0
  double precision                              , parameter :: wavelengthIonizationAtomicHydrogen =  +plancksConstant                &
       &                                                                                             *speedLight                     &
       &                                                                                             *metersToAngstroms              &
       &                                                                                             /electronVolt                   &
       &                                                                                             /energyIonizationAtomicHydrogen
  ! Pointers to cross-section functions.
  procedure       (crossSectionFunctionTemplate), pointer   :: crossSectionPhotoionization_       => crossSectionPhotoionization

  ! Sub-module scope pointers used in integrations.
  class           (nodeOperatorCGMChemistry    ), pointer   :: self_                              => null()
  !$omp threadprivate(self_)
  
  ! Minimum hydrogen number density for which we perform calculations.
  double precision                              , parameter :: numberDensityHydrogenMinimum       =  1.0d-20
  
contains

  function cgmChemistryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMChemistry} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Radiation_Fields, only : radiationFieldNull
    implicit none
    type            (nodeOperatorCGMChemistry              )                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (atomicIonizationRateCollisionalClass  ), pointer       :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateRadiativeClass ), pointer       :: atomicRecombinationRateRadiative_
    class           (atomicCrossSectionIonizationPhotoClass), pointer       :: atomicCrossSectionIonizationPhoto_
    class           (chemicalReactionRateClass             ), pointer       :: chemicalReactionRate_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (radiationFieldClass                   ), pointer       :: radiation_
    double precision                                                        :: fractionTimescaleEquilibrium

    !![
    <inputParameter>
      <name>fractionTimescaleEquilibrium</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The fraction of the halo dynamical time which, if atomic chemistry timescales are smaller than, switch to an equilibrium calculation of atomic abundances.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="chemicalReactionRate"              name="chemicalReactionRate_"              source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"               name="darkMatterHaloScale_"               source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                name="cosmologyFunctions_"                source="parameters"/>
    <objectBuilder class="atomicIonizationRateCollisional"   name="atomicIonizationRateCollisional_"   source="parameters"/>
    <objectBuilder class="atomicRecombinationRateRadiative"  name="atomicRecombinationRateRadiative_"  source="parameters"/>
    <objectBuilder class="atomicCrossSectionIonizationPhoto" name="atomicCrossSectionIonizationPhoto_" source="parameters"/>
    !!]
    radiation_ => null()
    if (parameters%isPresent('radiationFieldIntergalacticBackground',searchInParents=.true.)) then
       !![
       <objectBuilder class="radiationField" name="radiation_" parameterName="radiationFieldIntergalacticBackground" source="parameters"/>
       !!]
    else
       allocate(radiationFieldNull :: radiation_)
       select type (radiation_)
       type is (radiationFieldNull)
          !![
	  <referenceConstruct object="radiation_" constructor="radiationFieldNull()"/>
          !!]
       end select
    end if    
    self=nodeOperatorCGMChemistry(fractionTimescaleEquilibrium,atomicIonizationRateCollisional_,atomicRecombinationRateRadiative_,atomicCrossSectionIonizationPhoto_,chemicalReactionRate_,darkMatterHaloScale_,cosmologyFunctions_,radiation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="chemicalReactionRate_"             />
    <objectDestructor name="darkMatterHaloScale_"              />
    <objectDestructor name="cosmologyFunctions_"               />
    <objectDestructor name="radiation_"                        />
    <objectDestructor name="atomicIonizationRateCollisional_"  />
    <objectDestructor name="atomicRecombinationRateRadiative_" />
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"/>
    !!]
    return
  end function cgmChemistryConstructorParameters

  function cgmChemistryConstructorInternal(fractionTimescaleEquilibrium,atomicIonizationRateCollisional_,atomicRecombinationRateRadiative_,atomicCrossSectionIonizationPhoto_,chemicalReactionRate_,darkMatterHaloScale_,cosmologyFunctions_,radiation_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMChemistry} node operator class.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index, Chemicals_Property_Count
    implicit none
    type            (nodeOperatorCGMChemistry              )                         :: self
    double precision                                        , intent(in   )          :: fractionTimescaleEquilibrium
    class           (atomicIonizationRateCollisionalClass  ), intent(in   ), target  :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateRadiativeClass ), intent(in   ), target  :: atomicRecombinationRateRadiative_
    class           (atomicCrossSectionIonizationPhotoClass), intent(in   ), target  :: atomicCrossSectionIonizationPhoto_
    class           (chemicalReactionRateClass             ), intent(in   ), target  :: chemicalReactionRate_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target  :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass               ), intent(in   ), target  :: cosmologyFunctions_
    class           (radiationFieldClass                   ), intent(in   ), pointer :: radiation_
    !![
    <constructorAssign variables="fractionTimescaleEquilibrium, *atomicIonizationRateCollisional_, *atomicRecombinationRateRadiative_, *atomicCrossSectionIonizationPhoto_, *chemicalReactionRate_, *darkMatterHaloScale_, *cosmologyFunctions_, *radiation_"/>
    !!]

    ! Determine if chemicals are being solved for.
    self%chemicalsPresent=Chemicals_Property_Count() > 0
    if (.not.self%chemicalsPresent) return
    ! Get indices of chemicals needed for equilibrium hydrogen anion calculation.
    self%atomicHydrogenIndex      =Chemicals_Index("AtomicHydrogen"      )
    self%atomicHydrogenCationIndex=Chemicals_Index("AtomicHydrogenCation")
    self%electronIndex            =Chemicals_Index("Electron"            )
    ! Create a mask for properties that may be solved for analytically.
    allocate(self%maskAnalytic(Chemicals_Property_Count()))
    self%maskAnalytic                                =.false.
    self%maskAnalytic(self%atomicHydrogenIndex      )=.true.
    self%maskAnalytic(self%atomicHydrogenCationIndex)=.true.
    self%maskAnalytic(self%electronIndex            )=.true.
    return
  end function cgmChemistryConstructorInternal
  
  subroutine cgmChemistryDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMChemistry} node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMChemistry), intent(inout) :: self

    !![
    <objectDestructor name="self%chemicalReactionRate_"             />
    <objectDestructor name="self%darkMatterHaloScale_"              />
    <objectDestructor name="self%cosmologyFunctions_"               />
    <objectDestructor name="self%radiation_"                        />
    <objectDestructor name="self%atomicIonizationRateCollisional_"  />
    <objectDestructor name="self%atomicRecombinationRateRadiative_" />
    <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"/>
    !!]
    return
  end subroutine cgmChemistryDestructor

  subroutine cgmChemistryDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances 
    implicit none
    class           (nodeOperatorCGMChemistry), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                                          :: timescaleCollisionalIonization, timescaleRadiativeRecombination, &
         &                                                       timescalePhotoionization      , timescaleDynamical             , &
         &                                                       timescaleMinimum 

    ! Return instantly if no chemicals are tracked.
    if (.not.self%chemicalsPresent) return
    ! Compute atomic timescales, and the halo dynamical timescale. Test if timescales are short enough that we can assume
    ! equilibrium behavior.
    call self%atomicEquilibrium(node,timescaleCollisionalIonization,timescaleRadiativeRecombination,timescalePhotoionization)
    timescaleMinimum      =min(min(timescaleCollisionalIonization,timescaleRadiativeRecombination),timescalePhotoionization)
    timescaleDynamical    =self%darkMatterHaloScale_%timescaleDynamical(node)
    self%assumeEquilibrium=timescaleMinimum < self%fractionTimescaleEquilibrium*timescaleDynamical
    return
  end subroutine cgmChemistryDifferentialEvolutionAnalytics

  subroutine cgmChemistryDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve \gls{cgm} chemistry.
    !!}
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo
    use :: Chemical_Abundances_Structure, only : chemicalAbundances 
    use :: Numerical_Constants_Atomic   , only : atomicMassUnit
    use :: Numerical_Constants_Physical , only : electronMass
    implicit none
    class           (nodeOperatorCGMChemistry), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: time
    class           (nodeComponentHotHalo    ), pointer       :: hotHalo
    type            (chemicalAbundances      ), save          :: chemicalMasses
    !$omp threadprivate(chemicalMasses)
    double precision                                          :: timescaleCollisionalIonization, timescaleRadiativeRecombination, &
         &                                                       timescalePhotoionization      , fractionAtomicHydrogen         , &
         &                                                       rateCollisionalIonization     , rateRadiativeRecombination     , &
         &                                                       ratePhotoionization           , massHydrogen

    if (self%chemicalsPresent.and.self%assumeEquilibrium) then
       ! Compute densities of H, H⁺, and e⁻ from photoionization equilibrium by solving the quadratic equation for the equilibrium
       ! neutral hydrogen fraction.
       call self%atomicEquilibrium(node,timescaleCollisionalIonization,timescaleRadiativeRecombination,timescalePhotoionization)
       rateCollisionalIonization  =  +1.0d0/timescaleCollisionalIonization
       rateRadiativeRecombination =  +1.0d0/timescaleRadiativeRecombination
       ratePhotoionization        =  +1.0d0/timescalePhotoionization
       fractionAtomicHydrogen     =  +(                                                                 &
            &                                +                                   ratePhotoionization    &
            &                                +      rateCollisionalIonization                           &
            &                                +2.0d0*rateRadiativeRecombination                          &
            &                          -sqrt(                                                           &
            &                                +                                   ratePhotoionization**2 &
            &                                +2.0d0*rateCollisionalIonization   *ratePhotoionization    &
            &                                +      rateCollisionalIonization**2                        &
            &                                +4.0d0*rateRadiativeRecombination  *ratePhotoionization    &
            &                               )                                                           &
            &                          )                                                                &
            &                         /       2.0d0                                                     &
            &                         /(                                                                &
            &                           +           rateCollisionalIonization                           &
            &                           +           rateRadiativeRecombination                          &
            &                          )
       hotHalo                    =>   node          %hotHalo  (                              )
       chemicalMasses             =    hotHalo       %chemicals(                              )
       massHydrogen               =   +chemicalMasses%abundance(self%atomicHydrogenIndex      ) &
            &                         +chemicalMasses%abundance(self%atomicHydrogenCationIndex)
       call chemicalMasses%abundanceSet(self%atomicHydrogenIndex      ,        +fractionAtomicHydrogen *massHydrogen                            )
       call chemicalMasses%abundanceSet(self%atomicHydrogenCationIndex,+(+1.0d0-fractionAtomicHydrogen)*massHydrogen                            )
       call chemicalMasses%abundanceSet(self%electronIndex            ,+(+1.0d0-fractionAtomicHydrogen)*massHydrogen*electronMass/atomicMassUnit)
       call hotHalo%chemicalsSet(chemicalMasses)
    end if
    return
  end subroutine cgmChemistryDifferentialEvolutionSolveAnalytics

  subroutine cgmChemistryDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform differential evolution.
    !!}
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances 
    use :: Galacticus_Nodes                 , only : nodeComponentHotHalo
    use :: Numerical_Constants_Astronomical , only : gigaYear             , megaParsec
    use :: Numerical_Constants_Prefixes     , only : centi
    use :: Numerical_Constants_Math         , only : Pi
    use :: Mass_Distributions               , only : massDistributionClass
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (nodeOperatorCGMChemistry), intent(inout), target    :: self
    type            (treeNode                ), intent(inout), target    :: node
    logical                                   , intent(inout)            :: interrupt
    procedure       (interruptTask           ), intent(inout), pointer   :: functionInterrupt
    integer                                   , intent(in   )            :: propertyType
    class           (nodeComponentHotHalo    )               , pointer   :: hotHalo
    class           (massDistributionClass   )               , pointer   :: massDistribution_
    double precision                                         , parameter :: massHotHaloTiny        =1.0d-6
    type            (chemicalAbundances      ), save                     :: chemicalDensitiesRates        , chemicalMassesRates, &
         &                                                                  chemicalDensities
    !$omp threadprivate(chemicalDensities,chemicalDensitiesRates,chemicalMassesRates)
    double precision                                                     :: temperature                   , lengthColumn       , &
         &                                                                  massToDensityConversion       , radiusOuter        , &
         &                                                                  factorBoostColumn             , massHotHalo        , &
         &                                                                  factorClumping

    ! Return instantly if no chemicals are tracked.
    if (.not.self%chemicalsPresent) return
    ! Compute the state of the chemical system.
    call self%computeState(node,temperature,chemicalDensities,massToDensityConversion)
    ! Return if unphysical.
    if (massToDensityConversion <= 0.0d0) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Compute the column length through the halo (in cm).
    massHotHalo=hotHalo%mass()
    if (massHotHalo > massHotHaloTiny) then
       massDistribution_ =>  node             %massDistribution     (componentTypeHotHalo,massTypeGaseous)
       radiusOuter       =  +hotHalo          %outerRadius          (                                    )
       factorBoostColumn =  +massDistribution_%densityRadialMoment  (0.0d0,radiusOuter                   ) &
            &               *4.0d0                                                                         &
            &               *Pi                                                                            &
            &               *radiusOuter  **2                                                              &
            &               /3.0d0                                                                         &
            &               /massHotHalo
       lengthColumn      =  +radiusOuter                                                                   &
            &               *megaParsec                                                                    &
            &               /centi
       factorClumping    =  +massDistribution_%densitySquareIntegral(0.0d0,radiusOuter                   ) &
            &               *4.0d0                                                                         &
            &               /3.0d0                                                                         &
            &               *Pi                                                                            &
            &               *radiusOuter**3                                                                &
            &               /massHotHalo**2
       !![
       <objectDestructor name="massDistribution_"/>
       !!]          
     else
       lengthColumn     =+0.0d0
       factorClumping   =+1.0d0
    end if
    ! Compute the chemical reaction rates.
    call self%chemicalReactionRate_%rates(lengthColumn,temperature,chemicalDensities,factorClumping,self%radiation_,chemicalDensitiesRates,node)
    ! Convert to mass change rates.
    call chemicalDensitiesRates%numberToMass(chemicalMassesRates)
    call chemicalMassesRates   %scale       (                         &
         &                                   +gigaYear                &
         &                                   /massToDensityConversion &
         &                                  )
    ! Zero rates of analytically-solved properties.
    if (self%assumeEquilibrium) then
       call chemicalMassesRates%abundanceSet(self%atomicHydrogenIndex      ,0.0d0)
       call chemicalMassesRates%abundanceSet(self%atomicHydrogenCationIndex,0.0d0)
       call chemicalMassesRates%abundanceSet(self%electronIndex            ,0.0d0)
    end if
    ! Adjust rates appropriately.
    call hotHalo%chemicalsRate(chemicalMassesRates,interrupt,functionInterrupt)
    return
  end subroutine cgmChemistryDifferentialEvolution

  subroutine cgmChemistryComputeState(self,node,temperature,chemicalDensities,massToDensityConversion)
    !!{
    Compute the state of the chemical system.
    !!}
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances 
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                   , nodeComponentHotHalo
    implicit none
    class           (nodeOperatorCGMChemistry), intent(inout)           :: self
    type            (treeNode                ), intent(inout)           :: node
    double precision                          , intent(  out)           :: temperature
    type            (chemicalAbundances      ), intent(inout)           :: chemicalDensities
    double precision                          , intent(  out), optional :: massToDensityConversion
    class           (nodeComponentBasic      ), pointer                 :: basic
    class           (nodeComponentHotHalo    ), pointer                 :: hotHalo
    type            (chemicalAbundances      ), save                    :: chemicalMasses
    double precision                                                    :: massToDensityConversion_, massChemicals, &
         &          radiusOuter
    !$omp threadprivate(chemicalMasses)

    ! Get required components.
    basic   => node%basic  ()
    hotHalo => node%hotHalo()
    ! Get the temperature of the CGM.
    temperature=self%darkMatterHaloScale_%temperatureVirial(node)
    ! Set the radiation background.
    call self%radiation_%timeSet(basic%time())
    ! Get the masses of chemicals.
    chemicalMasses=hotHalo%chemicals()
    ! Truncate masses to zero to avoid unphysical behavior.
    call chemicalMasses%enforcePositive()
    massChemicals=chemicalMasses%sumOver()
    if     (                                                                        &
         &             massChemicals >         0.0d0                                &
         &  .and.                                                                   &
         &             massChemicals >         hotHalo%mass()                       &
         &  .and.                                                                   &
         &   -exponent(massChemicals)+exponent(hotHalo%mass()) < maxExponent(0.0d0) &
         & ) call chemicalMasses%scale(                                             &
         &                             +hotHalo%mass          ()                    &
         &                             /        massChemicals                       &
         &                            )
    ! Scale all chemical masses by their mass in atomic mass units to get a number density.
    call chemicalMasses%massToNumber(chemicalDensities)
    ! Compute factor converting mass of chemicals in (M☉/mᵤ) to number density in cm⁻³.
    radiusOuter=hotHalo%outerRadius()
    if (radiusOuter > 0.0d0) then
       massToDensityConversion_=Chemicals_Mass_To_Density_Conversion(radiusOuter)
    else
       massToDensityConversion_=0.0d0
    end if
    if (present(massToDensityConversion)) massToDensityConversion=massToDensityConversion_
    ! Convert to number density.
    call chemicalDensities%scale(massToDensityConversion_)
    return
  end subroutine cgmChemistryComputeState

  subroutine cgmChemistryAtomicEquilibrium(self,node,timescaleCollisionalIonization,timescaleRadiativeRecombination,timescalePhotoionization)
    !!{
    Determine if equilibrium should be assumed for atomic abundances.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class           (nodeOperatorCGMChemistry), intent(inout), target :: self
    type            (treeNode                ), intent(inout)         :: node
    double precision                          , intent(  out)         :: timescaleCollisionalIonization       , timescaleRadiativeRecombination     , &
         &                                                               timescalePhotoionization
    type            (chemicalAbundances      ), save                  :: chemicalDensities
    double precision                                                  :: temperature                          , rateCoefficientCollisionalIonization, &
         &                                                               rateCoefficientRadiativeRecombination, rateCoefficientPhotoionization      , &
         &                                                               numberDensityHydrogen
    !$omp threadprivate(chemicalDensities)

    call self%computeState(node,temperature,chemicalDensities)
    self_                                 =>  self
    rateCoefficientCollisionalIonization  =   self             %atomicIonizationRateCollisional_ %rate                     (1,1,temperature)
    rateCoefficientRadiativeRecombination =   self             %atomicRecombinationRateRadiative_%rate                     (1,1,temperature)
    rateCoefficientPhotoionization        =   self             %radiation_                       %integrateOverCrossSection([0.0d0,wavelengthIonizationAtomicHydrogen],crossSectionPhotoionization_,node)
    numberDensityHydrogen                 =  +chemicalDensities                                  %abundance                (self%atomicHydrogenIndex      ) &
         &                                   +chemicalDensities                                  %abundance                (self%atomicHydrogenCationIndex)
    if (numberDensityHydrogen > numberDensityHydrogenMinimum) then
       if (rateCoefficientCollisionalIonization  > 0.0d0) then
          timescaleCollisionalIonization =+1.0d0/rateCoefficientCollisionalIonization /numberDensityHydrogen/gigaYear
       else
          timescaleCollisionalIonization =+huge(0.0d0)
       end if
       if (rateCoefficientRadiativeRecombination > 0.0d0) then
          timescaleRadiativeRecombination=+1.0d0/rateCoefficientRadiativeRecombination/numberDensityHydrogen/gigaYear
       else
          timescaleRadiativeRecombination=+huge(0.0d0)
       end if
       if (rateCoefficientPhotoionization        > 0.0d0) then
          timescalePhotoionization       =+1.0d0/rateCoefficientPhotoionization                             /gigaYear
       else
          timescalePhotoionization       =+huge(0.0d0)
       end if
    else
       timescaleCollisionalIonization    =+huge(0.0d0)
       timescaleRadiativeRecombination   =+huge(0.0d0)
       timescalePhotoionization          =+huge(0.0d0)
    end if
    return
  end subroutine cgmChemistryAtomicEquilibrium

  double precision function crossSectionPhotoionization(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$ as given by
    \cite{abel_modeling_1997}.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    use :: Tables                    , only : table1DLogarithmicLinear
    use :: Table_Labels              , only : extrapolationTypeZero
    implicit none
    double precision                          , intent(in   ) :: wavelength
    double precision                          , parameter     :: wavelengthFactor=1.0d-3
    integer                                   , parameter     :: wavelengthCount =100
    type            (table1DLogarithmicLinear), save          :: interpolator_
    logical                                   , save          :: initialized     =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                          :: crossSection
    integer                                                   :: i

    if (.not.initialized) then
       call interpolator_%create(wavelengthFactor*lymanSeriesLimitWavelengthHydrogen_atomic,lymanSeriesLimitWavelengthHydrogen_atomic,wavelengthCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,wavelengthCount
          ! Use the hydrogen photoionization cross section method.
          crossSection=self_%atomicCrossSectionIonizationPhoto_%crossSection(1,1,1,interpolator_%x(i))
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Evaluate the cross-section.
    crossSectionPhotoionization=interpolator_%interpolate(wavelength)
    return
  end function crossSectionPhotoionization
