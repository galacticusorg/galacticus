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

  use :: Atomic_Cross_Sections_Ionization_Photo      , only : atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Ionization_Potentials                , only : atomicIonizationPotentialClass
  use :: Atomic_Radiation_Gaunt_Factors              , only : gauntFactorClass
  use :: Atomic_Rates_Excitation_Collisional         , only : atomicExcitationRateCollisionalClass
  use :: Atomic_Rates_Ionization_Collisional         , only : atomicIonizationRateCollisionalClass
  use :: Atomic_Rates_Recombination_Dielectronic     , only : atomicRecombinationRateDielectronicClass
  use :: Atomic_Rates_Recombination_Radiative        , only : atomicRecombinationRateRadiativeClass       , enumerationRecombinationCaseType
  use :: Atomic_Rates_Recombination_Radiative_Cooling, only : atomicRecombinationRateRadiativeCoolingClass
  use :: Mass_Distributions                          , only : massDistributionClass
  use :: Root_Finder                                 , only : rootFinder

  !![
  <radiativeTransferMatter name="radiativeTransferMatterAtomic">
   <description>A task which performs radiative transfer.</description>
  </radiativeTransferMatter>
  !!]
  type, extends(radiativeTransferMatterClass) :: radiativeTransferMatterAtomic
     !!{
     Implementation of a radiative transfer matter class for atomic matter.
     !!}
     private
     class           (massDistributionClass                       ), pointer                     :: massDistribution_                        => null()
     class           (atomicCrossSectionIonizationPhotoClass      ), pointer                     :: atomicCrossSectionIonizationPhoto_       => null()
     class           (atomicRecombinationRateRadiativeClass       ), pointer                     :: atomicRecombinationRateRadiative_        => null()
     class           (atomicRecombinationRateRadiativeCoolingClass), pointer                     :: atomicRecombinationRateRadiativeCooling_ => null()
     class           (atomicIonizationRateCollisionalClass        ), pointer                     :: atomicIonizationRateCollisional_         => null()
     class           (atomicRecombinationRateDielectronicClass    ), pointer                     :: atomicRecombinationRateDielectronic_     => null()
     class           (atomicIonizationPotentialClass              ), pointer                     :: atomicIonizationPotential_               => null()
     class           (atomicExcitationRateCollisionalClass        ), pointer                     :: atomicExcitationRateCollisional_         => null()
     class           (gauntFactorClass                            ), pointer                     :: gauntFactor_                             => null()
     type            (rootFinder                                  )                              :: finder
     integer                                                                                     :: indexAbundancePattern                             , iterationAverageCount       , &
          &                                                                                         countElements                                     , indexHydrogen
     integer         (c_size_t                                    )                              :: countOutputs_                                     , countStatesConvergence
     logical                                                                                     :: outputRates                                       , outputAbsorptionCoefficients
     double precision                                                                            :: metallicity                                       , temperatureMinimum          , &
          &                                                                                         convergencePercentile                             , wavelengthPrevious
     double precision                                              , allocatable, dimension(:  ) :: numberDensityMassDensityRatio                     , elementAtomicMasses
     type            (varying_string                              )                              :: abundancePattern
     character       (len=2                                       ), allocatable, dimension(:  ) :: elements
     integer                                                       , allocatable, dimension(:  ) :: elementAtomicNumbers
     double precision                                              , allocatable, dimension(:,:) :: crossSectionPhotoIonizationPrevious
   contains
     !![
     <methods>
       <method description="Return the total rate of recombinations (in units of s$^{-1}$)." method="recombinationRateHydrogen" />
       <method description="Return the total rate of recombinations (in units of s$^{-1}$)." method="absorptionCoefficientSpecies" />
       <method description="Update the ionization state hsitory in a properties object." method="historyUpdate" />
       <method description="Compute the total photoionization cross section for the given element and ionization state." method="crossSectionPhotoIonization" />
     </methods>
     !!]
     final     ::                                 atomicDestructor
     procedure :: propertyClass                => atomicPropertyClass
     procedure :: populateDomain               => atomicPopulateDomain
     procedure :: reset                        => atomicReset
     procedure :: absorptionCoefficient        => atomicAbsorptionCoefficient
     procedure :: absorptionCoefficientSpecies => atomicAbsorptionCoefficientSpecies
     procedure :: accumulatePhotonPacket       => atomicAccumulatePhotonPacket
     procedure :: interactWithPhotonPacket     => atomicInteractWithPhotonPacket
     procedure :: stateSolve                   => atomicStateSolve
     procedure :: convergenceMeasure           => atomicConvergenceMeasure
     procedure :: outputProperty               => atomicOutputProperty
     procedure :: countOutputs                 => atomicCountOutputs
     procedure :: outputName                   => atomicOutputName
     procedure :: recombinationRateHydrogen    => atomicRecombinationRateHydrogen
     procedure :: historyUpdate                => atomicHistoryUpdate
     procedure :: crossSectionPhotoIonization  => atomicCrossSectionPhotoionization
#ifdef USEMPI
     procedure :: accumulationReduction        => atomicAccumulationReduction
     procedure :: broadcastDomain              => atomicBroadcastDomain
     procedure :: broadcastState               => atomicBroadcastState
#endif
  end type radiativeTransferMatterAtomic
  
  interface radiativeTransferMatterAtomic
     !!{
     Constructors for the \refClass{radiativeTransferMatterAtomic} radiative transfer matter class.
     !!}
     module procedure atomicConstructorParameters
     module procedure atomicConstructorInternal
  end interface radiativeTransferMatterAtomic

  type, public :: element
     !!{
     Type used to store elemental states.
     !!}
     double precision                              :: densityNumber
     double precision, dimension(:,:), allocatable :: photoIonizationRateHistory    , photoHeatingRateHistory , &
          &                                           ionizationStateFractionHistory
     double precision, dimension(:  ), allocatable :: photoIonizationRate           , photoHeatingRate        , &
          &                                           photoIonizationRatePrevious   , photoHeatingRatePrevious, &
          &                                           ionizationStateFraction
  end type element
  
  type, extends(radiativeTransferPropertiesMatter), public :: radiativeTransferPropertiesMatterAtomic
     !!{
     Radiative transfer matter properties class for atomic matter.
     !!}
     integer                                              :: iterationCount
     double precision                                     :: volume        , temperature
     type            (element), dimension(:), allocatable :: elements
  end type radiativeTransferPropertiesMatterAtomic

  ! Module-scope objects needed for solving for thermal equilibrium.
  class           (radiativeTransferMatterAtomic          )                , pointer     :: self_
  type            (element                                ), dimension( : ), allocatable :: elementsPhotoRate
  type            (radiativeTransferPropertiesMatterAtomic)                , pointer     :: properties_
  double precision                                                                       :: densityNumberElectrons
  type            (enumerationRecombinationCaseType      )                               :: recombinationCase
  !$omp threadprivate(self_,elementsPhotoRate,properties_,densityNumberElectrons,recombinationCase)

  ! Tolerance parameters.
  double precision                                         , parameter                   :: ionizationStateFractionToleranceAbsolute=1.0d-12, ionizationStateFractionToleranceRelative=1.0d-2
  double precision                                         , parameter                   :: temperatureToleranceRelative            =1.0d-03, temperatureToleranceAbsolute            =1.0d+0
  double precision                                         , parameter                   :: temperatureMaximum                      =1.0d+10

contains

  function atomicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferMatterAtomic} radiative transfer matter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters                , only : inputParameter  , inputParameters
    use :: ISO_Varying_String              , only : var_str
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (radiativeTransferMatterAtomic               )                             :: self
    type            (inputParameters                             ), intent(inout)              :: parameters
    class           (massDistributionClass                       ), pointer                    :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass      ), pointer                    :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass       ), pointer                    :: atomicRecombinationRateRadiative_
    class           (atomicRecombinationRateRadiativeCoolingClass), pointer                    :: atomicRecombinationRateRadiativeCooling_
    class           (atomicIonizationRateCollisionalClass        ), pointer                    :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass    ), pointer                    :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass              ), pointer                    :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass        ), pointer                    :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                            ), pointer                    :: gauntFactor_
    character       (len=2                                       ), dimension(:) , allocatable :: elements
    integer                                                                                    :: iterationAverageCount    
    double precision                                                                           :: temperatureMinimum                      , metallicity                 , &
         &                                                                                        convergencePercentile
    type            (varying_string                              )                             :: abundancePattern
    logical                                                                                    :: outputRates                             , outputAbsorptionCoefficients
    
    !![
    <inputParameter>
      <name>iterationAverageCount</name>
      <defaultValue>5</defaultValue>
      <description>The number of iterations over which to average the photoionization rate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>temperatureMinimum</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The minimum temperature that matter is allowed to reach in the case of zero photoheating.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>abundancePattern</name>
      <defaultValue>var_str('solar')</defaultValue>
      <description>The abundance pattern to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>metallicity</name>
      <defaultValue>metallicitySolar</defaultValue>
      <description>The metallicity to use.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('elements')) then
       allocate(elements(parameters%count('elements')))
    else
       allocate(elements(1                           ))
    end if
    !![
    <inputParameter>
      <name>elements</name>
      <defaultValue>['H']</defaultValue>
      <description>The elements to include.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputRates</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, output photoionization and heating rates.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputAbsorptionCoefficients</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, output absorption coefficients of each species (at the hydrogen Lyman continuum edge).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>convergencePercentile</name>
      <defaultValue>0.90d0</defaultValue>
      <description>The percentile used in the convergence criterion.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution"                        name="massDistribution_"                        source="parameters"/>
    <objectBuilder class="atomicCrossSectionIonizationPhoto"       name="atomicCrossSectionIonizationPhoto_"       source="parameters"/>
    <objectBuilder class="atomicRecombinationRateRadiative"        name="atomicRecombinationRateRadiative_"        source="parameters"/>
    <objectBuilder class="atomicRecombinationRateRadiativeCooling" name="atomicRecombinationRateRadiativeCooling_" source="parameters"/>
    <objectBuilder class="atomicIonizationRateCollisional"         name="atomicIonizationRateCollisional_"         source="parameters"/>
    <objectBuilder class="atomicRecombinationRateDielectronic"     name="atomicRecombinationRateDielectronic_"     source="parameters"/>
    <objectBuilder class="atomicIonizationPotential"               name="atomicIonizationPotential_"               source="parameters"/>
    <objectBuilder class="atomicExcitationRateCollisional"         name="atomicExcitationRateCollisional_"         source="parameters"/>
    <objectBuilder class="gauntFactor"                             name="gauntFactor_"                             source="parameters"/>
    !!]
    self=radiativeTransferMatterAtomic(abundancePattern,metallicity,elements,iterationAverageCount,temperatureMinimum,outputRates,outputAbsorptionCoefficients,convergencePercentile,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicRecombinationRateRadiativeCooling_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"                       />
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"      />
    <objectDestructor name="atomicRecombinationRateRadiative_"       />
    <objectDestructor name="atomicRecombinationRateRadiativeCooling_"/>
    <objectDestructor name="atomicIonizationRateCollisional_"        />
    <objectDestructor name="atomicRecombinationRateDielectronic_"    />
    <objectDestructor name="atomicIonizationPotential_"              />
    <objectDestructor name="atomicExcitationRateCollisional_"        />
    <objectDestructor name="gauntFactor_"                            />
    !!]
    return
  end function atomicConstructorParameters

  function atomicConstructorInternal(abundancePattern,metallicity,elements,iterationAverageCount,temperatureMinimum,outputRates,outputAbsorptionCoefficients,convergencePercentile,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicRecombinationRateRadiativeCooling_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferMatterAtomic} radiative transfer matter class.
    !!}
    use :: Abundances_Structure            , only : Abundances_Index_From_Name, abundances                   , adjustElementsReset          , metallicityTypeLinearByMassSolar
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup  , Atomic_Abundance             , Atomic_Mass                  , Atomic_Number
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : char
    use :: Numerical_Constants_Astronomical, only : massSolar                 , megaParsec                   , metallicitySolar
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Root_Finder                     , only : rangeExpandMultiplicative , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: String_Handling                 , only : String_Count_Words        , String_Split_Words
    implicit none
    type            (radiativeTransferMatterAtomic               )                              :: self
    integer                                                       , intent(in   )               :: iterationAverageCount
    double precision                                              , intent(in   )               :: temperatureMinimum                      , metallicity                 , &
         &                                                                                         convergencePercentile
    type            (varying_string                              ), intent(in   )               :: abundancePattern
    logical                                                       , intent(in   )               :: outputRates                             , outputAbsorptionCoefficients
    character       (len=2                                       ), intent(in   ), dimension(:) :: elements
    class           (massDistributionClass                       ), intent(in   ), target       :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass      ), intent(in   ), target       :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass       ), intent(in   ), target       :: atomicRecombinationRateRadiative_
    class           (atomicRecombinationRateRadiativeCoolingClass), intent(in   ), target       :: atomicRecombinationRateRadiativeCooling_
    class           (atomicIonizationRateCollisionalClass        ), intent(in   ), target       :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass    ), intent(in   ), target       :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass              ), intent(in   ), target       :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass        ), intent(in   ), target       :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                            ), intent(in   ), target       :: gauntFactor_
    double precision                                              , allocatable  , dimension(:) :: abundancesRelative
    double precision                                                                            :: numberDensityMassDensityRatioHydrogen    , numberDensityMassDensityRatioHelium
    type            (abundances                                  )                              :: abundances_
    integer                                                                                     :: i                                        , countIonizationStates              , &
         &                                                                                         indexElement
    !![
    <constructorAssign variables="abundancePattern, metallicity, elements, iterationAverageCount, temperatureMinimum, outputRates, outputAbsorptionCoefficients, convergencePercentile, *massDistribution_, *atomicCrossSectionIonizationPhoto_, *atomicRecombinationRateRadiative_, *atomicRecombinationRateRadiativeCooling_, *atomicIonizationRateCollisional_, *atomicRecombinationRateDielectronic_, *atomicIonizationPotential_, *atomicExcitationRateCollisional_, *gauntFactor_"/>
    !!]

    ! Initialize count of outputs. (Just 1, for temperature.)
    self%countOutputs_=1_c_size_t
    ! Extract elements to compute.
    self%countElements=size(elements)
    self%indexHydrogen=-1
    allocate(self%elementAtomicNumbers         (self%countElements))
    allocate(self%elementAtomicMasses          (self%countElements))
    allocate(self%numberDensityMassDensityRatio(self%countElements))
    ! Get an abundance pattern and compute conversion factors from total gas-phase mass density to atomic number densities (in
    ! cm⁻³) for each element.    
    self%indexAbundancePattern= Abundance_Pattern_Lookup(abundanceName=char(abundancePattern))
    call abundances_%metallicitySet(metallicity,metallicityType=metallicityTypeLinearByMassSolar,adjustElements=adjustElementsReset,abundanceIndex=self%indexAbundancePattern)
    numberDensityMassDensityRatioHydrogen=+abundances_   %hydrogenMassFraction(               )    &
         &                                /Atomic_Mass                        (shortLabel='H' )    &
         &                                /atomicMassUnit                                          &
         &                                *massSolar                                               &
         &                                /megaParsec                                          **3 &
         &                                *centi                                               **3
    numberDensityMassDensityRatioHelium  =+abundances_   %heliumMassFraction  (               )    &
         &                                /Atomic_Mass                        (shortLabel='He')    &
         &                                /atomicMassUnit                                          &
         &                                *massSolar                                               &
         &                                /megaParsec                                          **3 &
         &                                *centi                                               **3
    allocate(abundancesRelative(abundances_%serializeCount()))
    call abundances_%serialize(abundancesRelative)
    countIonizationStates=0
    do i=1,self%countElements
       self%elementAtomicNumbers(i)=Atomic_Number(shortLabel=trim(self%elements(i)))
       self%elementAtomicMasses (i)=Atomic_Mass  (shortLabel=trim(self%elements(i)))
       self%countOutputs_=self%countOutputs_+self%elementAtomicNumbers(i)+1_c_size_t
       if (outputRates                 ) self%countOutputs_=self%countOutputs_+2*self%elementAtomicNumbers(i)
       if (outputAbsorptionCoefficients) self%countOutputs_=self%countOutputs_+  self%elementAtomicNumbers(i)
       select case (self%elementAtomicNumbers(i))
       case (1)     ! Hydrogen
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHydrogen
          self%indexHydrogen                   = i
       case (2)     ! Helium
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHelium
       case default ! Metals
          indexElement=Abundances_Index_From_Name(trim(self%elements(i)))
          if (indexElement <= 0) call Error_Report('unable to find element "'//trim(self%elements(i))//'"'//{introspection:location})
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHydrogen                          &
               &                                *abundancesRelative                   (           indexElement) &
               &                                *Atomic_Mass                          (shortLabel='H'         ) &
               &                                /self%elementAtomicMasses             (           i           )
       end select
       countIonizationStates=countIonizationStates+self%elementAtomicNumbers(i)
    end do
    self%countStatesConvergence=min(int(dble(countIonizationStates)*(1.0d0-convergencePercentile),c_size_t)+1_c_size_t,countIonizationStates)
    ! Initialize photoionization cross section memoization.
    allocate(self%crossSectionPhotoIonizationPrevious(self%countElements,0:maxval(self%elementAtomicNumbers)-1))
    self%wavelengthPrevious                 =-1.0d0
    self%crossSectionPhotoIonizationPrevious=-1.0d0
    ! Build a root finder for thermal state.
    self%finder=rootFinder(                                                                  &
         &                 rootFunction                 =     atomicStateThermalBalance    , &
         &                 toleranceRelative            =     temperatureToleranceRelative , &
         &                 toleranceAbsolute            =     temperatureToleranceAbsolute , &
         &                 rangeExpandUpward            =     2.0d0                        , &
         &                 rangeExpandDownward          =     0.5d0                        , &
         &                 rangeExpandUpwardSignExpect  =     rangeExpandSignExpectNegative, &
         &                 rangeExpandDownwardSignExpect=     rangeExpandSignExpectPositive, &
         &                 rangeUpwardLimit             =     temperatureMaximum           , &
         &                 rangeDownwardLimit           =self%temperatureMinimum           , &
         &                 rangeExpandType              =     rangeExpandMultiplicative      &
         &                )
    return
  end function atomicConstructorInternal

  subroutine atomicDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferMatterAtomic} radiative transfer matter class.
    !!}
    implicit none
    type(radiativeTransferMatterAtomic), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"                       />
    <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"      />
    <objectDestructor name="self%atomicRecombinationRateRadiative_"       />
    <objectDestructor name="self%atomicRecombinationRateRadiativeCooling_"/>
    <objectDestructor name="self%atomicIonizationRateCollisional_"        />
    <objectDestructor name="self%atomicRecombinationRateDielectronic_"    />
    <objectDestructor name="self%atomicIonizationPotential_"              />
    <objectDestructor name="self%atomicExcitationRateCollisional_"        />
    <objectDestructor name="self%gauntFactor_"                            />
    !!]
    return
  end subroutine atomicDestructor

  subroutine atomicPropertyClass(self,properties)
    !!{
    Return the property class to use for the computational domain.
    !!}
    implicit none
    class(radiativeTransferMatterAtomic    ), intent(inout)              :: self
    class(radiativeTransferPropertiesMatter), intent(inout), allocatable :: properties
    !$GLC attributes unused :: self
    
    allocate(radiativeTransferPropertiesMatterAtomic :: properties)
    return
  end subroutine atomicPropertyClass
  
  subroutine atomicPopulateDomain(self,properties,integrator,onProcess)
    !!{
    Populate a computational domain cell with atomic matter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiativeTransferMatterAtomic           ), intent(inout) :: self
    class           (radiativeTransferPropertiesMatter       ), intent(inout) :: properties
    class           (computationalDomainVolumeIntegratorClass), intent(inout) :: integrator
    logical                                                   , intent(in   ) :: onProcess
    double precision                                                          :: densityCell
    integer                                                                   :: i

    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Initialize properties:
       ! + Assume that the atomic species are fully neutral initially.
       ! + Assign an initial temperature.
       ! + Initialize photoionization and photo heating rates to impossible values. These will be used as the "previous" value on
       !   the first iteration, guaranteeing that our first iteration is never judged to be converged.
       properties%iterationCount=+0
       properties%temperature   =self%temperatureMinimum
       properties%volume        =+integrator%volume()
       allocate(properties%elements(self%countElements))
       do i=1,self%countElements
          allocate(properties%elements(i)%photoIonizationRateHistory    (self%iterationAverageCount,0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%photoHeatingRateHistory       (self%iterationAverageCount,0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%ionizationStateFractionHistory(self%iterationAverageCount,0:self%elementAtomicNumbers(i)  ))
          allocate(properties%elements(i)%ionizationStateFraction       (                           0:self%elementAtomicNumbers(i)  ))
          allocate(properties%elements(i)%photoIonizationRate           (                           0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%photoHeatingRate              (                           0:self%elementAtomicNumbers(i)-1))
          properties%elements(i)%ionizationStateFraction       =      0.0d0
          properties%elements(i)%ionizationStateFraction    (0)=      1.0d0
          properties%elements(i)%photoIonizationRate           =-huge(0.0d0)
          properties%elements(i)%photoHeatingRate              =-huge(0.0d0)
          properties%elements(i)%photoIonizationRateHistory    =-huge(0.0d0)
          properties%elements(i)%photoHeatingRateHistory       =-huge(0.0d0)
          properties%elements(i)%ionizationStateFractionHistory=-huge(0.0d0)
       end do
       ! Compute the number density of this atomic species.
       if (onProcess) then
          densityCell                      =+integrator %integrate                    (atomicDensityIntegrand) &
               &                            /properties %volume
          properties%elements%densityNumber=+densityCell                                                       &
               &                            *self       %numberDensityMassDensityRatio
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return

  contains

    double precision function atomicDensityIntegrand(coordinates)
      !!{
      Integrand of atomic matter density.
      !!}
      use :: Coordinates, only : coordinate
      implicit none
      class(coordinate), intent(in   ) :: coordinates
      
      atomicDensityIntegrand= self%massDistribution_%density(coordinates)
      return
    end function atomicDensityIntegrand

  end subroutine atomicPopulateDomain

#ifdef USEMPI
  subroutine atomicBroadcastDomain(self,sendFromProcess,properties)
    !!{
    Broadcast populated computational domain properties to other MPI processes.
    !!}
    use :: Error        , only : Error_Report
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    integer                                   , intent(in   ) :: sendFromProcess
    class  (radiativeTransferPropertiesMatter), intent(inout) :: properties
    !$GLC attributes unused :: self

    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       call mpiSelf%broadcastData(sendFromProcess,properties%elements%densityNumber)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicBroadcastDomain
#endif
  
  subroutine atomicReset(self,properties)
    !!{
    Reset a computational domain cell prior to a new iteration.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferPropertiesMatter), intent(inout) :: properties
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Store the current photoionization and photoheating rates as the "previous" values, and then reset these rates to zero.
       do i=1,self%countElements
          properties%elements(i)%photoIonizationRatePrevious=properties%elements(i)%photoIonizationRate
          properties%elements(i)%photoHeatingRatePrevious   =properties%elements(i)%photoHeatingRate
          properties%elements(i)%photoIonizationRate        =0.0d0
          properties%elements(i)%photoHeatingRate           =0.0d0
       end do
       properties               %iterationCount             =properties            %iterationCount     +1
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicReset

  double precision function atomicAbsorptionCoefficient(self,properties,photonPacket)
    !!{
    Return the absorption coefficient for the given photon packet and matter properties.
    !!}
    implicit none
    class  (radiativeTransferMatterAtomic     ), intent(inout) :: self
    class  (radiativeTransferPropertiesMatter ), intent(inout) :: properties
    class  (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    integer                                                    :: i           , j
    
    atomicAbsorptionCoefficient=0.0d0
    do i=1,self%countElements
       do j=0,self%elementAtomicNumbers(i)-1 ! j=0 is neutral atom; j=1 is first ionized state, etc.
          atomicAbsorptionCoefficient=+     atomicAbsorptionCoefficient                                            &
               &                      +self%absorptionCoefficientSpecies(i,j,photonPacket%wavelength(),properties)
       end do
    end do
    return
  end function atomicAbsorptionCoefficient

  double precision function atomicCrossSectionPhotoionization(self,elementIndex,ionizationState,wavelength)
    !!{
    Compute the photoionization cross section for the given element and ionization state.
    !!}
    implicit none
    implicit none
    class           (radiativeTransferMatterAtomic), intent(inout) :: self
    integer                                        , intent(in   ) :: elementIndex  , ionizationState
    double precision                               , intent(in   ) :: wavelength
    integer                                                        :: k             , n              , &
         &                                                            l             , m              , &
         &                                                            countElectrons

    ! Reset memoized values if the wavelength has changed.
    if (wavelength /= self%wavelengthPrevious) then
       self%wavelengthPrevious                 =wavelength
       self%crossSectionPhotoIonizationPrevious=-1.0d0
    end if
    ! Use the memoized value if possible.
    if (self%crossSectionPhotoIonizationPrevious(elementIndex,ionizationState) < 0.0d0) then
       ! Determine the maximum sub-shell occupied by electrons. Sub-shell number is 1s=1, 2s=2, 2p=3, 3s=4, etc.
       countElectrons=0
       n             =0
       m             =0
       do while (countElectrons < self%elementAtomicNumbers(elementIndex)-ionizationState)
          n=n+1
          l=-1
          do while (l < n-1)
             l             =l             +1
             m             =m             +1
             countElectrons=countElectrons+4*l+2
             if (countElectrons >= self%elementAtomicNumbers(elementIndex)-ionizationState) exit
          end do
       end do
       self%crossSectionPhotoIonizationPrevious(elementIndex,ionizationState)=0.0d0
       do k=1,m
          self%crossSectionPhotoIonizationPrevious(elementIndex,ionizationState)=+self%crossSectionPhotoIonizationPrevious             (                                                              &
               &                                                                                                                                                                  elementIndex      , &
               &                                                                                                                                                                  ionizationState     &
               &                                                                                                                       )                                                              &
               &                                                                 +self%atomicCrossSectionIonizationPhoto_ %crossSection(                                                              &
               &                                                                                                                        atomicNumber   =self%elementAtomicNumbers(elementIndex     ), &
               &                                                                                                                        ionizationState=                          ionizationState+1 , &
               &                                                                                                                        shellNumber    =                          k                 , &
               &                                                                                                                        wavelength     =     wavelength                               &
               &                                                                                                                       )
       end do
    end if
    atomicCrossSectionPhotoionization=self%crossSectionPhotoIonizationPrevious(elementIndex,ionizationState)
    return
  end function atomicCrossSectionPhotoionization
  
  double precision function atomicAbsorptionCoefficientSpecies(self,elementIndex,ionizationState,wavelength,properties)
    !!{
    Compute the absorption coefficient for the given species and wavelength.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout) :: self
    integer                                            , intent(in   ) :: elementIndex               , ionizationState
    double precision                                   , intent(in   ) :: wavelength
    class           (radiativeTransferPropertiesMatter), intent(inout) :: properties
    double precision                                                   :: crossSectionPhotoIonization

    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Get the photoionization cross section.
       crossSectionPhotoIonization=self%crossSectionPhotoionization(elementIndex,ionizationState,wavelength)
       ! Compute the absorption coefficient.
       atomicAbsorptionCoefficientSpecies=+crossSectionPhotoIonization                                                &
            &                             *properties%elements(elementIndex)%densityNumber                            &
            &                             *properties%elements(elementIndex)%ionizationStateFraction(ionizationState) &
            &                             /centi                                                                      &
            &                             *megaParsec
    class default
       atomicAbsorptionCoefficientSpecies=0.0d0
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicAbsorptionCoefficientSpecies
  
  subroutine atomicAccumulatePhotonPacket(self,properties,photonPacket,absorptionCoefficient,lengthTraversed)
    !!{
    Accumulate a photon packet.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : luminositySolar  , megaParsec
    use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : metersToAngstroms, electronVolt
    implicit none
    class           (radiativeTransferMatterAtomic     ), intent(inout) :: self
    class           (radiativeTransferPropertiesMatter ), intent(inout) :: properties
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    double precision                                    , intent(in   ) :: absorptionCoefficient      , lengthTraversed
    double precision                                                    :: energyPhoton               , rateIonization             , &
         &                                                                 crossSectionPhotoIonization, atomicAbsorptionCoefficient
    integer                                                             :: i                          , j
    !$GLC attributes unused :: absorptionCoefficient

    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Accumulate the photoionization rate per unit volume (cm⁻³), and photoheating rate per unit volume (J cm⁻³) using the Lucy
       ! (1999; A&A; 344; 282; https://ui.adsabs.harvard.edu/abs/1999A%26A...344..282L; see section 3.4) methodology.
       energyPhoton=+plancksConstant               &
            &       *speedLight                    &
            &       *metersToAngstroms             &
            &       /photonPacket     %wavelength()
       do i=1,self%countElements
          do j=0,self%elementAtomicNumbers(i)-1 ! j=0 is neutral atom; j=1 is first ionized state, etc.
             crossSectionPhotoIonization                  =+self%crossSectionPhotoionization(i,j,photonPacket%wavelength())
             atomicAbsorptionCoefficient                  =+crossSectionPhotoIonization                                                                                                                       &
                  &                                        *properties                 %elements                  (i)%densityNumber                                                                           &
                  &                                        *properties                 %elements                  (i)%ionizationStateFraction(                                                          j)    &
                  &                                        /centi                                                                                                                                             &
                  &                                        *megaParsec 
             rateIonization                               =+photonPacket                                             %luminosity             (                                                           )    &
                  &                                        *luminositySolar                                                                                                                                   &
                  &                                        /energyPhoton                                                                                                                                      &
                  &                                        *atomicAbsorptionCoefficient                                                                                                                       &
                  &                                        *lengthTraversed                                                                                                                                   &
                  &                                        /properties                                               %volume                                                                                  &
                  &                                        *centi                                                                                                                                         **3 &
                  &                                        /megaParsec                                                                                                                                    **3
             properties%elements(i)%photoIonizationRate(j)=+properties                 %elements                  (i)%photoIonizationRate    (                                                          j)    &
                  &                                        +rateIonization
             properties%elements(i)%photoHeatingRate   (j)=+properties                 %elements                  (i)%photoHeatingRate       (                                                          j)    &
                  &                                        +rateIonization                                                                                                                                    &
                  &                                        *(                                                                                                                                                 &
                  &                                          +energyPhoton                                                                                                                                    &
                  &                                          -self                     %atomicIonizationPotential_   %potential              (self%elementAtomicNumbers(i),self%elementAtomicNumbers(i)-j)    &
                  &                                          *electronVolt                                                                                                                                    &
                  &                                         )
          end do
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicAccumulatePhotonPacket

  logical function atomicInteractWithPhotonPacket(self,properties,photonPacket)
    !!{
    Interact with a photon packet. In this case the photon packet is always absorbed.
    !!}
    implicit none
    class(radiativeTransferMatterAtomic     ), intent(inout) :: self
    class(radiativeTransferPropertiesMatter ), intent(inout) :: properties
    class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    !$GLC attributes unused :: self, properties, photonPacket
    
    atomicInteractWithPhotonPacket=.false.
    return
  end function atomicInteractWithPhotonPacket

#ifdef USEMPI
  subroutine atomicAccumulationReduction(self,properties)
    !!{
    Perform reduction of accumulated properties across MPI processes.
    !!}
    use :: Error        , only : Error_Report
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout)  :: self
    class  (radiativeTransferPropertiesMatter), intent(inout)  :: properties
    integer                                                    :: i
    
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       do i=1,self%countElements
          properties%elements(i)%photoIonizationRate=mpiSelf%sum(properties%elements(i)%photoIonizationRate)
          properties%elements(i)%photoHeatingRate   =mpiSelf%sum(properties%elements(i)%photoHeatingRate   )
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicAccumulationReduction
#endif

  double precision function atomicStateThermalBalance(temperature)
    !!{
    Root function used in finding the equilibrium temperature for thermal balance.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : boltzmannsConstant, electronMass, electronRadius, fineStructure, &
          &                                     speedLight
    use :: Numerical_Constants_Prefixes, only : centi
    implicit none
    double precision, intent(in   ) :: temperature
    double precision                :: chargeIonic
    integer                         :: i          , j

    atomicStateThermalBalance=0.0d0
    ! Iterate over all states of all elements, accumulating heating/cooling rates and computing the equilibrium ionization fractions.
    do i=1,self_%countElements
       do j=0,self_%elementAtomicNumbers(i)
          ! Accumulate the photo-heating rate.
          if (j < self_%elementAtomicNumbers(i)) atomicStateThermalBalance=+atomicStateThermalBalance                        &
               &                                                           +elementsPhotoRate        (i)%photoHeatingRate(j)
          !! Accumulate rate of collisional excitation cooling.
          if  (j < self_%elementAtomicNumbers(i))                                                                                                                  &
               & atomicStateThermalBalance=+atomicStateThermalBalance                                                                                              &
               &                           -densityNumberElectrons                                                                                                 &
               &                           *properties_%elements                        (i)%densityNumber                                                          &
               &                           *properties_%elements                        (i)%ionizationStateFraction(                              j              ) &
               &                           *self_      %atomicExcitationRateCollisional_   %coolingRate            (self_%elementAtomicNumbers(i),j+1,temperature) &
               &                           /centi**3
          !! Accumulate rate of Bremsstrahlung cooling.
          if (j > 0) then
             chargeIonic              =dble(j)
             atomicStateThermalBalance=+atomicStateThermalBalance                                                                          &
                  &                    -16.0d0                                                                                             &
                  &                    / 3.0d0                                                                                             &
                  &                    *sqrt(                                                                                              &
                  &                          +2.0d0                                                                                        &
                  &                          *Pi                                                                                           &
                  &                          /3.0d0                                                                                        &
                  &                         )                                                                                              &
                  &                    *chargeIonic**2                                                                                     &
                  &                    *properties_%elements(i)%densityNumber                                                              &
                  &                    *properties_%elements(i)%ionizationStateFraction      (                              j            ) &
                  &                    *densityNumberElectrons                                                                             &
                  &                    *electronRadius          **3                                                                        &
                  &                    *speedLight                                                                                         &
                  &                    /electronRadius                                                                                     &
                  &                    *sqrt(                                                                                              &
                  &                          +electronMass                                                                                 &
                  &                          *speedLight        **2                                                                        &
                  &                          *boltzmannsConstant                                                                           &
                  &                          *temperature                                                                                  &
                  &                         )                                                                                              &
                  &                    *fineStructure                                                                                      &
                  &                    *self_                  %gauntFactor_           %total(self_%elementAtomicNumbers(i),j,temperature) &
                  &                    /centi**3
             !! Accumulate the rates of cooling due to recombinations.
             atomicStateThermalBalance=+                                                       atomicStateThermalBalance                                                                     &
                  &                    -self_      %atomicRecombinationRateRadiativeCooling_   %rate                   (self_%elementAtomicNumbers(i),j,temperature,level=recombinationCase) &
                  &                    *boltzmannsConstant                                                                                                                                   &
                  &                    *temperature                                                                                                                                          &
                  &                    *densityNumberElectrons                                                                                                                               &
                  &                    *properties_%elements                                (i)%densityNumber                                                                                &
                  &                    *properties_%elements                                (i)%ionizationStateFraction(                              j                                    )
             atomicStateThermalBalance=+atomicStateThermalBalance                                                                                                                            &
                  &                    -0.0d0 !! TO DO - get dielectronic cooling rates.     
          end if
       end do
    end do
    return
  end function atomicStateThermalBalance
  
  subroutine atomicStateSolve(self,properties,status)
    !!{
    Solve for the state of the matter.
    !!}
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseA, recombinationCaseB
    use :: Display                             , only : displayIndent     , displayMessage    , displayUnindent      , verbosityLevelStandard
    use :: Error                               , only : Error_Report      , errorStatusFail   , errorStatusOutOfRange, errorStatusSuccess
    use :: Numerical_Roman_Numerals            , only : Roman_Numerals
    use :: ISO_Varying_String                  , only : operator(//)
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout) , target      :: self
    class           (radiativeTransferPropertiesMatter), intent(inout) , target      :: properties
    integer                                            , intent(  out) , optional    :: status
    double precision                                   , dimension(2)                :: temperaturePrevious                             , densityElectronsPrevious
    integer                                            , parameter                   :: countIterationMaximum                   =1000   , temperatureOscillatingCounts =30    , &
         &                                                                              electronsOscillatingCounts              =30
    double precision                                   , parameter                   :: degreesOfFreedom                        =3.0d+00
    type            (element                          ), dimension( : ), allocatable :: elementsPrevious                                , elementsReference
    double precision                                                                 :: rateRecombinationRadiative                      , rateRecombinationDielectronic       , &
         &                                                                              temperatureChangePrevious                       , temperatureEquilibrium              , &
         &                                                                              rateUpward                                      , rateDownward                        , &
         &                                                                              temperatureReference
    logical                                                                          :: converged                                       , temperatureOscillating              , &
         &                                                                              electronsOscillating                            , heatingNonZero
    integer                                                                          :: countIteration                                  , i                                   , &
         &                                                                              temperatureOscillatingCount                     , j                                   , &
         &                                                                              electronsOscillatingCount                       , statusTemperature                   , &
         &                                                                              exponentIonizationStateFractionMaximum          , exponentIonizationState
#ifdef RADTRANSDEBUG
    integer                                                                          :: handlerPrevious
#endif
    character       (len=256                         )                               :: message
    
#ifdef RADTRANSDEBUG
    handlerPrevious=Signal(8,debugAbort)
#endif
    if (present(status)) status=errorStatusSuccess
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Set module-scope pointers to self and properties.
       self_       => self
       properties_ => properties
       ! Append the accumulated photoionization/photoheating rates to the history arrays, compute the average over those
       ! histories, and replace the photoionization/photoheating rates with those averages.
       if (self%iterationAverageCount > 1) then
          do i=1,self%countElements
             properties%elements(i)%photoIonizationRateHistory(2:self%iterationAverageCount,:)=properties%elements(i)%photoIonizationRateHistory(1:self%iterationAverageCount-1,:)
             properties%elements(i)%photoHeatingRateHistory   (2:self%iterationAverageCount,:)=properties%elements(i)%photoHeatingRateHistory   (1:self%iterationAverageCount-1,:)
             properties%elements(i)%photoIonizationRateHistory(1                           ,:)=properties%elements(i)%photoIonizationRate       (                               :)
             properties%elements(i)%photoHeatingRateHistory   (1                           ,:)=properties%elements(i)%photoHeatingRate          (                               :)
             do j=0,self%elementAtomicNumbers(i)-1
                properties%elements(i)%photoIonizationRate(j)=+sum (properties%elements(i)%photoIonizationRateHistory(1:min(properties%iterationCount,self%iterationAverageCount),j)) &
                     &                                        /dble(                                                    min(properties%iterationCount,self%iterationAverageCount)   )
                properties%elements(i)%photoHeatingRate   (j)=+sum (properties%elements(i)%photoHeatingRateHistory   (1:min(properties%iterationCount,self%iterationAverageCount),j)) &
                     &                                        /dble(                                                    min(properties%iterationCount,self%iterationAverageCount)   )
             end do
          end do
       end if
       ! Create array of elements matched to that in the domain cell.
       allocate(elementsPrevious       (self%countElements))
       allocate(elementsReference      (self%countElements))
       allocate(elementsPhotoRate(self%countElements))
       do i=1,self%countElements
          allocate(elementsPrevious (i)%ionizationStateFraction       (                           lbound(properties%elements(i)%ionizationStateFraction       ,dim=1):ubound(properties%elements(i)%ionizationStateFraction       ,dim=1)))
          allocate(elementsReference(i)%ionizationStateFraction       (                           lbound(properties%elements(i)%ionizationStateFraction       ,dim=1):ubound(properties%elements(i)%ionizationStateFraction       ,dim=1)))
          allocate(elementsPhotoRate(i)%ionizationStateFraction       (                           lbound(properties%elements(i)%ionizationStateFraction       ,dim=1):ubound(properties%elements(i)%ionizationStateFraction       ,dim=1)))
          allocate(elementsPrevious (i)%photoIonizationRate           (                           lbound(properties%elements(i)%photoIonizationRate           ,dim=1):ubound(properties%elements(i)%photoIonizationRate           ,dim=1)))
          allocate(elementsReference(i)%photoIonizationRate           (                           lbound(properties%elements(i)%photoIonizationRate           ,dim=1):ubound(properties%elements(i)%photoIonizationRate           ,dim=1)))
          allocate(elementsPhotoRate(i)%photoIonizationRate           (                           lbound(properties%elements(i)%photoIonizationRate           ,dim=1):ubound(properties%elements(i)%photoIonizationRate           ,dim=1)))
          allocate(elementsPrevious (i)%photoHeatingRate              (                           lbound(properties%elements(i)%photoHeatingRate              ,dim=1):ubound(properties%elements(i)%photoHeatingRate              ,dim=1)))
          allocate(elementsReference(i)%photoHeatingRate              (                           lbound(properties%elements(i)%photoHeatingRate              ,dim=1):ubound(properties%elements(i)%photoHeatingRate              ,dim=1)))
          allocate(elementsPhotoRate(i)%photoHeatingRate              (                           lbound(properties%elements(i)%photoHeatingRate              ,dim=1):ubound(properties%elements(i)%photoHeatingRate              ,dim=1)))
          allocate(elementsPrevious (i)%photoIonizationRatePrevious   (                           lbound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1):ubound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1)))
          allocate(elementsReference(i)%photoIonizationRatePrevious   (                           lbound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1):ubound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1)))
          allocate(elementsPhotoRate(i)%photoIonizationRatePrevious   (                           lbound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1):ubound(properties%elements(i)%photoIonizationRatePrevious   ,dim=1)))
          allocate(elementsPrevious (i)%photoHeatingRatePrevious      (                           lbound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1):ubound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1)))
          allocate(elementsReference(i)%photoHeatingRatePrevious      (                           lbound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1):ubound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1)))
          allocate(elementsPhotoRate(i)%photoHeatingRatePrevious      (                           lbound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1):ubound(properties%elements(i)%photoHeatingRatePrevious      ,dim=1)))
          allocate(elementsPrevious (i)%photoIonizationRateHistory    (self%iterationAverageCount,lbound(properties%elements(i)%photoIonizationRateHistory    ,dim=2):ubound(properties%elements(i)%photoIonizationRateHistory    ,dim=2)))
          allocate(elementsReference(i)%photoIonizationRateHistory    (self%iterationAverageCount,lbound(properties%elements(i)%photoIonizationRateHistory    ,dim=2):ubound(properties%elements(i)%photoIonizationRateHistory    ,dim=2)))
          allocate(elementsPhotoRate(i)%photoIonizationRateHistory    (self%iterationAverageCount,lbound(properties%elements(i)%photoIonizationRateHistory    ,dim=2):ubound(properties%elements(i)%photoIonizationRateHistory    ,dim=2)))
          allocate(elementsPrevious (i)%photoHeatingRateHistory       (self%iterationAverageCount,lbound(properties%elements(i)%photoHeatingRateHistory       ,dim=2):ubound(properties%elements(i)%photoHeatingRateHistory       ,dim=2)))
          allocate(elementsReference(i)%photoHeatingRateHistory       (self%iterationAverageCount,lbound(properties%elements(i)%photoHeatingRateHistory       ,dim=2):ubound(properties%elements(i)%photoHeatingRateHistory       ,dim=2)))
          allocate(elementsPhotoRate(i)%photoHeatingRateHistory       (self%iterationAverageCount,lbound(properties%elements(i)%photoHeatingRateHistory       ,dim=2):ubound(properties%elements(i)%photoHeatingRateHistory       ,dim=2)))
          allocate(elementsPrevious (i)%ionizationStateFractionHistory(self%iterationAverageCount,lbound(properties%elements(i)%ionizationStateFractionHistory,dim=2):ubound(properties%elements(i)%ionizationStateFractionHistory,dim=2)))
          allocate(elementsReference(i)%ionizationStateFractionHistory(self%iterationAverageCount,lbound(properties%elements(i)%ionizationStateFractionHistory,dim=2):ubound(properties%elements(i)%ionizationStateFractionHistory,dim=2)))
          allocate(elementsPhotoRate(i)%ionizationStateFractionHistory(self%iterationAverageCount,lbound(properties%elements(i)%ionizationStateFractionHistory,dim=2):ubound(properties%elements(i)%ionizationStateFractionHistory,dim=2)))
       end do
       elementsPrevious =properties%elements
       elementsReference=properties%elements
       elementsPhotoRate=properties%elements
       ! Iterate until convergence is achieved.
       converged                  =all(properties%elements(:)%densityNumber <= 0.0d0)
       countIteration             =0
       temperaturePrevious        =properties%temperature
       temperatureReference       =properties%temperature
       densityElectronsPrevious   =0.0d0
       temperatureChangePrevious  =0.0d0
       temperatureOscillating     =.false.
       electronsOscillating       =.false.
       temperatureOscillatingCount=0
       electronsOscillatingCount  =0
       do i=1,self%countElements
          elementsPrevious (i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          elementsReference(i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          if (.not.converged) then
             elementsPrevious(i)%photoHeatingRate   =+properties%elements(i)%photoHeatingRate
             if (properties%elements(i)%densityNumber > 0.0d0) then
                elementsPrevious(i)%photoIonizationRate=+properties%elements(i)%photoIonizationRate &
                     &                                  /properties%elements(i)%densityNumber
             else
                elementsPrevious(i)%photoIonizationRate=+0.0d0
             end if
          end if
       end do
       do while (.not.converged)
          ! Compute total electron density and particle density.
          densityNumberElectrons=0.0d0
          do i=1,self%countElements
             do j=1,self%elementAtomicNumbers(i)
                densityNumberElectrons=+                       densityNumberElectrons     &
                     &                 +properties%elements(i)%densityNumber              &
                     &                 *properties%elements(i)%ionizationStateFraction(j) &
                     &                 *dble                                          (j)
             end do
          end do
          if (electronsOscillatingCount > 0)                                                       &
               & densityNumberElectrons=0.5d0*(densityNumberElectrons+densityElectronsPrevious(1))
          ! Compute rates of change of the ionization states and temperature.
          !! Initialize rates and computed photoionization and photoheating rates.
          heatingNonZero=.false.
          do i=1,self%countElements
             do j=0,self%elementAtomicNumbers(i)-1
                ! Compute the photoionization rate.
                if (properties%elements(i)%densityNumber > 0.0d0) then
                   elementsPhotoRate (i)%photoIonizationRate(j)=+properties%elements(i)%photoIonizationRate(j) &
                        &                                       /properties%elements(i)%densityNumber
                else
                   elementsPhotoRate (i)%photoIonizationRate(j)=+0.0d0
                end if
                ! Compute photoheating rates (in W cm⁻³).
                elementsPhotoRate    (i)%photoHeatingRate   (j)=+properties%elements(i)%photoHeatingRate   (j)
                if (elementsPhotoRate(i)%photoHeatingRate   (j) > 0.0d0) heatingNonZero=.true.
             end do
          end do
          ! If there is photoheating occurring, solve for the state of matter.
          if (heatingNonZero) then
             ! Iterate over all states of all elements, computing the equilibrium ionization fractions.
             do i=1,self%countElements
                do j=0,self%elementAtomicNumbers(i)
                   if (j == 0) then
                      ! Set the neutral ion fraction to an arbitrary value - we will compute ratios of ionization fractions and
                      ! renormalize at the end.
                      properties%elements(i)%ionizationStateFraction(0)=1.0d0
                   else
                      ! Compute the abundance of this ionization state relative to the lower state by requiring balance between
                      ! ionization and recombination rates.
                      !! Rate of upward transitions from photoionization.
                      if     (                                                                                    &
                           &              elementsReference(i)%ionizationStateFraction(j-1)  > 0.0d0              &
                           &  .and.                                                                               &
                           &   properties%elements         (i)%densityNumber                 > 0.0d0              &
                           &  .and.                                                                               &
                           &   +exponent(elementsPhotoRate (i)%photoIonizationRate    (j-1))                      &
                           &   -exponent(elementsReference (i)%ionizationStateFraction(j-1)) < maxExponent(0.0d0) &
                           & ) then
                         rateUpward=+elementsPhotoRate(i)%photoIonizationRate    (j-1) &
                              &     /elementsReference(i)%ionizationStateFraction(j-1)
                      else
                         rateUpward=+0.0d0
                      end if
                      ! Rate of upward transitions from collisional ionizations.
                      rateUpward=+                                      rateUpward                                                                    &
                           &     +                                      densityNumberElectrons                                                        &
                           &     *self%atomicIonizationRateCollisional_%rate                  (self%elementAtomicNumbers(i),j,properties%temperature)
                      !! Rates of downward transitions from radiative and dielectronic recombintions.
                      select case (self%elementAtomicNumbers(i))
                      case (1,2)
                         ! Hydrogen, helium.
                         !! TO DO - hard-coded for case-B recombination
                         recombinationCase=recombinationCaseB
                      case default
                         ! Metals.
                         recombinationCase=recombinationCaseA
                      end select
                      rateRecombinationRadiative   =+                                          densityNumberElectrons                                                                                &
                           &                        *self%atomicRecombinationRateRadiative_   %rate                  (self%elementAtomicNumbers(i),j,properties%temperature,level=recombinationCase)
                      rateRecombinationDielectronic=+                                          densityNumberElectrons                                                                                &
                           &                        *self%atomicRecombinationRateDielectronic_%rate                  (self%elementAtomicNumbers(i),j,properties%temperature                        )
                      !! Find the net rate of downward transitions.
                      rateDownward                 =+rateRecombinationRadiative                                                                                                                      &
                           &                        +rateRecombinationDielectronic
                      !! Compute the relative abundance required to balance upward and downward transitions.
                      if (rateDownward == 0.0d0) then
                         if (rateUpward > 0.0d0) then
                            ! If there are no downward transitions, but some upward transitions no solution is available. So simply
                            ! assume that the relative abundance is 1 - this should get corrected in subsequent iterations.
                            properties%elements(i)%ionizationStateFraction(j)=+properties%elements(i)%ionizationStateFraction(j-1)
                         else
                            ! If there are no downward transitions, and no upward transitions then the upper state should have zero
                            ! fraction.
                            properties%elements(i)%ionizationStateFraction(j)=+0.0d0
                         end if
                      else
                         ! Solve for the relative abundance.
                         !! Ensure that the relative abundances are in a range such that the new abundance will not overflow.
                         exponentIonizationState               =+     exponent(properties%elements(i)%ionizationStateFraction(j-1))
                         exponentIonizationStateFractionMaximum=+  maxExponent(0.0d0                                              ) &
                              &                                 -(                                                                  &
                              &                                   +   exponent(rateUpward                                         ) &
                              &                                   -   exponent(rateDownward                                       ) &
                              &                                  )
                         if (exponentIonizationState >= exponentIonizationStateFractionMaximum) &
                              & properties%elements(i)%ionizationStateFraction(0:j-1)=properties%elements(i)%ionizationStateFraction(0:j-1)/dble(radix(0.0d0))**(exponentIonizationState-exponentIonizationStateFractionMaximum+1)
                         !! Compute the relative abundance in this ionization state.
                         properties%elements(i)%ionizationStateFraction   (j)=+properties%elements(i)%ionizationStateFraction(j-1) &
                              &                                               *rateUpward                                          &
                              &                                               /rateDownward
                      end if
                   end if
                   ! Normalize all ionization states computed so far.
                   properties%elements(i)%ionizationStateFraction(0:j)=+    properties%elements(i)%ionizationStateFraction(0:j)  &
                        &                                              /sum(properties%elements(i)%ionizationStateFraction(0:j))
                end do
                ! Renormalize ionization states so that the summed fraction is 1.
                properties%elements(i)%ionizationStateFraction=+    properties%elements(i)%ionizationStateFraction  &
                     &                                         /sum(properties%elements(i)%ionizationStateFraction)
                ! If electron density is oscillating then bisect the current and previous ionization state fractions.
                if (electronsOscillatingCount > 0) &
                     & properties%elements(i)%ionizationStateFraction=0.5d0*(properties%elements(i)%ionizationStateFraction+elementsPrevious(i)%ionizationStateFraction)
             end do
             ! Decrement electron oscillation count.
             if (electronsOscillatingCount > 0) &
                  & electronsOscillatingCount=electronsOscillatingCount-1
             ! Scale heating rates by the relevant ionization state fraction.
             do i=1,self%countElements
                do j=0,self%elementAtomicNumbers(i)-1
                   if (elementsReference   (i)%ionizationStateFraction(j) > 0.0d0)                                                    &
                        & elementsPhotoRate(i)%photoHeatingRate       (j)=+properties%elements         (i)%photoHeatingRate       (j) &
                        &                                                 *properties%elements         (i)%ionizationStateFraction(j) &
                        &                                                 /           elementsReference(i)%ionizationStateFraction(j)
                end do
             end do
             ! Solve for temperature - updating only if we find a solution in range.
             temperatureEquilibrium=self%finder%find(rootGuess=properties%temperature,status=statusTemperature)
             if (statusTemperature == errorStatusSuccess) then
                if (countIteration > 1 .and. abs(log10(properties%temperature/temperatureEquilibrium)) > 3.0d0) then
                   ! If a large change in temperature occurs, geometrically bisect with the previous temperature.
                   properties%temperature=sqrt(temperatureEquilibrium*properties%temperature)
                else
                   properties%temperature=     temperatureEquilibrium
                end if
             else
                ! Temperature solution failed. Attempt to adjust the temperature slightly to allow us to find a solution on
                ! subsequent iterations.
                if      (atomicStateThermalBalance(self%temperatureMinimum) < 0.0d0) then
                   properties%temperature=max(properties%temperature/1.1d0,self%temperatureMinimum)
                else if (atomicStateThermalBalance(     temperatureMaximum) > 0.0d0) then
                   properties%temperature=min(properties%temperature*1.1d0,     temperatureMaximum)
                end if
             end if
             ! If temperature is oscillating bisect the temperature, and then decrease the count of steps for which we attempted
             ! to break out of oscillation.
             if (temperatureOscillatingCount > 0) then
                properties%temperature=0.5d0*(properties%temperature+temperaturePrevious(1))
                if (temperatureOscillatingCount > 0) temperatureOscillatingCount=+temperatureOscillatingCount &
                     &                                                           -1
             end if
          else
             ! No photo-heating - set to minimum temperature and fully neutral.
             statusTemperature     =errorStatusSuccess
             properties%temperature=self%temperatureMinimum
             do i=1,self%countElements
                properties%elements(i)%ionizationStateFraction   =0.0d0
                properties%elements(i)%ionizationStateFraction(0)=1.0d0
             end do
          end if
          ! Check for convergence.
          converged=       statusTemperature                              ==     errorStatusSuccess                                                                &
               &    .and.                                                                                                                                          &
               &     (                                                                                                                                             &
               &       abs(properties%temperature-temperaturePrevious(1)) <  max(temperatureToleranceAbsolute,temperatureToleranceRelative*properties%temperature) &
               &      .or.                                                                                                                                         &
               &       abs(properties%temperature-temperaturePrevious(1)) <      temperatureToleranceAbsolute                                                      &
               &     )
          do i=1,self%countElements
             if (.not.converged) exit
             converged=all(                                                                                                                                       &
                  &         abs(properties%elements(i)%ionizationStateFraction-elementsPrevious(i)%ionizationStateFraction)                                       &
                  &        <                                                                                                                                      &
                  &         max(ionizationStateFractionToleranceAbsolute,ionizationStateFractionToleranceRelative*properties%elements(i)%ionizationStateFraction) &
                  &       )
          end do
          ! Check for exceeding maximum iterations.
          countIteration=countIteration+1
          if (countIteration > countIterationMaximum .and. .not.converged) then
             if (present(status)) then
                status=errorStatusFail
                call self%historyUpdate(properties)
#ifdef RADTRANSDEBUG
                call Signal(8,handlerPrevious)
#endif
                deallocate(elementsPhotoRate)                
                return
             else
                ! No solution was found - report on the state of this domain cell.
                call displayIndent('failed domain cell state report',verbosityLevelStandard)
                write (message,'(a,e23.16         )') 'volume                                    = ',properties%volume
                call displayMessage(message,verbosityLevelStandard)
                write (message,'(a,e23.16         )') 'initial temperature                       = ',temperatureReference
                call displayMessage(message,verbosityLevelStandard)
                do i=1,self%countElements
                   call displayIndent('element: '//trim(adjustl(self%elements(i))),verbosityLevelStandard)
                   write (message,'(a,e23.16         )') 'density                                   = ',properties%elements(i)%densityNumber
                   call displayMessage(message,verbosityLevelStandard)
                   do j=0,self%elementAtomicNumbers(i)
                      call displayIndent('ion: '//trim(adjustl(self%elements(i)))//Roman_Numerals(j+1),verbosityLevelStandard)
                      if (j < self%elementAtomicNumbers(i)) then
                         write (message,'(a,e23.16,a,e23.16)') 'photoionization rate (current : previous) = ',properties%elements         (i)%photoIonizationRate    (j),' : ',properties%elements(i)%photoIonizationRatePrevious(j)
                         call displayMessage(message,verbosityLevelStandard)
                         write (message,'(a,e23.16,a,e23.16)') 'photoheating rate    (current : previous) = ',properties%elements         (i)%photoHeatingRate       (j),' : ',properties%elements(i)%photoHeatingRatePrevious   (j)
                         call displayMessage(message,verbosityLevelStandard)
                      end if
                      write    (message,'(a,e23.16,a,e23.16)') 'initial ionization state                  = ',           elementsReference(i)%ionizationStateFraction(j)
                      call displayMessage(message,verbosityLevelStandard)
                      call displayUnindent('done',verbosityLevelStandard)
                   end do
                   call displayUnindent('done',verbosityLevelStandard)
                end do
                call displayUnindent('done',verbosityLevelStandard)
#ifdef RADTRANSDEBUG
                call debugReport()                
#endif
                call Error_Report('solution not found'//{introspection:location})
             end if
          end if
          ! Check for oscillating temperature or electron density.
          if (countIteration >= temperatureOscillatingCounts) then
             if (temperatureOscillatingCount == 0) then
                temperatureOscillating=(properties%temperature-temperaturePrevious     (1))*(temperaturePrevious     (1)-temperaturePrevious     (2)) < 0.0d0
                ! If oscillation is detected set the number of steps for which we will attempt to break out of it.
                if (temperatureOscillating) temperatureOscillatingCount=temperatureOscillatingCounts
             end if
             if (electronsOscillatingCount == 0) then
                electronsOscillating  =(densityNumberElectrons-densityElectronsPrevious(1))*(densityElectronsPrevious(1)-densityElectronsPrevious(2)) < 0.0d0
                ! If oscillation is detected set the number of steps for which we will attempt to break out of it.
                if (electronsOscillating  ) electronsOscillatingCount  =electronsOscillatingCounts
             end if
          end if
          ! Store current state to previous state.
          do i=1,self%countElements
             elementsPrevious(i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          end do
          temperatureChangePrevious   =abs(properties%temperature                 -temperaturePrevious(1))
          temperaturePrevious      (2)=               temperaturePrevious     (1)
          temperaturePrevious      (1)=    properties%temperature
          densityElectronsPrevious (2)=               densityElectronsPrevious(1)
          densityElectronsPrevious (1)=               densityNumberElectrons
       end do
       call self%historyUpdate(properties)       
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
#ifdef RADTRANSDEBUG
    call Signal(8,handlerPrevious)
#endif
    deallocate(elementsPhotoRate)
    return

  contains
    
#ifdef RADTRANSDEBUG
    subroutine debugAbort()
      !!{
      Perform debug reporting and then exit.
      !!}
      implicit none
      
      call debugReport()
      call Abort      ()
      return
    end subroutine debugAbort
    
    subroutine debugReport()
      !!{
      Report debugging information.
      !!}
      use, intrinsic :: ISO_Fortran_Env, only : output_unit
      implicit none
      
      ! Debugging output. Writes the initial state for this failed case in a format that be directly pasted into the
      ! relevant test suite code.
      select type (properties)
      type is (radiativeTransferPropertiesMatterAtomic)
         write    (output_unit,'(a,i1)'         ) "properties            %iterationCount             = "    ,1
         write    (output_unit,'(a,e23.16)'     ) "properties            %volume                     = "    ,properties%volume
         write    (output_unit,'(a,e23.16)'     ) "properties            %temperature                = "    ,temperatureReference
         do i=1,self%countElements
            write (output_unit,'(a,i1,a,e23.16)') "properties%elements(",i,")%densityNumber              = ",properties%elements(i)%densityNumber     
            write (output_unit,'(a,i1,a,$)'     ) "properties%elements(",i,")%ionizationStateFraction    =["
            do j=0,self%elementAtomicNumbers(i)
               if (j < self%elementAtomicNumbers(i)) then
                  write (output_unit,'(e23.16,a,$)') elementsReference(i)%ionizationStateFraction    (j),','
               else
                  write (output_unit,'(e23.16,a)'  ) elementsReference(i)%ionizationStateFraction    (j),']'
               end if
            end do
            write (output_unit,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoIonizationRate        =["
            do j=0,self%elementAtomicNumbers(i)-1
               if (j < self%elementAtomicNumbers(i)-1) then
                  write (output_unit,'(e23.16,a,$)') elementsReference(i)%photoIonizationRate        (j),','
               else
                  write (output_unit,'(e23.16,a)'  ) elementsReference(i)%photoIonizationRate        (j),']'
               end if
            end do
            write (output_unit,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoHeatingRate           =["
            do j=0,self%elementAtomicNumbers(i)-1
               if (j < self%elementAtomicNumbers(i)-1) then
                  write (output_unit,'(e23.16,a,$)') elementsReference(i)%photoHeatingRate           (j),','
               else
                  write (output_unit,'(e23.16,a)'  ) elementsReference(i)%photoHeatingRate           (j),']'
               end if
            end do
            write (output_unit,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoIonizationRatePrevious=["
            do j=0,self%elementAtomicNumbers(i)-1
               if (j < self%elementAtomicNumbers(i)-1) then
                  write (output_unit,'(e23.16,a,$)') elementsReference(i)%photoIonizationRatePrevious(j),','
               else
                  write (output_unit,'(e23.16,a)'  ) elementsReference(i)%photoIonizationRatePrevious(j),']'
               end if
            end do
            write (output_unit,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoHeatingRatePrevious   =["
            do j=0,self%elementAtomicNumbers(i)-1
               if (j < self%elementAtomicNumbers(i)-1) then
                  write (output_unit,'(e23.16,a,$)') elementsReference(i)%photoHeatingRatePrevious   (j),','
               else
                  write (output_unit,'(e23.16,a)'  ) elementsReference(i)%photoHeatingRatePrevious   (j),']'
               end if
            end do
         end do
      end select
      return
    end subroutine debugReport
#endif

  end subroutine atomicStateSolve
  
  subroutine atomicHistoryUpdate(self,properties)
    !!{
    Update the ionization state history in a properties object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferPropertiesMatter), intent(inout) :: properties
    integer                                                   :: i         , j

    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       ! Append the ionization fractions to the history arrays, compute the average over those
       ! histories, and replace the ionization fractions with those averages.
       if (self%iterationAverageCount > 1) then
          do i=1,self%countElements
             properties%elements(i)%ionizationStateFractionHistory(2:self%iterationAverageCount,:)=properties%elements(i)%ionizationStateFractionHistory(1:self%iterationAverageCount-1,:)
             properties%elements(i)%ionizationStateFractionHistory(1                           ,:)=properties%elements(i)%ionizationStateFraction       (                               :)
             do j=0,self%elementAtomicNumbers(i)
                properties%elements(i)%ionizationStateFraction(j)=+sum (properties%elements(i)%ionizationStateFractionHistory(1:min(properties%iterationCount,self%iterationAverageCount),j)) &
                     &                                            /dble(                                                        min(properties%iterationCount,self%iterationAverageCount)   )
             end do
          end do
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicHistoryUpdate
  
#ifdef USEMPI
  subroutine atomicBroadcastState(self,sendFromProcess,properties)
    !!{
    Broadcast computational domain state to other MPI processes.
    !!}
    use :: Error        , only : Error_Report
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    integer                                   , intent(in   ) :: sendFromProcess
    class  (radiativeTransferPropertiesMatter), intent(inout) :: properties
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       do i=1,self%countElements
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoIonizationRate        )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoIonizationRatePrevious)
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoHeatingRate           )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoHeatingRatePrevious   )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%ionizationStateFraction    )
       end do
       call    mpiSelf%broadcastData(sendFromProcess,properties            %temperature                )
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicBroadcastState
#endif
  
  double precision function atomicConvergenceMeasure(self,properties)
    !!{
    Return a convergence measure for the atomic matter.
    !!}
    use :: Arrays_Search   , only : searchArray
    use :: Disparity_Ratios, only : Disparity_Ratio
    use :: Error           , only : Error_Report
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout)                          :: self
    class           (radiativeTransferPropertiesMatter), intent(inout)                          :: properties
    double precision                                   , dimension(self%countStatesConvergence) :: convergenceMeasures
    integer                                                                                     :: i                  , j
    integer         (c_size_t                         )                                         :: l                  , m
    double precision                                                                            :: convergenceMeasure
    
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       do l=1,self%countStatesConvergence
          convergenceMeasures(l)=-dble(self%countStatesConvergence+1_c_size_t-l)
       end do
       do i=1,self%countElements
          do j=0,self%elementAtomicNumbers(i)-1
             ! Only consider states with some significant population when computing convergence measure.
             if     (                                                                                               &
                  &   properties%elements(i)%ionizationStateFraction (j) > ionizationStateFractionToleranceAbsolute &
                  &  .and.                                                                                          &
                  &   properties%elements(i)%photoHeatingRatePrevious(j) > 0.0d0                                    &
                  & ) then
                convergenceMeasure=max(                                                                                                                      &
                     &                 Disparity_Ratio(properties%elements(i)%photoIonizationRate(j),properties%elements(i)%photoIonizationRatePrevious(j)), & 
                     &                 Disparity_Ratio(properties%elements(i)%photoHeatingRate   (j),properties%elements(i)%photoHeatingRatePrevious   (j))  &
                     &                )
                if (convergenceMeasure > convergenceMeasures(1)) then
                   ! A new highest-k convergence measure has been found. Insert it into our list.
                   if (convergenceMeasure > convergenceMeasures(self%countStatesConvergence)) then
                      l=self%countStatesConvergence
                   else
                      l=searchArray(convergenceMeasures,convergenceMeasure)
                   end if
                   if (l > 1) then
                      do m=1,l-1
                         convergenceMeasures(m)=convergenceMeasures(m+1)
                      end do
                   end if
                   convergenceMeasures(l)=convergenceMeasure
                end if

             end if
          end do
       end do
       atomicConvergenceMeasure=max(convergenceMeasures(1),1.0d0)
       class default
       atomicConvergenceMeasure=0.0d0
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicConvergenceMeasure
  
  double precision function atomicOutputProperty(self,properties,output)
    !!{
    Return a scalar property to be output.
    !!}
    use :: Error                     , only : Error_Report
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferPropertiesMatter), intent(inout) :: properties
    integer(c_size_t                         ), intent(in   ) :: output
    integer(c_size_t                         )                :: output_                     , outputCountElement     , &
         &                                                       offsetPhotoionizationRates  , offsetPhotoheatingRates, &
         &                                                       offsetAbsorptionCoefficients
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferPropertiesMatterAtomic)
       if      (output == 1_c_size_t                                   ) then
          atomicOutputProperty=properties%temperature
       else if (output >= 2_c_size_t .and. output <= self%countOutputs_) then
          output_                  =output                                                    -1_c_size_t
          i                        =                                                           1
          outputCountElement       =                              self%elementAtomicNumbers(i)+1_c_size_t
          if (self%outputRates                 )                                                                    &
               & outputCountElement=outputCountElement+2_c_size_t*self%elementAtomicNumbers(i)
          if (self%outputAbsorptionCoefficients)                                                                    &
               & outputCountElement=outputCountElement+           self%elementAtomicNumbers(i)
          do while (output_ > outputCountElement)
             output_                  =output_                                                  -outputCountElement
             i                        =i                                                         +1
             outputCountElement       =                              self%elementAtomicNumbers(i)+1_c_size_t
             if (self%outputRates                )                                                                  &
                  & outputCountElement=outputCountElement+2_c_size_t*self%elementAtomicNumbers(i)
             if (self%outputAbsorptionCoefficients)                                                                 &
                  & outputCountElement=outputCountElement+           self%elementAtomicNumbers(i)
          end do
          offsetPhotoionizationRates  =1_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
          offsetPhotoheatingRates     =2_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
          offsetAbsorptionCoefficients=1_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
          if (self%outputRates                 )                                                                    &
               offsetAbsorptionCoefficients=offsetAbsorptionCoefficients+2_c_size_t*self%elementAtomicNumbers(i)
          if      (                                         output_                             == 1_c_size_t                  ) then
             atomicOutputProperty=properties%elements(i)%densityNumber
          else if (                                         output_-1_c_size_t                  <= self%elementAtomicNumbers(i)) then
             atomicOutputProperty=properties%elements(i)%ionizationStateFraction(output_-1_c_size_t                )
          else if (self%outputRates                   .and. output_-offsetPhotoionizationRates  <= self%elementAtomicNumbers(i)) then
             atomicOutputProperty=properties%elements(i)%photoIonizationRate    (output_-offsetPhotoionizationRates)
          else if (self%outputRates                  .and. output_-offsetPhotoheatingRates      <= self%elementAtomicNumbers(i)) then
             atomicOutputProperty=properties%elements(i)%photoHeatingRate       (output_-offsetPhotoheatingRates   )
          else if (self%outputAbsorptionCoefficients .and. output_-offsetAbsorptionCoefficients <= self%elementAtomicNumbers(i)) then
             atomicOutputProperty=self%absorptionCoefficientSpecies(i,int(output_-offsetAbsorptionCoefficients-1_c_size_t),(1.0d0-1.0d-3)*lymanSeriesLimitWavelengthHydrogen_atomic,properties)
          else
             atomicOutputProperty=0.0d0
             call Error_Report('output is out of range (this should not happen)'//{introspection:location})
          end if
       else
          atomicOutputProperty=0.0d0
          call Error_Report('output is out of range'//{introspection:location})
       end if
    class default
       atomicOutputProperty=0.0d0
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicOutputProperty
  
  function atomicCountOutputs(self)
    !!{
    Return the number of scalar properties to output.
    !!}
    implicit none
    integer(c_size_t                     )                :: atomicCountOutputs
    class  (radiativeTransferMatterAtomic), intent(inout) :: self
    !$GLC attributes unused :: self

    atomicCountOutputs=self%countOutputs_
    return
  end function atomicCountOutputs

  function atomicOutputName(self,output)
    !!{
    Return the name of the scalar property to be output.
    !!}
    use :: Error                   , only : Error_Report
    use :: ISO_Varying_String      , only : operator(//)  , var_str
    use :: Numerical_Roman_Numerals, only : Roman_Numerals
    implicit none
    type   (varying_string               )                :: atomicOutputName
    class  (radiativeTransferMatterAtomic), intent(inout) :: self
    integer(c_size_t                     ), intent(in   ) :: output
    integer(c_size_t                     )                :: output_                     , outputCountElement     , &
         &                                                   offsetPhotoionizationRates  , offsetPhotoheatingRates, &
         &                                                   offsetAbsorptionCoefficients
    integer                                               :: i
    
    if      (output == 1_c_size_t                                   ) then
       atomicOutputName=var_str('temperature')
    else if (output >= 2_c_size_t .and. output <= self%countOutputs_) then
       output_                  =output                                                    -1_c_size_t
       i                        =                                                           1
       outputCountElement       =                              self%elementAtomicNumbers(i)+1_c_size_t
       if (self%outputRates                 )                                                                     &
            & outputCountElement=outputCountElement+2_c_size_t*self%elementAtomicNumbers(i)
       if (self%outputAbsorptionCoefficients)                                                                     &
            & outputCountElement=outputCountElement+           self%elementAtomicNumbers(i)
       do while (output_ > outputCountElement)
          output_                  =output_                                                  -outputCountElement
          i                        =i                                                         +1
          outputCountElement       =                              self%elementAtomicNumbers(i)+1_c_size_t
          if (self%outputRates                 )                                                                  &
               & outputCountElement=outputCountElement+2_c_size_t*self%elementAtomicNumbers(i)
          if (self%outputAbsorptionCoefficients)                                                                  &
               & outputCountElement=outputCountElement+           self%elementAtomicNumbers(i)
       end do
       offsetPhotoionizationRates  =1_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
       offsetPhotoheatingRates     =2_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
       offsetAbsorptionCoefficients=1_c_size_t*self%elementAtomicNumbers(i)+1_c_size_t
       if (self%outputRates                 )                                                                    &
            offsetAbsorptionCoefficients=offsetAbsorptionCoefficients+2_c_size_t*self%elementAtomicNumbers(i)
       if      (                                         output_                             == 1_c_size_t                  ) then
          atomicOutputName=var_str('densityNumber'        )//trim(adjustl(self%elements(i)))
       else if (                                         output_-1_c_size_t                  <= self%elementAtomicNumbers(i)) then
          atomicOutputName=var_str('fraction'             )//trim(adjustl(self%elements(i)))//Roman_Numerals(int(output_                             ))
       else if (self%outputRates                   .and. output_-offsetPhotoionizationRates  <= self%elementAtomicNumbers(i)) then
          atomicOutputName=var_str('photoIonizationRate'  )//trim(adjustl(self%elements(i)))//Roman_Numerals(int(output_-offsetPhotoionizationRates  ))
       else if (self%outputRates                  .and. output_-offsetPhotoheatingRates      <= self%elementAtomicNumbers(i)) then
          atomicOutputName=var_str('photoHeatingRate'     )//trim(adjustl(self%elements(i)))//Roman_Numerals(int(output_-offsetPhotoheatingRates     ))
       else if (self%outputAbsorptionCoefficients .and. output_-offsetAbsorptionCoefficients <= self%elementAtomicNumbers(i)) then
          atomicOutputName=var_str('absorptionCoefficient')//trim(adjustl(self%elements(i)))//Roman_Numerals(int(output_-offsetAbsorptionCoefficients))
       else
          atomicOutputName=var_str(''                   )
          call Error_Report('output is out of range (this should not happen)'//{introspection:location})
       end if
    else
       atomicOutputName=var_str(''                )
       call Error_Report('output is out of range'//{introspection:location})
    end if
    return
  end function atomicOutputName

  double precision function atomicRecombinationRateHydrogen(self,properties)
    !!{
    Return the total recombination rate for the atomic matter.
    !!}
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use :: Error                               , only : Error_Report
    use :: Numerical_Constants_Astronomical    , only : megaParsec
    use :: Numerical_Constants_Prefixes        , only : centi
    implicit none
    class           (radiativeTransferMatterAtomic          ), intent(inout) :: self
    type            (radiativeTransferPropertiesMatterAtomic), intent(inout) :: properties
    double precision                                                         :: densityNumberElectrons
    integer                                                                  :: i                     , j

    densityNumberElectrons=0.0d0
    do i=1,self%countElements
       do j=1,self%elementAtomicNumbers(i)
          densityNumberElectrons=+                       densityNumberElectrons     &
               &                 +properties%elements(i)%densityNumber              &
               &                 *properties%elements(i)%ionizationStateFraction(j) &
               &                 *dble                                          (j)
       end do
    end do
    if (self%indexHydrogen > 0) then
       atomicRecombinationRateHydrogen=+                                                                 densityNumberElectrons                                                       &
            &                          *properties%elements                         (self%indexHydrogen)%densityNumber                                                                &
            &                          *properties%elements                         (self%indexHydrogen)%ionizationStateFraction(  1                                                ) &
            &                          *self      %atomicRecombinationRateRadiative_                    %rate                   (1,1,properties%temperature,level=recombinationCaseB) &
            &                          *properties%volume                                                                                                                             &
            &                          *(                                                                                                                                             &
            &                            +megaParsec                                                                                                                                  &
            &                            /centi                                                                                                                                       &
            &                           )**3
    else
       atomicRecombinationRateHydrogen=+0.0d0
       call Error_Report('hydrogren is not present'//{introspection:location})
    end if
    return
  end function atomicRecombinationRateHydrogen
