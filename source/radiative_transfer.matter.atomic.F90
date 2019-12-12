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

  use :: Atomic_Cross_Sections_Ionization_Photo , only : atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Ionization_Potentials           , only : atomicIonizationPotentialClass
  use :: Atomic_Radiation_Gaunt_Factors         , only : gauntFactorClass
  use :: Atomic_Rates_Excitation_Collisional    , only : atomicExcitationRateCollisionalClass
  use :: Atomic_Rates_Ionization_Collisional    , only : atomicIonizationRateCollisionalClass
  use :: Atomic_Rates_Recombination_Dielectronic, only : atomicRecombinationRateDielectronicClass
  use :: Atomic_Rates_Recombination_Radiative   , only : atomicRecombinationRateRadiativeClass
  use :: Mass_Distributions                     , only : massDistributionClass
  
  !# <radiativeTransferMatter name="radiativeTransferMatterAtomic">
  !#  <description>A task which performs radiative transfer.</description>
  !# </radiativeTransferMatter>
  type, extends(radiativeTransferMatterClass) :: radiativeTransferMatterAtomic
     !% Implementation of a radiative transfer matter class for atomic matter.
     private
     class           (massDistributionClass                   ), pointer                   :: massDistribution_                    => null()
     class           (atomicCrossSectionIonizationPhotoClass  ), pointer                   :: atomicCrossSectionIonizationPhoto_   => null()
     class           (atomicRecombinationRateRadiativeClass   ), pointer                   :: atomicRecombinationRateRadiative_    => null()
     class           (atomicIonizationRateCollisionalClass    ), pointer                   :: atomicIonizationRateCollisional_     => null()
     class           (atomicRecombinationRateDielectronicClass), pointer                   :: atomicRecombinationRateDielectronic_ => null()
     class           (atomicIonizationPotentialClass          ), pointer                   :: atomicIonizationPotential_           => null()
     class           (atomicExcitationRateCollisionalClass    ), pointer                   :: atomicExcitationRateCollisional_     => null()
     class           (gauntFactorClass                        ), pointer                   :: gauntFactor_                         => null()
     integer                                                                               :: indexAbundancePattern                         , iterationAverageCount, &
          &                                                                                   countElements                                 , indexHydrogen
     integer         (c_size_t                                )                            :: countOutputs_
     double precision                                                                      :: metallicity                                   , temperatureMinimum
     double precision                                          , allocatable, dimension(:) :: numberDensityMassDensityRatio                 , elementAtomicMasses
     type            (varying_string                          )                            :: abundancePattern
     character       (len=2                                   ), allocatable, dimension(:) :: elements
     integer                                                   , allocatable, dimension(:) :: elementAtomicNumbers
   contains
     !@ <objectMethods>
     !@   <object>radiativeTransferMatterAtomic</object>
     !@   <objectMethod>
     !@     <method>recombinationRateHydrogen</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless type(radiativeTransferMatterPropertiesAtomic)\textgreater} properties</arguments>
     !@     <description>Return the total rate of recombinations (in units of s$^{-1}$).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                              atomicDestructor
     procedure :: propertyClass             => atomicPropertyClass
     procedure :: populateDomain            => atomicPopulateDomain
     procedure :: reset                     => atomicReset
     procedure :: absorptionCoefficient     => atomicAbsorptionCoefficient
     procedure :: accumulatePhotonPacket    => atomicAccumulatePhotonPacket
     procedure :: interactWithPhotonPacket  => atomicInteractWithPhotonPacket
     procedure :: stateSolve                => atomicStateSolve
     procedure :: convergenceMeasure        => atomicConvergenceMeasure
     procedure :: outputProperty            => atomicOutputProperty
     procedure :: countOutputs              => atomicCountOutputs
     procedure :: outputName                => atomicOutputName
     procedure :: recombinationRateHydrogen => atomicRecombinationRateHydrogen
#ifdef USEMPI
     procedure :: accumulationReduction    => atomicAccumulationReduction
     procedure :: broadcastDomain          => atomicBroadcastDomain
     procedure :: broadcastState           => atomicBroadcastState
#endif
  end type radiativeTransferMatterAtomic
  
  interface radiativeTransferMatterAtomic
     !% Constructors for the {\normalfont \ttfamily atomic} radiative transfer matter class.
     module procedure atomicConstructorParameters
     module procedure atomicConstructorInternal
  end interface radiativeTransferMatterAtomic

  type, public :: element
     !% Type used to store elemental states.
     double precision                              :: densityNumber
     double precision, dimension(:,:), allocatable :: photoIonizationRateHistory , photoHeatingRateHistory
     double precision, dimension(:  ), allocatable :: photoIonizationRate        , photoHeatingRate        , &
          &                                           photoIonizationRatePrevious, photoHeatingRatePrevious, &
          &                                           ionizationStateFraction
  end type element
  
  type, extends(radiativeTransferMatterProperties), public :: radiativeTransferMatterPropertiesAtomic
     !% Radiative transfer matter properties class for atomic matter.
     integer                                              :: iterationCount
     double precision                                     :: volume        , temperature
     type            (element), dimension(:), allocatable :: elements
  end type radiativeTransferMatterPropertiesAtomic
  
contains

  function atomicConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily atomic} radiative transfer matter class which takes a parameter set as input.
    use :: Input_Parameters                , only : inputParameters , inputParameter
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (radiativeTransferMatterAtomic           )                             :: self
    type            (inputParameters                         ), intent(inout)              :: parameters
    class           (massDistributionClass                   ), pointer                    :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass  ), pointer                    :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass   ), pointer                    :: atomicRecombinationRateRadiative_
    class           (atomicIonizationRateCollisionalClass    ), pointer                    :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass), pointer                    :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass          ), pointer                    :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass    ), pointer                    :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                        ), pointer                    :: gauntFactor_
    character       (len=2                                   ), dimension(:) , allocatable :: elements
    integer                                                                                :: iterationAverageCount    
    double precision                                                                       :: temperatureMinimum                  , metallicity
    type            (varying_string                          )                             :: abundancePattern
    
    !# <inputParameter>
    !#   <name>iterationAverageCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>5</defaultValue>
    !#   <description>The number of iterations over which to average the photoionization rate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>temperatureMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>3.0d0</defaultValue>
    !#   <description>The minimum temperature that matter is allowed to reach in the case of zero photoheating.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>abundancePattern</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('solar')</defaultValue>
    !#   <description>The abundance pattern to use.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>metallicity</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>metallicitySolar</defaultValue>
    !#   <description>The metallicity to use.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    if (parameters%isPresent('elements')) then
       allocate(elements(parameters%count('elements')))
    else
       allocate(elements(1                           ))
    end if
    !# <inputParameter>
    !#   <name>elements</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>['H']</defaultValue>
    !#   <description>The elements to include.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="massDistribution"                    name="massDistribution_"                    source="parameters"/>
    !# <objectBuilder class="atomicCrossSectionIonizationPhoto"   name="atomicCrossSectionIonizationPhoto_"   source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateRadiative"    name="atomicRecombinationRateRadiative_"    source="parameters"/>
    !# <objectBuilder class="atomicIonizationRateCollisional"     name="atomicIonizationRateCollisional_"     source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateDielectronic" name="atomicRecombinationRateDielectronic_" source="parameters"/>
    !# <objectBuilder class="atomicIonizationPotential"           name="atomicIonizationPotential_"           source="parameters"/>
    !# <objectBuilder class="atomicExcitationRateCollisional"     name="atomicExcitationRateCollisional_"     source="parameters"/>
    !# <objectBuilder class="gauntFactor"                         name="gauntFactor_"                         source="parameters"/>
    self=radiativeTransferMatterAtomic(abundancePattern,metallicity,elements,iterationAverageCount,temperatureMinimum,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_)
    !# <objectDestructor name="massDistribution_"                   />
    !# <objectDestructor name="atomicCrossSectionIonizationPhoto_"  />
    !# <objectDestructor name="atomicRecombinationRateRadiative_"   />
    !# <objectDestructor name="atomicIonizationRateCollisional_"    />
    !# <objectDestructor name="atomicRecombinationRateDielectronic_"/>
    !# <objectDestructor name="atomicIonizationPotential_"          />
    !# <objectDestructor name="atomicExcitationRateCollisional_"    />
    !# <objectDestructor name="gauntFactor_"                        />
    return
  end function atomicConstructorParameters

  function atomicConstructorInternal(abundancePattern,metallicity,elements,iterationAverageCount,temperatureMinimum,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily atomic} radiative transfer matter class.
    use :: Abundances_Structure            , only : abundances              , metallicityTypeLinearByMassSolar, adjustElementsReset, Abundances_Index_From_Name
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup, Atomic_Abundance                , Atomic_Mass        , Atomic_Number
    use :: ISO_Varying_String              , only : char
    use :: Numerical_Constants_Astronomical, only : massSolar               , megaParsec                      , metallicitySolar
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: String_Handling                 , only : String_Count_Words      , String_Split_Words
    implicit none
    type            (radiativeTransferMatterAtomic           )                              :: self
    integer                                                   , intent(in   )               :: iterationAverageCount
    double precision                                          , intent(in   )               :: temperatureMinimum                  , metallicity
    type            (varying_string                          ), intent(in   )               :: abundancePattern
    character       (len=2                                   ), intent(in   ), dimension(:) :: elements
    class           (massDistributionClass                   ), intent(in   ), target       :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass  ), intent(in   ), target       :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass   ), intent(in   ), target       :: atomicRecombinationRateRadiative_
    class           (atomicIonizationRateCollisionalClass    ), intent(in   ), target       :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass), intent(in   ), target       :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass          ), intent(in   ), target       :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass    ), intent(in   ), target       :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                        ), intent(in   ), target       :: gauntFactor_
    double precision                                          , allocatable  , dimension(:) :: abundancesRelative
    double precision                                                                        :: numberDensityMassDensityRatioHydrogen, numberDensityMassDensityRatioHelium
    type            (abundances                              )                              :: abundances_
    integer                                                                                 :: i
    !# <constructorAssign variables="abundancePattern, metallicity, elements, iterationAverageCount, temperatureMinimum, *massDistribution_, *atomicCrossSectionIonizationPhoto_, *atomicRecombinationRateRadiative_, *atomicIonizationRateCollisional_, *atomicRecombinationRateDielectronic_, *atomicIonizationPotential_, *atomicExcitationRateCollisional_, *gauntFactor_"/>

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
    allocate(abundancesRelative(self%countElements))
    call abundances_%serialize(abundancesRelative)
    do i=1,self%countElements
       self%elementAtomicNumbers(i)=Atomic_Number(shortLabel=trim(self%elements(i)))
       self%elementAtomicMasses (i)=Atomic_Mass  (shortLabel=trim(self%elements(i)))
       self%countOutputs_=self%countOutputs_+self%elementAtomicNumbers(i)+1_c_size_t
       select case (self%elementAtomicNumbers(i))
       case (1)     ! Hydrogen
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHydrogen
          self%indexHydrogen                   = i
       case (2)     ! Helium
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHelium
       case default ! Metals
          self%numberDensityMassDensityRatio(i)=+numberDensityMassDensityRatioHydrogen                                                                &
               &                                *abundancesRelative                   (           Abundances_Index_From_Name(trim(self%elements(i)))) &
               &                                *Atomic_Mass                          (shortLabel='H'                                               ) &
               &                                /self%elementAtomicMasses             (                                                         i   )
       end select
    end do
    return
  end function atomicConstructorInternal

  subroutine atomicDestructor(self)
    !% Destructor for the {\normalfont \ttfamily atomic} radiative transfer matter class.
    implicit none
    type(radiativeTransferMatterAtomic), intent(inout) :: self

    !# <objectDestructor name="self%massDistribution_"                   />
    !# <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"  />
    !# <objectDestructor name="self%atomicRecombinationRateRadiative_"   />
    !# <objectDestructor name="self%atomicIonizationRateCollisional_"    />
    !# <objectDestructor name="self%atomicRecombinationRateDielectronic_"/>
    !# <objectDestructor name="self%atomicIonizationPotential_"          />
    !# <objectDestructor name="self%atomicExcitationRateCollisional_"    />
    !# <objectDestructor name="self%gauntFactor_"                        />
    return
  end subroutine atomicDestructor

  subroutine atomicPropertyClass(self,properties)
    !% Return the property class to use for the computational domain.
    implicit none
    class(radiativeTransferMatterAtomic    ), intent(inout)              :: self
    class(radiativeTransferMatterProperties), intent(inout), allocatable :: properties
    !GCC$ attributes unused :: self
    
    allocate(radiativeTransferMatterPropertiesAtomic :: properties)
    return
  end subroutine atomicPropertyClass
  
  subroutine atomicPopulateDomain(self,properties,integrator,onProcess)
    !% Populate a computational domain cell with atomic matter.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (radiativeTransferMatterAtomic           ), intent(inout) :: self
    class           (radiativeTransferMatterProperties       ), intent(inout) :: properties
    class           (computationalDomainVolumeIntegratorClass), intent(inout) :: integrator
    logical                                                   , intent(in   ) :: onProcess
    double precision                                                          :: densityCell
    integer                                                                   :: i

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Initialize properties:
       ! + Assume that the atomic species are fully neutral initially.
       ! + Assign an initial temperature.
       ! + Initialize photoionization and photo heating rates to impossible values. These will be used as the "previous" value on
       !   the first iteration, guaranteeing that our first iteration is never judged to be converged.
       properties%iterationCount=+0
       properties%temperature   =1.0d4
       properties%volume        =+integrator%volume()
       allocate(properties%elements(self%countElements))
       do i=1,self%countElements
          allocate(properties%elements(i)%photoIonizationRateHistory(self%iterationAverageCount,0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%photoHeatingRateHistory   (self%iterationAverageCount,0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%ionizationStateFraction   (                           0:self%elementAtomicNumbers(i)  ))
          allocate(properties%elements(i)%photoIonizationRate       (                           0:self%elementAtomicNumbers(i)-1))
          allocate(properties%elements(i)%photoHeatingRate          (                           0:self%elementAtomicNumbers(i)-1))
          properties%elements(i)%ionizationStateFraction      =      0.0d0
          properties%elements(i)%ionizationStateFraction   (0)=      1.0d0
          properties%elements(i)%photoIonizationRate          =-huge(0.0d0)
          properties%elements(i)%photoHeatingRate             =-huge(0.0d0)
          properties%elements(i)%photoIonizationRateHistory   =-huge(0.0d0)
          properties%elements(i)%photoHeatingRateHistory      =-huge(0.0d0)
       end do
       ! Compute the number density of this atomic species.
       if (onProcess) then
          densityCell                      =+integrator %integrate                    (atomicDensityIntegrand) &
               &                            /properties %volume
          properties%elements%densityNumber=+densityCell                                                       &
               &                            *self       %numberDensityMassDensityRatio
       end if
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return

  contains

    double precision function atomicDensityIntegrand(coordinates)
      !% Integrand of atomic matter density.
      use :: Coordinates, only : coordinate
      implicit none
      class(coordinate), intent(in   ) :: coordinates
      
      atomicDensityIntegrand= self%massDistribution_%density(coordinates)
      return
    end function atomicDensityIntegrand

  end subroutine atomicPopulateDomain

#ifdef USEMPI
  subroutine atomicBroadcastDomain(self,sendFromProcess,properties)
    !% Broadcast populated computational domain properties to other MPI processes.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: MPI_Utilities   , only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    integer                                   , intent(in   ) :: sendFromProcess
    class  (radiativeTransferMatterProperties), intent(inout) :: properties
    !GCC$ attributes unused :: self

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       call mpiSelf%broadcastData(sendFromProcess,properties%elements%densityNumber)
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicBroadcastDomain
#endif
  
  subroutine atomicReset(self,properties)
    !% Reset a computational domain cell prior to a new iteration.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferMatterProperties), intent(inout) :: properties
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Store the current photoionization and photoheating rates as the "previous" values, and then reset these rates to zero.
       do i=1,self%countElements
          properties%elements(i)%photoIonizationRatePrevious=properties%elements(i)%photoIonizationRate
          properties%elements(i)%photoHeatingRatePrevious   =properties%elements(i)%photoHeatingRate
          properties%elements(i)%photoIonizationRate        =0.0d0
          properties%elements(i)%photoHeatingRate           =0.0d0
       end do
       properties               %iterationCount             =properties            %iterationCount     +1
       class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicReset

  double precision function atomicAbsorptionCoefficient(self,properties,photonPacket)
    !% Return the absorption coefficient for the given photon packet and matter properties.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : angstromsPerMeter
    implicit none
    class           (radiativeTransferMatterAtomic     ), intent(inout) :: self
    class           (radiativeTransferMatterProperties ), intent(inout) :: properties
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    double precision                                                    :: crossSectionPhotoIonization
    integer                                                             :: i                          , j             , &
         &                                                                 k                          , n             , &
         &                                                                 l                          , m             , &
         &                                                                 countElectrons
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       atomicAbsorptionCoefficient=0.0d0
       do i=1,self%countElements
          do j=0,self%elementAtomicNumbers(i)-1 ! j=0 is neutral atom; j=1 is first ionized state, etc.
             ! Determine the maximum sub-shell occupied by electrons. Sub-shell number is 1s=1, 2s=2, 2p=3, 3s=4, etc.
             countElectrons=0
             n             =0
             m             =0
             do while (countElectrons < self%elementAtomicNumbers(i)-j)
                n=n+1
                l=-1
                do while (l < n-1)
                   l             =l             +1
                   m             =m             +1
                   countElectrons=countElectrons+4*l+2
                   if (countElectrons >= self%elementAtomicNumbers(i)-j) exit
                end do
             end do
             crossSectionPhotoIonization=0.0d0
             do k=1,m
                crossSectionPhotoIonization=+crossSectionPhotoIonization                                                                                  &
                     &                      +self%atomicCrossSectionIonizationPhoto_%crossSection(                                                        &
                     &                                                                            atomicNumber   =self        %elementAtomicNumbers(i  ), &
                     &                                                                            ionizationState=                                  j+1 , &
                     &                                                                            shellNumber    =                                  k   , &
                     &                                                                            wavelength     =photonPacket%wavelength          (   )  &
                     &                                                                           )
             end do
             atomicAbsorptionCoefficient=+atomicAbsorptionCoefficient                       &
                  &                      +crossSectionPhotoIonization                       &
                  &                      *properties%elements(i)%densityNumber              &
                  &                      *properties%elements(i)%ionizationStateFraction(j) &
                  &                      /centi                                             &
                  &                      *megaParsec
          end do
       end do
    class default
       atomicAbsorptionCoefficient=0.0d0
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicAbsorptionCoefficient

  subroutine atomicAccumulatePhotonPacket(self,properties,photonPacket,absorptionCoefficient,lengthTraversed)
    !% Accumulate a photon packet.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : luminositySolar        , megaParsec
    use :: Numerical_Constants_Physical    , only : plancksConstant        , speedLight
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : angstromsPerMeter      , electronVolt
    implicit none
    class           (radiativeTransferMatterAtomic     ), intent(inout) :: self
    class           (radiativeTransferMatterProperties ), intent(inout) :: properties
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    double precision                                    , intent(in   ) :: absorptionCoefficient      , lengthTraversed
    double precision                                                    :: energyPhoton               , rateIonization             , &
         &                                                                 crossSectionPhotoIonization, atomicAbsorptionCoefficient
    integer                                                             :: i                          , j                          , &
         &                                                                 k                          , n                          , &
         &                                                                 l                          , m                          , &
         &                                                                 countElectrons
    !GCC$ attributes unused :: absorptionCoefficient

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Accumulate the photoionization rate per unit volume (cm⁻³), and photoheating rate per unit volume (J cm⁻³) using the Lucy
       ! (1999; A&A; 344; 282; https://ui.adsabs.harvard.edu/abs/1999A%26A...344..282L; see section 3.4) methodology.
       energyPhoton=+plancksConstant               &
            &       *speedLight                    &
            &       *angstromsPerMeter             &
            &       /photonPacket     %wavelength()
       do i=1,self%countElements
          do j=0,self%elementAtomicNumbers(i)-1 ! j=0 is neutral atom; j=1 is first ionized state, etc.
             ! Determine the maximum sub-shell occupied by electrons. Sub-shell number is 1s=1, 2s=2, 2p=3, 3s=4, etc.
             countElectrons=0
             n             =0
             m             =0
             do while (countElectrons < self%elementAtomicNumbers(i)-j)
                n=n+1
                l=-1
                do while (l < n-1)
                   l             =l             +1
                   m             =m             +1
                   countElectrons=countElectrons+4*l+2
                   if (countElectrons >= self%elementAtomicNumbers(i)-j) exit
                end do
             end do
             crossSectionPhotoIonization=0.0d0
             do k=1,m
                crossSectionPhotoIonization=+crossSectionPhotoIonization                                                                                  &
                     &                      +self%atomicCrossSectionIonizationPhoto_%crossSection(                                                        &
                     &                                                                            atomicNumber   =self        %elementAtomicNumbers(i  ), &
                     &                                                                            ionizationState=                                  j+1 , &
                     &                                                                            shellNumber    =                                  k   , &
                     &                                                                            wavelength     =photonPacket%wavelength          (   )  &
                     &                                                                           )
             end do
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
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicAccumulatePhotonPacket

  logical function atomicInteractWithPhotonPacket(self,properties,photonPacket)
    !% Interact with a photon packet. In this case the photon packet is always absorbed.
    implicit none
    class(radiativeTransferMatterAtomic     ), intent(inout) :: self
    class(radiativeTransferMatterProperties ), intent(inout) :: properties
    class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    !GCC$ attributes unused :: self, properties, photonPacket
    
    atomicInteractWithPhotonPacket=.false.
    return
  end function atomicInteractWithPhotonPacket

#ifdef USEMPI
  subroutine atomicAccumulationReduction(self,properties)
    !% Perform reduction of accumulated properaties across MPI processes.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: MPI_Utilities   , only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout)  :: self
    class  (radiativeTransferMatterProperties), intent(inout)  :: properties
    integer                                                    :: i
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       do i=1,self%countElements
          properties%elements(i)%photoIonizationRate=mpiSelf%sum(properties%elements(i)%photoIonizationRate)
          properties%elements(i)%photoHeatingRate   =mpiSelf%sum(properties%elements(i)%photoHeatingRate   )
       end do
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicAccumulationReduction
#endif
  
  subroutine atomicStateSolve(self,properties,status)
    !% Solve for the state of the matter.
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseA       , recombinationCaseB
    use :: Galacticus_Display                  , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    use :: Galacticus_Error                    , only : Galacticus_Error_Report  , errorStatusSuccess         , errorStatusFail
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Constants_Physical        , only : boltzmannsConstant       , speedLight                 , electronMass              , fineStructure    , &
         &                                              electronRadius
    use :: Numerical_Constants_Prefixes        , only : centi
    use :: Numerical_Roman_Numerals            , only : Roman_Numerals
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout)               :: self
    class           (radiativeTransferMatterProperties), intent(inout)               :: properties
    integer                                            , intent(  out) , optional    :: status
    double precision                                   , dimension(2)                :: temperaturePrevious                            , densityElectronsPrevious
    double precision                                   , parameter                   :: ionizationStateFractionToleranceAbsolute=1.0d-6, ionizationStateFractionToleranceRelative=1.0d-3, &
         &                                                                              ionizationStateFractionStepRelative     =5.0d-1, temperatureToleranceRelative            =1.0d-2, &
         &                                                                              temperatureStepRelative                 =5.0d-1, temperatureToleranceAbsolute            =1.0d+0
    integer                                            , parameter                   :: countIterationMaximum                   =100000, temperatureOscillatingCounts            =30    , &
         &                                                                              electronsOscillatingCounts              =30
    double precision                                   , parameter                   :: degreesOfFreedom                        =3.0d+0
    type            (element                          ), dimension( : ), allocatable :: elementsPrevious                               , elementsReference                              , &
         &                                                                              elementsPhotoRate
    double precision                                                                 :: rateRecombinationRadiative                     , rateRecombinationDielectronic                  , &
         &                                                                              densityNumberTotal                             , rateImbalanceFactor                            , &
         &                                                                              densityNumberElectrons                         ,   temperatureReference                         , &
         &                                                                              rateCoolingRecombinationRadiative              , rateCoolingRecombinationDielectronic           , &
         &                                                                              rateCoolingCollisionalExcitation               , rateCoolingBremsstrahlung                      , &
         &                                                                              rateHeating                                    , ratePhotoHeating                               , &
         &                                                                              energyTemperatureRatio                         , timeStepTemperature                            , &
         &                                                                              temperatureChangePrevious                      , chargeIonic                                    , &
         &                                                                              rateUpward                                     , rateDownward
    logical                                                                          :: converged                                      , temperatureOscillating                         , &
         &                                                                              electronsOscillating
    integer                                                                          :: countIteration                                 , i                                              , &
         &                                                                              temperatureOscillatingCount                    , j                                              , &
         &                                                                              recombinationCase                              , electronsOscillatingCount
    character       (len=256                         )                               :: message
    
    if (present(status)) status=errorStatusSuccess
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
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
       allocate(elementsPrevious (self%countElements))
       allocate(elementsReference(self%countElements))
       allocate(elementsPhotoRate(self%countElements))
       elementsPrevious =properties%elements
       elementsReference=properties%elements
       elementsPhotoRate=properties%elements
       ! Iterate until convergence is achieved.
       converged                       =all(properties%elements(:)%densityNumber <= 0.0d0)
       countIteration                  =0
       temperaturePrevious             =properties%temperature
       temperatureReference            =properties%temperature
       densityElectronsPrevious        =0.0d0
       temperatureChangePrevious       =0.0d0
       temperatureOscillating          =.false.
       electronsOscillating            =.false.
       temperatureOscillatingCount     =0
       electronsOscillatingCount       =0
       do i=1,self%countElements
          elementsPrevious (i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          elementsReference(i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          if (.not.converged) then
             elementsPrevious(i)%photoHeatingRate   =+properties%elements(i)%photoHeatingRate   
             elementsPrevious(i)%photoIonizationRate=+properties%elements(i)%photoIonizationRate &
                  &                                  /properties%elements(i)%densityNumber
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
          if (electronsOscillatingCount > 0) then
             densityNumberElectrons=0.5d0*(densityNumberElectrons+densityElectronsPrevious(1))
             electronsOscillatingCount=electronsOscillatingCount-1
          end if
          densityNumberTotal=+sum(properties%elements%densityNumber         ) &
               &             +                        densityNumberElectrons
          ! Compute rates of change of the ionization states and temperature.
          !! Initialize rates and computed photoionization and photoheating rates.
          ratePhotoHeating=0.0d0
          do i=1,self%countElements
             do j=0,self%elementAtomicNumbers(i)-1
                ! Compute the photoionization rate.
                elementsPhotoRate(i)%photoIonizationRate    (j)=+properties%elements(i)%photoIonizationRate(j) &
                     &                                          /properties%elements(i)%densityNumber
                ! Compute heating/cooling rates (in W cm⁻³).
                elementsPhotoRate(i)%photoHeatingRate       (j)=+properties%elements(i)%photoHeatingRate   (j)
                ! Scale the photoionization and photoheating rates with the lower-level state fraction to better approximate the
                ! correct solution in the optically thin limit. If the lower-level state fraction becomes zero, then we instead
                ! set the photoionization rate to half of the previous value (otherwise the gas will immediately become fully
                ! lower-level again).
                if    (           elementsReference(i)%ionizationStateFraction(j) > 0.0d0) then
                   if (properties%elements         (i)%ionizationStateFraction(j) > 0.0d0) then
                      elementsPhotoRate(i)%photoIonizationRate(j)=+           elementsPhotoRate(i)%photoIonizationRate    (j) &
                           &                                      *properties%elements         (i)%ionizationStateFraction(j) &
                           &                                      /           elementsReference(i)%ionizationStateFraction(j)
                      elementsPhotoRate(i)%photoHeatingRate   (j)=+           elementsPhotoRate(i)%photoHeatingRate       (j) &
                           &                                      *properties%elements         (i)%ionizationStateFraction(j) &
                           &                                      /           elementsReference(i)%ionizationStateFraction(j) 
                   else
                      elementsPhotoRate(i)%photoIonizationRate(j)=+           elementsPrevious (i)%photoIonizationRate    (j) &
                           &                                      /2.0d0
                      elementsPhotoRate(i)%photoHeatingRate   (j)=+           elementsPrevious (i)%photoHeatingRate       (j) &
                           &                                      /2.0d0
                   end if
                   elementsPrevious(i)%photoIonizationRate(j)=elementsPhotoRate(i)%photoIonizationRate(j)
                   elementsPrevious(i)%photoHeatingRate   (j)=elementsPhotoRate(i)%photoHeatingRate   (j)
                end if
                ! Accumulate the photo-heating rate.
                ratePhotoHeating=+ratePhotoHeating                         &
                     &           +elementsPhotoRate(i)%photoHeatingRate(j)
             end do
             ! For the fully-ionized state, use a photoionization rate from the lower state - this is just used for scaling timescales.
             elementsPhotoRate(i)%photoIonizationRate(self%elementAtomicNumbers(i))=elementsPhotoRate(i)%photoIonizationRate(self%elementAtomicNumbers(i)-1)
             elementsPhotoRate(i)%photoHeatingRate   (self%elementAtomicNumbers(i))=elementsPhotoRate(i)%photoHeatingRate   (self%elementAtomicNumbers(i)-1)
          end do
          ! Initialize heating/cooling rates to zero.
          rateCoolingRecombinationRadiative=0.0d0
          rateCoolingCollisionalExcitation=0.0d0
          rateCoolingBremsstrahlung=0.0d0
          ! Iterate over all states of all elements, accumulating heating/cooling rates and computing the equilibrium ionization fractions.
          do i=1,self%countElements
             do j=0,self%elementAtomicNumbers(i)
                !! Accumulate rate of collisional excitation cooling.
                if  (j < self%elementAtomicNumbers(i))                                                                                                                                   &
                     & rateCoolingCollisionalExcitation=+rateCoolingCollisionalExcitation                                                                                                &
                     &                                  +densityNumberElectrons                                                                                                          &
                     &                                  *properties%elements                        (i)%densityNumber                                                                    &
                     &                                  *properties%elements                        (i)%ionizationStateFraction(                             j                         ) &
                     &                                  *self      %atomicExcitationRateCollisional_   %coolingRate            (self%elementAtomicNumbers(i),j+1,properties%temperature) &
                     &                                  /centi**3
                !! Accumulate rate of Bremsstrahlung cooling.
                if (j > 0) then
                   chargeIonic              =dble(j)
                   rateCoolingBremsstrahlung=+rateCoolingBremsstrahlung                                                                                     &
                        &                    +16.0d0                                                                                                        &
                        &                    / 3.0d0                                                                                                        &
                        &                    *sqrt(                                                                                                         &
                        &                          +2.0d0                                                                                                   &
                        &                          *Pi                                                                                                      &
                        &                          /3.0d0                                                                                                   &
                        &                         )                                                                                                         &
                        &                    *chargeIonic**2                                                                                                &
                        &                    *properties%elements(i)%densityNumber                                                                          &
                        &                    *properties%elements(i)%ionizationStateFraction      (                             j                         ) &
                        &                    *densityNumberElectrons                                                                                        &
                        &                    *electronRadius          **3                                                                                   &
                        &                    *speedLight                                                                                                    &
                        &                    /electronRadius                                                                                                &
                        &                    *sqrt(                                                                                                         &
                        &                          +electronMass                                                                                            &
                        &                          *speedLight        **2                                                                                   &
                        &                          *boltzmannsConstant                                                                                      &
                        &                          *properties%temperature                                                                                  &
                        &                         )                                                                                                         &
                        &                    *fineStructure                                                                                                 &
                        &                    *self                  %gauntFactor_           %total(self%elementAtomicNumbers(i),j,properties%temperature)   &
                        &                    /centi**3
                end if
                if (j == 0) then
                   ! Set the neutral ion fraction to an arbitrary value - we will compute ratios of ionization fractions and
                   ! renormalize at the end.
                   properties%elements(i)%ionizationStateFraction(0)=1.0d0
                else
                   ! Compute the abundance of this ionization state relative to the lower state by requiring balance between
                   ! ionization and recombination rates.
                   !! Rate of upward transitions from photoionization.
                   if    (elementsReference(i)%ionizationStateFraction(j-1) > 0.0d0) then
                      rateUpward=+properties%elements         (i)%photoIonizationRate    (j-1) &
                           &     /properties%elements         (i)%densityNumber                &
                           &     /           elementsReference(i)%ionizationStateFraction(j-1)
                   else
                      rateUpward=+0.0d0
                   end if
                   ! Rate of upward transitions from collisional ionizations.
                   rateUpward                          =+                                          rateUpward                                                                                             &
                        &                               +                                          densityNumberElectrons                                                                                 &
                        &                               *self%atomicIonizationRateCollisional_    %rate                   (self%elementAtomicNumbers(i),j,properties%temperature                        )
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
                   rateRecombinationRadiative          =+                                          densityNumberElectrons                                                                                 &
                        &                               *self%atomicRecombinationRateRadiative_   %rate                   (self%elementAtomicNumbers(i),j,properties%temperature,level=recombinationCase)
                   rateRecombinationDielectronic       =+                                          densityNumberElectrons                                                                                 &
                        &                               *self%atomicRecombinationRateDielectronic_%rate                   (self%elementAtomicNumbers(i),j,properties%temperature                        )
                   !! Accumulate the rates of cooling due to recombinations.
                   rateCoolingRecombinationRadiative   =+rateCoolingRecombinationRadiative                    &
                        &                               +rateRecombinationRadiative                           &
                        &                               *0.75d0                                               & !! TO DO use the correct factors.
                        &                               *boltzmannsConstant                                   &
                        &                               *properties%temperature                               &
                        &                               *properties%elements   (i)%ionizationStateFraction(j) 
                   rateCoolingRecombinationDielectronic=+0.0d0 !! TO DO - get dielectronic cooling rates.
                   !! Find the net rate of downward transitions.
                   rateDownward                        =+rateRecombinationRadiative                           &
                        &                               +rateRecombinationDielectronic
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
                      properties%elements(i)%ionizationStateFraction   (j)=+properties%elements(i)%ionizationStateFraction(j-1) &
                           &                                               *rateUpward                                          &
                           &                                               /rateDownward
                   end if
                end if
             end do
             ! Renormalize ionization states so that the summed fraction is 1.
             properties%elements(i)%ionizationStateFraction=+    properties%elements(i)%ionizationStateFraction  &
                  &                                         /sum(properties%elements(i)%ionizationStateFraction)
          end do
          ! Sum cooling and heating rates.
          rateHeating=+ratePhotoHeating                     &
               &      -rateCoolingRecombinationRadiative    &
               &      -rateCoolingRecombinationDielectronic &
               &      -rateCoolingCollisionalExcitation     &
               &      -rateCoolingBremsstrahlung
          ! Solve for temperature.         
          if (ratePhotoHeating > 0.0d0) then
             ! Compute timestep for temperature.
             energyTemperatureRatio=+0.5d0*degreesOfFreedom*densityNumberTotal*boltzmannsConstant
             if (abs(rateHeating) > 0.0d0) then
                rateImbalanceFactor=min(1.0d0,abs(rateHeating)/ratePhotoHeating)
                timeStepTemperature=+               temperatureStepRelative &
                     &              *               rateImbalanceFactor     &
                     &              *               energyTemperatureRatio  &
                     &              *    properties%temperature             &
                     &              /abs(           rateHeating           )
                ! If temperature was oscillating, then limit the timestep such that the temperature can not change by more than
                ! half of the change in the previous timestep - this helps to avoid further oscillation.
                if (temperatureOscillatingCount > 0 .and. temperatureChangePrevious > 0.0d0) &
                     & timeStepTemperature=min(timeStepTemperature,abs(0.5d0*temperatureChangePrevious/(rateHeating/energyTemperatureRatio)))
                ! Adjust the temperature.
                properties%temperature=+properties%temperature &
                     &                 +timeStepTemperature    &
                     &                 *rateHeating            &
                     &                 /energyTemperatureRatio
                if (properties%temperature < 0.0d0) properties%temperature=self%temperatureMinimum        
             end if
             ! If temperature is oscillating decrease the count of steps for which we attempted to break out of oscillation.
             if (temperatureOscillatingCount > 0) temperatureOscillatingCount=+temperatureOscillatingCount &
                  &                                                           -1
          else
             ! No photo-heating - set to minimum temperature.
             properties%temperature                                    =self%temperatureMinimum
          end if
          ! Check for convergence.
          converged=   abs(properties%temperature-temperaturePrevious(1)) <  max(temperatureToleranceAbsolute,temperatureToleranceRelative*properties%temperature     ) &
               &    .and.                                                                                                                                               &
               &     (                                                                                                                                                  &
               &       abs(rateHeating                                  ) <=                                  temperatureToleranceRelative*           ratePhotoHeating  &
               &      .or.                                                                                                                                              &
               &       abs(properties%temperature-temperaturePrevious(1)) <      temperatureToleranceAbsolute                                                           &
               &     )
          do i=1,self%countElements
             if (.not.converged) exit
             converged=all(abs(properties%elements(i)%ionizationStateFraction-elementsPrevious(i)%ionizationStateFraction) < max(ionizationStateFractionToleranceAbsolute,ionizationStateFractionToleranceRelative*properties%elements(i)%ionizationStateFraction))
          end do
          ! Check for exceeding maximum iterations.
          countIteration=countIteration+1
          if (countIteration > countIterationMaximum .and. .not.converged) then
             ! No solution was found - report on the state of this domain cell.
             call Galacticus_Display_Indent('failed domain cell state report',verbosityStandard)
             write (message,'(a,e23.16         )') 'volume                                    = ',properties%volume
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16         )') 'initial temperature                       = ',temperatureReference
             call Galacticus_Display_Message(message,verbosityStandard)
             do i=1,self%countElements
                call Galacticus_Display_Indent('element: '//trim(adjustl(self%elements(i))),verbosityStandard)
                write (message,'(a,e23.16         )') 'density                                   = ',properties%elements(i)%densityNumber
                call Galacticus_Display_Message(message,verbosityStandard)
                do j=0,self%elementAtomicNumbers(i)
                   call Galacticus_Display_Indent('ion: '//trim(adjustl(self%elements(i)))//Roman_Numerals(j+1),verbosityStandard)
                   if (j < self%elementAtomicNumbers(i)) then
                      write (message,'(a,e23.16,a,e23.16)') 'photoionization rate (current : previous) = ',properties%elements         (i)%photoIonizationRate    (j),' : ',properties%elements(i)%photoIonizationRatePrevious(j)
                      call Galacticus_Display_Message(message,verbosityStandard)
                      write (message,'(a,e23.16,a,e23.16)') 'photoheating rate    (current : previous) = ',properties%elements         (i)%photoHeatingRate       (j),' : ',properties%elements(i)%photoHeatingRatePrevious   (j)
                      call Galacticus_Display_Message(message,verbosityStandard)
                   end if
                   write    (message,'(a,e23.16,a,e23.16)') 'initial ionization state                  = ',           elementsReference(i)%ionizationStateFraction(j)
                   call Galacticus_Display_Message(message,verbosityStandard)
                   call Galacticus_Display_Unindent('done',verbosityStandard)
                end do
                call Galacticus_Display_Unindent('done',verbosityStandard)
             end do
#ifdef RADTRANSDEBUG
             ! Debugging output. Writes the initial state for this failed case in a format that be directly pasted into the
             ! relevant test suite code.
             write    (0,'(a,i1)'         ) "properties            %iterationCount             = "    ,1
             write    (0,'(a,e23.16)'     ) "properties            %volume                     = "    ,properties%volume
             write    (0,'(a,e23.16)'     ) "properties            %temperature                = "    ,temperatureReference
             do i=1,self%countElements
                write (0,'(a,i1,a,e23.16)') "properties%elements(",i,")%densityNumber              = ",properties%elements(i)%densityNumber     
                write (0,'(a,i1,a,$)'     ) "properties%elements(",i,")%ionizationStateFraction    =["
                do j=0,self%elementAtomicNumbers(i)
                   if (j < self%elementAtomicNumbers(i)) then
                      write (0,'(e23.16,a,$)') elementsReference(i)%ionizationStateFraction    (j),','
                   else
                      write (0,'(e23.16,a)'  ) elementsReference(i)%ionizationStateFraction    (j),']'
                   end if
                end do
                write (0,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoIonizationRate        =["
                do j=0,self%elementAtomicNumbers(i)-1
                   if (j < self%elementAtomicNumbers(i)-1) then
                      write (0,'(e23.16,a,$)') elementsReference(i)%photoIonizationRate        (j),','
                   else
                      write (0,'(e23.16,a)'  ) elementsReference(i)%photoIonizationRate        (j),']'
                   end if
                end do
                write (0,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoHeatingRate           =["
                do j=0,self%elementAtomicNumbers(i)-1
                   if (j < self%elementAtomicNumbers(i)-1) then
                      write (0,'(e23.16,a,$)') elementsReference(i)%photoHeatingRate           (j),','
                   else
                      write (0,'(e23.16,a)'  ) elementsReference(i)%photoHeatingRate           (j),']'
                   end if
                end do
                write (0,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoIonizationRatePrevious=["
                do j=0,self%elementAtomicNumbers(i)-1
                   if (j < self%elementAtomicNumbers(i)-1) then
                      write (0,'(e23.16,a,$)') elementsReference(i)%photoIonizationRatePrevious(j),','
                   else
                      write (0,'(e23.16,a)'  ) elementsReference(i)%photoIonizationRatePrevious(j),']'
                   end if
                end do
                write (0,'(a,i1,a,$)'     ) "properties%elements(",i,")%photoHeatingRatePrevious   =["
                do j=0,self%elementAtomicNumbers(i)-1
                   if (j < self%elementAtomicNumbers(i)-1) then
                      write (0,'(e23.16,a,$)') elementsReference(i)%photoHeatingRatePrevious   (j),','
                   else
                      write (0,'(e23.16,a)'  ) elementsReference(i)%photoHeatingRatePrevious   (j),']'
                   end if
                end do
             end do
#endif
             call Galacticus_Display_Unindent('done',verbosityStandard)
             if (present(status)) then
                status=errorStatusFail
                return
             else
                call Galacticus_Error_Report('solution not found'//{introspection:location})
             end if
          end if
          ! Check for oscillating temperature or electron density.
          if (countIteration >= 3) then
             temperatureOscillating=(properties%temperature-temperaturePrevious     (1))*(temperaturePrevious     (1)-temperaturePrevious     (2)) < 0.0d0
             electronsOscillating  =(densityNumberElectrons-densityElectronsPrevious(1))*(densityElectronsPrevious(1)-densityElectronsPrevious(2)) < 0.0d0
             ! If oscillation is detected set the number of steps for which we will attempt to break out of it.
             if (temperatureOscillating) temperatureOscillatingCount=temperatureOscillatingCounts
             if (electronsOscillating  ) electronsOscillatingCount  =electronsOscillatingCounts
          end if
          ! Store current state to previous state.
          do i=1,self%countElements
             elementsPrevious(i)%ionizationStateFraction=properties%elements(i)%ionizationStateFraction
          end do
          temperatureChangePrevious        =abs(properties%temperature                 -temperaturePrevious     (1))
          temperaturePrevious           (2)=               temperaturePrevious     (1)
          temperaturePrevious           (1)=    properties%temperature
          densityElectronsPrevious      (2)=               densityElectronsPrevious(1)
          densityElectronsPrevious      (1)=               densityNumberElectrons
       end do
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicStateSolve

#ifdef USEMPI
  subroutine atomicBroadcastState(self,sendFromProcess,properties)
    !% Broadcast computational domain state to other MPI processes.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: MPI_Utilities   , only : mpiSelf
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    integer                                   , intent(in   ) :: sendFromProcess
    class  (radiativeTransferMatterProperties), intent(inout) :: properties
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       do i=1,self%countElements
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoIonizationRate        )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoIonizationRatePrevious)
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoHeatingRate           )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%photoHeatingRatePrevious   )
          call mpiSelf%broadcastData(sendFromProcess,properties%elements(i)%ionizationStateFraction    )
       end do
       call    mpiSelf%broadcastData(sendFromProcess,properties            %temperature                )
       class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicBroadcastState
#endif
  
  double precision function atomicConvergenceMeasure(self,properties)
    !% Return a convergence measure for the atomic matter.
    use :: Disparity_Ratios, only : Disparity_Ratio
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferMatterProperties), intent(inout) :: properties
    integer                                                   :: i
    !GCC$ attributes unused :: self
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       atomicConvergenceMeasure=-huge(0.0d0)
       do i=1,self%countElements
          atomicConvergenceMeasure=max(                                                                                                                          &
               &                           atomicConvergenceMeasure                                                                                            , &
               &                       max(                                                                                                                      &
               &                           Disparity_Ratio(properties%elements(i)%photoIonizationRate(0),properties%elements(i)%photoIonizationRatePrevious(0)), &
               &                           Disparity_Ratio(properties%elements(i)%photoHeatingRate   (0),properties%elements(i)%photoHeatingRatePrevious   (0))  &
               &                          )                                                                                                                      &
               &                      )
       end do
    class default
       atomicConvergenceMeasure=0.0d0
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicConvergenceMeasure
  
  double precision function atomicOutputProperty(self,properties,output)
    !% Return a scalar property to be output.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (radiativeTransferMatterAtomic    ), intent(inout) :: self
    class  (radiativeTransferMatterProperties), intent(inout) :: properties
    integer(c_size_t                         ), intent(in   ) :: output
    integer(c_size_t                         )                :: output_
    integer                                                   :: i
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       if (output == 1_c_size_t) then
          atomicOutputProperty=properties%temperature
       else if (output >= 2_c_size_t .and. output <= self%countOutputs_) then
          output_=output-1_c_size_t
          i      =1
          do while (output_ > self%elementAtomicNumbers(i)+1_c_size_t)
             output_=output_-self%elementAtomicNumbers(i)-1_c_size_t
             i      =i      +1
          end do
          select case (output_)
          case (1)
             atomicOutputProperty=properties%elements(i)%densityNumber
          case default
             atomicOutputProperty=properties%elements(i)%ionizationStateFraction(output_-1)
          end select
       else
          atomicOutputProperty=0.0d0
          call Galacticus_Error_Report('output is out of range'//{introspection:location})
       end if
    class default
       atomicOutputProperty=0.0d0
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end function atomicOutputProperty
  
  function atomicCountOutputs(self)
    !% Return the number of scalar properties to output.
    implicit none
    integer(c_size_t                     )                :: atomicCountOutputs
    class  (radiativeTransferMatterAtomic), intent(inout) :: self
    !GCC$ attributes unused :: self

    atomicCountOutputs=self%countOutputs_
    return
  end function atomicCountOutputs

  function atomicOutputName(self,output)
    !% Return the name of the scalar property to be output.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: ISO_Varying_String      , only : operator(//)
    use :: Numerical_Roman_Numerals, only : Roman_Numerals
    implicit none
    type   (varying_string               )                :: atomicOutputName
    class  (radiativeTransferMatterAtomic), intent(inout) :: self
    integer(c_size_t                     ), intent(in   ) :: output
    integer(c_size_t                     )                :: output_
    integer                                               :: i
    
    if (output == 1_c_size_t) then
       atomicOutputName=var_str('temperature'     )
    else if (output >= 2_c_size_t .and. output <= self%countOutputs_) then
       output_=output-1_c_size_t
       i      =1
       do while (output_ > self%elementAtomicNumbers(i)+1_c_size_t)
          output_=output_-self%elementAtomicNumbers(i)-1_c_size_t
          i      =i      +1
       end do
       select case (output_)
       case (1)
          atomicOutputName=var_str('densityNumber')//trim(adjustl(self%elements(i)))
       case default
          atomicOutputName=var_str('fraction'     )//trim(adjustl(self%elements(i)))//Roman_Numerals(int(output_))
       end select
    else
       atomicOutputName=var_str(''                )
       call Galacticus_Error_Report('output is out of range'//{introspection:location})
    end if
    return
  end function atomicOutputName

  double precision function atomicRecombinationRateHydrogen(self,properties)
    !% Return the total recombination rate for the atomic matter.
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical    , only : megaParsec
    use :: Numerical_Constants_Prefixes        , only : centi
    implicit none
    class           (radiativeTransferMatterAtomic          ), intent(inout) :: self
    type            (radiativeTransferMatterPropertiesAtomic), intent(inout) :: properties
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
       call Galacticus_Error_Report('hydrogren is not present'//{introspection:location})
    end if
    return
  end function atomicRecombinationRateHydrogen

