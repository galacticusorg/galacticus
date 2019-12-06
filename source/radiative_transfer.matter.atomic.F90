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
     class           (massDistributionClass                   ), pointer :: massDistribution_                    => null()
     class           (atomicCrossSectionIonizationPhotoClass  ), pointer :: atomicCrossSectionIonizationPhoto_   => null()
     class           (atomicRecombinationRateRadiativeClass   ), pointer :: atomicRecombinationRateRadiative_    => null()
     class           (atomicIonizationRateCollisionalClass    ), pointer :: atomicIonizationRateCollisional_     => null()
     class           (atomicRecombinationRateDielectronicClass), pointer :: atomicRecombinationRateDielectronic_ => null()
     class           (atomicIonizationPotentialClass          ), pointer :: atomicIonizationPotential_           => null()
     class           (atomicExcitationRateCollisionalClass    ), pointer :: atomicExcitationRateCollisional_     => null()
     class           (gauntFactorClass                        ), pointer :: gauntFactor_                         => null()
     integer                                                             :: indexAbundancePattern                         , iterationAverageCount
     double precision                                                    :: numberDensityMassDensityRatio                 , temperatureMinimum
   contains
     !@ <objectMethods>
     !@   <object>radiativeTransferMatterAtomic</object>
     !@   <objectMethod>
     !@     <method>recombinationRate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless type(radiativeTransferMatterPropertiesAtomic)\textgreater} properties</arguments>
     !@     <description>Return the total rate of recombinations (in units of s$^{-1}$).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                             atomicDestructor
     procedure :: propertyClass            => atomicPropertyClass
     procedure :: populateDomain           => atomicPopulateDomain
     procedure :: reset                    => atomicReset
     procedure :: absorptionCoefficient    => atomicAbsorptionCoefficient
     procedure :: accumulatePhotonPacket   => atomicAccumulatePhotonPacket
     procedure :: interactWithPhotonPacket => atomicInteractWithPhotonPacket
     procedure :: stateSolve               => atomicStateSolve
     procedure :: convergenceMeasure       => atomicConvergenceMeasure
     procedure :: outputProperty           => atomicOutputProperty
     procedure :: countOutputs             => atomicCountOutputs
     procedure :: outputName               => atomicOutputName
     procedure :: recombinationRate        => atomicRecombinationRate
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

  type, extends(radiativeTransferMatterProperties), public :: radiativeTransferMatterPropertiesAtomic
     !% Radiative transfer matter properties class for atomic matter.
     integer                                       :: iterationCount
     double precision                              :: densityNumber              , volume                  , &
          &                                           temperature
     double precision, dimension(0:1)              :: ionizationStateFraction 
     double precision, dimension(:,:), allocatable :: photoIonizationRateHistory , photoHeatingRateHistory
     double precision, dimension(0:0)              :: photoIonizationRate        , photoHeatingRate        , &
          &                                           photoIonizationRatePrevious, photoHeatingRatePrevious
  end type radiativeTransferMatterPropertiesAtomic
  
contains

  function atomicConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily atomic} radiative transfer matter class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferMatterAtomic           )                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (massDistributionClass                   ), pointer       :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass  ), pointer       :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass   ), pointer       :: atomicRecombinationRateRadiative_
    class           (atomicIonizationRateCollisionalClass    ), pointer       :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass), pointer       :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass          ), pointer       :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass    ), pointer       :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                        ), pointer       :: gauntFactor_
    integer                                                                   :: iterationAverageCount    
    double precision                                                          :: temperatureMinimum
    
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
    !# <objectBuilder class="massDistribution"                    name="massDistribution_"                    source="parameters"/>
    !# <objectBuilder class="atomicCrossSectionIonizationPhoto"   name="atomicCrossSectionIonizationPhoto_"   source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateRadiative"    name="atomicRecombinationRateRadiative_"    source="parameters"/>
    !# <objectBuilder class="atomicIonizationRateCollisional"     name="atomicIonizationRateCollisional_"     source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateDielectronic" name="atomicRecombinationRateDielectronic_" source="parameters"/>
    !# <objectBuilder class="atomicIonizationPotential"           name="atomicIonizationPotential_"           source="parameters"/>
    !# <objectBuilder class="atomicExcitationRateCollisional"     name="atomicExcitationRateCollisional_"     source="parameters"/>
    !# <objectBuilder class="gauntFactor"                         name="gauntFactor_"                         source="parameters"/>
    self=radiativeTransferMatterAtomic(iterationAverageCount,temperatureMinimum,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_)
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

  function atomicConstructorInternal(iterationAverageCount,temperatureMinimum,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_,atomicExcitationRateCollisional_,gauntFactor_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily atomic} radiative transfer matter class.
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup, Atomic_Abundance, Atomic_Mass
    use :: Numerical_Constants_Astronomical, only : massSolar               , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Prefixes    , only : centi
    implicit none
    type            (radiativeTransferMatterAtomic           )                        :: self
    integer                                                   , intent(in   )         :: iterationAverageCount
    double precision                                          , intent(in   )         :: temperatureMinimum
    class           (massDistributionClass                   ), intent(in   ), target :: massDistribution_
    class           (atomicCrossSectionIonizationPhotoClass  ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class           (atomicRecombinationRateRadiativeClass   ), intent(in   ), target :: atomicRecombinationRateRadiative_
    class           (atomicIonizationRateCollisionalClass    ), intent(in   ), target :: atomicIonizationRateCollisional_
    class           (atomicRecombinationRateDielectronicClass), intent(in   ), target :: atomicRecombinationRateDielectronic_
    class           (atomicIonizationPotentialClass          ), intent(in   ), target :: atomicIonizationPotential_
    class           (atomicExcitationRateCollisionalClass    ), intent(in   ), target :: atomicExcitationRateCollisional_
    class           (gauntFactorClass                        ), intent(in   ), target :: gauntFactor_
    !# <constructorAssign variables="iterationAverageCount, temperatureMinimum, *massDistribution_, *atomicCrossSectionIonizationPhoto_, *atomicRecombinationRateRadiative_, *atomicIonizationRateCollisional_, *atomicRecombinationRateDielectronic_, *atomicIonizationPotential_, *atomicExcitationRateCollisional_, *gauntFactor_"/>

    ! Get an abundance pattern and compute conversion factors from total gas-phase mass density to atomic number densities (in
    ! cm⁻³) for each element.    
    !! TO DO - for now we're just doing hydrogen.
    self%indexAbundancePattern        = Abundance_Pattern_Lookup(                           abundanceName="solar")
    self%numberDensityMassDensityRatio=+Atomic_Abundance        (self%indexAbundancePattern,shortLabel   ="H"    ) &
         &                             /Atomic_Mass             (                           shortLabel   ="H"    ) &
         &                             /atomicMassUnit                                                             &
         &                             *massSolar                                                                  &
         &                             /megaParsec    **3                                                          &
         &                             *centi         **3
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
    class  (radiativeTransferMatterAtomic           ), intent(inout) :: self
    class  (radiativeTransferMatterProperties       ), intent(inout) :: properties
    class  (computationalDomainVolumeIntegratorClass), intent(inout) :: integrator
    logical                                          , intent(in   ) :: onProcess

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Initialize properties:
       ! + Assume that the atomic species are fully neutral initially.
       ! + Assign an initial temperature.
       ! + Initialize photoionization and photo heating rates to impossible values. These will be used as the "previous" value on
       !   the first iteration, guaranteeing that our first iteration is never judged to be converged.
       allocate(properties%photoIonizationRateHistory(self%iterationAverageCount,0:0))
       allocate(properties%photoHeatingRateHistory   (self%iterationAverageCount,0:0))
       properties%iterationCount            =+0
       properties%ionizationStateFraction   =[1.0d0,0.0d0]
       properties%temperature               =1.0d4
       properties%photoIonizationRate       =-huge(0.0d0)
       properties%photoHeatingRate          =-huge(0.0d0)
       properties%photoIonizationRateHistory=-huge(0.0d0)
       properties%photoHeatingRateHistory   =-huge(0.0d0)
       properties%volume                    =+integrator%volume()
       ! Compute the number density of this atomic species.
       if (onProcess) then
          properties%densityNumber=+integrator%integrate                    (atomicDensityIntegrand) &
               &                   /properties%volume                                                &
               &                   *self      %numberDensityMassDensityRatio
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
       call mpiSelf%broadcastData(sendFromProcess,properties%densityNumber)
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
    class(radiativeTransferMatterAtomic    ), intent(inout) :: self
    class(radiativeTransferMatterProperties), intent(inout) :: properties
    !GCC$ attributes unused :: self
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Store the current photoionization and photoheating rates as the "previous" values, and then reset these rates to zero.
       properties%photoIonizationRatePrevious=properties%photoIonizationRate
       properties%photoHeatingRatePrevious   =properties%photoHeatingRate
       properties%photoIonizationRate        =0.0d0
       properties%photoHeatingRate           =0.0d0
       properties%iterationCount             =properties%iterationCount     +1
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
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       crossSectionPhotoIonization=self%atomicCrossSectionIonizationPhoto_%crossSection(                                           &
                                                                                        atomicNumber   =1                        , &
                                                                                        ionizationState=1                        , &
                                                                                        shellNumber    =1                        , &
                                                                                        wavelength     =photonPacket%wavelength()  &
                                                                                       ) !! TO DO - currently hard-coded for hydrogen
       atomicAbsorptionCoefficient=+crossSectionPhotoIonization           &
            &                      *properties%densityNumber              &
            &                      *properties%ionizationStateFraction(0) &
            &                      /centi                                 &
            &                      *megaParsec     
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
    double precision                                    , intent(in   ) :: absorptionCoefficient, lengthTraversed
    double precision                                                    :: energyPhoton         , rateIonization
    !GCC$ attributes unused :: self

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Accumulate the photoionization rate per unit volume (cm⁻³), and photoheating rate per unit volume (J cm⁻³) using the Lucy
       ! (1999; A&A; 344; 282; https://ui.adsabs.harvard.edu/abs/1999A%26A...344..282L; see section 3.4) methodology.
       energyPhoton                     =+plancksConstant                                                    &
            &                            *speedLight                                                         &
            &                            *angstromsPerMeter                                                  &
            &                            /photonPacket         %wavelength                            ( )
       rateIonization                   =+photonPacket         %luminosity                            ( )    &
            &                            *luminositySolar                                                    &
            &                            /energyPhoton                                                       &
            &                            *absorptionCoefficient                                              &
            &                            *lengthTraversed                                                    &
            &                            /properties           %volume                                       &
            &                            *centi                                                          **3 &
            &                            /megaParsec                                                     **3
       properties%photoIonizationRate(0)=+properties           %photoIonizationRate                 (0  )    &
            &                            +rateIonization
       properties%photoHeatingRate   (0)=+properties           %photoHeatingRate                    (0  )    &
            &                            +rateIonization                                                     &
            &                            *(                                                                  &
            &                              +energyPhoton                                                     &
            &                              -self               %atomicIonizationPotential_%potential(1,1)    & ! TO DO - hard coded for hydrogen
            &                              *electronVolt                                                     &
            &                             )
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
    class(radiativeTransferMatterAtomic    ), intent(inout)  :: self
    class(radiativeTransferMatterProperties), intent(inout)  :: properties
    !GCC$ attributes unused :: self
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       properties%photoIonizationRate=mpiSelf%sum(properties%photoIonizationRate)
       properties%photoHeatingRate   =mpiSelf%sum(properties%photoHeatingRate   )
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicAccumulationReduction
#endif
  
  subroutine atomicStateSolve(self,properties,status)
    !% Solve for the state of the matter.
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use :: Galacticus_Display                  , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    use :: Galacticus_Error                    , only : Galacticus_Error_Report  , errorStatusSuccess         , errorStatusFail
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Constants_Physical        , only : boltzmannsConstant       , speedLight                 , electronMass              , fineStructure    , &
         &                                              electronRadius
    use :: Numerical_Constants_Prefixes        , only : centi
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout)           :: self
    class           (radiativeTransferMatterProperties), intent(inout)           :: properties
    integer                                            , intent(  out), optional :: status
    double precision                                   , dimension(0:1)          :: rateIonizationState                            , ionizationStateFractionPrevious                , &
         &                                                                          ionizationStateFractionReference
    double precision                                   , dimension(2)            :: temperaturePrevious
    double precision                                   , parameter               :: ionizationStateFractionToleranceAbsolute=1.0d-6, ionizationStateFractionToleranceRelative=1.0d-3, &
         &                                                                          ionizationStateFractionStepRelative     =5.0d-1, temperatureToleranceRelative            =1.0d-2, &
         &                                                                          temperatureStepRelative                 =5.0d-1, temperatureToleranceAbsolute            =1.0d+0
    integer                                            , parameter               :: countIterationMaximum                   =10000 , temperatureOscillatingCounts            =30
    double precision                                   , parameter               :: degreesOfFreedom                        =3.0d+0
    double precision                                                             :: ratePhotoIonization                            , rateRecombination                              , &
         &                                                                          densityNumberTotal                             , rateImbalanceFactor                            , &
         &                                                                          densityNumberElectrons                         , rateNonPhotoIonization                         , &
         &                                                                          rateCoolingRecombinationRadiative              , rateCoolingRecombinationDielectronic           , &
         &                                                                          rateCoolingCollisionalExcitation               , rateCoolingBremsstrahlung                      , &
         &                                                                          ratePhotoHeating                               , rateHeating                                    , &
         &                                                                          energyTemperatureRatio                         , timeStepTemperature                            , &
         &                                                                          timeStepIonization                             , temperatureReference                           , &
         &                                                                          ratePhotoHeatingPrevious                       , ratePhotoIonizationPrevious                    , &
         &                                                                          temperatureChangePrevious
    logical                                                                      :: converged                                      , temperatureOscillating
    integer                                                                      :: countIteration                                 , i                                              , &
         &                                                                          temperatureOscillatingCount
    character       (len=256                         )                           :: message

    if (present(status)) status=errorStatusSuccess
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Append the accumulated photoionization/photoheating rates to the history arrays, compute the average over those histories, and replace the photoionization/photoheating rates with those averages.
       properties%photoIonizationRateHistory(2:self%iterationAverageCount,:)=properties%photoIonizationRateHistory(1:self%iterationAverageCount-1,:)
       properties%photoHeatingRateHistory   (2:self%iterationAverageCount,:)=properties%photoHeatingRateHistory   (1:self%iterationAverageCount-1,:)
       properties%photoIonizationRateHistory(1                           ,:)=properties%photoIonizationRate       (                               :)
       properties%photoHeatingRateHistory   (1                           ,:)=properties%photoHeatingRate          (                               :)
       do i=0,0
          properties%photoIonizationRate(i)=sum(properties%photoIonizationRateHistory(1:min(properties%iterationCount,self%iterationAverageCount),i))/dble(min(properties%iterationCount,self%iterationAverageCount))
          properties%photoHeatingRate   (i)=sum(properties%photoHeatingRateHistory   (1:min(properties%iterationCount,self%iterationAverageCount),i))/dble(min(properties%iterationCount,self%iterationAverageCount))
       end do
       ! Iterate until convergence is achieved.
       converged                       =properties%densityNumber <= 0.0d0
       countIteration                  =0
       ionizationStateFractionPrevious =properties%ionizationStateFraction
       ionizationStateFractionReference=properties%ionizationStateFraction
       temperaturePrevious             =properties%temperature
       temperatureReference            =properties%temperature
       temperatureChangePrevious       =0.0d0
       temperatureOscillating          =.false.
       temperatureOscillatingCount     =0
       if (.not.converged) then
          ratePhotoHeatingPrevious        =properties%photoHeatingRate   (0)
          ratePhotoIonizationPrevious     =properties%photoIonizationRate(0)/properties%densityNumber
       end if
       do while (.not.converged)
          ! Compute total electron density.
          densityNumberElectrons=+properties%densityNumber              &
               &                 *properties%ionizationStateFraction(1)
          ! Compute total particle density.
          densityNumberTotal    =+           densityNumberElectrons     &
               &                 +properties%densityNumber
          ! Compute rates of change of the ionization states.       
          rateIonizationState   =+0.0d0
          ratePhotoIonization   =+properties%photoIonizationRate                      (  0                                                ) &
               &                 /properties%densityNumber
          ! Scale the photoionization rate with the neutral fraction to better approximate the correct solution in the optically
          ! thin limit. If the neutral fraction becomes zero, then we instead set the photoionization rate to half of the previous value
          ! (otherwise the gas will immediately become fully neutral again).
          if (ionizationStateFractionReference(0) > 0.0d0) then
             if (properties%ionizationStateFraction(0) > 0.0d0) then
                ratePhotoIonization=+ratePhotoIonization                   &
                     &              *properties%ionizationStateFraction(0) &
                     &              /ionizationStateFractionReference  (0)
             else
                ratePhotoIonization=+ratePhotoIonizationPrevious           &
                     &              /2.0d0
             end if
             ratePhotoIonizationPrevious=ratePhotoIonization
          end if
          rateNonPhotoIonization=+           densityNumberElectrons                                                                         &
               &                 *properties%ionizationStateFraction                  (  0                                                ) &
               &                 *  self    %atomicIonizationRateCollisional_    %rate(1,1,properties%temperature                         )   !! TO DO - hard-coded for hydrogen, fixed temperature
          rateRecombination     =+           densityNumberElectrons                                                                         &
               &                 *properties%ionizationStateFraction                  (  1                                                ) &
               &                 *(                                                                                                         &
               &                   +self    %atomicRecombinationRateRadiative_   %rate(1,1,properties%temperature,level=recombinationCaseB) & !! TO DO - hard-coded for hydrogen, fixed temperature, and case-B recombination
               &                   +self    %atomicRecombinationRateDielectronic_%rate(1,1,properties%temperature                         ) & !! TO DO - hard-coded for hydrogen, fixed temperature
               &                  )
          rateIonizationState(0)=-ratePhotoIonization-rateNonPhotoIonization+rateRecombination
          rateIonizationState(1)=+ratePhotoIonization+rateNonPhotoIonization-rateRecombination
          ! Compute heating/cooling rates (in W cm⁻³).
          ratePhotoHeating                 =+properties%photoHeatingRate(0)
          ! Scale the photoheating rate with the neutral fraction to better approximate the correct solution in the optically thin
          ! limit. If the neutral fraction becomes zero, then we instead set the photoheating rate to half of the previous value
          ! (otherwise the temperature will drop to the minimum).
          if (ionizationStateFractionReference(0) > 0.0d0) then
             if (properties%ionizationStateFraction(0) > 0.0d0) then
                ratePhotoHeating=+ratePhotoHeating                      &
                     &           *properties%ionizationStateFraction(0) &
                     &           /ionizationStateFractionReference  (0)
             else
                ratePhotoHeating=+ratePhotoHeatingPrevious              &
                     &           /2.0d0
             end if
             ratePhotoHeatingPrevious=ratePhotoHeating
          end if
          rateCoolingRecombinationRadiative   =+rateRecombination              &
               &                               *0.75d0                         & !! TO DO use the correct factors (from Spizter or something more recent?) here
               &                               *boltzmannsConstant             &
               &                               *properties%temperature           !! TO DO temperature fixed at 10^4 K here.
          if (.false.) then !! TO DO - should implement helium cooling rate here
             rateCoolingRecombinationDielectronic=+0.0d0 !! TO DO - add in Helium dielectronic cooling rate here
          else
             rateCoolingRecombinationDielectronic=+0.0d0 !! TO DO - get dielectronic cooling rates for other elements
          end if
          rateCoolingCollisionalExcitation    =+           densityNumberElectrons                                                   &
               &                               *properties%densityNumber                                                            &
               &                               *properties%ionizationStateFraction                     (0                         ) &
               &                               *self      %atomicExcitationRateCollisional_%coolingRate(1,1,properties%temperature) &
               &                               /centi**3
          rateCoolingBremsstrahlung           =+16.0d0                                                                              &
               &                               / 3.0d0                                                                              &
               &                               *sqrt(                                                                               &
               &                                     +2.0d0                                                                         &
               &                                     *Pi                                                                            &
               &                                     /3.0d0                                                                         &
               &                                    )                                                                               &
               &                               *dble(1) **2                                                                         & !! Ionic charge
               &                               *properties%densityNumber                                                            &
               &                               *properties%ionizationStateFraction                     (1                         ) &
               &                               *densityNumberElectrons                                                              &
               &                               *electronRadius          **3                                                         &
               &                               *speedLight                                                                          &
               &                               /electronRadius                                                                      &
               &                               *sqrt(                                                                               &
               &                                     +electronMass                                                                  &
               &                                     *speedLight        **2                                                         &
               &                                     *boltzmannsConstant                                                            &
               &                                     *properties%temperature                                                        &
               &                                    )                                                                               &
               &                               *fineStructure                                                                       &
               &                               *self      %gauntFactor_                    %total      (1,1,properties%temperature) &
               &                               /centi**3
          rateHeating                         =+ratePhotoHeating                     &
               &                               -rateCoolingRecombinationRadiative    &
               &                               -rateCoolingRecombinationDielectronic &
               &                               -rateCoolingCollisionalExcitation     &
               &                               -rateCoolingBremsstrahlung
          ! Check for zero photoionization.
          if (ratePhotoionization <= 0.0d0) then
             properties%ionizationStateFraction(0  )=1.0d0
             properties%ionizationStateFraction(1:1)=0.0d0
             properties%temperature                 =self%temperatureMinimum
          else
             ! Compute the equilibrium ionization states and temperature.             
             !! Determine a time-step. We seek a solution where the net rate of change for each ionization state and temperature
             !! is zero. Therefore, we compute a "rate imbalance factor" which is the absolute value of the rate of change divided
             !! by the photoionization/photoheating rate (which is fixed). The larger the factor, the bigger the step we should take to
             !! approach equilibrium.
             !!! Limit timestep for temperature.
             energyTemperatureRatio=+0.5d0*degreesOfFreedom*densityNumberTotal*boltzmannsConstant
             timeStepTemperature=huge(0.0d0)
             if (abs(rateHeating) > 0.0d0) then
                rateImbalanceFactor=min(1.0d0,abs(rateHeating)/ratePhotoHeating)
                timeStepTemperature=temperatureStepRelative*rateImbalanceFactor*energyTemperatureRatio*properties%temperature/abs(rateHeating)
                ! If temperature was oscillating, then limit the timestep such that the temperature can not change by more than
                ! half of the change in the previous timestep - this helps to avoid further oscillation.
                if (temperatureOscillatingCount > 0 .and. temperatureChangePrevious > 0.0d0) &
                     & timeStepTemperature=min(timeStepTemperature,abs(0.5d0*temperatureChangePrevious/(rateHeating/energyTemperatureRatio)))
             end if
             ! If temperature is oscillating, force the timestep to be determined by temperature, and decrease the count of
             ! steps for which we attempted to break out of oscillation.
             if (temperatureOscillatingCount > 0) &
                  & temperatureOscillatingCount=temperatureOscillatingCount-1
             !!! Limit timestep for ionization states.
             timeStepIonization=huge(0.0d0)
             do i=0,1
                if (abs(rateIonizationState(i)) > 0.0d0 .and. properties%ionizationStateFraction(i) > 0.0d0) then
                   rateImbalanceFactor=min(1.0d0,abs(rateIonizationState(i))/ratePhotoIonization)
                   timeStepIonization=min(timeStepIonization,ionizationStateFractionStepRelative*rateImbalanceFactor*properties%ionizationStateFraction(i)/abs(rateIonizationState(i)))
                end if
             end do
             ! Adjust the temperature.
             properties%temperature=+properties%temperature &
                  &                 +timeStepTemperature    &
                  &                 *rateHeating            &
                  &                 /energyTemperatureRatio
             if (properties%temperature < 0.0d0) properties%temperature=self%temperatureMinimum
             ! Adjust ionization state.
             properties%ionizationStateFraction=properties%ionizationStateFraction+timeStepIonization*rateIonizationState
             if (any(properties%ionizationStateFraction < 0.0d0)) then
                where (properties%ionizationStateFraction < 0.0d0)
                   properties%ionizationStateFraction=0.0d0
                end where
                properties%ionizationStateFraction=properties%ionizationStateFraction/sum(properties%ionizationStateFraction)
             end if
          end if
          ! Check for convergence.
          converged= all(abs(properties%ionizationStateFraction-ionizationStateFractionPrevious) <  max(ionizationStateFractionToleranceAbsolute,ionizationStateFractionToleranceRelative*properties%ionizationStateFraction)) &
               &    .and.                                                                                                                                                                                                      &
               &     all(abs(rateIonizationState                                               ) <=                                              ionizationStateFractionToleranceRelative*           ratePhotoIonization     ) &
               &    .and.                                                                                                                                                                                                      &
               &         abs(properties%temperature            -temperaturePrevious(1)         ) <  max(temperatureToleranceAbsolute            ,temperatureToleranceRelative            *properties%temperature            )  &
               &    .and.                                                                                                                                                                                                      &
               &     (                                                                                                                                                                                                         &
               &         abs(rateHeating                                                       ) <=                                              temperatureToleranceRelative            *           ratePhotoHeating          &
               &      .or.                                                                                                                                                                                                     &
               &         abs(properties%temperature            -temperaturePrevious(1)         ) < temperatureToleranceAbsolute                                                                                                &
               &     )
          ! Check for exceeding maximum iterations.
          countIteration=countIteration+1
          if (countIteration > countIterationMaximum .and. .not.converged) then
             ! No solution was found - report on the state of this domain cell.
             call Galacticus_Display_Indent('failed domain cell state report',verbosityStandard)
             write (message,'(a,e23.16         )') 'volume                                    = ',properties%volume
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16         )') 'density                                   = ',properties%densityNumber
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16,a,e23.16)') 'photoionization rate (current : previous) = ',properties%photoIonizationRate     ,' : ',properties%photoIonizationRatePrevious
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16,a,e23.16)') 'photoheating rate    (current : previous) = ',properties%photoHeatingRate        ,' : ',properties%photoHeatingRatePrevious
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16,a,e23.16)') 'initial ionization state                  = ',ionizationStateFractionReference(0),';  ',ionizationStateFractionReference      (1)
             call Galacticus_Display_Message(message,verbosityStandard)
             write (message,'(a,e23.16         )') 'initial temperature                       = ',temperatureReference
             call Galacticus_Display_Message(message,verbosityStandard)
             call Galacticus_Display_Unindent('done',verbosityStandard)
             if (present(status)) then
                status=errorStatusFail
                return
             else
                call Galacticus_Error_Report('solution not found'//{introspection:location})
             end if
          end if
          ! Check for oscillating temperature.
          if (countIteration >= 3) then
             temperatureOscillating=(properties%temperature-temperaturePrevious(1))*(temperaturePrevious(1)-temperaturePrevious(2)) < 0.0d0
             ! If oscillation is detected set the number of steps for which we will attempt to break out of it.
             if (temperatureOscillating) temperatureOscillatingCount=temperatureOscillatingCounts
          end if
          ! Store current state to previous state.
          ionizationStateFractionPrevious=    properties%ionizationStateFraction
          temperatureChangePrevious      =abs(properties%temperature            -temperaturePrevious(1))
          temperaturePrevious(2)         =    temperaturePrevious(1)
          temperaturePrevious(1)         =    properties%temperature
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
    !GCC$ attributes unused :: self

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       call mpiSelf%broadcastData(sendFromProcess,properties%photoIonizationRate        )
       call mpiSelf%broadcastData(sendFromProcess,properties%photoIonizationRatePrevious)
       call mpiSelf%broadcastData(sendFromProcess,properties%photoHeatingRate           )
       call mpiSelf%broadcastData(sendFromProcess,properties%photoHeatingRatePrevious   )
       call mpiSelf%broadcastData(sendFromProcess,properties%ionizationStateFraction    )
       call mpiSelf%broadcastData(sendFromProcess,properties%temperature                )
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
    class(radiativeTransferMatterAtomic    ), intent(inout) :: self
    class(radiativeTransferMatterProperties), intent(inout) :: properties
    !GCC$ attributes unused :: self
    
    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       atomicConvergenceMeasure=max(                                                                                              &
            &                       Disparity_Ratio(properties%photoIonizationRate(0),properties%photoIonizationRatePrevious(0)), &
            &                       Disparity_Ratio(properties%photoHeatingRate   (0),properties%photoHeatingRatePrevious   (0))  &
            &                      )
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
    !GCC$ attributes unused :: self

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       select case (output)
       case (1_c_size_t)
          atomicOutputProperty=properties%temperature
       case (2_c_size_t)
          atomicOutputProperty=properties%densityNumber
       case (3_c_size_t)
          atomicOutputProperty=properties%ionizationStateFraction(1)
       case default
          atomicOutputProperty=0.0d0
          call Galacticus_Error_Report('output is out of range'//{introspection:location})
       end select
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

    atomicCountOutputs=3_c_size_t
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
    !GCC$ attributes unused :: self
    
    select case (output)
    case (1_c_size_t)
       atomicOutputName=var_str('temperature'     )
    case (2_c_size_t)
       atomicOutputName=var_str('densityNumber'   )
    case (3_c_size_t)
       atomicOutputName=var_str('fractionHydrogen')//Roman_Numerals(1+1)
    case default
       atomicOutputName=var_str(''                )
       call Galacticus_Error_Report('output is out of range'//{introspection:location})
    end select
    return
  end function atomicOutputName

  double precision function atomicRecombinationRate(self,properties)
    !% Return the total recombination rate for the atomic matter.
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use :: Numerical_Constants_Astronomical    , only : megaParsec
    use :: Numerical_Constants_Prefixes        , only : centi
    implicit none
    class           (radiativeTransferMatterAtomic          ), intent(inout) :: self
    type            (radiativeTransferMatterPropertiesAtomic), intent(inout) :: properties
    double precision                                                         :: densityNumberElectrons
    
    densityNumberElectrons =+properties%densityNumber                                                                               &
         &                  *properties%ionizationStateFraction               (  1)
    atomicRecombinationRate=+           densityNumberElectrons                                                                      &
         &                  *properties%densityNumber                                                                               &
         &                  *properties%ionizationStateFraction               (  1)                                                 &
         &                  *self      %atomicRecombinationRateRadiative_%rate(1,1,properties%temperature,level=recombinationCaseB) &
         &                  *properties%volume                                                                                      &
         &                  *(                                                                                                      &
         &                    +megaParsec                                                                                           &
         &                    /centi                                                                                                &
         &                   )**3
    return
  end function atomicRecombinationRate

