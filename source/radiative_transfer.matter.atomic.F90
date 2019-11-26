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
     integer                                                             :: indexAbundancePattern                         , iterationAverageCount
     double precision                                                    :: numberDensityMassDensityRatio
   contains
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
  end type radiativeTransferMatterAtomic
  
  interface radiativeTransferMatterAtomic
     !% Constructors for the {\normalfont \ttfamily atomic} radiative transfer matter class.
     module procedure atomicConstructorParameters
     module procedure atomicConstructorInternal
  end interface radiativeTransferMatterAtomic

  type, extends(radiativeTransferMatterProperties) :: radiativeTransferMatterPropertiesAtomic
     !% Radiative transfer matter properties class for atomic matter.
     integer                                       :: iterationCount
     double precision                              :: densityNumber             , volume                   , &
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
    type   (radiativeTransferMatterAtomic           )                :: self
    type   (inputParameters                         ), intent(inout) :: parameters
    class  (massDistributionClass                   ), pointer       :: massDistribution_
    class  (atomicCrossSectionIonizationPhotoClass  ), pointer       :: atomicCrossSectionIonizationPhoto_
    class  (atomicRecombinationRateRadiativeClass   ), pointer       :: atomicRecombinationRateRadiative_
    class  (atomicIonizationRateCollisionalClass    ), pointer       :: atomicIonizationRateCollisional_
    class  (atomicRecombinationRateDielectronicClass), pointer       :: atomicRecombinationRateDielectronic_
    class  (atomicIonizationPotentialClass          ), pointer       :: atomicIonizationPotential_
    integer                                                          :: iterationAverageCount    

    !# <inputParameter>
    !#   <name>iterationAverageCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>5</defaultValue>
    !#   <description>The number of iterations over which to average the photoionization rate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="massDistribution"                    name="massDistribution_"                    source="parameters"/>
    !# <objectBuilder class="atomicCrossSectionIonizationPhoto"   name="atomicCrossSectionIonizationPhoto_"   source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateRadiative"    name="atomicRecombinationRateRadiative_"    source="parameters"/>
    !# <objectBuilder class="atomicIonizationRateCollisional"     name="atomicIonizationRateCollisional_"     source="parameters"/>
    !# <objectBuilder class="atomicRecombinationRateDielectronic" name="atomicRecombinationRateDielectronic_" source="parameters"/>
    !# <objectBuilder class="atomicIonizationPotential"           name="atomicIonizationPotential_"           source="parameters"/>
    self=radiativeTransferMatterAtomic(iterationAverageCount,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_)
    !# <objectDestructor name="massDistribution_"                   />
    !# <objectDestructor name="atomicCrossSectionIonizationPhoto_"  />
    !# <objectDestructor name="atomicRecombinationRateRadiative_"   />
    !# <objectDestructor name="atomicIonizationRateCollisional_"    />
    !# <objectDestructor name="atomicRecombinationRateDielectronic_"/>
    !# <objectDestructor name="atomicIonizationPotential_"          />
    return
  end function atomicConstructorParameters

  function atomicConstructorInternal(iterationAverageCount,massDistribution_,atomicCrossSectionIonizationPhoto_,atomicRecombinationRateRadiative_,atomicIonizationRateCollisional_,atomicRecombinationRateDielectronic_,atomicIonizationPotential_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily atomic} radiative transfer matter class.
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup, Atomic_Abundance, Atomic_Mass
    use :: Numerical_Constants_Astronomical, only : massSolar               , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Prefixes    , only : centi
    implicit none
    type   (radiativeTransferMatterAtomic           )                        :: self
    integer                                          , intent(in   )         :: iterationAverageCount
    class  (massDistributionClass                   ), intent(in   ), target :: massDistribution_
    class  (atomicCrossSectionIonizationPhotoClass  ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class  (atomicRecombinationRateRadiativeClass   ), intent(in   ), target :: atomicRecombinationRateRadiative_
    class  (atomicIonizationRateCollisionalClass    ), intent(in   ), target :: atomicIonizationRateCollisional_
    class  (atomicRecombinationRateDielectronicClass), intent(in   ), target :: atomicRecombinationRateDielectronic_
    class  (atomicIonizationPotentialClass          ), intent(in   ), target :: atomicIonizationPotential_
    !# <constructorAssign variables="iterationAverageCount, *massDistribution_, *atomicCrossSectionIonizationPhoto_, *atomicRecombinationRateRadiative_, *atomicIonizationRateCollisional_, *atomicRecombinationRateDielectronic_, *atomicIonizationPotential_"/>

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
  
  subroutine atomicPopulateDomain(self,properties,integrator)
    !% Populate a computational domain cell with atomic matter.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(radiativeTransferMatterAtomic           ), intent(inout), target :: self
    class(radiativeTransferMatterProperties       ), intent(inout)         :: properties
    class(computationalDomainVolumeIntegratorClass), intent(inout)         :: integrator

    select type (properties)
    type is (radiativeTransferMatterPropertiesAtomic)
       ! Compute the number density of this atomic species.
       properties%iterationCount         =+0
       properties%volume                 =+integrator%volume                       (                      )
       properties%densityNumber          =+integrator%integrate                    (atomicDensityIntegrand) &
            &                             /properties%volume                                                &
            &                             *self      %numberDensityMassDensityRatio
       ! Assume that the atomic species is fully neutral initially.
       properties%ionizationStateFraction=[1.0d0,0.0d0]
       ! Assign an initial temperature.
       properties%temperature            =1.0d4
       ! Initialize photoionization and photo heating rates to impossible values. These will be used as the "previous" value on
       ! the first iteration, guaranteeing that our first iteration is never judged to be converged.
       allocate(properties%photoIonizationRateHistory(self%iterationAverageCount,0:0))
       allocate(properties%photoHeatingRateHistory   (self%iterationAverageCount,0:0))
       properties%photoIonizationRate       =-huge(0.0d0)
       properties%photoHeatingRate          =-huge(0.0d0)
       properties%photoIonizationRateHistory=-huge(0.0d0)
       properties%photoHeatingRateHistory   =-huge(0.0d0)
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
  
  subroutine atomicStateSolve(self,properties)
    !% Solve for the state of the matter.
    use :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Numerical_Constants_Physical        , only : boltzmannsConstant
    implicit none
    class           (radiativeTransferMatterAtomic    ), intent(inout)  :: self
    class           (radiativeTransferMatterProperties), intent(inout)  :: properties
    double precision                                   , dimension(0:1) :: rateIonizationState                            , ionizationStateFractionPrevious
    double precision                                   , parameter      :: ionizationStateFractionToleranceAbsolute=1.0d-6, ionizationStateFractionToleranceRelative=1.0d-3, &
         &                                                                 ionizationStateFractionStepRelative     =5.0d-1
    integer                                            , parameter      :: countIterationMaximum                   =10000
    double precision                                                    :: ratePhotoIonization                            , rateRecombination                              , &
         &                                                                 timeStep                                       , rateImbalanceFactor                            , &
         &                                                                 densityNumberElectrons                         , rateNonPhotoIonization                         , &
         &                                                                 rateCoolingRecombinationRadiative              , rateCoolingRecombinationDielectronic           , &
         &                                                                 rateCoolingCollisionalExcitation               , rateCoolingBremsstrahlung                      , &
         &                                                                 ratePhotoHeating
    logical                                                             :: converged
    integer                                                             :: countIteration                                 , i
    
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
       ionizationStateFractionPrevious=properties%ionizationStateFraction
       countIteration                 =0
       converged                      =properties%densityNumber <= 0.0d0
       do while (.not.converged)
          ! Compute total electron density.
          densityNumberElectrons=+properties%densityNumber              &
               &                 *properties%ionizationStateFraction(1)
          ! Compute rates of change of the ionization states.       
          rateIonizationState   =+0.0d0
          ratePhotoIonization   =+properties%photoIonizationRate                      (  0                                                ) &
               &                 /properties%densityNumber
          rateNonPhotoIonization=+           densityNumberElectrons                                                                         &
               &                 *properties%ionizationStateFraction                  (  0                                                ) &
               &                 *  self    %atomicIonizationRateCollisional_    %rate(1,0,properties%temperature                         )   !! TO DO - hard-coded for hydrogen, fixed temperature
          rateRecombination     =+           densityNumberElectrons                                                                         &
               &                 *properties%ionizationStateFraction                  (  1                                                ) &
               &                 *(                                                                                                         &
               &                   +self    %atomicRecombinationRateRadiative_   %rate(1,1,properties%temperature,level=recombinationCaseB) & !! TO DO - hard-coded for hydrogen, fixed temperature, and case-B recombination
               &                   +self    %atomicRecombinationRateDielectronic_%rate(1,1,properties%temperature                         ) & !! TO DO - hard-coded for hydrogen, fixed temperature
               &                  )
          rateIonizationState(0)=-ratePhotoIonization-rateNonPhotoIonization+rateRecombination
          rateIonizationState(1)=+ratePhotoIonization+rateNonPhotoIonization-rateRecombination
          ! Compute heating/cooling rates (in J cm⁻³).
          rateNonPhotoIonization           =+properties%photoHeatingRate(0)
          rateCoolingRecombinationRadiative=+rateRecombination              &
               &                            *0.75d0                         &
               &                            *boltzmannsConstant             &
               &                            *properties%temperature           !! TO DO temperature fixed at 10^4 K here.
          ! Check for zero photoionization.
          if (ratePhotoionization <= 0.0d0) then
             properties%ionizationStateFraction(0  )=1.0d0
             properties%ionizationStateFraction(1:1)=0.0d0
          else
             ! Compute the equilibrium ionization state.
             !! Determine a time-step. We seek a solution where the net rate of change for each ionization state is zero. Therefore,
             !! we compute a "rate imbalance factor" which is the absolute value of the rate of change divided by the photoionization
             !! rate (which is fixed). The larger the factor, the bigger the step we should take to approach equilibrium.
             timeStep=huge(0.0d0)
             do i=0,1
                if (abs(rateIonizationState(i)) > 0.0d0 .and. properties%ionizationStateFraction(i) > 0.0d0) then
                   rateImbalanceFactor=min(1.0d0,abs(rateIonizationState(i))/ratePhotoIonization)
                   timeStep=min(timeStep,ionizationStateFractionStepRelative*rateImbalanceFactor*properties%ionizationStateFraction(i)/abs(rateIonizationState(i)))
                end if
             end do
             ! Adjust ionization state.
             properties%ionizationStateFraction=properties%ionizationStateFraction+timeStep*rateIonizationState
             if (any(properties%ionizationStateFraction < 0.0d0)) then
                where (properties%ionizationStateFraction < 0.0d0)
                   properties%ionizationStateFraction=0.0d0
                end where
                properties%ionizationStateFraction=properties%ionizationStateFraction/sum(properties%ionizationStateFraction)
             end if
          end if
          ! Check for convergence.
          converged=all(abs(properties%ionizationStateFraction-ionizationStateFractionPrevious) < max(ionizationStateFractionToleranceAbsolute,ionizationStateFractionToleranceRelative*properties%ionizationStateFraction))
          ionizationStateFractionPrevious=properties%ionizationStateFraction
          ! Check for exceeding maximum iterations.
          countIteration=countIteration+1
          if (countIteration > countIterationMaximum .and. .not.converged) call Galacticus_Error_Report('solution not found'//{introspection:location})
       end do
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine atomicStateSolve

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
          atomicOutputProperty=properties%densityNumber
       case (2_c_size_t)
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

    atomicCountOutputs=2_c_size_t
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
       atomicOutputName=var_str('densityNumber'   )
    case (2_c_size_t)
       atomicOutputName=var_str('fractionHydrogen')//Roman_Numerals(1+1)
    case default
       atomicOutputName=var_str(''                )
       call Galacticus_Error_Report('output is out of range'//{introspection:location})
    end select
    return
  end function atomicOutputName
