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

  !!{
  Implements a node operator class that performs star formation in disks.
  !!}

  use :: Star_Formation_Rates_Disks    , only : starFormationRateDisksClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass
  use :: Star_Formation_Histories      , only : starFormationHistoryClass

  !![
  <nodeOperator name="nodeOperatorStarFormationDisks">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationDisks
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class           (starFormationRateDisksClass     ), pointer :: starFormationRateDisks_      => null()
     class           (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class           (starFormationHistoryClass       ), pointer :: starFormationHistory_        => null()
     logical                                                     :: luminositiesStellarInactive
     double precision                                            :: fractionMassRetainedInitial           , fractionMassRetainedFinal
   contains
     final     ::                                        starFormationDisksDestructor
     procedure :: differentialEvolutionAnalytics      => starFormationDisksDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionPre            => starFormationDisksDifferentialEvolutionPre
     procedure :: differentialEvolution               => starFormationDisksDifferentialEvolution
     procedure :: differentialEvolutionStepFinalState => starFormationDisksDifferentialEvolutionStepFinalState
  end type nodeOperatorStarFormationDisks
  
  interface nodeOperatorStarFormationDisks
     !!{
     Constructors for the \refClass{nodeOperatorStarFormationDisks} node operator class.
     !!}
     module procedure starFormationDisksConstructorParameters
     module procedure starFormationDisksConstructorInternal
  end interface nodeOperatorStarFormationDisks
  
contains

  function starFormationDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorStarFormationDisks} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorStarFormationDisks  )                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (starFormationRateDisksClass     ), pointer       :: starFormationRateDisks_
    class  (stellarPopulationPropertiesClass), pointer       :: stellarPopulationProperties_
    class  (starFormationHistoryClass       ), pointer       :: starFormationHistory_
    logical                                                  :: luminositiesStellarInactive
    
    !![
    <inputParameter>
      <name>luminositiesStellarInactive</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, stellar luminosities will be treated as inactive properties.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="starFormationHistory"        name="starFormationHistory_"        source="parameters"/>
    !!]
    self=nodeOperatorStarFormationDisks(luminositiesStellarInactive,starFormationRateDisks_,stellarPopulationProperties_,starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"     />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="starFormationHistory_"       />
    !!]
    return
  end function starFormationDisksConstructorParameters

  function starFormationDisksConstructorInternal(luminositiesStellarInactive,starFormationRateDisks_,stellarPopulationProperties_,starFormationHistory_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorStarFormationDisks} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStarFormationDisks  )                        :: self
    class  (starFormationRateDisksClass     ), intent(in   ), target :: starFormationRateDisks_
    class  (stellarPopulationPropertiesClass), intent(in   ), target :: stellarPopulationProperties_
    class(starFormationHistoryClass         ), intent(in   ), target :: starFormationHistory_
    logical                                  , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *starFormationRateDisks_, *stellarPopulationProperties_, *starFormationHistory_"/>
    !!]

    ! Initialize values.
    self%fractionMassRetainedInitial=1.0d0
    self%fractionMassRetainedFinal  =1.0d0
    return
  end function starFormationDisksConstructorInternal

  subroutine starFormationDisksDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorStarFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"     />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%starFormationHistory_"       />
    !!]
    return
  end subroutine starFormationDisksDestructor
  
  subroutine starFormationDisksDifferentialEvolutionAnalytics(self,node)
    !!{
    Initialize the mass transfer fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorStarFormationDisks), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDisk             ), pointer       :: disk

    ! Mark the formed stellar mass as analytically-solvable (it is always zero) if we are not solving for luminosities as inactive
    ! properties.
    if (.not.self%luminositiesStellarInactive) then
       disk => node%disk()
       select type (disk)
       type is (nodeComponentDisk)
          ! Disk does not yet exist.
       class default
          call disk%massStellarFormedAnalytic()
       end select
    end if
    return
  end subroutine starFormationDisksDifferentialEvolutionAnalytics

  subroutine starFormationDisksDifferentialEvolutionPre(self,node)
    !!{
    Initialize the mass transfer fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorStarFormationDisks), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDisk             ), pointer       :: disk

    ! Initialize the mass transferred fraction to unity. The value is arbitrary as only ratios of this quantity are used, but
    ! must be non-zero.
    disk => node%disk()
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not yet exist.
    class default
       self%fractionMassRetainedInitial=1.0d0
       self%fractionMassRetainedFinal  =1.0d0
       if (disk%fractionMassRetainedIsSettable()) call disk%fractionMassRetainedSet(1.0d0)
    end select
    return
  end subroutine starFormationDisksDifferentialEvolutionPre

  subroutine starFormationDisksDifferentialEvolutionStepFinalState(self,node)
    !!{
    Record the final state of the disk at the end of the timestep prior to begin evaluation of integrals for inactive
    properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorStarFormationDisks), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDisk             ), pointer       :: disk

    ! The retained mass fraction at the start of this step is just the fraction at the end of the previous step. Then update the
    ! retained fraction at the end of the current step.
    disk                             =>           node%disk                     ()
    self%fractionMassRetainedInitial =            self%fractionMassRetainedFinal
    self%fractionMassRetainedFinal   =  max(0.0d0,disk%fractionMassRetained     ())
    return
  end subroutine starFormationDisksDifferentialEvolutionStepFinalState

  subroutine starFormationDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a disk.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : propertyInactive     , propertyTypeActive, propertyEvaluate, nodeComponentDisk, &
         &                                        nodeComponentSpheroid
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStarFormationDisks), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDisk             )               , pointer :: disk
    class           (nodeComponentSpheroid         )               , pointer :: spheroid
    double precision                                                         :: rateStarFormation       , massFuel                , &
         &                                                                      rateMassStellar         , rateEnergyInput         , &
         &                                                                      rateMassFuel            , fractionMassRetained    , &
         &                                                                      fractionMassRetainedRate
    logical                                                                  :: luminositiesCompute
    type            (abundances                    )                         :: abundancesFuel          , rateAbundancesFuels     , &
         &                                                                      rateAbundancesStellar
    type            (history                       )                         :: rateHistoryStarFormation, ratePropertiesStellar
    type            (stellarLuminosities           )                         :: rateLuminositiesStellar , rateLuminositiesTransfer
    
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (     disk%angularMomentum() < 0.0d0      &
         &  .or. disk%radius         () < 0.0d0      &
         &  .or. disk%massGas        () < 0.0d0      &
         &  .or.                                     &
         &   (                                       &
         &          propertyInactive(propertyType)   &
         &    .and.                                  & 
         &     .not.self%luminositiesStellarInactive &
         &   )                                       &
         & ) return
    if (propertyInactive(propertyType)) then
       ! For inactive property solution make use of the "massStellarFormed" property to determine the star formation rate.
       rateStarFormation=disk%massStellarFormedRateGet()
    else
       ! During active property solution, integrate the star formation rate so that we will have a solution for the total mass
       ! of stars formed as a function of time. This differs from the stellar mass due to recycling, and possibly transfer of
       ! stellar mass to other components.
       rateStarFormation=self%starFormationRateDisks_%rate(node)   
       if (self%luminositiesStellarInactive) call disk%massStellarFormedRate(rateStarFormation)
    end if
    ! Compute abundances of star forming gas.
    massFuel      =disk%massGas      ()
    abundancesFuel=disk%abundancesGas()
    call abundancesFuel%massToMassFraction(massFuel)
    ! Determine if luminosities must be computed.
    luminositiesCompute=propertyEvaluate(propertyType,self%luminositiesStellarInactive)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=disk%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                         &
         &                                       rateStarFormation      , &
         &                                       abundancesFuel         , &
         &                                       disk                   , &
         &                                       node                   , &
         &                                       ratePropertiesStellar  , &
         &                                       rateMassStellar        , &
         &                                       rateMassFuel           , &
         &                                       rateEnergyInput        , &
         &                                       rateAbundancesFuels    , &
         &                                       rateAbundancesStellar  , &
         &                                       rateLuminositiesStellar, &
         &                                       luminositiesCompute      &
         &                                      )
    ! Adjust rates.
    if (propertyEvaluate(propertyTypeActive,propertyIsInactive=.false.)) then
       rateHistoryStarFormation=disk%starFormationHistory()
       call        rateHistoryStarFormation%reset                       (                                                              )
       call self  %starFormationHistory_   %                        rate(node,rateHistoryStarFormation,abundancesFuel,rateStarFormation)
       call        disk                    %             massStellarRate(     rateMassStellar                                          )
       call        disk                    %                 massGasRate(     rateMassFuel                                             )
       call        disk                    %       abundancesStellarRate(     rateAbundancesStellar                                    )
       call        disk                    %           abundancesGasRate(     rateAbundancesFuels                                      )
       if (ratePropertiesStellar   %exists())                                                                                            &
            & call disk                    %stellarPropertiesHistoryRate(     ratePropertiesStellar                                    )
       if (rateHistoryStarFormation%exists())                                                                                            &
            & call disk                    %    starFormationHistoryRate(     rateHistoryStarFormation                                 )
    end if
    if (luminositiesCompute) then
       ! For inactive property calculations we must check if any mass (and, therefore, light) is being transferred to the
       ! spheroid component. If it is, our integrand must account for this mass transfer. The fractions of mass retained and
       ! transferred are determined from the "fractionMassRetained" property which is computed during differential evolution.
       if (propertyInactive(propertyType) .and. self%fractionMassRetainedFinal < self%fractionMassRetainedInitial) then
          spheroid => node%spheroid()
          ! Determine the fraction of mass (and light) formed at this time which will be retained in the disk at the final time in the step.
          if      (self%fractionMassRetainedFinal   == 0.0d0                      ) then
             ! The retained fraction reached zero by the end of the step, so no mass is retained.
             fractionMassRetained    =                                                                          0.0d0
          else if (self%fractionMassRetainedFinal   >  disk%fractionMassRetained()) then
             fractionMassRetained    =                                                                          1.0d0
          else
             ! Limit the retained fraction to unity (to avoid any rounding errors).
             fractionMassRetained    =    self%fractionMassRetainedFinal    /disk%fractionMassRetained       ()
          end if
          ! Determine the rate at which mass (and light) that was pre-existing at the start of this timestep is being transferred.
          if      (self%fractionMassRetainedInitial == 0.0d0                      ) then
             ! The initial retained fraction was zero, so there should be no light to transfer - set a transfer rate of zero.
             fractionMassRetainedRate=                                                                        0.0d0
          else
             ! Limit the transfer rate of pre-existing light to be negative - it is not possible to transfer light *to* the
             ! disk, so any positive value here can arise only via rounding errors.
             fractionMassRetainedRate=min(disk%fractionMassRetainedRateGet()/self%fractionMassRetainedInitial  ,0.0d0)
          end if
          ! Find the rate of transfer of pre-existing light.
          rateLuminositiesTransfer=+disk%luminositiesStellar() &
               &                   *fractionMassRetainedRate
          ! Evaluate the integrand for the disk, and the corresponding one for the spheroid to account for the transfer of light.
          call    disk %luminositiesStellarRate(rateLuminositiesStellar*       fractionMassRetained +rateLuminositiesTransfer                            )
          call spheroid%luminositiesStellarRate(rateLuminositiesStellar*(1.0d0-fractionMassRetained)-rateLuminositiesTransfer,interrupt,functionInterrupt)
       else
          ! In this case we do not need to account for transfer of light to the spheroid because either:
          !  a) there is none, or;
          !  b) we are solving for luminosities as active properties in which case transfer to the spheroid is handled directly in the ODE.
          call     disk%luminositiesStellarRate(rateLuminositiesStellar                                                                                  )
       end if
    end if
    return
  end subroutine starFormationDisksDifferentialEvolution
  
