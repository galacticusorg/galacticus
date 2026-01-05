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
  Implements a node operator class that accumulates an estimate of the energy injected into the dark matter halo due to impulsive outflows.
  !!}

  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Stellar_Feedback_Outflows       , only : stellarFeedbackOutflowsClass
  use :: Star_Formation_Rates_Disks      , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids  , only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorImpulsiveOutflowEnergy">
    <description>
      A node operator class that accumulates an estimate of the energy injected into the dark matter halo due to impulsive
      outflows. The model assumed is that the energy injection is given by
      \begin{equation}
       \dot{\epsilon}(r) = \alpha \frac{\mathrm{G} \dot{M}_\mathrm{outflow}(r)}{r} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
      \end{equation}
      where $\alpha$ is a normalization factor, $t_\phi = M_\mathrm{gas}/\dot{M}_\mathrm{outflow}$ is the timescale for the
      outflow, and $t_\mathrm{dyn} = r_{1/2}/v_{1/2}$ is the dynamical time at the half-mass radius. The function $f(x)$ accounts
      for the fact that only impulsive changes in the potential should be accounted for, and is defined as
      \begin{equation}
       f(x) = ( 1 + \beta x )^\gamma,
      \end{equation}
      where $\beta=${\normalfont \ttfamily [impulsiveCorrectionScale]} and $\gamma=${\normalfont \ttfamily
      [impulsiveCorrectionExponent]}.

      In practice, this operator accumulates just
      \begin{equation}
      \dot{\epsilon}^\prime = \dot{M}_\mathrm{outflow} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
      \end{equation}
      allowing the radial dependence to be inserted later.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorImpulsiveOutflowEnergy
     !!{
     A node operator class that accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
     !!}
     private
     class           (darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_         => null()
     class           (stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflowsDisks_ => null(), stellarFeedbackOutflowsSpheroids_ => null()
     class           (starFormationRateDisksClass     ), pointer :: starFormationRateDisks_       => null()
     class           (starFormationRateSpheroidsClass ), pointer :: starFormationRateSpheroids_   => null()
     class           (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_  => null()
     class           (mergerMassMovementsClass        ), pointer :: mergerMassMovements_          => null()
     integer                                                     :: energyImpulsiveOutflowDiskID           , energyImpulsiveOutflowSpheroidID
     double precision                                            :: impulsiveCorrectionScale               , impulsiveCorrectionExponent
   contains
     final     ::                                impulsiveOutflowEnergyDestructor
     procedure :: differentialEvolutionScales => impulsiveOutflowEnergyDifferentialEvolutionScales
     procedure :: differentialEvolution       => impulsiveOutflowEnergyDifferentialEvolution
     procedure :: nodePromote                 => impulsiveOutflowEnergyNodePromote
     procedure :: galaxiesMerge               => impulsiveOutflowEnergyGalaxiesMerge
  end type nodeOperatorImpulsiveOutflowEnergy
  
  interface nodeOperatorImpulsiveOutflowEnergy
     !!{
     Constructors for the \refClass{nodeOperatorImpulsiveOutflowEnergy} node operator class.
     !!}
     module procedure impulsiveOutflowEnergyConstructorParameters
     module procedure impulsiveOutflowEnergyConstructorInternal
  end interface nodeOperatorImpulsiveOutflowEnergy
  
contains

  function impulsiveOutflowEnergyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorImpulsiveOutflowEnergy} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorImpulsiveOutflowEnergy)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class           (stellarFeedbackOutflowsClass      ), pointer       :: stellarFeedbackOutflowsDisks_, stellarFeedbackOutflowsSpheroids_
    class           (starFormationRateDisksClass       ), pointer       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass   ), pointer       :: starFormationRateSpheroids_
    class           (stellarPopulationPropertiesClass  ), pointer       :: stellarPopulationProperties_
    class           (mergerMassMovementsClass          ), pointer       :: mergerMassMovements_
    type            (inputParameters                   )                :: parametersDisks              , parametersSpheroids
    double precision                                                    :: impulsiveCorrectionScale     , impulsiveCorrectionExponent

    parametersDisks    =parameters%subParameters('disks'    ,requireValue=.false.)
    parametersSpheroids=parameters%subParameters('spheroids',requireValue=.false.)
    !![
    <inputParameter>
      <name>impulsiveCorrectionScale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter $\beta$ appearing in the impulsive correction function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>impulsiveCorrectionExponent</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The parameter $\gamma$ appearing in the impulsive correction function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"        name="darkMatterProfileDMO_"             source="parameters"         />
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflowsDisks_"     source="parametersDisks"    />
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflowsSpheroids_" source="parametersSpheroids"/>
    <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"           source="parameters"         />
    <objectBuilder class="starFormationRateSpheroids"  name="starFormationRateSpheroids_"       source="parameters"         />
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_"      source="parameters"         />
    <objectBuilder class="mergerMassMovements"         name="mergerMassMovements_"              source="parameters"         />
    !!]
     self=nodeOperatorImpulsiveOutflowEnergy(impulsiveCorrectionScale,impulsiveCorrectionExponent,darkMatterProfileDMO_,stellarFeedbackOutflowsDisks_,stellarFeedbackOutflowsSpheroids_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationProperties_,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"            />
    <objectDestructor name="stellarFeedbackOutflowsDisks_"    />
    <objectDestructor name="stellarFeedbackOutflowsSpheroids_"/>
    <objectDestructor name="starFormationRateDisks_"          />
    <objectDestructor name="starFormationRateSpheroids_"      />
    <objectDestructor name="stellarPopulationProperties_"     />
    <objectDestructor name="mergerMassMovements_"             />
    !!]
    return
  end function impulsiveOutflowEnergyConstructorParameters

  function impulsiveOutflowEnergyConstructorInternal(impulsiveCorrectionScale,impulsiveCorrectionExponent,darkMatterProfileDMO_,stellarFeedbackOutflowsDisks_,stellarFeedbackOutflowsSpheroids_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationProperties_,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorImpulsiveOutflowEnergy} node operator class.
    !!}
    implicit none
    type            (nodeOperatorImpulsiveOutflowEnergy)                        :: self
    class           (darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    class           (stellarFeedbackOutflowsClass      ), intent(in   ), target :: stellarFeedbackOutflowsDisks_, stellarFeedbackOutflowsSpheroids_
    class           (starFormationRateDisksClass       ), intent(in   ), target :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass   ), intent(in   ), target :: starFormationRateSpheroids_
    class           (stellarPopulationPropertiesClass  ), intent(in   ), target :: stellarPopulationProperties_
    class           (mergerMassMovementsClass          ), intent(in   ), target :: mergerMassMovements_
    double precision                                    , intent(in   )         :: impulsiveCorrectionScale     , impulsiveCorrectionExponent
    !![
    <constructorAssign variables="impulsiveCorrectionScale, impulsiveCorrectionExponent, *darkMatterProfileDMO_, *stellarFeedbackOutflowsDisks_, *stellarFeedbackOutflowsSpheroids_, *starFormationRateDisks_, *starFormationRateSpheroids_, *stellarPopulationProperties_, *mergerMassMovements_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowDisk"     id="self%energyImpulsiveOutflowDiskID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowSpheroid" id="self%energyImpulsiveOutflowSpheroidID" isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function impulsiveOutflowEnergyConstructorInternal

  subroutine impulsiveOutflowEnergyDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorImpulsiveOutflowEnergy} node operator class.
    !!}
    implicit none
    type(nodeOperatorImpulsiveOutflowEnergy), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"            />
    <objectDestructor name="self%stellarFeedbackOutflowsDisks_"    />
    <objectDestructor name="self%stellarFeedbackOutflowsSpheroids_"/>
    <objectDestructor name="self%starFormationRateDisks_"          />
    <objectDestructor name="self%starFormationRateSpheroids_"      />
    <objectDestructor name="self%stellarPopulationProperties_"     />
    <objectDestructor name="self%mergerMassMovements_"             />
    !!]
    return
  end subroutine impulsiveOutflowEnergyDestructor

  subroutine impulsiveOutflowEnergyDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the impulsive mass outflowed.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorImpulsiveOutflowEnergy), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , parameter     :: scaleRelative=1.0d-6
    class           (nodeComponentBasic                ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile    ), pointer       :: darkMatterProfile
    
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertyScale(self%energyImpulsiveOutflowDiskID    ,basic%mass()*scaleRelative)
    call darkMatterProfile%floatRank0MetaPropertyScale(self%energyImpulsiveOutflowSpheroidID,basic%mass()*scaleRelative)
    return
  end subroutine impulsiveOutflowEnergyDifferentialEvolutionScales
  
  subroutine impulsiveOutflowEnergyDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulates an estimate of the impulsive mass outflowed.
    !!}
    use :: Abundances_Structure            , only : abundances
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, nodeComponentDisk, nodeComponentSpheroid
    use :: Histories                       , only : history
    use :: Stellar_Luminosities_Structure  , only : stellarLuminosities
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (nodeOperatorImpulsiveOutflowEnergy), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout), target  :: node
    logical                                             , intent(inout)          :: interrupt
    procedure       (interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                             , intent(in   )          :: propertyType
    class           (nodeComponentDarkMatterProfile    )               , pointer :: darkMatterProfile
    class           (nodeComponentDisk                 )               , pointer :: disk
    class           (nodeComponentSpheroid             )               , pointer :: spheroid
    double precision                                                             :: rateStarFormationDisk       , rateStarFormationSpheroid       , &
         &                                                                          massGasDisk                 , massGasSpheroid                 , &
         &                                                                          rateEnergyInputDisk         , rateEnergyInputSpheroid         , &
         &                                                                          rateMassStellarDisk         , rateMassStellarSpheroid         , &
         &                                                                          rateMassFuelDisk            , rateMassFuelSpheroid            , &
         &                                                                          rateMassOutflowEjectiveDisk , rateMassOutflowEjectiveSpheroid , &
         &                                                                          rateMassOutflowExpulsiveDisk, rateMassOutflowExpulsiveSpheroid, &
         &                                                                          rateMassOutflowTotalDisk    , rateMassOutflowTotalSpheroid    , &
         &                                                                          timescaleOutflowDisk        , timescaleOutflowSpheroid        , &
         &                                                                          timescaleDynamicalDisk      , timescaleDynamicalSpheroid
    type            (abundances                      )                           :: abundancesGasDisk           , abundancesGasSpheroid           , &
         &                                                                          rateAbundancesFuelsDisk     , rateAbundancesFuelsSpheroid     , &
         &                                                                          rateAbundancesStellarDisk   , rateAbundancesStellarSpheroid
    type            (history                         )                           :: ratePropertiesStellarDisk   , ratePropertiesStellarSpheroid
    type            (stellarLuminosities             )                           :: rateLuminositiesStellarDisk , rateLuminositiesStellarSpheroid
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType
    
    ! Evaluate star formation rates and outflow rates.
    darkMatterProfile             => node                                %darkMatterProfile       (    )
    disk                          => node                                %disk                    (    )
    spheroid                      => node                                %spheroid                (    )
    rateStarFormationDisk         =  self    %starFormationRateDisks_    %rate                    (node)   
    rateStarFormationSpheroid     =  self    %starFormationRateSpheroids_%rate                    (node)   
    massGasDisk                   =  disk                                %massGas                 (    )
    massGasSpheroid               =  spheroid                            %massGas                 (    )
    abundancesGasDisk             =  disk                                %abundancesGas           (    )
    abundancesGasSpheroid         =  spheroid                            %abundancesGas           (    )
    ratePropertiesStellarDisk     =  disk                                %stellarPropertiesHistory(    )
    ratePropertiesStellarSpheroid =  spheroid                            %stellarPropertiesHistory(    )
    call abundancesGasDisk                                      %massToMassFraction(                                                               &
         &                                                                                                        massGasDisk                      &
         &                                                                         )
    call abundancesGasSpheroid                                  %massToMassFraction(                                                               &
         &                                                                                                        massGasSpheroid                  &
         &                                                                         )
    call self                 %stellarPopulationProperties_     %rates             (                                                               &
         &                                                                                                        rateStarFormationDisk          , &
         &                                                                                                        abundancesGasDisk              , &
         &                                                                                                        disk                           , &
         &                                                                                                        node                           , &
         &                                                                                                        ratePropertiesStellarDisk      , &
         &                                                                                                        rateMassStellarDisk            , &
         &                                                                                                        rateMassFuelDisk               , &
         &                                                                                                        rateEnergyInputDisk            , &
         &                                                                                                        rateAbundancesFuelsDisk        , &
         &                                                                                                        rateAbundancesStellarDisk      , &
         &                                                                                                        rateLuminositiesStellarDisk    , &
         &                                                                           computeRateLuminosityStellar=.false.                          &
         &                                                                          )
    call self                  %stellarPopulationProperties_     %rates              (                                                              &
         &                                                                                                        rateStarFormationSpheroid      , &
         &                                                                                                        abundancesGasSpheroid          , &
         &                                                                                                        spheroid                       , &
         &                                                                                                        node                           , &
         &                                                                                                        ratePropertiesStellarSpheroid  , &
         &                                                                                                        rateMassStellarSpheroid        , &
         &                                                                                                        rateMassFuelSpheroid           , &
         &                                                                                                        rateEnergyInputSpheroid        , &
         &                                                                                                        rateAbundancesFuelsSpheroid    , &
         &                                                                                                        rateAbundancesStellarSpheroid  , &
         &                                                                                                        rateLuminositiesStellarSpheroid, &
         &                                                                           computeRateLuminosityStellar=.false.                          &
         &                                                                          )
    call self                  %stellarFeedbackOutflowsDisks_    %outflowRate        (                                                              &
         &                                                                                                        disk                            , &
         &                                                                                                        rateStarFormationDisk           , &
         &                                                                                                        rateEnergyInputDisk             , &
         &                                                                                                        rateMassOutflowEjectiveDisk     , &
         &                                                                                                        rateMassOutflowExpulsiveDisk      &
         &                                                                           )
    call self                  %stellarFeedbackOutflowsSpheroids_%outflowRate        (                                                              &
         &                                                                                                        spheroid                        , &
         &                                                                                                        rateStarFormationSpheroid       , &
         &                                                                                                        rateEnergyInputSpheroid         , &
         &                                                                                                        rateMassOutflowEjectiveSpheroid , &
         &                                                                                                        rateMassOutflowExpulsiveSpheroid  &
         &                                                                           )
    ! Compute total outflow rates and relevant timescales.
    rateMassOutflowTotalDisk    =+rateMassOutflowEjectiveDisk    +rateMassOutflowExpulsiveDisk
    rateMassOutflowTotalSpheroid=+rateMassOutflowEjectiveSpheroid+rateMassOutflowExpulsiveSpheroid
    if (rateMassOutflowTotalDisk     > 0.0d0) then
       timescaleOutflowDisk      =+massGasDisk    /rateMassOutflowTotalDisk
    else
       timescaleOutflowDisk      =+0.0d0
    end if
    if (rateMassOutflowTotalSpheroid > 0.0d0) then
       timescaleOutflowSpheroid  =+massGasSpheroid/rateMassOutflowTotalSpheroid
    else
       timescaleOutflowSpheroid  =+0.0d0
    end if
    if (disk    %velocity() > 0.0d0) then
       timescaleDynamicalDisk    =+disk    %radius()/disk    %velocity()*MpcPerKmPerSToGyr
    else
       timescaleDynamicalDisk    =+0.0d0
    end if 
    if (spheroid%velocity() > 0.0d0) then
       timescaleDynamicalSpheroid=+spheroid%radius()/spheroid%velocity()*MpcPerKmPerSToGyr
    else
       timescaleDynamicalSpheroid=+0.0d0
    end if 
    ! Set the rate of change for the impulsive energy input.
    call darkMatterProfile%floatRank0MetaPropertyRate(self%energyImpulsiveOutflowDiskID    ,rateMassOutflowTotalDisk    *impulsiveCorrection(timescaleOutflowDisk    ,timescaleDynamicalDisk    ))
    call darkMatterProfile%floatRank0MetaPropertyRate(self%energyImpulsiveOutflowSpheroidID,rateMassOutflowTotalSpheroid*impulsiveCorrection(timescaleOutflowSpheroid,timescaleDynamicalSpheroid))
    return
    
  contains
    
    double precision function impulsiveCorrection(timescaleOutflow,timescaleDynamical)
      !!{
      Compute the correction factor for impulsive outflows.
      !!}
      implicit none
      double precision, intent(in   ) :: timescaleOutflow, timescaleDynamical

      if (timescaleDynamical > 0.0d0) then
         impulsiveCorrection=+1.0d0                               &
              &              /(                                   &
              &                +1.0d0                             &
              &                +self%impulsiveCorrectionScale     &
              &                *     timescaleOutflow             &
              &                /     timescaleDynamical           &
              &               )**self%impulsiveCorrectionExponent
      else
         impulsiveCorrection=+0.0d0
      end if
      return
    end function impulsiveCorrection
    
  end subroutine impulsiveOutflowEnergyDifferentialEvolution

  subroutine impulsiveOutflowEnergyNodePromote(self,node)
    !!{
    Sum any impulsive energy inputs in the parent and child halos prior to promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorImpulsiveOutflowEnergy), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile    ), pointer       :: darkMatterProfile, darkMatterProfileParent

    darkMatterProfile       => node       %darkMatterProfile(autoCreate=.true.)
    darkMatterProfileParent => node%parent%darkMatterProfile(autoCreate=.true.)
    call darkMatterProfile%floatRank0MetaPropertySet(                                                                                           &
         &                                                                                              self%energyImpulsiveOutflowDiskID     , &
         &                                           +darkMatterProfile      %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
         &                                           +darkMatterProfileParent%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
         &                                          )
    call darkMatterProfile%floatRank0MetaPropertySet(                                                                                           &
         &                                                                                              self%energyImpulsiveOutflowSpheroidID , &
         &                                           +darkMatterProfile      %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
         &                                           +darkMatterProfileParent%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
         &                                          )
    return
  end subroutine impulsiveOutflowEnergyNodePromote

  subroutine impulsiveOutflowEnergyGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk         , destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatorImpulsiveOutflowEnergy), intent(inout) :: self
    type   (treeNode                          ), intent(inout) :: node
    type   (treeNode                          ), pointer       :: nodeHost
    class  (nodeComponentDarkMatterProfile    ), pointer       :: darkMatterProfile      , darkMatterProfileHost
    type   (enumerationDestinationMergerType  )                :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                        destinationGasHost     , destinationStarsHost
    logical                                                    :: mergerIsMajor

    ! Find the node to merge with.
    nodeHost              => node    %mergesWith       ()
    darkMatterProfile     => node    %darkMatterProfile()
    darkMatterProfileHost => nodeHost%darkMatterProfile()
    ! Get mass movement descriptors.
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Move impulsive energy within the host if necessary.
    select case (destinationGasHost     %ID)
    case (destinationMergerDisk    %ID)
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowDiskID     , &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
            &                                              )
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowSpheroidID , &
            &                                               +0.0d0                                                                                   &
            &                                              )
     case (destinationMergerSpheroid%ID)
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowSpheroidID , &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
            &                                              )
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowDiskID     , &
            &                                               +0.0d0                                                                                   &
            &                                              )
    case (destinationMergerUnmoved%ID)
       ! Do nothing.
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Move impulsive energy from secondary to primary.
    select case (destinationGasSatellite%ID)
    case (destinationMergerDisk    %ID)
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowDiskID     , &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
            &                                               +darkMatterProfile    %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
            &                                               +darkMatterProfile    %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
            &                                              )
    case (destinationMergerSpheroid%ID)
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowSpheroidID , &
            &                                               +darkMatterProfileHost%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
            &                                               +darkMatterProfile    %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    )  &
            &                                               +darkMatterProfile    %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID)  &
            &                                              )
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Zero rates in the secondary,
    call    darkMatterProfile    %floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowDiskID     , &
         &                                                  +0.0d0                                                                                   &
         &                                                 )
    call    darkMatterProfile    %floatRank0MetaPropertySet(                                                 self%energyImpulsiveOutflowSpheroidID , &
         &                                                  +0.0d0                                                                                   &
         &                                                 )
    return
  end subroutine impulsiveOutflowEnergyGalaxiesMerge
  
