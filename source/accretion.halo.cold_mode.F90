!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  !% mimic the effects of reionization and accounting for cold mode accretion.

  !# <accretionHalo name="accretionHaloColdMode">
  !#  <description>Accretion onto halos using simple truncation to mimic the effects of reionization and accounting for cold mode accretion.</description>
  !# </accretionHalo>

  type, extends(accretionHaloSimple) :: accretionHaloColdMode
     !% A halo accretion class using simple truncation to mimic the effects of reionization and accounting for cold mode accretion.
     private
     double precision :: shockStabilityThreshold, shockStabilityTransitionWidth
   contains
     procedure :: accretionRate          => coldModeAccretionRate
     procedure :: accretedMass           => coldModeAccretedMass
     procedure :: failedAccretionRate    => coldModeFailedAccretionRate
     procedure :: failedAccretedMass     => coldModeFailedAccretedMass
     procedure :: accretionRateMetals    => coldModeAccretionRateMetals
     procedure :: accretedMassMetals     => coldModeAccretedMassMetals
     procedure :: accretionRateChemicals => coldModeAccretionRateChemicals
     procedure :: accretedMassChemicals  => coldModeAccretedMassChemicals
     procedure :: chemicalMasses         => coldModeChemicalMasses
     procedure :: coldModeFraction       => coldModeColdModeFraction
  end type accretionHaloColdMode

  interface accretionHaloColdMode
     !% Constructors for the {\tt coldMode} halo accretion class.
     module procedure coldModeConstructor
     module procedure coldModeDefaultConstructor
  end interface accretionHaloColdMode

  interface assignment(=)
     module procedure coldModeFromSimple
  end interface assignment(=)

  ! Default parameters.
  double precision :: coldModeShockStabilityThreshold, coldModeShockStabilityTransitionWidth

  ! Initialization state.
  logical          :: coldModeDefaultInitialized=.false.

contains
  
  subroutine coldModeFromSimple(coldMode,simple)
    !% Assign a {\tt simple} halo accretion object to a {\tt coldMode} halo accretion object.
    implicit none
    type(accretionHaloColdMode), intent(inout) :: coldMode
    type(accretionHaloSimple  ), intent(in   ) :: simple
    
    coldMode%reionizationSuppressionTime    =simple%reionizationSuppressionTime
    coldMode%reionizationSuppressionVelocity=simple%reionizationSuppressionVelocity
    coldMode%negativeAccretionAllowed       =simple%negativeAccretionAllowed
    coldMode%accreteNewGrowthOnly           =simple%accreteNewGrowthOnly
    coldMode%radiation                      =simple%radiation
    return
  end subroutine coldModeFromSimple
  
  function coldModeDefaultConstructor()
    !% Default constructor for the {\tt coldMode} halo accretion class.
    use Intergalactic_Medium_State
    use Cosmology_Functions
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type(accretionHaloColdMode), target  :: coldModeDefaultConstructor

    ! Get default parameters.
    if (.not.coldModeDefaultInitialized) then
       !$omp critical(accretionHaloColdModeDefaultInitialize)
       if (.not.coldModeDefaultInitialized) then
          ! Read parameters controlling the cold mode.
          !@ <inputParameter>
          !@   <name>accretionColdModeShockStabilityThreshold</name>
          !@   <defaultValue>0.0126 \citep{birnboim_virial_2003}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The threshold value, $\epsilon_{\rm s,crit}$, for shock stability in the model of \cite{birnboim_virial_2003}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("accretionColdModeShockStabilityThreshold",coldModeShockStabilityThreshold,defaultValue=0.0126d0)
          !@ <inputParameter>
          !@   <name>accretionColdModeShockStabilityTransitionWidth</name>
          !@   <defaultValue>0.01 \citep{benson_cold_2010}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The width of the transition from stability to instability for cold mode accretion \citep{benson_cold_2010}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("accretionColdModeShockStabilityTransitionWidth",coldModeShockStabilityTransitionWidth,defaultValue=0.01d0)
          ! Record that class is now initialized.
          coldModeDefaultInitialized=.true.
       end if
       !$omp end critical(accretionHaloColdModeDefaultInitialize)
    end if
    coldModeDefaultConstructor=accretionHaloSimple()
    coldModeDefaultConstructor%shockStabilityThreshold      =coldModeShockStabilityThreshold
    coldModeDefaultConstructor%shockStabilityTransitionWidth=coldModeShockStabilityTransitionWidth
    return
  end function coldModeDefaultConstructor
       
  function coldModeConstructor(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly,shockStabilityThreshold,shockStabilityTransitionWidth)
    !% Default constructor for the {\tt coldMode} halo accretion class.
    implicit none
    type            (accretionHaloColdMode), target        :: coldModeConstructor
    double precision                       , intent(in   ) :: reionizationSuppressionTime, reionizationSuppressionVelocity, &
         &                                                    shockStabilityThreshold    , shockStabilityTransitionWidth
    logical                                , intent(in   ) :: negativeAccretionAllowed   , accreteNewGrowthOnly

    coldModeConstructor=accretionHaloSimple(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    coldModeConstructor%shockStabilityThreshold      =shockStabilityThreshold
    coldModeConstructor%shockStabilityTransitionWidth=shockStabilityTransitionWidth
    return
  end function coldModeConstructor

  double precision function coldModeAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt node}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    class  (accretionHaloColdMode), intent(inout)          :: self
    type   (treeNode             ), intent(inout), pointer :: node
    integer                       , intent(in   )          :: accretionMode

    coldModeAccretionRate=                                                     &
         &  self%accretionHaloSimple%accretionRate   (node,accretionModeTotal) &
         & *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretionRate
  
  double precision function coldModeAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\tt node}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    class  (accretionHaloColdMode), intent(inout)          :: self
    type   (treeNode             ), intent(inout), pointer :: node
    integer                       , intent(in   )          :: accretionMode

    coldModeAccretedMass=                                                      &
         &  self%accretionHaloSimple%accretedMass    (node,accretionModeTotal) &
         & *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretedMass

  double precision function coldModeFailedAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt node}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    class  (accretionHaloColdMode), intent(inout)          :: self
    type   (treeNode             ), intent(inout), pointer :: node
    integer                       , intent(in   )          :: accretionMode

    coldModeFailedAccretionRate=                                                  &
         &  self%accretionHaloSimple%failedAccretionRate(node,accretionModeTotal) &
         & *self                    %coldModeFraction   (node,accretionMode     )
    return
  end function coldModeFailedAccretionRate

  double precision function coldModeFailedAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\tt node}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    class  (accretionHaloColdMode), intent(inout)          :: self
    type   (treeNode             ), intent(inout), pointer :: node
    integer                       , intent(in   )          :: accretionMode

    coldModeFailedAccretedMass=                                                  &
         &  self%accretionHaloSimple%failedAccretedMass(node,accretionModeTotal) &
         & *self                    %coldModeFraction  (node,accretionMode     )
    return
  end function coldModeFailedAccretedMass

  function coldModeAccretionRateMetals(self,node,accretionMode)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type  (abundances           )                         :: coldModeAccretionRateMetals
    class (accretionHaloColdMode), intent(inout)          :: self
    type  (treeNode             ), intent(inout), pointer :: node
    integer                      , intent(in   )          :: accretionMode
    
    coldModeAccretionRateMetals=                                             &
         &  self%accretionHaloSimple%accretionRateMetals(node,accretionMode) &
         & *self                    %coldModeFraction   (node,accretionMode)
    return
  end function coldModeAccretionRateMetals

  function coldModeAccretedMassMetals(self,node,accretionMode)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type   (abundances           )                         :: coldModeAccretedMassMetals
    class  (accretionHaloColdMode), intent(inout)          :: self
    type   (treeNode             ), intent(inout), pointer :: node
    integer                       , intent(in   )          :: accretionMode

    coldModeAccretedMassMetals=                                             &
         &  self%accretionHaloSimple%accretedMassMetals(node,accretionMode) &
         & *self                    %coldModeFraction  (node,accretionMode)
    return
  end function coldModeAccretedMassMetals

  function coldModeAccretionRateChemicals(self,node,accretionMode)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\tt node} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances   )                         :: coldModeAccretionRateChemicals
    class           (accretionHaloColdMode), intent(inout)          :: self
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(in   )          :: accretionMode
    double precision                                                :: massAccretionRate

    ! Return immediately if no chemicals are being tracked.
    if (simpleChemicalsCount == 0) return
    ! Ensure that chemicals are reset to zero.
    call coldModeAccretionRateChemicals%reset()
    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=self%accretionRate(node,accretionMode)
    ! Get the mass accretion rates.
    coldModeAccretionRateChemicals=self%chemicalMasses(node,massAccretionRate,accretionMode)
    return
  end function coldModeAccretionRateChemicals

  function coldModeAccretedMassChemicals(self,node,accretionMode)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances   )                         :: coldModeAccretedMassChemicals
    class           (accretionHaloColdMode), intent(inout)          :: self
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(in   )          :: accretionMode
    double precision                                                :: massAccreted

    ! Return immediately if no chemicals are being tracked.
    if (simpleChemicalsCount == 0) return
    ! Ensure that chemicals are reset to zero.
    call coldModeAccretedMassChemicals%reset()
    ! Total mass of material accreted.
    massAccreted=self%accretedMass(node,accretionMode)
    ! Get the masses of chemicals accreted.
    coldModeAccretedMassChemicals=self%chemicalMasses(node,massAccreted,accretionMode)
    return
  end function coldModeAccretedMassChemicals

  function coldModeChemicalMasses(self,node,massAccreted,accretionMode)
    !% Compute the masses of chemicals accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Intergalactic_Medium_State
    implicit none
    type            (chemicalAbundances           )                         :: coldModeChemicalMasses
    class           (accretionHaloColdMode        ), intent(inout)          :: self
    type            (treeNode                     ), intent(inout), pointer :: node
    double precision                               , intent(in   )          :: massAccreted
    integer                                        , intent(in   )          :: accretionMode
    class           (nodeComponentBasic           )               , pointer :: thisBasicComponent
    class           (cosmologyParametersClass     )               , pointer :: thisCosmologyParameters
    class           (intergalacticMediumStateClass)               , pointer :: intergalacticMediumState_
    class           (darkMatterHaloScaleClass     )               , pointer :: darkMatterHaloScale_
    type            (chemicalAbundances           ), save                   :: chemicalDensities      , chemicalDensitiesHot , &
         &                                                                     chemicalDensitiesCold
    !$omp threadprivate(chemicalDensities,chemicalDensitiesCold,chemicalDensitiesHot)
    double precision                                                        :: massToDensityConversion, numberDensityHydrogen, &
         &                                                                     temperature            , temperatureHot       , &
         &                                                                     temperatureCold        , fractionCold         , &
         &                                                                     fractionHot

    ! Get required objects.
    thisCosmologyParameters => cosmologyParameters()
    darkMatterHaloScale_    => darkMatterHaloScale()
    ! Get the basic component.
    thisBasicComponent   => node%basic()
    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(darkMatterHaloScale_%virialRadius(node))/3.0d0
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperatureHot            =  darkMatterHaloScale_%virialTemperature(node                 )
    intergalacticMediumState_ => intergalacticMediumState              (                         )
    temperature               =  intergalacticMediumState_%temperature (thisBasicComponent%time())
    numberDensityHydrogen     =  hydrogenByMassPrimordial*(thisCosmologyParameters%omegaBaryon()/thisCosmologyParameters%omegaMatter())*thisBasicComponent%mass()*massToDensityConversion&
         &/atomicMassHydrogen
    ! Set the radiation field.
    call self%radiation%set(node)
    ! Get hot and cold mode fractions.
    fractionHot =self%coldModeFraction(node,accretionModeHot )
    fractionCold=self%coldModeFraction(node,accretionModeCold)
    ! Get the chemical densities.
    call Chemical_Densities(chemicalDensitiesHot ,temperatureHot ,numberDensityHydrogen,zeroAbundances,self%radiation)
    call Chemical_Densities(chemicalDensitiesCold,temperatureCold,numberDensityHydrogen,zeroAbundances,self%radiation)
    select case (accretionMode)
    case (accretionModeTotal)
       chemicalDensities=chemicalDensitiesHot*fractionHot+chemicalDensitiesCold*fractionCold
    case (accretionModeHot  )
       chemicalDensities=chemicalDensitiesHot*fractionHot
    case (accretionModeCold)
       chemicalDensities=                                 chemicalDensitiesCold*fractionCold
    end select
    ! Convert from densities to masses.
    call chemicalDensities%numberToMass(coldModeChemicalMasses)
    coldModeChemicalMasses=coldModeChemicalMasses*massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen
    return
  end function coldModeChemicalMasses

  double precision function coldModeColdModeFraction(self,node,accretionMode)
    !% Computes the fraction of accretion occuring in the specified mode.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Shocks_1D
    use Numerical_Constants_Atomic
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Cooling_Functions
    implicit none
    class           (accretionHaloColdMode   ), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: node
    integer                                   , intent(in   )          :: accretionMode
    double precision                          , parameter              :: adiabaticIndex             =5.0d0/3.0d0  
    double precision                          , parameter              :: perturbationInitialExponent=0.0d0
    double precision                          , parameter              :: logStabilityRatioMaximum   =60.0d0
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    class           (nodeComponentBasic      )               , pointer :: thisBasic
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    type            (chemicalAbundances      ), save                   :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                                   :: shockStability       , coldFraction        , &
         &                                                                radiusShock          , coolingFunction     , &
         &                                                                densityPreShock      , densityPostShock    , &
         &                                                                numberDensityHydrogen, temperaturePostShock, &
         &                                                                velocityPreShock     , stabilityRatio

    select case (accretionMode)
    case (accretionModeTotal)
       coldModeColdModeFraction=1.0d0
    case (accretionModeHot,accretionModeCold)
       ! Get required objects.
       thisCosmologyParameters => cosmologyParameters()
       darkMatterHaloScale_    => darkMatterHaloScale()
       ! Set the radiation field.
       call self%radiation%set(node)
       ! Get the basic component.
       thisBasic => node%basic()
       ! Compute factors required for stability analysis.
       radiusShock          =darkMatterHaloScale_%virialRadius  (node)
       velocityPreShock     =darkMatterHaloScale_%virialVelocity(node)
       temperaturePostShock =                                            &
            &                 (3.0d0/16.0d0)                             &
            &                *atomicMassUnit                             &
            &                *meanAtomicMassPrimordial                   &
            &                *(kilo*velocityPreShock)**2                 &
            &                /boltzmannsConstant
       densityPreShock      =                                            &
            &                 (adiabaticIndex-1.0d0)                     &
            &                /(adiabaticIndex+1.0d0)                     &
            &                *(3.0d0/4.0d0/Pi)                           &
            &                *thisBasic%mass()                           &
            &                *thisCosmologyParameters%omegaBaryon()      &
            &                /thisCosmologyParameters%omegaMatter()      &
            &                /radiusShock**3                             &
            &                /(                                          &
            &                   1.0d0                                    &
            &                  +(perturbationInitialExponent+3.0d0)      &
            &                  *(10.0d0+9.0d0*Pi)                        &
            &                  /4.0d0                                    &
            &                 )
       densityPostShock     =                                            &
            &                 densityPreShock                            &
            &                *Shocks_1D_Density_Jump(                    &
            &                                        adiabaticIndex    , &
            &                                        machNumberInfinite  &
            &                                       )
       numberDensityHydrogen=                                            &
            &                 massSolar                                  &
            &                /megaParsec              **3                &
            &                *centi                   **3                &
            &                *densityPreShock                            &
            &                *hydrogenByMassPrimordial                   &
            &                /atomicMassUnit                             &
            &                /atomicMassHydrogen
       call Chemical_Densities(                                          &
            &                  chemicalDensities    ,                    &
            &                  temperaturePostShock ,                    &
            &                  numberDensityHydrogen,                    &
            &                  zeroAbundances       ,                    &
            &                  self%radiation                            &
            &                 )
       coolingFunction     =                                             &
            &               Cooling_Function(                            &
            &                                temperaturePostShock ,      &
            &                                numberDensityHydrogen,      &
            &                                zeroAbundances       ,      &
            &                                chemicalDensities    ,      &
            &                                self%radiation              &
            &                               )
       ! Compute the shock stability parameter from Birnboim & Dekel (2003).
       shockStability=                      &
            &          megaParsec       **4 &
            &          /massSolar           &
            &          *ergs                &
            &          /centi           **3 &
            &          /kilo            **3 &
            &          /densityPreShock     &
            &          *radiusShock         &
            &          *coolingFunction     &
            &          /velocityPreShock**3
       ! Compute the cold fraction using the model from eqn. (2) of Benson & Bower (2011). The original form doesn't allow the
       ! cold fraction to go to zero in high mass halos, since "shockStability" can never be less than zero. This form is
       ! basically the equivalent functional form, but defined in terms of ln(epsilon) rather than epsilon.
       stabilityRatio=self%shockStabilityThreshold/shockStability
       if (log(stabilityRatio) > self%shockStabilityTransitionWidth*logStabilityRatioMaximum) then
          coldFraction=0.0d0
       else
          coldFraction=                                                               &
               &        1.0d0                                                         &
               &        /(                                                            &
               &           1.0d0                                                      &
               &          +stabilityRatio**(1.0d0/self%shockStabilityTransitionWidth) &
               &         )
       end if
       ! Return the appropriate fraction.
       select case (accretionMode)
       case (accretionModeHot )
          coldModeColdModeFraction=1.0d0-coldFraction
       case (accretionModeCold)
          coldModeColdModeFraction=     +coldFraction
       end select
    end select
    return
  end function coldModeColdModeFraction
