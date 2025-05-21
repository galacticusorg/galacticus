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

  !+ Contributions to this file made by: Sachi Weerasooriya

  !!{
  An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  mimic the effects of reionization and accounting for cold mode accretion.
  !!}

  use :: Cooling_Functions, only : coolingFunction, coolingFunctionClass
  use :: Kind_Numbers     , only : kind_int8

  !![
  <accretionHalo name="accretionHaloColdMode">
   <description>
    Accretion onto halos using simple truncation to mimic the effects of reionization and accounting for cold mode
    accretion. This class extends the {\normalfont \ttfamily simple} class by dividing the accretion into hot and cold mode
    components. The cold mode fraction follows the approximation introduced by \cite{benson_cold_2010}, namely:
    \begin{equation}
    f_\mathrm{cold}=(1+r^{1/\delta})^{-1},
    \end{equation}
    where $\delta=${\normalfont \ttfamily [accretionColdModeShockStabilityTransitionWidth]}, $r =
    \epsilon_\mathrm{crit}/\epsilon$, and
    \begin{equation}
    \epsilon = r_\mathrm{s} \Lambda / \rho_\mathrm{s} v_\mathrm{s}^3,
    \end{equation}
    where $r_\mathrm{s}$ is the accretion shock radius, $\Lambda$ is the post-shock cooling function, $\rho_\mathrm{s}$ is the
    pre-shock density, $v_\mathrm{s}$ is the pre-shock velocity, and $\epsilon_\mathrm{crit}=${\normalfont \ttfamily
    [accretionColdModeShockStabilityThreshold]}. The pre-shock radius is set equal to the halo virial radius, the pre-shock
    velocity is set equal to the halo virial velocity, while the pre-shock density is given by
    \begin{equation}
    \rho_\mathrm{s} = {\gamma - 1 \over \gamma + 1} { 3 \over 4 \pi } { \Omega_\mathrm{b} \over \Omega_\mathrm{m} } {M \over
    r_\mathrm{s}^3} \left[ 1 + {(\alpha + 3) (10 + 9 \pi) \over 4} \right]^{-1},
    \end{equation}
    where $M$ is the total halo mass, $\gamma(=5/3)$ is the adiabatic index of the gas, and $\alpha$ is the exponent of the
    initial power-law density perturbation ($\alpha=0$ is assumed). The post-shock density and temperature are found assuming
    the strong-shock limit.
   </description>
  </accretionHalo>
  !!]
  type, extends(accretionHaloSimple) :: accretionHaloColdMode
     !!{
     A halo accretion class using simple truncation to mimic the effects of reionization and accounting for cold mode accretion.
     !!}
     private
     class           (coolingFunctionClass), pointer :: coolingFunction_        => null()
     double precision                                :: thresholdStabilityShock          , widthTransitionStabilityShock
     integer         (kind=kind_int8      )          :: lastUniqueID
     double precision                                :: coldFractionStored
     logical                                         :: coldFractionComputed
   contains
     !![
     <methods>
       <method description="Initialize after construction."                                                           method="initialize"      />
       <method description="Reset memoized calculations."                                                             method="calculationReset"/>
       <method description="Returns the total accretion rate from the \gls{igm} onto a halo (including dark matter)." method="chemicalMasses"  />
       <method description="Returns the total accretion rate from the \gls{igm} onto a halo (including dark matter)." method="coldModeFraction"/>
     </methods>
     !!]
     final     ::                              coldModeDestructor
     procedure :: initialize                => coldModeInitialize
     procedure :: autoHook                  => coldModeAutoHook
     procedure :: calculationReset          => coldModeCalculationReset
     procedure :: accretionRate             => coldModeAccretionRate
     procedure :: accretedMass              => coldModeAccretedMass
     procedure :: failedAccretionRate       => coldModeFailedAccretionRate
     procedure :: failedAccretedMass        => coldModeFailedAccretedMass
     procedure :: accretionRateMetals       => coldModeAccretionRateMetals
     procedure :: accretedMassMetals        => coldModeAccretedMassMetals
     procedure :: failedAccretionRateMetals => coldModeAccretionRateMetals
     procedure :: failedAccretedMassMetals  => coldModeAccretedMassMetals
     procedure :: accretionRateChemicals    => coldModeAccretionRateChemicals
     procedure :: accretedMassChemicals     => coldModeAccretedMassChemicals
     procedure :: chemicalMasses            => coldModeChemicalMasses
     procedure :: coldModeFraction          => coldModeColdModeFraction
  end type accretionHaloColdMode

  interface accretionHaloColdMode
     !!{
     Constructors for the \refClass{accretionHaloColdMode} halo accretion class.
     !!}
     module procedure coldModeConstructorParameters
     module procedure coldModeConstructorInternal
  end interface accretionHaloColdMode

contains

  function coldModeConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily coldMode} halo accretion class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (accretionHaloColdMode)                :: self
    type (inputParameters      ), intent(inout) :: parameters

    self%accretionHaloSimple=accretionHaloSimple(parameters)
    !![
    <inputParameter>
      <name>thresholdStabilityShock</name>
      <defaultSource>\citep{birnboim_virial_2003}</defaultSource>
      <defaultValue>0.0126d0</defaultValue>
      <description>The threshold value, $\epsilon_\mathrm{s,crit}$, for shock stability in the model of \cite{birnboim_virial_2003}.</description>
      <source>parameters</source>
      <variable>self%thresholdStabilityShock</variable>
    </inputParameter>
    <inputParameter>
      <name>widthTransitionStabilityShock</name>
      <defaultSource>\citep{benson_cold_2010}</defaultSource>
      <defaultValue>0.01d0</defaultValue>
      <description>The width of the transition from stability to instability for cold mode accretion \citep{benson_cold_2010}.</description>
      <source>parameters</source>
      <variable>self%widthTransitionStabilityShock</variable>
    </inputParameter>
    <objectBuilder class="coolingFunction" name="self%coolingFunction_" source="parameters"/>
    <inputParametersValidate source="parameters"/>
    !!]
    call self%initialize()
    return
  end function coldModeConstructorParameters

  function coldModeConstructorInternal(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,thresholdStabilityShock,widthTransitionStabilityShock,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_,coolingFunction_) result(self)
    !!{
    Internal constructor for the \refClass{accretionHaloColdMode} halo accretion class.
    !!}
    implicit none
    type            (accretionHaloColdMode        )                        :: self
    double precision                               , intent(in   )         :: timeReionization        , velocitySuppressionReionization, &
         &                                                                    thresholdStabilityShock , widthTransitionStabilityShock
    logical                                        , intent(in   )         :: accretionNegativeAllowed, accretionNewGrowthOnly
    class           (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (accretionHaloTotalClass      ), intent(in   ), target :: accretionHaloTotal_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    class           (chemicalStateClass           ), intent(in   ), target :: chemicalState_
    class           (coolingFunctionClass         ), intent(in   ), target :: coolingFunction_
    class           (intergalacticMediumStateClass), intent(in   ), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="thresholdStabilityShock, widthTransitionStabilityShock, *coolingFunction_"/>
    !!]

    self%accretionHaloSimple=accretionHaloSimple(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_)
    call self%initialize()
    return
  end function coldModeConstructorInternal

  subroutine coldModeInitialize(self)
    !!{
    Initialize the object after construction.
    !!}
    implicit none
    class(accretionHaloColdMode), intent(inout) :: self

    self%coldFractionComputed=.false.
    self%lastUniqueID        =-1_kind_int8
    return
  end subroutine coldModeInitialize

  subroutine coldModeAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(accretionHaloColdMode), intent(inout) :: self

    call calculationResetEvent%attach(self,coldModeCalculationReset,openMPThreadBindingAllLevels,label='accretionHaloColdMode')
    return
  end subroutine coldModeAutoHook

  subroutine coldModeDestructor(self)
    !!{
    Destructor for the \refClass{accretionHaloColdMode} halo accretion class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(accretionHaloColdMode), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingFunction_"/>
    !!]
    if (calculationResetEvent%isAttached(self,coldModeCalculationReset)) call calculationResetEvent%detach(self,coldModeCalculationReset)
    return
  end subroutine coldModeDestructor

  subroutine coldModeCalculationReset(self,node,uniqueID)
    !!{
    Reset the accretion rate calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer(kind_int8            ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%coldFractionComputed=.false.
    self%lastUniqueID        =uniqueID
    return
  end subroutine coldModeCalculationReset

  double precision function coldModeAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeAccretionRate=+self%accretionHaloSimple%accretionRate   (node,accretionModeTotal) &
         &                *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretionRate
  
  double precision function coldModeAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeAccretedMass=+self%accretionHaloSimple%accretedMass    (node,accretionModeTotal) &
         &               *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretedMass

  double precision function coldModeFailedAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeFailedAccretionRate=+self%accretionHaloSimple%failedAccretionRate(node,accretionModeTotal) &
         &                      *self                    %coldModeFraction   (node,accretionMode     )
    return
  end function coldModeFailedAccretionRate
  
  double precision function coldModeFailedAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeFailedAccretedMass=+self%accretionHaloSimple%failedAccretedMass(node,accretionModeTotal) &
         &                     *self                    %coldModeFraction  (node,accretionMode     )
    return
  end function coldModeFailedAccretedMass

  function coldModeAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    implicit none
    type (abundances                  )                :: coldModeAccretionRateMetals
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeAccretionRateMetals=+self%accretionHaloSimple%accretionRateMetals(node,accretionMode) &
         &                      *self                    %coldModeFraction   (node,accretionMode)
    return
  end function coldModeAccretionRateMetals

  function coldModeAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    implicit none
    type (abundances                  )                :: coldModeAccretedMassMetals
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeAccretedMassMetals=+self%accretionHaloSimple%accretedMassMetals(node,accretionMode) &
         &                     *self                    %coldModeFraction  (node,accretionMode)
    return
  end function coldModeAccretedMassMetals

  function coldModeFailedAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of failed mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    implicit none
    type (abundances                  )                :: coldModeFailedAccretionRateMetals
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeFailedAccretionRateMetals=+self%accretionHaloSimple%failedAccretionRateMetals(node,accretionMode) &
         &                            *self                    %coldModeFraction         (node,accretionMode)
    return
  end function coldModeFailedAccretionRateMetals

  function coldModeFailedAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances that failed to accrete (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    implicit none
    type (abundances                  )                :: coldModeFailedAccretedMassMetals
    class(accretionHaloColdMode       ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    coldModeFailedAccretedMassMetals=+self%accretionHaloSimple%failedAccretedMassMetals(node,accretionMode) &
         &                           *self                    %coldModeFraction        (node,accretionMode)
    return
  end function coldModeFailedAccretedMassMetals
  
  function coldModeAccretionRateChemicals(self,node,accretionMode)
    !!{
    Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium. Assumes a
    primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    temperature.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    type            (chemicalAbundances          )                :: coldModeAccretionRateChemicals
    class           (accretionHaloColdMode       ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    double precision                                              :: massAccretionRate

    ! Ensure that chemicals are reset to zero.
    call coldModeAccretionRateChemicals%reset()
    ! Return immediately if no chemicals are being tracked.
    if (self%countChemicals == 0) return
    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=self%accretionRate(node,accretionMode)
    ! Get the mass accretion rates.
    coldModeAccretionRateChemicals=self%chemicalMasses(node,massAccretionRate,accretionMode)
    return
  end function coldModeAccretionRateChemicals

  function coldModeAccretedMassChemicals(self,node,accretionMode)
    !!{
    Computes the mass of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    type            (chemicalAbundances          )                :: coldModeAccretedMassChemicals
    class           (accretionHaloColdMode       ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    double precision                                              :: massAccreted

    ! Ensure that chemicals are reset to zero.
    call coldModeAccretedMassChemicals%reset()
    ! Return if no chemicals are being tracked.
    if (self%countChemicals == 0) return
    ! Total mass of material accreted.
    massAccreted=self%accretedMass(node,accretionMode)
    ! Get the masses of chemicals accreted.
    coldModeAccretedMassChemicals=self%chemicalMasses(node,massAccreted,accretionMode)
    return
  end function coldModeAccretedMassChemicals

  function coldModeChemicalMasses(self,node,massAccreted,accretionMode)
    !!{
    Compute the masses of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure             , only : zeroAbundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                  , treeNode
    use :: Numerical_Constants_Astronomical , only : hydrogenByMassPrimordial
    use :: Numerical_Constants_Atomic       , only : atomicMassHydrogen
    implicit none
    type            (chemicalAbundances          )                :: coldModeChemicalMasses
    class           (accretionHaloColdMode       ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: massAccreted
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    type            (chemicalAbundances          ), save          :: chemicalDensities      , chemicalDensitiesHot , &
         &                                                           chemicalDensitiesCold
    !$omp threadprivate(chemicalDensities,chemicalDensitiesCold,chemicalDensitiesHot)
    double precision                                              :: massToDensityConversion, numberDensityHydrogen, &
         &                                                           temperatureHot         , temperatureCold      , &
         &                                                           fractionHot            , fractionCold

    ! Get the basic component.
    basic                => node%basic()
    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(self%darkMatterHaloScale_%radiusVirial(node))/3.0d0
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperatureHot            =  self%darkMatterHaloScale_     %temperatureVirial(node        )
    temperatureCold           =  self%intergalacticMediumState_%temperature      (basic%time())
    numberDensityHydrogen     =  hydrogenByMassPrimordial                  &
         &                       /atomicMassHydrogen                       &
         &                       *self %cosmologyParameters_%omegaBaryon() &
         &                       /self %cosmologyParameters_%omegaMatter() &
         &                       *basic                     %mass       () &
         &                       *massToDensityConversion
    ! Set the radiation field.
    call self%radiation%timeSet(basic%time())
    ! Get hot and cold mode fractions.
    fractionHot =self%coldModeFraction(node,accretionModeHot )
    fractionCold=self%coldModeFraction(node,accretionModeCold)
    ! Get the chemical densities.
    call self%chemicalState_%chemicalDensities(chemicalDensitiesHot ,numberDensityHydrogen,temperatureHot ,zeroAbundances,self%radiation)
    call self%chemicalState_%chemicalDensities(chemicalDensitiesCold,numberDensityHydrogen,temperatureCold,zeroAbundances,self%radiation)
    select case (accretionMode%ID)
    case (accretionModeTotal%ID)
       chemicalDensities=chemicalDensitiesHot*fractionHot+chemicalDensitiesCold*fractionCold
    case (accretionModeHot  %ID)
       chemicalDensities=chemicalDensitiesHot*fractionHot
    case (accretionModeCold %ID)
       chemicalDensities=                                 chemicalDensitiesCold*fractionCold
    end select
    ! Convert from densities to masses.
    call chemicalDensities%numberToMass(coldModeChemicalMasses)
    coldModeChemicalMasses=coldModeChemicalMasses*massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen
    return
  end function coldModeChemicalMasses

  double precision function coldModeColdModeFraction(self,node,accretionMode)
    !!{
    Computes the fraction of accretion occurring in the specified mode.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Chemical_Abundances_Structure   , only : chemicalAbundances
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentBasic      , treeNode
    use :: Numerical_Constants_Astronomical, only : hydrogenByMassPrimordial, massSolar         , meanAtomicMassPrimordial, megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassHydrogen      , atomicMassUnit
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : centi                   , kilo
    use :: Numerical_Constants_Units       , only : ergs
    use :: Shocks_1D                       , only : Shocks_1D_Density_Jump  , machNumberInfinite
    implicit none
    class           (accretionHaloColdMode       ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    double precision                              , parameter     :: adiabaticIndex             =5.0d0/3.0d0
    double precision                              , parameter     :: perturbationInitialExponent=0.0d0
    double precision                              , parameter     :: logStabilityRatioMaximum   =60.0d0
    class           (nodeComponentBasic          ), pointer       :: basic
    type            (chemicalAbundances          ), save          :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                              :: shockStability                         , stabilityRatio      , &
         &                                                           radiusShock                            , coolingFunctionValue, &
         &                                                           densityPreShock                        , densityPostShock    , &
         &                                                           numberDensityHydrogen                  , temperaturePostShock, &
         &                                                           velocityPreShock

    select case (accretionMode%ID)
    case (accretionModeTotal%ID)
       coldModeColdModeFraction=1.0d0
    case (accretionModeHot%ID,accretionModeCold%ID)
       ! Reset calculations if necessary.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
       ! Compute cold fraction if not already computed.
       if (.not.self%coldFractionComputed) then
          ! Get the basic component.
          basic => node%basic()
          ! Set the radiation field.
          call self%radiation%timeSet(basic%time())
          ! Compute factors required for stability analysis.
          radiusShock          =self%darkMatterHaloScale_%radiusVirial  (node)
          velocityPreShock     =self%darkMatterHaloScale_%velocityVirial(node)
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
               &                *basic                    %mass       ()    &
               &                *self%cosmologyParameters_%omegaBaryon()    &
               &                /self%cosmologyParameters_%omegaMatter()    &
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
          call self%chemicalState_%chemicalDensities(                            &
               &                                     chemicalDensities    ,      &
               &                                     numberDensityHydrogen,      &
               &                                     temperaturePostShock ,      &
               &                                     zeroAbundances       ,      &
               &                                     self%radiation              &
               &                                    )
          coolingFunctionValue=                                                             &
               &               self%coolingFunction_%coolingFunction(                       &
               &                                                     node                 , &
               &                                                     numberDensityHydrogen, &
               &                                                     temperaturePostShock , &
               &                                                     zeroAbundances       , &
               &                                                     chemicalDensities    , &
               &                                                     self%radiation         &
               &                                                    )
          ! Compute the shock stability parameter from Birnboim & Dekel (2003).
          shockStability=                           &
               &          megaParsec            **4 &
               &          /massSolar                &
               &          *ergs                     &
               &          /centi                **3 &
               &          /kilo                 **3 &
               &          /densityPreShock          &
               &          *radiusShock              &
               &          *coolingFunctionValue     &
               &          /velocityPreShock     **3
          ! Compute the cold fraction using the model from eqn. (2) of Benson & Bower (2011). The original form does not allow the
          ! cold fraction to go to zero in high mass halos, since "shockStability" can never be less than zero. This form is
          ! basically the equivalent functional form, but defined in terms of ln(ε) rather than ε.
          if (shockStability <= 0.0d0) then
             self%coldFractionStored=+0.0d0
          else
             stabilityRatio=self%thresholdStabilityShock/shockStability
             if (log(stabilityRatio) > self%widthTransitionStabilityShock*logStabilityRatioMaximum) then
                self%coldFractionStored=+0.0d0
             else
                self%coldFractionStored=+1.0d0                                                        &
                     &                  /(                                                            &
                     &                    +1.0d0                                                      &
                     &                    +stabilityRatio**(1.0d0/self%widthTransitionStabilityShock) &
                     &                   )
             end if
          end if
          ! Mark cold fraction as computed.
          self%coldFractionComputed=.true.
       end if
       ! Return the appropriate fraction.
       select case (accretionMode%ID)
       case (accretionModeHot %ID)
          coldModeColdModeFraction=1.0d0-self%coldFractionStored
       case (accretionModeCold%ID)
          coldModeColdModeFraction=     +self%coldFractionStored
       case default
          coldModeColdModeFraction=1.0d0
          call Error_Report('unknown accretion mode - this should not happen'//{introspection:location})
       end select
    case default
       coldModeColdModeFraction=1.0d0
       call Error_Report('unknown accretion mode - this should not happen'//{introspection:location})
    end select
    return
  end function coldModeColdModeFraction
