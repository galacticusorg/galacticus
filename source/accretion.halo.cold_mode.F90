!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  !% mimic the effects of reionization and accounting for cold mode accretion.

  use Cooling_Functions

  !# <accretionHalo name="accretionHaloColdMode">
  !#  <description>Accretion onto halos using simple truncation to mimic the effects of reionization and accounting for cold mode accretion.</description>
  !# </accretionHalo>
  type, extends(accretionHaloSimple) :: accretionHaloColdMode
     !% A halo accretion class using simple truncation to mimic the effects of reionization and accounting for cold mode accretion.
     private
     class           (coolingFunctionClass), pointer :: coolingFunction_
     double precision                                :: thresholdStabilityShock, widthTransitionStabilityShock
     integer         (kind=kind_int8      )          :: lastUniqueID
     double precision                                :: coldFractionStored
     logical                                         :: coldFractionComputed
   contains
     !@ <objectMethods>
     !@   <object>accretionHaloColdMode</object>
     !@   <objectMethod>
     !@     <method>chemicalMasses</method>
     !@     <type>\textcolor{red}{\textless type(chemicalAbundances)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout, \doublezero massAccreted\argin, \intzero accretionMode\argin</arguments>
     !@     <description>Returns the total accretion rate from the \gls{igm} onto a halo (including dark matter).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>coldModeFraction</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout, \intzero accretionMode\argin</arguments>
     !@     <description>Returns the total accretion rate from the \gls{igm} onto a halo (including dark matter).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                           coldModeDestructor
     procedure :: calculationReset       => coldModeCalculationReset
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
     !% Constructors for the {\normalfont \ttfamily coldMode} halo accretion class.
     module procedure coldModeConstructorParameters
     module procedure coldModeConstructorInternal
  end interface accretionHaloColdMode

contains
  
  function coldModeConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily coldMode} halo accretion class.
    use Input_Parameters
    implicit none
    type (accretionHaloColdMode)                :: self
    type (inputParameters      ), intent(inout) :: parameters

    self%accretionHaloSimple=accretionHaloSimple(parameters)
    !# <inputParameter>
    !#   <name>thresholdStabilityShock</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>\citep{birnboim_virial_2003}</defaultSource>
    !#   <defaultValue>0.0126d0</defaultValue>
    !#   <description>The threshold value, $\epsilon_\mathrm{s,crit}$, for shock stability in the model of \cite{birnboim_virial_2003}.</description>
    !#   <source>globalParameters</source>
    !#   <type>real</type>
    !#   <variable>self%thresholdStabilityShock</variable>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>widthTransitionStabilityShock</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>\citep{benson_cold_2010}</defaultSource>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The width of the transition from stability to instability for cold mode accretion \citep{benson_cold_2010}.</description>
    !#   <source>globalParameters</source>
    !#   <type>real</type>
    !#   <variable>self%widthTransitionStabilityShock</variable>
    !# </inputParameter>
    !# <objectBuilder class="coolingFunction" name="self%coolingFunction_" source="parameters"/>
    !# <inputParametersValidate source="parameters"/>
    self%coldFractionComputed=.false.
    return
  end function coldModeConstructorParameters

  function coldModeConstructorInternal(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,thresholdStabilityShock,widthTransitionStabilityShock,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_,coolingFunction_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily coldMode} halo accretion class.
    use Galacticus_Error
    use Atomic_Data
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
    !# <constructorAssign variables="thresholdStabilityShock, widthTransitionStabilityShock, *coolingFunction_"/>

    self%accretionHaloSimple=accretionHaloSimple(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_)
    self%coldFractionComputed=.false.
    return
  end function coldModeConstructorInternal

  subroutine coldModeDestructor(self)
    !% Destructor for the {\normalfont \ttfamily coldMode} halo accretion class.
    implicit none
    type(accretionHaloColdMode), intent(inout) :: self

    !# <objectDestructor name="self%coolingFunction_"/>
    return
  end subroutine coldModeDestructor

  subroutine coldModeCalculationReset(self,node)
    !% Reset the accretion rate calculation.
    implicit none
    class(accretionHaloColdMode), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node

    self%coldFractionComputed=.false.
    self%lastUniqueID        =node%uniqueID()
    return
  end subroutine coldModeCalculationReset

  double precision function coldModeAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer                       , intent(in   ) :: accretionMode

    coldModeAccretionRate=+self%accretionHaloSimple%accretionRate   (node,accretionModeTotal) &
         &                *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretionRate
  
  double precision function coldModeAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer                       , intent(in   ) :: accretionMode

    coldModeAccretedMass=+self%accretionHaloSimple%accretedMass    (node,accretionModeTotal) &
         &               *self                    %coldModeFraction(node,accretionMode     )
    return
  end function coldModeAccretedMass

  double precision function coldModeFailedAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer                       , intent(in   ) :: accretionMode

    coldModeFailedAccretionRate=+self%accretionHaloSimple%failedAccretionRate(node,accretionModeTotal) &
         &                      *self                    %coldModeFraction   (node,accretionMode     )
    return
  end function coldModeFailedAccretionRate

  double precision function coldModeFailedAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer                       , intent(in   ) :: accretionMode

    coldModeFailedAccretedMass=+self%accretionHaloSimple%failedAccretedMass(node,accretionModeTotal) &
         &                     *self                    %coldModeFraction  (node,accretionMode     )
    return
  end function coldModeFailedAccretedMass

  function coldModeAccretionRateMetals(self,node,accretionMode)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type  (abundances           )                :: coldModeAccretionRateMetals
    class (accretionHaloColdMode), intent(inout) :: self
    type  (treeNode             ), intent(inout) :: node
    integer                      , intent(in   ) :: accretionMode
    
    coldModeAccretionRateMetals=+self%accretionHaloSimple%accretionRateMetals(node,accretionMode) &
         &                      *self                    %coldModeFraction   (node,accretionMode)
    return
  end function coldModeAccretionRateMetals

  function coldModeAccretedMassMetals(self,node,accretionMode)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type   (abundances           )                :: coldModeAccretedMassMetals
    class  (accretionHaloColdMode), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer                       , intent(in   ) :: accretionMode

    coldModeAccretedMassMetals=+self%accretionHaloSimple%accretedMassMetals(node,accretionMode) &
         &                     *self                    %coldModeFraction  (node,accretionMode)
    return
  end function coldModeAccretedMassMetals

  function coldModeAccretionRateChemicals(self,node,accretionMode)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances   )                :: coldModeAccretionRateChemicals
    class           (accretionHaloColdMode), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    integer                                , intent(in   ) :: accretionMode
    double precision                                       :: massAccretionRate

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
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances   )                :: coldModeAccretedMassChemicals
    class           (accretionHaloColdMode), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    integer                                , intent(in   ) :: accretionMode
    double precision                                       :: massAccreted

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
    !% Compute the masses of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Numerical_Constants_Astronomical
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    implicit none
    type            (chemicalAbundances   )                :: coldModeChemicalMasses
    class           (accretionHaloColdMode), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: massAccreted
    integer                                , intent(in   ) :: accretionMode
    class           (nodeComponentBasic   ), pointer       :: basic
    type            (chemicalAbundances   ), save          :: chemicalDensities      , chemicalDensitiesHot , &
         &                                                    chemicalDensitiesCold
    !$omp threadprivate(chemicalDensities,chemicalDensitiesCold,chemicalDensitiesHot)
    double precision                                       :: massToDensityConversion, numberDensityHydrogen, &
         &                                                    temperature            , temperatureHot       , &
         &                                                    temperatureCold        , fractionCold         , &
         &                                                    fractionHot

    ! Get the basic component.
    basic                => node%basic()
    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(self%darkMatterHaloScale_%virialRadius(node))/3.0d0
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperatureHot            =  self%darkMatterHaloScale_     %virialTemperature(node        )
    temperature               =  self%intergalacticMediumState_%temperature      (basic%time())
    numberDensityHydrogen     =  hydrogenByMassPrimordial*(self%cosmologyParameters_%omegaBaryon()/self%cosmologyParameters_%omegaMatter())*basic%mass()*massToDensityConversion&
         &/atomicMassHydrogen
    ! Set the radiation field.
    call self%radiation%set(node)
    ! Get hot and cold mode fractions.
    fractionHot =self%coldModeFraction(node,accretionModeHot )
    fractionCold=self%coldModeFraction(node,accretionModeCold)
    ! Get the chemical densities.
    call self%chemicalState_%chemicalDensities(chemicalDensitiesHot ,numberDensityHydrogen,temperatureHot ,zeroAbundances,self%radiation)
    call self%chemicalState_%chemicalDensities(chemicalDensitiesCold,numberDensityHydrogen,temperatureCold,zeroAbundances,self%radiation)
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
    use Galacticus_Error
    use Shocks_1D
    use Numerical_Constants_Atomic
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    implicit none
    class           (accretionHaloColdMode   ), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    integer                                   , intent(in   ) :: accretionMode
    double precision                          , parameter     :: adiabaticIndex             =5.0d0/3.0d0  
    double precision                          , parameter     :: perturbationInitialExponent=0.0d0
    double precision                          , parameter     :: logStabilityRatioMaximum   =60.0d0
    class           (nodeComponentBasic      ), pointer       :: basic
    type            (chemicalAbundances      ), save          :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                          :: shockStability                         , stabilityRatio      , &
         &                                                       radiusShock                            , coolingFunctionValue, &
         &                                                       densityPreShock                        , densityPostShock    , &
         &                                                       numberDensityHydrogen                  , temperaturePostShock, &
         &                                                       velocityPreShock

    select case (accretionMode)
    case (accretionModeTotal)
       coldModeColdModeFraction=1.0d0
    case (accretionModeHot,accretionModeCold)
       ! Reset calculations if necessary.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       ! Compute cold fraction if not already computed.
       if (.not.self%coldFractionComputed) then
          ! Set the radiation field.
          call self%radiation%set(node)
          ! Get the basic component.
          basic => node%basic()
          ! Compute factors required for stability analysis.
          radiusShock          =self%darkMatterHaloScale_%virialRadius  (node)
          velocityPreShock     =self%darkMatterHaloScale_%virialVelocity(node)
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
               &          /velocityPreShock**3
          ! Compute the cold fraction using the model from eqn. (2) of Benson & Bower (2011). The original form doesn't allow the
          ! cold fraction to go to zero in high mass halos, since "shockStability" can never be less than zero. This form is
          ! basically the equivalent functional form, but defined in terms of ln(epsilon) rather than epsilon.
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
          ! Mark cold fraction as computed.
          self%coldFractionComputed=.true.
       end if
       ! Return the appropriate fraction.
       select case (accretionMode)
       case (accretionModeHot )
          coldModeColdModeFraction=1.0d0-self%coldFractionStored
       case (accretionModeCold)
          coldModeColdModeFraction=     +self%coldFractionStored
       case default
          coldModeColdModeFraction=1.0d0
          call Galacticus_Error_Report('unknown accretion mode - this should not happen'//{introspection:location})
       end select
    case default
       coldModeColdModeFraction=1.0d0
       call Galacticus_Error_Report('unknown accretion mode - this should not happen'//{introspection:location})  
    end select
    return
  end function coldModeColdModeFraction
