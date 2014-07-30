!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of baryonic accretion onto halos using a simple truncation to mimic
!% reionization.

module Accretion_Halos_Cold_Mode
  !% Implements calculations of baryonic accretion onto halos using a simple truncation to mimic reionization and accounting for
  !% cold mode accretion.
  use Radiation_Structure
  use Abundances_Structure
  use Accretion_Halos_Options
  implicit none
  private
  public :: Accretion_Halos_Cold_Mode_Initialize

  ! Internal record of the number of chemicals being tracked.
  integer                              :: chemicalsCount
  ! Radiation structure.
  type            (radiationStructure) :: radiation
  !$omp threadprivate(radiation)
  ! Parameters controlling onset of cold mode.
  double precision                     :: accretionColdModeShockStabilityThreshold
  double precision                     :: accretionColdModeShockStabilityTransitionWidth

contains

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_Cold_Mode_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_Cold_Mode_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
       &,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get &
       &,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accreted_Abundances_Get&
       &,Halo_Baryonic_Accretion_Rate_Chemicals_Get,Halo_Baryonic_Accreted_Chemicals_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmology_Functions
    use Atomic_Data
    use Galacticus_Error
    use Chemical_Abundances_Structure
    use Intergalactic_Medium_State
    implicit none
    type            (varying_string                                       ), intent(in   )          :: accretionHalosMethod
    procedure       (Halo_Baryonic_Accretion_Rate_Cold_Mode_Get           ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Get
    procedure       (Halo_Baryonic_Accreted_Mass_Cold_Mode_Get            ), intent(inout), pointer :: Halo_Baryonic_Accreted_Mass_Get
    procedure       (Halo_Baryonic_Failed_Accretion_Rate_Cold_Mode_Get    ), intent(inout), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get
    procedure       (Halo_Baryonic_Failed_Accreted_Mass_Cold_Mode_Get     ), intent(inout), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get
    procedure       (Halo_Baryonic_Accretion_Rate_Abundances_Cold_Mode_Get), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Abundances_Get
    procedure       (Halo_Baryonic_Accreted_Abundances_Cold_Mode_Get      ), intent(inout), pointer :: Halo_Baryonic_Accreted_Abundances_Get
    procedure       (Halo_Baryonic_Accretion_Rate_Chemicals_Cold_Mode_Get ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Chemicals_Get
    procedure       (Halo_Baryonic_Accreted_Chemicals_Cold_Mode_Get       ), intent(inout), pointer :: Halo_Baryonic_Accreted_Chemicals_Get

    if (accretionHalosMethod == 'coldMode') then
       ! Set pointers to our implementations of accretion functions.
       Halo_Baryonic_Accretion_Rate_Get            => Halo_Baryonic_Accretion_Rate_Cold_Mode_Get
       Halo_Baryonic_Accreted_Mass_Get             => Halo_Baryonic_Accreted_Mass_Cold_Mode_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get     => Halo_Baryonic_Failed_Accretion_Rate_Cold_Mode_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get      => Halo_Baryonic_Failed_Accreted_Mass_Cold_Mode_Get
       Halo_Baryonic_Accretion_Rate_Abundances_Get => Halo_Baryonic_Accretion_Rate_Abundances_Cold_Mode_Get
       Halo_Baryonic_Accreted_Abundances_Get       => Halo_Baryonic_Accreted_Abundances_Cold_Mode_Get
       Halo_Baryonic_Accretion_Rate_Chemicals_Get  => Halo_Baryonic_Accretion_Rate_Chemicals_Cold_Mode_Get
       Halo_Baryonic_Accreted_Chemicals_Get        => Halo_Baryonic_Accreted_Chemicals_Cold_Mode_Get
       ! Get a count of the number of chemicals being tracked.
       chemicalsCount=Chemicals_Property_Count()
       ! Define the radiation structure.
       call radiation%define([radiationTypeCMB])
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
       call Get_Input_Parameter("accretionColdModeShockStabilityThreshold",accretionColdModeShockStabilityThreshold,defaultValue=0.0126d0)
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
       call Get_Input_Parameter("accretionColdModeShockStabilityTransitionWidth",accretionColdModeShockStabilityTransitionWidth,defaultValue=0.01d0)
    end if
    return
  end subroutine Accretion_Halos_Cold_Mode_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_Cold_Mode_Get(thisNode,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type   (treeNode), intent(inout), pointer :: thisNode
    integer          , intent(in   )          :: accretionMode

    Halo_Baryonic_Accretion_Rate_Cold_Mode_Get=                                     &
         &  Halo_Baryonic_Accretion_Rate_Simple_Get   (thisNode,accretionModeTotal) &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionMode     )
    return
  end function Halo_Baryonic_Accretion_Rate_Cold_Mode_Get
  
  double precision function Halo_Baryonic_Accreted_Mass_Cold_Mode_Get(thisNode,accretionMode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type   (treeNode          ), intent(inout), pointer :: thisNode
    integer                    , intent(in   )          :: accretionMode

    Halo_Baryonic_Accreted_Mass_Cold_Mode_Get=                                      &
         &  Halo_Baryonic_Accreted_Mass_Simple_Get    (thisNode,accretionModeTotal) &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionMode     )
    return
  end function Halo_Baryonic_Accreted_Mass_Cold_Mode_Get

  double precision function Halo_Baryonic_Failed_Accretion_Rate_Cold_Mode_Get(thisNode,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    integer                               , intent(in   )          :: accretionMode

    Halo_Baryonic_Failed_Accretion_Rate_Cold_Mode_Get=                                  &
         &  Halo_Baryonic_Failed_Accretion_Rate_Simple_Get(thisNode,accretionModeTotal) &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction    (thisNode,accretionMode     )
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_Cold_Mode_Get

  double precision function Halo_Baryonic_Failed_Accreted_Mass_Cold_Mode_Get(thisNode,accretionMode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type   (treeNode          ), intent(inout), pointer :: thisNode
    integer                    , intent(in   )          :: accretionMode

    Halo_Baryonic_Failed_Accreted_Mass_Cold_Mode_Get=                                  &
         &  Halo_Baryonic_Failed_Accreted_Mass_Simple_Get(thisNode,accretionModeTotal) &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction   (thisNode,accretionMode     )
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_Cold_Mode_Get

  subroutine Halo_Baryonic_Accretion_Rate_Abundances_Cold_Mode_Get(thisNode,accretionRateAbundances,accretionMode)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type  (treeNode  ), intent(inout), pointer :: thisNode
    type  (abundances), intent(inout)          :: accretionRateAbundances
    integer           , intent(in   )          :: accretionMode
    
    call Halo_Baryonic_Accretion_Rate_Abundances_Simple_Get(thisNode,accretionRateAbundances,accretionModeTotal)
    accretionRateAbundances=                                                   &
         &  accretionRateAbundances                                            &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionMode)
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances_Cold_Mode_Get

  subroutine Halo_Baryonic_Accreted_Abundances_Cold_Mode_Get(thisNode,accretedAbundances,accretionMode)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Accretion_Halos_Simple
    implicit none
    type   (treeNode  ), intent(inout), pointer :: thisNode
    type   (abundances), intent(inout)          :: accretedAbundances
    integer            , intent(in   )          :: accretionMode

    call Halo_Baryonic_Accreted_Abundances_Simple_Get(thisNode,accretedAbundances,accretionModeTotal)
    accretedAbundances=                                                        &
         &  accretedAbundances                                                 &
         & *Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionMode)
    return
  end subroutine Halo_Baryonic_Accreted_Abundances_Cold_Mode_Get

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Cold_Mode_Get(thisNode,accretionRateChemicals,accretionMode)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (chemicalAbundances), intent(inout)          :: accretionRateChemicals
    integer                             , intent(in   )          :: accretionMode
    double precision                                             :: massAccretionRate

    ! Return immediately if no chemicals are being tracked.
    if (chemicalsCount == 0) return
    ! Ensure that chemicals are reset to zero.
    call accretionRateChemicals%reset()
    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=Halo_Baryonic_Accretion_Rate_Cold_Mode_Get(thisNode,accretionMode)
    ! Get the mass accretion rates.
    call Get_Chemical_Masses(thisNode,massAccretionRate,accretionRateChemicals,accretionMode)
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Cold_Mode_Get

  subroutine Halo_Baryonic_Accreted_Chemicals_Cold_Mode_Get(thisNode,accretedChemicals,accretionMode)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (chemicalAbundances), intent(inout)          :: accretedChemicals
    integer                             , intent(in   )          :: accretionMode
    double precision                                             :: massAccreted

    ! Return immediately if no chemicals are being tracked.
    if (chemicalsCount == 0) return
    ! Ensure that chemicals are reset to zero.
    call accretedChemicals%reset()
    ! Total mass of material accreted.
    massAccreted=Halo_Baryonic_Accreted_Mass_Cold_Mode_Get(thisNode,accretionMode)
    ! Get the masses of chemicals accreted.
    call Get_Chemical_Masses(thisNode,massAccreted,accretedChemicals,accretionMode)
    return
  end subroutine Halo_Baryonic_Accreted_Chemicals_Cold_Mode_Get

  subroutine Get_Chemical_Masses(thisNode,massAccreted,chemicalMasses,accretionMode)
    !% Compute the masses of chemicals accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Intergalactic_Medium_State
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    double precision                    , intent(in   )          :: massAccreted
    type            (chemicalAbundances), intent(  out)          :: chemicalMasses
    integer                             , intent(in   )          :: accretionMode
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    class           (cosmologyParametersClass)         , pointer :: thisCosmologyParameters
    class           (intergalacticMediumStateClass)    , pointer :: intergalacticMediumState_
    type            (chemicalAbundances), save                   :: chemicalDensities      , chemicalDensitiesHot , &
         &                                                          chemicalDensitiesCold
    !$omp threadprivate(chemicalDensities,chemicalDensitiesCold,chemicalDensitiesHot)
    double precision                                             :: massToDensityConversion, numberDensityHydrogen, &
         &                                                          temperature            , temperatureHot       , &
         &                                                          temperatureCold        , fractionCold         , &
         &                                                          fractionHot

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Get the basic component.
    thisBasicComponent   => thisNode%basic()
    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(thisNode))/3.0d0
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperatureHot            =  Dark_Matter_Halo_Virial_Temperature  (thisNode                 )
    intergalacticMediumState_ => intergalacticMediumState             (                         )
    temperature               =  intergalacticMediumState_%temperature(thisBasicComponent%time())
    numberDensityHydrogen     =  hydrogenByMassPrimordial*(thisCosmologyParameters%omegaBaryon()/thisCosmologyParameters%omegaMatter())*thisBasicComponent%mass()*massToDensityConversion&
         &/atomicMassHydrogen
    ! Set the radiation field.
    call radiation%set(thisNode)
    ! Get hot and cold mode fractions.
    fractionHot =Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionModeHot )
    fractionCold=Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionModeCold)
    ! Get the chemical densities.
    call Chemical_Densities(chemicalDensitiesHot ,temperatureHot ,numberDensityHydrogen,zeroAbundances,radiation)
    call Chemical_Densities(chemicalDensitiesCold,temperatureCold,numberDensityHydrogen,zeroAbundances,radiation)
    select case (accretionMode)
    case (accretionModeTotal)
       chemicalDensities=chemicalDensitiesHot*fractionHot+chemicalDensitiesCold*fractionCold
    case (accretionModeHot  )
       chemicalDensities=chemicalDensitiesHot*fractionHot
    case (accretionModeCold)
       chemicalDensities=                                 chemicalDensitiesCold*fractionCold
    end select
    ! Convert from densities to masses.
    call chemicalDensities%numberToMass(chemicalMasses)
    chemicalMasses=chemicalMasses*massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen
    return
  end subroutine Get_Chemical_Masses

  double precision function Halo_Baryonic_Accretion_Cold_Mode_Fraction(thisNode,accretionMode)
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
    type            (treeNode          ), intent(inout), pointer :: thisNode
    integer                             , intent(in   )          :: accretionMode
    double precision                    , parameter              :: adiabaticIndex             =5.0d0/3.0d0  
    double precision                    , parameter              :: perturbationInitialExponent=0.0d0
    double precision                    , parameter              :: logStabilityRatioMaximum   =60.0d0
    class           (cosmologyParametersClass)         , pointer :: thisCosmologyParameters
    class           (nodeComponentBasic)               , pointer :: thisBasic
    type            (chemicalAbundances), save                   :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                             :: shockStability       , coldFraction        , &
         &                                                          radiusShock          , coolingFunction     , &
         &                                                          densityPreShock      , densityPostShock    , &
         &                                                          numberDensityHydrogen, temperaturePostShock, &
         &                                                          velocityPreShock     , stabilityRatio

    select case (accretionMode)
    case (accretionModeTotal)
       Halo_Baryonic_Accretion_Cold_Mode_Fraction=1.0d0
    case (accretionModeHot,accretionModeCold)
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       ! Set the radiation field.
       call radiation%set(thisNode)
       ! Get the basic component.
       thisBasic => thisNode%basic()
       ! Compute factors required for stability analysis.
       radiusShock          =Dark_Matter_Halo_Virial_Radius  (thisNode)
       velocityPreShock     =Dark_Matter_Halo_Virial_Velocity(thisNode)
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
            &                  radiation                                 &
            &                 )
       coolingFunction     =                                             &
            &               Cooling_Function(                            &
            &                                temperaturePostShock ,      &
            &                                numberDensityHydrogen,      &
            &                                zeroAbundances       ,      &
            &                                chemicalDensities    ,      &
            &                                radiation                   &
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
       ! Compute the cold fraction using the model from eqn. (2) of Benson & Bower (2011).
!! AJB HACK: This original form doesn't allow the cold fraction to go to zero in high mass halos, since "shockStability" can never be less than zero.
       ! coldFraction=                                                        &
       !      &        1.0d0                                                  &
       !      &        /(                                                     &
       !      &           1.0d0                                               &
       !      &          +exp(                                                &
       !      &               (                                               &
       !      &                 accretionColdModeShockStabilityThreshold      &
       !      &                -shockStability                                &
       !      &               )                                               &
       !      &               /accretionColdModeShockStabilityTransitionWidth &
       !      &              )                                                &
       !      &         )
!! AJB HACK: This form is basically the equivalent functional form, but defined in terms of ln(epsilon) rather than epsilon.
       stabilityRatio=accretionColdModeShockStabilityThreshold/shockStability
       if (log(stabilityRatio) > accretionColdModeShockStabilityTransitionWidth*logStabilityRatioMaximum) then
          coldFraction=0.0d0
       else
          coldFraction=                                                                           &
               &        1.0d0                                                                     &
               &        /(                                                                        &
               &           1.0d0                                                                  &
               &          +stabilityRatio**(1.0d0/accretionColdModeShockStabilityTransitionWidth) &
               &         )
       end if
      ! Return the appropriate fraction.
       select case (accretionMode)
       case (accretionModeHot )
          Halo_Baryonic_Accretion_Cold_Mode_Fraction=1.0d0-coldFraction
       case (accretionModeCold)
          Halo_Baryonic_Accretion_Cold_Mode_Fraction=     +coldFraction
       end select
    end select
    return
  end function Halo_Baryonic_Accretion_Cold_Mode_Fraction

end module Accretion_Halos_Cold_Mode
