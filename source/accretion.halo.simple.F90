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
  !% mimic the effects of reionization.

  use Radiation_Structure

  !# <accretionHalo name="accretionHaloSimple">
  !#  <description>Accretion onto halos using simple truncation to mimic the effects of reionization.</description>
  !# </accretionHalo>

  type, extends(accretionHaloClass) :: accretionHaloSimple
     !% A halo accretion class using simple truncation to mimic the effects of reionization.
     private
     double precision                     :: reionizationSuppressionTime, reionizationSuppressionVelocity
     logical                              :: negativeAccretionAllowed   , accreteNewGrowthOnly
     type            (radiationStructure) :: radiation
   contains
     procedure :: accretionRate          => simpleAccretionRate
     procedure :: accretedMass           => simpleAccretedMass
     procedure :: failedAccretionRate    => simpleFailedAccretionRate
     procedure :: failedAccretedMass     => simpleFailedAccretedMass
     procedure :: accretionRateMetals    => simpleAccretionRateMetals
     procedure :: accretedMassMetals     => simpleAccretedMassMetals
     procedure :: accretionRateChemicals => simpleAccretionRateChemicals
     procedure :: accretedMassChemicals  => simpleAccretedMassChemicals
     procedure :: velocityScale          => simpleVelocityScale
     procedure :: accretionRateTotal     => simpleAccretionRateTotal
     procedure :: massTotal              => simpleMassTotal
  end type accretionHaloSimple

  interface accretionHaloSimple
     !% Constructors for the {\tt simple} halo accretion class.
     module procedure simpleConstructor
     module procedure simpleDefaultConstructor
  end interface accretionHaloSimple

  ! Parameters controlling when accretion is suppressed.
  double precision                     :: simpleReionizationSuppressionRedshift, simpleReionizationSuppressionTime, &
       &                                  simpleReionizationSuppressionVelocity

  ! Options controlling accretion.
  logical                              :: simpleNegativeAccretionAllowed
  logical                              :: simpleAccreteNewGrowthOnly

  ! Index of Solar abundance pattern.
  integer                              :: simpleAbundanceIndexSolar

  ! Internal record of the number of chemicals being tracked.
  integer                              :: simpleChemicalsCount

  ! Initialization state.
  logical                              :: simpleInitialized=.false., simpleDefaultInitialized=.false.

contains

  function simpleDefaultConstructor()
    !% Default constructor for the {\tt simple} halo accretion class.
    use Intergalactic_Medium_State
    use Cosmology_Functions
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type            (accretionHaloSimple          ), target  :: simpleDefaultConstructor
    class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctionsDefault
    class           (intergalacticMediumStateClass), pointer :: intergalacticMediumState_ 
    double precision                                         :: reionizationSuppressionOpticalDepth

    ! Get default parameters.
    if (.not.simpleDefaultInitialized) then
       !$omp critical(accretionHaloSimpleDefaultInitialize)
       if (.not.simpleDefaultInitialized) then
          if (Input_Parameter_Is_Present("simpleReionizationSuppressionOpticalDepth")) then
             if (Input_Parameter_Is_Present("simpleReionizationSuppressionRedshift")) call Galacticus_Error_Report("simpleDefaultConstructor","only one of [simpleReionizationSuppressionOpticalDepth] and [simpleReionizationSuppressionRedshift] should be specified")
             !@ <inputParameter>
             !@   <name>reionizationSuppressionOpticalDepth</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@    The optical depth to electron scattering below which baryonic accretion is suppressed.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("reionizationSuppressionOpticalDepth",reionizationSuppressionOpticalDepth)
             intergalacticMediumState_        => intergalacticMediumState()
             simpleReionizationSuppressionTime=intergalacticMediumState_%electronScatteringTime(reionizationSuppressionOpticalDepth&
                  &,assumeFullyIonized=.true.)
          else
             !@ <inputParameter>
             !@   <name>reionizationSuppressionRedshift</name>
             !@   <defaultValue>9.97 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@    The redshift below which baryonic accretion is suppressed.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("reionizationSuppressionRedshift",simpleReionizationSuppressionRedshift,defaultValue=9.97d0)
             cosmologyFunctionsDefault        => cosmologyFunctions()
             simpleReionizationSuppressionTime=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(simpleReionizationSuppressionRedshift))
          end if
          !@ <inputParameter>
          !@   <name>reionizationSuppressionVelocity</name>
          !@   <defaultValue>35.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The velocity scale below which baryonic accretion is suppressed.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("reionizationSuppressionVelocity",simpleReionizationSuppressionVelocity,defaultValue=35.0d0)
          !@ <inputParameter>
          !@   <name>accretionHalosSimpleNegativeAccretionAllowed</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether negative accretion (mass loss) is allowed in the simple halo accretion model.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("accretionHalosSimpleNegativeAccretionAllowed",simpleNegativeAccretionAllowed,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>accretionHalosSimpleAccreteNewGrowthOnly</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether accretion from the \gls{igm} is allowed only when a halo is growing past its previous greatest mass.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("accretionHalosSimpleAccreteNewGrowthOnly",simpleAccreteNewGrowthOnly,defaultValue=.false.)
          ! Record that class is now initialized.
          simpleDefaultInitialized=.true.
       end if
       !$omp end critical(accretionHaloSimpleDefaultInitialize)
    end if
    simpleDefaultConstructor=simpleConstructor(simpleReionizationSuppressionTime,simpleReionizationSuppressionVelocity,simpleNegativeAccretionAllowed,simpleAccreteNewGrowthOnly)
    return
  end function simpleDefaultConstructor
       
  function simpleConstructor(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    !% Default constructor for the {\tt simple} halo accretion class.
    use Input_Parameters
    use Galacticus_Nodes
    use Galacticus_Error
    use Atomic_Data
    use Chemical_Abundances_Structure
    implicit none
    type            (accretionHaloSimple), target        :: simpleConstructor
    double precision                     , intent(in   ) :: reionizationSuppressionTime, reionizationSuppressionVelocity
    logical                              , intent(in   ) :: negativeAccretionAllowed   , accreteNewGrowthOnly

    simpleConstructor%reionizationSuppressionTime    =reionizationSuppressionTime
    simpleConstructor%reionizationSuppressionVelocity=reionizationSuppressionVelocity
    simpleConstructor%negativeAccretionAllowed       =negativeAccretionAllowed
    simpleConstructor%accreteNewGrowthOnly           =accreteNewGrowthOnly
    call simpleConstructor%radiation%define([radiationTypeCMB])
    ! Check that required properties have required attributes.       
    if     (                                                                                                      &
         &   accreteNewGrowthOnly                                                                                 &
         &  .and.                                                                                                 &
         &   .not.defaultBasicComponent%massMaximumIsGettable()                                                   &
         & ) call Galacticus_Error_Report                                                                         &
         &   (                                                                                                    &
         &    'simpleConstructor'                                                                               , &
         &    'accreteNewGrowthOnly=true requires that the "massMaximum" '//                                      &
         &    'property of the basic component be gettable.'              //                                      &
         &    Galacticus_Component_List(                                                                          &
         &                              'basic'                                                                 , &
         &                               defaultBasicComponent%massMaximumAttributeMatch(requireGettable=.true.)  &
         &                             )                                                                          &
         &   )
    ! Perform class global initialization.
    if (.not.simpleInitialized) then
       !$omp critical(accretionHaloSimpleInitialize)
       if (.not.simpleInitialized) then
          simpleAbundanceIndexSolar=Abundance_Pattern_Lookup(abundanceName="solar")
          ! Get a count of the number of chemicals being tracked.
          simpleChemicalsCount     =Chemicals_Property_Count(                     )
          ! Record that initialization is completed.
          simpleInitialized=.true.
       end if
       !$omp end critical(accretionHaloSimpleInitialize)
    end if
    return
  end function simpleConstructor

  double precision function simpleAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    implicit none
    class           (accretionHaloSimple     ), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: node
    integer                                   , intent(in   )          :: accretionMode
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (nodeComponentHotHalo    )               , pointer :: thisHotHaloComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    double precision                                                   :: growthRate             , unaccretedMass

    simpleAccretionRate=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    thisBasicComponent   => node%basic()
    if (thisBasicComponent%time() > self%reionizationSuppressionTime .and. self%velocityScale(node) <&
         & self%reionizationSuppressionVelocity) then
       simpleAccretionRate=0.0d0
    else
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       thisHotHaloComponent => node%hotHalo()
       simpleAccretionRate=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%accretionRateTotal(node)
       ! Test for negative accretion.
       if (.not.self%negativeAccretionAllowed.and.self%accretionRateTotal(node) < 0.0d0) then
          ! Accretion rate is negative, and not allowed. Return zero accretion rate.
          simpleAccretionRate=0.0d0
       else
          ! Return the standard accretion rate.
          unaccretedMass=thisHotHaloComponent%unaccretedMass()
          growthRate=self%accretionRateTotal(node)/self%massTotal(node)
          simpleAccretionRate=simpleAccretionRate+unaccretedMass*growthRate
       end if
       ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
       if (self%accreteNewGrowthOnly) then
          if (self%massTotal(node) < thisBasicComponent%massMaximum()) simpleAccretionRate=0.0d0
       end if
    end if
    return
  end function simpleAccretionRate

  double precision function simpleAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\tt node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    implicit none
    class (accretionHaloSimple    ), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node
    integer                        , intent(in   )          :: accretionMode
    class(nodeComponentBasic      )               , pointer :: thisBasicComponent
    class(cosmologyParametersClass)               , pointer :: thisCosmologyParameters

    simpleAccretedMass=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    thisBasicComponent   => node%basic     ()
    if (thisBasicComponent%time() > self%reionizationSuppressionTime .and. self%velocityScale(node) <&
         & self%reionizationSuppressionVelocity) then
       simpleAccretedMass=0.0d0
    else
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       simpleAccretedMass=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%massTotal(node)
    end if
    return
  end function simpleAccretedMass

  double precision function simpleFailedAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\tt node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    class           (accretionHaloSimple     ), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: node
    integer                                   , intent(in   )          :: accretionMode
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (nodeComponentHotHalo    )               , pointer :: thisHotHaloComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    double precision                                                   :: growthRate             , unaccretedMass

    simpleFailedAccretionRate=0.0d0
    if (accretionMode          == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
       thisBasicComponent      => node%basic     ()
       if (thisBasicComponent%time() > self%reionizationSuppressionTime .and. self%velocityScale(node) <&
         & self%reionizationSuppressionVelocity) then
          simpleFailedAccretionRate=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%accretionRateTotal(node)
    else
       thisHotHaloComponent => node%hotHalo()
       ! Test for negative accretion.
       if (.not.self%negativeAccretionAllowed.and.self%accretionRateTotal(node) < 0.0d0) then
          simpleFailedAccretionRate=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%accretionRateTotal(node)
       else
          unaccretedMass=thisHotHaloComponent%unaccretedMass()
          growthRate=self%accretionRateTotal(node)/self%massTotal(node)
          simpleFailedAccretionRate=-unaccretedMass*growthRate
       end if
       ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
       if (self%accreteNewGrowthOnly) then
          if (self%massTotal(node) < thisBasicComponent%massMaximum()) simpleFailedAccretionRate=0.0d0
       end if
    end if
    return
  end function simpleFailedAccretionRate

  double precision function simpleFailedAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\tt node}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    class  (accretionHaloSimple     ), intent(inout)          :: self
    type   (treeNode                ), intent(inout), pointer :: node
    integer                          , intent(in   )          :: accretionMode
    class  (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class  (cosmologyParametersClass)               , pointer :: thisCosmologyParameters

    simpleFailedAccretedMass=0.0d0
    if (accretionMode          == accretionModeCold) return
    if (node%isSatellite()                     ) return
    thisBasicComponent   => node%basic     ()
    if (thisBasicComponent%time() > self%reionizationSuppressionTime .and. self%velocityScale(node) <&
         & self%reionizationSuppressionVelocity) then
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       simpleFailedAccretedMass=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%massTotal(node)
    end if
    return
  end function simpleFailedAccretedMass
  
  function simpleAccretionRateMetals(self,node,accretionMode)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type  (abundances         )                         :: simpleAccretionRateMetals
    class (accretionHaloSimple), intent(inout)          :: self
    type  (treeNode           ), intent(inout), pointer :: node
    integer                    , intent(in   )          :: accretionMode

    ! Assume zero metallicity.
    simpleAccretionRateMetals=zeroAbundances
    return
  end function simpleAccretionRateMetals

  function simpleAccretedMassMetals(self,node,accretionMode)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type   (abundances         )                         :: simpleAccretedMassMetals
    class  (accretionHaloSimple), intent(inout)          :: self
    type   (treeNode           ), intent(inout), pointer :: node
    integer                     , intent(in   )          :: accretionMode
    
    ! Assume zero metallicity.
    simpleAccretedMassMetals=zeroAbundances
    return
  end function simpleAccretedMassMetals

  function simpleAccretionRateChemicals(self,node,accretionMode)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\tt node} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances )                         :: simpleAccretionRateChemicals
    class           (accretionHaloSimple), intent(inout)          :: self
    type            (treeNode           ), intent(inout), pointer :: node
    integer                              , intent(in   )          :: accretionMode
    double precision                                              :: massAccretionRate

    ! Return immediately if no chemicals are being tracked.
    if (simpleChemicalsCount == 0) return

    ! Ensure that chemicals are reset to zero.
    call simpleAccretionRateChemicals%reset()

    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=simpleAccretionRate(self,node,accretionMode)

    ! Get the mass accretion rates.
    simpleAccretionRateChemicals=simpleChemicalMasses(self,node,massAccretionRate)

    return
  end function simpleAccretionRateChemicals

  function simpleAccretedMassChemicals(self,node,accretionMode)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances )                         :: simpleAccretedMassChemicals
    class           (accretionHaloSimple), intent(inout)          :: self
    type            (treeNode           ), intent(inout), pointer :: node
    integer                              , intent(in   )          :: accretionMode
    double precision                                              :: massAccreted

    ! Return immediately if no chemicals are being tracked.
    if (simpleChemicalsCount == 0) return

    ! Ensure that chemicals are reset to zero.
    call simpleAccretedMassChemicals%reset()

    ! Total mass of material accreted.
    massAccreted=simpleAccretedMass(self,node,accretionMode)

    ! Get the masses of chemicals accreted.
    simpleAccretedMassChemicals=simpleChemicalMasses(self,node,massAccreted)

    return
  end function simpleAccretedMassChemicals

  function simpleChemicalMasses(self,node,massAccreted)
    !% Compute the masses of chemicals accreted (in $M_\odot$) onto {\tt node} from the intergalactic medium.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    implicit none
    class           (accretionHaloSimple     ), intent(inout)          :: self
    type            (chemicalAbundances      )                         :: simpleChemicalMasses
    type            (treeNode                ), intent(inout), pointer :: node
    double precision                          , intent(in   )          :: massAccreted
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    type            (chemicalAbundances      ), save                   :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                                   :: massToDensityConversion, numberDensityHydrogen, &
         &                                                                temperature

    ! Compute coefficient in conversion of mass to density for this node.
    darkMatterHaloScale_ => darkMatterHaloScale()
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(darkMatterHaloScale_%virialRadius(node))/3.0d0
    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperature          =darkMatterHaloScale_%virialTemperature(node)
    thisBasicComponent   => node%basic()
    numberDensityHydrogen=hydrogenByMassPrimordial*(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*self%massTotal(node)*massToDensityConversion/atomicMassHydrogen
    ! Set the radiation field.
    call self%radiation%set(node)
    ! Get the chemical densities.
    call Chemical_Densities(chemicalDensities,temperature,numberDensityHydrogen,zeroAbundances,self%radiation)
    ! Convert from densities to masses.
    call chemicalDensities%numberToMass(simpleChemicalMasses)
    simpleChemicalMasses=simpleChemicalMasses*massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen
    return
  end function simpleChemicalMasses

  double precision function simpleVelocityScale(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class(accretionHaloSimple     ), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node
    class(darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_

    darkMatterHaloScale_ => darkMatterHaloScale                (    )
    simpleVelocityScale  =  darkMatterHaloScale_%virialVelocity(node)
    return
  end function simpleVelocityScale

  double precision function simpleAccretionRateTotal(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloSimple     ), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node
    class(nodeComponentBasic      )               , pointer :: basic

    basic                    => node %basic        ()
    simpleAccretionRateTotal =  basic%accretionRate()
    return
  end function simpleAccretionRateTotal

  double precision function simpleMassTotal(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloSimple     ), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node
    class(nodeComponentBasic      )               , pointer :: basic

    basic           => node %basic()
    simpleMassTotal =  basic%mass ()
    return
  end function simpleMassTotal
