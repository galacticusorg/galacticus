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

!% Contains a module which implements calculations of baryonic accretion onto halos using a simple truncation to mimic
!% reionization.

module Accretion_Halos_Simple
  !% Implements calculations of baryonic accretion onto halos using a simple truncation to mimic reionization.
  use Radiation_Structure
  use Abundances_Structure
  implicit none
  private
  public :: Accretion_Halos_Simple_Initialize

  ! Parameters controlling when accretion is suppressed.
  double precision                     :: reionizationSuppressionRedshift             , reionizationSuppressionTime, &
       &                                  reionizationSuppressionVelocity

  ! Options controlling accretion.
  logical                              :: accretionHalosSimpleNegativeAccretionAllowed
  logical                              :: accretionHalosSimpleAccreteNewGrowthOnly

  ! Index of Solar abundance pattern.
  integer                              :: abundanceIndexSolar

  ! Internal record of the number of chemicals being tracked.
  integer                              :: chemicalsCount

  ! Radiation structure.
  type            (radiationStructure) :: radiation
  !$omp threadprivate(radiation)
contains

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_Simple_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_Simple_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
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
    use Galacticus_Nodes
    implicit none
    type            (varying_string                                    ), intent(in   )          :: accretionHalosMethod
    procedure       (Halo_Baryonic_Accretion_Rate_Simple_Get           ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Get
    procedure       (Halo_Baryonic_Accreted_Mass_Simple_Get            ), intent(inout), pointer :: Halo_Baryonic_Accreted_Mass_Get
    procedure       (Halo_Baryonic_Failed_Accretion_Rate_Simple_Get    ), intent(inout), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get
    procedure       (Halo_Baryonic_Failed_Accreted_Mass_Simple_Get     ), intent(inout), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get
    procedure       (Halo_Baryonic_Accretion_Rate_Abundances_Simple_Get), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Abundances_Get
    procedure       (Halo_Baryonic_Accreted_Abundances_Simple_Get      ), intent(inout), pointer :: Halo_Baryonic_Accreted_Abundances_Get
    procedure       (Halo_Baryonic_Accretion_Rate_Chemicals_Simple_Get ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Chemicals_Get
    procedure       (Halo_Baryonic_Accreted_Chemicals_Simple_Get       ), intent(inout), pointer :: Halo_Baryonic_Accreted_Chemicals_Get
    class           (cosmologyFunctionsClass                           )               , pointer :: cosmologyFunctionsDefault
    double precision                                                                             :: reionizationSuppressionOpticalDepth

    if (accretionHalosMethod == 'simple') then
       ! Set pointers to our implementations of accretion functions.
       Halo_Baryonic_Accretion_Rate_Get            => Halo_Baryonic_Accretion_Rate_Simple_Get
       Halo_Baryonic_Accreted_Mass_Get             => Halo_Baryonic_Accreted_Mass_Simple_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get     => Halo_Baryonic_Failed_Accretion_Rate_Simple_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get      => Halo_Baryonic_Failed_Accreted_Mass_Simple_Get
       Halo_Baryonic_Accretion_Rate_Abundances_Get => Halo_Baryonic_Accretion_Rate_Abundances_Simple_Get
       Halo_Baryonic_Accreted_Abundances_Get       => Halo_Baryonic_Accreted_Abundances_Simple_Get
       Halo_Baryonic_Accretion_Rate_Chemicals_Get  => Halo_Baryonic_Accretion_Rate_Chemicals_Simple_Get
       Halo_Baryonic_Accreted_Chemicals_Get        => Halo_Baryonic_Accreted_Chemicals_Simple_Get
       ! Get the index of the solar composition abundance pattern.
       abundanceIndexSolar=Abundance_Pattern_Lookup(abundanceName="solar")
       ! Get a count of the number of chemicals being tracked.
       chemicalsCount=Chemicals_Property_Count()
       ! Read parameters.
       if (Input_Parameter_Is_Present("reionizationSuppressionOpticalDepth")) then
          if (Input_Parameter_Is_Present("reionizationSuppressionRedshift")) call Galacticus_Error_Report("Accretion_Halos_Simple_Initialize","only one of [reionizationSuppressionOpticalDepth] and [reionizationSuppressionRedshift] should be specified")
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
          reionizationSuppressionTime=Intergalactic_Medium_Electron_Scattering_Time(reionizationSuppressionOpticalDepth&
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
          call Get_Input_Parameter("reionizationSuppressionRedshift",reionizationSuppressionRedshift,defaultValue=9.97d0)
          cosmologyFunctionsDefault => cosmologyFunctions()
          reionizationSuppressionTime=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(reionizationSuppressionRedshift))
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
       call Get_Input_Parameter("reionizationSuppressionVelocity",reionizationSuppressionVelocity,defaultValue=35.0d0)
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
       call Get_Input_Parameter("accretionHalosSimpleNegativeAccretionAllowed",accretionHalosSimpleNegativeAccretionAllowed,defaultValue=.true.)
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
       call Get_Input_Parameter("accretionHalosSimpleAccreteNewGrowthOnly",accretionHalosSimpleAccreteNewGrowthOnly,defaultValue=.false.)
       ! Define the radiation structure.
       call radiation%define([radiationTypeCMB])
       ! Check that required properties have required attributes.       
       if     (                                                                                                      &
            &   accretionHalosSimpleAccreteNewGrowthOnly                                                             &
            &  .and.                                                                                                 &
            &   .not.defaultBasicComponent%massMaximumIsGettable()                                                   &
            & ) call Galacticus_Error_Report                                                                         &
            &   (                                                                                                    &
            &    'Accretion_Halos_Simple_Initialize'                                                               , &
            &    'AccretionHalosSimpleAccreteNewGrowthOnly=true requires that the "massMaximum" '//                  &
            &    'property of the basic component be gettable.'                                  //                  &
            &    Galacticus_Component_List(                                                                          &
            &                              'basic'                                                                 , &
            &                               defaultBasicComponent%massMaximumAttributeMatch(requireGettable=.true.)  &
            &                             )                                                                          &
            &   )
    end if
    return
  end subroutine Accretion_Halos_Simple_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                ), intent(inout), pointer :: thisNode
    double precision                                                   :: growthRate             , unaccretedMass
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (nodeComponentHotHalo    )               , pointer :: thisHotHaloComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
    else
       thisBasicComponent => thisNode%basic()
       if (thisBasicComponent%time() > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
       else
          ! Get the default cosmology.
          thisCosmologyParameters => cosmologyParameters()
          thisHotHaloComponent => thisNode%hotHalo()
          Halo_Baryonic_Accretion_Rate_Simple_Get=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%accretionRate()
          ! Test for negative accretion.
          if (.not.accretionHalosSimpleNegativeAccretionAllowed.and.thisBasicComponent%accretionRate() < 0.0d0) then
             ! Accretion rate is negative, and not allowed. Return zero accretion rate.
             Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
          else
             ! Return the standard accretion rate.
             unaccretedMass=thisHotHaloComponent%unaccretedMass()
             growthRate=thisBasicComponent%accretionRate()/thisBasicComponent%mass()
             Halo_Baryonic_Accretion_Rate_Simple_Get=Halo_Baryonic_Accretion_Rate_Simple_Get+unaccretedMass*growthRate
          end if
       end if
       ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
       if (accretionHalosSimpleAccreteNewGrowthOnly) then
          if (thisBasicComponent%mass() < thisBasicComponent%massMaximum()) Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
       end if
    end if
    return
  end function Halo_Baryonic_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode                ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic      )               , pointer :: thisBasicComponent
    class(cosmologyParametersClass)               , pointer :: thisCosmologyParameters

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Accreted_Mass_Simple_Get=0.0d0
    else
       thisBasicComponent => thisNode%basic()
       if (thisBasicComponent%time() > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Accreted_Mass_Simple_Get=0.0d0
       else
          ! Get the default cosmology.
          thisCosmologyParameters => cosmologyParameters()
          Halo_Baryonic_Accreted_Mass_Simple_Get=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%mass()
       end if
    end if
    return
  end function Halo_Baryonic_Accreted_Mass_Simple_Get

  double precision function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (nodeComponentHotHalo    )               , pointer :: thisHotHaloComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    double precision                                                   :: growthRate             , unaccretedMass

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
    else
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       thisBasicComponent => thisNode%basic()
       if (thisBasicComponent%time() > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%accretionRate()
       else
          thisHotHaloComponent => thisNode%hotHalo()
          ! Test for negative accretion.
          if (.not.accretionHalosSimpleNegativeAccretionAllowed.and.thisBasicComponent%accretionRate() < 0.0d0) then
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%accretionRate()
          else
             unaccretedMass=thisHotHaloComponent%unaccretedMass()
             growthRate=thisBasicComponent%accretionRate()/thisBasicComponent%mass()
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=-unaccretedMass*growthRate
          end if
       end if
       ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
       if (accretionHalosSimpleAccreteNewGrowthOnly) then
          if (thisBasicComponent%mass() < thisBasicComponent%massMaximum()) Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
       end if
    end if
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode                ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic      )               , pointer :: thisBasicComponent
    class(cosmologyParametersClass)               , pointer :: thisCosmologyParameters

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=0.0d0
    else
       thisBasicComponent => thisNode%basic()
       if (thisBasicComponent%time() > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          ! Get the default cosmology.
          thisCosmologyParameters => cosmologyParameters()
         Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%mass()
       else
          Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=0.0d0
       end if
    end if
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get

  subroutine Halo_Baryonic_Accretion_Rate_Abundances_Simple_Get(thisNode,accretionRateAbundances)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretionRateAbundances

    ! Assume zero metallicity.
    accretionRateAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances_Simple_Get

  subroutine Halo_Baryonic_Accreted_Abundances_Simple_Get(thisNode,accretedAbundances)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretedAbundances

    ! Assume zero metallicity.
    accretedAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accreted_Abundances_Simple_Get

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Simple_Get(thisNode,accretionRateChemicals)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (chemicalAbundances), intent(inout)          :: accretionRateChemicals
    double precision                                             :: massAccretionRate

    ! Return immediately if no chemicals are being tracked.
    if (chemicalsCount == 0) return

    ! Ensure that chemicals are reset to zero.
    call accretionRateChemicals%reset()

    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=Halo_Baryonic_Accretion_Rate_Simple_Get(thisNode)

    ! Get the mass accretion rates.
    call Get_Chemical_Masses(thisNode,massAccretionRate,accretionRateChemicals)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Simple_Get

  subroutine Halo_Baryonic_Accreted_Chemicals_Simple_Get(thisNode,accretedChemicals)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (chemicalAbundances), intent(inout)          :: accretedChemicals
    double precision                                             :: massAccreted

    ! Return immediately if no chemicals are being tracked.
    if (chemicalsCount == 0) return

    ! Ensure that chemicals are reset to zero.
    call accretedChemicals%reset()

    ! Total mass of material accreted.
    massAccreted=Halo_Baryonic_Accreted_Mass_Simple_Get(thisNode)

    ! Get the masses of chemicals accreted.
    call Get_Chemical_Masses(thisNode,massAccreted,accretedChemicals)

    return
  end subroutine Halo_Baryonic_Accreted_Chemicals_Simple_Get

  subroutine Get_Chemical_Masses(thisNode,massAccreted,chemicalMasses)
    !% Compute the masses of chemicals accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    implicit none
    type            (treeNode                ), intent(inout), pointer :: thisNode
    double precision                          , intent(in   )          :: massAccreted
    type            (chemicalAbundances      ), intent(  out)          :: chemicalMasses
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    type            (chemicalAbundances      ), save                   :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                                   :: massToDensityConversion, numberDensityHydrogen, &
         &                                                                temperature

    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(thisNode))/3.0d0

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperature          =Dark_Matter_Halo_Virial_Temperature(thisNode)
    thisBasicComponent   => thisNode%basic()
    numberDensityHydrogen=hydrogenByMassPrimordial*(thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter())*thisBasicComponent%mass()*massToDensityConversion/atomicMassHydrogen

    ! Set the radiation field.
    call radiation%set(thisNode)

    ! Get the chemical densities.
    call Chemical_Densities(chemicalDensities,temperature,numberDensityHydrogen,zeroAbundances,radiation)

    ! Convert from densities to masses.
    call chemicalDensities%numberToMass(chemicalMasses)
    chemicalMasses=chemicalMasses*massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen

    return
  end subroutine Get_Chemical_Masses

end module Accretion_Halos_Simple
