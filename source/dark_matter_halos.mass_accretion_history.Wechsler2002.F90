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

!% Contains a module which implements the \cite{wechsler_concentrations_2002} halo mass accretion algorithm.

module Dark_Matter_Halo_Mass_Accretion_Histories_Wechsler2002
  !% Implements the \cite{wechsler_concentrations_2002} halo mass accretion algorithm.
  use Galacticus_Nodes
  implicit none
  private
  public :: Dark_Matter_Mass_Accretion_Wechsler2002_Initialize

  ! Parameters controlling the calculation of formation histories.
  logical          :: accretionHistoryWechslerFormationRedshiftCompute
  double precision :: accretionHistoryWechslerFormationRedshift

contains

  !# <darkMatterAccretionHistoryMethod>
  !#  <unitName>Dark_Matter_Mass_Accretion_Wechsler2002_Initialize</unitName>
  !# </darkMatterAccretionHistoryMethod>
  subroutine Dark_Matter_Mass_Accretion_Wechsler2002_Initialize(darkMatterAccretionHistoryMethod,Dark_Matter_Halo_Mass_Accretion_Time_Get)
    !% Initializes the ``Wechsler2002'' mass accretion history module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterAccretionHistoryMethod
    procedure(Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002), pointer, intent(inout) :: Dark_Matter_Halo_Mass_Accretion_Time_Get
    
    if (darkMatterAccretionHistoryMethod == 'Wechsler2002') then
       ! Set procedure pointers.
       Dark_Matter_Halo_Mass_Accretion_Time_Get => Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002

       !@ <inputParameter>
       !@   <name>accretionHistoryWechslerFormationRedshiftCompute</name>
       !@   <defaultValue>true</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Compute formation redshift automatically for \cite{wechsler_concentrations_2002} halo mass accretion histories?
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('accretionHistoryWechslerFormationRedshiftCompute',accretionHistoryWechslerFormationRedshiftCompute,defaultValue=.true.)
       if (.not.accretionHistoryWechslerFormationRedshiftCompute) then
          ! In this case, read the formation redshift.
          !@ <inputParameter>
          !@   <name>accretionHistoryWechslerFormationRedshift</name>
          !@   <defaultValue>0.4</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The formation redshift to use in \cite{wechsler_concentrations_2002} halo mass accretion histories.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('accretionHistoryWechslerFormationRedshift',accretionHistoryWechslerFormationRedshift,defaultValue=0.4d0)
        end if
       
    end if

    return
  end subroutine Dark_Matter_Mass_Accretion_Wechsler2002_Initialize

  double precision function Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002(baseNode,nodeMass)
    !% Compute the time corresponding to {\tt nodeMass} in the mass accretion history of {\tt thisNode} using the algorithm of
    !% \cite{wechsler_concentrations_2002}.
    use Cosmology_Functions
    implicit none
    type (treeNode          ), intent(inout), pointer :: baseNode
    double precision,          intent(in   )          :: nodeMass
    class(nodeComponentBasic),                pointer :: baseBasicComponent
    double precision                                  :: expansionFactorBase,expansionFactor,mergerTreeFormationExpansionFactor

    baseBasicComponent => baseNode%basic()
    select case (accretionHistoryWechslerFormationRedshiftCompute)
    case (.true.)
       ! Compute the expansion factor at formation.
       mergerTreeFormationExpansionFactor=Expansion_Factor_At_Formation (baseBasicComponent%mass()                )
    case (.false.)
       ! Use the specified formation redshift.
       mergerTreeFormationExpansionFactor=Expansion_Factor_from_Redshift(accretionHistoryWechslerFormationRedshift)
    end select
    
    ! Get the expansion factor at the tree base.
    expansionFactorBase=Expansion_Factor(baseBasicComponent%time())

    ! Compute the expansion factor for the current node.
    expansionFactor    =expansionFactorBase/(1.0d0-0.5d0*dlog(nodeMass/baseBasicComponent%mass())&
         &/mergerTreeFormationExpansionFactor)

    ! Find the time corresponding to this expansion factor.
    Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002=Cosmology_Age(expansionFactor)

   return
  end function Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002
  
  double precision function Expansion_Factor_At_Formation(haloMass)
    !% Computes the expansion factor at formation using the simple model of \cite{bullock_profiles_2001}.
    use Cosmology_Functions
    use Power_Spectra
    use Critical_Overdensity
    implicit none
    double precision, intent(in) :: haloMass
    double precision, parameter  :: haloMassFraction=0.015d0 ! Wechsler et al. (2002;  Astrophysical Journal, 568:52-70).
    double precision             :: haloMassCharacteristic,sigmaCharacteristic,formationTime

    ! Compute the characteristic mass at formation time.    
    haloMassCharacteristic=haloMassFraction*haloMass

    ! Compute the corresponding rms fluctuation in the density field (i.e. sigma(M)).
    sigmaCharacteristic=Cosmological_Mass_Root_Variance(haloMassCharacteristic)

    ! Get the time at which this equals the critical overdensity for collapse.
    formationTime=Time_of_Collapse(criticalOverdensity=sigmaCharacteristic,mass=haloMass)
    
    ! Get the corresponding expansion factor.
    Expansion_Factor_At_Formation=Expansion_Factor(formationTime)
    
    return
  end function Expansion_Factor_At_Formation

end module Dark_Matter_Halo_Mass_Accretion_Histories_Wechsler2002
