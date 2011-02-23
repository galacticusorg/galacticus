!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements the \cite{wechsler_concentrations_2002} halo mass accretion algorithm.

module Dark_Matter_Halo_Mass_Accretion_Histories_Wechsler2002
  !% Implements the \cite{wechsler_concentrations_2002} halo mass accretion algorithm.
  use Tree_Nodes
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
    !% Initializes the ``Wechsler 2002'' mass accretion history module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterAccretionHistoryMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Halo_Mass_Accretion_Time_Get
    
    if (darkMatterAccretionHistoryMethod == 'Wechsler 2002') then
       ! Set procedure pointers.
       Dark_Matter_Halo_Mass_Accretion_Time_Get => Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002

       !@ <inputParameter>
       !@   <name>accretionHistoryWechslerFormationRedshiftCompute</name>
       !@   <defaultValue>true</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Compute formation redshift automatically for \cite{wechsler_concentrations_2002} halo mass accretion histories?
       !@   </description>
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
          !@ </inputParameter>
          call Get_Input_Parameter('accretionHistoryWechslerFormationRedshift',accretionHistoryWechslerFormationRedshift,defaultValue=0.4d0)
        end if
       
    end if

    return
  end subroutine Dark_Matter_Mass_Accretion_Wechsler2002_Initialize

  double precision function Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002(baseNode,nodeMass)
    !% Compute the time corresponding to {\tt nodeMass} in the mass accretion history of {\tt thisNode} using the algorithm of
    !% \cite{wechsler_concentrations_2002}.
    use Tree_Nodes
    use Cosmology_Functions
    implicit none
    type(treeNode),   intent(inout), pointer :: baseNode
    double precision, intent(in)             :: nodeMass
    double precision                         :: expansionFactorBase,expansionFactor,mergerTreeFormationExpansionFactor

    select case (accretionHistoryWechslerFormationRedshiftCompute)
    case (.true.)
       ! Compute the expansion factor at formation.
       mergerTreeFormationExpansionFactor=Expansion_Factor_At_Formation(Tree_Node_Mass(baseNode))
    case (.false.)
       ! Use the specified formation redshift.
       mergerTreeFormationExpansionFactor=Expansion_Factor_from_Redshift(accretionHistoryWechslerFormationRedshift)
    end select
    
    ! Get the expansion factor at the tree base.
    expansionFactorBase=Expansion_Factor(Tree_Node_Time(baseNode))

    ! Compute the expansion factor for the current node.
    expansionFactor    =expansionFactorBase/(1.0d0-0.5d0*dlog(nodeMass/Tree_Node_Mass(baseNode))&
         &/mergerTreeFormationExpansionFactor)

    ! Find the time corresponding to this expansion factor.
    Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002=Cosmology_Age(expansionFactor)

   return
  end function Dark_Matter_Halo_Mass_Accretion_Time_Wechsler2002
  
  double precision function Expansion_Factor_At_Formation(haloMass)
    !% Computes the expansion factor at formation using the simple model of \cite{bullock_profiles_2001}.
    use Cosmology_Functions
    use CDM_Power_Spectrum
    use Critical_Overdensity
    implicit none
    double precision, intent(in) :: haloMass
    double precision, parameter  :: haloMassFraction=0.015d0 ! Wechsler et al. (2002;  Astrophysical Journal, 568:52-70).
    double precision             :: haloMassCharacteristic,sigmaCharacteristic,formationTime

    ! Compute the characteristic mass at formation time.    
    haloMassCharacteristic=haloMassFraction*haloMass

    ! Compute the corresponding rms fluctuation in the density field (i.e. sigma(M)).
    sigmaCharacteristic=sigma_CDM(haloMassCharacteristic)

    ! Get the time at which this equals the critical overdensity for collapse.
    formationTime=Time_of_Collapse(sigmaCharacteristic)
    
    ! Get the corresponding expansion factor.
    Expansion_Factor_At_Formation=Expansion_Factor(formationTime)
    
    return
  end function Expansion_Factor_At_Formation

end module Dark_Matter_Halo_Mass_Accretion_Histories_Wechsler2002
