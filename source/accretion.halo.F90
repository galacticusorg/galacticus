!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of baryonic accretion into halos.

module Accretion_Halos
  !% Implements calculations of baryonic accretion into halos.
  use ISO_Varying_String
  use Tree_Nodes
  use Abundances_Structure
  use Chemical_Abundances_Structure
  implicit none
  private
  public :: Halo_Baryonic_Accretion_Rate, Halo_Baryonic_Accreted_Mass, Halo_Baryonic_Failed_Accretion_Rate,&
       & Halo_Baryonic_Failed_Accreted_Mass, Halo_Baryonic_Accretion_Rate_Abundances, Halo_Baryonic_Accreted_Abundances,&
       & Halo_Baryonic_Accretion_Rate_Chemicals, Halo_Baryonic_Accreted_Chemicals

  ! Flag to indicate if this module has been initialized.  
  logical                                                       :: accretionHalosInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                                          :: accretionHalosMethod

  ! Pointers to functions that return baryonic mass accretion rates/masses.
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Accretion_Rate_Get            => null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Accreted_Mass_Get             => null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get     => null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get      => null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer :: Halo_Baryonic_Accreted_Abundances_Get       => null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer :: Halo_Baryonic_Accreted_Chemicals_Get        => null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer :: Halo_Baryonic_Accretion_Rate_Abundances_Get => null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer :: Halo_Baryonic_Accretion_Rate_Chemicals_Get  => null()
  abstract interface
     double precision function Halo_Baryonic_Accretion_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Baryonic_Accretion_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Abundances_Get_Template(thisNode,accretedAbundances)
       import treeNode
       import abundancesStructure
       type(treeNode),            intent(inout), pointer :: thisNode
       type(abundancesStructure), intent(inout)          :: accretedAbundances
     end subroutine Halo_Baryonic_Abundances_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Chemicals_Get_Template(thisNode,accretedChemicals)
       import treeNode
       import chemicalAbundancesStructure
       type(treeNode),                    intent(inout), pointer :: thisNode
       type(chemicalAbundancesStructure), intent(inout)          :: accretedChemicals
     end subroutine Halo_Baryonic_Chemicals_Get_Template
  end interface

contains

  subroutine Accretion_Halos_Initialize
    !% Initalize the accretion disk module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="accretionHalosMethod" type="moduleUse">
    include 'accretion.halos.modules.inc'
    !# </include>
    implicit none

    if (.not.accretionHalosInitialized) then
       !$omp critical(accretionHalosInitialize)
       if (.not.accretionHalosInitialized) then
          ! Get the halo accretion method parameter.
          !@ <inputParameter>
          !@   <name>accretionHalosMethod</name>
          !@   <defaultValue>simple</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Selects which method should be used for accretion onto halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('accretionHalosMethod',accretionHalosMethod,defaultValue='simple')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="accretionHalosMethod" type="code" action="subroutine">
          !#  <subroutineArgs>accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get,Halo_Baryonic_Accreted_Abundances_Get,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accretion_Rate_Chemicals_Get,Halo_Baryonic_Accreted_Chemicals_Get</subroutineArgs>
          include 'accretion.halos.inc'
          !# </include>
          if     (.not.(     associated(Halo_Baryonic_Accretion_Rate_Get           ) &
               &        .and.associated(Halo_Baryonic_Accreted_Mass_Get            ) &
               &        .and.associated(Halo_Baryonic_Failed_Accretion_Rate_Get    ) &
               &        .and.associated(Halo_Baryonic_Failed_Accreted_Mass_Get     ) &
               &        .and.associated(Halo_Baryonic_Accretion_Rate_Abundances_Get) &
               &        .and.associated(Halo_Baryonic_Accreted_Abundances_Get      ) &
               &        .and.associated(Halo_Baryonic_Accretion_Rate_Chemicals_Get ) &
               &        .and.associated(Halo_Baryonic_Accreted_Chemicals_Get       ) &
               &       )                                                             &
               & ) call Galacticus_Error_Report('Accretion_Halos_Initialize','method ' //char(accretionHalosMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          accretionHalosInitialized=.true.
       end if
       !$omp end critical(accretionHalosInitialize)
    end if
    return
  end subroutine Accretion_Halos_Initialize
  
  double precision function Halo_Baryonic_Accretion_Rate(thisNode)
    !% Computes the rate of baryonic mass accretion (in $M_\odot$/Gyr) onto {\tt thisNode} from the intergalactic medium.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Accretion_Rate=Halo_Baryonic_Accretion_Rate_Get(thisNode)

    return
  end function Halo_Baryonic_Accretion_Rate

  double precision function Halo_Baryonic_Accreted_Mass(thisNode)
    !% Computes the mass of baryons accreted (in $M_\odot$) into {\tt thisNode} from the intergalactic medium. Used to initialize
    !% nodes.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Accreted_Mass=Halo_Baryonic_Accreted_Mass_Get(thisNode)

    return
  end function Halo_Baryonic_Accreted_Mass

  double precision function Halo_Baryonic_Failed_Accretion_Rate(thisNode)
    !% Computes the rate of failed baryonic mass accretion (in $M_\odot$/Gyr) onto {\tt thisNode} from the intergalactic medium.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accretion_Rate=Halo_Baryonic_Failed_Accretion_Rate_Get(thisNode)

    return
  end function Halo_Baryonic_Failed_Accretion_Rate

  double precision function Halo_Baryonic_Failed_Accreted_Mass(thisNode)
    !% Computes the mass of baryons that failed to accrete (in $M_\odot$) into {\tt thisNode} from the intergalactic medium. Used to initialize
    !% nodes.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accreted_Mass=Halo_Baryonic_Failed_Accreted_Mass_Get(thisNode)

    return
  end function Halo_Baryonic_Failed_Accreted_Mass

  subroutine Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances)
    !% Compute the rate of mass accretion of abundances (in $M_\odot/$Gyr) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(abundancesStructure), intent(inout)          :: accretionRateAbundances

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Abundances_Get(thisNode,accretionRateAbundances)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances

  subroutine Halo_Baryonic_Accreted_Abundances(thisNode,accretedAbundances)
    !% Compute the mass of abundances (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(abundancesStructure), intent(inout)          :: accretedAbundances

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Abundances_Get(thisNode,accretedAbundances)

    return
  end subroutine Halo_Baryonic_Accreted_Abundances

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals(thisNode,accretionRateChemicals)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode),                    intent(inout), pointer :: thisNode
    type(chemicalAbundancesStructure), intent(inout)          :: accretionRateChemicals

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Chemicals_Get(thisNode,accretionRateChemicals)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals

  subroutine Halo_Baryonic_Accreted_Chemicals(thisNode,accretedChemicals)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode),                    intent(inout), pointer :: thisNode
    type(chemicalAbundancesStructure), intent(inout)          :: accretedChemicals

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Chemicals_Get(thisNode,accretedChemicals)

    return
  end subroutine Halo_Baryonic_Accreted_Chemicals

end module Accretion_Halos
