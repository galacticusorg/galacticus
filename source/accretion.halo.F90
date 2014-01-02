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

!% Contains a module which implements calculations of baryonic accretion into halos.

module Accretion_Halos
  !% Implements calculations of baryonic accretion into halos.
  use ISO_Varying_String
  use Galacticus_Nodes
  use Abundances_Structure
  use Chemical_Abundances_Structure
  implicit none
  private
  public :: Halo_Baryonic_Accretion_Rate, Halo_Baryonic_Accreted_Mass, Halo_Baryonic_Failed_Accretion_Rate,&
       & Halo_Baryonic_Failed_Accreted_Mass, Halo_Baryonic_Accretion_Rate_Abundances, Halo_Baryonic_Accreted_Abundances,&
       & Halo_Baryonic_Accretion_Rate_Chemicals, Halo_Baryonic_Accreted_Chemicals

  ! Flag to indicate if this module has been initialized.
  logical                                                   :: accretionHalosInitialized                  =.false.

  ! Name of mass movement method used.
  type     (varying_string                       )          :: accretionHalosMethod

  ! Pointers to functions that return baryonic mass accretion rates/masses.
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Accretion_Rate_Get           =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Accreted_Mass_Get            =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get    =>null()
  procedure(Halo_Baryonic_Accretion_Get_Template ), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get     =>null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer :: Halo_Baryonic_Accreted_Abundances_Get      =>null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer :: Halo_Baryonic_Accreted_Chemicals_Get       =>null()
  procedure(Halo_Baryonic_Abundances_Get_Template), pointer :: Halo_Baryonic_Accretion_Rate_Abundances_Get=>null()
  procedure(Halo_Baryonic_Chemicals_Get_Template ), pointer :: Halo_Baryonic_Accretion_Rate_Chemicals_Get =>null()
  abstract interface
     double precision function Halo_Baryonic_Accretion_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Baryonic_Accretion_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Abundances_Get_Template(thisNode,accretedAbundances)
       import treeNode
       import abundances
       type(treeNode  ), intent(inout), pointer :: thisNode
       type(abundances), intent(inout)          :: accretedAbundances
     end subroutine Halo_Baryonic_Abundances_Get_Template
  end interface
  abstract interface
     subroutine Halo_Baryonic_Chemicals_Get_Template(thisNode,accretedChemicals)
       import treeNode
       import chemicalAbundances
       type(treeNode          ), intent(inout), pointer :: thisNode
       type(chemicalAbundances), intent(inout)          :: accretedChemicals
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
          !# <include directive="accretionHalosMethod" type="functionCall" functionType="void">
          !#  <functionArgs>accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get,Halo_Baryonic_Accreted_Abundances_Get,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accretion_Rate_Chemicals_Get,Halo_Baryonic_Accreted_Chemicals_Get</functionArgs>
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
    type(treeNode), intent(inout), pointer :: thisNode

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
    type(treeNode), intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Accreted_Mass=Halo_Baryonic_Accreted_Mass_Get(thisNode)

    return
  end function Halo_Baryonic_Accreted_Mass

  double precision function Halo_Baryonic_Failed_Accretion_Rate(thisNode)
    !% Computes the rate of failed baryonic mass accretion (in $M_\odot$/Gyr) onto {\tt thisNode} from the intergalactic medium.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

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
    type(treeNode), intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accreted_Mass=Halo_Baryonic_Failed_Accreted_Mass_Get(thisNode)

    return
  end function Halo_Baryonic_Failed_Accreted_Mass

  subroutine Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances)
    !% Compute the rate of mass accretion of abundances (in $M_\odot/$Gyr) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretionRateAbundances

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Abundances_Get(thisNode,accretionRateAbundances)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances

  subroutine Halo_Baryonic_Accreted_Abundances(thisNode,accretedAbundances)
    !% Compute the mass of abundances (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretedAbundances

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Abundances_Get(thisNode,accretedAbundances)

    return
  end subroutine Halo_Baryonic_Accreted_Abundances

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals(thisNode,accretionRateChemicals)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretionRateChemicals

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accretion_Rate_Chemicals_Get(thisNode,accretionRateChemicals)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals

  subroutine Halo_Baryonic_Accreted_Chemicals(thisNode,accretedChemicals)
    !% Compute the mass of chemicals (in $M_\odot$) accreted onto {\tt thisNode}.
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretedChemicals

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    call Halo_Baryonic_Accreted_Chemicals_Get(thisNode,accretedChemicals)

    return
  end subroutine Halo_Baryonic_Accreted_Chemicals

end module Accretion_Halos
