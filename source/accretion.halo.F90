!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  use Tree_Nodes
  private
  public :: Halo_Baryonic_Accretion_Rate, Halo_Baryonic_Accreted_Mass, Halo_Baryonic_Failed_Accretion_Rate,&
       & Halo_Baryonic_Failed_Accreted_Mass

  ! Flag to indicate if this module has been initialized.  
  logical                                                       :: accretionHalosInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                                          :: accretionHalosMethod

  ! Pointers to functions that return baryonic mass accretion rates/masses.
  procedure(Halo_Baryonic_Accretion_Get_Template), pointer :: Halo_Baryonic_Accretion_Rate_Get        => null()
  procedure(Halo_Baryonic_Accretion_Get_Template), pointer :: Halo_Baryonic_Accreted_Mass_Get         => null()
  procedure(Halo_Baryonic_Accretion_Get_Template), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get => null()
  procedure(Halo_Baryonic_Accretion_Get_Template), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get  => null()

  abstract interface
     double precision function Halo_Baryonic_Accretion_Get_Template(thisNode)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
     end function Halo_Baryonic_Accretion_Get_Template
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

    !$omp critical(accretionHalosInitialize)
    if (.not.accretionHalosInitialized) then
       ! Do the binary black hole merger method parameter.
       !@ <inputParameter>
       !@   <name>accretionHalosMethod</name>
       !@   <defaultValue>simple</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Selects which method should be used for accretion onto halos.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('accretionHalosMethod',accretionHalosMethod,defaultValue='simple')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="accretionHalosMethod" type="code" action="subroutine">
       !#  <subroutineArgs>accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get</subroutineArgs>
       include 'accretion.halos.inc'
       !# </include>
       if (.not.(associated(Halo_Baryonic_Accretion_Rate_Get).and.associated(Halo_Baryonic_Accreted_Mass_Get) &
            .and.associated(Halo_Baryonic_Failed_Accretion_Rate_Get).and.associated(Halo_Baryonic_Failed_Accreted_Mass_Get))) &
            & call Galacticus_Error_Report('Accretion_Halos_Initialize','method ' //char(accretionHalosMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       accretionHalosInitialized=.true.
    end if
    !$omp end critical(accretionHalosInitialized)

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
    !% Computes the mass of baryons that failed to accre (in $M_\odot$) into {\tt thisNode} from the intergalactic medium. Used to initialize
    !% nodes.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Accretion_Halos_Initialize

    ! Get the accretion rate.
    Halo_Baryonic_Failed_Accreted_Mass=Halo_Baryonic_Failed_Accreted_Mass_Get(thisNode)

    return
  end function Halo_Baryonic_Failed_Accreted_Mass

end module Accretion_Halos
