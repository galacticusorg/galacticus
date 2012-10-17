!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of mass loss rates from dark matter halos.

module Dark_Matter_Halos_Mass_Loss_Rates
  !% Implements calculations of mass loss rates from dark matter halos.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Dark_Matter_Halos_Mass_Loss_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterHaloMassLossRateInitialized=.false.

  ! Name of mass loss rate method used.
  type(varying_string) :: darkMatterHaloMassLossRateMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Halo_Mass_Loss_Rate_Template), pointer :: Dark_Matter_Halos_Mass_Loss_Rate_Get => null()
  abstract interface
     double precision function Dark_Matter_Halo_Mass_Loss_Rate_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Dark_Matter_Halo_Mass_Loss_Rate_Template
  end interface

contains

  subroutine Dark_Matter_Halo_Mass_Loss_Rates_Initialize
    !% Initialize the dark matter halos mass loss rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterHaloMassLossRateMethod" type="moduleUse">
    include 'dark_matter_halos.mass_loss_rates.modules.inc'
    !# </include>
    implicit none
    
    ! Initialize if necessary.
    if (.not.darkMatterHaloMassLossRateInitialized) then
       !$omp critical(Dark_Matter_Halo_Mass_Loss_Rates_Initialization) 
       if (.not.darkMatterHaloMassLossRateInitialized) then
          ! Get the dark matter halo mass loss rate method parameter.
          !@ <inputParameter>
          !@   <name>darkMatterHaloMassLossRateMethod</name>
          !@   <defaultValue>dynamicalTime</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing mass loss rates from dark matter halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterHaloMassLossRateMethod',darkMatterHaloMassLossRateMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="darkMatterHaloMassLossRateMethod" type="code" action="subroutine">
          !#  <subroutineArgs>darkMatterHaloMassLossRateMethod,Dark_Matter_Halos_Mass_Loss_Rate_Get</subroutineArgs>
          include 'dark_matter_halos.mass_loss_rates.inc'
          !# </include>
          if (.not.associated(Dark_Matter_Halos_Mass_Loss_Rate_Get)) call Galacticus_Error_Report('Dark_Matter_Halo_Mass_Loss_Rates_Initialize'&
               &,'method ' //char(darkMatterHaloMassLossRateMethod)//' is unrecognized')
          darkMatterHaloMassLossRateInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Halo_Mass_Loss_Rates_Initialization) 
    end if
    return
  end subroutine Dark_Matter_Halo_Mass_Loss_Rates_Initialize

  double precision function Dark_Matter_Halos_Mass_Loss_Rate(thisNode)
    !% Returns the rate of mass loss (in $M_\odot$/Gyr) from {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Halo_Mass_Loss_Rates_Initialize

    ! Get the energy using the selected method.
    Dark_Matter_Halos_Mass_Loss_Rate=Dark_Matter_Halos_Mass_Loss_Rate_Get(thisNode)

    return
  end function Dark_Matter_Halos_Mass_Loss_Rate
  
end module Dark_Matter_Halos_Mass_Loss_Rates
