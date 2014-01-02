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

!% Contains a module that implements calculations of mass loss rates from disks due to tidal
!% stripping.

module Tidal_Stripping_Mass_Loss_Rate_Disks
  !% Implements calculations of mass loss rates from disks due to tidal stripping.
  implicit none
  private
  public :: Tidal_Stripping_Mass_Loss_Rate_Disk

  ! Flag to indicate if this module has been initialized.
  logical                                                 :: moduleInitialized                      =.false.

  ! Pointer to the function that actually does the calculation.
  procedure(Tidal_Stripping_Mass_Loss_Rate_Disk), pointer :: Tidal_Stripping_Mass_Loss_Rate_Disk_Get=>null()

contains

  double precision function Tidal_Stripping_Mass_Loss_Rate_Disk(thisNode)
    !% Return the tidal force for the hot halo of {\tt thisNode}.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    use ISO_Varying_String
    !# <include directive="tidalStrippingMassLossRateDisksMethod" type="moduleUse">
    include 'tidal_stripping.mass_loss_rate.disks.modules.inc'
    !# </include>
    implicit none
    type(treeNode      ), intent(inout), pointer :: thisNode
    type(varying_string)                         :: tidalStrippingMassLossRateDisksMethod

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Tidal_Stripping_Mass_Loss_Rate_Disk_Init)
       if (.not.moduleInitialized) then
          ! Get the tidal stripping mass loss rate method parameter.
          !@ <inputParameter>
          !@   <name>tidalStrippingMassLossRateDisksMethod</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing mass loss rates from disks due to tidal stripping.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('tidalStrippingMassLossRateDisksMethod',tidalStrippingMassLossRateDisksMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="tidalStrippingMassLossRateDisksMethod" type="functionCall" functionType="void">
          !#  <functionArgs>tidalStrippingMassLossRateDisksMethod,Tidal_Stripping_Mass_Loss_Rate_Disk_Get</functionArgs>
          include 'tidal_stripping.mass_loss_rate.disks.inc'
          !# </include>
          if (.not.associated(Tidal_Stripping_Mass_Loss_Rate_Disk_Get)) call Galacticus_Error_Report('Tidal_Stripping_Mass_Loss_Rate_Disk','method ' &
               &//char(tidalStrippingMassLossRateDisksMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Tidal_Stripping_Mass_Loss_Rate_Disk_Init)
    end if

    ! Get the mass loss rate using the selected method.
    Tidal_Stripping_Mass_Loss_Rate_Disk=Tidal_Stripping_Mass_Loss_Rate_Disk_Get(thisNode)

    return
  end function Tidal_Stripping_Mass_Loss_Rate_Disk

end module Tidal_Stripping_Mass_Loss_Rate_Disks
