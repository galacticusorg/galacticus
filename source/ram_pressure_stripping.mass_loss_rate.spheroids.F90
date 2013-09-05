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

!% Contains a module that implements calculations of mass loss rates from spheroids due to ram pressure
!% stripping.

module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids
  !% Implements calculations of mass loss rates from spheroids due to ram pressure stripping.
  implicit none
  private
  public :: Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid

  ! Flag to indicate if this module has been initialized.  
  logical                                                            :: moduleInitialized=.false.

  ! Pointer to the function that actually does the calculation.
  procedure(Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid), pointer :: Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get => null()

contains

  double precision function Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid(thisNode)
    !% Return the ram pressure force for the hot halo of {\tt thisNode}.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    use ISO_Varying_String
    !# <include directive="ramPressureStrippingMassLossRateSpheroidsMethod" type="moduleUse">
    include 'ram_pressure_stripping.mass_loss_rate.spheroids.modules.inc'
    !# </include>
    implicit none
    type(treeNode      ), intent(inout), pointer :: thisNode
    type(varying_string)                         :: ramPressureStrippingMassLossRateSpheroidsMethod

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Init) 
       if (.not.moduleInitialized) then
          ! Get the ram pressure stripping mass loss rate method parameter.
          !@ <inputParameter>
          !@   <name>ramPressureStrippingMassLossRateSpheroidsMethod</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing mass loss rates from spheroids due to ram pressure stripping.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('ramPressureStrippingMassLossRateSpheroidsMethod',ramPressureStrippingMassLossRateSpheroidsMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="ramPressureStrippingMassLossRateSpheroidsMethod" type="functionCall" functionType="void">
          !#  <functionArgs>ramPressureStrippingMassLossRateSpheroidsMethod,Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get</functionArgs>
          include 'ram_pressure_stripping.mass_loss_rate.spheroids.inc'
          !# </include>
          if (.not.associated(Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get)) call Galacticus_Error_Report('Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid','method ' &
               &//char(ramPressureStrippingMassLossRateSpheroidsMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Init) 
    end if

    ! Get the mass loss rate using the selected method.
    Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid=Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get(thisNode)

    return
  end function Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid

end module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids
