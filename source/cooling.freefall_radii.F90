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

!% Contains a module that implements calculations of the freefall radius in cooling calculations.

module Freefall_Radii
  !% Implements calculations of the freefall radius in cooling calculations.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="freefallRadiusMethod" type="moduleUse">
  include 'cooling.freefall_radius.modules.inc'
  !# </include>
  implicit none
  private
  public :: Freefall_Radius, Freefall_Radius_Growth_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: freefallRadiusInitialized=.false.

  ! Name of cooling radius available method used.
  type(varying_string) :: freefallRadiusMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Freefall_Radius_Get_Template), pointer :: Freefall_Radius_Get => null()
  procedure(Freefall_Radius_Get_Template), pointer :: Freefall_Radius_Growth_Rate_Get => null()
  abstract interface
     double precision function Freefall_Radius_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Freefall_Radius_Get_Template
  end interface
  
contains

  subroutine Freefall_Radius_Initialize
    !% Initialize the cooling radius module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.freefallRadiusInitialized) then
       !$omp critical(Freefall_Radius_Initialization) 
       if (.not.freefallRadiusInitialized) then
          ! Get the cooling radius method parameter.
          !@ <inputParameter>
          !@   <name>freefallRadiusMethod</name>
          !@   <defaultValue>darkMatterHalo</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the freefall radius in cooling calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('freefallRadiusMethod',freefallRadiusMethod,defaultValue='darkMatterHalo')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="freefallRadiusMethod" type="functionCall" functionType="void">
          !#  <functionArgs>freefallRadiusMethod,Freefall_Radius_Get,Freefall_Radius_Growth_Rate_Get</functionArgs>
          include 'cooling.freefall_radius.inc'
          !# </include>
          if (.not.(associated(Freefall_Radius_Get).and.associated(Freefall_Radius_Growth_Rate_Get))) call&
               & Galacticus_Error_Report('Freefall_Radius','method ' //char(freefallRadiusMethod)//' is unrecognized')
          freefallRadiusInitialized=.true.
       end if
       !$omp end critical(Freefall_Radius_Initialization) 
    end if
    return
  end subroutine Freefall_Radius_Initialize

  double precision function Freefall_Radius(thisNode)
    !% Return the freefall radius for cooling calculations for {\tt thisNode} (in units of Mpc).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Freefall_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Freefall_Radius=Freefall_Radius_Get(thisNode)

    return
  end function Freefall_Radius

  double precision function Freefall_Radius_Growth_Rate(thisNode)
    !% Return the rate at which the freefall radius for cooling calculations grows for {\tt thisNode} (in units of Mpc/Gyr).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Freefall_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Freefall_Radius_Growth_Rate=Freefall_Radius_Growth_Rate_Get(thisNode)

    return
  end function Freefall_Radius_Growth_Rate

end module Freefall_Radii
