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

!% Contains a module that implements calculations of the cooling radius.

module Cooling_Radii
  !% Implements calculations of the cooling radius.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="coolingRadiusMethod" type="moduleUse">
  include 'cooling.cooling_radius.modules.inc'
  !# </include>
  implicit none
  private
  public :: Cooling_Radius, Cooling_Radius_Growth_Rate, Cooling_Radius_Hot_Halo_Output, Cooling_Radius_Hot_Halo_Output_Names,&
       & Cooling_Radius_Hot_Halo_Output_Count

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingRadiusInitialized      =.false.
  logical              :: coolingRadiusOutputInitialized=.false.

  ! Name of cooling radius available method used.
  type(varying_string) :: coolingRadiusMethod

  ! Option controlling whether cooling radii are output.
  logical              :: outputHotHaloCoolingRadii

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Radius_Get_Template), pointer :: Cooling_Radius_Get => null()
  procedure(Cooling_Radius_Get_Template), pointer :: Cooling_Radius_Growth_Rate_Get => null()
  abstract interface
     double precision function Cooling_Radius_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Radius_Get_Template
  end interface
  
contains

  subroutine Cooling_Radius_Initialize
    !% Initialize the cooling radius module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingRadiusInitialized) then
       !$omp critical(Cooling_Radius_Initialization) 
       if (.not.coolingRadiusInitialized) then
          ! Get the cooling radius method parameter.
          !@ <inputParameter>
          !@   <name>coolingRadiusMethod</name>
          !@   <defaultValue>simple</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the cooling radius.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('coolingRadiusMethod',coolingRadiusMethod,defaultValue='simple')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="coolingRadiusMethod" type="functionCall" functionType="void">
          !#  <functionArgs>coolingRadiusMethod,Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get</functionArgs>
          include 'cooling.cooling_radius.inc'
          !# </include>
          if (.not.(associated(Cooling_Radius_Get).and.associated(Cooling_Radius_Growth_Rate_Get))) call&
               & Galacticus_Error_Report('Cooling_Radius','method ' //char(coolingRadiusMethod)//' is unrecognized')

          coolingRadiusInitialized=.true.
       end if
       !$omp end critical(Cooling_Radius_Initialization) 
    end if
    return
  end subroutine Cooling_Radius_Initialize

  subroutine Cooling_Radius_Output_Initialize()
    !% Initialize output in the cooling radius module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingRadiusOutputInitialized) then
       !$omp critical(Cooling_Radius_Output_Initialization) 
       if (.not.coolingRadiusOutputInitialized) then
          ! Get options controlling output.
          !@ <inputParameter>
          !@   <name>outputHotHaloCoolingRadii</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Determines whether or not cooling radii are output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHotHaloCoolingRadii',outputHotHaloCoolingRadii,defaultValue=.false.)

          coolingRadiusOutputInitialized=.true.
       end if
       !$omp end critical(Cooling_Radius_Output_Initialization) 
    end if
    return
  end subroutine Cooling_Radius_Output_Initialize

  double precision function Cooling_Radius(thisNode)
    !% Return the cooling radius for {\tt thisNode} (in units of Mpc).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Radius=Cooling_Radius_Get(thisNode)

    return
  end function Cooling_Radius

  double precision function Cooling_Radius_Growth_Rate(thisNode)
    !% Return the rate at which the cooling radius grows for {\tt thisNode} (in units of Mpc/Gyr).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Radius_Growth_Rate=Cooling_Radius_Growth_Rate_Get(thisNode)

    return
  end function Cooling_Radius_Growth_Rate

  !# <mergerTreeOutputNames>
  !#  <unitName>Cooling_Radius_Hot_Halo_Output_Names</unitName>
  !#  <sortName>Cooling_Radius_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Cooling_Radius_Hot_Halo_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Abundances_Structure
    use ISO_Varying_String
    implicit none
    type (treeNode            ), intent(inout), pointer      :: thisNode
    double precision           , intent(in   )               :: time
    integer                    , intent(inout)               :: integerProperty,doubleProperty
    character(len=*)           , intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision           , intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
    ! Initialize the module.
    call Cooling_Radius_Output_Initialize()

    if (outputHotHaloCoolingRadii) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>hotHaloCoolingRadius</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Cooling radius in the hot halo.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>hotHalo</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='hotHaloCoolingRadius'
       doublePropertyComments(doubleProperty)='Cooling radius in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
    end if
    return
  end subroutine Cooling_Radius_Hot_Halo_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Cooling_Radius_Hot_Halo_Output_Count</unitName>
  !#  <sortName>Cooling_Radius_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Cooling_Radius_Hot_Halo_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of hot halo cooling properties to be written to the the \glc\ output file.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    double precision           , intent(in   )          :: time
    integer                    , intent(inout)          :: integerPropertyCount,doublePropertyCount
    integer                    , parameter              :: propertyCount=1

    ! Initialize the module.
    call Cooling_Radius_Output_Initialize()

    if (outputHotHaloCoolingRadii) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Cooling_Radius_Hot_Halo_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Cooling_Radius_Hot_Halo_Output</unitName>
  !#  <sortName>Cooling_Radius_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Cooling_Radius_Hot_Halo_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )              :: time
    type            (treeNode      ), intent(inout), pointer     :: thisNode
    integer                         , intent(inout)              :: integerProperty,integerBufferCount,doubleProperty&
         &,doubleBufferCount
    integer         (kind=kind_int8), intent(inout)              :: integerBuffer(:,:)
    double precision                , intent(inout)              :: doubleBuffer (:,:)

    ! Initialize the module.
    call Cooling_Radius_Output_Initialize()

    if (outputHotHaloCoolingRadii) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Cooling_Radius(thisNode)
    end if
    return
  end subroutine Cooling_Radius_Hot_Halo_Output

end module Cooling_Radii
