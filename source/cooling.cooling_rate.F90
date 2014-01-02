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

!% Contains a module that implements calculations of the cooling rate.

module Cooling_Rates
  !% Implements calculations of the cooling rate.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Rate, Cooling_Rate_Hot_Halo_Output_Names, Cooling_Rate_Hot_Halo_Output_Count, Cooling_Rate_Hot_Halo_Output

  ! Flag to indicate if this module has been initialized.
  logical                                       :: coolingRateInitialized      =.false.
  logical                                       :: coolingRateOutputInitialized=.false.

  ! Name of cooling rate available method used.
  type     (varying_string           )          :: coolingRateMethod

  ! Option controlling whether cooling rates are output.
  logical                                       :: outputHotHaloCoolingRates

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Rate_Get_Template), pointer :: Cooling_Rate_Get            =>null()
  abstract interface
     double precision function Cooling_Rate_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Rate_Get_Template
  end interface

contains

  subroutine Cooling_Rate_Initialize()
    !% Initialize the cooling rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="coolingRateMethod" type="moduleUse">
    include 'cooling.cooling_rate.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.coolingRateInitialized) then
       !$omp critical(Cooling_Rate_Initialization)
       if (.not.coolingRateInitialized) then
          ! Get the cooling rate method parameter.
          !@ <inputParameter>
          !@   <name>coolingRateMethod</name>
          !@   <defaultValue>White-Frenk1991</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing the cooling rate.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('coolingRateMethod',coolingRateMethod,defaultValue='White-Frenk1991')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="coolingRateMethod" type="functionCall" functionType="void">
          !#  <functionArgs>coolingRateMethod,Cooling_Rate_Get</functionArgs>
          include 'cooling.cooling_rate.inc'
          !# </include>
          if (.not.associated(Cooling_Rate_Get)) call Galacticus_Error_Report('Cooling_Rate','method ' &
               &//char(coolingRateMethod)//' is unrecognized')
          coolingRateInitialized=.true.
       end if
       !$omp end critical(Cooling_Rate_Initialization)
    end if
    return
  end subroutine Cooling_Rate_Initialize

  subroutine Cooling_Rate_Output_Initialize()
    !% Initialize output in the cooling rate module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingRateOutputInitialized) then
       !$omp critical(Cooling_Rate_Output_Initialization)
       if (.not.coolingRateOutputInitialized) then
          ! Get options controlling output.
          !@ <inputParameter>
          !@   <name>outputHotHaloCoolingRates</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Determines whether or not cooling rates and radii are output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHotHaloCoolingRates',outputHotHaloCoolingRates,defaultValue=.false.)
          coolingRateOutputInitialized=.true.
       end if
       !$omp end critical(Cooling_Rate_Output_Initialization)
    end if
    return
  end subroutine Cooling_Rate_Output_Initialize

  double precision function Cooling_Rate(thisNode)
    !% Return the cooling rate for {\tt thisNode} (in units of $M_\odot$ Gyr$^{-1}$).
    !# <include directive="coolingRateModifierMethod" type="moduleUse">
    include 'cooling.cooling_rate.modifier.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Rate_Initialize()

    ! Get the cooling rate using the selected method.
    Cooling_Rate=Cooling_Rate_Get(thisNode)

    ! Call routines that modify the cooling rate.
    !# <include directive="coolingRateModifierMethod" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,Cooling_Rate</functionArgs>
    include 'cooling.cooling_rate.modifier.tasks.inc'
    !# </include>

    return
  end function Cooling_Rate

  !# <mergerTreeOutputNames>
  !#  <unitName>Cooling_Rate_Hot_Halo_Output_Names</unitName>
  !#  <sortName>Cooling_Rate_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Cooling_Rate_Hot_Halo_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode            )              , intent(inout), pointer :: thisNode
    double precision                                    , intent(in   )          :: time
    integer                                             , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*               ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                                          integerPropertyComments, integerPropertyNames
    double precision                      , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Cooling_Rate_Output_Initialize()

    if (outputHotHaloCoolingRates) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>hotHaloCoolingRate</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Rate of mass cooling in the hot halo.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>hotHalo</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='hotHaloCoolingRate'
       doublePropertyComments(doubleProperty)='Rate of mass cooling in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
    end if
    return
  end subroutine Cooling_Rate_Hot_Halo_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Cooling_Rate_Hot_Halo_Output_Count</unitName>
  !#  <sortName>Cooling_Rate_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Cooling_Rate_Hot_Halo_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of hot halo cooling properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    double precision                      , intent(in   )          :: time
    integer                               , intent(inout)          :: doublePropertyCount  , integerPropertyCount
    integer                               , parameter              :: propertyCount      =1

    ! Initialize the module.
    call Cooling_Rate_Output_Initialize()

    if (outputHotHaloCoolingRates) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Cooling_Rate_Hot_Halo_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Cooling_Rate_Hot_Halo_Output</unitName>
  !#  <sortName>Cooling_Rate_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Cooling_Rate_Hot_Halo_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)

    ! Initialize the module.
    call Cooling_Rate_Output_Initialize()

    if (outputHotHaloCoolingRates) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Cooling_Rate(thisNode)
    end if
    return
  end subroutine Cooling_Rate_Hot_Halo_Output

end module Cooling_Rates
