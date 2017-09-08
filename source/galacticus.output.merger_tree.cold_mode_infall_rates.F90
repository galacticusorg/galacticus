!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module that implements output of cold mode infall.

module Cold_Mode_Infall_Rates_Output
  !% Implements output of cold mode infall rates.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Cold_Mode_Infall_Output_Names, Cooling_Cold_Mode_Infall_Output_Count, Cooling_Cold_Mode_Infall_Output

  ! Flag to indicate if this module has been initialized.
  logical :: coldModeInfallRateOutputInitialized=.false.

  ! Option controlling whether cold mode infall rates are output.
  logical :: outputColdModeInfallRate

contains

  subroutine Cooling_Cold_Mode_Infall_Output_Initialize()
    !% Initialize output in the cold mode infall rate module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coldModeInfallRateOutputInitialized) then
       !$omp critical(Cooling_Cold_Mode_Infall_Output_Initialization)
       if (.not.coldModeInfallRateOutputInitialized) then
          ! Get options controlling output.
          !@ <inputParameter>
          !@   <name>outputColdModeInfallRate</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Determines whether or not cold mode infall rates are output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('outputColdModeInfallRate',outputColdModeInfallRate,defaultValue=.false.)
          coldModeInfallRateOutputInitialized=.true.
       end if
       !$omp end critical(Cooling_Cold_Mode_Infall_Output_Initialization)
    end if
    return
  end subroutine Cooling_Cold_Mode_Infall_Output_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Cooling_Cold_Mode_Infall_Output_Names</unitName>
  !#  <sortName>Cooling_Cold_Mode_Infall_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Cooling_Cold_Mode_Infall_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of cold mode infall rate properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode            )              , intent(inout) :: thisNode
    double precision                                    , intent(in   ) :: time
    integer                                             , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*               ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                                 integerPropertyComments, integerPropertyNames
    double precision                      , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: thisNode, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, time
    
    ! Initialize the module.
    call Cooling_Cold_Mode_Infall_Output_Initialize()

    if (outputColdModeInfallRate) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>coldModeInfallRate</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Rate of infall of cold mode material onto the galaxy.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>hotHalo</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='coldModeInfallRate'
       doublePropertyComments(doubleProperty)='Rate of infall of cold mode material onto the galaxy.'
       doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
    end if
    return
  end subroutine Cooling_Cold_Mode_Infall_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Cooling_Cold_Mode_Infall_Output_Count</unitName>
  !#  <sortName>Cooling_Cold_Mode_Infall_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Cooling_Cold_Mode_Infall_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of cold mode infall rate properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode            ), intent(inout) :: thisNode
    double precision                      , intent(in   ) :: time
    integer                               , intent(inout) :: doublePropertyCount  , integerPropertyCount
    integer                               , parameter     :: propertyCount      =1
    !GCC$ attributes unused :: thisNode, integerPropertyCount, time
    
    ! Initialize the module.
    call Cooling_Cold_Mode_Infall_Output_Initialize()

    if (outputColdModeInfallRate) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Cooling_Cold_Mode_Infall_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Cooling_Cold_Mode_Infall_Output</unitName>
  !#  <sortName>Cooling_Cold_Mode_Infall_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Cooling_Cold_Mode_Infall_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store cold mode infall rate properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Multi_Counters
    use Cooling_Cold_Mode_Infall_Rates
    implicit none
    double precision                         , intent(in   ) :: time
    type            (treeNode               ), intent(inout) :: thisNode
    integer                                  , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8         ), intent(inout) :: integerBuffer    (:,:)
    double precision                         , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter           ), intent(inout) :: instance
    class           (coldModeInfallRateClass), pointer       :: coldModeInfallRate_
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Cooling_Cold_Mode_Infall_Output_Initialize()

    if (outputColdModeInfallRate) then
       coldModeInfallRate_ => coldModeInfallRate()
       doubleProperty      =  doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=coldModeInfallRate_%infallRate(thisNode)
    end if
    return
  end subroutine Cooling_Cold_Mode_Infall_Output

end module Cold_Mode_Infall_Rates_Output
