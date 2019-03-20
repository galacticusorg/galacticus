!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module that implements output of cooling rates.

module Cooling_Rates_Output
  !% Implements output of cooling rates.
  use ISO_Varying_String
  use Galacticus_Nodes, only : treeNode
  implicit none
  private
  public :: Cooling_Rate_Hot_Halo_Output_Names, Cooling_Rate_Hot_Halo_Output_Count, Cooling_Rate_Hot_Halo_Output

  ! Flag to indicate if this module has been initialized.
  logical :: coolingRateOutputInitialized=.false.

  ! Option controlling whether cooling rates are output.
  logical :: outputHotHaloCoolingRates

contains

  subroutine Cooling_Rate_Output_Initialize()
    !% Initialize output in the cooling rate module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingRateOutputInitialized) then
       !$omp critical(Cooling_Rate_Output_Initialization)
       if (.not.coolingRateOutputInitialized) then
          ! Get options controlling output.
          !# <inputParameter>
          !#   <name>outputHotHaloCoolingRates</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Determines whether or not cooling rates and radii are output.</description>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          coolingRateOutputInitialized=.true.
       end if
       !$omp end critical(Cooling_Rate_Output_Initialization)
    end if
    return
  end subroutine Cooling_Rate_Output_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Cooling_Rate_Hot_Halo_Output_Names</unitName>
  !#  <sortName>Cooling_Rate_Hot_Halo_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Cooling_Rate_Hot_Halo_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
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
    call Cooling_Rate_Output_Initialize()

    if (outputHotHaloCoolingRates) then
       doubleProperty=doubleProperty+1
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
    type            (treeNode            ), intent(inout) :: thisNode
    double precision                      , intent(in   ) :: time
    integer                               , intent(inout) :: doublePropertyCount  , integerPropertyCount
    integer                               , parameter     :: propertyCount      =1
    !GCC$ attributes unused :: thisNode, integerPropertyCount, time
    
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
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Multi_Counters
    use Cooling_Rates
    implicit none
    double precision                  , intent(in   ) :: time
    type            (treeNode        ), intent(inout) :: thisNode
    integer                           , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                               integerProperty
    integer         (kind=kind_int8  ), intent(inout) :: integerBuffer    (:,:)
    double precision                  , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter    ), intent(inout) :: instance
    class           (coolingRateClass), pointer       :: coolingRate_
    !GCC$ attributes unused :: integerProperty, integerBufferCount, integerBuffer, time, instance
    
    ! Initialize the module.
    call Cooling_Rate_Output_Initialize()
    if (outputHotHaloCoolingRates) then
       coolingRate_ => coolingRate()
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=coolingRate_%rate(thisNode)
    end if
    return
  end subroutine Cooling_Rate_Hot_Halo_Output

end module Cooling_Rates_Output
