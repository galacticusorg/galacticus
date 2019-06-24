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

!% Contains a module which handles outputting of node virial data to the \glc\ output file.

module Galacticus_Output_Trees_Halo_Environment
  !% Handles outputting of halo environemtn data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Halo_Environment, Galacticus_Output_Tree_Halo_Environment_Property_Count, Galacticus_Output_Tree_Halo_Environment_Names

  ! Number of environment properties.
  integer, parameter :: haloEnvironmentPropertyCount        =2

  ! Flag indicating whether or not haloEnvironment information is to be output.
  logical            :: outputHaloEnvironmentData

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputHaloEnvironmentDataInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Halo_Environment_Initialize
    !% Initializes the module by determining whether or not halo environment data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputHaloEnvironmentDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Halo_Environment_Initialize)
       if (.not.outputHaloEnvironmentDataInitialized) then
          !# <inputParameter>
          !#   <name>outputHaloEnvironmentData</name>
          !#   <defaultValue>.false.</defaultValue>
          !#   <source>globalParameters</source>
          !#   <description>Specifies whether or not halo environment data (overdensity) should be included in the output.</description>
          !#   <type>boolean</type>
          !#   <cardinality>1</cardinality>
          !#   <group>output</group>
          !# </inputParameter>
          ! Flag that module is now initialized.
          outputHaloEnvironmentDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Halo_Environment_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Halo_Environment_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Halo_Environment_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Halo_Environment</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Halo_Environment_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of haloEnvironment properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: thisNode
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: thisNode, time, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Halo_Environment_Initialize()
    ! Return property names if we are outputting haloEnvironment data.
    if (outputHaloEnvironmentData) then
       !@ <outputProperty>
       !@   <name>haloEnvironmentOverdensityLinear</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Environmental linear overdensity of the halo [].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doubleProperty                        =doubleProperty+1
       doublePropertyNames   (doubleProperty)='haloEnvironmentOverdensityLinear'
       doublePropertyComments(doubleProperty)='Environmental linear overdensity of the halo [].'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
       !@ <outputProperty>
       !@   <name>haloEnvironmentOverdensityNonLinear</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Environmental linear overdensity of the halo [].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doubleProperty                        =doubleProperty+1
       doublePropertyNames   (doubleProperty)='haloEnvironmentOverdensityNonLinear'
       doublePropertyComments(doubleProperty)='Environmental non-linear overdensity of the halo [].'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Halo_Environment_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Halo_Environment_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Halo_Environment</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Halo_Environment_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of haloEnvironment properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: thisNode, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Halo_Environment_Initialize()
    ! Increment property count if we are outputting haloEnvironment data.
    if (outputHaloEnvironmentData) doublePropertyCount=doublePropertyCount+haloEnvironmentPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Halo_Environment_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Halo_Environment</unitName>
  !#  <sortName>Galacticus_Output_Tree_Halo_Environment</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Halo_Environment(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store haloEnvironment properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Cosmological_Density_Field
    use Kind_Numbers
    use Multi_Counters
    implicit none
    double precision                      , intent(in   ) :: time
    type            (treeNode            ), intent(inout) :: thisNode
    integer                               , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                   integerProperty
    integer         (kind=kind_int8      ), intent(inout) :: integerBuffer    (:,:)
    double precision                      , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter        ), intent(inout) :: instance
    class           (haloEnvironmentClass), pointer       :: haloEnvironment_
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Halo_Environment_Initialize()
    ! Store property data if we are outputting haloEnvironment data.
    if (outputHaloEnvironmentData) then
       haloEnvironment_ => haloEnvironment()
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=haloEnvironment_%overdensityLinear   (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=haloEnvironment_%overdensityNonLinear(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Halo_Environment

end module Galacticus_Output_Trees_Halo_Environment
