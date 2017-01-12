! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which filters output for lightcone geometry.

module Galacticus_Merger_Tree_Output_Filter_Lightcones
  !% Filters output for lightcone geometry.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Galacticus_Output_Tree_Lightcone      , Galacticus_Output_Tree_Lightcone_Property_Count, &
       &    Galacticus_Output_Tree_Lightcone_Names, Galacticus_Output_Tree_Lightcone_Instances

  ! Number of lightcone properties.
  integer          , parameter :: lightconePropertyCount    =8

  ! Flags indicating if the module has been initialized and if this output is active.
  logical                      :: lightconeOutputInitialized=.false.
  logical                      :: outputLightconeData

  ! Index of lightcone replicant counter in the output instance counter.
  integer(c_size_t)            :: instanceIndex
  
contains

  subroutine Galacticus_Output_Lightcone_Initialize()
    !% Initializes the module by determining whether or not lightcone data should be output.
    use Input_Parameters
    implicit none

    if (.not.lightconeOutputInitialized) then
       !$omp critical(Galacticus_Output_Lightcone_Initialize)
       if (.not.lightconeOutputInitialized) then
          !@ <inputParameter>
          !@   <name>outputLightconeData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not satellite host data (node mass) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputLightconeData',outputLightconeData,defaultValue=.false.)
          ! Flag that module is now initialized.
          lightconeOutputInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Lightcone_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Lightcone_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Lightcone_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: node
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: node, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, time

    call Galacticus_Output_Lightcone_Initialize()
    if (outputLightconeData) then
       !@ <outputPropertyGroup>
       !@   <name>lightconePosition</name>
       !@   <description>Lightcone X-Y-Z coordinates</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       !@ <outputPropertyGroup>
       !@   <name>lightconeVelocity</name>
       !@   <description>Lightcone X-Y-Z velocity components</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionX</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (radial axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionX'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (radial axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionY</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (1st angular axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionY'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (1st angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionZ</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (2nd angular axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionZ'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (2nd angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityX</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (radial axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityX'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (radial axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityY</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (1st angular axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityY'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (1st angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityZ</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (2nd angular axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityZ'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (2nd angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeRedshift</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Reshift of galaxy in lightcone.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeRedshift'
       doublePropertyComments(doubleProperty)='Reshift of galaxy in lightcone.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>angularWeight</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Number of such galaxies per unit area [degrees^-2].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='angularWeight'
       doublePropertyComments(doubleProperty)='Number of such galaxies per unit area [degrees⁻²].'
       doublePropertyUnitsSI (doubleProperty)=1.0d0/degreesToRadians**2
    end if
    return
  end subroutine Galacticus_Output_Tree_Lightcone_Names

  !# <mergerTreeOutputInstances>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Instances</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputInstances>
  subroutine Galacticus_Output_Tree_Lightcone_Instances(node,output,instance)
    !% Set the number of instances for the given {\normalfont \ttfamily node} to be output in the lightcone.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Nodes
    use               Multi_Counters
    use               Geometry_Lightcones
    use               Galacticus_Error
    use               String_Handling
    use               ISO_Varying_String
    implicit none
    type     (treeNode              ), intent(inout), pointer :: node
    integer  (c_size_t              ), intent(in   )          :: output
    type     (multiCounter          ), intent(inout)          :: instance
    class    (geometryLightconeClass)               , pointer :: geometryLightcone_
    integer  (c_size_t              )                         :: replicationCount
    type     (varying_string        )                         :: message
    character(len=5                 )                         :: label
    
    call Galacticus_Output_Lightcone_Initialize()
    if (outputLightconeData) then
       geometryLightcone_  => geometryLightcone ()
       replicationCount=geometryLightcone_%replicationCount(node,output)
       if (replicationCount < 1_c_size_t) then
          if (geometryLightcone_%isInLightcone(node)) then
             label="true"
          else
             label="false"
          end if
          message=var_str("Node ")//node%index()//" appears in "//replicationCount//"(<1) replicants - this should not happen - lightcone intersection reports '"//trim(label)//"'"
          call Galacticus_Error_Report('Galacticus_Output_Tree_Lightcone_Instances',message)
       end if
       call instance%append(replicationCount)
       instanceIndex=instance%count()
    end if
    return
  end subroutine Galacticus_Output_Tree_Lightcone_Instances
  
  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Lightcone_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: node, time, integerPropertyCount
    
    call Galacticus_Output_Lightcone_Initialize()
    if (outputLightconeData) doublePropertyCount=doublePropertyCount+lightconePropertyCount
    return
  end subroutine Galacticus_Output_Tree_Lightcone_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Lightcone</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Lightcone(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store link properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Geometry_Lightcones
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Vectors
    use Multi_Counters
    implicit none
    double precision                         , intent(in   )                 :: time
    type            (treeNode               ), intent(inout)                 :: node
    integer                                  , intent(inout)                 :: doubleBufferCount , doubleProperty , &
         &                                                                      integerBufferCount, integerProperty
    integer         (kind=kind_int8         ), intent(inout), dimension(:,:) :: integerBuffer
    double precision                         , intent(inout), dimension(:,:) :: doubleBuffer
    type            (multiCounter           ), intent(inout)                 :: instance
    class           (geometryLightconeClass )               , pointer        :: geometryLightcone_
    class           (cosmologyFunctionsClass)               , pointer        :: cosmologyFunctions_
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer
    
    call Galacticus_Output_Lightcone_Initialize()
    if (outputLightconeData) then
       geometryLightcone_  => geometryLightcone ()
       cosmologyFunctions_ => cosmologyFunctions()
       doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+3)=geometryLightcone_%position(node,instance%state(instanceIndex))
       doubleProperty=doubleProperty+3
       doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+3)=geometryLightcone_%velocity(node,instance%state(instanceIndex))
       doubleProperty=doubleProperty+3
       doubleBuffer(doubleBufferCount,doubleProperty+1                 )=cosmologyFunctions_   %redshiftFromExpansionFactor(     &
            &                                                             cosmologyFunctions_  %expansionFactor             (    &
            &                                                              cosmologyFunctions_ %timeAtDistanceComoving       (   &
            &                                                               Vector_Magnitude                                  (  &
            &                                                                geometryLightcone_%position                       ( &
            &                                                                 node                         ,                     &
            &                                                                 instance%state(instanceIndex)                      &
            &                                                                                                                  ) &
            &                                                                                                                 )  &
            &                                                                                                                )   &
            &                                                                                                               )    &
            &                                                                                                              )
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty+1                 )=geometryLightcone_%solidAngle()/degreesToRadians**2
       doubleProperty=doubleProperty+1
    end if
    return
  end subroutine Galacticus_Output_Tree_Lightcone

end module Galacticus_Merger_Tree_Output_Filter_Lightcones
