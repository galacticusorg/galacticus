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

!% Contains a module which handles outputting of node dark matter profile rotation velocity maximum to the \glc\ output file.

module Galacticus_Output_Trees_Velocity_Maximum
  !% Handles outputting of node dark matter profile rotation velocity maxima to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Velocity_Maximum, Galacticus_Output_Tree_Velocity_Maximum_Property_Count, Galacticus_Output_Tree_Velocity_Maximum_Names

  ! Number of velocity maximum properties.
  integer, parameter :: velocityMaximumPropertyCount        =1

  ! Flag indicating whether or not velocityMaximum information is to be output.
  logical            :: outputVelocityMaximumData

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputVelocityMaximumDataInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Velocity_Maximum_Initialize
    !% Initializes the module by determining whether or not velocity maximum data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputVelocityMaximumDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Velocity_Maximum_Initialize)
       if (.not.outputVelocityMaximumDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputVelocityMaximumData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not dark matter profile rotation velocity maximum data should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputVelocityMaximumData',outputVelocityMaximumData,defaultValue=.false.)
          ! Flag that module is now initialized.
          outputVelocityMaximumDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Velocity_Maximum_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Maximum_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Maximum_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Maximum</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Velocity_Maximum_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of velocityMaximum properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: node
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: node, time, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Maximum_Initialize()
    ! Return property names if we are outputting velocityMaximum data.
    if (outputVelocityMaximumData) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>darkMatterProfileVelocityMaximum</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Maximum rotation velocity of the dark matter profile [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='darkMatterProfileVelocityMaximum'
       doublePropertyComments(doubleProperty)='Maximum rotation velocity of the dark matter profile [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Maximum_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Maximum_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Maximum</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Velocity_Maximum_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of velocityMaximum properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: node, integerPropertyCount, time
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Maximum_Initialize()
    ! Increment property count if we are outputting velocityMaximum data.
    if (outputVelocityMaximumData) doublePropertyCount=doublePropertyCount+velocityMaximumPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Velocity_Maximum_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Maximum</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Maximum</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Velocity_Maximum(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store velocityMaximum properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Dark_Matter_Profiles
    use Kind_Numbers
    use Multi_Counters
    implicit none
    double precision                        , intent(in   )                 :: time
    type            (treeNode              ), intent(inout)                 :: node
    integer                                 , intent(inout)                 :: doubleBufferCount , doubleProperty , &
         &                                                                     integerBufferCount, integerProperty
    integer         (kind=kind_int8        ), intent(inout), dimension(:,:) :: integerBuffer    
    double precision                        , intent(inout), dimension(:,:) :: doubleBuffer     
    type            (multiCounter          ), intent(inout)                 :: instance
    class           (darkMatterProfileClass)               , pointer        :: darkMatterProfile_
    !GCC$ attributes unused :: time, integerProperty, integerBuffer, integerBufferCount, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Maximum_Initialize()
    ! Store property data if we are outputting velocity maximum data.
    if (outputVelocityMaximumData) then
       darkMatterProfile_  => darkMatterProfile()
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=darkMatterProfile_%circularVelocityMaximum(node)
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Maximum

end module Galacticus_Output_Trees_Velocity_Maximum
