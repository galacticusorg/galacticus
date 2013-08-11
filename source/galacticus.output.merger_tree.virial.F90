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

!% Contains a module which handles outputting of node virial data to the \glc\ output file.

module Galacticus_Output_Trees_Virial
  !% Handles outputting of node virial data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Virial, Galacticus_Output_Tree_Virial_Property_Count, Galacticus_Output_Tree_Virial_Names

  ! Number of virial properties.
  integer, parameter :: virialPropertyCount        =2

  ! Flag indicating whether or not virial information is to be output.
  logical            :: outputVirialData

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputVirialDataInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Virial_Initialize
    !% Initializes the module by determining whether or not virial data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputVirialDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Virial_Initialize)
       if (.not.outputVirialDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputVirialData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not virial data (radius, velocity) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputVirialData',outputVirialData,defaultValue=.false.)

          ! Flag that module is now initialized.
          outputVirialDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Virial_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Virial_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Virial_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Virial</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Virial_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of virial properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Virial_Initialize

    ! Return property names if we are outputting virial data.
    if (outputVirialData) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeVirialRadius</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Virial radius of the node [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeVirialRadius'
       doublePropertyComments(doubleProperty)='Virial radius of the node [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeVirialVelocity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Virial velocity of the node [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeVirialVelocity'
       doublePropertyComments(doubleProperty)='Virial velocity of the node [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
    end if
    return
  end subroutine Galacticus_Output_Tree_Virial_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Virial_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Virial</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Virial_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of virial properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Virial_Initialize

    ! Increment property count if we are outputting virial data.
    if (outputVirialData) doublePropertyCount=doublePropertyCount+virialPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Virial_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Virial</unitName>
  !#  <sortName>Galacticus_Output_Tree_Virial</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Virial(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store virial properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)

    ! Initialize the module.
    call Galacticus_Output_Tree_Virial_Initialize

    ! Store property data if we are outputting virial data.
    if (outputVirialData) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Dark_Matter_Halo_Virial_Radius  (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Dark_Matter_Halo_Virial_Velocity(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Virial

end module Galacticus_Output_Trees_Virial
