!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which handles outputting of tree descendent data to the \glc\ output file.

module Galacticus_Output_Trees_Descendents
  !% Handles outputting of tree descendent data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Descendents, Galacticus_Output_Tree_Descendents_Property_Count, Galacticus_Output_Tree_Descendents_Names

  ! Number of descendent properties.
  integer, parameter :: descendentPropertyCount=1

  ! Flag indicating whether or not descendent index information is to be output.
  logical            :: outputDescendentIndices

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputDescendentsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Descendents_Initialize
    !% Initializes the module by determining whether or not descendent data should be output.
    use Input_Parameters
    implicit none

    !$omp critical(Galacticus_Output_Tree_Descendents_Initialize)
    if (.not.outputDescendentsInitialized) then
       !@ <inputParameter>
       !@   <name>outputDescendentIndices</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not descendent indices (i.e. index of the node at the next output) should be included in the output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('outputDescendentIndices',outputDescendentIndices,defaultValue=.false.)

       ! Flag that module is now initialized.
       outputDescendentsInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Tree_Descendents_Initialize)
    return
  end subroutine Galacticus_Output_Tree_Descendents_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Descendents_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Descendents_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of descendent properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>descendentIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>ID of the node which this node will have descended into by the next timestep.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='descendentIndex'
       integerPropertyComments(integerProperty)='ID of the node which this node will have descended into by the next timestep.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Descendents_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Descendents_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Descendents_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of descendent properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) integerPropertyCount=integerPropertyCount+descendentPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Descendents_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Descendents</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Descendents(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store descendent properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    type(treeNode),                         pointer :: descendentNode
    double precision                                :: outputTimeNext
    logical                                         :: foundDescendent

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) then
       ! Get the time of the next output.
       outputTimeNext=Galacticus_Next_Output_Time(time)
       
       foundDescendent=.false.
       integerProperty=integerProperty+1
       if (outputTimeNext < 0.0d0) then
          ! There is no next output time.
          integerBuffer(integerBufferCount,integerProperty)=thisNode%index()
          foundDescendent=.true.
       else if (thisNode%isSatellite()) then
          ! Node is a satellite, so it's node index will remain unchanged.
          if (Tree_Node_Satellite_Time_Of_Merging(thisNode) > outputTimeNext) then
             ! Satellite will not have merged prior to the next output time, so retains its own index.
             integerBuffer(integerBufferCount,integerProperty)=thisNode%index()
             foundDescendent=.true.
          else
             ! Satellite will merge prior to the next output time - find the node it merges with.
             call thisNode%mergesWith(descendentNode)
          end if
       else
          ! Node is not a satellite, so set the initial descendent to itself.
          descendentNode => thisNode
       end if
       ! Check if we still need to find the descendent.
       if (.not.foundDescendent) then
          ! No descendent has yet been found, so trace forward in time until we find one.
          ! Continue until the tree base is reached, or the next output time is reached.
          do while (.not.foundDescendent)
             ! If the next output time has been surpassed, then we are finished.
             if (Tree_Node_Time(descendentNode) >= outputTimeNext) then
                foundDescendent=.true.
             else
                ! Test whether this node is the primary progenitor.
                if (descendentNode%isPrimaryProgenitor()) then
                   ! It is, so simply move to the parent node.
                   descendentNode => descendentNode%parentNode
                else
                   ! It is not, so it becomes a satellite. Test whether it has a merge target associated with it.
                   if (associated(descendentNode%mergeNode)) then
                      ! It does. If merging occurs before the next output time, jump to that node. Otherwise we are finished.
                      if (Tree_Node_Satellite_Time_Of_Merging(descendentNode) <= outputTimeNext) then
                         descendentNode => descendentNode%mergeNode
                      else
                         foundDescendent=.true.
                      end if
                   else
                      ! We no longer can tell if this node will exist as a separate entity at the next output time. Assume that it
                      ! will, and therefore we are finished.
                      foundDescendent=.true.
                   end if
                end if
             end if
          end do
          ! If the descendent exists after the next output time, then we've gone one step too far - back up a step.
          if (Tree_Node_Time(descendentNode) > outputTimeNext .and. associated(descendentNode%childNode)) descendentNode => descendentNode%childNode
          ! Store the descendent index.
          integerBuffer(integerBufferCount,integerProperty)=descendentNode%index()
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Descendents

end module Galacticus_Output_Trees_Descendents
