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


!% Contains a module which handles outputting of tree main branch data to the \glc\ output file.

module Galacticus_Output_Trees_Main_Branch
  !% Handles outputting of tree main branch data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Main_Branch, Galacticus_Output_Tree_Main_Branch_Property_Count,&
       & Galacticus_Output_Tree_Main_Branch_Names

  ! Number of main branch properties.
  integer, parameter   :: mainBranchPropertyCount=1

  ! Flag indicating whether main branch output is required.
  logical              :: outputMainBranchStatus

  ! Flag indicating if module is initialized.
  logical              :: mainBranchOutputIsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Main_Branch_Initalize
    !% Intialize the ``main branch status'' output module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    !$omp critical(Galacticus_Output_Tree_Main_Branch_Initalization)
    if (.not.mainBranchOutputIsInitialized) then
       ! Read parameter controlling whether or not this module should output.       
       !@ <inputParameter>
       !@   <name>outputMainBranchStatus</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Controls whether or not the main branch status of each node will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('outputMainBranchStatus',outputMainBranchStatus,defaultValue=.false.)
       ! Flag that the module is now initialized.
       mainBranchOutputIsInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Tree_Main_Branch_Initalization)

    return
  end subroutine Galacticus_Output_Tree_Main_Branch_Initalize
    
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Main_Branch_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Main_Branch</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Main_Branch_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of main branch properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    ! Ensure the module is initialized.
    call Galacticus_Output_Tree_Main_Branch_Initalize

    if (outputMainBranchStatus) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>nodeIsOnMainBranch</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Indicates if the node is on the main branch of the merger tree (0|1).</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='nodeIsOnMainBranch'
       integerPropertyComments(integerProperty)='Indicates if the node is on the main branch of the merger tree (0|1).'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Main_Branch_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Main_Branch_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Main_Branch</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Main_Branch_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of main branch properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    ! Ensure the module is initialized.
    call Galacticus_Output_Tree_Main_Branch_Initalize

    if (outputMainBranchStatus) integerPropertyCount=integerPropertyCount+mainBranchPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Main_Branch_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Main_Branch</unitName>
  !#  <sortName>Galacticus_Output_Tree_Main_Branch</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Main_Branch(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store mainBranch properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

    ! Ensure the module is initialized.
    call Galacticus_Output_Tree_Main_Branch_Initalize

    if (outputMainBranchStatus) then
       integerProperty=integerProperty+1
       select case (thisNode%isOnMainBranch())
       case (.true.)
          integerBuffer(integerBufferCount,integerProperty)=1
       case (.false.)
          integerBuffer(integerBufferCount,integerProperty)=0
       end select
    end if
    return
  end subroutine Galacticus_Output_Tree_Main_Branch

end module Galacticus_Output_Trees_Main_Branch
