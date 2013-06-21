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

!% Contains a module which handles outputting of tree main branch data to the \glc\ output file.

module Galacticus_Output_Trees_Main_Branch
  !% Handles outputting of tree main branch data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Main_Branch, Galacticus_Output_Tree_Main_Branch_Property_Count,&
       & Galacticus_Output_Tree_Main_Branch_Names

  ! Number of main branch properties.
  integer, parameter :: mainBranchPropertyCount      =1

  ! Flag indicating whether main branch output is required.
  logical            :: outputMainBranchStatus

  ! Flag indicating if module is initialized.
  logical            :: mainBranchOutputIsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Main_Branch_Initalize
    !% Intialize the ``main branch status'' output module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.mainBranchOutputIsInitialized) then
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
    end if

    return
  end subroutine Galacticus_Output_Tree_Main_Branch_Initalize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Main_Branch_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Main_Branch</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Main_Branch_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of main branch properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

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
  subroutine Galacticus_Output_Tree_Main_Branch_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of main branch properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

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
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)

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
