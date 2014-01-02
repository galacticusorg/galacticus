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

!% Contains a module which handles outputting a flag for the most massive progenitor in a tree.

module Galacticus_Output_Most_Massive_Progenitors
  !% Handles outputting a flag for the most massive progenitor in a tree.
  implicit none
  private
  public :: Galacticus_Output_Most_Massive_Progenitor, Galacticus_Output_Most_Massive_Progenitor_Property_Count, Galacticus_Output_Most_Massive_Progenitor_Names

  ! Number of properties.
  integer, parameter :: propertyCount              =1

  ! Record of module initialization.
  logical            :: moduleIsInitialized        =.false.

  ! Output option.
  logical            :: outputMostMassiveProgenitor

contains

  subroutine Galacticus_Output_Most_Massive_Progenitor_Initialize()
    !% Initialize the module that outputs flags for the most massive progenitor in a tree.
    use Input_Parameters
    implicit none

    if (.not.moduleIsInitialized) then
       !$omp critical(Galacticus_Output_Most_Massive_Progenitor_Initialize)
       if (.not.moduleIsInitialized) then
          !@ <inputParameter>
          !@   <name>outputMostMassiveProgenitor</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not most massive progenitor status should be output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputMostMassiveProgenitor',outputMostMassiveProgenitor,defaultValue=.false.)
          ! Record that the module is initialized.
          moduleIsInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Most_Massive_Progenitor_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Most_Massive_Progenitor_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Most_Massive_Progenitor_Names</unitName>
  !#  <sortName>Galacticus_Output_Most_Massive_Progenitor</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Most_Massive_Progenitor_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames    , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Ensure the module is initialized.
    call Galacticus_Output_Most_Massive_Progenitor_Initialize()

    if (outputMostMassiveProgenitor) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>isMostMassiveProgenitor</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Flag indicating if this node is the most massive progenitor in its tree at this time.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='isMostMassiveProgenitor'
       integerPropertyComments(integerProperty)='Flag indicating if this node is the most massive progenitor in its tree at this time.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Most_Massive_Progenitor_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Most_Massive_Progenitor_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Most_Massive_Progenitor</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Most_Massive_Progenitor_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Ensure the module is initialized.
    call Galacticus_Output_Most_Massive_Progenitor_Initialize()

    if (outputMostMassiveProgenitor) integerPropertyCount=integerPropertyCount+propertyCount
    return
  end subroutine Galacticus_Output_Most_Massive_Progenitor_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Most_Massive_Progenitor</unitName>
  !#  <sortName>Galacticus_Output_Most_Massive_Progenitor</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Most_Massive_Progenitor(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store link properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    double precision                          , intent(in   )          :: time
    type            (treeNode          )      , intent(inout), pointer :: thisNode
    integer                                   , intent(inout)          :: doubleBufferCount                    , doubleProperty      , &
         &                                                                integerBufferCount                   , integerProperty
    integer         (kind=kind_int8    )      , intent(inout)          :: integerBuffer            (:,:)
    double precision                          , intent(inout)          :: doubleBuffer             (:,:)
    double precision                    , save                         :: timePrevious                  =-1.0d0
    integer         (kind=kind_int8    ), save                         :: uniqueIdMatched                      , uniqueIdPrevious=-1
    !$omp threadprivate(timePrevious,uniqueIdPrevious,uniqueIdMatched)
    type            (treeNode          )                     , pointer :: currentNode
    class           (nodeComponentBasic)                     , pointer :: currentBasicComponent
    double precision                                                   :: mostMassiveProgenitorMass

    ! Ensure the module is initialized.
    call Galacticus_Output_Most_Massive_Progenitor_Initialize()

    if (outputMostMassiveProgenitor) then
       ! Find the root node in the tree.
       currentNode => thisNode
       do while (associated(currentNode%parent))
          currentNode => currentNode%parent
       end do
       ! Check if this is the same tree, at the same time as on the previous call.
       if (time /= timePrevious .or. currentNode%uniqueId() /= uniqueIdPrevious) then
          ! It is not, so record the new tree root unique ID and the new time.
          timePrevious    =time
          uniqueIdPrevious=currentNode%uniqueId()
          ! Find the most massive progenitor in the tree at this time.
          mostMassiveProgenitorMass=0.0d0
          do while (associated(currentNode))
             currentBasicComponent => currentNode%basic()
             if (currentBasicComponent%time() == time .and. currentBasicComponent%mass() > mostMassiveProgenitorMass) then
                uniqueIdMatched          =currentNode          %uniqueId()
                mostMassiveProgenitorMass=currentBasicComponent%mass    ()
             end if
             call currentNode%walkTree(currentNode)
          end do
       end if
       integerProperty=integerProperty+1
       if (thisNode%uniqueId() == uniqueIdMatched) then
          integerBuffer(integerBufferCount,integerProperty)=1
       else
          integerBuffer(integerBufferCount,integerProperty)=0
       end if
    end if
    return
  end subroutine Galacticus_Output_Most_Massive_Progenitor

end module Galacticus_Output_Most_Massive_Progenitors
