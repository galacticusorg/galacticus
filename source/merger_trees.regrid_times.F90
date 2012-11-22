!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which prunes branches below a given mass threshold from merger trees.

module Merger_Trees_Regrid_Times
  !% Forces a merger tree onto a specified time grid.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Regrid_Time

  ! Flag indicating if module is initialized.
  logical                                         :: regridTimeModuleInitialized=.false.

  ! Flag indicating if regridding is required.
  logical                                         :: mergerTreeRegridTimes

  ! Flag indicating if dumping of merger trees is required.
  logical                                         :: mergerTreeRegridDumpTrees

  ! Variables expressing the distribution of grid points.
  integer                                         :: mergerTreeRegridCount,mergerTreeRegridSpacing
  type(varying_string)                            :: mergerTreeRegridSpacingAsText
  double precision                                :: mergerTreeRegridStartExpansionFactor,mergerTreeRegridEndExpansionFactor
  double precision,     allocatable, dimension(:) :: mergerTreeRegridTimeGrid
  integer,              parameter                 :: mergerTreeRegridSpacingLinear                =0
  integer,              parameter                 :: mergerTreeRegridSpacingLogarithmic           =1
  integer,              parameter                 :: mergerTreeRegridSpacingLogCriticalOverdensity=2

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Regrid_Time</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Regrid_Time(thisTree)
    !% Regrid times of halos in {\tt thisTree}.
    use Merger_Trees
    use Galacticus_Nodes
    use Input_Parameters
    use Numerical_Ranges
    use Memory_Management
    use Cosmology_Functions
    use Critical_Overdensity
    use Galacticus_Error
    use FGSL
    use Numerical_Interpolation
    use Kind_Numbers
    use Merger_Trees_Dump
    implicit none
    type(mergerTree),          intent(inout), target       :: thisTree
    type(treeNode),            pointer                     :: thisNode,childNode,siblingNode,nextNode
    type(treeNode),            pointer,       dimension(:) :: newNodes
    integer(kind=kind_int8),   allocatable,   dimension(:) :: highlightNodes
    class(nodeComponentBasic), pointer                     :: thisBasicComponent,parentBasicComponent,childBasicComponent
    type (mergerTree        ), pointer                     :: currentTree
    type(fgsl_interp_accel)                                :: interpolationAccelerator
    logical                                                :: interpolationReset
    integer                                                :: iNow,iParent,iTime,allocErr
    double precision                                       :: timeNow,timeParent,massNow,massParent
    integer(kind=kind_int8)                                :: nodeIndex,firstNewNode

    ! Check if module is initialized.
    if (.not.regridTimeModuleInitialized) then
       !$omp critical (Merger_Tree_Regrid_Time_Initialize)
       if (.not.regridTimeModuleInitialized) then
          ! Get parameter specifying if regridding is required.
          !@ <inputParameter>
          !@   <name>mergerTreeRegridTimes</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to regrid merger tree times.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeRegridTimes',mergerTreeRegridTimes,defaultValue=.false.)
          if (mergerTreeRegridTimes) then
             !@ <inputParameter>
             !@   <name>mergerTreeRegridDumpTrees</name>
             !@   <defaultValue>false</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Specifies whether or not to dump merger trees as they are regridded.
             !@   </description>
             !@   <type>boolean</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeRegridDumpTrees',mergerTreeRegridDumpTrees,defaultValue=.false.)
             !@ <inputParameter>
             !@   <name>mergerTreeRegridCount</name>
             !@   <defaultValue>false</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Number of points in time to use when regridding merger trees.
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeRegridCount',mergerTreeRegridCount,defaultValue=100)
             !@ <inputParameter>
             !@   <name>mergerTreeRegridStartExpansionFactor</name>
             !@   <defaultValue>false</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Starting expansion factor to use when regridding merger trees.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             if (mergerTreeRegridCount < 2) call Galacticus_Error_Report('Merger_Tree_Regrid_Time','mergerTreeRegridCount > 2 is required')
             call Get_Input_Parameter('mergerTreeRegridStartExpansionFactor',mergerTreeRegridStartExpansionFactor,defaultValue=0.1d0)
             !@ <inputParameter>
             !@   <name>mergerTreeRegridEndExpansionFactor</name>
             !@   <defaultValue>false</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Ending expansion factor to use when regridding merger trees.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeRegridEndExpansionFactor',mergerTreeRegridEndExpansionFactor,defaultValue=1.0d0)
             !@ <inputParameter>
             !@   <name>mergerTreeRegridSpacing</name>
             !@   <defaultValue>false</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Type of spacing to use in merger tree regridding (linear or logarithmic).
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeRegridSpacing',mergerTreeRegridSpacingAsText,defaultValue='logarithmic')
             select case (char(mergerTreeRegridSpacingAsText))
             case ("linear")
                mergerTreeRegridSpacing=mergerTreeRegridSpacingLinear
             case ("logarithmic")
                mergerTreeRegridSpacing=mergerTreeRegridSpacingLogarithmic
             case ("log critical density")
                mergerTreeRegridSpacing=mergerTreeRegridSpacingLogCriticalOverdensity
             case default
                call Galacticus_Error_Report('Merger_Tree_Regrid_Time','unrecognized spacing type: '//mergerTreeRegridSpacingAsText)
             end select

             ! Construct array of grid expansion factors.
             call Alloc_Array(mergerTreeRegridTimeGrid,[mergerTreeRegridCount])
             select case (mergerTreeRegridSpacing)
             case (mergerTreeRegridSpacingLinear                )
                mergerTreeRegridTimeGrid=Make_Range(mergerTreeRegridStartExpansionFactor,mergerTreeRegridEndExpansionFactor&
                     &,mergerTreeRegridCount,rangeTypeLinear     )
                ! Convert expansion factors to time.
                do iTime=1,mergerTreeRegridCount
                   mergerTreeRegridTimeGrid(iTime)=Cosmology_Age(mergerTreeRegridTimeGrid(iTime))
                end do
             case (mergerTreeRegridSpacingLogarithmic           )
                mergerTreeRegridTimeGrid=Make_Range(mergerTreeRegridStartExpansionFactor,mergerTreeRegridEndExpansionFactor&
                     &,mergerTreeRegridCount,rangeTypeLogarithmic)
                ! Convert expansion factors to time.
                do iTime=1,mergerTreeRegridCount
                   mergerTreeRegridTimeGrid(iTime)=Cosmology_Age(mergerTreeRegridTimeGrid(iTime))
                end do
             case (mergerTreeRegridSpacingLogCriticalOverdensity)
                ! Build a logarithmic grid in critical overdensity.
                mergerTreeRegridTimeGrid&
                     & =Make_Range(                                                                                        &
                     &              Critical_Overdensity_for_Collapse(Cosmology_Age(mergerTreeRegridStartExpansionFactor)) &
                     &             ,Critical_Overdensity_for_Collapse(Cosmology_Age(mergerTreeRegridEndExpansionFactor  )) &
                     &             ,mergerTreeRegridCount                                                                  &
                     &             ,rangeTypeLogarithmic                                                                   &
                     &            )
                ! Convert critical overdensity to time.
                do iTime=1,mergerTreeRegridCount
                   mergerTreeRegridTimeGrid(iTime)=Time_of_Collapse(mergerTreeRegridTimeGrid(iTime))
                end do
             end select

          end if

          ! Flag that module is initialized.
          regridTimeModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Regrid_Time_Initialize)
    end if

    ! Prune tree if necessary.
    if (mergerTreeRegridTimes) then

       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))

          ! Dump the unprocessed tree if required.
          if (mergerTreeRegridDumpTrees) call Merger_Tree_Dump(                              &
               &                                               currentTree%index,            &
               &                                               currentTree%baseNode        , &
               &                                               backgroundColor    ='white' , &
               &                                               nodeColor          ='black' , &
               &                                               highlightColor     ='black' , &
               &                                               edgeColor          ='black' , &
               &                                               nodeStyle          ='solid' , &
               &                                               highlightStyle     ='filled', &
               &                                               edgeStyle          ='solid' , &
               &                                               labelNodes         =.false. , &
               &                                               scaleNodesByLogMass=.true.  , &
               &                                               edgeLengthsToTimes =.true.    &
               &                                              )

          ! Ensure interpolation accelerator gets reset.
          interpolationReset=.true.

          ! Find the current maximum node index in the tree.
          nodeIndex=0_kind_int8
          thisNode => currentTree%baseNode
          do while (associated(thisNode))
             nodeIndex=max(nodeIndex,thisNode%index())
             call thisNode%walkTree(thisNode)
          end do
          firstNewNode=nodeIndex+1

          ! Walk the tree, locating branches which intersect grid times.
          thisNode => currentTree%baseNode
          do while (associated(thisNode))
             thisBasicComponent => thisNode%basic()

             ! Skip this node if it is the root node.
             if (associated(thisNode%parent)) then
                parentBasicComponent => thisNode%parent%basic()

                ! Get the time of this node and its parent.
                timeNow   =thisBasicComponent  %time()
                timeParent=parentBasicComponent%time()

                ! Get masses of these halos.
                massNow   =thisBasicComponent  %mass()
                massParent=parentBasicComponent%mass()
                if (thisNode%isPrimaryProgenitor()) then
                   ! Remove the mass in any non-primary progenitors - we don't want to include their mass in the estimated mass
                   ! growth rate of this node.
                   childNode => thisNode%parent%firstChild%sibling
                   do while (associated(childNode))
                      childBasicComponent => childNode%basic()
                      massParent          =  massParent-childBasicComponent%mass()
                      childNode           => childNode%sibling
                   end do
                end if

                ! Locate these times in the list of grid times.
                iNow   =Interpolate_Locate(mergerTreeRegridCount,mergerTreeRegridTimeGrid,interpolationAccelerator,timeNow   ,reset=interpolationReset)
                iParent=Interpolate_Locate(mergerTreeRegridCount,mergerTreeRegridTimeGrid,interpolationAccelerator,timeParent,reset=interpolationReset)

                ! For nodes existing precisely at a grid time, ignore this grid point. (These are, typically, nodes which have been created at these points.)
                if (timeNow == mergerTreeRegridTimeGrid(iNow)) iNow=iNow+1

                ! If the branch from node to parent spans one or more grid times, insert new nodes at those points.
                if (iParent > iNow) then
                   ! Create new nodes.
                   allocate(newNodes(iParent-iNow),stat=allocErr)
                   if (allocErr/=0) call Galacticus_Error_Report('Merger_Tree_Regrid_Time','unable to allocate new nodes')
                   do iTime=iNow+1,iParent
                      nodeIndex=nodeIndex+1_kind_int8
                      call newNodes(iTime-iNow)%initialize(nodeIndex)
                   end do
                   ! Assign node properties and build links.
                   do iTime=iNow+1,iParent
                      ! Assign a time and a mass
                      thisBasicComponent => newNodes(iTime-iNow)%basic()
                      call thisBasicComponent%timeSet(                              mergerTreeRegridTimeGrid(iTime)                              )
                      call thisBasicComponent%massSet(massNow+(massParent-massNow)*(mergerTreeRegridTimeGrid(iTime)-timeNow)/(timeParent-timeNow))
                      ! Link to child node.
                      if (iTime > iNow+1 ) newNodes(iTime-iNow)%firstChild => newNodes(iTime-iNow-1)
                      ! Link to parent node.
                      if (iTime < iParent) newNodes(iTime-iNow)%parent     => newNodes(iTime-iNow+1)
                   end do
                   ! Link final node to the parent.
                   newNodes(iParent-iNow)%parent  => thisNode%parent
                   ! Link final node sibling to current node sibling.
                   newNodes(iParent-iNow)%sibling => thisNode%sibling
                   ! Link the parent to the final node.
                   if (thisNode%isPrimaryProgenitor()) then
                      ! Node is the main progenitor of its parent, so simply replace it with the final node in our list.
                      thisNode%parent%firstChild  => newNodes(iParent-iNow)
                   else
                      ! Node is not the main progenitor of its parent, so find the child node that has it as a sibling.
                      childNode => thisNode%parent%firstChild
                      do while (.not.associated(childNode%sibling,thisNode))
                         childNode => childNode%sibling
                      end do
                      childNode%sibling => newNodes(iParent-iNow)
                   end if
                   ! Link the child of the first node to the node being processed.
                   newNodes(1)%firstChild  => thisNode
                   ! Nullify any sibling of the node being processed.
                   thisNode%sibling => null()
                   ! Link the parent of the node being processed to the first node of the list.
                   thisNode%parent  => newNodes(1)
                   ! Erase the node list.
                   deallocate(newNodes)
                end if

             end if

             ! Step to the next node.
             call thisNode%walkTree(thisNode)

          end do

          ! Dump the intermediate tree if required.
          if (mergerTreeRegridDumpTrees) then
             allocate(highlightNodes(nodeIndex-firstNewNode+2))
             highlightNodes(1)=currentTree%baseNode%index()
             do nodeIndex=1,nodeIndex-firstNewNode+1
                highlightNodes(nodeIndex+1)=firstNewNode+nodeIndex-1
             end do
             call Merger_Tree_Dump(                                    &
                  &                currentTree%index,                  &
                  &                currentTree%baseNode              , &
                  &                highlightNodes     =highlightNodes, &
                  &                backgroundColor    ='white'       , &
                  &                nodeColor          ='black'       , &
                  &                highlightColor     ='black'       , &
                  &                edgeColor          ='black'       , &
                  &                nodeStyle          ='solid'       , &
                  &                highlightStyle     ='filled'      , &
                  &                edgeStyle          ='dotted'      , &
                  &                labelNodes         =.false.       , &
                  &                scaleNodesByLogMass=.true.        , &
                  &                edgeLengthsToTimes =.true.          &
                  &               )
             deallocate(highlightNodes)
          end if

          ! Walk the tree removing nodes not at grid times.
          thisNode => currentTree%baseNode
          do while (associated(thisNode))
             thisBasicComponent => thisNode%basic()

             ! Record the next node to walk to.
             call thisNode%walkTree(nextNode)

             ! Get the time for this node.
             timeNow=thisBasicComponent%time()

             ! Find the closest time in the new time grid.
             iNow   =Interpolate_Locate(mergerTreeRegridCount,mergerTreeRegridTimeGrid,interpolationAccelerator,timeNow,reset=interpolationReset,closest=.true.)

             ! If this node does not lie precisely on the grid then remove it.
             if (associated(thisNode%parent) .and. timeNow /= mergerTreeRegridTimeGrid(iNow)) then
                if (thisNode%isPrimaryProgenitor()) then
                   ! Handle primary progenitor nodes.
                   if (associated(thisNode%firstChild)) then
                      ! Handle primary progenitors with children
                      childNode => thisNode%firstChild
                      ! Assign all children a parent that is the parent of the current node.
                      do while (associated(childNode))
                         childNode%parent => thisNode %parent
                         if (.not.associated(childNode%sibling)) then
                            childNode%sibling => thisNode%sibling
                            childNode             => null()
                         else
                            childNode             => childNode%sibling
                         end if
                      end do
                      ! Assign the current node's parent a child that is the child of the current node.
                      thisNode%parent%firstChild => thisNode%firstChild
                   else
                      ! Handle primary nodes with no children - simply make the parents main progenitor the sibling of the current node.
                      thisNode%parent%firstChild => thisNode%sibling
                   end if
                else
                   ! Handle non-primary nodes.
                   if (associated(thisNode%firstChild)) then
                      ! Handle nod-primary nodes with children.
                      ! Assign all children a parent that is the parent of the current node.
                      childNode => thisNode%firstChild
                      do while (associated(childNode))
                         childNode%parent => thisNode %parent
                         if (.not.associated(childNode%sibling)) then
                            childNode%sibling => thisNode%sibling
                            childNode => null()
                         else
                            childNode            => childNode%sibling
                         end if
                      end do
                      ! Find which sibling points the current node and link in the children of the current node.
                      siblingNode => thisNode%parent%firstChild
                      do while (.not.associated(siblingNode%sibling,thisNode))
                         siblingNode => siblingNode%sibling
                      end do
                      siblingNode%sibling => thisNode%firstChild
                   else
                      ! Handle non-primary nodes with no children - just snip it out of the sibling list.
                      siblingNode => thisNode%parent%firstChild
                      do while (.not.associated(siblingNode%sibling,thisNode))
                         siblingNode => siblingNode%sibling
                      end do
                      siblingNode%sibling => thisNode%sibling
                   end if
                end if

                ! Destroy the node.
                call thisNode%destroy()

             end if

             ! Step to the next node.
             thisNode => nextNode

          end do

          ! Clean up interpolation objects.
          call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)

          ! Dump the processed tree if required.
          if (mergerTreeRegridDumpTrees) call Merger_Tree_Dump(                                     &
               &                                             currentTree%index,                  &
               &                                             currentTree%baseNode              , &
               &                                             backgroundColor    ='white'       , &
               &                                             nodeColor          ='black'       , &
               &                                             highlightColor     ='black'       , &
               &                                             edgeColor          ='black'       , &
               &                                             nodeStyle          ='solid'       , &
               &                                             highlightStyle     ='filled'      , &
               &                                             edgeStyle          ='solid'       , &
               &                                             labelNodes         =.false.       , &
               &                                             scaleNodesByLogMass=.true.        , &
               &                                             edgeLengthsToTimes =.true.          &
               &                                            )
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if

    return
  end subroutine Merger_Tree_Regrid_Time

end module Merger_Trees_Regrid_Times
