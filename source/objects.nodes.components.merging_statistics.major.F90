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

!% Contains a module which implements the major merging statistics component.

module Node_Component_Merging_Statistics_Major
  !% Implements the major merging statistics component.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Major_Satellite_Merging, Node_Component_Merging_Statistics_Major_Node_Promotion, &
      &     Node_Component_Merging_Statistics_Major_Output

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>major</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>majorMergerTime</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains
  
  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Major_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Merging_Statistics_Major_Satellite_Merging(node)
    !% Record any major merger of {\normalfont \ttfamily node}.
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type            (treeNode                      ), intent(inout), pointer      :: node
    class           (nodeComponentMergingStatistics)               , pointer      :: mergingStatistics
    class           (nodeComponentBasic            )               , pointer      :: basicHost          , basic
    type            (treeNode                      )               , pointer      :: nodeHost
    double precision                                , allocatable  , dimension(:) :: majorMergerTimesNew, majorMergerTimes

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%majorIsActive()) return
    ! Record the merger time if this is a major merger.
    if (thisMergerIsMajor) then
       ! Get required components    
       nodeHost          => node    %mergesWith       ()
       basic             => node    %basic            ()
       basicHost         => nodeHost%basic            ()
       mergingStatistics => nodeHost%mergingStatistics(autoCreate=.true.)
       ! Expand the merger time array and record the new time.
       majorMergerTimes=mergingStatistics%majorMergerTime()
       if (allocated(majorMergerTimes)) then
          allocate(majorMergerTimesNew(size(majorMergerTimes)+1))
          majorMergerTimesNew(1:size(majorMergerTimes))=majorMergerTimes
          deallocate(majorMergerTimes)
       else
          allocate(majorMergerTimesNew(1))
       end if
       majorMergerTimesNew(size(majorMergerTimesNew))=basicHost%time()
       call mergingStatistics%majorMergerTimeSet(majorMergerTimesNew)
       deallocate(majorMergerTimesNew)
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Major_Satellite_Merging

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Merging_Statistics_Major_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Merging_Statistics_Major_Node_Promotion(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the node merger time.
    implicit none
    type (treeNode                      ), intent(inout), pointer      :: node
    class(nodeComponentMergingStatistics)               , pointer      :: mergingStatisticsParent, mergingStatistics
    double precision                     , allocatable  , dimension(:) :: timeMajorMerger        , timeMajorMergerParent, &
         &                                                                timeMajorMergerNew

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%majorIsActive()) return
    ! Get the merging statistics components.
    mergingStatisticsParent => node%parent            %mergingStatistics(autoCreate=.true.)
    mergingStatistics       => node                   %mergingStatistics(autoCreate=.true.)
    timeMajorMerger         =  mergingStatistics      %majorMergerTime  (                 )
    timeMajorMergerParent   =  mergingStatisticsParent%majorMergerTime  (                 )
    allocate(timeMajorMergerNew(size(timeMajorMerger)+size(timeMajorMergerParent)))
    timeMajorMergerNew(                      1:size(timeMajorMerger)                            )=timeMajorMerger
    timeMajorMergerNew(size(timeMajorMerger)+1:size(timeMajorMerger)+size(timeMajorMergerParent))=timeMajorMergerParent    
    call mergingStatistics%majorMergerTimeSet(timeMajorMergerNew)
    return
  end subroutine Node_Component_Merging_Statistics_Major_Node_Promotion

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Merging_Statistics_Major_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Merging_Statistics_Major_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Output properties for all black holes in {\normalfont \ttfamily node}.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_HDF5
    use Kind_Numbers
    use ISO_Varying_String
    use String_Handling
    implicit none
    type            (treeNode                      ), intent(inout), pointer      :: node
    integer         (kind=kind_int8                ), intent(in   )               :: treeIndex
    integer         (c_size_t                      ), intent(in   )               :: iOutput
    logical                                         , intent(in   )               :: nodePassesFilter
    class           (nodeComponentMergingStatistics)               , pointer      :: mergingStatistics
    double precision                                , allocatable  , dimension(:) :: majorMergerTimes
    type            (hdf5Object                    )                              :: majorMergersGroup, outputGroup, &
         &                                                                           treeGroup
    type            (varying_string                )                              :: groupName

    ! Return immediately if this class is not active or the filter has not passed.
    if (.not.(defaultMergingStatisticsComponent%majorIsActive().and.nodePassesFilter)) return
    ! Get the major merger times.
    mergingStatistics => node%mergingStatistics()
    majorMergerTimes=mergingStatistics%majorMergerTime()
    ! Return if no major mergers occurred.
    if (.not.allocated(majorMergerTimes).or.size(majorMergerTimes) == 0) return
    ! Open the output group.
    !$omp critical (HDF5_Access)
    majorMergersGroup=galacticusOutputFile%openGroup("majorMergers","Major merger times.")
    groupName="Output"
    groupName=groupName//iOutput
    outputGroup=majorMergersGroup%openGroup(char(groupName),"Major merger times for all trees at each output.")
    groupName="mergerTree"
    groupName=groupName//treeIndex
    treeGroup=outputGroup%openGroup(char(groupName),"Major merger times for all nodes in this tree")
    groupName="node"
    groupName=groupName//node%index()
    call treeGroup%writeDataset(majorMergerTimes,char(groupName),"Major merger times.")
    call treeGroup        %close()
    call outputGroup      %close()
    call majorMergersGroup%close()
    !$omp end critical (HDF5_Access)
    return
  end subroutine Node_Component_Merging_Statistics_Major_Output
  
end module Node_Component_Merging_Statistics_Major
