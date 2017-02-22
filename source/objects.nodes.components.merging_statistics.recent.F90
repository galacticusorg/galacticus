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

!% Contains a module which implements the recent merging statistics component.

module Node_Component_Merging_Statistics_Recent
  !% Implements the recent merging statistics component.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Recent_Merger_Tree_Init, Node_Component_Merging_Statistics_Recent_Node_Merger , &
       &    Node_Component_Merging_Statistics_Recent_Node_Promotion  , Node_Component_Merging_Statistics_Recent_Output_Names, &
       &    Node_Component_Merging_Statistics_Recent_Output_Count    , Node_Component_Merging_Statistics_Recent_Output

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>recent</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>recentMajorMergerCount</name>
  !#     <type>integer</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault modules="Galacticus_Output_Times" count="Galacticus_Output_Time_Count()">0</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Parameters controlling the statistics gathered.
  double precision                             :: nodeMajorMergerFraction                           , nodeRecentMajorMergerInterval
  integer                                      :: nodeRecentMajorMergerIntervalType
  integer          , parameter                 :: nodeRecentMajorMergerIntervalTypeAbsolute =0
  integer          , parameter                 :: nodeRecentMajorMergerIntervalTypeDynamical=1
  logical                                      :: nodeRecentMajorMergerFromInfall

  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized                         =.false.

  ! Initialization array for count of recent mergers.
  integer(c_size_t)                            :: outputCount
  integer          , allocatable, dimension(:) :: zeroCount

contains

  subroutine Node_Component_Merging_Statistics_Recent_Initialize()
    !% Initializes the recent merging statistics component.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Memory_Management
    implicit none
    type(varying_string) :: nodeRecentMajorMergerIntervalTypeText

    ! Test whether module is already initialize.
    !$omp critical (Node_Component_Merging_Statistics_Recent_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>nodeMajorMergerFraction</name>
       !@   <defaultValue>0.25</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeMajorMergerFraction',nodeMajorMergerFraction,defaultValue=0.25d0)
       !@ <inputParameter>
       !@   <name>nodeRecentMajorMergerInterval</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The time interval used to define ``recent'' mergers in the {\normalfont \ttfamily recent} merging statistics component. This parameter is in units of Gyr if {\normalfont \ttfamily [nodeRecentMajorMergerIntervalType]}$=${\normalfont \ttfamily absolute}, or in units of the halo dynamical time if {\normalfont \ttfamily [nodeRecentMajorMergerIntervalType]}$=${\normalfont \ttfamily dynmical}.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeRecentMajorMergerInterval',nodeRecentMajorMergerInterval,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>nodeRecentMajorMergerIntervalType</name>
       !@   <defaultValue>dynamical</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies the units for the {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} parameter. If set to {\normalfont \ttfamily absolute} then {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} is given in Gyr, while if set to {\normalfont \ttfamily dynamical} {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} is given in units of the halo dynamical time.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeRecentMajorMergerIntervalType',nodeRecentMajorMergerIntervalTypeText,defaultValue="dynamical")
       select case (char(nodeRecentMajorMergerIntervalTypeText))
       case ("absolute" )
          nodeRecentMajorMergerIntervalType=nodeRecentMajorMergerIntervalTypeAbsolute
       case ("dynamical")
          nodeRecentMajorMergerIntervalType=nodeRecentMajorMergerIntervalTypeDynamical
       case default
          call Galacticus_Error_Report('Node_Component_Merging_Statistics_Recent_Initialize','[nodeRecentMajorMergerIntervalType] has unrecognized value')
       end select
       !@ <inputParameter>
       !@   <name>nodeRecentMajorMergerFromInfall</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether ``recent'' for satellite galaxies is measured from the current time, or from the time at which they were last isolated.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeRecentMajorMergerFromInfall',nodeRecentMajorMergerFromInfall,defaultValue=.false.)
       ! Determine the number of output times.
       outputCount=Galacticus_Output_Time_Count()
       call allocateArray(zeroCount,[outputCount])
       zeroCount=0
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Merging_Statistics_Recent_Initialize)
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init(node)
    !% Initialize the merging statistics component by creating components in nodes.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    ! Ensure that the module is initialized.
    call Node_Component_Merging_Statistics_Recent_Initialize()

    ! Create a merger statistics component and initialize it.
    mergingStatistics => node%mergingStatistics(autoCreate=.true.)
    select type (mergingStatistics)
    class is (nodeComponentMergingStatisticsRecent)
       call mergingStatistics%recentMajorMergerCountSet(zeroCount)
    end select
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Merging_Statistics_Recent_Node_Merger(node)
    !% Record any major merger of {\normalfont \ttfamily node}.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    use Galacticus_Output_Times
    implicit none
    type            (treeNode                      ), intent(inout)         , pointer :: node
    type            (treeNode                      )                        , pointer :: nodeDescendent
    class           (nodeComponentMergingStatistics)                        , pointer :: mergingStatisticsParent
    class           (nodeComponentBasic            )                        , pointer :: basicDescendentParent  , basicParent, &
         &                                                                               basic
    class           (darkMatterHaloScaleClass      )                        , pointer :: darkMatterHaloScale_
    integer                                         , dimension(outputCount)          :: mergerIncrement
    integer         (c_size_t                      )                                  :: i
    double precision                                                                  :: recentTimeInterval     , timeBase

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    basic                   => node       %basic            ()
    basicParent             => node%parent%basic            ()
    mergingStatisticsParent => node%parent%mergingStatistics()
    ! Record the merger time if this is a major merger.
    if (basic%mass() >= nodeMajorMergerFraction*basicParent%mass()) then
       ! Iterate over output times and check if this merger is sufficient close to them to be counted.
       mergerIncrement=0
       do i=1,outputCount
          select case (nodeRecentMajorMergerIntervalType)
          case (nodeRecentMajorMergerIntervalTypeAbsolute)
             recentTimeInterval=nodeRecentMajorMergerInterval
          case (nodeRecentMajorMergerIntervalTypeDynamical)
             darkMatterHaloScale_ => darkMatterHaloScale()
             recentTimeInterval=nodeRecentMajorMergerInterval*darkMatterHaloScale_%dynamicalTimescale(node)
          case default
             call Galacticus_Error_Report('Node_Component_Merging_Statistics_Recent_Node_Merger','unrecognized time interval type')
          end select
          if (nodeRecentMajorMergerFromInfall) then
             if (node%parent%isSatellite()) then
                timeBase=basicParent%timeLastIsolated()
             else
                timeBase=Galacticus_Output_Time(i)
                nodeDescendent => node%parent
                do while (associated(nodeDescendent))
                   if (nodeDescendent%isPrimaryProgenitor()) then
                      nodeDescendent => nodeDescendent%parent
                   else
                      if (associated(nodeDescendent%parent)) then
                         basicDescendentParent => nodeDescendent%parent%basic()
                         timeBase=min(timeBase,basicDescendentParent%time())
                      end if
                      exit
                   end if
                end do
             end if
          else
             timeBase=Galacticus_Output_Time(i)
          end if
          if     (                                                        &
               &   basic%time() <= timeBase                               &
               &  .and.                                                   &
               &   basic%time() >  timeBase-nodeRecentMajorMergerInterval &
               & ) mergerIncrement(i)=1
       end do
       call mergingStatisticsParent%recentMajorMergerCountSet(mergingStatisticsParent%recentMajorMergerCount()+mergerIncrement)
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Node_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Merging_Statistics_Recent_Node_Promotion(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the node merger time.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent, mergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    ! Get the merging statistics components.
    mergingStatisticsParent => node%parent%mergingStatistics()
    mergingStatistics       => node       %mergingStatistics()
    call      mergingStatistics%recentMajorMergerCountSet        &
         & (                                                     &
         &   mergingStatistics      %recentMajorMergerCount   () &
         &  +mergingStatisticsParent%recentMajorMergerCount   () &
         & )
    call    mergingStatisticsParent%recentMajorMergerCountSet    &
         & (                                                     &
         &  zeroCount                                            &
         & )
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Node_Promotion

  !# <mergerTreeOutputNames>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Output_Names</unitName>
  !#  <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Node_Component_Merging_Statistics_Recent_Output_Names(node,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI ,doubleProperty,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: node
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: time, doubleProperty, doublePropertyComments, doublePropertyNames, doublePropertyUnitsSI
    
    if (Node_Component_Merging_Statistics_Recent_Matches(node)) then
       !@ <outputPropertyGroup>
       !@   <name>mergingStatistics</name>
       !@   <description>Statistics on mergers</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>mergingStatisticsRecentMajorMergerCount</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Number of major mergers occuring in a recent time interval.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>mergingStatistics</group>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='mergingStatisticsRecentMajorMergerCount'
       integerPropertyComments(integerProperty)='Number of major mergers occuring in a recent time interval.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Output_Count</unitName>
  !#  <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Node_Component_Merging_Statistics_Recent_Output_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: node
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: doublePropertyCount, time
    
    if (Node_Component_Merging_Statistics_Recent_Matches(node)) integerPropertyCount=integerPropertyCount+1
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Output</unitName>
  !#  <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Node_Component_Merging_Statistics_Recent_Output(node,integerProperty,integerBufferCount,integerBuffer&
       &,doubleProperty ,doubleBufferCount,doubleBuffer,time,instance)
    !% Store black hole properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galacticus_Output_Times
    use Multi_Counters
    implicit none
    double precision                                , intent(in   )                   :: time
    type            (treeNode                      ), intent(inout)         , pointer :: node
    integer                                         , intent(inout)                   :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                                               integerProperty
    integer         (kind=kind_int8                ), intent(inout)                   :: integerBuffer        (:,:)
    double precision                                , intent(inout)                   :: doubleBuffer         (:,:)
    type            (multiCounter                  ), intent(inout)                   :: instance
    class           (nodeComponentMergingStatistics)                        , pointer :: mergingStatistics
    integer                                         , dimension(outputCount)          :: mergerIncrement
    !GCC$ attributes unused :: doubleBufferCount, doubleProperty, doubleBuffer, instance

    if (Node_Component_Merging_Statistics_Recent_Matches(node)) then
       ! Store the properties.
       mergingStatistics => node             %mergingStatistics     ()
       mergerIncrement       =  mergingStatistics%recentMajorMergerCount()
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=mergerIncrement(Galacticus_Output_Time_Index(time,findClosest=.true.))
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output

  logical function Node_Component_Merging_Statistics_Recent_Matches(node)
    !% Return true if the black hole component of {\normalfont \ttfamily node} is a match to the standard implementation.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatistics

    ! Get the merging statistics component.
    mergingStatistics => node%mergingStatistics()
    ! Ensure that it is of the recent class.
    Node_Component_Merging_Statistics_Recent_Matches=.false.
    select type (mergingStatistics)
    class is (nodeComponentMergingStatisticsRecent)
       Node_Component_Merging_Statistics_Recent_Matches=.true.
    type  is (nodeComponentMergingStatistics        )
       Node_Component_Merging_Statistics_Recent_Matches=defaultMergingStatisticsComponent%recentIsActive()
    end select
    return
  end function Node_Component_Merging_Statistics_Recent_Matches

end module Node_Component_Merging_Statistics_Recent
