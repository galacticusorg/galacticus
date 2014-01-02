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

!% Contains a module which implements the recent merging statistics component.

module Node_Component_Merging_Statistics_Recent
  !% Implements the recent merging statistics component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Recent_Merger_Tree_Init, Node_Component_Merging_Statistics_Recent_Node_Merger , &
       &    Node_Component_Merging_Statistics_Recent_Node_Promotion  , Node_Component_Merging_Statistics_Recent_Output_Names, &
       &    Node_Component_Merging_Statistics_Recent_Output_Count    , Node_Component_Merging_Statistics_Recent_Output

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>recent</name>
  !#  <isDefault>no</isDefault>
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
  double precision                            :: nodeMajorMergerFraction                           , nodeRecentMajorMergerInterval
  integer                                     :: nodeRecentMajorMergerIntervalType
  integer         , parameter                 :: nodeRecentMajorMergerIntervalTypeAbsolute =0
  integer         , parameter                 :: nodeRecentMajorMergerIntervalTypeDynamical=1
  logical                                     :: nodeRecentMajorMergerFromInfall

  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized                         =.false.

  ! Initialization array for count of recent mergers.
  integer                                     :: outputCount
  integer         , allocatable, dimension(:) :: zeroCount

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
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeMajorMergerFraction',nodeMajorMergerFraction,defaultValue=0.25d0)
       !@ <inputParameter>
       !@   <name>nodeRecentMajorMergerInterval</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The time interval used to define ``recent'' mergers in the {\tt recent} merging statistics component. This parameter is in units of Gyr if {\tt [nodeRecentMajorMergerIntervalType]}$={\tt absolute}, or in units of the halo dynamical time if {\tt [nodeRecentMajorMergerIntervalType]}$={\tt dynmical}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeRecentMajorMergerInterval',nodeRecentMajorMergerInterval,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>nodeRecentMajorMergerIntervalType</name>
       !@   <defaultValue>dynamical</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies the units for the {\tt [nodeRecentMajorMergerInterval]} parameter. If set to {\tt absolute} then {\tt [nodeRecentMajorMergerInterval]} is given in Gyr, while if set to {\tt dynamical} {\tt [nodeRecentMajorMergerInterval]} is given in units of the halo dynamical time.
       !@   </description>
       !@   <type>real</type>
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
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeRecentMajorMergerFromInfall',nodeRecentMajorMergerFromInfall,defaultValue=.false.)
       ! Determine the number of output times.
       outputCount=Galacticus_Output_Time_Count()
       call Alloc_Array(zeroCount,[outputCount])
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
  subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init(thisNode)
    !% Initialize the merging statistics component by creating components in nodes.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentMergingStatistics)               , pointer :: thisMergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    ! Ensure that the module is initialized.
    call Node_Component_Merging_Statistics_Recent_Initialize()

    ! Create a merger statistics component and initialize it.
    thisMergingStatistics => thisNode%mergingStatistics(autoCreate=.true.)
    select type (thisMergingStatistics)
    class is (nodeComponentMergingStatisticsRecent)
       call thisMergingStatistics%recentMajorMergerCountSet(zeroCount)
    end select
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Merging_Statistics_Recent_Node_Merger(thisNode)
    !% Record any major merger of {\tt thisNode}.
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    use Galacticus_Output_Times
    implicit none
    type            (treeNode                      ), intent(inout)         , pointer :: thisNode
    type            (treeNode                      )                        , pointer :: descendentNode
    class           (nodeComponentMergingStatistics)                        , pointer :: parentMergingStatistics
    class           (nodeComponentBasic            )                        , pointer :: descendentParentBasic  , parentBasic, &
         &                                                                               thisBasic
    integer                                         , dimension(outputCount)          :: mergerIncrement
    integer                                                                           :: i
    double precision                                                                  :: recentTimeInterval     , timeBase

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    thisBasic               => thisNode       %basic            ()
    parentBasic             => thisNode%parent%basic            ()
    parentMergingStatistics => thisNode%parent%mergingStatistics()
    ! Record the merger time if this is a major merger.
    if (thisBasic%mass() >= nodeMajorMergerFraction*parentBasic%mass()) then
       ! Iterate over output times and check if this merger is sufficient close to them to be counted.
       mergerIncrement=0
       do i=1,outputCount
          select case (nodeRecentMajorMergerIntervalType)
          case (nodeRecentMajorMergerIntervalTypeAbsolute)
             recentTimeInterval=nodeRecentMajorMergerInterval
          case (nodeRecentMajorMergerIntervalTypeDynamical)
             recentTimeInterval=nodeRecentMajorMergerInterval*Dark_Matter_Halo_Dynamical_Timescale(thisNode)
          case default
             call Galacticus_Error_Report('Node_Component_Merging_Statistics_Recent_Node_Merger','unrecognized time interval type')
          end select
          if (nodeRecentMajorMergerFromInfall) then
             if (thisNode%parent%isSatellite()) then
                timeBase=parentBasic%timeLastIsolated()
             else
                timeBase=Galacticus_Output_Time(i)
                descendentNode => thisNode%parent
                do while (associated(descendentNode))
                   if (descendentNode%isPrimaryProgenitor()) then
                      descendentNode => descendentNode%parent
                   else
                      if (associated(descendentNode%parent)) then
                         descendentParentBasic => descendentNode%parent%basic()
                         timeBase=min(timeBase,descendentParentBasic%time())
                      end if
                      exit
                   end if
                end do
             end if
          else
             timeBase=Galacticus_Output_Time(i)
          end if
          if     (                                                            &
               &   thisBasic%time() <= timeBase                               &
               &  .and.                                                       &
               &   thisBasic%time() >  timeBase-nodeRecentMajorMergerInterval &
               & ) mergerIncrement(i)=1
       end do
       call parentMergingStatistics%recentMajorMergerCountSet(parentMergingStatistics%recentMajorMergerCount()+mergerIncrement)
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Node_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Merging_Statistics_Recent_Node_Promotion(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the node merger time.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentMergingStatistics)               , pointer :: parentMergingStatistics, thisMergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return

    ! Get the merging statistics components.
    parentMergingStatistics => thisNode%parent%mergingStatistics()
    thisMergingStatistics   => thisNode       %mergingStatistics()
    call      thisMergingStatistics%recentMajorMergerCountSet    &
         & (                                                     &
         &     thisMergingStatistics%recentMajorMergerCount   () &
         &  +parentMergingStatistics%recentMajorMergerCount   () &
         & )
    call    parentMergingStatistics%recentMajorMergerCountSet    &
         & (                                                     &
         &  zeroCount                                            &
         & )
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Node_Promotion

  !# <mergerTreeOutputNames>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Output_Names</unitName>
  !#  <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Node_Component_Merging_Statistics_Recent_Output_Names(thisNode,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI ,doubleProperty,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    if (Node_Component_Merging_Statistics_Recent_Matches(thisNode)) then
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
  subroutine Node_Component_Merging_Statistics_Recent_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    if (Node_Component_Merging_Statistics_Recent_Matches(thisNode)) integerPropertyCount=integerPropertyCount+1
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Node_Component_Merging_Statistics_Recent_Output</unitName>
  !#  <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Node_Component_Merging_Statistics_Recent_Output(thisNode,integerProperty,integerBufferCount,integerBuffer&
       &,doubleProperty ,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    double precision                                , intent(in   )                   :: time
    type            (treeNode                      ), intent(inout)         , pointer :: thisNode
    integer                                         , intent(inout)                   :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                                               integerProperty
    integer         (kind=kind_int8                ), intent(inout)                   :: integerBuffer        (:,:)
    double precision                                , intent(inout)                   :: doubleBuffer         (:,:)
    class           (nodeComponentMergingStatistics)                        , pointer :: thisMergingStatistics
    integer                                         , dimension(outputCount)          :: mergerIncrement


    if (Node_Component_Merging_Statistics_Recent_Matches(thisNode)) then
       ! Store the properties.
       thisMergingStatistics => thisNode             %mergingStatistics     ()
       mergerIncrement       =  thisMergingStatistics%recentMajorMergerCount()
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=mergerIncrement(Galacticus_Output_Time_Index(time))
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output

  logical function Node_Component_Merging_Statistics_Recent_Matches(thisNode)
    !% Return true if the black hole component of {\tt thisNode} is a match to the standard implementation.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentMergingStatistics)               , pointer :: thisMergingStatistics

    ! Get the merging statistics component.
    thisMergingStatistics => thisNode%mergingStatistics()
    ! Ensure that it is of the recent class.
    Node_Component_Merging_Statistics_Recent_Matches=.false.
    select type (thisMergingStatistics)
    class is (nodeComponentMergingStatisticsRecent)
       Node_Component_Merging_Statistics_Recent_Matches=.true.
    type  is (nodeComponentMergingStatistics        )
       Node_Component_Merging_Statistics_Recent_Matches=defaultMergingStatisticsComponent%recentIsActive()
    end select
    return
  end function Node_Component_Merging_Statistics_Recent_Matches

end module Node_Component_Merging_Statistics_Recent
