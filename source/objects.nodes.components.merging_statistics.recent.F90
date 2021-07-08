!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!!{
Contains a module which implements the recent merging statistics component.
!!}

module Node_Component_Merging_Statistics_Recent
  !!{
  Implements the recent merging statistics component.
  !!}
  use            :: Dark_Matter_Halo_Scales                      , only : darkMatterHaloScaleClass
  use, intrinsic :: ISO_C_Binding                                , only : c_size_t
  use            :: Node_Component_Merging_Statistics_Recent_Data, only : Node_Component_Merging_Statistics_Recent_Count
  use            :: Output_Times                                 , only : outputTimesClass
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Recent_Merger_Tree_Init , Node_Component_Merging_Statistics_Recent_Node_Merger        , &
       &    Node_Component_Merging_Statistics_Recent_Output_Count     , Node_Component_Merging_Statistics_Recent_Output             , &
       &    Node_Component_Merging_Statistics_Recent_Thread_Initialize, Node_Component_Merging_Statistics_Recent_Thread_Uninitialize, &
       &    Node_Component_Merging_Statistics_Recent_Initialize       , Node_Component_Merging_Statistics_Recent_Output_Names 

  !![
  <component>
   <class>mergingStatistics</class>
   <name>recent</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>recentMajorMergerCount</name>
      <type>integer</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault modules="Node_Component_Merging_Statistics_Recent_Data" count="Node_Component_Merging_Statistics_Recent_Count()">0</classDefault>
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  class(outputTimesClass        ), pointer :: outputTimes_
  !$omp threadprivate(darkMatterHaloScale_,outputTimes_)

  ! Parameters controlling the statistics gathered.
  double precision                             :: nodeMajorMergerFraction                     , nodeRecentMajorMergerInterval
  integer                                      :: nodeRecentMajorMergerIntervalType
  integer          , parameter                 :: nodeRecentMajorMergerIntervalTypeAbsolute =0
  integer          , parameter                 :: nodeRecentMajorMergerIntervalTypeDynamical=1
  logical                                      :: nodeRecentMajorMergerFromInfall

  ! Initialization array for count of recent mergers.
  integer(c_size_t)                            :: outputCount
  integer          , allocatable, dimension(:) :: zeroCount

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Initialize(parameters_)
    !!{
    Initializes the recent merging statistics component.
    !!}
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : char                   , var_str        , varying_string
    use :: Input_Parameters  , only : inputParameter         , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(varying_string )                :: nodeRecentMajorMergerIntervalTypeText

    !![
    <inputParameter>
      <name>nodeMajorMergerFraction</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>nodeRecentMajorMergerInterval</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The time interval used to define ``recent'' mergers in the {\normalfont \ttfamily recent} merging statistics component. This parameter is in units of Gyr if {\normalfont \ttfamily [nodeRecentMajorMergerIntervalType]}$=${\normalfont \ttfamily absolute}, or in units of the halo dynamical time if {\normalfont \ttfamily [nodeRecentMajorMergerIntervalType]}$=${\normalfont \ttfamily dynmical}.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>nodeRecentMajorMergerIntervalType</name>
      <defaultValue>var_str('dynamical')</defaultValue>
      <description>Specifies the units for the {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} parameter. If set to {\normalfont \ttfamily absolute} then {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} is given in Gyr, while if set to {\normalfont \ttfamily dynamical} {\normalfont \ttfamily [nodeRecentMajorMergerInterval]} is given in units of the halo dynamical time.</description>
      <source>parameters_</source>
      <variable>nodeRecentMajorMergerIntervalTypeText</variable>
    </inputParameter>
    !!]
    select case (char(nodeRecentMajorMergerIntervalTypeText))
    case ("absolute" )
       nodeRecentMajorMergerIntervalType=nodeRecentMajorMergerIntervalTypeAbsolute
    case ("dynamical")
       nodeRecentMajorMergerIntervalType=nodeRecentMajorMergerIntervalTypeDynamical
    case default
       call Galacticus_Error_Report('[nodeRecentMajorMergerIntervalType] has unrecognized value'//{introspection:location})
    end select
    !![
    <inputParameter>
      <name>nodeRecentMajorMergerFromInfall</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether ``recent'' for satellite galaxies is measured from the current time, or from the time at which they were last isolated.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node recent merging flow statistics module.
    !!}
    use :: Events_Hooks                                 , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes                             , only : defaultMergingStatisticsComponent
    use :: Input_Parameters                             , only : inputParameter                   , inputParameters
    use :: Memory_Management                            , only : allocateArray
    use :: Node_Component_Merging_Statistics_Recent_Data, only : mergingStatisticsRecentCount
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultMergingStatisticsComponent%recentIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters_"/>
       !!]
       call nodePromotionEvent%attach(defaultMergingStatisticsComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentMergingStatisticsRecent')
       !$omp critical (Node_Component_Merging_Statistics_Recent_Thread_Initialize)
       if (.not.allocated(zeroCount)) then
          mergingStatisticsRecentCount=outputTimes_%count()
          call allocateArray(zeroCount,[mergingStatisticsRecentCount])
          zeroCount=0
       end if
       !$omp end critical (Node_Component_Merging_Statistics_Recent_Thread_Initialize)
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Thread_Uninitialize()
    !!{
    Uninitializes the tree node recent merging flow statistics module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent
    implicit none

    if (defaultMergingStatisticsComponent%recentIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       <objectDestructor name="outputTimes_"        />
       !!]
       call nodePromotionEvent%detach(defaultMergingStatisticsComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Thread_Uninitialize

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Merger_Tree_Init</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init(node)
    !!{
    Initialize the merging statistics component by creating components in nodes.
    !!}
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent, nodeComponentMergingStatistics, nodeComponentMergingStatisticsRecent, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%recentIsActive()) return
    ! Create a merger statistics component and initialize it.
    mergingStatistics => node%mergingStatistics(autoCreate=.true.)
    select type (mergingStatistics)
    class is (nodeComponentMergingStatisticsRecent)
       call mergingStatistics%recentMajorMergerCountSet(zeroCount)
    end select
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Merger_Tree_Init

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Node_Merger</unitName>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Node_Merger(node)
    !!{
    Record any major merger of {\normalfont \ttfamily node}.
    !!}
    use            :: Galacticus_Error, only : Galacticus_Error_Report
    use            :: Galacticus_Nodes, only : defaultMergingStatisticsComponent, nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    type            (treeNode                      ), intent(inout)          :: node
    type            (treeNode                      ), pointer                :: nodeDescendent
    class           (nodeComponentMergingStatistics), pointer                :: mergingStatisticsParent
    class           (nodeComponentBasic            ), pointer                :: basicDescendentParent  , basicParent, &
         &                                                                      basic
    integer                                         , dimension(outputCount) :: mergerIncrement
    integer         (c_size_t                      )                         :: i
    double precision                                                         :: recentTimeInterval     , timeBase

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
             recentTimeInterval=nodeRecentMajorMergerInterval*darkMatterHaloScale_%dynamicalTimescale(node)
          case default
             call Galacticus_Error_Report('unrecognized time interval type'//{introspection:location})
          end select
          if (nodeRecentMajorMergerFromInfall) then
             if (node%parent%isSatellite()) then
                timeBase=basicParent%timeLastIsolated()
             else
                timeBase=outputTimes_%time(i)
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
             timeBase=outputTimes_%time(i)
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

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the node merger time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentMergingStatistics, treeNode
    implicit none
    class(*                             ), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent, mergingStatistics
    !$GLC attributes unused :: self
    
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
  end subroutine nodePromotion

  !![
  <mergerTreeOutputNames>
   <unitName>Node_Component_Merging_Statistics_Recent_Output_Names</unitName>
   <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  </mergerTreeOutputNames>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Output_Names(node,integerProperty,integerProperties,doubleProperty,doubleProperties,time)
    !!{
    Set names of black hole properties to be written to the \glc\ output file.
    !!}
    use :: Galacticus_Nodes                  , only : treeNode
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    type            (treeNode             )              , intent(inout) :: node
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty   , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    !$GLC attributes unused :: time, doubleProperty, doubleProperties

    if (Node_Component_Merging_Statistics_Recent_Matches(node)) then
       integerProperty=integerProperty+1
       integerProperties(integerProperty)%name     ='mergingStatisticsRecentMajorMergerCount'
       integerProperties(integerProperty)%comment  ='Number of major mergers occuring in a recent time interval.'
       integerProperties(integerProperty)%unitsInSI=0.0d0
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output_Names

  !![
  <mergerTreeOutputPropertyCount>
   <unitName>Node_Component_Merging_Statistics_Recent_Output_Count</unitName>
   <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  </mergerTreeOutputPropertyCount>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Output_Count(node,integerPropertyCount,doublePropertyCount,time)
    !!{
    Account for the number of black hole properties to be written to the the \glc\ output file.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !$GLC attributes unused :: doublePropertyCount, time

    if (Node_Component_Merging_Statistics_Recent_Matches(node)) integerPropertyCount=integerPropertyCount+1
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output_Count

  !![
  <mergerTreeOutputTask>
   <unitName>Node_Component_Merging_Statistics_Recent_Output</unitName>
   <sortName>Node_Component_Merging_Statistics_Recent_Output</sortName>
  </mergerTreeOutputTask>
  !!]
  subroutine Node_Component_Merging_Statistics_Recent_Output(node,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,instance)
    !!{
    Store black hole properties in the \glc\ output file buffers.
    !!}
    use :: Galacticus_Nodes                  , only : nodeComponentMergingStatistics, treeNode
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger         , outputPropertyDouble
    implicit none
    double precision                                                        , intent(in   ) :: time
    type            (treeNode                      )                        , intent(inout) :: node
    integer                                                                 , intent(inout) :: doubleBufferCount , doubleProperty , &
         &                                                                                     integerBufferCount, integerProperty
    type            (outputPropertyInteger         ), dimension(:          ), intent(inout) :: integerProperties
    type            (outputPropertyDouble          ), dimension(:          ), intent(inout) :: doubleProperties
    type            (multiCounter                  )                        , intent(inout) :: instance
    class           (nodeComponentMergingStatistics)                        , pointer       :: mergingStatistics
    integer                                         , dimension(outputCount)                :: mergerIncrement
    !$GLC attributes unused :: doubleBufferCount, doubleProperty, doubleProperties, instance

    if (Node_Component_Merging_Statistics_Recent_Matches(node)) then
       ! Store the properties.
       mergingStatistics => node             %mergingStatistics     ()
       mergerIncrement   =  mergingStatistics%recentMajorMergerCount()
       integerProperty=integerProperty+1
       integerProperties(integerProperty)%scalar(integerBufferCount)=mergerIncrement(outputTimes_%index(time,findClosest=.true.))
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Recent_Output

  logical function Node_Component_Merging_Statistics_Recent_Matches(node)
    !!{
    Return true if the black hole component of {\normalfont \ttfamily node} is a match to the standard implementation.
    !!}
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent, nodeComponentMergingStatistics, nodeComponentMergingStatisticsRecent, treeNode
    implicit none
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentMergingStatistics), pointer       :: mergingStatistics

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
