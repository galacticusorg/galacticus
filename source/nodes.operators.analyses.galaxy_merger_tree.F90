!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node operator class that records properties of galaxies in every galaxy's merger tree.
  !!}

  use :: Node_Property_Extractors, only : multiExtractorList
  
  !![
  <nodeOperator name="nodeOperatorGalaxyMergerTree">
   <description>A node operator class that records properties of galaxies in every galaxy's merger tree.</description>
   <linkedList type="multiExtractorList" variable="extractors" next="next" object="extractor_" objectType="nodePropertyExtractorClass" module="Node_Property_Extractors"/>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGalaxyMergerTree
     !!{
     A node operator class that records properties of galaxies in every galaxy's merger tree.
     !!}
     private
     type   (multiExtractorList), pointer :: extractors          => null()
     integer                              :: uniqueIDBranchTipID          , nodeIndexID , &
          &                                  countID                      , timeID      , &
          &                                  propertyID                   , mergeIndexID, &
          &                                  mergeTargetID                , timeMergeID , &
          &                                  countExtractors
     double precision                     :: timeStep
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize galaxy merger tree meta-properties."/>
       <method method="record"     description="Record a step in the merger tree."             />
     </methods>
     !!]
     final     ::                              galaxyMergerTreeDestructor
     procedure :: autoHook                  => galaxyMergerTreeAutoHook
     procedure :: initialize                => galaxyMergerTreeInitialize
     procedure :: record                    => galaxyMergerTreeRecord
     procedure :: differentialEvolutionPost => galaxyMergerTreeDifferentialEvolutionPost
     procedure :: deepCopy                  => galaxyMergerTreeDeepCopy
  end type nodeOperatorGalaxyMergerTree
  
  interface nodeOperatorGalaxyMergerTree
     !!{
     Constructors for the \refClass{nodeOperatorGalaxyMergerTree} node operator class.
     !!}
     module procedure galaxyMergerTreeConstructorParameters
     module procedure galaxyMergerTreeConstructorInternal
  end interface nodeOperatorGalaxyMergerTree
  
contains

  function galaxyMergerTreeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorGalaxyMergerTree} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorGalaxyMergerTree)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    type   (multiExtractorList          ), pointer       :: extractor_
    integer                                              :: i

    !![
    <inputParameter>
      <name>timeStep</name>
      <variable>self%timestep</variable>
      <description>The minimum timestep at which to record galaxy properties.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self      %extractors => null()
    extractor_            => null()
    do i=1,parameters%copiesCount('nodePropertyExtractor')
       if (associated(extractor_)) then
          allocate(extractor_%next)
          extractor_ => extractor_%next
       else
          allocate(self%extractors)
          extractor_ => self%extractors
       end if
       !![
       <objectBuilder class="nodePropertyExtractor" name="extractor_%extractor_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    call self%initialize()
    return
  end function galaxyMergerTreeConstructorParameters

  function galaxyMergerTreeConstructorInternal(timeStep,extractors) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorGalaxyMergerTree} node operator class.
    !!}
    implicit none
    type            (nodeOperatorGalaxyMergerTree)                         :: self
    double precision                              , intent(in   )          :: timeStep
    type            (multiExtractorList          ), intent(in   ), target  :: extractors
    type            (multiExtractorList          )               , pointer :: extractor_
    !![
    <constructorAssign variables="timeStep"/>
    !!]

    self      %extractors => extractors
    extractor_            => extractors
    do while (associated(extractor_))
       !![
       <referenceCountIncrement owner="extractor_" object="extractor_"/>
       !!]
       extractor_ => extractor_%next
    end do
    call self%initialize()
    return
  end function galaxyMergerTreeConstructorInternal

  subroutine galaxyMergerTreeInitialize(self)
    !!{
    Initialize meta-properties needed for galaxy merger trees.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    class(nodeOperatorGalaxyMergerTree), intent(inout):: self
    type (multiExtractorList          ), pointer      :: extractor_

    !![
    <addMetaProperty component="basic" name="nodeUniqueIDBranchTip"       type="longInteger" id="self%uniqueIDBranchTipID"                           isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeProperty"                       id="self%propertyID"          rank="1" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeTime"                           id="self%timeID"              rank="1" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeTimeMerge"                      id="self%timeMergeID"         rank="1" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeNodeIndex"   type="longInteger" id="self%nodeIndexID"         rank="1"                  isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeCount"       type="longInteger" id="self%countID"             rank="1"                  isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeMergeIndex"  type="longInteger" id="self%mergeIndexID"        rank="1"                  isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeMergeTarget" type="longInteger" id="self%mergeTargetID"       rank="1"                  isCreator="yes"/>
    !!]
    ! Count extractors.
    self%countExtractors =  0
    extractor_           => self%extractors
    do while (associated(extractor_))
       self%countExtractors =  self      %countExtractors+1
       extractor_           => extractor_%next
    end do
    return
  end subroutine galaxyMergerTreeInitialize

  subroutine galaxyMergerTreeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorGalaxyMergerTree), intent(inout) :: self

    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='galaxyMergerTree')
    return
  end subroutine galaxyMergerTreeAutoHook
  
  subroutine galaxyMergerTreeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorGalaxyMergerTree} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorGalaxyMergerTree), intent(inout) :: self
    type(multiExtractorList          ), pointer       :: extractor_, extractorNext

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    if (associated(self%extractors)) then
       extractor_ => self%extractors
       do while (associated(extractor_))
          extractorNext => extractor_%next
          !![
          <objectDestructor name="extractor_%extractor_"/>
          !!]
          deallocate(extractor_)
          extractor_ => extractorNext
       end do
    end if
    return
  end subroutine galaxyMergerTreeDestructor

  subroutine galaxyMergerTreeDifferentialEvolutionPost(self,node)
    !!{
    Operate on the node after differential evolution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorGalaxyMergerTree), intent(inout)             :: self
    type            (treeNode                    ), intent(inout)             :: node
    double precision                              , allocatable, dimension(:) :: times
    class           (nodeComponentBasic          ), pointer                   :: basic
    logical                                                                   :: record, isNew

    basic  => node %basic                    (           )
    times  =  basic%floatRank1MetaPropertyGet(self%timeID)
    isNew  =  size(times) == 0
    record =  isNew .or. basic%time() >= times(size(times))+self%timeStep
    if (record) call self%record(node)
    return
  end subroutine galaxyMergerTreeDifferentialEvolutionPost
  
  subroutine galaxyMergerTreeRecord(self,node)
    !!{
    Record a new time in the merger tree.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    implicit none
    class           (nodeOperatorGalaxyMergerTree), intent(inout)               :: self
    type            (treeNode                    ), intent(inout)               :: node
    double precision                              , allocatable  , dimension(:) :: times      , timesNew     , &
         &                                                                         properties , propertiesNew
    integer         (c_size_t                    ), allocatable  , dimension(:) :: nodeIndices, counts
    class           (nodeComponentBasic          ), pointer                     :: basic
    type            (multiExtractorList          ), pointer                     :: extractor_
    logical                                                                     :: isNew
    integer                                                                     :: i
 
    basic => node %basic                    (           )
    times =  basic%floatRank1MetaPropertyGet(self%timeID)
    isNew =  size(times) == 0
    nodeIndices=basic%longIntegerRank1MetaPropertyGet(self%nodeIndexID)
    counts     =basic%longIntegerRank1MetaPropertyGet(self%countID    )
    properties =basic%      floatRank1MetaPropertyGet(self%propertyID )
    if (isNew) then
       deallocate(nodeIndices   )
       deallocate(counts        )
       allocate  (nodeIndices(1))
       allocate  (counts     (1))
       counts(1)=0_c_size_t
    end if
    nodeIndices(size(counts))=basic%longIntegerRank0MetaPropertyGet(self%uniqueIDBranchTipID)
    counts     (size(counts))=counts(size(counts))+1_c_size_t
    allocate(timesNew     (size(times     )+1                   ))
    allocate(propertiesNew(size(properties)+self%countExtractors))
    if (size(times     ) > 0) timesNew     (1:size(times     ))=times
    if (size(properties) > 0) propertiesNew(1:size(properties))=properties
    timesNew(size(times)+1) =  basic%time      ()
    extractor_              => self %extractors
    i                       =  0
    do while (associated(extractor_))
       i=i+1
       select type (nodePropertyExtractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          propertiesNew(size(properties)+i)=nodePropertyExtractor_%extract(node)
       end select
       extractor_=> extractor_%next
    end do
    call basic%longIntegerRank1MetaPropertySet(self%nodeIndexID,nodeIndices  )
    call basic%longIntegerRank1MetaPropertySet(self%countID    ,counts       )
    call basic%floatRank1MetaPropertySet      (self%propertyID ,propertiesNew)
    call basic%floatRank1MetaPropertySet      (self%timeID     ,timesNew     )
    return
  end subroutine galaxyMergerTreeRecord
  
  subroutine satelliteMerger(self,node)
    !!{
    Record galaxy-galaxy merges and combine trees.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (*                               ), intent(inout)               :: self
    type            (treeNode                        ), intent(inout), target       :: node
    type            (treeNode                        )               , pointer      :: nodeTarget
    class           (nodeComponentBasic              )               , pointer      :: basic             , basicTarget
    double precision                                  , allocatable  , dimension(:) :: times             , timesTarget       , &
         &                                                                             timesNew          , properties        , &
         &                                                                             propertiesTarget  , propertiesNew     , &
         &                                                                             timesMerge        , timesMergeTarget  , &
         &                                                                             timesMergeNew
    integer         (c_size_t                        ), allocatable  , dimension(:) :: nodeIndices       , nodeIndicesTarget , &
         &                                                                             nodeIndicesNew    , counts            , &
         &                                                                             countsTarget      , countsNew         , &
         &                                                                             mergeIndices      , mergeIndicesTarget, &
         &                                                                             mergeIndicesNew   , mergeTargets      , &
         &                                                                             mergeTargetsTarget, mergeTargetsNew

    select type (self)
    class is (nodeOperatorGalaxyMergerTree)
       nodeTarget  => node      %mergesWith()
       basic       => node      %basic     ()
       basicTarget => nodeTarget%basic     ()
       ! Add new records of the galaxies immediately prior to merging.
       call self%record(node      )
       call self%record(nodeTarget)
       ! Combine the records from both galaxies.
       times             =basic      %      floatRank1MetaPropertyGet(self%       timeID)
       timesMerge        =basic      %      floatRank1MetaPropertyGet(self%  timeMergeID)
       nodeIndices       =basic      %longIntegerRank1MetaPropertyGet(self%  nodeIndexID)
       mergeIndices      =basic      %longIntegerRank1MetaPropertyGet(self% mergeIndexID)
       mergeTargets      =basic      %longIntegerRank1MetaPropertyGet(self%mergeTargetID)
       counts            =basic      %longIntegerRank1MetaPropertyGet(self%      countID)
       properties        =basic      %      floatRank1MetaPropertyGet(self%   propertyID)
       timesTarget       =basicTarget%      floatRank1MetaPropertyGet(self%       timeID)
       timesMergeTarget  =basicTarget%      floatRank1MetaPropertyGet(self%  timeMergeID)
       nodeIndicesTarget =basicTarget%longIntegerRank1MetaPropertyGet(self%  nodeIndexID)
       mergeIndicesTarget=basicTarget%longIntegerRank1MetaPropertyGet(self% mergeIndexID)
       mergeTargetsTarget=basicTarget%longIntegerRank1MetaPropertyGet(self%mergeTargetID)
       countsTarget      =basicTarget%longIntegerRank1MetaPropertyGet(self      %countID)
       propertiesTarget  =basicTarget%      floatRank1MetaPropertyGet(self%   propertyID)
       allocate(timesNew       (size(times       )+size(timesTarget       )  ))
       allocate(nodeIndicesNew (size(nodeIndices )+size(nodeIndicesTarget )  ))
       allocate(countsNew      (size(counts      )+size(countsTarget      )  ))
       allocate(propertiesNew  (size(properties  )+size(propertiesTarget  )  ))
       allocate(mergeIndicesNew(size(mergeIndices)+size(mergeIndicesTarget)+1))
       allocate(mergeTargetsNew(size(mergeTargets)+size(mergeTargetsTarget)+1))
       allocate(timesMergeNew  (size(timesMerge  )+size(timesMergeTarget  )+1))
       timesNew       (1:size(times       ))=times       ; timesNew       (size(times       )+1:size(timesNew       )  )=timesTarget
       propertiesNew  (1:size(properties  ))=properties  ; propertiesNew  (size(properties  )+1:size(propertiesNew  )  )=propertiesTarget
       nodeIndicesNew (1:size(nodeIndices ))=nodeIndices ; nodeIndicesNew (size(nodeIndices )+1:size(nodeIndicesNew )  )=nodeIndicesTarget
       countsNew      (1:size(counts      ))=counts      ; countsNew      (size(counts      )+1:size(countsNew      )  )=countsTarget
       timesMergeNew  (1:size(timesMerge  ))=timesMerge  ; timesMergeNew  (size(timesMerge  )+1:size(timesMergeNew  )-1)=timesMergeTarget  ; timesMergeNew  (size(timesMergeNew  ))=basicTarget%time                           (                        )
       mergeIndicesNew(1:size(mergeIndices))=mergeIndices; mergeIndicesNew(size(mergeIndices)+1:size(mergeIndicesNew)-1)=mergeIndicesTarget; mergeIndicesNew(size(mergeIndicesNew))=basic      %longIntegerRank0MetaPropertyGet(self%uniqueIDBranchTipID)
       mergeTargetsNew(1:size(mergeTargets))=mergeTargets; mergeTargetsNew(size(mergeTargets)+1:size(mergeTargetsNew)-1)=mergeTargetsTarget; mergeTargetsNew(size(mergeTargetsNew))=basicTarget%longIntegerRank0MetaPropertyGet(self%uniqueIDBranchTipID)
       call basicTarget%longIntegerRank1MetaPropertySet(self%nodeIndexID  ,nodeIndicesNew )
       call basicTarget%longIntegerRank1MetaPropertySet(self%countID      ,countsNew      )
       call basicTarget%      floatRank1MetaPropertySet(self%propertyID   ,propertiesNew  )
       call basicTarget%      floatRank1MetaPropertySet(self%timeID       ,timesNew       )
       call basicTarget%      floatRank1MetaPropertySet(self%timeMergeID  ,timesMergeNew  )
       call basicTarget%longIntegerRank1MetaPropertySet(self%mergeIndexID ,mergeIndicesNew)
       call basicTarget%longIntegerRank1MetaPropertySet(self%mergeTargetID,mergeTargetsNew)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger

  subroutine galaxyMergerTreeDeepCopy(self,destination)
    !!{
    Deep copy a {\normalfont \ttfamily nodeOperatorGalaxyMergerTree} object.
    !!}
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount, nodePropertyExtractorGalaxyMergerTree_
    implicit none
    class  (nodeOperatorGalaxyMergerTree), intent(inout), target :: self
    class  (nodeOperatorClass           ), intent(inout)         :: destination
    type   (multiExtractorList          ), pointer               :: extractor_
    integer                                                      :: i

    call self%deepCopy_(destination)
    select type (destination)
    class is (nodeOperatorGalaxyMergerTree)
       nodePropertyExtractorGalaxyMergerTreeCount=destination%countExtractors
       allocate(nodePropertyExtractorGalaxyMergerTree_(destination%countExtractors))
       extractor_ => self %extractors
       i          =  0
       do while (associated(extractor_))
          i                                                    =  i                    +1
          nodePropertyExtractorGalaxyMergerTree_(i)%extractor_ => extractor_%extractor_
          extractor_                                           => extractor_%next
       end do
    end select
    return
  end subroutine galaxyMergerTreeDeepCopy
