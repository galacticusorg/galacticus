!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{RST
  Implements a node operator class that records the complete evolutionary trajectory of selected nodes.
  !!}

  use            :: Node_Trajectory_Events, only : enumerationNodeTrajectoryEventType
  use            :: Kind_Numbers          , only : kind_int8
  use, intrinsic :: ISO_C_Binding         , only : c_size_t

  !![
  <nodeOperator name="nodeOperatorOutputTrajectory" docformat="rst">
   <description>
   A node operator class which records the complete evolutionary trajectory of one or more selected nodes to the Galacticus
   output file. This is intended as a debugging tool---it makes it possible to inspect how the properties of a specific galaxy
   evolve *between* the regular output times, without having to add temporary output statements to the code.

   The nodes to follow are specified either by matching pairs of ``[treeIndex]`` and ``[nodeIndex]`` values, or by ``[uniqueID]``
   values. Note that ``uniqueID`` values are assigned from a global counter which is incremented as nodes are created, and which
   is strided by MPI rank---consequently they are *not* reproducible between runs which use differing numbers of OpenMP threads
   or MPI processes, and selection by ``[treeIndex]``/``[nodeIndex]`` should normally be preferred.

   Records are made at whichever of the available events are enabled via the ``[recordNodeInitialize]``,
   ``[recordDifferentialEvolutionPostStep]``, ``[recordDifferentialEvolutionPost]``, ``[recordNodePromote]``,
   ``[recordNodesMerge]``, and ``[recordGalaxiesMerge]`` parameters. ``[recordDifferentialEvolutionPostStep]`` records after every
   step accepted by the ODE solver and so gives the finest-grained view of the evolution, while
   ``[recordDifferentialEvolutionPost]`` records once per call to the node evolver (i.e. at the end of each timestep). The event
   responsible for each record is available in the output via the
   :galacticus-class:`nodePropertyExtractorTrajectoryEvent` property extractor.

   The properties recorded are entirely determined by the :galacticus-class:`mergerTreeOutputterClass` object supplied as a
   sub-parameter of this operator---in particular, by the :galacticus-class:`nodePropertyExtractorClass` objects given to a
   :galacticus-class:`mergerTreeOutputterStandard` outputter. That outputter's ``[outputsGroupName]`` parameter should be set to
   something other than ``Outputs`` so that trajectory data are kept separate from the regular outputs. If no outputter is
   given, a :galacticus-class:`mergerTreeOutputterStandard` writing the ``time``, ``indicesTree``, ``nodeIndices``, and
   ``trajectoryEvent`` properties to the ``nodeTrajectories`` group is used. Note that the outputter is looked for only among
   this operator's own sub-parameters---the outputter used for regular outputs is deliberately *not* inherited, as that would
   cause trajectory records to be written into the regular ``Outputs`` group. Records are written to
   ``&lt;outputsGroupName&gt;/Trajectory&lt;N&gt;/nodeData``, where ``N`` is the position of the node in the list of tracked nodes if
   ``[groupPerNode]`` is true, or ``1`` (so that all tracked nodes share a single table) if it is false.

   Because each record is appended to the output file immediately there is a significant cost per record. This class is intended
   for use with a small number of nodes; the ``[countRecordsMaximum]`` parameter provides a safety limit on the number of records
   made for each tracked node.

   When a node is promoted its components are moved to its parent and the node itself is destroyed. This class detects such
   promotions and re-targets its tracking to the parent node, so that the trajectory of a galaxy continues uninterrupted across
   promotions.

   Note that property extractors which modify the state of the node, or which assume that they are being evaluated at a regular
   output time, are not safe to use with this class.
   </description>
   <deepCopy>
    <ignore variables="mergerTreeOutputter_"/>
   </deepCopy>
   <assignment>
    <functionClass variables="mergerTreeOutputter_"/>
   </assignment>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorOutputTrajectory
     !!{RST
     A node operator class that records the complete evolutionary trajectory of selected nodes.
     !!}
     private
     ! Held as an unlimited polymorphic pointer, and manipulated via the globally-callable functions in
     ! `Merger_Tree_Outputters_Utilities`. Naming `mergerTreeOutputterClass` here would make `Nodes_Operators` depend on
     ! `Merger_Tree_Outputters`, which depends (via `mergerTreeOutputterFullState`) on `Merger_Tree_Construction`, which depends
     ! back on `Nodes_Operators` - a circular module dependency.
     class  (*                       ), pointer                   :: mergerTreeOutputter_            => null()
     integer(kind_int8               ), allocatable, dimension(:) :: indexTree                                , indexNode                          , &
          &                                                          idUnique                                 , idTracked
     integer(c_size_t                ), allocatable, dimension(:) :: countRecords
     integer(c_size_t                )                            :: countRecordsMaximum
     logical                                                      :: selectByUniqueID                         , groupPerNode
     logical                                                      :: recordNodeInitialize                     , recordNodePromote                  , &
          &                                                          recordNodesMerge                         , recordGalaxiesMerge                , &
          &                                                          recordDifferentialEvolutionPost          , recordDifferentialEvolutionPostStep
   contains
     !![
     <methods docformat="rst">
       <method description="Return the position of ``node`` in the list of tracked nodes, or zero if it is not tracked." method="slot"  />
       <method description="Record the current properties of ``node`` if it is being tracked."                           method="record"/>
     </methods>
     !!]
     final     ::                                  outputTrajectoryDestructor
     procedure :: deepCopy                      => outputTrajectoryDeepCopy
     procedure :: deepCopyReset                 => outputTrajectoryDeepCopyReset
     procedure :: deepCopyFinalize              => outputTrajectoryDeepCopyFinalize
     procedure :: slot                          => outputTrajectorySlot
     procedure :: record                        => outputTrajectoryRecord
     procedure :: nodeInitialize                => outputTrajectoryNodeInitialize
     procedure :: differentialEvolutionPostStep => outputTrajectoryDifferentialEvolutionPostStep
     procedure :: differentialEvolutionPost     => outputTrajectoryDifferentialEvolutionPost
     procedure :: nodePromote                   => outputTrajectoryNodePromote
     procedure :: nodesMerge                    => outputTrajectoryNodesMerge
     procedure :: galaxiesMerge                 => outputTrajectoryGalaxiesMerge
  end type nodeOperatorOutputTrajectory

  interface nodeOperatorOutputTrajectory
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorOutputTrajectory` node operator class.
     !!}
     module procedure outputTrajectoryConstructorParameters
     module procedure outputTrajectoryConstructorInternal
  end interface nodeOperatorOutputTrajectory

contains

  function outputTrajectoryConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorOutputTrajectory` node operator class which takes a parameter set as
    input.
    !!}
    use :: Error           , only : Error_Report
    use :: Functions_Global, only : mergerTreeOutputterConstruct_, mergerTreeOutputterDestruct_
    use :: Input_Parameters, only : inputParameter               , inputParameters
    implicit none
    type   (nodeOperatorOutputTrajectory)                              :: self
    type   (inputParameters             ), intent(inout)               :: parameters
    class  (*                           ), pointer                     :: mergerTreeOutputter_
    integer(kind_int8                   ), allocatable  , dimension(:) :: indexTree                      , indexNode                          , &
         &                                                                idUnique
    integer(c_size_t                    )                              :: countNodes                     , countRecordsMaximum
    logical                                                            :: selectByUniqueID               , groupPerNode                       , &
         &                                                                recordNodeInitialize           , recordNodePromote                  , &
         &                                                                recordNodesMerge               , recordGalaxiesMerge                , &
         &                                                                recordDifferentialEvolutionPost, recordDifferentialEvolutionPostStep

    ! Determine how nodes are to be selected.
    selectByUniqueID=parameters%isPresent('uniqueID')
    if (selectByUniqueID) then
       if (parameters%isPresent('treeIndex').or.parameters%isPresent('nodeIndex'))                                      &
            & call Error_Report(                                                                                        &
            &                   "specify nodes using either [uniqueID], or [treeIndex] and [nodeIndex], but not both"// &
            &                   {introspection:location}                                                                &
            &                  )
       countNodes=parameters%count('uniqueID')
       allocate(idUnique (countNodes))
       allocate(indexTree(0         ))
       allocate(indexNode(0         ))
       !![
       <inputParameter docformat="rst">
         <name>uniqueID</name>
         <source>parameters</source>
         <variable>idUnique</variable>
         <type>integer</type>
         <cardinality>1..*</cardinality>
         <description>
         The unique IDs of the nodes whose evolutionary trajectories are to be recorded. Note that unique IDs are not
         reproducible between runs which use differing numbers of OpenMP threads or MPI processes.
         </description>
       </inputParameter>
       !!]
    else
       if (.not.(parameters%isPresent('treeIndex').and.parameters%isPresent('nodeIndex')))                &
            & call Error_Report(                                                                          &
            &                   "specify nodes using either [uniqueID], or [treeIndex] and [nodeIndex]"// &
            &                   {introspection:location}                                                  &
            &                  )
       countNodes=parameters%count('treeIndex')
       if (parameters%count('nodeIndex') /= countNodes)                                                   &
            & call Error_Report(                                                                          &
            &                   "[treeIndex] and [nodeIndex] must contain the same number of values"   // &
            &                   {introspection:location}                                                  &
            &                  )
       allocate(indexTree(countNodes))
       allocate(indexNode(countNodes))
       allocate(idUnique (0         ))
       !![
       <inputParameter docformat="rst">
         <name>treeIndex</name>
         <source>parameters</source>
         <variable>indexTree</variable>
         <type>integer</type>
         <cardinality>1..*</cardinality>
         <description>
         The indices of the merger trees containing the nodes whose evolutionary trajectories are to be recorded.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>nodeIndex</name>
         <source>parameters</source>
         <variable>indexNode</variable>
         <type>integer</type>
         <cardinality>1..*</cardinality>
         <description>
         The indices, within the corresponding merger tree, of the nodes whose evolutionary trajectories are to be recorded.
         </description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter docformat="rst">
      <name>groupPerNode</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, the trajectory of each tracked node is written to its own HDF5 group. If false, the trajectories of all tracked
      nodes are written to a single group (in which case the nodes may be distinguished using, for example, the ``indicesTree``
      and ``nodeIndices`` property extractors).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countRecordsMaximum</name>
      <source>parameters</source>
      <defaultValue>100000_c_size_t</defaultValue>
      <description>
      The maximum number of records to make for each tracked node. Once this many records have been made for a node no further
      records are made for it. Set to zero for no limit.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordNodeInitialize</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made when a tracked node is initialized prior to evolution.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordDifferentialEvolutionPostStep</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made after every step of the ODE solver which is accepted for a tracked node. This gives the
      finest-grained view of the evolution, but also generates the largest volume of output.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordDifferentialEvolutionPost</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made each time a tracked node has been evolved to the end of its current timestep.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordNodePromote</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made immediately before a tracked node is promoted to its parent.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordNodesMerge</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made when a tracked node merges with another node.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordGalaxiesMerge</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a record is made when the galaxy in a tracked node merges with another galaxy.
      </description>
    </inputParameter>
    !!]
    call mergerTreeOutputterConstruct_(parameters,mergerTreeOutputter_)
    self=nodeOperatorOutputTrajectory(                                     &
         &                            indexTree                          , &
         &                            indexNode                          , &
         &                            idUnique                           , &
         &                            selectByUniqueID                   , &
         &                            groupPerNode                       , &
         &                            countRecordsMaximum                , &
         &                            recordNodeInitialize               , &
         &                            recordDifferentialEvolutionPostStep, &
         &                            recordDifferentialEvolutionPost    , &
         &                            recordNodePromote                  , &
         &                            recordNodesMerge                   , &
         &                            recordGalaxiesMerge                , &
         &                            mergerTreeOutputter_                 &
         &                           )
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="mergerTreeOutputter"/>
    !!]
    call mergerTreeOutputterDestruct_(mergerTreeOutputter_)
    return
  end function outputTrajectoryConstructorParameters

  function outputTrajectoryConstructorInternal(indexTree,indexNode,idUnique,selectByUniqueID,groupPerNode,countRecordsMaximum,recordNodeInitialize,recordDifferentialEvolutionPostStep,recordDifferentialEvolutionPost,recordNodePromote,recordNodesMerge,recordGalaxiesMerge,mergerTreeOutputter_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorOutputTrajectory` node operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (nodeOperatorOutputTrajectory)                              :: self
    integer(kind_int8                   ), intent(in   ), dimension(:) :: indexTree                      , indexNode                          , &
         &                                                                idUnique
    integer(c_size_t                    ), intent(in   )               :: countRecordsMaximum
    logical                              , intent(in   )               :: selectByUniqueID               , groupPerNode                       , &
         &                                                                recordNodeInitialize           , recordNodePromote                  , &
         &                                                                recordNodesMerge               , recordGalaxiesMerge                , &
         &                                                                recordDifferentialEvolutionPost, recordDifferentialEvolutionPostStep
    class  (*                           ), intent(in   ), target       :: mergerTreeOutputter_
    integer(c_size_t                    )                              :: countNodes
    !![
    <constructorAssign variables="indexTree, indexNode, idUnique, selectByUniqueID, groupPerNode, countRecordsMaximum, recordNodeInitialize, recordDifferentialEvolutionPostStep, recordDifferentialEvolutionPost, recordNodePromote, recordNodesMerge, recordGalaxiesMerge, *mergerTreeOutputter_"/>
    !!]

    if (selectByUniqueID) then
       countNodes=size(idUnique ,kind=c_size_t)
    else
       countNodes=size(indexTree,kind=c_size_t)
    end if
    if (countNodes < 1_c_size_t) call Error_Report('at least one node must be specified'//{introspection:location})
    allocate(self%idTracked   (countNodes))
    allocate(self%countRecords(countNodes))
    self%idTracked   =-1_kind_int8
    self%countRecords= 0_c_size_t
    return
  end function outputTrajectoryConstructorInternal

  subroutine outputTrajectoryDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorOutputTrajectory` node operator class.
    !!}
    use :: Functions_Global, only : mergerTreeOutputterDestruct_
    implicit none
    type(nodeOperatorOutputTrajectory), intent(inout) :: self

    if (associated(self%mergerTreeOutputter_)) call mergerTreeOutputterDestruct_(self%mergerTreeOutputter_)
    return
  end subroutine outputTrajectoryDestructor

  function outputTrajectorySlot(self,node) result(slot)
    !!{RST
    Return the position of ``node`` in the list of tracked nodes, or zero if it is not being tracked.
    !!}
    implicit none
    integer(c_size_t                    )                :: slot
    class  (nodeOperatorOutputTrajectory), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(c_size_t                    )                :: i
    integer(kind_int8                   )                :: idUnique_

    slot     =0_c_size_t
    idUnique_=node%uniqueID()
    ! First look for a tracked node which we have already resolved to a unique ID.
    do i=1_c_size_t,size(self%idTracked,kind=c_size_t)
       if (self%idTracked(i) == idUnique_) then
          slot=i
          return
       end if
    end do
    ! Otherwise, attempt to resolve any as-yet-unresolved tracked node. Note that resolution is attempted only once for each
    ! tracked node. This matters because, once our tracked node has been promoted (at which point we re-target to its parent),
    ! some other node may acquire the original indices of our tracked node - this happens, for example, if the
    ! `nodeOperatorIndexShift` node operator is in use. Resolving only once ensures that such a node can not be mistaken for the
    ! node that we are tracking.
    do i=1_c_size_t,size(self%idTracked,kind=c_size_t)
       if (self%idTracked(i) >= 0_kind_int8) cycle
       if (self%selectByUniqueID) then
          if (idUnique_ /= self%idUnique(i)) cycle
       else
          if (.not.associated(node%hostTree)) cycle
          if (node%hostTree%index /= self%indexTree(i) .or. node%index() /= self%indexNode(i)) cycle
       end if
       self%idTracked(i)=idUnique_
       slot             =i
       return
    end do
    return
  end function outputTrajectorySlot

  subroutine outputTrajectoryRecord(self,node,event)
    !!{RST
    Record the current properties of ``node``, if it is one of the nodes that we are tracking.
    !!}
    use :: Functions_Global      , only : mergerTreeOutputterOutputTrajectory_
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventCurrent
    implicit none
    class  (nodeOperatorOutputTrajectory      ), intent(inout) :: self
    type   (treeNode                          ), intent(inout) :: node
    type   (enumerationNodeTrajectoryEventType), intent(in   ) :: event
    integer(c_size_t                          )                :: slot_, indexOutput

    slot_=self%slot(node)
    if (slot_ == 0_c_size_t) return
    ! Impose the limit on the number of records made for this node.
    if     (                                                             &
         &   self%countRecordsMaximum        >  0_c_size_t               &
         &  .and.                                                        &
         &   self%countRecords       (slot_) >= self%countRecordsMaximum &
         & ) return
    ! Record the event responsible for this record so that it can be extracted by the `trajectoryEvent` property extractor during
    ! the output which we trigger below.
    nodeTrajectoryEventCurrent=event%ID
    ! Output the node. Each tracked node is given its own output group (i.e. its own index) unless all nodes are to share a
    ! single group.
    if (self%groupPerNode) then
       indexOutput=slot_
    else
       indexOutput=1_c_size_t
    end if
    call mergerTreeOutputterOutputTrajectory_(self%mergerTreeOutputter_,node,indexOutput)
    self%countRecords(slot_)=self%countRecords(slot_)+1_c_size_t
    return
  end subroutine outputTrajectoryRecord

  subroutine outputTrajectoryNodeInitialize(self,node)
    !!{RST
    Record the properties of a tracked node as it is initialized.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventNodeInitialize
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout), target :: self
    type (treeNode                    ), intent(inout), target :: node

    if (self%recordNodeInitialize) call self%record(node,nodeTrajectoryEventNodeInitialize)
    return
  end subroutine outputTrajectoryNodeInitialize

  subroutine outputTrajectoryDifferentialEvolutionPostStep(self,node,status)
    !!{RST
    Record the properties of a tracked node after each accepted step of the ODE solver.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventDifferentialEvolutionPostStep
    implicit none
    class  (nodeOperatorOutputTrajectory), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer                              , intent(inout) :: status
    !$GLC attributes unused :: status

    ! Note that `status` is deliberately left unchanged - we make no modification to the state of the node, so the ODE solver
    ! must not be asked to reset.
    if (self%recordDifferentialEvolutionPostStep) call self%record(node,nodeTrajectoryEventDifferentialEvolutionPostStep)
    return
  end subroutine outputTrajectoryDifferentialEvolutionPostStep

  subroutine outputTrajectoryDifferentialEvolutionPost(self,node)
    !!{RST
    Record the properties of a tracked node once it has been evolved to the end of its timestep.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventDifferentialEvolutionPost
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (self%recordDifferentialEvolutionPost) call self%record(node,nodeTrajectoryEventDifferentialEvolutionPost)
    return
  end subroutine outputTrajectoryDifferentialEvolutionPost

  subroutine outputTrajectoryNodePromote(self,node)
    !!{RST
    Record the properties of a tracked node as it is promoted, and re-target our tracking to its parent.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventNodePromote
    implicit none
    class  (nodeOperatorOutputTrajectory), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(c_size_t                    )                :: slot_

    slot_=self%slot(node)
    if (slot_ == 0_c_size_t) return
    if (self%recordNodePromote) call self%record(node,nodeTrajectoryEventNodePromote)
    ! The components of this node are about to be moved to its parent, after which this node will be destroyed. The parent has a
    ! different unique ID (and, unless the `nodeOperatorIndexShift` operator is in use, a different index), so re-target our
    ! tracking to the parent to allow the trajectory to continue uninterrupted.
    if (associated(node%parent)) self%idTracked(slot_)=node%parent%uniqueID()
    return
  end subroutine outputTrajectoryNodePromote

  subroutine outputTrajectoryNodesMerge(self,node)
    !!{RST
    Record the properties of a tracked node as it merges with another node.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventNodesMerge
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (self%recordNodesMerge) call self%record(node,nodeTrajectoryEventNodesMerge)
    return
  end subroutine outputTrajectoryNodesMerge

  subroutine outputTrajectoryGalaxiesMerge(self,node)
    !!{RST
    Record the properties of a tracked node as its galaxy merges with another galaxy.
    !!}
    use :: Node_Trajectory_Events, only : nodeTrajectoryEventGalaxiesMerge
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (self%recordGalaxiesMerge) call self%record(node,nodeTrajectoryEventGalaxiesMerge)
    return
  end subroutine outputTrajectoryGalaxiesMerge

  subroutine outputTrajectoryDeepCopyReset(self)
    !!{RST
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : mergerTreeOutputterDeepCopyReset_
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout) :: self

    self%copiedSelf => null()
    if (associated(self%mergerTreeOutputter_)) call mergerTreeOutputterDeepCopyReset_(self%mergerTreeOutputter_)
    return
  end subroutine outputTrajectoryDeepCopyReset

  subroutine outputTrajectoryDeepCopyFinalize(self)
    !!{RST
    Finalize a deep copy of the object.
    !!}
    use :: Functions_Global, only : mergerTreeOutputterDeepCopyFinalize_
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout) :: self

    if (associated(self%mergerTreeOutputter_)) call mergerTreeOutputterDeepCopyFinalize_(self%mergerTreeOutputter_)
    return
  end subroutine outputTrajectoryDeepCopyFinalize

  subroutine outputTrajectoryDeepCopy(self,destination)
    !!{RST
    Perform a deep copy of the object.
    !!}
    use :: Error           , only : Error_Report
    use :: Functions_Global, only : mergerTreeOutputterDeepCopy_
    implicit none
    class(nodeOperatorOutputTrajectory), intent(inout), target :: self
    class(nodeOperatorClass           ), intent(inout)         :: destination

    call self%deepCopy_(destination)
    select type (destination)
    type is (nodeOperatorOutputTrajectory)
       nullify(destination%mergerTreeOutputter_)
       if (associated(self%mergerTreeOutputter_)) then
          allocate(destination%mergerTreeOutputter_,mold=self%mergerTreeOutputter_)
          call mergerTreeOutputterDeepCopy_(self%mergerTreeOutputter_,destination%mergerTreeOutputter_)
       end if
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine outputTrajectoryDeepCopy
