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
Contains a module which implements an object hierarchy for nodes in merger trees and all of their constituent physical
components.
!!}

module Galacticus_Nodes
  !!{
  Implements an object hierarchy for nodes in merger trees and all of their constituent physical components.
  !!}
  use            :: Abundances_Structure               , only : abundances
  use            :: Chemical_Abundances_Structure      , only : chemicalAbundances
  use            :: Galactic_Structure_Options         , only : enumerationComponentTypeType , enumerationMassTypeType       , enumerationWeightByType
  use            :: Hashes                             , only : doubleHash                   , genericHash
  use            :: Histories                          , only : history                      , longIntegerHistory
  use            :: IO_HDF5                            , only : hdf5Object
  use, intrinsic :: ISO_C_Binding                      , only : c_size_t
  use            :: ISO_Varying_String                 , only : varying_string
  use            :: Kepler_Orbits                      , only : keplerOrbit
  use            :: Kind_Numbers                       , only : kind_int8
  use            :: Locks                              , only : ompLock
  use            :: Mass_Distributions                 , only : massDistributionClass
  use            :: Merger_Trees_Evolve_Deadlock_Status, only : enumerationDeadlockStatusType
  use            :: Numerical_Constants_Astronomical   , only : gigaYear                     , luminosityZeroPointAB         , massSolar             , megaParsec
  use            :: Numerical_Constants_Prefixes       , only : kilo
  use            :: Numerical_Random_Numbers           , only : randomNumberGeneratorClass
  use            :: Stellar_Luminosities_Structure     , only : stellarLuminosities
  use            :: Tensors                            , only : tensorNullR2D3Sym            , tensorRank2Dimension3Symmetric
  private
  public :: nodeClassHierarchyInitialize    , nodeClassHierarchyFinalize, Galacticus_Nodes_Unique_ID_Set, interruptTask   , &
       &    nodeEventBuildFromRaw           , propertyEvaluate          , propertyActive                , propertyInactive, &
       &    massDistributionCalculationReset, massDistributionsLast     , massDistributionsDestroy
  
  type, public :: treeNodeList
     !!{
     Type to give a list of treeNodes.
     !!}
     type(treeNode), pointer :: node => null()
  end type treeNodeList
  
  type, public :: treeNodeLinkedList
     !!{
     Type to give a linked list of treeNodes.
     !!}
     type(treeNode          ), pointer :: node => null()
     type(treeNodeLinkedList), pointer :: next => null()
  end type treeNodeLinkedList

  type, public :: mergerTree
     !!{
     The merger tree object type.
     !!}
     integer         (kind=kind_int8            )                  :: index
     type            (hdf5Object                )                  :: hdf5Group
     double precision                                              :: volumeWeight                    , initializedUntil
     type            (treeNode                  ), pointer         :: nodeBase               => null()
     type            (mergerTree                ), pointer         :: nextTree               => null(), firstTree        => null()
     type            (universe                  ), pointer         :: hostUniverse           => null()
     type            (treeEvent                 ), pointer, public :: event                  => null()
     class           (randomNumberGeneratorClass), pointer         :: randomNumberGenerator_ => null()
     type            (doubleHash                )                  :: properties
   contains
     ! Tree creation/destruction.
     !![
     <methods>
       <method description="Destroys the merger tree, including all nodes and their components." method="destroy" />
       <method description="Returns a pointer to the node with given index in the merger tree, or a null pointer if no such node exists." method="getNode" />
       <method description="Create a {\normalfont \ttfamily treeEvent} object in this tree." method="createEvent" />
       <method description="Remove a {\normalfont \ttfamily treeEvent} from this tree." method="removeEvent" />
       <method description="Return the earliest time in a merger tree." method="earliestTime" />
       <method description="Return the earliest time in an evolving merger tree." method="earliestTimeEvolving" />
       <method description="Return the latest time in a merger tree." method="latestTime" />
       <method description="Return the size (in bytes) of a merger tree." method="sizeOf" />
     </methods>
     !!]
     procedure :: destroy              => Merger_Tree_Destroy
     procedure :: getNode              => Merger_Tree_Node_Get
     procedure :: createEvent          => Merger_Tree_Create_Event
     procedure :: removeEvent          => Merger_Tree_Remove_Event
     procedure :: earliestTime         => Merger_Tree_Earliest_Time
     procedure :: earliestTimeEvolving => Merger_Tree_Earliest_Time_Evolving
     procedure :: latestTime           => Merger_Tree_Latest_Time
     procedure :: sizeOf               => Merger_Tree_Size_Of
  end type mergerTree

  interface mergerTree
     !!{
     Interface to merger tree constructors.
     !!}
     module procedure mergerTreeConstructor
  end interface mergerTree
    
  type, public :: treeEvent
     !!{
     Type for events attached to trees.
     !!}
     private
     integer         (kind=kind_int8)         , public :: ID
     type            (mergerTree    ), pointer, public :: tree => null()
     double precision                         , public :: time
     type            (treeEvent     ), pointer, public :: next => null()
     procedure       (treeEventTask ), pointer, public :: task => null()
  end type treeEvent
  
  ! Interface for tree event tasks.
  abstract interface
     logical function treeEventTask(event,tree,deadlockStatus)
       import treeEvent, mergerTree, enumerationDeadlockStatusType
       class(treeEvent                    ), intent(in   ) :: event
       type (mergerTree                   ), intent(inout) :: tree
       type (enumerationDeadlockStatusType), intent(inout) :: deadlockStatus
     end function treeEventTask
  end interface
  
  type, public :: mergerTreeList
     !!{
     A class used for building linked lists of merger trees.
     !!}
     type(mergerTreeList), pointer :: next => null()
     type(mergerTree    ), pointer :: tree => null()
  end type mergerTreeList

  type, public :: universe
     !!{
     The universe object class.
     !!}
     type   (mergerTreeList), pointer         :: trees         => null()
     logical                                  :: allTreesBuilt =  .false.
     type   (universeEvent ), pointer, public :: event         => null()
     type   (genericHash   )                  :: attributes
     integer(kind_int8     )                  :: uniqueID
     type   (ompLock       )                  :: lock
   contains
     !![
     <methods>
       <method description="Create a {\normalfont \ttfamily treeEvent} object in this universe." method="createEvent"/>
       <method description="Remove a {\normalfont \ttfamily treeEvent} from this universe."      method="removeEvent"/>
       <method description="Pop a {\normalfont \ttfamily mergerTree} from this universe."        method="popTree"    />
       <method description="Pop a {\normalfont \ttfamily mergerTree} from this universe."        method="pushTree"   />
     </methods>
     !!]
     procedure :: createEvent => universeCreateEvent
     procedure :: removeEvent => universeRemoveEvent
     procedure :: popTree     => universePopTree
     procedure :: pushTree    => universePushTree
  end type universe

  interface universe
     !!{
     Interface to universe constructors.
     !!}
     module procedure universeConstructor
  end interface universe
    
  type, public :: universeEvent
     !!{
     Type for events attached to universes.
     !!}
     private
     integer         (kind=kind_int8   )         , public :: ID
     type            (universe         ), pointer, public :: universe => null()
     double precision                            , public :: time
     type            (universeEvent    ), pointer, public :: next     => null()
     procedure       (universeEventTask), pointer, public :: task     => null()
     class           (*                ), pointer, public :: creator  => null()
  end type universeEvent

  ! Interface for universe event tasks.
  abstract interface
     logical function universeEventTask(event,universe_)
       import universeEvent, universe
       class  (universeEvent), intent(in   ) :: event
       type   (universe     ), intent(inout) :: universe_
     end function universeEventTask
  end interface

  ! Meta-property types.
  type :: floatRank1MetaProperty
     !!{
     Type used to store float rank-1 meta-properties.
     !!}
     double precision, allocatable, dimension(:) :: values
  end type floatRank1MetaProperty

  type :: integerRank1MetaProperty
     !!{
     Type used to store integer rank-1 meta-properties.
     !!}
     integer, allocatable, dimension(:) :: values
  end type integerRank1MetaProperty
  
  type :: longIntegerRank1MetaProperty
     !!{
     Type used to store long integer rank-1 meta-properties.
     !!}
     integer(kind_int8), allocatable, dimension(:) :: values
  end type longIntegerRank1MetaProperty
  
  ! Zero dimension arrays to be returned as defaults.
  double precision                                   , dimension(0) :: nullDouble1d

  ! Labels for function mapping reduction types.
  integer                         , parameter, public               :: reductionSummation  =1
  integer                         , parameter, public               :: reductionProduct    =2

  ! Unique ID counter.
  integer         (kind=kind_int8)                                  :: uniqueIdCount       =0

  ! Event ID counter.
  integer         (kind=kind_int8)                                  :: eventID             =0

  ! Universe unique ID counter.
  integer         (kind=kind_int8)                                  :: universeUniqueIdCount=0

  ! Enumeration for active/inactive properties.
  integer, parameter, public :: propertyTypeAll       =0
  integer, parameter, public :: propertyTypeActive    =1
  integer, parameter, public :: propertyTypeInactive  =2
  integer, parameter, public :: propertyTypeNumerics  =3
  integer, parameter, public :: propertyTypeNone      =4

  ! Enumeration for analytically/numerically-solved properties.
  integer, parameter, public :: solutionTypeNumerical =0
  integer, parameter, public :: solutionTypeAnalytical=1
  
  ! State for rate computations.
  integer           , public :: rateComputeState    =propertyTypeActive
  !$omp threadprivate(rateComputeState)

  ! Memoized massDistributions
  type :: massDistributionArray
     private
     integer(kind_int8                   )          :: uniqueID          =  -huge(kind_int8)
     type   (enumerationComponentTypeType)          :: componentType
     type   (enumerationMassTypeType     )          :: massType
     type   (enumerationWeightByType     )          :: weightBy
     integer                                        :: weightIndex
     class  (massDistributionClass       ), pointer :: massDistribution_ =>  null(         )
  end type massDistributionArray
  integer                       , parameter                         :: massDistributionsCount=20
  integer                                                           :: massDistributionsLast = 0
  type   (massDistributionArray), dimension(massDistributionsCount) :: massDistributions__
  !$omp threadprivate(massDistributions__,massDistributionsLast)
  
  ! Define a constructor for treeNodes.
  interface treeNode
     module procedure Tree_Node_Constructor
  end interface treeNode

  ! Include node methods.
  !![
  <include directive="component" type="component">
  !!]
  include 'objects.nodes.components.inc'
  !![
  </include>
  !!]

  !
  ! Nodes functions.
  subroutine Galacticus_Nodes_Unique_ID_Set(uniqueID)
    !!{
    Resets the global unique ID number.
    !!}
    implicit none
    integer(kind=kind_int8), intent(in   ) :: uniqueID

    uniqueIdCount=uniqueID
    return
  end subroutine Galacticus_Nodes_Unique_ID_Set

  !
  ! Functions for treeNode class.
  function Tree_Node_Constructor(index,hostTree) result(self)
    !!{
    Return a pointer to a newly created and initialized {\normalfont \ttfamily treeNode}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (treeNode      ), pointer                         :: self
    integer(kind=kind_int8), intent(in   ), optional         :: index
    type   (mergerTree    ), intent(in   ), optional, target :: hostTree
    integer                                                  :: allocErr

    ! Allocate the object.
    allocate(self,stat=allocErr)
    if (allocErr/=0) call Error_Report('unable to allocate node'//{introspection:location})
    ! Initialize the node.
    call self%initialize(index,hostTree)
    return
  end function Tree_Node_Constructor

  function Tree_Node_Type(self)
    !!{
    Returns the name of a {\normalfont \ttfamily treeNode} object.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class(treeNode      ), intent(in   ) :: self
    type (varying_string)                :: Tree_Node_Type
    !$GLC attributes unused :: self

    Tree_Node_Type="treeNode"
    return
  end function Tree_Node_Type

  function Tree_Node_Index(self)
    !!{
    Returns the index of a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    class  (treeNode      ), pointer               :: workNode
    integer(kind=kind_int8)                        :: Tree_Node_Index

    workNode => self
    if (associated(workNode)) then
       Tree_Node_Index=workNode%indexValue
    else
       Tree_Node_Index=-1
    end if
    return
  end function Tree_Node_Index

  subroutine Tree_Node_Index_Set(self,index)
    !!{
    Sets the index of a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class  (treeNode      ), intent(inout) :: self
    integer(kind=kind_int8), intent(in   ) :: index

    self%indexValue=index
    return
  end subroutine Tree_Node_Index_Set

  function Tree_Node_Unique_ID(self) result(uniqueID)
    !!{
    Returns the unique ID of a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class  (treeNode      ), intent(in   ) :: self
    integer(kind=kind_int8)                        :: uniqueID

    uniqueID=self%uniqueIdValue
    return
  end function Tree_Node_Unique_ID

  subroutine Tree_Node_Unique_ID_Set(self,uniqueID)
    !!{
    Sets the index of a {\normalfont \ttfamily treeNode}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (treeNode      ), intent(inout)           :: self
    integer(kind=kind_int8), intent(in   ), optional :: uniqueID

    if (present(uniqueID)) then
       self%uniqueIdValue=uniqueID
    else
       !$omp critical(UniqueID_Assign)
       uniqueIDCount=uniqueIDCount+1
       if (uniqueIDCount <= 0) call Error_Report('ran out of unique ID numbers'//{introspection:location})
       self%uniqueIdValue=uniqueIDCount
       !$omp end critical(UniqueID_Assign)
    end if
    return
  end subroutine Tree_Node_Unique_ID_Set

  double precision function Tree_Node_Time_Step(self)
    !!{
    Returns the time-step last used by a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class(treeNode), intent(in   ) :: self

    Tree_Node_Time_Step=self%timeStepValue
    return
  end function Tree_Node_Time_Step

  subroutine Tree_Node_Time_Step_Set(self,timeStep)
    !!{
    Sets the time-step used by a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class           (treeNode      ), intent(inout) :: self
    double precision                , intent(in   ) :: timeStep

    self%timeStepValue=timeStep
    return
  end subroutine Tree_Node_Time_Step_Set

  double precision function Tree_Node_Subsampling_Weight(self)
    !!{
    Returns the subsampling weight of a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class(treeNode), intent(in   ) :: self

    Tree_Node_Subsampling_Weight=self%subsamplingWeightValue
    return
  end function Tree_Node_Subsampling_Weight

  subroutine Tree_Node_Subsampling_Weight_Set(self,subsamplingWeight)
    !!{
    Sets the time-step used by a {\normalfont \ttfamily treeNode}.
    !!}
    implicit none
    class           (treeNode      ), intent(inout) :: self
    double precision                , intent(in   ) :: subsamplingWeight

    self%subsamplingWeightValue=subsamplingWeight
    return
  end subroutine Tree_Node_Subsampling_Weight_Set

  subroutine Tree_Node_Attach_Event(self,newEvent)
    !!{
    Create a new event in a tree node.
    !!}
    implicit none
    class(treeNode ), intent(inout)          :: self
    class(nodeEvent), intent(inout), pointer :: newEvent
    class(nodeEvent)               , pointer :: event

    !$omp critical (nodeEventIncrement)
    eventID       =  eventID+1
    newEvent%ID   =  eventID
    !$omp end critical (nodeEventIncrement)
    newEvent%next => null()
    if (associated(self%event)) then
       event => self%event
       do while (associated(event%next))
          event => event%next
       end do
       event%next => newEvent
    else
       self%event => newEvent
    end if
    return
  end subroutine Tree_Node_Attach_Event

  subroutine Tree_Node_Remove_Paired_Event(self,event)
    !!{
    Removed a paired event from {\normalfont \ttfamily self}. Matching is done on the basis of event ID.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (treeNode ), intent(inout) :: self
    class  (nodeEvent), intent(in   ) :: event
    class  (nodeEvent), pointer       :: eventLast  , eventNext, eventPaired
    logical                           :: pairMatched

    ! Locate the paired event in self and remove it.
    eventPaired => self%event
    eventLast   => self%event
    ! Iterate over all events.
    pairMatched=.false.
    do while (associated(eventPaired).and..not.pairMatched)
       ! Match the paired event ID with the current event ID.
       if (eventPaired%ID == event%ID) then
          pairMatched=.true.
          if (associated(eventPaired,self%event)) then
             self     %event => eventPaired%next
             eventLast       => self       %event
          else
             eventLast%next  => eventPaired%next
          end if
          eventNext => eventPaired%next
          deallocate(eventPaired)
          eventPaired => eventNext
       else
          eventLast   => eventPaired
          eventPaired => eventPaired%next
       end if
    end do
    if (.not.pairMatched) call Error_Report('unable to find paired event'//{introspection:location})
    return
  end subroutine Tree_Node_Remove_Paired_Event

  logical function treeNodeIsPrimaryProgenitor(self) result(isPrimaryProgenitor)
    !!{
    Returns true if {\normalfont \ttfamily self} is the primary progenitor of its parent node.
    !!}
    implicit none
    class(treeNode), intent(inout), target:: self

    if (associated(self%parent)) then
       isPrimaryProgenitor=associated(self%parent%firstChild,self)
    else
       isPrimaryProgenitor=.false.
    end if
    return
  end function treeNodeIsPrimaryProgenitor

  logical function Tree_Node_Is_Primary_Progenitor_Of_Index(self,targetNodeIndex)
    !!{
    Return true if {\normalfont \ttfamily self} is a progenitor of the node with index {\normalfont \ttfamily targetNodeIndex}.
    !!}
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    integer(kind=kind_int8), intent(in   )         :: targetNodeIndex
    class  (treeNode      ), pointer               :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Index=.false.
    workNode => self
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Primary_Progenitor_Of_Index=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Index

  logical function Tree_Node_Is_Primary_Progenitor_Of_Node(self,targetNode)
    !!{
    Return true if {\normalfont \ttfamily self} is a progenitor of {\normalfont \ttfamily targetNode}.
    !!}
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode), intent(in   ), pointer :: targetNode
    class(treeNode)               , pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Node=.false.
    workNode => self
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Primary_Progenitor_Of_Node=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Node

  logical function Tree_Node_Is_Progenitor_Of_Index(self,targetNodeIndex)
    !!{
    Return true if {\normalfont \ttfamily self} is a progenitor of the node with index {\normalfont \ttfamily targetNodeIndex}.
    !!}
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    integer(kind=kind_int8), intent(in   )         :: targetNodeIndex
    class  (treeNode      ), pointer               :: workNode

    Tree_Node_Is_Progenitor_Of_Index=.false.
    workNode => self
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Progenitor_Of_Index=.true.
          return
       end if
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Index

  logical function Tree_Node_Is_Progenitor_Of_Node(self,targetNode)
    !!{
    Return true if {\normalfont \ttfamily self} is a progenitor of {\normalfont \ttfamily targetNode}.
    !!}
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode), intent(in   ), pointer :: targetNode
    class(treeNode)               , pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Node=.false.
    workNode => self
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Progenitor_Of_Node=.true.
          return
       end if
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Node

  logical function Tree_Node_Is_On_Main_Branch(self)
    !!{
    Returns true if {\normalfont \ttfamily self} is on the main branch.
    !!}
    implicit none
    class(treeNode), intent(inout), target :: self
    class(treeNode), pointer               :: workNode

    Tree_Node_Is_On_Main_Branch=.not.associated(self%parent)
    workNode => self
    do while (associated(workNode%parent))
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    Tree_Node_Is_On_Main_Branch=.true.
    return
  end function Tree_Node_Is_On_Main_Branch

  logical function Tree_Node_Is_Satellite(self)
    !!{
    Returns true if {\normalfont \ttfamily self} is a satellite.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(treeNode), intent(in   ), target :: self
    type (treeNode), pointer               :: childNode, parentNode, selfActual

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Error_Report('treeNode of unknown class - this should not happen'//{introspection:location})
    end select
    parentNode => selfActual%parent
    select case (associated(parentNode))
    case (.false.)
       Tree_Node_Is_Satellite=.false.
       return
    case (.true.)
       childNode => parentNode%firstChild
       Tree_Node_Is_Satellite=.true.
       do while (associated(childNode))
          if (associated(childNode,selfActual)) then
             Tree_Node_Is_Satellite=.false.
             exit
          end if
          childNode => childNode%sibling
       end do
    end select
    return
  end function Tree_Node_Is_Satellite

  function Tree_Node_Get_Isolated_Parent(self) result (isolatedNode)
    !!{
    Returns a pointer to the isolated parent node of {\normalfont \ttfamily self}.
    !!}
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode)               , pointer :: isolatedNode

    isolatedNode => self
    do while (isolatedNode%isSatellite())
       isolatedNode => isolatedNode%parent
    end do
    return
  end function Tree_Node_Get_Isolated_Parent

  function Tree_Node_Get_Last_Satellite(self) result (satelliteNode)
    !!{
    Returns a pointer to the final satellite node associated with {\normalfont \ttfamily self}.
    !!}
    implicit none
    class(treeNode), intent(in   ) :: self
    type (treeNode), pointer       :: satelliteNode

    satelliteNode => self%firstSatellite
    do while (associated(satelliteNode%sibling))
       satelliteNode => satelliteNode%sibling
    end do
    return
  end function Tree_Node_Get_Last_Satellite

  function Tree_Node_Get_Earliest_Progenitor(self) result (progenitorNode)
    !!{
    Returns a pointer to the earliest progenitor of {\normalfont \ttfamily self}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (treeNode), pointer       :: progenitorNode
    class(treeNode), intent(inout) :: self

    select type (self)
    type is (treeNode)
       progenitorNode => self
       do while (associated(progenitorNode%firstChild))
          progenitorNode => progenitorNode%firstChild
       end do
    class default
       progenitorNode => null()
       call Error_Report('treeNode of unknown class - this should not happen'//{introspection:location})
    end select
    return
  end function Tree_Node_Get_Earliest_Progenitor

  function Tree_Node_Merges_With_Node(node)
    !!{
    Returns a pointer to the node with which {\normalfont \ttfamily node} will merge.
    !!}
    implicit none
    class(treeNode), intent(in   ) :: node
    type (treeNode), pointer       :: Tree_Node_Merges_With_Node

    ! Check if a specific merge node has been set.
    if (associated(node%mergeTarget)) then
       ! One has, so simply return it.
       Tree_Node_Merges_With_Node => node%mergeTarget
    else
       ! No specific merge node has been set, assume merging with the parent node.
       Tree_Node_Merges_With_Node => node%parent
    end if
    return
  end function Tree_Node_Merges_With_Node

  subroutine Tree_Node_Remove_From_Host(self)
    !!{
    Remove {\normalfont \ttfamily self} from the linked list of its host node's satellites.
    !!}
    use :: Display           , only : displayMessage, displayVerbosity, verbosityLevelInfo
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=) , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: nodeHost, nodePrevious, selfActual, node
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Error_Report('treeNode of unknown class'//{introspection:location})
    end select

    ! Remove from the parent node satellite list.
    nodeHost => selfActual%parent
    if (displayVerbosity() >= verbosityLevelInfo) then
       message='Satellite node ['
       message=message//selfActual%index()//'] being removed from host node ['//nodeHost%index()//']'
       call displayMessage(message,verbosityLevelInfo)
    end if
    if (associated(nodeHost%firstSatellite,selfActual)) then
       ! This is the first satellite, unlink it, and link to any sibling.
       nodeHost%firstSatellite => selfActual%sibling
    else
       node         => nodeHost%firstSatellite
       nodePrevious => null()
       do while (associated(node))
          if (associated(node,selfActual)) then
             ! Found our node, link its older sibling to its younger sibling.
             nodePrevious%sibling => node%sibling
             exit
          end if
          nodePrevious => node
          node     => node%sibling
       end do
    end if
    return
  end subroutine Tree_Node_Remove_From_Host

  subroutine Tree_Node_Remove_from_Mergee(self)
    !!{
    Remove {\normalfont \ttfamily self} from the linked list of its host node's satellites.
    !!}
    use :: Display           , only : displayMessage, verbosityLevelInfo
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=) , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: nodeHost, nodePrevious, selfActual, node
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Error_Report('treeNode of unknown class'//{introspection:location})
    end select

    ! Remove from the mergee list of any merge target.
    if (associated(selfActual%mergeTarget)) then
       nodeHost => selfActual%mergeTarget
       message='Mergee node ['
       message=message//selfActual%index()//'] being removed from merge target ['//nodeHost%index()//']'
       call displayMessage(message,verbosityLevelInfo)
       if (associated(nodeHost%firstMergee,selfActual)) then
          ! This is the first mergee, unlink it, and link to any sibling.
          nodeHost%firstMergee => selfActual%siblingMergee
       else
          node         => nodeHost%firstMergee
          nodePrevious => null()
          do while (associated(node))
             if (associated(node,selfActual)) then
                ! Found our node, link its older sibling to its younger sibling.
                nodePrevious%siblingMergee => node%siblingMergee
                exit
             end if
             nodePrevious => node
             node         => node%siblingMergee
          end do
       end if
    end if
    return
  end subroutine Tree_Node_Remove_from_Mergee

  function treeNodeWalkTreeWithSatellites(self)
    !!{
    Merger tree walk function which also descends through satellite nodes. Note that it is
    important that the walk descends to satellites before descending to children: the
    routines that destroy merger tree branches rely on this since child nodes are used in
    testing whether a node is a satellite---if they are destroyed prior to the test being
    made then problems with dangling pointers will occur.
    !!}
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkTreeWithSatellites
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: workNode

    workNode => self
    if (.not.associated(workNode%parent)) then
       ! This is the base of the merger tree.
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,self)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(workNode%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parent%firstChild
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no children, so move to the parent.
                workNode => workNode%parent
             end if
          else
             ! It is not a satellite, so all satellites and children have been processed.
             workNode => workNode%parent
          end if
          ! Terminate when back at tree base.
          if (.not.associated(workNode%parent)) workNode => null()
       end if
    end if
    treeNodeWalkTreeWithSatellites => workNode
    return
  end function treeNodeWalkTreeWithSatellites

  function treeNodeWalkBranchWithSatellites(self,startNode)
    !!{
    This function provides a mechanism for walking through the branches of the merger
    tree. Given a pointer {\normalfont \ttfamily self} to a branch of the tree, it will
    return the next node that should be visited in the tree. Thus, if {\normalfont \ttfamily
    self} is initially set to the base of the merger tree and {\normalfont \ttfamily
    Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every node of the
    branch. Once the entire branch has been walked, a {\normalfont \ttfamily null()} pointer
    will be returned, indicating that there are no more nodes to walk. Each node will be
    visited once and once only if the branch is walked in this way. Note that it is important
    that the walk descends to satellites before descending to children: the routines that
    destroy merger tree branches rely on this since child nodes are used in testing whether a
    node is a satellite---if they are destroyed prior to the test being made then problems
    with dangling pointers will occur.
    !!}
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkBranchWithSatellites
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: startNode
    type (treeNode)               , pointer :: workNode                        , selfNode

    selfNode => self
    workNode => self
    if (associated(selfNode,startNode)) then
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,selfNode)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(workNode%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parent%firstChild
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no satellites, so move to the parent.
                workNode => workNode%parent
             end if
          else
             ! It is not a satellite, so all satellites and children of the parent must have
             ! been processed. Therefore, move to the parent.
             workNode => workNode%parent
          end if
          ! Terminate when back at starting node.
          if (associated(workNode,startNode)) workNode => null()
       end if
    end if
    treeNodeWalkBranchWithSatellites => workNode
    return
  end function treeNodeWalkBranchWithSatellites

  function Merger_Tree_Walk_Descend_to_Progenitors(self) result (progenitorNode)
    !!{
    Descend to the deepest progenitor (satellites and children) of {\normalfont \ttfamily self}.
    !!}
    implicit none
    type(treeNode), intent(in   ), target  :: self
    type(treeNode)               , pointer :: progenitorNode

    ! Begin at the input node.
    progenitorNode => self

    ! Descend through satellites and children.
    do while (associated(progenitorNode%firstSatellite).or.associated(progenitorNode%firstChild))
       if (associated(progenitorNode%firstSatellite)) then
          progenitorNode => progenitorNode%firstSatellite
       else
          progenitorNode => progenitorNode%firstChild
       end if
  end do
    return
  end function Merger_Tree_Walk_Descend_to_Progenitors

  subroutine treeNodeDestroyBranch(self)
    !!{
    Destroy the tree branch rooted at this given node.
    !!}
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: nodeDestroy, nodeNext, &
         &                                     branchTip

    ! Descend to the tip of the branch.
    branchTip => self
    nodeNext  => branchTip%walkBranchWithSatellites(branchTip)
    ! Loop over all tree nodes.
    do while (associated(nodeNext))
       ! Keep of a record of the current node, so that we can destroy it.
       nodeDestroy => nodeNext
       ! Walk to the next node in the tree.
       nodeNext => nodeDestroy%walkBranchWithSatellites(branchTip)
       ! If the node about to be destroyed is the primary progenitor of its parent we must move the child pointer of the parent to
       ! point to the node's sibling. This is necessary as parent-child pointers are used to establish satellite status and so
       ! will be utilized when walking the tree. Failure to do this can result in attempts to use dangling pointers.
       if (associated(nodeDestroy%parent).and.associated(nodeDestroy%parent%firstChild,nodeDestroy)) &
            & nodeDestroy%parent%firstChild => nodeDestroy%sibling
       ! Destroy the current node.
       call nodeDestroy%destroy()
       deallocate(nodeDestroy)
    end do
    ! Destroy the base node of the branch.
    if (associated(self%parent)) then
       if (associated(self%parent%firstChild,self)) self%parent%firstChild => self%sibling
    end if
    call self%destroy()
    return
  end subroutine treeNodeDestroyBranch

  !
  ! Functions for nodeComponent class.
  function Node_Component_Generic_Type(self)
    !!{
    Returns the name of a generic tree node component.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class(nodeComponent ), intent(in   ) :: self
    type (varying_string)                :: Node_Component_Generic_Type
    !$GLC attributes unused :: self

    Node_Component_Generic_Type="nodeComponent"
    return
  end function Node_Component_Generic_Type

  subroutine Node_Component_Generic_Destroy(self)
    !!{
    Destroy a generic tree node component.
    !!}
    implicit none
    class(nodeComponent), intent(inout) :: self
    !$GLC attributes unused :: self

    ! Do nothing.
    return
  end subroutine Node_Component_Generic_Destroy

  integer function Node_Component_Generic_Add_Float_Rank0_Meta_Property(self,label,name,isEvolvable,isCreator)
    !!{
    Add a meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isEvolvable, isCreator
    !$GLC attributes unused :: self, label, name, isEvolvable, isCreator

    Node_Component_Generic_Add_Float_Rank0_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_Float_Rank0_Meta_Property

  integer function Node_Component_Generic_Add_Float_Rank1_Meta_Property(self,label,name,isCreator)
    !!{
    Add a rank-1 meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isCreator
    !$GLC attributes unused :: self, label, name, isCreator

    Node_Component_Generic_Add_Float_Rank1_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_Float_Rank1_Meta_Property

  integer function Node_Component_Generic_Add_LongInteger_Rank0_Meta_Property(self,label,name,isCreator)
    !!{
    Add a long integer meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isCreator
    !$GLC attributes unused :: self, label, name, isCreator

    Node_Component_Generic_Add_LongInteger_Rank0_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_LongInteger_Rank0_Meta_Property
  
  integer function Node_Component_Generic_Add_Integer_Rank0_Meta_Property(self,label,name,isCreator)
    !!{
    Add a long integer meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isCreator
    !$GLC attributes unused :: self, label, name, isCreator

    Node_Component_Generic_Add_Integer_Rank0_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_Integer_Rank0_Meta_Property
  
  integer function Node_Component_Generic_Add_Integer_Rank1_Meta_Property(self,label,name,isCreator)
    !!{
    Add an integer meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isCreator
    !$GLC attributes unused :: self, label, name, isCreator

    Node_Component_Generic_Add_Integer_Rank1_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_Integer_Rank1_Meta_Property
  
  integer function Node_Component_Generic_Add_LongInteger_Rank1_Meta_Property(self,label,name,isCreator)
    !!{
    Add a longinteger meta-property to a node component.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (nodeComponent ), intent(inout)           :: self
    type     (varying_string), intent(in   )           :: label
    character(len=*         ), intent(in   )           :: name
    logical                  , intent(in   ), optional :: isCreator
    !$GLC attributes unused :: self, label, name, isCreator

    Node_Component_Generic_Add_LongInteger_Rank1_Meta_Property=-1
    call Error_Report('can not add meta-properties to a generic nodeComponent'//{introspection:location})
    return
  end function Node_Component_Generic_Add_LongInteger_Rank1_Meta_Property
  
  !![
  <metaPropertyDatabase/>
  !!]
  
  subroutine Node_Component_ODE_Step_Initialize_Null(self)
    !!{
    Initialize a generic tree node component for an ODE solver step.
    !!}
    implicit none
    class(nodeComponent), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine Node_Component_ODE_Step_Initialize_Null

  subroutine Node_Component_Dump_Null(self,verbosityLevel)
    !!{
    Dump a generic tree node component.
    !!}
    use :: Display, only : enumerationVerbosityLevelType
    implicit none
    class(nodeComponent                ), intent(in   ) :: self
    type (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel 
    !$GLC attributes unused :: self, verbosityLevel

    return
  end subroutine Node_Component_Dump_Null

  subroutine Node_Component_Dump_XML_Null(self,fileHandle)
    !!{
    Dump a generic tree node component to XML.
    !!}
    implicit none
    class  (nodeComponent), intent(inout) :: self
    integer               , intent(in   ) :: fileHandle
    !$GLC attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Dump_XML_Null

  subroutine Node_Component_Dump_Raw_Null(self,fileHandle)
    !!{
    Dump a generic tree node component in binary.
    !!}
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(in   ) :: fileHandle
    !$GLC attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Dump_Raw_Null

  subroutine Node_Component_Read_Raw_Null(self,fileHandle)
    !!{
    Read a generic tree node component in binary.
    !!}
    implicit none
    class  (nodeComponent), intent(inout) :: self
    integer               , intent(in   ) :: fileHandle
    !$GLC attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Read_Raw_Null

  subroutine Node_Component_Output_Count_Null(self,integerPropertyCount,doublePropertyCount,time,instance)
    !!{
    Dump a generic tree node component.
    !!}
    implicit none
    class           (nodeComponent), intent(inout) :: self
    integer                        , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision               , intent(in   ) :: time
    integer                        , intent(in   ) :: instance
    !$GLC attributes unused :: self, integerPropertyCount, doublePropertyCount, time, instance

    return
  end subroutine Node_Component_Output_Count_Null

  subroutine Node_Component_Output_Names_Null(self,integerProperty,integerProperties,doubleProperty,doubleProperties,time,instance)
    !!{
    Dump a generic tree node component.
    !!}
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    class           (nodeComponent        )              , intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty   , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    integer                                              , intent(in   ) :: instance
    !$GLC attributes unused :: self, integerProperty, integerProperties, doubleProperty, doubleProperties, time, instance

    return
  end subroutine Node_Component_Output_Names_Null

  subroutine Node_Component_Output_Null(self,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,outputInstance,instance)
    !!{
    Dump a generic tree node component.
    !!}
    use :: Multi_Counters, only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    class           (nodeComponent        )              , intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleBufferCount , doubleProperty ,  &
         &                                                                  integerBufferCount, integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    type            (multiCounter         )              , intent(in   ) :: outputInstance
    integer                                              , intent(in   ) :: instance
    !$GLC attributes unused :: self, integerProperty, integerBufferCount, integerProperties, doubleProperty, doubleBufferCount, doubleProperties, time, outputInstance, instance

    return
  end subroutine Node_Component_Output_Null

  integer function Node_Component_Serialize_Count_Zero(self,propertyType)
    !!{
    Return the serialization count of a generic tree node component.
    !!}
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(in   ) :: propertyType
    !$GLC attributes unused :: self, propertyType

    Node_Component_Serialize_Count_Zero=0
    return
  end function Node_Component_Serialize_Count_Zero

  subroutine Node_Component_Serialization_Offsets(self,count,countSubset,propertyType)
    !!{
    Return the serialization count of a generic tree node component.
    !!}
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(inout) :: count       , countSubset
    integer               , intent(in   ) :: propertyType
    !$GLC attributes unused :: self, count, countSubset, propertyType

    return
  end subroutine Node_Component_Serialization_Offsets

  subroutine Node_Component_Serialize_Null(self,array,propertyType)
    !!{
    Serialize a generic tree node component.
    !!}
    implicit none
    class           (nodeComponent)              , intent(in   ) :: self
    double precision               , dimension(:), intent(  out) :: array
    integer                                      , intent(in   ) :: propertyType
    !$GLC attributes unused :: self, array, propertyType

    return
  end subroutine Node_Component_Serialize_Null

  subroutine Node_Component_Serialize_NonNegative_Null(self,array)
    !!{
    Serialize the non-negative status for a generic tree node component.
    !!}
    implicit none
    class  (nodeComponent)              , intent(in   ) :: self
    logical               , dimension(:), intent(  out) :: array
    !$GLC attributes unused :: self, array

    return
  end subroutine Node_Component_Serialize_NonNegative_Null

  subroutine Node_Component_Deserialize_Null(self,array,propertyType)
    !!{
    Deserialize a generic tree node component.
    !!}
    implicit none
    class           (nodeComponent)              , intent(inout) :: self
    double precision               , dimension(:), intent(in   ) :: array
    integer                                      , intent(in   ) :: propertyType
    !$GLC attributes unused :: self, array, propertyType

    return
  end subroutine Node_Component_Deserialize_Null

  function Node_Component_Host_Node(self)
    !!{
    Return the host tree node of a tree node component.
    !!}
    implicit none
    class(nodeComponent), intent(in   ) :: self
    type (treeNode     ), pointer       :: Node_Component_Host_Node

    Node_Component_Host_Node => self%hostNode
    return
  end function Node_Component_Host_Node

  subroutine Node_Component_Null_Void0_InOut(self)
    !!{
    A null {\normalfont \ttfamily void} function for rank 0 {\normalfont \ttfamily nodeComponent} arrays.
    !!}
    implicit none
    class(nodeComponent), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine Node_Component_Null_Void0_InOut

  double precision function Node_Component_Null_Double0_InOut(self)
    !!{
    A null {\normalfont \ttfamily double} function for rank 0 {\normalfont \ttfamily nodeComponent} arrays.
    !!}
    implicit none
    class(nodeComponent), intent(inout) :: self
    !$GLC attributes unused :: self

    Node_Component_Null_Double0_InOut=0.0d0
    return
  end function Node_Component_Null_Double0_InOut

  function Node_Component_Null_Double1_InOut(self,resultSize)
    !!{
    A null {\normalfont \ttfamily double} function for rank 1 {\normalfont \ttfamily nodeComponent} arrays.
    !!}
    implicit none
    class           (nodeComponent), intent(inout)         :: self
    integer                        , intent(in   )         :: resultSize
    double precision               , dimension(resultSize) :: Node_Component_Null_Double1_InOut
    !$GLC attributes unused :: self

    Node_Component_Null_Double1_InOut=0.0d0
    return
  end function Node_Component_Null_Double1_InOut

  function Node_Component_Null_TensorR2D3_InOut(self)
    !!{
    A null {\normalfont \ttfamily tensorRank2Dimension3Symmetric} function for {\normalfont \ttfamily nodeComponent}s.
    !!}
    implicit none
    type (tensorRank2Dimension3Symmetric)                :: Node_Component_Null_TensorR2D3_InOut
    class(nodeComponent                 ), intent(inout) :: self
    !$GLC attributes unused :: self

    Node_Component_Null_TensorR2D3_InOut=tensorNullR2D3Sym
    return
  end function Node_Component_Null_TensorR2D3_InOut

  function Node_Component_Mass_Distribution_Null(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
    !!{
    A null implementation of the mass distribution factory for a component. Always returns null.
    !!}
    use :: Galactic_Structure_Options, only : enumerationWeightByType, enumerationComponentTypeType, enumerationMassTypeType
    implicit none
    class  (massDistributionClass       ), pointer                 :: massDistribution_
    class  (nodeComponent               ), intent(inout)           :: self
    type   (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                              , intent(in   ), optional :: weightIndex
    !$GLC attributes unused :: self, componentType, massType, weightBy, weightIndex

    massDistribution_ => null()
    return
  end function Node_Component_Mass_Distribution_Null

  double precision function Node_Component_Mass_Baryonic_Null(self)
    !!{
    A null implementation of the total baryonic mass distribution. Always returns zero.
    !!}
    implicit none
    class(nodeComponent), intent(inout):: self
    !$GLC attributes unused :: self

    Node_Component_Mass_Baryonic_Null=0.0d0
    return
  end function Node_Component_Mass_Baryonic_Null

  double precision function Node_Component_Density_Null(self,positionSpherical,componentType,massType,weightBy,weightIndex)
    !!{
    A null implementation of the density in a component. Always returns zero.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
    implicit none
    class           (nodeComponent               )              , intent(inout) :: self
    type            (enumerationComponentTypeType)              , intent(in   ) :: componentType
    type            (enumerationMassTypeType     )              , intent(in   ) :: massType
    type            (enumerationWeightByType     )              , intent(in   ) :: weightBy
    integer                                                     , intent(in   ) :: weightIndex
    double precision                              , dimension(3), intent(in   ) :: positionSpherical
    !$GLC attributes unused :: self, positionSpherical, componentType, massType, weightBy, weightIndex

    Node_Component_Density_Null=0.0d0
    return
  end function Node_Component_Density_Null

  double precision function Node_Component_Density_Spherical_Average_Null(self,radius,componentType,massType,weightBy,weightIndex)
    !!{
    A null implementation of the spherically-averaged density in a component. Always returns zero.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
    implicit none
    class           (nodeComponent               ), intent(inout) :: self
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    type            (enumerationWeightByType     ), intent(in   ) :: weightBy
    integer                                       , intent(in   ) :: weightIndex
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius, componentType, massType, weightBy, weightIndex

    Node_Component_Density_Spherical_Average_Null=0.0d0
    return
  end function Node_Component_Density_Spherical_Average_Null

  double precision function Node_Component_Surface_Density_Null(self,positionCylindrical,componentType,massType,weightBy,weightIndex)
    !!{
    A null implementation of the surface density in a component. Always returns zero.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
    implicit none
    class           (nodeComponent               )              , intent(inout) :: self
    type            (enumerationComponentTypeType)              , intent(in   ) :: componentType
    type            (enumerationMassTypeType     )              , intent(in   ) :: massType
    type            (enumerationWeightByType     )              , intent(in   ) :: weightBy
    integer                                                     , intent(in   ) :: weightIndex
    double precision                              , dimension(3), intent(in   ) :: positionCylindrical
    !$GLC attributes unused :: self, positionCylindrical, componentType, massType, weightBy, weightIndex

    Node_Component_Surface_Density_Null=0.0d0
    return
  end function Node_Component_Surface_Density_Null

  ! Simple Boolean functions.
  logical function Boolean_False()
    !!{
    Returns Boolean false always.
    !!}
    implicit none

    Boolean_False=.false.
    return
  end function Boolean_False

  logical function Boolean_True()
    !!{
    Returns Boolean true always.
    !!}
    implicit none

    Boolean_True=.true.
    return
  end function Boolean_True

  subroutine nodeComponentGetError(nameComponent,indexNode)
    !!{
    Function used to report errors when attempting to get a component from a node. Error reporting is handled here to avoid
    having the relatively expensive creation/destruction of a varying string object in the actual get functions (which are
    called a very large number of times).
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=), operator(//), varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    character(len=*         ), intent(in   ) :: nameComponent
    integer  (kind_int8     ), intent(in   ) :: indexNode
    type     (varying_string)                :: message

    message='"'//nameComponent//'" component is not allocated in node '
    message=message//indexNode
    call Error_Report(message//{introspection:location})
    return
  end subroutine nodeComponentGetError
  
  function universeConstructor() result(self)
    !!{
    Constructor for the universe class.
    !!}
    implicit none
    type(universe) :: self

    !$omp critical(universeUniqueIDAssign)
    universeUniqueIdCount=universeUniqueIdCount+1
    self%uniqueID        =universeUniqueIdCount
    !$omp end critical(universeUniqueIDAssign)
    self%trees         => null   ()
    self%event         => null   ()
    self%allTreesBuilt =  .false.
    self%lock          =  ompLock()
    call self%attributes%initialize()
    return
  end function universeConstructor

  function mergerTreeConstructor() result(self)
    !!{
    Constructor for the merger tree class. Currently does nothing.
    !!}
    implicit none
    type(mergerTree) :: self
    
    return
  end function mergerTreeConstructor

  subroutine Merger_Tree_Destroy(tree)
    !!{
    Destroys the entire merger tree.
    !!}
    implicit none
    class(mergerTree), intent(inout) :: tree
    type (treeEvent ), pointer       :: event, eventNext

    select type (tree)
    type is (mergerTree)
       ! Destroy all nodes.
       if (associated(tree%nodeBase)) then
          call tree%nodeBase%destroyBranch()
          deallocate(tree%nodeBase)
       end if
       ! Destroy any events attached to this tree.
       event => tree%event
       do while (associated(event))
          eventNext => event%next
          deallocate(event)
          event => eventNext
       end do
       ! Destroy the property hash.
       call tree%properties%destroy()
       ! Destroy the random number generator.
       !![
       <objectDestructor name="tree%randomNumberGenerator_"/>
       !!]
    end select
    return
  end subroutine Merger_Tree_Destroy

  function Merger_Tree_Node_Get(tree,nodeIndex)
    !!{
    Return a pointer to a node in {\normalfont \ttfamily tree} given the index of the node.
    !!}
    implicit none
    class  (mergerTree    ), intent(in   ), target :: tree
    integer(kind=kind_int8), intent(in   )         :: nodeIndex
    type   (mergerTree    ), pointer               :: treeCurrent
    type   (treeNode      ), pointer               :: Merger_Tree_Node_Get, node

    Merger_Tree_Node_Get => null()
    treeCurrent => tree
    do while (associated(treeCurrent))
       node => treeCurrent%nodeBase
       do while (associated(node))
          if (node%index() == nodeIndex) then
             Merger_Tree_Node_Get => node
             return
          end if
          node => node%walkTreeWithSatellites()
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end function Merger_Tree_Node_Get

  function Merger_Tree_Create_Event(self) result (eventNew)
    !!{
    Create a new event in a merger tree.
    !!}
    implicit none
    class(mergerTree), intent(inout) :: self
    type (treeEvent ), pointer       :: eventNew, event

    allocate(eventNew)
    nullify(eventNew%next)
    !$omp critical (nodeEventIncrement)
    eventID=eventID+1
    eventNew%ID=eventID
    !$omp end critical (nodeEventIncrement)
    if (associated(self%event)) then
       event => self%event
       do while (associated(event%next))
          event => event%next
       end do
       event%next => eventNew
    else
       self%event => eventNew
    end if
    return
  end function Merger_Tree_Create_Event

  subroutine Merger_Tree_Remove_Event(self,event)
    !!{
    Removed an event from {\normalfont \ttfamily self}.
    !!}
    implicit none
    class  (mergerTree), intent(inout) :: self
    type   (treeEvent ), intent(in   ) :: event
    type   (treeEvent ), pointer       :: eventLast, eventNext, eventCurrent

    ! Destroy the event.
    eventLast    => null()
    eventCurrent => self%event
    do while (associated(eventCurrent))
       ! Match the event.
       if (eventCurrent%ID == event%ID) then
          if (associated(eventCurrent,self%event)) then
             self     %event => eventCurrent%next
             eventLast       => self %event
          else
             eventLast%next  => eventCurrent%next
          end if
          eventNext => eventCurrent%next
          deallocate(eventCurrent)
          eventCurrent => eventNext
       else
          eventLast    => eventCurrent
          eventCurrent => eventCurrent%next
       end if
    end do
    return
  end subroutine Merger_Tree_Remove_Event

  double precision function Merger_Tree_Earliest_Time(self)
    !!{
    Return the earliest time in a merger tree.
    !!}
    implicit none
    class           (mergerTree        ), intent(inout), target :: self
    double precision                    , parameter             :: timeInfinity=huge(1.0d0)
    type            (mergerTree        ), pointer               :: treeCurrent
    type            (treeNode          ), pointer               :: node
    class           (nodeComponentBasic), pointer               :: basic

    Merger_Tree_Earliest_Time =  timeInfinity
    treeCurrent               => self
    do while (associated(treeCurrent))
       node => treeCurrent%nodeBase
       do while (associated(node))
          if (.not.associated(node%firstChild)) then
             basic                     =>                               node %basic()
             Merger_Tree_Earliest_Time =  min(Merger_Tree_Earliest_Time,basic%time ())
          end if
          node => node%walkTreeWithSatellites()
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end function Merger_Tree_Earliest_Time

  double precision function Merger_Tree_Earliest_Time_Evolving(self)
    !!{
    Return the earliest time in a merger tree.
    !!}
    implicit none
    class           (mergerTree        ), intent(inout), target :: self
    double precision                    , parameter             :: timeInfinity=huge(1.0d0)
    type            (mergerTree        ), pointer               :: treeCurrent
    type            (treeNode          ), pointer               :: node
    class           (nodeComponentBasic), pointer               :: basic

    Merger_Tree_Earliest_Time_Evolving =  timeInfinity
    treeCurrent                        => self
    do while (associated(treeCurrent))
       node => treeCurrent%nodeBase
       do while (associated(node))
          if (.not.associated(node%firstChild).and.(associated(node%parent).or..not.associated(node%firstSatellite))) then
             basic                              =>                                        node %basic()
             Merger_Tree_Earliest_Time_Evolving =  min(Merger_Tree_Earliest_Time_Evolving,basic%time ())
          end if
          node => node%walkTreeWithSatellites()
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end function Merger_Tree_Earliest_Time_Evolving

  double precision function Merger_Tree_Latest_Time(self)
    !!{
    Return the latest time in a merger tree.
    !!}
    implicit none
    class           (mergerTree        ), intent(inout), target :: self
    type            (mergerTree        ), pointer               :: treeCurrent
    class           (nodeComponentBasic), pointer               :: basicBase

    Merger_Tree_Latest_Time =  -1.0d0
    treeCurrent             => self
    do while (associated(treeCurrent))
       if (associated(treeCurrent%nodeBase)) then
          basicBase               =>                             treeCurrent%nodeBase%basic()
          Merger_Tree_Latest_Time =  max(Merger_Tree_Latest_Time,basicBase  %time          ())
       end if
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end function Merger_Tree_Latest_Time

  integer(c_size_t) function Merger_Tree_Size_Of(self)
    !!{
    Return the size (in bytes) of a merger tree.
    !!}
    implicit none
    class(mergerTree), intent(in   ) :: self
    type (treeEvent ), pointer       :: event
    
    !![
    <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
     <description>
      Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
     </description>
    </workaround>
    !!]
    !Merger_Tree_Size_Of=sizeof(self)
    Merger_Tree_Size_Of=0_c_size_t
    event => self%event
    do while (associated(event))
       !![
       <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
        <description>
         Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
        </description>
       </workaround>
       !!]
       ! Merger_Tree_Size_Of=Merger_Tree_Size_Of+sizeof(event)
       event => event%next
    end do
    return
  end function Merger_Tree_Size_Of


  function universeCreateEvent(self) result (newEvent)
    !!{
    Create a new event in a universe.
    !!}
    implicit none
    class(universe     ), intent(inout) :: self
    type (universeEvent), pointer       :: newEvent, event

    allocate(newEvent)
    nullify(newEvent%next)
    !$omp critical (nodeEventIncrement)
    eventID=eventID+1
    newEvent%ID=eventID
    !$omp end critical (nodeEventIncrement)
    if (associated(self%event)) then
       event => self%event
       do while (associated(event%next))
          event => event%next
       end do
       event%next => newEvent
    else
       self%event => newEvent
    end if
    return
  end function universeCreateEvent

  subroutine universeRemoveEvent(self,event)
    !!{
    Remove an event from {\normalfont \ttfamily self}.
    !!}
    implicit none
    class  (universe     ), intent(inout) :: self
    type   (universeEvent), intent(in   ) :: event
    type   (universeEvent), pointer       :: eventLast, eventNext, event_
    integer(kind_int8    )                :: eventID
    
    ! Destroy the event.
    eventID   =  event%ID
    eventLast => null()
    event_    => self%event
    do while (associated(event_))
       ! Match the event.
       if (event_%ID == eventID) then
          if (associated(event_,self%event)) then
             self     %event => event_%next
             eventLast       => self  %event
          else
             eventLast%next  => event_%next
          end if
          eventNext => event_   %next
          deallocate(event_)
          event_    => eventNext
       else
          eventLast => event_
          event_    => event_%next
       end if
    end do
    return
  end subroutine universeRemoveEvent
  
  function universePopTree(self)
    !!{
    Pop a tree from the universe.
    !!}
    implicit none
    type (mergerTree    ), pointer       :: universePopTree
    class(universe      ), intent(inout) :: self
    type (mergerTreeList), pointer       :: nextTree

    if (associated(self%trees)) then
       universePopTree => self%trees%tree
       nextTree        => self%trees%next
       deallocate(self%trees)
       self%trees      => nextTree
    else
       universePopTree => null()
    end if
    return
  end function universePopTree

  subroutine universePushTree(self,tree)
    !!{
    Pop a tree from the universe.
    !!}
    implicit none
    class(universe      ), intent(inout)          :: self
    type (mergerTree    ), intent(in   ), pointer :: tree
    type (mergerTreeList)               , pointer :: treeNext, treeNew

    allocate(treeNew)
    treeNew%tree => tree
    treeNew%next => null()
    if (associated(self%trees)) then
       treeNext => self%trees
       do while (associated(treeNext%next))
          treeNext => treeNext%next
       end do
       treeNext%next => treeNew
    else
       self%trees => treeNew
    end if
    return
  end subroutine universePushTree

  logical function propertyActive(propertyType)
    !!{
    Returns true if active property evaluate is underway.
    !!}
    implicit none
    integer, intent(in   ) :: propertyType

    propertyActive=propertyType == propertyTypeActive
    return
  end function propertyActive
  
  logical function propertyInactive(propertyType)
    !!{
    Returns true if inactive property evaluate is underway.
    !!}
    implicit none
    integer, intent(in   ) :: propertyType

    propertyInactive=propertyType == propertyTypeInactive
    return
  end function propertyInactive
  
  logical function propertyEvaluate(propertyType,propertyIsInactive)
    !!{
    Returns true if a property should be evaluated during the current stage of evolution.
    !!}
    implicit none
    integer, intent(in   ) :: propertyType
    logical, intent(in   ) :: propertyIsInactive

    propertyEvaluate= (propertyType == propertyTypeAll                                   ) &
         &           .or.                                                                  &
         &            (propertyType == propertyTypeActive   .and. .not.propertyIsInactive) &
         &           .or.                                                                  &
         &             propertyType == propertyTypeNumerics                                &
         &           .or.                                                                  &
         &            (propertyType == propertyTypeInactive .and.      propertyIsInactive)
    return
  end function propertyEvaluate

  subroutine massDistributionCalculationReset(massDistributionsLast,node,uniqueID)
    !!{
    Reset the memoized {\normalfont \ttfamily massDistribution} due to a {\normalfont \ttfamily calculationReset} event.
    !!}
    implicit none
    integer           , intent(inout) :: massDistributionsLast
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    integer                           :: i
    !$GLC attributes unused :: massDistributionsLast, node, uniqueID
    
    do i=1,massDistributionsCount
       !![
       <objectDestructor name="massDistributions__(i)%massDistribution_"/>
       !!]
       massDistributions__(i)%uniqueID=-huge(kind_int8)
    end do
    return
  end subroutine massDistributionCalculationReset
  
  subroutine massDistributionsDestroy()
    !!{
    Destroy memoized {\normalfont \ttfamily massDistributions}.
    !!}
    implicit none    
    integer :: i

    do i=1,massDistributionsCount
       !![
       <objectDestructor name="massDistributions__(i)%massDistribution_"/>
       !!]
    end do
    return
  end subroutine massDistributionsDestroy
  
end module Galacticus_Nodes
