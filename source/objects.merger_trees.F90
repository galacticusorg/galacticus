!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module defining the merger tree object type.

module Merger_Trees
  !% Defines the merger tree object type.
  use, intrinsic :: ISO_C_Binding
  use               Tree_Nodes
  use               Components
  use               FGSL
  use               Events_Interrupts
  use               IO_HDF5
  use               ISO_Varying_String
  use               Kind_Numbers
  !# <include directive="treeNodeCreateInitialize" type="moduleUse">
  include 'objects.tree_node.create.modules.inc'
  !# </include>
  implicit none
  private
  public :: mergerTree, Tree_Node_Is_Accurate, Tree_Node_Get

  type mergerTree
     !% The merger tree object type.
     integer(kind=kind_int8)   :: index
     type(hdf5Object)          :: hdf5Group
     double precision          :: volumeWeight
     logical                   :: initialized
     type(treeNode),   pointer :: baseNode => null()
     type(mergerTree), pointer :: nextTree => null()
   contains
     ! Tree creation/destruction.
     !@ <objectMethods>
     !@   <object>mergerTree</object>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroys the merger tree, including all nodes and their components.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroyBranch</method>
     !@     <description>
     !@       Destroys a branch of a merger tree starting from the supplied node. All nodes and their components are destroyed.
     !@     </description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: destroy       => Merger_Tree_Destroy
     procedure                 :: destroyBranch => Merger_Tree_Destroy_Branch
     ! Node creation.
     !@ <objectMethod>
     !@   <object>mergerTree</object>
     !@   <method>createNode(thisNode)</method>
     !@   <description>Creates {\tt thisNode}.</description>
     !@ </objectMethod>
     procedure                 :: createNode    => Tree_Node_Create
     ! Node promotion.
     !@ <objectMethod>
     !@   <object>mergerTree</object>
     !@   <method>promoteNode(thisNode)</method>
     !@   <description>Promotes {\tt thisNode} to take the place of its parent node.</description>
     !@ </objectMethod>
     procedure                 :: promoteNode   => Tree_Node_Promote
     ! Node evolution.
     !@ <objectMethod>
     !@   <object>mergerTree</object>
     !@   <method>evolveNode(thisNode)</method>
     !@   <description>Evolves {\tt thisNode} over a suitable timestep.</description>
     !@ </objectMethod>
     procedure                 :: evolveNode    => Tree_Node_Evolve
     ! Node merging.
     !@ <objectMethod>
     !@   <object>mergerTree</object>
     !@   <method>mergeNode(thisNode)</method>
     !@   <description>Causes {\tt thisNode} to undergo a node merger event with its parent node.</description>
     !@ </objectMethod>
     procedure                 :: mergeNode     => Events_Node_Merger
     ! Node locating.
     !@ <objectMethod>
     !@   <object>mergerTree</object>
     !@   <method>getNode(thisNode)</method>
     !@   <description>Returns a pointer to the node with given index in the merger tree, or a null pointer if no such node exists.</description>
     !@ </objectMethod>
     procedure                 :: getNode       => Tree_Node_Get
  end type mergerTree

  ! Flag to indicate if the tree node create routines have been initialized.
  logical :: treeNodeCreateInitialized=.false.

  ! Count of number of different possible component types.
  integer :: componentTypesCount=0
  
  ! Flag to indicate if the node evolver method has been initialized.
  logical :: evolverInitialized=.false.

  ! Parameters controlling the accuracy of ODE solving.
  double precision :: odeToleranceAbsolute, odeToleranceRelative

  ! Arrays that point to node properties and their derivatives.
  integer                                     :: nPropertiesMax=0, nProperties
  double precision, allocatable, dimension(:) :: propertyValues,propertyScales
  !$omp threadprivate(nPropertiesMax,nProperties,propertyValues,propertyScales)
#ifdef PROFILE
  logical                                     :: profileOdeEvolver
  integer,          allocatable, dimension(:) :: propertyComponent,propertyObject,propertyIndex
  !$omp threadprivate(propertyComponent,propertyObject,propertyIndex)
#endif

  ! Module global pointer to the node being processed.
  integer(kind=kind_int8)                     :: activeTreeIndex
  type(treeNode),   pointer                   :: activeNode
  !$omp threadprivate(activeNode,activeTreeIndex)

  ! Variables to track interrupt events.
  logical                                          :: firstInterruptFound
  double precision                                 :: firstInterruptTime
  procedure(Interrupt_Procedure_Template), pointer :: firstInterruptProcedure
  !$omp threadprivate(firstInterruptFound,firstInterruptTime,firstInterruptProcedure)

  ! Flag to indicate if node merging event method has been initialized.  
  logical              :: nodeMergersInitialized=.false.

  ! Name of node mergers method used.
  type(varying_string) :: nodeMergersMethod

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure(Node_Mergers_Template), pointer :: Events_Node_Merger_Do => null()
  abstract interface
     subroutine Node_Mergers_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end subroutine Node_Mergers_Template
  end interface

  ! Pointer to an error handler for failures in the ODE solver.
  procedure(), pointer :: Galacticus_ODE_Error_Handler => Tree_Node_ODEs_Error_Handler

contains

  subroutine Merger_Tree_Destroy(thisTree)
    !% Destroys the entire merger tree.
    implicit none
    class(mergerTree), intent(inout) :: thisTree

    select type (thisTree)
    type is (mergerTree)
       ! Destroy all nodes.
       if (associated(thisTree%baseNode)) call thisTree%destroyBranch(thisTree%baseNode)
       ! Destroy the HDF5 group associated with this tree.
       call thisTree%hdf5Group%destroy()
    end select
    return
  end subroutine Merger_Tree_Destroy

  subroutine Merger_Tree_Destroy_Branch(thisTree,thisNode)
    !% Destroy a branch of a tree which begins at {\tt thisNode}.
    use Tree_Nodes
    implicit none
    class(mergerTree),          intent(inout) :: thisTree
    type(treeNode),    pointer, intent(inout) :: thisNode
    type(treeNode),    pointer                :: destroyNode,nextNode

    ! Descend to the tip of the branch.
    call thisNode%walkBranchWithSatellites(thisNode,nextNode)
    ! Loop over all tree nodes.
    do while (associated(nextNode))
       ! Keep of a record of the current node, so that we can destroy it.
       destroyNode => nextNode

       ! Walk to the next node in the tree.
       call destroyNode%walkBranchWithSatellites(thisNode,nextNode)

       ! If the node about to be destroyed is the primary progenitor of its parent we must move the child pointer of the parent to
       ! point to the node's sibling. This is necessary as parent-child pointers are used to establish satellite status and so
       ! will be utilized when walking the tree. Failure to do this can result in attempts to use dangling pointers.
       if (associated(destroyNode%parentNode)) then
          if (associated(destroyNode%parentNode%childNode,destroyNode)) destroyNode%parentNode%childNode => destroyNode%siblingNode
       end if

       ! Destroy the current node.
       call destroyNode%destroy
    end do
    ! Destroy the base node of the branch.
    if (associated(thisNode%parentNode)) then
       if (associated(thisNode%parentNode%childNode,thisNode)) thisNode%parentNode%childNode => thisNode%siblingNode
    end if
    call thisNode%destroy
    return
  end subroutine Merger_Tree_Destroy_Branch

  subroutine Tree_Node_Create(thisTree,thisNode,index)
    !% Return a pointer to a newly created and initialized tree node.
    use Galacticus_Error
    use Memory_Management
    implicit none
    class(mergerTree),                intent(inout)        :: thisTree
    type(treeNode),          pointer, intent(inout)        :: thisNode
    integer(kind=kind_int8),          intent(in), optional :: index
    integer                                                :: allocErr

    ! Initialize tree node methods if necessary.
    call Tree_Node_Create_Initialize

    ! Allocate the object.
    allocate(thisNode,stat=allocErr)
    if (allocErr/=0) call Galacticus_Error_Report('Tree_Node_Create','unable to allocate node')
    call Memory_Usage_Record(sizeof(thisNode),memoryType=memoryTypeNodes)

    ! Initialize list of component indices.
    allocate(thisNode%componentIndex(componentTypesCount))
    call Memory_Usage_Record(sizeof(thisNode%componentIndex),memoryType=memoryTypeNodes)
    thisNode%componentIndex=-1

    ! Ensure pointers are nullified.
    nullify(thisNode%parentNode,thisNode%childNode,thisNode%siblingNode,thisNode%satelliteNode,thisNode%mergeNode&
         &,thisNode%mergeeNode,thisNode%nextMergee,thisNode%formationNode)

    ! Assign index if supplied.
    if (present(index)) call thisNode%indexSet(index)

    ! Assign a unique ID.
    call thisNode%uniqueIDSet(Tree_Nodes_New_Unique_ID())

    return
  end subroutine Tree_Node_Create

  subroutine Tree_Node_Create_Initialize
    !% Initializes tree node create by calling all relevant initialization routines.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    !# <include directive="treeNodeCreateInitialize" type="optionDefinitions">
    include 'objects.tree_node.create.definitions.inc'
    !# </include>
    
    if (.not.treeNodeCreateInitialized) then
       !$omp critical (Tree_Node_Create_Initialize)
       if (.not.treeNodeCreateInitialized) then
          
          ! Read all parameters needed by methods.
          !# <include directive="treeNodeCreateInitialize" type="optionNames">
          include 'objects.tree_node.create.parameters.inc'
          !# </include>
          
          ! Initialize rate adjust and compute pointers to dummy implementations.
          !# <include directive="treeNodeMethodsPointer" type="initializeMethods">
          include 'objects.tree_node.initializeMethods.inc'
          !# </include>
          
          ! Call all routines to initialize tree node create.
          !# <include directive="treeNodeCreateInitialize" type="code" action="subroutine">
          !#  <subroutineArgs>componentTypesCount</subroutineArgs>
          include 'objects.tree_node.create.initialize.inc'
          !# </include>
          
          ! Flag that tree node methods are now initialized.
          treeNodeCreateInitialized=.true.
          
       end if
       !$omp end critical (Tree_Node_Create_Initialize)
    end if
    return
  end subroutine Tree_Node_Create_Initialize

  subroutine Tree_Node_Promote(thisTree,thisNode)
    !% Transfer the properties of {\tt thisNode} to its parent node, then destroy it.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    !# <include directive="nodePromotionTask" type="moduleUse">
    include 'objects.tree_node.promote.modules.inc'
    !# </include>
    implicit none
    class(mergerTree),            intent(inout) :: thisTree
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(treeNode),      pointer                :: parentNode,satelliteNode
    type(varying_string)                        :: message

    ! Get pointer to parent node.
    parentNode => thisNode%parentNode

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Promoting node '
       message=message//thisNode%index()//' to '//parentNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if

    ! Perform any processing necessary before this halo is promoted.
    !# <include directive="nodePromotionTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.promote.inc'
    !# </include>

    ! Move the components of thisNode to the parent.
    call parentNode%destroyAllComponents
    call Move_Alloc(thisNode%componentIndex,parentNode%componentIndex)
    call Move_Alloc(thisNode%components    ,parentNode%components    )
    
    ! Copy any formation node data to the parent, and update the formation node's parentNode pointer to point to the new parent.
    if (associated(thisNode%formationNode)) then
       if (associated(parentNode%formationNode)) then
          call parentNode%formationNode%destroy()
          deallocate(parentNode%formationNode)
       end if
       allocate(parentNode%formationNode)
       call parentNode%formationNode%copy(thisNode%formationNode)
       parentNode%formationNode%parentNode => parentNode
    end if

    ! Transfer any satellite nodes to the parent.
    if (associated(thisNode%satelliteNode)) then
       ! Attach the satellite nodes to the parent.
       if (associated(parentNode%satelliteNode)) then
          ! Find the last satellite of the parent node.
          call parentNode%lastSatellite(satelliteNode)
          satelliteNode%siblingNode => thisNode%satelliteNode
       else
          parentNode%satelliteNode  => thisNode%satelliteNode
       end if
       ! Get the first satellite of thisNode.
       satelliteNode => thisNode%satelliteNode
       do while (associated(satelliteNode))
          ! Set the parent node for this satellite to the parent.
          satelliteNode%parentNode => parentNode
          satelliteNode => satelliteNode%siblingNode
       end do
    end if

    ! Nullify the child pointer for the parent.
    parentNode%childNode => null()

    ! Destroy the node.
    call thisNode%destroy

    return
  end subroutine Tree_Node_Promote

  subroutine Tree_Node_Evolve_Initialize
    !% Initializes the tree evolving routines by reading in parameters
    use Input_Parameters

    ! Initialize if necessary.
    !$omp critical (Tree_Node_Evolve_Initialize)
    if (.not.evolverInitialized) then
       ! Get tolerance values for the ODE solver.
       !@ <inputParameter>
       !@   <name>odeToleranceAbsolute</name>
       !@   <defaultValue>$0.01$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The absolute tolerance used in solving differential equations for node evolution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>timeStepping</group>
       !@ </inputParameter>
       call Get_Input_Parameter('odeToleranceAbsolute',odeToleranceAbsolute,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>odeToleranceRelative</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The relative tolerance used in solving differential equations for node evolution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>timeStepping</group>
       !@ </inputParameter>
       call Get_Input_Parameter('odeToleranceRelative',odeToleranceRelative,defaultValue=1.0d-2)
#ifdef PROFILE
       !@ <inputParameter>
       !@   <name>profileOdeEvolver</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to profile the ODE evolver.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>timeStepping</group>
       !@ </inputParameter>
       call Get_Input_Parameter('profileOdeEvolver',profileOdeEvolver,defaultValue=.false.)
#endif
       ! Flag that the module is now initialized.
       evolverInitialized=.true.
    end if
    !$omp end critical (Tree_Node_Evolve_Initialize)
    return
  end subroutine Tree_Node_Evolve_Initialize

  subroutine Tree_Node_Evolve(thisTree,thisNode,endTime,interrupted,interruptProcedure)
    !% Evolves {\tt thisNode} to time {\tt endTime}, or until evolution is interrupted.
    use ODEIV2_Solver
    use FODEIV2
    use Memory_Management
    !# <include directive="postEvolveTask" type="moduleUse">
    include 'objects.tree_node.post_evolve.modules.inc'
    !# </include>
    implicit none
    class(mergerTree),                       intent(inout)          :: thisTree
    type(treeNode),                          intent(inout), pointer :: thisNode
    double precision,                        intent(in)             :: endTime
    logical,                                 intent(out)            :: interrupted
    procedure(Interrupt_Procedure_Template), intent(out),   pointer :: interruptProcedure
    integer,                                 save                   :: nPropertiesPrevious=-1
    !$omp threadprivate(nPropertiesPrevious)
    double precision                                                :: startTimeThisNode
    ! Variables used in the ODE solver.
    type(fodeiv2_system    ), save                                  :: ode2System
    type(fodeiv2_driver    ), save                                  :: ode2Driver
    logical                 , save                                  :: odeReset
    !$omp threadprivate(ode2System,ode2Driver,odeReset)
    type(c_ptr)                                                     :: parameterPointer
#ifdef PROFILE
    type(c_funptr)                                                  :: Error_Analyzer
#endif

    ! Initialize.
    call Tree_Node_Evolve_Initialize()

    ! Determine the end time for this node - either the specified end time, or the time associated with the parent node, whichever
    ! occurs first.
    startTimeThisNode=Tree_Node_Time(thisNode)

    ! Find number of variables to evolve for this node.
    nProperties=Count_Properties(thisNode)

    ! Allocate pointer arrays if necessary.
    if (nProperties > nPropertiesMax) then
       if (allocated(propertyValues)) then
          call Memory_Usage_Record(sizeof(propertyValues),addRemove=-1)
          deallocate(propertyValues)
          call Memory_Usage_Record(sizeof(propertyScales),addRemove=-1)
          deallocate(propertyScales)
#ifdef PROFILE
          if (profileOdeEvolver) then
             call Memory_Usage_Record(sizeof(propertyComponent),addRemove=-1)
             deallocate(propertyComponent)
             call Memory_Usage_Record(sizeof(propertyObject   ),addRemove=-1)
             deallocate(propertyObject   )
             call Memory_Usage_Record(sizeof(propertyIndex    ),addRemove=-1)
             deallocate(propertyIndex    )
          end if
#endif
       end if
       allocate(propertyValues(nProperties))
       call Memory_Usage_Record(sizeof(propertyValues))
       allocate(propertyScales(nProperties))
       call Memory_Usage_Record(sizeof(propertyScales))
#ifdef PROFILE
       if (profileOdeEvolver) then
          allocate(propertyComponent(nProperties))
          call Memory_Usage_Record(sizeof(propertyComponent))
          allocate(propertyObject   (nProperties))
          call Memory_Usage_Record(sizeof(propertyObject   ))
          allocate(propertyIndex    (nProperties))
          call Memory_Usage_Record(sizeof(propertyIndex    ))
       end if
#endif
       nPropertiesMax=nProperties
    end if

    ! Assign pointers to node variables.
    call Map_Properties_To_ODE_Array(thisNode,propertyValues,labelValue)

    ! Compute scales for all properties and extract from the node.
    call Initialize_Scales_to_Unity(thisNode)
    !# <include directive="scaleSetTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.set_scale.inc'
    !# </include>
    call Map_Properties_To_ODE_Array(thisNode,propertyScales,labelScale)

    ! Assign module global pointer to this node.
    activeTreeIndex=  thisTree%index
    activeNode     => thisNode

    ! Reset interrupt variables.
    firstInterruptFound     =  .false.
    firstInterruptTime      =  0.0d0
    firstInterruptProcedure => null()

    ! Call ODE solver routines.
    startTimeThisNode=Tree_Node_Time(thisNode)
#ifdef PROFILE
    if (profileOdeEvolver) then
       Error_Analyzer=c_funloc(Tree_Node_Evolve_Error_Analyzer)
    else
       Error_Analyzer=C_NULL_FUNPTR
    end if
#endif
    if (nPropertiesPrevious > 0 .and. .not.odeReset) call ODEIV2_Solver_Free(ode2Driver,ode2System)
    odeReset=.true.
    nPropertiesPrevious=nProperties
    if (startTimeThisNode /= endTime)                                   &
         & call ODEIV2_Solve(                                           &
         &                   ode2Driver,ode2System                    , &
         &                   startTimeThisNode,endTime                , &
         &                   nProperties                              , &
         &                   propertyValues                           , &
         &                   Tree_Node_ODEs                           , &
         &                   parameterPointer                         , &
         &                   odeToleranceAbsolute,odeToleranceRelative, &
#ifdef PROFILE
         &                   Error_Analyzer                           , &
#endif
         &                   propertyScales                           , &
         &                   reset=odeReset                           , &
         &                   errorHandler=Galacticus_ODE_Error_Handler  &
         &                  )

    ! Extract values.
    call Map_Properties_From_ODE_Array(thisNode,propertyValues)

    ! Ensure that the maximum time has not been exceed (can happen due to rounding errors).
    if (Tree_Node_Time(thisNode) > endTime) call Tree_Node_Time_Set(thisNode,endTime)

    ! Flag interruption if one occurred.
    if (firstInterruptTime /= 0.0d0) then
       interrupted=.true.
       interruptProcedure => firstInterruptProcedure
    else
       interrupted=.false.
    end if

    ! Call routines to perform any post-evolution tasks.
    !# <include directive="postEvolveTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.post_evolve.inc'
    !# </include>
 
    return
  end subroutine Tree_Node_Evolve

  logical function Tree_Node_Is_Accurate(valueNode,valueExpected)  
    !% Return true if a tree node property is within expected accuracy of a given value.
    use Numerical_Comparison
    implicit none
    double precision, intent(in) :: valueNode,valueExpected

    ! Initialize.
    call Tree_Node_Evolve_Initialize()

    ! Compare values.
    Tree_Node_Is_Accurate=Values_Agree(valueNode,valueExpected,odeToleranceAbsolute,odeToleranceRelative)
    return
  end function Tree_Node_Is_Accurate

  function Tree_Node_ODEs(time,y,dydt,parameterPointer) bind(c)
    !% Function which evaluates the set of ODEs for the evolution of a specific node.
    use ODE_Solver_Error_Codes
    implicit none
    integer(c_int)                                       :: Tree_Node_ODEs
    real(c_double),                          value       :: time
    real(c_double),                          intent(in)  :: y(nProperties)
    real(c_double),                          intent(out) :: dydt(nProperties)
    type(c_ptr),                             value       :: parameterPointer
    logical                                              :: interrupt
    procedure(Interrupt_Procedure_Template), pointer     :: interruptProcedure

    ! Extract values.
    call Map_Properties_From_ODE_Array(activeNode,y)

    ! Set derivatives to zero initially.
    call Initialize_Derivatives_to_Zero(activeNode)

    if (firstInterruptFound .and. time >= firstInterruptTime) then
       ! Already beyond the location of the first interrupt, simply return zero derivatives.
       dydt=0.0d0
    else
       ! Compute derivatives.
       call Tree_Node_Compute_Derivatives(activeNode,interrupt,interruptProcedure)
       
       ! Check whether an interrupt has been requested.
       select case (interrupt)
       case (.false.)
          ! No interrupt - place derivatives into ODE arrays.
          call Map_Properties_To_ODE_Array(activeNode,dydt,labelDerivative)
       case (.true.)
          ! Interrupt requested - freeze evolution and store the interrupt if it is the earliest one to occur.
          dydt=0.0d0
          if (time < firstInterruptTime .or. .not.firstInterruptFound) then
             firstInterruptFound     =  .true.
             firstInterruptTime      =  time
             firstInterruptProcedure => interruptProcedure
             ! Let the ODE solver know that an interrupt occured, and when it happened.
             Tree_Node_ODEs          =  odeSolverInterrupt
             interruptedAtX          =  time
             return
          end if
       end select
    end if

    ! Return success.
    Tree_Node_ODEs=FGSL_Success
    return
  end function Tree_Node_ODEs

  subroutine Initialize_Derivatives_to_Zero(thisNode)
    !% Initialize the derivatives of all properties of {\tt thisNode} to zero.
    implicit none
    type(treeNode), intent(inout) :: thisNode
    integer                       :: iComponent,iHistory,iInstance

    do    iComponent=1,size(thisNode%components                     )
       do iInstance =1,size(thisNode%components(iComponent)%instance)
          if     (allocated(thisNode &
               &                   %components(iComponent          )       &
               &                   %instance  (iInstance           )       &
               &                   %properties                             &
               &           )                                               &
               & )          thisNode                                       &
               &                   %components(iComponent          )       &
               &                   %instance  (iInstance           )       &
               &                   %properties(:,propertyDerivative)=0.0d0
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories )) then
             do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time))&
                     & thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%rates(:,:)=0.0d0
             end do
          end if
       end do
    end do
    return
  end subroutine Initialize_Derivatives_to_Zero

  subroutine Initialize_Scales_to_Unity(thisNode)
    !% Initialize the scale of all properties of {\tt thisNode} to unity.
    implicit none
    type(treeNode),  intent(inout) :: thisNode
    integer                        :: iComponent,iHistory,iInstance

    do    iComponent=1,size(thisNode%components                     )
       do iInstance =1,size(thisNode%components(iComponent)%instance)
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%properties)) thisNode%components(iComponent)%instance(iInstance)%properties(:&
               &,propertyScale)=1.0d0
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories )) then
             do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time))&
                     & thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%scales(:,:)=1.0d0
             end do
          end if
       end do
    end do
    return
  end subroutine Initialize_Scales_to_Unity

  integer function Count_Properties(thisNode)
    !% Determine the number of active properties for {\tt thisNode}.
    implicit none
    type(treeNode),  intent(in) :: thisNode
    integer                     :: iComponent,iHistory,iInstance

    Count_Properties=0
    do    iComponent=1,size(thisNode%components                     )
       do iInstance =1,size(thisNode%components(iComponent)%instance)
          if     (allocated(               thisNode                       &
               &                                  %components(iComponent) &
               &                                  %instance  (iInstance ) &
               &                                  %properties             &
               &           )                                              &
               & ) Count_Properties= Count_Properties                     &
               &                    +size( thisNode                       &
               &                                  %components(iComponent) &
               &                                  %instance  (iInstance ) &
               &                                  %properties             &
               &                          ,dim=1                          &
               &                         )
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories)) then
             do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)) Count_Properties=Count_Properties&
                     &+size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data,dim=1) &
                     &*size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data,dim=2)
             end do
          end if
       end do
    end do
    return
  end function Count_Properties

  subroutine Map_Properties_To_ODE_Array(thisNode,propertyValuesODE,labelType)
    !% Copies values of properties from {\tt thisNode} to a single array for use in the ODE solver routines.
    implicit none
    double precision, intent(out)          :: propertyValuesODE(:)
    type(treeNode),   intent(in)           :: thisNode
    integer,          intent(in), optional :: labelType
    integer                                :: propertyCounter,iComponent,iProperty,labelTypeActual,propertyTypeActual,iHistory,iTime,iValue,iInstance
    
    ! Find which type of property (value or derivative) is required.
    if (present(labelType)) then
       select case (labelType)
       case (labelValue     )
          propertyTypeActual=propertyValue
       case (labelDerivative)
          propertyTypeActual=propertyDerivative
       case (labelScale     )
          propertyTypeActual=propertyScale
       end select
       labelTypeActual   =labelType
    else
       propertyTypeActual=propertyValue
       labelTypeActual   =labelValue
    end if

    propertyCounter=0
    do    iComponent=1,size(thisNode%components                     )
       do iInstance =1,size(thisNode%components(iComponent)%instance)
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%properties)) then
             do iProperty=1,size(thisNode%components(iComponent)%instance(iInstance)%properties,dim=1)
                propertyCounter=propertyCounter+1
                propertyValuesODE(propertyCounter)=thisNode%components(iComponent)%instance(iInstance)%properties(iProperty,propertyTypeActual)
#ifdef PROFILE
                if (profileOdeEvolver) then
                   propertyComponent(propertyCounter)=iComponent
                   propertyObject   (propertyCounter)=objectTypeProperty
                   propertyIndex    (propertyCounter)=iProperty
                end if
#endif
             end do
          end if
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories)) then
             do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)) then
                   do iValue=1,size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data,dim=2)
                      do iTime=1,size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)
                         propertyCounter=propertyCounter+1
                         select case (labelTypeActual)
                         case (labelValue     )
                            propertyValuesODE(propertyCounter)=thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data  (iTime,iValue)
                         case (labelDerivative)
                            propertyValuesODE(propertyCounter)=thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%rates (iTime,iValue)
                         case (labelScale     )
                            propertyValuesODE(propertyCounter)=thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%scales(iTime,iValue)
                         end select
#ifdef PROFILE
                         if (profileOdeEvolver) then
                            propertyComponent(propertyCounter)=iComponent
                            propertyObject   (propertyCounter)=objectTypeHistory
                            propertyIndex    (propertyCounter)=iHistory
                         end if
#endif
                      end do
                   end do
                end if
             end do
          end if
       end do
    end do

#ifdef PROFILE
    ! Invert the component index mapping.
    if (profileOdeEvolver) then
       do while (propertyCounter > 0)
          do iProperty=1,size(thisNode%componentIndex)
             if (thisNode%componentIndex(iProperty) == propertyComponent(propertyCounter)) then
                propertyComponent(propertyCounter)=iProperty
                exit
             end if
          end do
          propertyCounter=propertyCounter-1
       end do
    end if
#endif
    return
  end subroutine Map_Properties_To_ODE_Array

  subroutine Map_Properties_From_ODE_Array(thisNode,propertyValuesODE)
    !% Copies values of properties from {\tt thisNode} to a single array for use in the ODE solver routines.
    implicit none
    double precision, intent(in)    :: propertyValuesODE(:)
    type(treeNode),   intent(inout) :: thisNode
    integer                         :: propertyCounter,iComponent,iProperty,iHistory,iValue,iInstance,entryCount
    
    propertyCounter=0
    do    iComponent=1,size(thisNode%components                     )
       do iInstance =1,size(thisNode%components(iComponent)%instance)
          if     (allocated(thisNode                       &
               &                   %components(iComponent) &
               &                   %instance  (iInstance ) &
               &                   %properties             &
               &           )                               &
               & ) then
             do iProperty=1,size(thisNode%components(iComponent)%instance(iInstance)%properties,dim=1)
                propertyCounter=propertyCounter+1
                thisNode%components(iComponent)%instance(iInstance)%properties(iProperty,propertyValue)=propertyValuesODE(propertyCounter)
             end do
          end if
          if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories)) then
             do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)) then
                   do iValue=1,size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data,dim=2)
                      entryCount=size(thisNode                       &
                           &                 %components(iComponent) &
                           &                 %instance  (iInstance ) &
                           &                 %histories (iHistory  ) &
                           &                 %time                   &
                           &         )
                      thisNode                                                                                         &
                           & %components(  iComponent)                                                                 &
                           & %instance  (  iInstance )                                                                 &
                           & %histories (  iHistory  )                                                                 &
                           & %data      (:,iValue    )=propertyValuesODE(                                              &
                           &                                              propertyCounter+1                            &
                           &                                             :propertyCounter+entryCount                   &
                           &                                            )
                      propertyCounter=propertyCounter+entryCount
                   end do
                end if
             end do
          end if
       end do
    end do
    return
  end subroutine Map_Properties_From_ODE_Array

  subroutine Tree_Node_Compute_Derivatives(thisNode,interrupt,interruptProcedureReturn)
    !% Call routines to set alls derivatives for {\tt thisNode}.
    use Galacticus_Calculations_Resets
    !# <include directive="preDerivativeTask" type="moduleUse">
    include 'objects.merger_trees.prederivative.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(out)   :: interrupt
    procedure(),    pointer, intent(out)   :: interruptProcedureReturn
    procedure(),    pointer                :: interruptProcedure

    ! Initialize interrupt status.
    interrupt=.false.
    interruptProcedure => null()

    ! Call component routines to indicate that derivative calculation is commencing.
    call Galacticus_Calculations_Reset(thisNode)

    ! Call routines to perform any pre-derivative calculations.
    !# <include directive="preDerivativeTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.merger_trees.prederivative.tasks.inc'
    !# </include>
    
    ! Call component routines to compute derivatives.
    !# <include directive="treeNodeMethodsPointer" type="derivatives">
    include 'objects.tree_node.derivatives.inc'
    !# </include>

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure

    return
  end subroutine Tree_Node_Compute_Derivatives

  subroutine Tree_Node_ODEs_Error_Handler()
    !% Handles errors in the ODE solver when evolving \glc\ nodes. Dumps the content of the node.
    use Tree_Nodes_Dump
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string) :: message
    
    message="ODE solver failed in tree #"
    message=message//activeTreeIndex
    call Galacticus_Display_Message(message)
    call Node_Dump(activeNode)
    return
  end subroutine Tree_Node_ODEs_Error_Handler
  
#ifdef PROFILE
  subroutine Tree_Node_Evolve_Error_Analyzer(currentPropertyValue,currentPropertyError,timeStep,stepStatus) bind(c)
    !% Profiles ODE solver step sizes and errors.
    use FGSL
    use ISO_Varying_String
    use Galacticus_Meta_Evolver_Profiler
    !# <include directive="decodePropertyIdentifiersTask" type="moduleUse">
    include 'objects.merger_trees.decode_property_identifiers.modules.inc'
    !# </include>
    implicit none
    real(c_double),   intent(in), dimension(nProperties) :: currentPropertyValue
    real(c_double),   intent(in), dimension(nProperties) :: currentPropertyError
    real(c_double),   intent(in), value                  :: timeStep
    integer(c_int),   intent(in), value                  :: stepStatus
    double precision                                     :: scaledErrorMaximum,scaledError,scale
    integer                                              :: iProperty,limitingProperty
    logical                                              :: matchedProperty
    type(varying_string)                                 :: propertyName

    ! If the step was not good, return immediately.
    if (stepStatus /= FGSL_SUCCESS) return

    ! Find the property with the largest error (i.e. that which is limiting the step).
    scaledErrorMaximum=0.0d0
    do iProperty=1,nProperties
       scale=odeToleranceAbsolute*propertyScales(iProperty)+odeToleranceRelative*dabs(currentPropertyValue(iProperty))
       scaledError=dabs(currentPropertyError(iProperty))/scale
       if (scaledError > scaledErrorMaximum) then
          scaledErrorMaximum=scaledError
          limitingProperty=iProperty
       end if
    end do
    ! Decode the step limiting property.
    matchedProperty=.false.
    propertyName="unknown"
    !# <include directive="decodePropertyIdentifiersTask" type="code" action="subroutine">
    !#  <subroutineArgs>propertyComponent(limitingProperty),propertyObject(limitingProperty),propertyIndex(limitingProperty),matchedProperty,propertyName</subroutineArgs>
    include 'objects.merger_trees.decode_property_identifiers.inc'
    !# </include>

    ! Record this information.
    call Galacticus_Meta_Evolver_Profile(timeStep,propertyName)
    return
  end subroutine Tree_Node_Evolve_Error_Analyzer
#endif

  subroutine Events_Node_Merger(thisTree,thisNode)
    !% Handles instances where {\tt thisNode} is about to merge with its parent node.
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
    !# <include directive="nodeMergersMethod" type="moduleUse">
    include 'events.node_mergers.moduleUse.inc'
    !# </include>
    !# <include directive="nodeMergerTask" type="moduleUse">
    include 'events.node_mergers.process.modules.inc'
    !# </include>
    implicit none
    class(mergerTree),   intent(inout)          :: thisTree
    type(treeNode),      intent(inout), pointer :: thisNode
    type(varying_string)                        :: message
    
    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Making node '
       message=message//thisNode%index()//' a satellite in '//thisNode%parentNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if

    ! Call subroutines to perform any necessary processing prior to this node merger event.
    !# <include directive="nodeMergerTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'events.node_mergers.process.inc'
    !# </include>
 
    !$omp critical (Events_Node_Merger_Initialize)
    if (.not.nodeMergersInitialized) then
       ! Get the node mergers method parameter.
       !@ <inputParameter>
       !@   <name>nodeMergersMethod</name>
       !@   <defaultValue>singleLevelHierarchy</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used for handling node merger events.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeMergersMethod',nodeMergersMethod,defaultValue='singleLevelHierarchy')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="nodeMergersMethod" type="code" action="subroutine">
       !#  <subroutineArgs>nodeMergersMethod,Events_Node_Merger_Do</subroutineArgs>
       include 'events.node_mergers.inc'
       !# </include>
       if (.not.associated(Events_Node_Merger_Do)) call Galacticus_Error_Report('Events_Node_Merger','method '&
            &//char(nodeMergersMethod)//' is unrecognized')
       nodeMergersInitialized=.true.
    end if
    !$omp end critical (Events_Node_Merger_Initialize)

    ! Call the routine to perform the merger.
    call Events_Node_Merger_Do(thisNode)

    return
  end subroutine Events_Node_Merger

  
  function Tree_Node_Get(thisTree,nodeIndex)
    !% Return a pointer to a node in {\tt thisTree} given the index of the node.
    use Tree_Nodes
    implicit none
    class(mergerTree),       intent(in) :: thisTree
    integer(kind=kind_int8), intent(in) :: nodeIndex
    type(treeNode),          pointer    :: Tree_Node_Get,thisNode

    Tree_Node_Get => null()
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       if (thisNode%index() == nodeIndex) then
          Tree_Node_Get => thisNode
          return
       end if
       ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
       call thisNode%walkTreeWithSatellites(thisNode)
    end do
    return
  end function Tree_Node_Get

end module Merger_Trees
