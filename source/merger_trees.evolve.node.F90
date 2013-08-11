!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements evolution of a single node in a merger tree.

module Merger_Trees_Evolve_Node
  !% Implements evolution of a single node in a merger tree.
  use Galacticus_Nodes
  use Kind_Numbers
  use ISO_Varying_String
  use FODEIV2
  private
  public :: Tree_Node_Evolve, Tree_Node_Promote, Events_Node_Merger, Tree_Node_Is_Accurate

  ! Flag to indicate if the node evolver method has been initialized.
  logical                                     :: evolverInitialized       =.false.

  ! Parameters controlling the accuracy of ODE solving.
  double precision                            :: odeToleranceAbsolute             , odeToleranceRelative

  ! Arrays that point to node properties and their derivatives.
  integer                                     :: nProperties                      , nPropertiesMax      =0
  double precision, allocatable, dimension(:) :: propertyScales                   , propertyValues
  !$omp threadprivate(nPropertiesMax,nProperties,propertyValues,propertyScales)
#ifdef PROFILE
  logical :: profileOdeEvolver
#endif

  ! Module global pointer to the node being processed.
  integer         (kind=kind_int8              )          :: activeTreeIndex
  type            (treeNode                    ), pointer :: activeNode
  !$omp threadprivate(activeNode,activeTreeIndex)
  ! Variables to track interrupt events.
  logical                                                 :: firstInterruptFound
  double precision                                        :: firstInterruptTime
  procedure       (Interrupt_Procedure_Template), pointer :: firstInterruptProcedure
  !$omp threadprivate(firstInterruptFound,firstInterruptTime,firstInterruptProcedure)
  ! Flag to indicate if node merging event method has been initialized.
  logical                                                 :: nodeMergersInitialized  =.false.

  ! Name of node mergers method used.
  type            (varying_string              )          :: nodeMergersMethod

  ! The algorithm to use for ODE solving.
  type            (fodeiv2_step_type           )          :: Galacticus_ODE_Algorithm

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure       (Node_Mergers_Template       ), pointer :: Events_Node_Merger_Do   =>null()
  abstract interface
     subroutine Node_Mergers_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end subroutine Node_Mergers_Template
  end interface

  ! Pointer to an error handler for failures in the ODE solver.
  procedure(), pointer :: Galacticus_ODE_Error_Handler=>Tree_Node_ODEs_Error_Handler

contains

  subroutine Tree_Node_Evolve_Initialize
    !% Initializes the tree evolving routines by reading in parameters
    use Input_Parameters
    use Galacticus_Error
    type(varying_string) :: odeAlgorithm

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
       !@ <inputParameter>
       !@   <name>odeAlgorithm</name>
       !@   <defaultValue>Runge-Kutta-Cash-Karp</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The algorithm to use in the ODE solver.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>timeStepping</group>
       !@ </inputParameter>
       call Get_Input_Parameter('odeAlgorithm',odeAlgorithm,defaultValue='Runge-Kutta-Cash-Karp')
       select case (char(odeAlgorithm))
       case ('Runge-Kutta-Cash-Karp')
          Galacticus_ODE_Algorithm=Fodeiv2_Step_RKCK
       case ('Runge-Kutta-Second-Order')
          Galacticus_ODE_Algorithm=Fodeiv2_Step_RK2
       case default
          call Galacticus_Error_Report('Tree_Node_Evolve_Initialize','odeAlgorithm is unrecognized')
       end select
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
    use Memory_Management
    use Galacticus_Calculations_Resets
    use, intrinsic :: ISO_C_Binding
    !# <include directive="preEvolveTask" type="moduleUse">
    include 'objects.tree_node.pre_evolve.modules.inc'
    !# </include>
    !# <include directive="postEvolveTask" type="moduleUse">
    include 'objects.tree_node.post_evolve.modules.inc'
    !# </include>
    !# <include directive="scaleSetTask" type="moduleUse">
    include 'objects.tree_node.set_scale.modules.inc'
    !# </include>
    implicit none
    class           (mergerTree                  )      , intent(inout)          :: thisTree
    type            (treeNode                    )      , intent(inout), pointer :: thisNode
    double precision                                    , intent(in   )          :: endTime
    logical                                             , intent(  out)          :: interrupted
    procedure       (Interrupt_Procedure_Template)      , intent(  out), pointer :: interruptProcedure
    class           (nodeComponentBasic          )                     , pointer :: basicComponent
    integer                                       , save                         :: nPropertiesPrevious=-1
    !$omp threadprivate(nPropertiesPrevious)
    double precision                                                             :: startTimeThisNode
    ! Variables used in the ODE solver.
    type            (fodeiv2_system              ), save                         :: ode2System
    type            (fodeiv2_driver              ), save                         :: ode2Driver
    logical                                       , save                         :: odeReset
    !$omp threadprivate(ode2System,ode2Driver,odeReset)
    type            (c_ptr                       )                               :: parameterPointer
#ifdef PROFILE
    type(c_funptr) :: Error_Analyzer
#endif

    ! Initialize.
    call Tree_Node_Evolve_Initialize()

    ! Call routines to perform any pre-evolution tasks.
    !# <include directive="preEvolveTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'objects.tree_node.pre_evolve.inc'
    !# </include>

    ! Determine the end time for this node - either the specified end time, or the time associated with the parent node, whichever
    ! occurs first.
    basicComponent => thisNode%basic()
    startTimeThisNode=basicComponent%time()

    ! Find number of variables to evolve for this node.
    nProperties=thisNode%serializeCount()

    ! Allocate pointer arrays if necessary.
    if (nProperties > nPropertiesMax) then
       if (allocated(propertyValues)) then
          call Memory_Usage_Record(sizeof(propertyValues),addRemove=-1)
          deallocate(propertyValues)
          call Memory_Usage_Record(sizeof(propertyScales),addRemove=-1)
          deallocate(propertyScales)
       end if
       allocate(propertyValues(nProperties))
       call Memory_Usage_Record(sizeof(propertyValues))
       allocate(propertyScales(nProperties))
       call Memory_Usage_Record(sizeof(propertyScales))
       nPropertiesMax=nProperties
    end if

    ! Assign pointers to node variables.
    call thisNode%serializeValues(propertyValues)

    ! Compute scales for all properties and extract from the node.
    call Galacticus_Calculations_Reset(thisNode)
    call thisNode%odeStepScalesInitialize()
    !# <include directive="scaleSetTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'objects.tree_node.set_scale.inc'
    !# </include>
    call thisNode%serializeScales(propertyScales)

    ! Assign module global pointer to this node.
    activeTreeIndex=  thisTree%index
    activeNode     => thisNode

    ! Reset interrupt variables.
    firstInterruptFound     =  .false.
    firstInterruptTime      =  0.0d0
    firstInterruptProcedure => null()

    ! Call ODE solver routines.
    startTimeThisNode=basicComponent%time()
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
         &                   errorHandler=Galacticus_ODE_Error_Handler, &
         &                   algorithm   =Galacticus_ODE_Algorithm      &
         &                  )

    ! Extract values.
    call thisNode%deserializeValues(propertyValues)

    ! Ensure that the maximum time has not been exceed (can happen due to rounding errors).
    if (basicComponent%time() > endTime) call basicComponent%timeSet(endTime)

    ! Flag interruption if one occurred.
    if (firstInterruptTime /= 0.0d0) then
       interrupted=.true.
       interruptProcedure => firstInterruptProcedure
    else
       interrupted=.false.
    end if

    ! Call routines to perform any post-evolution tasks.
    !# <include directive="postEvolveTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'objects.tree_node.post_evolve.inc'
    !# </include>

    return
  end subroutine Tree_Node_Evolve

  logical function Tree_Node_Is_Accurate(valueNode,valueExpected)
    !% Return true if a tree node property is within expected accuracy of a given value.
    use Numerical_Comparison
    implicit none
    double precision, intent(in   ) :: valueExpected, valueNode

    ! Initialize.
    call Tree_Node_Evolve_Initialize()

    ! Compare values.
    Tree_Node_Is_Accurate=Values_Agree(valueNode,valueExpected,odeToleranceAbsolute,odeToleranceRelative)
    return
  end function Tree_Node_Is_Accurate

  function Tree_Node_ODEs(time,y,dydt,parameterPointer) bind(c)
    !% Function which evaluates the set of ODEs for the evolution of a specific node.
    use ODE_Solver_Error_Codes
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer  (kind=c_int                  )                       :: Tree_Node_ODEs
    real     (kind=c_double               )               , value :: time
    real     (kind=c_double               ), intent(in   )        :: y                 (nProperties)
    real     (kind=c_double               ), intent(  out)        :: dydt              (nProperties)
    type     (c_ptr                       )               , value :: parameterPointer
    logical                                                       :: interrupt
    procedure(Interrupt_Procedure_Template), pointer              :: interruptProcedure

    ! Extract values.
    call activeNode%deserializeValues(y)

    ! Set derivatives to zero initially.
    call activeNode%odeStepRatesInitialize()

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
          call activeNode%serializeRates(dydt)
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

  subroutine Tree_Node_Compute_Derivatives(thisNode,interrupt,interruptProcedureReturn)
    !% Call routines to set alls derivatives for {\tt thisNode}.
    use Galacticus_Calculations_Resets
    !# <include directive="preDerivativeTask" type="moduleUse">
    include 'objects.merger_trees.prederivative.tasks.modules.inc'
    !# </include>
    !# <include directive="rateComputeTask" type="moduleUse">
    include 'objects.node.component.derivatives.modules.inc'
    !# </include>
    implicit none
    type     (treeNode), intent(inout), pointer :: thisNode
    logical            , intent(  out)          :: interrupt
    procedure(        ), intent(  out), pointer :: interruptProcedureReturn
    procedure(        )               , pointer :: interruptProcedure

    ! Initialize interrupt status.
    interrupt=.false.
    interruptProcedure => null()

    ! Call component routines to indicate that derivative calculation is commencing.
    call Galacticus_Calculations_Reset(thisNode)

    ! Call routines to perform any pre-derivative calculations.
    !# <include directive="preDerivativeTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'objects.merger_trees.prederivative.tasks.inc'
    !# </include>

    ! Call component routines to compute derivatives.
    !# <include directive="rateComputeTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,interrupt,interruptProcedure</functionArgs>
    include 'objects.node.component.derivatives.inc'
    !# </include>

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure

    return
  end subroutine Tree_Node_Compute_Derivatives

  subroutine Tree_Node_ODEs_Error_Handler()
    !% Handles errors in the ODE solver when evolving \glc\ nodes. Dumps the content of the node.
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string) :: message

    message="ODE solver failed in tree #"
    message=message//activeTreeIndex
    call Galacticus_Display_Message(message)
    call activeNode%dump()
    return
  end subroutine Tree_Node_ODEs_Error_Handler

#ifdef PROFILE
  subroutine Tree_Node_Evolve_Error_Analyzer(currentPropertyValue,currentPropertyError,timeStep,stepStatus) bind(c)
    !% Profiles ODE solver step sizes and errors.
    !# <include directive="decodePropertyIdentifiersTask" type="moduleUse">
    include 'objects.merger_trees.decode_property_identifiers.modules.inc'
    !# </include>
    implicit none
    real            (kind=c_double ), dimension(nProperties), intent(in   )        :: currentPropertyValue
    real            (kind=c_double ), dimension(nProperties), intent(in   )        :: currentPropertyError
    real            (kind=c_double )                        , intent(in   ), value :: timeStep
    integer         (kind=c_int    )                        , intent(in   ), value :: stepStatus
    double precision                                                               :: scale               , scaledError     , scaledErrorMaximum
    integer                                                                        :: iProperty           , limitingProperty
    type            (varying_string)                                               :: propertyName

    ! If the step was not good, return immediately.
    if (stepStatus /= FGSL_Success) return

    ! Find the property with the largest error (i.e. that which is limiting the step).
    scaledErrorMaximum=0.0d0
    do iProperty=1,nProperties
       scale=odeToleranceAbsolute*propertyScales(iProperty)+odeToleranceRelative*abs(currentPropertyValue(iProperty))
       scaledError=abs(currentPropertyError(iProperty))/scale
       if (scaledError > scaledErrorMaximum) then
          scaledErrorMaximum=scaledError
          limitingProperty=iProperty
       end if
    end do
    ! Decode the step limiting property.
    propertyName=activeNode%nameFromIndex(limitingProperty)
    ! Record this information.
    call Galacticus_Meta_Evolver_Profile(timeStep,propertyName)
    return
  end subroutine Tree_Node_Evolve_Error_Analyzer
#endif

  subroutine Tree_Node_Promote(thisTree,thisNode)
    !% Transfer the properties of {\tt thisNode} to its parent node, then destroy it.
    use String_Handling
    use Galacticus_Display
    !# <include directive="nodePromotionTask" type="moduleUse">
    include 'objects.tree_node.promote.modules.inc'
    !# </include>
    implicit none
    class(mergerTree    ), intent(inout)          :: thisTree
    type (treeNode      ), intent(inout), pointer :: thisNode
    type (treeNode      )               , pointer :: parentNode, satelliteNode
    type (varying_string)                         :: message

    ! Get pointer to parent node.
    parentNode => thisNode%parent

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Promoting node '
       message=message//thisNode%index()//' to '//parentNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if

    ! Perform any processing necessary before this halo is promoted.
    !# <include directive="nodePromotionTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'objects.tree_node.promote.inc'
    !# </include>

    ! Move the components of thisNode to the parent.
    call thisNode%moveComponentsTo(parentNode)

    ! Copy any formation node data to the parent, and update the formation node's parentNode pointer to point to the new parent.
    if (associated(thisNode%formationNode)) then
       if (associated(parentNode%formationNode)) then
          call parentNode%formationNode%destroy()
          deallocate(parentNode%formationNode)
       end if
       allocate(parentNode%formationNode)
       call thisNode%formationNode%copyNodeTo(parentNode%formationNode)
       parentNode%formationNode%parent => parentNode
    end if

    ! Transfer any satellite nodes to the parent.
    if (associated(thisNode%firstSatellite)) then
       ! Attach the satellite nodes to the parent.
       if (associated(parentNode%firstSatellite)) then
          ! Find the last satellite of the parent node.
          satelliteNode         => parentNode%lastSatellite()
          satelliteNode%sibling => thisNode  %firstSatellite
       else
          parentNode%firstSatellite => thisNode%firstSatellite
       end if
       ! Get the first satellite of thisNode.
       satelliteNode => thisNode%firstSatellite
       do while (associated(satelliteNode))
          ! Set the parent node for this satellite to the parent.
          satelliteNode%parent => parentNode
          satelliteNode        => satelliteNode%sibling
       end do
    end if

    ! Nullify the child pointer for the parent.
    parentNode%firstChild => null()

    ! Destroy the node.
    call thisNode%destroy()
    deallocate(thisNode)

    return
  end subroutine Tree_Node_Promote

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
    class(mergerTree    ), intent(in   )          :: thisTree
    type (treeNode      ), intent(inout), pointer :: thisNode
    type (varying_string)                         :: message

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Making node '
       message=message//thisNode%index()//' a satellite in '//thisNode%parent%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if

    ! Call subroutines to perform any necessary processing prior to this node merger event.
    !# <include directive="nodeMergerTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
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
       !# <include directive="nodeMergersMethod" type="functionCall" functionType="void">
       !#  <functionArgs>nodeMergersMethod,Events_Node_Merger_Do</functionArgs>
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

end module Merger_Trees_Evolve_Node
