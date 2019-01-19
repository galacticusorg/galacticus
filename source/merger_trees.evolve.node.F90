!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements evolution of a single node in a merger tree.

module Merger_Trees_Evolve_Node
  !% Implements evolution of a single node in a merger tree.
  use Galacticus_Nodes  , only : treeNode, interruptTask
  use Kind_Numbers
  use ISO_Varying_String
  use FODEIV2
  implicit none
  private
  public :: Tree_Node_Evolve, Tree_Node_Promote, Events_Node_Merger, Tree_Node_Is_Accurate
  
  ! Variables used in the ODE solver.
  type            (fodeiv2_driver)                            :: ode2Driver
  type            (fodeiv2_system)                            :: ode2System
  logical                                                     :: odeReset
  !$omp threadprivate(ode2System,ode2Driver,odeReset)
  
  ! Flag to indicate if the node evolver method has been initialized.
  logical                                                     :: evolverInitialized         =.false.

  ! Parameters controlling the accuracy of ODE solving.
  double precision                                            :: odeToleranceAbsolute               , odeToleranceRelative           , &
       &                                                         odeJacobianStepSizeRelative

  ! Options for latent variable integrator.
  !# <enumeration>
  !#  <name>latentIntegratorType</name>
  !#  <description>Used to specify the type of latent variable integrator to use.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="gaussKronrod"/>
  !#  <entry label="trapezoidal" />
  !# </enumeration>
  integer                                                     :: odeLatentIntegratorType_           , odeLatentIntegratorOrder       , &
       &                                                         odeLatentIntegratorIntervalsMaximum
  
  ! Arrays that point to node properties and their derivatives.
  integer                                                     :: propertyCountAll                   , propertyCountMaximum         =0, &
       &                                                         propertyCountInactive              , propertyCountActive
  double precision                , allocatable, dimension(:) :: propertyScalesActive               , propertyValuesActive           , &
       &                                                         propertyErrors                     , propertyTolerances             , &
       &                                                         propertyValuesActiveSaved          , propertyScalesInactive         , &
       &                                                         propertyValuesInactiveSaved        , propertyValuesInactive         , &
       &                                                         odeTolerancesInactiveRelative      , odeTolerancesInactiveAbsolute
  !$omp threadprivate(propertyCountMaximum,propertyCountAll,propertyCountActive,propertyCountInactive,propertyValuesActive,propertyValuesActiveSaved,propertyValuesInactive,propertyValuesInactiveSaved,propertyScalesActive,propertyScalesInactive,propertyErrors,propertyTolerances,odeTolerancesInactiveRelative,odeTolerancesInactiveAbsolute)
  logical                                                     :: profileOdeEvolver

  ! Module global pointer to the node being processed.
  integer         (kind=kind_int8              )            :: activeTreeIndex
  type            (treeNode                    ), pointer   :: activeNode
  integer                                                   :: trialCount                     , propertyTypeODE, &
       &                                                       propertyTypeIntegrator
  integer                                       , parameter :: trialCountMaximum=8
  !$omp threadprivate(activeNode,activeTreeIndex,trialCount,propertyTypeODE,propertyTypeIntegrator)
  ! Variables to track interrupt events.
  logical                                                 :: interruptFirstFound
  double precision                                        :: timeInterruptFirst
  procedure       (interruptTask               ), pointer :: functionInterruptFirst
  !$omp threadprivate(interruptFirstFound,timeInterruptFirst,functionInterruptFirst)
  ! Flag to indicate if node merging event method has been initialized.
  logical                                                 :: nodeMergersInitialized  =.false.

  ! Name of node mergers method used.
  type            (varying_string              )          :: nodeMergersMethod

  ! The algorithm to use for ODE solving.
  type            (fodeiv2_step_type           )          :: Galacticus_ODE_Algorithm        , Galacticus_ODE_Algorithm_Non_Jacobian
  logical                                                 :: useJacobian

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure       (Node_Mergers_Template       ), pointer :: Events_Node_Merger_Do   =>null()
  abstract interface
     subroutine Node_Mergers_Template(node)
       import treeNode
       type(treeNode), intent(inout), pointer :: node
     end subroutine Node_Mergers_Template
  end interface

  ! Pointer to an error handler for failures in the ODE solver.
  procedure(), pointer :: Galacticus_ODE_Error_Handler => Tree_Node_ODEs_Error_Handler

  ! Previous state of ODE system - used to restore this state if the ODE evaluation
  ! functions are called in succession without change.
  double precision                            :: timePrevious          =-1.0d0
  double precision, allocatable, dimension(:) :: propertyValuesPrevious       , propertyRatesPrevious
  !$omp threadprivate(timePrevious,propertyValuesPrevious,propertyRatesPrevious)
  
contains

  subroutine Tree_Node_Evolve_Initialize
    !% Initializes the tree evolving routines by reading in parameters
    use Input_Parameters
    use Galacticus_Error
    type(varying_string) :: odeAlgorithm, odeAlgorithmNonJacobian, odeLatentIntegratorType

    ! Initialize if necessary.
    if (.not.evolverInitialized) then
       !$omp critical (Tree_Node_Evolve_Initialize)
       if (.not.evolverInitialized) then
          ! Get tolerance values for the ODE solver.
          !# <inputParameter>
          !#   <name>odeToleranceAbsolute</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.01d0</defaultValue>
          !#   <description>The absolute tolerance used in solving differential equations for node evolution.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>odeToleranceRelative</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>1.0d-2</defaultValue>
          !#   <description>The relative tolerance used in solving differential equations for node evolution.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>odeJacobianStepSizeRelative</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.01d0</defaultValue>
          !#   <description>The relative step size to use when perturbing properties for purposes of computing a finite difference approximation to the ODE system Jacobian.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>odeAlgorithm</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('Runge-Kutta-Cash-Karp')</defaultValue>
          !#   <description>The algorithm to use in the ODE solver.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          useJacobian=.false.
          select case (char(odeAlgorithm))
          case ('Runge-Kutta-Cash-Karp')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_RKCK
          case ('Runge-Kutta-Second-Order')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_RK2
          case ('Runge-Kutta')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_RK4
          case ('Runge-Kutta-Fehlberg')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_RKF45
          case ('Runge-Kutta-Prince-Dormand')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_RK8PD
          case ('multistepAdams')
             Galacticus_ODE_Algorithm=Fodeiv2_Step_msAdams
          case ('Bulirsch-Stoer')
             Galacticus_ODE_Algorithm=Fodeiv2_step_BSimp
             useJacobian             =.true.
          case ('BDF')
             Galacticus_ODE_Algorithm=Fodeiv2_step_MSBDFActive
             useJacobian             =.true.
          case default
             call Galacticus_Error_Report('odeAlgorithm is unrecognized'//{introspection:location})
          end select
          !# <inputParameter>
          !#   <name>odeAlgorithmNonJacobian</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('Runge-Kutta-Cash-Karp')</defaultValue>
          !#   <description>The algorithm to use in the ODE solver.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          select case (char(odeAlgorithmNonJacobian))
          case ('Runge-Kutta-Cash-Karp')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_RKCK
          case ('Runge-Kutta-Second-Order')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_RK2
          case ('Runge-Kutta')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_RK4
          case ('Runge-Kutta-Fehlberg')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_RKF45
          case ('Runge-Kutta-Prince-Dormand')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_RK8PD
          case ('multistepAdams')
             Galacticus_ODE_Algorithm_Non_Jacobian=Fodeiv2_Step_msAdams
          case default
             call Galacticus_Error_Report('odeAlgorithmNonJacobian is unrecognized'//{introspection:location})
          end select
          if (useJacobian) then
             !# <inputParameter>
             !#   <name>odeLatentIntegratorType</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>var_str('trapezoidal')</defaultValue>
             !#   <description>The type of integrator to use for latent variables.</description>
             !#   <source>globalParameters</source>
             !#   <type>string</type>
             !# </inputParameter>
             odeLatentIntegratorType_=enumerationLatentIntegratorTypeEncode(char(odeLatentIntegratorType),includesPrefix=.false.)
             select case (odeLatentIntegratorType_)
             case (latentIntegratorTypeGaussKronrod)
                !# <inputParameter>
                !#   <name>odeLatentIntegratorOrder</name>
                !#   <cardinality>1</cardinality>
                !#   <defaultValue>15</defaultValue>
                !#   <description>The order of the integrator for latent variables.</description>
                !#   <source>globalParameters</source>
                !#   <type>integer</type>
                !# </inputParameter>
             end select
             !# <inputParameter>
             !#   <name>odeLatentIntegratorIntervalsMaximum</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>1000</defaultValue>
             !#   <description>The maxium number of intervals allowed in the integrator for latent variables.</description>
             !#   <source>globalParameters</source>
             !#   <type>integer</type>
             !# </inputParameter>
          end if
          !# <inputParameter>
          !#   <name>profileOdeEvolver</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not to profile the ODE evolver.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Flag that the module is now initialized.
          evolverInitialized=.true.
       end if
       !$omp end critical (Tree_Node_Evolve_Initialize)
    end if
    return
  end subroutine Tree_Node_Evolve_Initialize

  subroutine Tree_Node_Evolve(thisTree,node,timeEnd,interrupted,functionInterrupt)
    !% Evolves {\normalfont \ttfamily node} to time {\normalfont \ttfamily timeEnd}, or until evolution is interrupted.
    use Galacticus_Nodes               , only : nodeComponentBasic  , mergerTree      , propertyTypeAll , propertyTypeActive, &
         &                                      propertyTypeInactive, rateComputeState, propertyTypeNone
    use ODEIV2_Solver
    use Memory_Management
    use Galacticus_Calculations_Resets
    use Galacticus_Display
    use Galacticus_Error
    use ODE_Solver_Error_Codes
    use Numerical_Integration2
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
    !# <include directive="inactiveSetTask" type="moduleUse">
    include 'objects.tree_node.set_inactive.modules.inc'
    !# </include>
    !# <include directive="analyticSolverTask" type="moduleUse">
    include 'objects.tree_node.analytic_solver_task.modules.inc'
    !# </include>
    implicit none
    class           (mergerTree                  )      , intent(inout)          :: thisTree
    type            (treeNode                    )      , intent(inout), pointer :: node
    double precision                                    , intent(in   )          :: timeEnd
    logical                                             , intent(  out)          :: interrupted
    procedure       (interruptTask               )      , intent(  out), pointer :: functionInterrupt
    class           (nodeComponentBasic          )                     , pointer :: basicComponent
    integer                                       , save                         :: propertyCountPrevious=-1
    !$omp threadprivate(propertyCountPrevious)
    class          (integratorMultiVectorized1D  ), allocatable                  :: integrator_
    logical                                                                      :: solvedAnalytically      , solvedNumerically, &
         &                                                                          jacobianSolver
    double precision                                                             :: timeStart               , stepSize         , &
         &                                                                          timeStartSaved
    integer                                                                      :: lengthMaximum           , i                , &
         &                                                                          odeStatus
    type            (varying_string              )                               :: message                 , line
    character       (len =12                     )                               :: label
    type            (c_funptr                    )                               :: Error_Analyzer          , postStep
    
    ! Initialize.
    call Tree_Node_Evolve_Initialize()

    ! Call routines to perform any pre-evolution tasks.
    !# <include directive="preEvolveTask" type="functionCall" functionType="void">
    !#  <functionArgs>node</functionArgs>
    include 'objects.tree_node.pre_evolve.inc'
    !# </include>

    ! Determine the end time for this node - either the specified end time, or the time associated with the parent node, whichever
    ! occurs first.
    basicComponent => node%basic()
    timeStart=basicComponent%time()
    timeStartSaved   =timeStart

    ! Ensure calculations are reset for this new step.
    call Galacticus_Calculations_Reset(node)
    
    ! Attempt to find analytic solutions.
    solvedAnalytically=.false.     
    !# <include directive="analyticSolverTask" type="functionCall" functionType="void">
    !#  <functionArgs>node,timeStart,timeEnd,solvedAnalytically</functionArgs>
    include 'objects.tree_node.analytic_solver_task.inc'
    !# </include>
    ! Check if an analytic solution was available - use numerical solution if not.
    if (solvedAnalytically) then
       ! An analytic solution was available. Record that no interrupt therefore occurred.
       interrupted=.false.
    else       
       ! Compute offsets into serialization arrays for rates and scales.
       call node%serializationOffsets(propertyTypeAll)
       ! Find number of all evolvable variables (active and inactive) for this node.
       propertyCountAll=node%serializeCount(propertyTypeAll)
       ! Allocate pointer arrays if necessary.
       if (propertyCountAll > propertyCountMaximum) then
          if (allocated(propertyValuesActive)) then
             call Memory_Usage_Record(sizeof(propertyValuesActive         ),addRemove=-1)
             deallocate(propertyValuesActive         )
             call Memory_Usage_Record(sizeof(propertyValuesActiveSaved    ),addRemove=-1)
             deallocate(propertyValuesActiveSaved    )
             call Memory_Usage_Record(sizeof(propertyValuesInactive       ),addRemove=-1)
             deallocate(propertyValuesInactive       )
             call Memory_Usage_Record(sizeof(propertyValuesInactiveSaved  ),addRemove=-1)
             deallocate(propertyValuesInactiveSaved  )
             call Memory_Usage_Record(sizeof(propertyScalesActive         ),addRemove=-1)
             deallocate(propertyScalesActive         )
             call Memory_Usage_Record(sizeof(propertyScalesInactive       ),addRemove=-1)
             deallocate(propertyScalesInactive       )
             call Memory_Usage_Record(sizeof(propertyValuesPrevious       ),addRemove=-1)
             deallocate(propertyValuesPrevious       )
             call Memory_Usage_Record(sizeof(propertyRatesPrevious        ),addRemove=-1)
             deallocate(propertyRatesPrevious        )
             call Memory_Usage_Record(sizeof(propertyErrors               ),addRemove=-1)
             deallocate(propertyErrors               )
             call Memory_Usage_Record(sizeof(propertyTolerances           ),addRemove=-1)
             deallocate(propertyTolerances           )
             call Memory_Usage_Record(sizeof(odeTolerancesInactiveRelative),addRemove=-1)
             deallocate(odeTolerancesInactiveRelative)
             call Memory_Usage_Record(sizeof(odeTolerancesInactiveAbsolute),addRemove=-1)
             deallocate(odeTolerancesInactiveAbsolute)
          end if
          allocate(propertyValuesActive         (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyValuesActive         ))
          allocate(propertyValuesActiveSaved  (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyValuesActiveSaved    ))
          allocate(propertyValuesInactive       (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyValuesInactive       ))
          allocate(propertyValuesInactiveSaved  (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyValuesInactiveSaved  ))
          allocate(propertyScalesActive         (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyScalesActive         ))
          allocate(propertyScalesInactive       (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyScalesInactive       ))
          allocate(propertyValuesPrevious       (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyValuesPrevious       ))
          allocate(propertyRatesPrevious        (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyRatesPrevious        ))
          allocate(propertyErrors               (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyErrors               ))
          allocate(propertyTolerances           (propertyCountAll))
          call Memory_Usage_Record(sizeof(propertyTolerances           ))
          allocate(odeTolerancesInactiveRelative(propertyCountAll))
          call Memory_Usage_Record(sizeof(odeTolerancesInactiveRelative))
          allocate(odeTolerancesInactiveAbsolute(propertyCountAll))
          call Memory_Usage_Record(sizeof(odeTolerancesInactiveAbsolute))
          propertyCountMaximum  =propertyCountAll
          propertyValuesPrevious=0.0d0
       end if
       ! Iterate until the step has been solved numerically.
       solvedNumerically=.false.
       jacobianSolver   =useJacobian
       do while (.not.solvedNumerically)
          if (jacobianSolver) then
             propertyTypeODE       =propertyTypeActive
             propertyTypeIntegrator=propertyTypeInactive
          else
             propertyTypeODE       =propertyTypeAll
             propertyTypeIntegrator=propertyTypeNone
          end if
          rateComputeState=propertyTypeODE
          ! Determine active and inactive properties.
          call node%odeStepInactivesInitialize()       
          if (jacobianSolver) then
             !# <include directive="inactiveSetTask" type="functionCall" functionType="void">
             !#  <functionArgs>node</functionArgs>
             include 'objects.tree_node.set_inactive.inc'
             !# </include>
          end if
          ! Compute offsets into serialization arrays for rates and scales.
          call node%serializationOffsets(propertyTypeODE       )
          call node%serializationOffsets(propertyTypeIntegrator)
          ! Find number of active and inactive variables to evolve for this node.
          propertyCountActive  =node%serializeCount(propertyTypeODE       )
          propertyCountInactive=node%serializeCount(propertyTypeIntegrator)
          ! Serialize property values to array, and save a copy in case of a need to reset.
          call node%serializeValues(propertyValuesActive  ,propertyTypeODE       )
          call node%serializeValues(propertyValuesInactive,propertyTypeIntegrator)
          propertyValuesActiveSaved  (1:propertyCountActive  )=propertyValuesActive  (1:propertyCountActive  )
          propertyValuesInactiveSaved(1:propertyCountInactive)=propertyValuesInactive(1:propertyCountInactive)
          ! Compute scales for all properties and extract from the node.
          call node%odeStepScalesInitialize()
          !# <include directive="scaleSetTask" type="functionCall" functionType="void">
          !#  <functionArgs>node</functionArgs>
          include 'objects.tree_node.set_scale.inc'
          !# </include>
          call node%serializeScales(propertyScalesActive  ,propertyTypeODE       )
          call node%serializeScales(propertyScalesInactive,propertyTypeIntegrator)
          ! Check for zero property scales which will cause floating point overflow in the ODE solver.
          if (any(propertyScalesActive(1:propertyCountActive) == 0.0d0)) then
             message='WARNING: Zero entry in ODE system scales for node'
             call Galacticus_Warn          (message)
             call Galacticus_Display_Indent(message)
             lengthMaximum=0    
             do i=1,propertyCountActive
                lengthMaximum=max(lengthMaximum,len(node%nameFromIndex(i,propertyTypeODE)))
             end do
             line=repeat("―",lengthMaximum)//repeat("―――――――――――――――",2)
             call Galacticus_Display_Message(line)
             call Galacticus_Display_Message(repeat(" ",lengthMaximum)//' : y            : yScale')
             call Galacticus_Display_Message(line)
             do i=1,propertyCountActive
                message=node%nameFromIndex(i,propertyTypeODE)
                message=repeat(" ",lengthMaximum-len(message))//message
                write (label,'(e12.6)') propertyValuesActive(i)
                message=message//" : "//label
                write (label,'(e12.6)') propertyScalesActive(i)
                message=message//" : "//label
                call Galacticus_Display_Message(message)
             end do
             call Galacticus_Display_Message(line)
             call node%serializeASCII()
             call Galacticus_Display_Unindent('done')          
          end if
          ! Assign module global pointer to this node.
          activeTreeIndex=  thisTree%index
          activeNode     => node
          ! Reset interrupt variables.
          interruptFirstFound     =  .false.
          timeInterruptFirst      =  0.0d0
          functionInterruptFirst => null()
          ! Call ODE solver routines.
          if (profileOdeEvolver) then
             Error_Analyzer=c_funloc(Tree_Node_Evolve_Error_Analyzer)
          else
             Error_Analyzer=C_NULL_FUNPTR
          end if
          postStep=c_funloc(Tree_Node_Post_Step)
          solvedNumerically=.true.
          if (timeStart /= timeEnd) then          
             if (propertyCountPrevious > 0 .and. .not.odeReset) call ODEIV2_Solver_Free(ode2Driver,ode2System)
             propertyCountPrevious=propertyCountActive
             odeStatus            =errorStatusFail
             trialCount           =0
             do while (trialCount < trialCountMaximum .and. .not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt))
                ! Initialize state. Stepsize is reduced by half on each successive trial.
                odeReset    =.true.
                timePrevious=-1.0d0
                stepSize    =node%timeStep()/2.0d0**trialCount
                timeStart   =timeStartSaved
                if (jacobianSolver) then
                   ! Initialize integrator.
                   if (.not.allocated(integrator_)) then
                      select case (odeLatentIntegratorType_)
                      case (latentIntegratorTypeGaussKronrod)
                         allocate(integratorMultiVectorizedCompositeGaussKronrod1D :: integrator_)
                         select type (integrator_)
                         type is (integratorMultiVectorizedCompositeGaussKronrod1D)
                            call integrator_%initialize(odeLatentIntegratorIntervalsMaximum,odeLatentIntegratorOrder)
                         end select
                      case (latentIntegratorTypeTrapezoidal )
                         allocate(integratorMultiVectorizedCompositeTrapezoidal1D  :: integrator_)
                         select type (integrator_)
                         type is (integratorMultiVectorizedCompositeTrapezoidal1D )
                            call integrator_%initialize(odeLatentIntegratorIntervalsMaximum                         )
                         end select
                      end select
                   end if
                   odeTolerancesInactiveRelative(1:propertyCountInactive)=odeToleranceRelative
                   odeTolerancesInactiveAbsolute(1:propertyCountInactive)=odeToleranceAbsolute*propertyScalesInactive(1:propertyCountInactive)
                   call integrator_%tolerancesSet(odeTolerancesInactiveAbsolute(1:propertyCountInactive),odeTolerancesInactiveRelative(1:propertyCountInactive))
                   call ODEIV2_Solve(                                                                         &
                        &            ode2Driver                                                             , &
                        &            ode2System                                                             , &
                        &            timeStart                                                              , &
                        &            timeEnd                                                                , &
                        &            propertyCountActive                                                    , &
                        &            propertyValuesActive                                                   , &
                        &            Tree_Node_ODEs                                                         , &
                        &            odeToleranceAbsolute                                                   , &
                        &            odeToleranceRelative                                                   , &
                        &            postStep                                                               , &
                        &            Error_Analyzer                                                         , &
                        &            propertyScalesActive                                                   , &
                        &            reset                  =odeReset                                       , &
                        &            errorHandler           =Galacticus_ODE_Error_Handler                   , &
                        &            odeStatus              =odeStatus                                      , &
                        &            algorithm              =Galacticus_ODE_Algorithm                       , &
                        &            stepSize               =stepSize                                       , &
                        &            jacobian               =Tree_Node_ODEs_Jacobian                        , &
                        &            integrator_            =integrator_                                    , &
                        &            integratorErrorTolerate=.true.                                         , &
                        &            zCount                 =propertyCountInactive                          , &
                        &            z                      =propertyValuesInactive(1:propertyCountInactive), &
                        &            integrands             =Tree_Node_Integrands                             &
                        &           )
                else
                   call ODEIV2_Solve(                                                                         &
                        &            ode2Driver                                                             , &
                        &            ode2System                                                             , &
                        &            timeStart                                                              , &
                        &            timeEnd                                                                , &
                        &            propertyCountActive                                                    , &
                        &            propertyValuesActive                                                   , &
                        &            Tree_Node_ODEs                                                         , &
                        &            odeToleranceAbsolute                                                   , &
                        &            odeToleranceRelative                                                   , &
                        &            postStep                                                               , &
                        &            Error_Analyzer                                                         , &
                        &            propertyScalesActive                                                   , &
                        &            reset                  =odeReset                                       , &
                        &            errorHandler           =Galacticus_ODE_Error_Handler                   , &
                        &            odeStatus              =odeStatus                                      , &
                        &            algorithm              =Galacticus_ODE_Algorithm_Non_Jacobian          , &
                        &            stepSize               =stepSize                                         &
                        &           )
                end if
                ! Check for failure.
                if (.not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt)) then
                   ! Increment number of trials.
                   trialCount=trialCount+1
                   ! Restore state of the node.
                   propertyValuesActive  (1:propertyCountActive  )=propertyValuesActiveSaved  (1:propertyCountActive  )
                   propertyValuesInactive(1:propertyCountInactive)=propertyValuesInactiveSaved(1:propertyCountInactive)
                end if
             end do
             if (.not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt)) then
                if (jacobianSolver) then
                   ! Restore state of the node, switch to using a non-jacobian solver, and indicate that the step was not solved numerically.
                   solvedNumerically=.false.
                   jacobianSolver   =.false.
                   propertyValuesActive  (1:propertyCountActive  )=propertyValuesActiveSaved  (1:propertyCountActive  )
                   propertyValuesInactive(1:propertyCountInactive)=propertyValuesInactiveSaved(1:propertyCountInactive)
                else
                   call Galacticus_Error_Report('ODE integration failed '//{introspection:location})
                end if
             end if
          end if
          if (solvedNumerically) then
             call node%timeStepSet(stepSize)
             ! Extract values.
             call node%deserializeValues(propertyValuesActive  ,propertyTypeODE       )
             call node%deserializeValues(propertyValuesInactive,propertyTypeIntegrator)
             ! Flag interruption if one occurred, and ensure that the time is matched precisely to the end or interrupt time (can differ
             ! due to finite precision of the ODE integrator).
             if (timeInterruptFirst /= 0.0d0) then
                interrupted=.true.
                functionInterrupt => functionInterruptFirst
                call basicComponent%timeSet(timeInterruptFirst)          
             else
                interrupted=.false.
                call basicComponent%timeSet(           timeEnd)          
             end if
          end if
       end do
    endif
    ! Call routines to perform any post-evolution tasks.
    if (associated(node)) then
       !# <include directive="postEvolveTask" type="functionCall" functionType="void">
       !#  <functionArgs>node</functionArgs>
       include 'objects.tree_node.post_evolve.inc'
       !# </include>
    end if
    
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

  subroutine Tree_Node_Integrands(propertyCountActive,propertyCountInactive,time,propertyValues,evaluate,integrands)
    !% A set of integrands for unit tests.
    use Galacticus_Nodes, only : rateComputeState
    implicit none
    integer                        , intent(in   )                                              :: propertyCountActive       , propertyCountInactive
    double precision               , intent(in   ), dimension(                              : ) :: time
    double precision               , intent(in   ), dimension(propertyCountActive  ,size(time)) :: propertyValues
    logical                        , intent(inout), dimension(                              : ) :: evaluate
    double precision               , intent(  out), dimension(propertyCountInactive,size(time)) :: integrands
    logical                        , parameter                                                  :: odeConverged       =.true.
    procedure       (interruptTask), pointer                                                    :: functionInterrupt
    logical                                                                                     :: interrupt
    integer                                                                                     :: iTime
    ! "evaluate" array is currently not used. It indicates which integrands must be evaluated, and which can (optionally) be
    ! ignored as they have already converged to the required tolerance. It is currently not used because the potential for
    ! significant speed up appears to be small based on profiling. This will be model-depdendent though, so this decision can be
    ! revisited.    
    !GCC$ attributes unused :: evaluate
    
    ! Set state to indicate that rate calculations are for inactive variables.
    rateComputeState=propertyTypeIntegrator
    ! Iterate over times.
    do iTime=1,size(time)
       ! Deserialize the active values.
       call activeNode%deserializeValues(propertyValues(1:propertyCountActive,iTime),propertyTypeODE)
       ! If past the time of the first interrupt, integrands are set to zero. Otherwise, evaluate integrands.
       if (interruptFirstFound .and. time(iTime) >= timeInterruptFirst) then
          integrands(:,iTime)=0.0d0
       else
          ! Set derivatives to zero initially.
          call activeNode%odeStepRatesInitialize()
          ! Compute derivatives.
          call Tree_Node_Compute_Derivatives(activeNode,odeConverged,interrupt,functionInterrupt,propertyTypeIntegrator)      
          ! Serialize rates into integrand array.
          call activeNode%serializeRates(integrands(:,iTime),propertyTypeIntegrator)
       end if
    end do
    ! Set state to indicate that rate calculations are for active variables.
    rateComputeState=propertyTypeODE
    return
  end subroutine Tree_Node_Integrands

  integer function Tree_Node_ODEs(time,y,dydt)
    !% Function which evaluates the set of ODEs for the evolution of a specific node.
    use ODE_Solver_Error_Codes
    use Galacticus_Nodes      , only : nodeComponentBasic
    use FGSL                  , only : FGSL_Success
    implicit none
    double precision                     , intent(in   )               :: time
    double precision                     , intent(in   ), dimension(:) :: y
    double precision                     , intent(  out), dimension(:) :: dydt
    logical                                                            :: interrupt         , odeConverged
    procedure       (interruptTask     ), pointer                      :: functionInterrupt
    class           (nodeComponentBasic), pointer                      :: basic

    
    ! Return success by default.
    Tree_Node_ODEs=FGSL_Success
    ! Check if we can reuse the previous derivatives.
    if   (                                                                                              &
       &   time                                       == timePrevious                                   &
       &  .and.                                                                                         &
       &   all(y(             1:propertyCountActive)  == propertyValuesPrevious(1:propertyCountActive)) &
       & ) then
       dydt                  (1:propertyCountActive)  =  propertyRatesPrevious (1:propertyCountActive)
       return
    else
       timePrevious                                   =  time
       propertyValuesPrevious(1:propertyCountActive)  =  y                     (1:propertyCountActive)
    end if
    ! Extract values.
    call activeNode%deserializeValues(y(1:propertyCountActive),propertyTypeODE)
    ! If the node is significantly inaccurate (as judged by the node time being different from the system time), then set rates to
    ! zero, as the ODE solver is likely just taking a step which is too large.
    basic => activeNode%basic()
    if (.not.Tree_Node_Is_Accurate(basic%time(),time)) then
       dydt=0.0d0
       return
    end if
    ! Set derivatives to zero initially.
    call activeNode%odeStepRatesInitialize()
    ! Determine if the ODE evolver has reached sufficiently small errors for this step.
    call FODEIV2_Driver_Errors(ode2Driver,propertyErrors(1:propertyCountActive))
    propertyTolerances(1:propertyCountActive)=treeNodeODEStepTolerances(y)
    odeConverged=all(propertyErrors(1:propertyCountActive) <= propertyTolerances(1:propertyCountActive))
    if (interruptFirstFound .and. time >= timeInterruptFirst) then
       ! Already beyond the location of the first interrupt, simply return zero derivatives.
       dydt                 (1:propertyCountActive)=0.0d0
       propertyRatesPrevious(1:propertyCountActive)=0.0d0
    else
       ! Compute derivatives.
       call Tree_Node_Compute_Derivatives(activeNode,odeConverged,interrupt,functionInterrupt,propertyTypeODE)      
       ! Check whether an interrupt has been requested.
       select case (interrupt)
       case (.false.)
          ! No interrupt - place derivatives into ODE arrays.
          call activeNode%serializeRates(dydt(1:propertyCountActive),propertyTypeODE)
          propertyRatesPrevious(1:propertyCountActive)=dydt(1:propertyCountActive)
       case (.true.)
          ! Interrupt requested - freeze evolution and store the interrupt if it is the earliest one to occur.
          dydt                 (1:propertyCountActive)=0.0d0
          propertyRatesPrevious(1:propertyCountActive)=0.0d0
          if (time < timeInterruptFirst .or. .not.interruptFirstFound) then
             interruptFirstFound    =  .true.
             timeInterruptFirst     =  time
             functionInterruptFirst => functionInterrupt
             ! Let the ODE solver know that an interrupt occured, and when it happened.
             Tree_Node_ODEs         =  odeSolverInterrupt
             interruptedAtX         =  time
             return
          end if
       end select
    end if
    return
  end function Tree_Node_ODEs

  integer function Tree_Node_ODEs_Jacobian(time,propertyValues0,derivativeRatesValues,derivativeRatesTime)
    !% Function which evaluates the set of ODEs for the evolution of a specific node.
    use ODE_Solver_Error_Codes
    use Galactic_Structure_Radii
    use Numerical_Comparison
    use Galacticus_Error
    use FGSL                    , only : FGSL_Success
    implicit none
    double precision                                                                   , intent(in   ) :: time
    double precision               , dimension(:                                      ), intent(in   ) :: propertyValues0
    double precision               , dimension(:                                      ), intent(  out) :: derivativeRatesValues        , derivativeRatesTime
    double precision               , dimension(propertyCountActive                    )                :: propertyRates0               , propertyRates1     , &
         &                                                                                                propertyValues1
    double precision               , dimension(propertyCountActive,propertyCountActive)                :: jacobian
    procedure       (interruptTask), pointer                                                           :: functionInterrupt
    double precision               , parameter                                                         :: deltaTiny            =1.0d-10
    logical                                                                                            :: interrupt                    , odeConverged
    integer                                                                                            :: i
    double precision                                                                                   :: propertyValueDelta
    
    ! Return success by default.
    Tree_Node_ODEs_Jacobian=FGSL_Success
    ! No explicit time dependence.
    derivativeRatesTime=0.0d0
    ! Check for interrupts.
    if (interruptFirstFound .and. time >= timeInterruptFirst) then
       ! Already beyond the location of the first interrupt, simply return zero derivatives.
       jacobian(1:propertyCountActive,1:propertyCountActive)=0.0d0
    else
       ! Compute rates at current parameter values.
       call activeNode%deserializeValues     (propertyValues0(1:propertyCountActive)             ,propertyTypeODE)
       call activeNode%odeStepRatesInitialize(                                                                   )
       call Tree_Node_Compute_Derivatives    (activeNode,odeConverged,interrupt,functionInterrupt,propertyTypeODE)
       call activeNode%serializeRates        (propertyRates0                                     ,propertyTypeODE)
       ! If an interrupt was triggered, then derivatives will all be zero, so we set the Jacobian to zero here and exit.
       if (interrupt) then
          jacobian(1:propertyCountActive,1:propertyCountActive)=0.0d0
       else
          ! Iterate over parameters, computing Jacobian using finite differences.
          do i=1,propertyCountActive
             ! To compute the finite difference we make a small perturbation in one parameter. If the parameter is non-zero, use a
             ! small, fractional perturbation. For parameters with zero value, use a perturbation equal to the absolute tolerance
             ! supplied to the ODE solver.
             if (propertyValues0(i)==0.0d0) then
                propertyValueDelta       =+propertyScalesActive       (i)
             else
                propertyValueDelta       =+odeJacobianStepSizeRelative    &
                     &                    *propertyValues0            (i)
                if (abs(propertyValueDelta) < deltaTiny*propertyScalesActive(i)) &
                     & propertyValueDelta=+propertyScalesActive       (i)
             end if
             propertyValues1       =+propertyValues0
             propertyValues1(i)    =+propertyValues1            (i) &
                  &                 +propertyValueDelta
             call activeNode%deserializeValues     (propertyValues1,propertyTypeODE                                         )
             call activeNode%odeStepRatesInitialize(                                                                        )
             call Galactic_Structure_Radii_Revert  (activeNode                                                              )
             call Tree_Node_Compute_Derivatives    (activeNode     ,odeConverged,interrupt,functionInterrupt,propertyTypeODE)
             call activeNode%serializeRates        (propertyRates1                                          ,propertyTypeODE)
             jacobian(i,:)=+(                  &
                  &          +propertyRates1   &
                  &          -propertyRates0   &
                  &         )                  &
                  &        /propertyValueDelta
          end do
       end if
    end if
    ! Map Jacobian back to output array.
    derivativeRatesValues=reshape(jacobian,[propertyCountActive**2])
    return
  end function Tree_Node_ODEs_Jacobian
  
  subroutine Tree_Node_Compute_Derivatives(node,odeConverged,interrupt,functionInterruptReturn,propertyType)
    !% Call routines to set alls derivatives for {\normalfont \ttfamily node}.
    use Galacticus_Calculations_Resets
    !# <include directive="preDerivativeTask" type="moduleUse">
    include 'objects.merger_trees.prederivative.tasks.modules.inc'
    !# </include>
    !# <include directive="rateComputeTask" type="moduleUse">
    include 'objects.node.component.derivatives.modules.inc'
    !# </include>
    implicit none
    type     (treeNode), intent(inout), pointer :: node
    logical            , intent(in   )          :: odeConverged
    logical            , intent(  out)          :: interrupt
    procedure(        ), intent(  out), pointer :: functionInterruptReturn
    integer            , intent(in   )          :: propertyType
    procedure(        )               , pointer :: functionInterrupt

    ! Initialize interrupt status.
    interrupt         =  .false.
    functionInterrupt => null()
    ! Call component routines to indicate that derivative calculation is commencing.
    call Galacticus_Calculations_Reset(node)
    ! Call routines to perform any pre-derivative calculations.
    !# <include directive="preDerivativeTask" type="functionCall" functionType="void">
    !#  <functionArgs>node</functionArgs>
    include 'objects.merger_trees.prederivative.tasks.inc'
    !# </include>
    ! Do not attempt to compute derivatives for nodes which are not solvable.    
    if (.not.node%isSolvable) return   
    ! Call component routines to compute derivatives.
    !# <include directive="rateComputeTask" type="functionCall" functionType="void">
    !#  <functionArgs>node,odeConverged,interrupt,functionInterrupt,propertyType</functionArgs>
    include 'objects.node.component.derivatives.inc'
    !# </include>
    ! Return the procedure pointer.
    functionInterruptReturn => functionInterrupt
    return
  end subroutine Tree_Node_Compute_Derivatives

  subroutine Tree_Node_ODEs_Error_Handler(status,time,y)
    !% Handles errors in the ODE solver when evolving \glc\ nodes. Dumps the content of the node.
    use, intrinsic :: ISO_C_Binding
    use String_Handling
    use Galacticus_Display
    implicit none
    integer         (kind=c_int    ), intent(in   )                                 :: status
    real            (kind=c_double ), intent(in   )                                 :: time
    real            (kind=c_double ), intent(in   ), dimension(propertyCountActive) :: y
    real            (kind=c_double )               , dimension(propertyCountActive) :: dydt      , yError       , &
         &                                                                             yTolerance
    type            (varying_string)                                                :: message   , line
    integer                                                                         :: i         , lengthMaximum
    character       (len =12       )                                                :: label
    integer         (kind=c_int    )                                                :: odeStatus
    double precision                                                                :: stepFactor

    ! Check if this is the final trial for this node.
    if (trialCount == trialCountMaximum-1) then
       ! Get the current errors and tolerances in the ODE driver.
       call FODEIV2_Driver_Errors(ode2Driver,yError)
       yTolerance=treeNodeODEStepTolerances(y)
       ! Report the failure message.
       if (Galacticus_Verbosity_Level() < verbosityStandard) call Galacticus_Verbosity_Level_Set(verbosityStandard)
       message="ODE solver failed with error code "
       message=message//status//" in tree #"//activeTreeIndex
       call Galacticus_Display_Message(message)
       ! Dump all node properties.
       call activeNode%serializeASCII()
       ! Evaluate derivatives.
       odeStatus=Tree_Node_ODEs(time,y,dydt)
       call Galacticus_Display_Indent('ODE system parameters')
       lengthMaximum=0    
       do i=1,propertyCountActive
          lengthMaximum=max(lengthMaximum,len(activeNode%nameFromIndex(i,propertyTypeODE)))
       end do
       line=repeat("―",lengthMaximum)//repeat("―――――――――――――――",5)
       call Galacticus_Display_Message(line)
       call Galacticus_Display_Message(repeat(" ",lengthMaximum)//' : y            : dy/dt        : yScale       : yError       : yErrorScaled')
       call Galacticus_Display_Message(line)
       do i=1,propertyCountActive
          stepFactor=abs(yError(i))/yTolerance(i)
          message=activeNode%nameFromIndex(i,propertyTypeODE)
          message=repeat(" ",lengthMaximum-len(message))//message
          write (label,'(e12.6)') y                   (i)
          message=message//" : "//label
          write (label,'(e12.6)') dydt                (i)
          message=message//" : "//label
          write (label,'(e12.6)') propertyScalesActive(i)
          message=message//" : "//label
          write (label,'(e12.6)') yError              (i)
          message=message//" : "//label
          write (label,'(e12.6)') stepFactor
          message=message//" : "//label
          call Galacticus_Display_Message(message)
       end do
       call Galacticus_Display_Message(line)
       call Galacticus_Display_Unindent('done')
    end if
    return
  end subroutine Tree_Node_ODEs_Error_Handler

  function treeNodeODEStepTolerances(propertyValues)
    !% Compute the tolerances on each property being evolved in the ODE stystem at the current timestep.
    implicit none
    double precision               , dimension(propertyCountActive) :: treeNodeODEStepTolerances
    double precision, intent(in   ), dimension(propertyCountActive) :: propertyValues
    integer                                                         :: i
    
    forall(i=1:propertyCountActive)
       treeNodeODEStepTolerances(i)=+odeToleranceRelative         &
            &                       *abs(propertyValues      (i)) &
            &                       +odeToleranceAbsolute         &
            &                       *    propertyScalesActive(i)
    end forall
    return
  end function treeNodeODEStepTolerances
  
  subroutine Tree_Node_Post_Step(y,status) bind(c)
    !% Perform any post-step actions on the node.
    use, intrinsic :: ISO_C_Binding
    use               FGSL         , only : FGSL_Success
    !# <include directive="postStepTask" type="moduleUse">
    include 'objects.tree_node.post_step.modules.inc'
    !# </include>
    implicit none
    real   (kind=c_double), intent(inout), dimension(propertyCountActive) :: y
    integer(kind=c_int   ), intent(inout)                                 :: status

    call activeNode%deserializeValues(y,propertyTypeODE)
    !# <include directive="postStepTask" type="functionCall" functionType="void">
    !#  <functionArgs>activeNode,status</functionArgs>
    include 'objects.tree_node.post_step.inc'
    !# </include>
    if (status /= FGSL_Success) call activeNode%serializeValues(y,propertyTypeODE)
    return
  end subroutine Tree_Node_Post_Step
  
  subroutine Tree_Node_Evolve_Error_Analyzer(currentPropertyValue,currentPropertyError,timeStep,stepStatus) bind(c)
    !% Profiles ODE solver step sizes and errors.
    use, intrinsic :: ISO_C_Binding
    use               FGSL                            , only : FGSL_Success
    use               Galacticus_Meta_Evolver_Profiler
    implicit none
    real            (kind=c_double ), dimension(propertyCountActive), intent(in   )        :: currentPropertyValue
    real            (kind=c_double ), dimension(propertyCountActive), intent(in   )        :: currentPropertyError
    real            (kind=c_double )                                , intent(in   ), value :: timeStep
    integer         (kind=c_int    )                                , intent(in   ), value :: stepStatus
    double precision                                                                       :: scale               , scaledError     , scaledErrorMaximum
    integer                                                                                :: iProperty           , limitingProperty
    type            (varying_string)                                                       :: propertyName
   
    ! If the step was not good, return immediately.
    if (stepStatus /= FGSL_Success) return    
    ! Find the property with the largest error (i.e. that which is limiting the step).
    scaledErrorMaximum=0.0d0
    do iProperty=1,propertyCountActive
       scale=odeToleranceAbsolute*propertyScalesActive(iProperty)+odeToleranceRelative*abs(currentPropertyValue(iProperty))
       scaledError=abs(currentPropertyError(iProperty))/scale
       if (scaledError > scaledErrorMaximum) then
          scaledErrorMaximum=scaledError
          limitingProperty=iProperty
       end if
    end do
    ! Check that we found a limiting property.
    if (scaledErrorMaximum > 0.0d0) then
       ! Decode the step limiting property.
       propertyName=activeNode%nameFromIndex(limitingProperty,propertyTypeODE)
       ! Record this information.
       call Galacticus_Meta_Evolver_Profile(timeStep,propertyName)
    end if
    return
  end subroutine Tree_Node_Evolve_Error_Analyzer
  
  subroutine Tree_Node_Promote(node)
    !% Transfer the properties of {\normalfont \ttfamily node} to its parent node, then destroy it.
    use String_Handling
    use Galacticus_Display
    !# <include directive="nodePromotionTask" type="moduleUse">
    include 'objects.tree_node.promote.modules.inc'
    !# </include>
    implicit none
    type (treeNode      ), intent(inout), pointer :: node
    type (treeNode      )               , pointer :: parentNode, satelliteNode, &
         &                                           mergeeNode, hostNode
    type (varying_string)                         :: message
    
    ! Get pointer to parent node.
    parentNode => node%parent

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Promoting node '
       message=message//node%index()//' to '//parentNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if

    ! Perform any processing necessary before this halo is promoted.
    !# <include directive="nodePromotionTask" type="functionCall" functionType="void">
    !#  <functionArgs>node</functionArgs>
    include 'objects.tree_node.promote.inc'
    !# </include>
    ! Copy timestep to the parent.
    call parentNode%timeStepSet(node%timeStep())
    ! Move the components of node to the parent.
    call node%moveComponentsTo(parentNode)

    ! Copy any formation node data to the parent, and update the formation node's parentNode pointer to point to the new parent.
    if (associated(node%formationNode)) then
       if (associated(parentNode%formationNode)) then
          call parentNode%formationNode%destroy()
          deallocate(parentNode%formationNode)
       end if
       allocate(parentNode%formationNode)
       call node%formationNode%copyNodeTo(parentNode%formationNode)
       parentNode%formationNode%parent => parentNode
    end if
    ! Transfer any satellite nodes to the parent.
    if (associated(node%firstSatellite)) then
       hostNode => parentNode
       do while (hostNode%isSatellite())
          hostNode => hostNode%parent
       end do   
       ! Attach the satellite nodes to the parent.
       if (associated(hostNode%firstSatellite)) then
          ! Find the last satellite of the parent node.
          satelliteNode                 => hostNode% lastSatellite()
          satelliteNode%sibling         => node%firstSatellite
       else
          hostNode      %firstSatellite => node%firstSatellite
       end if
       ! Get the first satellite of node.
       satelliteNode => node%firstSatellite
       do while (associated(satelliteNode))
          ! Set the parent node for this satellite to the parent.
          satelliteNode%parent => hostNode
          satelliteNode        => satelliteNode%sibling
       end do
    end if

    ! Mergees of the node to be promoted must have their merge targets reset to the parent node.
    mergeeNode => node%firstMergee
    do while (associated(mergeeNode))
       mergeeNode%mergeTarget => parentNode
       mergeeNode => mergeeNode%siblingMergee
    end do
    if (associated(parentNode%firstMergee)) then
       mergeeNode => parentNode%firstMergee
       do while (associated(mergeeNode%siblingMergee))
          mergeeNode => mergeeNode%siblingMergee
       end do
       mergeeNode%siblingMergee => node%firstMergee
    else
       parentNode%firstMergee => node%firstMergee
    end if
    
    ! Nullify the child pointer for the parent.
    parentNode%firstChild => null()

    ! Destroy the node.
    call node%destroy()
    deallocate(node)

    return
  end subroutine Tree_Node_Promote

  subroutine Events_Node_Merger(node)
    !% Handles instances where {\normalfont \ttfamily node} is about to merge with its parent node.
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
    type (treeNode      ), intent(inout), pointer :: node
    type (varying_string)                         :: message

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Making node '
       message=message//node%index()//' a satellite in '//node%parent%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if    
    ! Call subroutines to perform any necessary processing prior to this node merger event.
    !# <include directive="nodeMergerTask" type="functionCall" functionType="void">
    !#  <functionArgs>node</functionArgs>
    include 'events.node_mergers.process.inc'
    !# </include>
    if (.not.nodeMergersInitialized) then
       !$omp critical (Events_Node_Merger_Initialize)
       if (.not.nodeMergersInitialized) then
          ! Get the node mergers method parameter.
          !# <inputParameter>
          !#   <name>nodeMergersMethod</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('singleLevelHierarchy')</defaultValue>
          !#   <description>Selects the method to be used for handling node merger events.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !# </inputParameter>
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="nodeMergersMethod" type="functionCall" functionType="void">
          !#  <functionArgs>nodeMergersMethod,Events_Node_Merger_Do</functionArgs>
          include 'events.node_mergers.inc'
          !# </include>
          if (.not.associated(Events_Node_Merger_Do)) call Galacticus_Error_Report('method '//char(nodeMergersMethod)//' is unrecognized'//{introspection:location})
          nodeMergersInitialized=.true.
       end if
       !$omp end critical (Events_Node_Merger_Initialize)
    end if
    ! Call the routine to perform the merger.
    call Events_Node_Merger_Do(node)
    return
  end subroutine Events_Node_Merger

end module Merger_Trees_Evolve_Node
