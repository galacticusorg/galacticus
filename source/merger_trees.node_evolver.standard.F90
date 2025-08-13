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
  Implements the standard class for evolving nodes in merger trees.
  !!}

  use :: Kind_Numbers                , only : kind_int8
  use :: Merger_Tree_Evolve_Profilers, only : mergerTreeEvolveProfilerClass
  use :: Merger_Trees_Merge_Node     , only : mergerTreeNodeMerger         , mergerTreeNodeMergerClass
  use :: Nodes_Operators             , only : nodeOperatorClass
  use :: Numerical_ODE_Solvers       , only : odeSolver
  use :: Timers                      , only : timer
  
  ! Enumeration of latent variable integrator.
  !![
  <enumeration>
   <name>latentIntegratorType</name>
   <description>Used to specify the type of latent variable integrator to use.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="gaussKronrod"/>
   <entry label="trapezoidal" />
  </enumeration>
  !!]

  ! Enumeration of ODE algorithms.
  !![
  <enumeration>
   <name>standardODEAlgorithm</name>
   <description>Used to specify the type of ODE algorithm to use.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="rungeKuttaCashKarp"     />
   <entry label="rungeKuttaSecondOrder"  />
   <entry label="rungeKutta"             />
   <entry label="rungeKuttaFehlberg"     />
   <entry label="rungeKuttaPrinceDormand"/>
   <entry label="multistepAdams"         />
   <entry label="bulirschStoer"          />
   <entry label="bdf"                    />
  </enumeration>
  !!]

  !![
  <mergerTreeNodeEvolver name="mergerTreeNodeEvolverStandard">
    <description>
      The standard merger tree node evolver.

      If the parameter {\normalfont \ttfamily [enforceNonNegativity] = true} then properties which are marked as being
      non-negative (e.g. masses) are evolved in such a way to ensure that they remain non-negative. This typically requires
      smaller time step size and so longer run times. In some cases it may be impossible to ensure non-negativity even for
      arbitrarily small timesteps\footnote{This can occur if a property as a non-zero, negative derivative as the property
      approaches zero. Such cases are quite likely unphysical, but are tolerated here.}. In such cases, if a property remains
      negative with the smallest possible time step, it will be zeroed and evolution continues.
    </description>
   <deepCopy>
    <ignore variables="galacticStructureSolver_"/>
   </deepCopy>
   <stateStorable>
    <exclude variables="galacticStructureSolver_"/>
   </stateStorable>
  </mergerTreeNodeEvolver>
  !!]
  type, extends(mergerTreeNodeEvolverClass) :: mergerTreeNodeEvolverStandard
     !!{
     Implementation of the standard merger tree node evolver.
     !!}
     private
     class           (mergerTreeNodeMergerClass          ), pointer                   :: mergerTreeNodeMerger_               => null()
     double precision                                                                 :: odeToleranceAbsolute                         , odeToleranceRelative         , &
          &                                                                              odeJacobianStepSizeRelative
     integer                                                                          :: odeAlgorithm                                 , odeAlgorithmNonJacobian
     class           (galacticStructureSolverClass       ), pointer                   :: galacticStructureSolver_            => null()
     class           (nodeOperatorClass                  ), pointer                   :: nodeOperator_                       => null()
     class           (mergerTreeEvolveProfilerClass      ), pointer                   :: mergerTreeEvolveProfiler_           => null()
     type            (enumerationLatentIntegratorTypeType)                            :: odeLatentIntegratorType
     integer                                                                          :: odeLatentIntegratorIntervalsMaximum          , odeLatentIntegratorOrder
     integer         (c_size_t                           )                            :: propertyCountAll                             , propertyCountMaximum         , &
          &                                                                              propertyCountInactive                        , propertyCountActive          , &
          &                                                                              propertyCountPrevious                        , countEvaluationsToSuccess
     double precision                                     , allocatable, dimension(:) :: propertyScalesActive                         , propertyValuesActive         , &
          &                                                                              propertyErrors                               , propertyTolerances           , &
          &                                                                              propertyValuesActiveSaved                    , propertyScalesInactive       , &
          &                                                                              propertyValuesInactiveSaved                  , propertyValuesInactive       , &
          &                                                                              odeTolerancesInactiveRelative                , odeTolerancesInactiveAbsolute
     logical                                              , allocatable, dimension(:) :: isNonNegative
     logical                                                                          :: profileOdeEvolver                            , reuseODEStepSize             , &
          &                                                                              enforceNonNegativity
     integer         (kind=kind_int8                     )                            :: activeTreeIndex
     type            (treeNode                           ), pointer                   :: activeNode                          => null()
     integer                                                                          :: trialCount                                   , propertyTypeODE              , &
          &                                                                              propertyTypeIntegrator
     logical                                                                          :: interruptFirstFound
     double precision                                                                 :: timeInterruptFirst
     procedure       (interruptTask                      ), pointer    , nopass       :: functionInterruptFirst
     logical                                                                          :: useJacobian
     double precision                                                                 :: timePrevious
     double precision                                     , allocatable, dimension(:) :: propertyValuesPrevious                       , propertyRatesPrevious
     integer         (kind_int8                          )                            :: systemClockMaximum
     type            (timer                              )                            :: stepTimer
   contains
     final     ::               standardDestructor
     procedure :: evolve     => standardEvolve
     procedure :: promote    => standardPromote
     procedure :: merge      => standardMerge
     procedure :: isAccurate => standardIsAccurate
     procedure :: autoHook   => standardAutoHook
  end type mergerTreeNodeEvolverStandard

  interface mergerTreeNodeEvolverStandard
     !!{
     Constructors for the \refClass{mergerTreeNodeEvolverStandard} merger tree node evolver.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeNodeEvolverStandard

  ! Maximum number of trials to solve ODEs.
  integer                                        , parameter :: trialCountMaximum=8

  ! Pointer to self and solver.
  class           (mergerTreeNodeEvolverStandard), pointer   :: self_
  type            (odeSolver                    ), pointer   :: solver_
  double precision                                           :: timeStartSaved
  !$omp threadprivate(self_,solver_,timeStartSaved)
  !$GLC ignore outlive :: solver_
  
contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeNodeEvolverStandard} merger tree node evolver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeNodeEvolverStandard)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (mergerTreeNodeMergerClass    ), pointer       :: mergerTreeNodeMerger_
    class           (nodeOperatorClass            ), pointer       :: nodeOperator_
    class           (mergerTreeEvolveProfilerClass), pointer       :: mergerTreeEvolveProfiler_
    type            (varying_string               )                :: odeAlgorithm               , odeAlgorithmNonJacobian            , &
         &                                                            odeLatentIntegratorType
    integer                                                        :: odeLatentIntegratorOrder   , odeLatentIntegratorIntervalsMaximum
    double precision                                               :: odeToleranceAbsolute       , odeToleranceRelative               , &
         &                                                            odeJacobianStepSizeRelative
    logical                                                        :: profileOdeEvolver          , reuseODEStepSize                   , &
         &                                                            enforceNonNegativity

    !![
    <inputParameter>
      <name>odeToleranceAbsolute</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The absolute tolerance used in solving differential equations for node evolution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeToleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>The relative tolerance used in solving differential equations for node evolution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeJacobianStepSizeRelative</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The relative step size to use when perturbing properties for purposes of computing a finite difference approximation to the ODE system Jacobian.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeAlgorithm</name>
      <defaultValue>var_str('rungeKuttaCashKarp')</defaultValue>
      <description>The algorithm to use in the ODE solver.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeAlgorithmNonJacobian</name>
      <defaultValue>var_str('rungeKuttaCashKarp')</defaultValue>
      <description>The algorithm to use in the ODE solver.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeLatentIntegratorType</name>
      <defaultValue>var_str('trapezoidal')</defaultValue>
      <description>The type of integrator to use for latent variables.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeLatentIntegratorOrder</name>
      <defaultValue>15</defaultValue>
      <description>The order of the integrator for latent variables.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>odeLatentIntegratorIntervalsMaximum</name>
      <defaultValue>1000</defaultValue>
      <description>The maximum number of intervals allowed in the integrator for latent variables.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>profileOdeEvolver</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not to profile the ODE evolver.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reuseODEStepSize</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, re-use the previous ODE step size when resuming the evolution of a node. Otherwise, the initial step size is not specified.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>enforceNonNegativity</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, properties that are marked as non-negative (e.g. masses) will be evolved in such a way as to enforce that non-negativity.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerTreeNodeMerger"     name="mergerTreeNodeMerger_"     source="parameters"/>
    <objectBuilder class="nodeOperator"             name="nodeOperator_"             source="parameters"/>
    <objectBuilder class="mergerTreeEvolveProfiler" name="mergerTreeEvolveProfiler_" source="parameters"/>
    !!]
    self=mergerTreeNodeEvolverStandard(                                                                                                         &
         &                                                                        odeToleranceAbsolute                                        , &
         &                                                                        odeToleranceRelative                                        , &
         &                             enumerationStandardODEAlgorithmEncode(char(odeAlgorithm                       ),includesPrefix=.false.), &
         &                             enumerationStandardODEAlgorithmEncode(char(odeAlgorithmNonJacobian            ),includesPrefix=.false.), &
         &                                                                        odeJacobianStepSizeRelative                                 , &
         &                             enumerationLatentIntegratorTypeEncode(char(odeLatentIntegratorType            ),includesPrefix=.false.), &
         &                                                                        odeLatentIntegratorOrder                                    , &
         &                                                                        odeLatentIntegratorIntervalsMaximum                         , &
         &                                                                        profileOdeEvolver                                           , &
         &                                                                        reuseODEStepSize                                            , &
         &                                                                        enforceNonNegativity                                        , &
         &                                                                        mergerTreeNodeMerger_                                       , &
         &                                                                        nodeOperator_                                               , &
         &                                                                        mergerTreeEvolveProfiler_                                     &
         &                            )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeNodeMerger_"    />
    <objectDestructor name="nodeOperator_"            />
    <objectDestructor name="mergerTreeEvolveProfiler_"/>
    !!]
    return
  end function standardConstructorParameters

   function standardConstructorInternal(odeToleranceAbsolute,odeToleranceRelative,odeAlgorithm,odeAlgorithmNonJacobian,odeJacobianStepSizeRelative,odeLatentIntegratorType,odeLatentIntegratorOrder,odeLatentIntegratorIntervalsMaximum,profileOdeEvolver,reuseODEStepSize,enforceNonNegativity,mergerTreeNodeMerger_,nodeOperator_,mergerTreeEvolveProfiler_) result(self)
     !!{
     Internal constructor for the \refClass{mergerTreeNodeEvolverStandard} merger tree node evolver class.
     !!}
     use :: Error                , only : Error_Report
     use :: Numerical_ODE_Solvers, only : GSL_ODEIV2_Step_RK2  , GSL_ODEIV2_Step_RK4    , GSL_ODEIV2_Step_RK8PD, GSL_ODEIV2_Step_RKCK       , &
          &                               GSL_ODEIV2_Step_RKF45, GSL_ODEIV2_Step_msAdams, GSL_ODEIV2_step_BSimp, GSL_ODEIV2_step_MSBDFActive
     implicit none
     type            (mergerTreeNodeEvolverStandard      )                        :: self
     type            (enumerationLatentIntegratorTypeType), intent(in   )         :: odeLatentIntegratorType
     type            (enumerationStandardODEAlgorithmType), intent(in   )         :: odeAlgorithm                       , odeAlgorithmNonJacobian
     integer                                              , intent(in   )         :: odeLatentIntegratorIntervalsMaximum, odeLatentIntegratorOrder
     double precision                                     , intent(in   )         :: odeToleranceAbsolute               , odeToleranceRelative    , &
          &                                                                          odeJacobianStepSizeRelative
     logical                                              , intent(in   )         :: profileOdeEvolver                  , reuseODEStepSize        , &
          &                                                                          enforceNonNegativity
     class           (mergerTreeNodeMergerClass          ), intent(in   ), target :: mergerTreeNodeMerger_
     class           (nodeOperatorClass                  ), intent(in   ), target :: nodeOperator_
     class           (mergerTreeEvolveProfilerClass      ), intent(in   ), target :: mergerTreeEvolveProfiler_
    !![
    <constructorAssign variables="odeToleranceAbsolute, odeToleranceRelative, odeJacobianStepSizeRelative, odeLatentIntegratorType, odeLatentIntegratorOrder, odeLatentIntegratorIntervalsMaximum, profileOdeEvolver, reuseODEStepSize, enforceNonNegativity, *mergerTreeNodeMerger_, *nodeOperator_, *mergerTreeEvolveProfiler_"/>
    !!]

     ! Construct ODE solver object.
     self%useJacobian=.false.
     select case (odeAlgorithm           %ID)
     case (standardODEAlgorithmRungeKuttaCashKarp     %ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_RKCK
     case (standardODEAlgorithmRungeKuttaSecondOrder  %ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_RK2
     case (standardODEAlgorithmRungeKutta             %ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_RK4
     case (standardODEAlgorithmRungeKuttaFehlberg     %ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_RKF45
     case (standardODEAlgorithmRungeKuttaPrinceDormand%ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_RK8PD
     case (standardODEAlgorithmMultistepAdams         %ID)
        self%odeAlgorithm=GSL_ODEIV2_Step_msAdams
     case (standardODEAlgorithmBulirschStoer          %ID)
        self%odeAlgorithm=GSL_ODEIV2_step_BSimp
        self%useJacobian             =.true.
     case (standardODEAlgorithmBDF                    %ID)
        self%odeAlgorithm=GSL_ODEIV2_step_MSBDFActive
        self%useJacobian             =.true.
     case default
        call Error_Report('odeAlgorithm is unrecognized'//{introspection:location})
     end select
     ! Construct non-Jacobian ODE solver object.
     select case (odeAlgorithmNonJacobian%ID)
     case (standardODEAlgorithmRungeKuttaCashKarp     %ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_RKCK
     case (standardODEAlgorithmRungeKuttaSecondOrder  %ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_RK2
     case (standardODEAlgorithmRungeKutta             %ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_RK4
     case (standardODEAlgorithmRungeKuttaFehlberg     %ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_RKF45
     case (standardODEAlgorithmRungeKuttaPrinceDormand%ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_RK8PD
     case (standardODEAlgorithmMultistepAdams         %ID)
        self%odeAlgorithmNonJacobian=GSL_ODEIV2_Step_msAdams
     case default
        call Error_Report('odeAlgorithm is unrecognized'//{introspection:location})
     end select
     ! Initialize
     self%propertyCountPrevious    =-1
     self%propertyCountMaximum     = 0
     self%timePrevious             =-1.0d0
     self%countEvaluationsToSuccess=0_c_size_t
     self%stepTimer                =timer()
     return
   end function standardConstructorInternal

  subroutine standardAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, subhaloPromotionEvent
    implicit none
    class(mergerTreeNodeEvolverStandard), intent(inout) :: self

    call subhaloPromotionEvent%attach(self,standardNodeSubhaloPromotion,openMPThreadBindingAtLevel,label='mergerTreeNodeEvolver')
    return
  end subroutine standardAutoHook

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeNodeEvolverStandard} merger tree node evolver class.
    !!}
    use :: Events_Hooks, only : subhaloPromotionEvent
    implicit none
    type(mergerTreeNodeEvolverStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeNodeMerger_"    />
    <objectDestructor name="self%nodeOperator_"            />
    <objectDestructor name="self%mergerTreeEvolveProfiler_"/>
    !!]
    if (subhaloPromotionEvent%isAttached(self,standardNodeSubhaloPromotion)) call subhaloPromotionEvent%detach(self,standardNodeSubhaloPromotion)
    return
  end subroutine standardDestructor

  subroutine standardEvolve(self,tree,node,timeEnd,interrupted,functionInterrupt,galacticStructureSolver__,treeLock,systemClockMaximum,status)
    !!{
    Evolves {\normalfont \ttfamily node} to time {\normalfont \ttfamily timeEnd}, or until evolution is interrupted.
    !!}
    use            :: Display               , only : displayIndent              , displayMessage                                  , displayUnindent                                , displayMagenta    , &
         &                                           displayReset
    use            :: Calculations_Resets   , only : Calculations_Reset
    use            :: Error                 , only : Error_Report               , Warn                                            , errorStatusFail                                , errorStatusSuccess, &
          &                                          errorStatusXCPU            , errorStatusUnderflow
    use            :: Galacticus_Nodes      , only : interruptTask              , mergerTree                                      , nodeComponentBasic                             , propertyTypeActive, &
          &                                          propertyTypeAll            , propertyTypeInactive                            , propertyTypeNone                               , rateComputeState  , &
          &                                          treeNode                   , propertyTypeNumerics
    use, intrinsic :: ISO_C_Binding         , only : c_funloc                   , c_funptr                                        , c_null_funptr
    use            :: Numerical_Integration2, only : integratorMultiVectorized1D, integratorMultiVectorizedCompositeGaussKronrod1D, integratorMultiVectorizedCompositeTrapezoidal1D
    use            :: ODE_Solver_Error_Codes, only : odeSolverInterrupt
    !![
    <include directive="preEvolveTask"      type="moduleUse">
    !!]
    include 'objects.tree_node.pre_evolve.modules.inc'
    !![
    </include>
    <include directive="scaleSetTask"       type="moduleUse">
    !!]
    include 'objects.tree_node.set_scale.modules.inc'
    !![
    </include>
    <include directive="inactiveSetTask"    type="moduleUse">
    !!]
    include 'objects.tree_node.set_inactive.modules.inc'
    !![
    </include>
    <include directive="analyticSolverTask" type="moduleUse">
    !!]
    include 'objects.tree_node.analytic_solver_task.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (mergerTreeNodeEvolverStandard      )             , intent(inout), target  :: self
    type            (mergerTree                         )             , intent(inout)          :: tree
    type            (treeNode                           )             , intent(inout), pointer :: node
    double precision                                                  , intent(in   )          :: timeEnd
    logical                                                           , intent(  out)          :: interrupted
    procedure       (interruptTask                      )             , intent(  out), pointer :: functionInterrupt
    class           (galacticStructureSolverClass       )             , intent(in   ), target  :: galacticStructureSolver__
    class           (ompLockClass                       )             , intent(inout)          :: treeLock
    integer         (kind_int8                          ), optional   , intent(in   )          :: systemClockMaximum
    integer                                              , optional   , intent(  out)          :: status
    procedure       (standardErrorHandler               )                            , pointer :: errorHandler
    class           (nodeComponentBasic                 )                            , pointer :: basic
    class           (integratorMultiVectorized1D        ), allocatable                         :: integrator_
    type            (odeSolver                          )                            , target  :: solver
    logical                                                                                    :: solvedAnalytically       , solvedNumerically, &
         &                                                                                        jacobianSolver           , solverInitialized
    double precision                                                                           :: timeStart                , stepSize
    integer                                                                                    :: lengthMaximum            , odeStatus        , &
         &                                                                                        odeAlgorithm
    integer         (c_size_t                           )                                      :: i
    integer         (kind_int8                          )                                      :: systemClockCount
    type            (varying_string                     )                                      :: message                  , line
    character       (len =12                            )                                      :: label

    ! Set status to success.
    if (present(status)) status=errorStatusSuccess
    ! Set time limit.
    if (present(systemClockMaximum)) then
       self%systemClockMaximum=systemClockMaximum
    else
       self%systemClockMaximum=-1_kind_int8
    end if
    ! Call functions to perform any pre-evolution tasks.
    call treeLock              %set                     (    )
    call self    %nodeOperator_%differentialEvolutionPre(node)
    !![
    <include directive="preEvolveTask" type="functionCall" functionType="void">
     <functionArgs>node</functionArgs>
    !!]
    include 'objects.tree_node.pre_evolve.inc'
    !![
    </include>
    !!]
    call treeLock              %unset                   (    )

    ! Determine the end time for this node - either the specified end time, or the time associated with the parent node, whichever
    ! occurs first.
    basic          => node %basic    ()
    timeStart      =  basic%time     ()
    timeStartSaved =        timeStart
    ! Set a pointer to the galactic structure solver. We do this to ensure that the solver we act on here is the same one as used
    ! by the calling function (which will therefore be called by various event hook triggers).
    self%galacticStructureSolver_ => galacticStructureSolver__
    ! Ensure calculations are reset for this new step.
    call Calculations_Reset(node)
    ! Attempt to find analytic solutions.
    solvedAnalytically=.false.
    !![
    <include directive="analyticSolverTask" type="functionCall" functionType="void">
     <functionArgs>node,timeStart,timeEnd,solvedAnalytically</functionArgs>
    !!]
    include 'objects.tree_node.analytic_solver_task.inc'
    !![
    </include>
    !!]
    ! Check if an analytic solution was available - use numerical solution if not.
    if (solvedAnalytically) then
       ! An analytic solution was available. Record that no interrupt therefore occurred.
       interrupted=.false.
    else
       ! Compute offsets into serialization arrays for rates and scales.
       call node%serializationOffsets(propertyTypeAll)
       ! Find number of all evolvable variables (active and inactive) for this node.
       self%propertyCountAll=node%serializeCount(propertyTypeAll)
       ! Allocate pointer arrays if necessary.
       if (self%propertyCountAll > self%propertyCountMaximum) then
          if (allocated(self%propertyValuesActive)) then
             deallocate(self%isNonNegative                )
             deallocate(self%propertyValuesActive         )
             deallocate(self%propertyValuesActiveSaved    )
             deallocate(self%propertyValuesInactive       )
             deallocate(self%propertyValuesInactiveSaved  )
             deallocate(self%propertyScalesActive         )
             deallocate(self%propertyScalesInactive       )
             deallocate(self%propertyValuesPrevious       )
             deallocate(self%propertyRatesPrevious        )
             deallocate(self%propertyErrors               )
             deallocate(self%propertyTolerances           )
             deallocate(self%odeTolerancesInactiveRelative)
             deallocate(self%odeTolerancesInactiveAbsolute)
          end if
          allocate(self%isNonNegative                (self%propertyCountAll))
          allocate(self%propertyValuesActive         (self%propertyCountAll))
          allocate(self%propertyValuesActiveSaved    (self%propertyCountAll))
          allocate(self%propertyValuesInactive       (self%propertyCountAll))
          allocate(self%propertyValuesInactiveSaved  (self%propertyCountAll))
          allocate(self%propertyScalesActive         (self%propertyCountAll))
          allocate(self%propertyScalesInactive       (self%propertyCountAll))
          allocate(self%propertyValuesPrevious       (self%propertyCountAll))
          allocate(self%propertyRatesPrevious        (self%propertyCountAll))
          allocate(self%propertyErrors               (self%propertyCountAll))
          allocate(self%propertyTolerances           (self%propertyCountAll))
          allocate(self%odeTolerancesInactiveRelative(self%propertyCountAll))
          allocate(self%odeTolerancesInactiveAbsolute(self%propertyCountAll))
          self%propertyCountMaximum  =self%propertyCountAll
          self%propertyValuesPrevious=0.0d0
          self%isNonNegative         =.false.
       end if
       ! Iterate until the step has been solved numerically.
       solvedNumerically=.false.
       jacobianSolver   =self%useJacobian
       do while (.not.solvedNumerically)
          if (jacobianSolver) then
             self%propertyTypeODE       =propertyTypeActive
             self%propertyTypeIntegrator=propertyTypeInactive
          else
             self%propertyTypeODE       =propertyTypeNumerics
             self%propertyTypeIntegrator=propertyTypeNone
          end if
          rateComputeState=self%propertyTypeODE
          ! Determine active and inactive properties.
          call node%odeStepInactivesInitialize()
          call node%odeStepAnalyticsInitialize()
          if (jacobianSolver) then
             !![
             <include directive="inactiveSetTask" type="functionCall" functionType="void">
              <functionArgs>node</functionArgs>
             !!]
             include 'objects.tree_node.set_inactive.inc'
             !![
             </include>
             !!]
             call self%nodeOperator_%differentialEvolutionInactives(node)
          end if
          ! Determine analytically-solvable properties.
          call self%nodeOperator_%differentialEvolutionAnalytics(node)
          ! Compute offsets into serialization arrays for rates and scales.
          call node%serializationOffsets(self%propertyTypeODE       )
          call node%serializationOffsets(self%propertyTypeIntegrator)
          ! Find number of active and inactive variables to evolve for this node.
          self%propertyCountActive  =node%serializeCount(self%propertyTypeODE       )
          self%propertyCountInactive=node%serializeCount(self%propertyTypeIntegrator)
          ! Serialize property values to array, and save a copy in case of a need to reset.
          call node%serializeValues(self%propertyValuesActive  ,self%propertyTypeODE       )
          call node%serializeValues(self%propertyValuesInactive,self%propertyTypeIntegrator)
          self%propertyValuesActiveSaved  (1:self%propertyCountActive  )=self%propertyValuesActive  (1:self%propertyCountActive  )
          self%propertyValuesInactiveSaved(1:self%propertyCountInactive)=self%propertyValuesInactive(1:self%propertyCountInactive)
          ! Find properties that must be non-negative.
          if (self%enforceNonNegativity) call node%serializeNonNegative(self%isNonNegative)
          ! Compute scales for all properties and extract from the node.
          call node%odeStepScalesInitialize()
          !![
          <include directive="scaleSetTask" type="functionCall" functionType="void">
           <functionArgs>node</functionArgs>
          !!]
          include 'objects.tree_node.set_scale.inc'
          !![
          </include>
          !!]
          call self%nodeOperator_%differentialEvolutionScales(node)
          call node%serializeScales(self%propertyScalesActive  ,self%propertyTypeODE       )
          call node%serializeScales(self%propertyScalesInactive,self%propertyTypeIntegrator)
          ! Check for zero property scales which will cause floating point overflow in the ODE solver.
          if (any(self%propertyScalesActive(1:self%propertyCountActive) == 0.0d0)) then
             message=displayMagenta()//'WARNING:'//displayReset()//' Zero entry in ODE system scales for node'
             call Warn         (message)
             call displayIndent(message)
             lengthMaximum=0
             do i=1,self%propertyCountActive
                lengthMaximum=max(lengthMaximum,len(node%nameFromIndex(int(i),self%propertyTypeODE)))
             end do
             line=repeat("―",lengthMaximum)//repeat("―――――――――――――――",2)
             call displayMessage(line)
             call displayMessage(repeat(" ",lengthMaximum)//' : y            : yScale')
             call displayMessage(line)
             do i=1,self%propertyCountActive
                message=node%nameFromIndex(int(i),self%propertyTypeODE)
                message=repeat(" ",lengthMaximum-len(message))//message
                write (label,'(e12.6)') self%propertyValuesActive(i)
                message=message//" : "//label
                write (label,'(e12.6)') self%propertyScalesActive(i)
                message=message//" : "//label
                call displayMessage(message)
             end do
             call displayMessage(line)
             call node%serializeASCII()
             call displayUnindent('done')
          end if          
          ! Assign module global pointer to this node.
          self_                   => self
          solver_                 => solver
          self   %activeTreeIndex =  tree%index
          self   %activeNode      => node
          ! Reset interrupt variables.
          self%interruptFirstFound    =  .false.
          self%timeInterruptFirst     =  0.0d0
          self%functionInterruptFirst => null()
          ! Call ODE solver routines.
          solvedNumerically=.true.
          stepSize         =-1.0d0
          if (timeStart /= timeEnd .and. self%propertyCountActive > 0) then
             self%propertyCountPrevious=self%propertyCountActive
             self%trialCount           =0
             solverInitialized         =.false.
             odeStatus                 =errorStatusFail
             do while (self%trialCount < trialCountMaximum .and. .not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt))
                ! Initialize state. Stepsize is reduced by half on each successive trial.
                if (.not.solverInitialized) then
                   if (jacobianSolver) then
                      ! Initialize integrator.
                      if (.not.allocated(integrator_)) then
                         select case (self%odeLatentIntegratorType%ID)
                         case (latentIntegratorTypeGaussKronrod%ID)
                            allocate(integratorMultiVectorizedCompositeGaussKronrod1D :: integrator_)
                            select type (integrator_)
                            type is (integratorMultiVectorizedCompositeGaussKronrod1D)
                               call integrator_%initialize(self%odeLatentIntegratorIntervalsMaximum,self%odeLatentIntegratorOrder)
                            end select
                         case (latentIntegratorTypeTrapezoidal %ID)
                            allocate(integratorMultiVectorizedCompositeTrapezoidal1D  :: integrator_)
                            select type (integrator_)
                            type is (integratorMultiVectorizedCompositeTrapezoidal1D )
                               call integrator_%initialize(self%odeLatentIntegratorIntervalsMaximum                              )
                            end select
                         end select
                      end if
                      odeAlgorithm=self%odeAlgorithm
                   else
                      odeAlgorithm=self%odeAlgorithmNonJacobian
                   end if
                   if (self%profileOdeEvolver) call self%stepTimer%start()
                   !![
                   <conditionalCall>
		    <call>
                     solver=odeSolver(                                                     &amp;
                      &amp;                             self%propertyCountActive         , &amp;
                      &amp;                                  standardODEs                , &amp;
                      &amp;           toleranceAbsolute=self%odeToleranceAbsolute        , &amp;
                      &amp;           toleranceRelative=self%odeToleranceRelative        , &amp;
                      &amp;           scale            =self%propertyScalesActive        , &amp;
                      &amp;           stepperType      =     odeAlgorithm                , &amp;
                      &amp;           postStep         =     standardPostStepProcessing  , &amp;
                      &amp;           finalState       =     standardFinalStateProcessing, &amp;
                      &amp;           errorHandler     =     standardErrorHandler          &amp;                                      
                      &amp;           {conditions}                                         &amp;
                      &amp;          )
		    </call>
                    <argument name="jacobian"                value="standardODEsJacobian"      condition="jacobianSolver"           />
                    <argument name="integrator"              value="integrator_"               condition="jacobianSolver"           />
                    <argument name="integratorErrorTolerant" value=".true."                    condition="jacobianSolver"           />
                    <argument name="integrands"              value="standardIntegrands"        condition="jacobianSolver"           />
                    <argument name="errorAnalyzer"           value="standardStepErrorAnalyzer" condition="self%profileOdeEvolver"   />
                    <argument name="isNonNegative"           value="self%isNonNegative"        condition="self%enforceNonNegativity"/>
                   </conditionalCall>
                   !!]
                   solverInitialized=.true.
                end if
                if (self%reuseODEStepSize) then
                   stepSize      =node%timeStep()/2.0d0**self%trialCount
                else
                   stepSize      =-1.0d0
                end if
                self%timePrevious=-1.0d0
                timeStart        =timeStartSaved
                if (jacobianSolver) then
                   self%odeTolerancesInactiveRelative(1:self%propertyCountInactive)=self%odeToleranceRelative
                   self%odeTolerancesInactiveAbsolute(1:self%propertyCountInactive)=self%odeToleranceAbsolute*self%propertyScalesInactive(1:self%propertyCountInactive)
                   call integrator_%tolerancesSet(self%odeTolerancesInactiveAbsolute(1:self%propertyCountInactive),self%odeTolerancesInactiveRelative(1:self%propertyCountInactive))
                   call solver     %solve(timeStart,timeEnd,self%propertyValuesActive,status=odeStatus,xStep=stepSize,z=self%propertyValuesInactive(1:self%propertyCountInactive))
                else
                   call solver     %solve(timeStart,timeEnd,self%propertyValuesActive,status=odeStatus,xStep=stepSize                                                            )
                end if
                ! If non-negativity is being enforced for some properties, check for failures in the step and, if these are due to such a property
                ! being negative, zero that property and ignore the failure. This is necessary as it may be that a non-negative property has
                ! non-zero derivative as it approaches zero, such that even an arbitrarily small step size still results in it becoming negative.                
                if     (                                                                    &
                     &   self%enforceNonNegativity                                          &
                     &  .and.                                                               &
                     &   .not.                                                              &
                     &    (                                                                 &
                     &      odeStatus == errorStatusSuccess                                 &
                     &     .or.                                                             &
                     &      odeStatus == odeSolverInterrupt                                 &
                     &    )                                                                 &
                     &  .and.                                                               &
                     &   any(                                                               &
                     &        self%propertyValuesActive(1:self%propertyCountActive) < 0.0d0 &
                     &       .and.                                                          &
                     &        self%isNonNegative       (1:self%propertyCountActive)         &
                     &      )                                                               &
                     & ) then
                   where  (                                                               &
                        &   self%propertyValuesActive(1:self%propertyCountActive) < 0.0d0 &
                        &  .and.                                                          &
                        &   self%isNonNegative       (1:self%propertyCountActive)         &
                        & )
                      self%propertyValuesActive(1:self%propertyCountActive)=0.0d0
                   end where
                   ! Ignore the returned step size in this case.
                   stepSize =node%timeStep()/2.0d0**self%trialCount
                   ! Reset status to success.
                   odeStatus=errorStatusSuccess
                end if
                ! Check for failure.
                if (.not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt)) then
                   ! Increment number of trials.
                   self%trialCount=self%trialCount+1
                   ! Restore state of the node.
                   self%propertyValuesActive  (1:self%propertyCountActive  )=self%propertyValuesActiveSaved  (1:self%propertyCountActive  )
                   self%propertyValuesInactive(1:self%propertyCountInactive)=self%propertyValuesInactiveSaved(1:self%propertyCountInactive)
                end if
                ! Check for exceeding wall time.
                if (self%systemClockMaximum > 0_kind_int8) then
                   call System_Clock(systemClockCount)
                   if (systemClockCount > self%systemClockMaximum) then
                      if (present(status)) then
                         call displayMessage('maximum wall time exceeded'                          )
                         status=errorStatusXCPU
                         return
                      else
                         call Error_Report  ('maximum wall time exceeded'//{introspection:location})
                      end if
                   end if
                end if
             end do
             if (.not.(odeStatus == errorStatusSuccess .or. odeStatus == odeSolverInterrupt)) then
                if (jacobianSolver) then
                   ! Restore state of the node, switch to using a non-Jacobian solver, and indicate that the step was not solved numerically.
                   solvedNumerically=.false.
                   jacobianSolver   =.false.
                   solverInitialized=.false.
                   self%propertyValuesActive  (1:self%propertyCountActive  )=self%propertyValuesActiveSaved  (1:self%propertyCountActive  )
                   self%propertyValuesInactive(1:self%propertyCountInactive)=self%propertyValuesInactiveSaved(1:self%propertyCountInactive)
                else if (present(status)) then
                   call displayMessage('ODE integration failed '//{introspection:location})
                   status=errorStatusFail
                   return
                else
                   call Error_Report  ('ODE integration failed '//{introspection:location})
                end if
             end if
          end if
          if (solvedNumerically) then
             ! Extract values.
             call node%deserializeValues(self%propertyValuesActive  ,self%propertyTypeODE       )
             call node%deserializeValues(self%propertyValuesInactive,self%propertyTypeIntegrator)
             call self%nodeOperator_%differentialEvolutionSolveAnalytics(node,timeEnd)
             ! Flag interruption if one occurred, and ensure that the time is matched precisely to the end or interrupt time (can differ
             ! due to finite precision of the ODE integrator).
             if (self%timeInterruptFirst /= 0.0d0) then
                interrupted=.true.
                functionInterrupt => self%functionInterruptFirst
                call basic%timeSet    (self%timeInterruptFirst)
                call node %timeStepSet(     -1.0d0            )
             else
                interrupted=.false.
                call basic%timeSet    (     timeEnd           )
                call node %timeStepSet(     stepSize          )
             end if
          end if
       end do
    end if
    ! Call routines to perform any post-evolution tasks.
    if (associated(node)) then
       call treeLock              %set                      (    )
       call self    %nodeOperator_%differentialEvolutionPost(node)
       !![
       <eventHook name="postEvolve">
        <callWith>node</callWith>
       </eventHook>
       !!]
       call treeLock              %unset                    (    )
    end if
    return
  end subroutine standardEvolve

  logical function standardIsAccurate(self,valueNode,valueExpected)
    !!{
    Return true if a tree node property is within expected accuracy of a given value.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (mergerTreeNodeEvolverStandard), intent(inout) :: self
    double precision                               , intent(in   ) :: valueNode, valueExpected

    standardIsAccurate=Values_Agree(valueNode,valueExpected,self%odeToleranceAbsolute,self%odeToleranceRelative)
    return
  end function standardIsAccurate

  subroutine standardIntegrands(time,propertyValues,propertyRates,inactivePropertyInitialValues,evaluate,integrands)
    !!{
    Evaluates integrands for node evolution.
    !!}
    use :: Galacticus_Nodes, only : interruptTask, propertyTypeInactive, rateComputeState
    implicit none
    double precision               , intent(in   ), dimension(:  ) :: time
    double precision               , intent(in   ), dimension(:,:) :: propertyValues               , propertyRates
    double precision               , intent(in   ), dimension(:  ) :: inactivePropertyInitialValues
    logical                        , intent(inout), dimension(:  ) :: evaluate
    double precision               , intent(  out), dimension(:,:) :: integrands
    procedure       (interruptTask), pointer                       :: functionInterrupt
    logical                                                        :: interrupt
    integer                                                        :: iTime
    ! "evaluate" array is currently not used. It indicates which integrands must be evaluated, and which can (optionally) be
    ! ignored as they have already converged to the required tolerance. It is currently not used because the potential for
    ! significant speed up appears to be small based on profiling. This will be model-dependent though, so this decision can be
    ! revisited.
    !$GLC attributes unused :: evaluate

    ! Set state to indicate that rate calculations are for inactive variables.
    rateComputeState=self_%propertyTypeIntegrator
    ! Iterate over times.
    do iTime=1,size(time)
       ! Deserialize the active values.
       call self_%activeNode%deserializeValues(propertyValues               (1:self_%propertyCountActive  ,iTime),self_%propertyTypeODE     )
       call self_%activeNode%deserializeRates (propertyRates                (1:self_%propertyCountActive  ,iTime),self_%propertyTypeODE     )
       call self_%activeNode%deserializeValues(inactivePropertyInitialValues(1:self_%propertyCountInactive      ),      propertyTypeInactive)
       call self_%nodeOperator_%differentialEvolutionSolveAnalytics(self_%activeNode,time(iTime))
       ! If past the time of the first interrupt, integrands are set to zero. Otherwise, evaluate integrands.
       if (self_%interruptFirstFound .and. time(iTime) >= self_%timeInterruptFirst) then
          integrands(:,iTime)=0.0d0
       else
          ! Set derivatives to zero initially.
          call self_%activeNode%odeStepRatesInitialize()
          ! Compute derivatives.
          call self_%galacticStructureSolver_%revert(self_%activeNode                                                                )
          call standardDerivativesCompute           (self_%activeNode,interrupt,functionInterrupt,self_%propertyTypeIntegrator)
          ! Serialize rates into integrand array.
          call self_%activeNode%serializeRates(integrands(:,iTime),self_%propertyTypeIntegrator)
       end if
    end do
    ! Set state to indicate that rate calculations are for active variables.
    rateComputeState=self_%propertyTypeODE
    return
  end subroutine standardIntegrands

  subroutine standardFinalStateProcessing(time,propertyValues)
    !!{
    Perform any actions based on the final state of the ODE step.
    !!}
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: propertyValues

    call self_%activeNode   %deserializeValues                  (propertyValues(1:self_%propertyCountActive),self_%propertyTypeODE)
    call self_%nodeOperator_%differentialEvolutionSolveAnalytics(                 self_%activeNode          ,time                 )
    call self_%nodeOperator_%differentialEvolutionStepFinalState(                 self_%activeNode                                )
    return
  end subroutine standardFinalStateProcessing

  integer function standardODEs(time,y,dydt)
    !!{
    Function which evaluates the set of ODEs for the evolution of a specific node.
    !!}
#ifdef DEBUGGING
    use :: Debugging             , only : isDebugging    , debugLog
    use :: ISO_Varying_String    , only : varying_string , assignment(=), var_str
    use :: Galacticus_Nodes      , only : propertyTypeAll
    use :: String_Handling       , only : operator(//)
#endif
    use :: Error                 , only : errorStatusXCPU
    use :: Galacticus_Nodes      , only : interruptTask  , nodeComponentBasic
    use :: Interface_GSL         , only : GSL_Success
    use :: ODE_Solver_Error_Codes, only : interruptedAtX , odeSolverInterrupt
    implicit none
    double precision                    , intent(in   )                                       :: time
    double precision                    , intent(in   ), dimension(:                        ) :: y
    double precision                    , intent(  out), dimension(:                        ) :: dydt
    logical                                                                                   :: interrupt
    procedure       (interruptTask     ), pointer                                             :: functionInterrupt
    class           (nodeComponentBasic), pointer                                             :: basic
    integer         (kind_int8         )                                                      :: systemClockCount
#ifdef DEBUGGING
    integer         (c_size_t          )                                                      :: iProperty
    character       (len=32            )                                                      :: label
    type            (varying_string    )                                                      :: message
    double precision                                   , dimension(self_%propertyCountActive) :: valuesDebug, ratesDebug
#endif

    ! Check for exceeding wall time.
    if (self_%systemClockMaximum > 0_kind_int8) then
       call System_Clock(systemClockCount)
       if (systemClockCount > self_%systemClockMaximum) then
          standardODEs=errorStatusXCPU
          return
       end if
    end if
    ! Return success by default.
    standardODEs=GSL_Success
    ! Check if we can reuse the previous derivatives.
    if   (                                                                                                                &
       &   time                                             == self_%timePrevious                                         &
       &  .and.                                                                                                           &
       &   all(y(             1:self_%propertyCountActive)  == self_%propertyValuesPrevious(1:self_%propertyCountActive)) &
       & ) then
       dydt                  (1:self_%propertyCountActive)  =  self_%propertyRatesPrevious (1:self_%propertyCountActive)
       return
    else
       self_%timePrevious                                       =time
       self_%propertyValuesPrevious(1:self_%propertyCountActive)=y                         (1:self_%propertyCountActive)
    end if
    ! Extract values.
    call self_%activeNode%deserializeValues(y(1:self_%propertyCountActive),self_%propertyTypeODE)
    call self_%nodeOperator_%differentialEvolutionSolveAnalytics(self_%activeNode,time)
    ! If the node is significantly inaccurate (as judged by the node time being different from the system time), then set rates to
    ! zero, as the ODE solver is likely just taking a step which is too large.
    basic => self_%activeNode%basic()
    if (.not.self_%isAccurate(basic%time(),time)) then
       dydt=0.0d0
       return
    end if
#ifdef DEBUGGING
    if (isDebugging()) then
       write (label,'(f9.4)') time
       message=var_str("step: ")//trim(adjustl(label))//' [treeIndex: '//self_%activeNode%hostTree%index//' ; nodeIndex '//self_%activeNode%index()//']'
       call debugLog(message)
    end if
#endif
    ! Set derivatives to zero initially.
    call self_%activeNode%odeStepRatesInitialize()
    if (self_%interruptFirstFound .and. time >= self_%timeInterruptFirst) then
       ! Already beyond the location of the first interrupt, simply set any analytic solutions to the interrupt time and return
       ! zero derivatives.
       dydt                              (1:self_%propertyCountActive)=0.0d0
       self_%propertyRatesPrevious(1:self_%propertyCountActive)=0.0d0
       call self_%nodeOperator_%differentialEvolutionSolveAnalytics(self_%activeNode,self_%timeInterruptFirst)
    else
       ! Compute derivatives.
       call standardDerivativesCompute(self_%activeNode,interrupt,functionInterrupt,self_%propertyTypeODE)
       ! Check whether an interrupt has been requested.
       select case (interrupt)
       case (.false.)
          ! No interrupt - place derivatives into ODE arrays.
          call self_%activeNode%serializeRates(dydt(1:self_%propertyCountActive),self_%propertyTypeODE)
          self_%propertyRatesPrevious(1:self_%propertyCountActive)=dydt(1:self_%propertyCountActive)             
       case (.true.)
          ! Interrupt requested - freeze evolution and store the interrupt if it is the earliest one to occur.
          dydt                              (1:self_%propertyCountActive)=0.0d0
          self_%propertyRatesPrevious(1:self_%propertyCountActive)=0.0d0
          if (time < self_%timeInterruptFirst .or. .not.self_%interruptFirstFound) then
             self_%interruptFirstFound    =  .true.
             self_%timeInterruptFirst     =  time
             self_%functionInterruptFirst => functionInterrupt
             ! Let the ODE solver know that an interrupt occurred, and when it happened.
             standardODEs                 =  odeSolverInterrupt
             interruptedAtX               =  time
             return
          end if
       end select
    end if
#ifdef DEBUGGING
    if (isDebugging()) then
       call self_%activeNode%serializeValues(valuesDebug,self_%propertyTypeODE)
       call self_%activeNode%serializeRates (ratesDebug ,self_%propertyTypeODE)
       do iProperty=1,self_%propertyCountActive
          write (label,'(e12.6)') valuesDebug(iProperty)
          message="  value: "        //self_%activeNode%nameFromIndex(int(iProperty),self_%propertyTypeODE)//" "//trim(adjustl(label))
          call debugLog(message)          
          write (label,'(e12.6)') ratesDebug(iProperty)
          message="   rate: (total) "//self_%activeNode%nameFromIndex(int(iProperty),self_%propertyTypeODE)//" "//trim(adjustl(label))
          call debugLog(message)          
       end do
    end if
#endif
    return
  end function standardODEs

  integer function standardODEsJacobian(time,propertyValues0,derivativeRatesValues,derivativeRatesTime)
    !!{
    Function which evaluates the set of ODEs for the evolution of a specific node.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision                                                                               , intent(in   ) :: time
    double precision               , dimension(:                                                  ), intent(in   ) :: propertyValues0
    double precision               , dimension(:                                                  ), intent(  out) :: derivativeRatesValues        , derivativeRatesTime
    double precision               , dimension(self_%propertyCountActive                          )                :: propertyRates0               , propertyRates1     , &
         &                                                                                                            propertyValues1
    double precision               , dimension(self_%propertyCountActive,self_%propertyCountActive)                :: jacobian
    procedure       (interruptTask), pointer                                                                       :: functionInterrupt
    double precision               , parameter                                                                     :: deltaTiny            =1.0d-10
    logical                                                                                                        :: interrupt
    integer   (c_size_t)                                                                                           :: i
    double precision                                                                                               :: propertyValueDelta

    ! Return success by default.
    standardODEsJacobian=GSL_Success
    ! No explicit time dependence.
    derivativeRatesTime=0.0d0
    ! Check for interrupts.
    if (self_%interruptFirstFound .and. time >= self_%timeInterruptFirst) then
       ! Already beyond the location of the first interrupt, simply return zero derivatives.
       jacobian(1:self_%propertyCountActive,1:self_%propertyCountActive)=0.0d0
    else
       ! Compute rates at current parameter values.
       call self_%activeNode%deserializeValues     (propertyValues0(1:self_%propertyCountActive),self_%propertyTypeODE)
       call self_%nodeOperator_%differentialEvolutionSolveAnalytics(self_%activeNode,time)
       call self_%activeNode%odeStepRatesInitialize(                                                                  )
       call standardDerivativesCompute             (self_%activeNode,interrupt,functionInterrupt,self_%propertyTypeODE)
       call self_%activeNode%serializeRates        (propertyRates0                              ,self_%propertyTypeODE)
       ! If an interrupt was triggered, then derivatives will all be zero, so we set the Jacobian to zero here and exit.
       if (interrupt) then
          jacobian(1:self_%propertyCountActive,1:self_%propertyCountActive)=0.0d0
       else
          ! Iterate over parameters, computing Jacobian using finite differences.
          do i=1,self_%propertyCountActive
             ! To compute the finite difference we make a small perturbation in one parameter. If the parameter is non-zero, use a
             ! small, fractional perturbation. For parameters with zero value, use a perturbation equal to the absolute tolerance
             ! supplied to the ODE solver.
             if (propertyValues0(i)==0.0d0) then
                propertyValueDelta       =+self_%propertyScalesActive       (i)
             else
                propertyValueDelta       =+self_%odeJacobianStepSizeRelative    &
                     &                    *propertyValues0                         (i)
                if (abs(propertyValueDelta) < deltaTiny*self_%propertyScalesActive(i)) &
                     & propertyValueDelta=+self_%propertyScalesActive       (i)
             end if
             propertyValues1   =+propertyValues0
             propertyValues1(i)=+propertyValues1           (i) &
                  &             +propertyValueDelta
             call self_%activeNode%deserializeValues     (propertyValues1                             ,self_%propertyTypeODE)
             call self_%activeNode%odeStepRatesInitialize(                                                                  )
             call self_%galacticStructureSolver_%revert  (self_%activeNode                                                  )
             call standardDerivativesCompute             (self_%activeNode,interrupt,functionInterrupt,self_%propertyTypeODE)
             call self_%activeNode%serializeRates        (propertyRates1                              ,self_%propertyTypeODE)
             jacobian(i,:)=+(                  &
                  &          +propertyRates1   &
                  &          -propertyRates0   &
                  &         )                  &
                  &        /propertyValueDelta
          end do
       end if
    end if
    ! Map Jacobian back to output array.
    derivativeRatesValues=reshape(jacobian,[self_%propertyCountActive**2])
    return
  end function standardODEsJacobian

  subroutine standardDerivativesCompute(node,interrupt,functionInterruptReturn,propertyType)
    !!{
    Call routines to set all derivatives for {\normalfont \ttfamily node}.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    !![
    <include directive="rateComputeTask" type="moduleUse">
    !!]
    include 'objects.node.component.derivatives.modules.inc'
    !![
    </include>
    !!]
    implicit none
    type     (treeNode), intent(inout)          :: node
    logical            , intent(  out)          :: interrupt
    procedure(        ), intent(  out), pointer :: functionInterruptReturn
    integer            , intent(in   )          :: propertyType
    procedure(        )               , pointer :: functionInterrupt

    ! Initialize interrupt status.
    interrupt         =  .false.
    functionInterrupt => null()
    ! Call component routines to indicate that derivative calculation is commencing.
    call Calculations_Reset(node)
    ! Trigger an event to perform any pre-derivative calculations.
    !![
    <eventHook name="preDerivative">
     <callWith>node,propertyType</callWith>
    </eventHook>
    !!]
    ! Do not attempt to compute derivatives for nodes which are not solvable.
    if (.not.node%isSolvable) return
    ! Call component routines to compute derivatives.
    !![
    <include directive="rateComputeTask" type="functionCall" functionType="void">
     <functionArgs>node,interrupt,functionInterrupt,propertyType</functionArgs>
    !!]
    include 'objects.node.component.derivatives.inc'
    !![
    </include>
    !!]
    call self_%nodeOperator_%differentialEvolution(node,interrupt,functionInterrupt,propertyType)
    ! Return the procedure pointer.
    functionInterruptReturn => functionInterrupt
    return
  end subroutine standardDerivativesCompute

  subroutine standardErrorHandler(status,time,timeStep,y)
    !!{
    Handles errors in the ODE solver when evolving \glc\ nodes. Dumps the content of the node.
    !!}
    use            :: Display        , only : displayIndent      , displayMessage        , displayUnindent              , displayVerbosity, &
          &                                   displayVerbositySet, verbosityLevelStandard, enumerationVerbosityLevelType
    use, intrinsic :: ISO_C_Binding  , only : c_double           , c_int
    use            :: String_Handling, only : operator(//)
    implicit none
    integer         (kind=c_int                   ), intent(in   )                                       :: status
    real            (kind=c_double                ), intent(in   )                                       :: time          , timeStep
    real            (kind=c_double                ), intent(in   ), dimension(:                        ) :: y
    real            (kind=c_double                )               , dimension(self_%propertyCountActive) :: dydt          , yError  , &
         &                                                                                                  yTolerance
    type            (varying_string               )                                                      :: message       , line
    integer                                                                                              :: lengthMaximum
    type            (enumerationVerbosityLevelType)                                                      :: verbosityLevel
    integer         (c_size_t                     )                                                      :: i
    character       (len =12                      )                                                      :: label
    integer         (kind=c_int                   )                                                      :: odeStatus
    double precision                                                                                     :: stepFactor

    ! Check if this is the final trial for this node.
    if (self_%trialCount == trialCountMaximum-1) then
       ! Get the current errors and tolerances in the ODE driver.
       call solver_%errors(yError)
       ! Report the failure message.
       verbosityLevel=displayVerbosity()
       if (verbosityLevel < verbosityLevelStandard) call displayVerbositySet(verbosityLevelStandard)
       message="ODE solver failed with error code "
       message=message//status//" in tree #"//self_%activeTreeIndex
       call displayMessage(message)
       ! Dump all node properties.
       call self_%activeNode%serializeASCII()
       ! Report timestep.
       write (label,'(e12.6)') time
       message="time, timeStep = "//label
       write (label,'(e12.6)') timeStep
       message=message//", "//label
       call displayMessage(message)       
       ! Evaluate derivatives.
       odeStatus =standardODEs             (time,y,dydt)
       yTolerance=standardODEStepTolerances(     y     )
       call displayIndent('ODE system parameters')
       lengthMaximum=0
       do i=1,self_%propertyCountActive
          lengthMaximum=max(lengthMaximum,len(self_%activeNode%nameFromIndex(int(i),self_%propertyTypeODE)))
       end do
       line=repeat("―",lengthMaximum)//repeat("―――――――――――――――",5)
       call displayMessage(line)
       call displayMessage(repeat(" ",lengthMaximum)//' : y            : dy/dt        : yScale       : yError       : yErrorScaled')
       call displayMessage(line)
       do i=1,self_%propertyCountActive
          message=self_%activeNode%nameFromIndex(int(i),self_%propertyTypeODE)
          message=repeat(" ",lengthMaximum-len(message))//message
          write (label,'(e12.6)') y                         (i)
          message=message//" : "//label
          write (label,'(e12.6)') dydt                      (i)
          message=message//" : "//label
          write (label,'(e12.6)') self_%propertyScalesActive(i)
          message=message//" : "//label
          write (label,'(e12.6)') yError                    (i)
          message=message//" : "//label
          if (exponent(yError(i))-exponent(yTolerance(i)) < maxExponent(0.0d0)) then
             stepFactor=abs(yError(i))/yTolerance(i)
             write (label,'(e12.6)') stepFactor
          else
             label="infinity"
          end if
          message=message//" : "//label
          call displayMessage(message)
       end do
       call displayMessage(line)
       call displayUnindent('done')
       call displayVerbositySet(verbosityLevel)
    end if
    return
  end subroutine standardErrorHandler

  function standardODEStepTolerances(propertyValues)
    !!{
    Compute the tolerances on each property being evolved in the ODE system at the current timestep.
    !!}
    implicit none
    double precision                         , dimension(self_%propertyCountActive) :: standardODEStepTolerances
    double precision          , intent(in   ), dimension(self_%propertyCountActive) :: propertyValues
    integer         (c_size_t)                                                      :: i

    forall(i=1:self_%propertyCountActive)
       standardODEStepTolerances(i)=+    self_%odeToleranceRelative     &
            &                       *abs(      propertyValues      (i)) &
            &                       +    self_%odeToleranceAbsolute     &
            &                       *    self_%propertyScalesActive(i)
    end forall
    return
  end function standardODEStepTolerances

  subroutine standardPostStepProcessing(time,y,postStepStatus) bind(c)
    !!{
    Perform any post-step actions on the node.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double   , c_int
    use            :: Interface_GSL, only : GSL_Success, GSL_Continue
    !![
    <include directive="postStepTask" type="moduleUse">
    !!]
    include 'objects.tree_node.post_step.modules.inc'
    !![
    </include>
    !!]
    implicit none
    real   (kind=c_double), intent(in   ), value        :: time
    real   (kind=c_double), intent(inout), dimension(*) :: y
    integer(kind=c_int   ), intent(inout)               :: postStepStatus

    call self_%activeNode   %deserializeValues                  (y(1:self_%propertyCountActive),self_%propertyTypeODE)
    call self_%nodeOperator_%differentialEvolutionSolveAnalytics(    self_%activeNode          ,time                 )
    call self_%nodeOperator_%differentialEvolutionPostStep      (    self_%activeNode          ,postStepStatus       )
    !![
    <include directive="postStepTask" type="functionCall" functionType="void">
     <functionArgs>self_%activeNode,postStepStatus</functionArgs>
    !!]
    include 'objects.tree_node.post_step.inc'
    !![
    </include>
    !!]
    ! If the post-step processing returned a non-success error code - indicating that the node state was changed - reserialize the
    ! node state to the ODE solver arrays.
    if (postStepStatus /= GSL_Success ) call self_%activeNode%serializeValues(y(1:self_%propertyCountActive),self_%propertyTypeODE)
    ! If post-step processing returned a "continue" code, we do not need to reset ODE evolution, so set status back to success.
    if (postStepStatus == GSL_Continue) postStepStatus=GSL_Success
    return
  end subroutine standardPostStepProcessing

  subroutine standardStepErrorAnalyzer(time,timeEnd,currentPropertyValue,currentPropertyError,timeStep,stepStatus) bind(c)
    !!{
    Profiles ODE solver step sizes and errors.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double   , c_int
    use            :: Interface_GSL, only : GSL_Success
    implicit none
    real            (kind=c_double ), intent(in   ), dimension(*                         ) :: currentPropertyValue
    real            (kind=c_double ), intent(in   ), dimension(*                         ) :: currentPropertyError
    real            (kind=c_double ), intent(in   ), value                                 :: time                , timeEnd         , &
         &                                                                                    timeStep
    integer         (kind=c_int    ), intent(in   ), value                                 :: stepStatus
    double precision                                                                       :: scale               , scaledError     , &
         &                                                                                    scaledErrorMaximum
    integer         (c_size_t      )                                                       :: iProperty           , limitingProperty
    type            (varying_string)                                                       :: propertyName
    double precision                                , dimension(self_%propertyCountActive) :: scale               , rate

    ! Count steps.
    self_%countEvaluationsToSuccess=self_%countEvaluationsToSuccess+1_c_size_t
    ! If the step was not good, return immediately.
    if (stepStatus /= GSL_Success) return
    ! Stop the timer.
    call self_%stepTimer%stop()
    ! Find the property with the largest error (i.e. that which is limiting the step).
    scaledErrorMaximum=0.0d0
    limitingProperty  =-1_c_size_t
    scale             =+self_%odeToleranceAbsolute*    self_%propertyScalesActive(1:self_%propertyCountActive)  &
         &             +self_%odeToleranceRelative*abs(      currentPropertyValue(1:self_%propertyCountActive))
    do iProperty=1,self_%propertyCountActive
       scaledError=abs(currentPropertyError(iProperty))/scale(iProperty)
       if (scaledError > scaledErrorMaximum) then
          scaledErrorMaximum=scaledError
          limitingProperty  =iProperty
       end if
    end do
    ! Check that we found a limiting property.
    if (scaledErrorMaximum > 0.0d0) then
       ! Decode the step limiting property.
       propertyName=self_%activeNode%nameFromIndex(int(limitingProperty),self_%propertyTypeODE)
    else
       propertyName="unknown"
    end if
    ! Serialize rates.
    call self_%activeNode%serializeRates(rate,self_%propertyTypeODE)
    ! Profile the step.
    call self_%mergerTreeEvolveProfiler_%profile(self_%activeNode,time,timeStartSaved,timeEnd,timeStep,self_%countEvaluationsToSuccess,self_%interruptFirstFound,limitingProperty,propertyName,currentPropertyValue(1:self_%propertyCountActive),rate,scale,currentPropertyError(1:self_%propertyCountActive),self_%stepTimer%report())
    ! Reset the count of steps to success.
    self_%countEvaluationsToSuccess=0_c_size_t
    ! Restart the timer.
    call self_%stepTimer%start()
    return
  end subroutine standardStepErrorAnalyzer

  subroutine standardPromote(self,node)
    !!{
    Transfer the properties of {\normalfont \ttfamily node} to its parent node, then destroy it.
    !!}
    use :: Display        , only : displayMessage, displayVerbosity, verbosityLevelInfo
    use :: String_Handling, only : operator(//)
    implicit none
    class(mergerTreeNodeEvolverStandard), intent(inout)          :: self
    type (treeNode                     ), intent(inout), pointer :: node
    type (treeNode                     )               , pointer :: parentNode, satelliteNode, &
         &                                                          mergeeNode, hostNode
    type (varying_string               )                         :: message
    !$GLC attributes unused :: self

    ! Get pointer to parent node.
    parentNode => node%parent
    ! Display a message.
    if (displayVerbosity() >= verbosityLevelInfo) then
       message='Promoting node '
       message=message//node%index()//' to '//parentNode%index()
       call displayMessage(message,verbosityLevelInfo)
    end if
    ! Perform any processing necessary before this halo is promoted.
    call self%nodeOperator_%nodePromote(node)
    !![
    <eventHook name="nodePromotion">
     <callWith>node</callWith>
    </eventHook>
    !!]
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
          satelliteNode                => hostNode%lastSatellite ()
          satelliteNode%sibling        => node    %firstSatellite
       else
          hostNode     %firstSatellite => node    %firstSatellite
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
       mergeeNode             => mergeeNode%siblingMergee
    end do
    if (associated(parentNode%firstMergee)) then
       mergeeNode => parentNode%firstMergee
       do while (associated(mergeeNode%siblingMergee))
          mergeeNode => mergeeNode%siblingMergee
       end do
       mergeeNode%siblingMergee => node%firstMergee
    else
       parentNode%firstMergee   => node%firstMergee
    end if
    ! Nullify the child pointer for the parent.
    parentNode%firstChild => null()
    ! Destroy the node.
    call node%destroy()
    deallocate(node)
    return
  end subroutine standardPromote

  subroutine standardMerge(self,node)
    !!{
    Handles instances where {\normalfont \ttfamily node} is about to merge with its parent node.
    !!}
    use :: Display         , only : displayMessage    , displayVerbosity, verbosityLevelInfo
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: String_Handling , only : operator(//)
    !![
    <include directive="nodeMergerTask" type="moduleUse">
    !!]
    include 'events.node_mergers.process.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class    (mergerTreeNodeEvolverStandard), intent(inout) :: self
    type     (treeNode                     ), intent(inout) :: node
    class    (nodeComponentBasic           ), pointer       :: basic
    type     (varying_string               )                :: message
    character(len=7                        )                :: label

    ! Display a message.
    if (displayVerbosity() >= verbosityLevelInfo) then
       basic => node%basic()
       write (label,'(f7.4)') basic%time()
       message='Making node '
       message=message//node%index()//' a satellite in '//node%parent%index()//' at time '//trim(adjustl(label))//' Gyr'
       call displayMessage(message,verbosityLevelInfo)
    end if
    ! Call subroutines to perform any necessary processing prior to this node merger event.
    !![
    <include directive="nodeMergerTask" type="functionCall" functionType="void">
     <functionArgs>node</functionArgs>
    !!]
    include 'events.node_mergers.process.inc'
    !![
    </include>
    !!]
    call self%nodeOperator_        %nodesMerge(node)
    ! Process the merger.
    call self%mergerTreeNodeMerger_%process   (node)
    return
  end subroutine standardMerge

  subroutine standardNodeSubhaloPromotion(self,node,nodePromotion)
    !!{
    Promote a recently promoted subhalo to its new parent.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)          :: self
    type (treeNode), intent(inout), pointer :: node, nodePromotion
    !$GLC attributes unused :: nodePromotion
    
    select type (self)
    class is (mergerTreeNodeEvolverStandard)
       call self%promote(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine standardNodeSubhaloPromotion
