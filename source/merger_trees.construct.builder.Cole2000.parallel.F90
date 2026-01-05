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

  !!{
  An implementation of a merger tree builder using the algorithm of \cite{cole_hierarchical_2000} utilizing a recursive
  construction approach.  
  !!}
  
  !![
  <mergerTreeBuilder name="mergerTreeBuilderCole2000Parallel">
   <description>
     A merger tree builder class derived fro the \refClass{mergerTreeBuilderCole2000} class which implements a recursive tree
     construction algorithm utilizing OpenMP task parallelism.
   </description>
  </mergerTreeBuilder>
  !!]
  type, extends(mergerTreeBuilderCole2000) :: mergerTreeBuilderCole2000Parallel
     !!{
     A merger tree builder class using the algorithm of \cite{cole_hierarchical_2000}.
     !!}
     private
   contains
     procedure         :: build           => cole2000ParallelBuild
     procedure, nopass :: onBranch        => cole2000ParallelOnBranch
  end type mergerTreeBuilderCole2000Parallel

  interface mergerTreeBuilderCole2000Parallel
     !!{
     Constructors for the \refClass{mergerTreeBuilderCole2000Parallel} merger tree builder class.
     !!}
     module procedure cole2000ParallelConstructorParameters
     module procedure cole2000ParallelConstructorInternal
  end interface mergerTreeBuilderCole2000Parallel

  ! Sub-module scope variables used in tree building.
  !![
  <workaround type="gfortran" PR="110547" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=110547">
    <description>
      We use a pointer to self here rather than have self be passed to various methods (which are defined as "nopass" in the
      parent class) because otherwise gfortran calls the destructor of self on exit from these functions when using OpenMP
      task-based parallelism. This may be a compiler bug. Note that we use a separate pointer here that is not threadprivate as
      threadprivate variables can not be passed as arguments in untied OpenMP tasks.
    </description>
  </workaround>
  !!]  
  class(mergerTreeBuilderCole2000Parallel), pointer :: self__
  
contains

  function cole2000ParallelConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} merger tree building class which reads parameters from a provided parameter list.
    !!}
    implicit none
    type            (mergerTreeBuilderCole2000Parallel  )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (mergerTreeMassResolutionClass      ), pointer       :: mergerTreeMassResolution_
    class           (criticalOverdensityClass           ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), pointer       :: cosmologicalMassVariance_
    class           (mergerTreeBuildControllerClass     ), pointer       :: mergerTreeBuildController_
    double precision                                                     :: mergeProbability          , accretionLimit         , &
         &                                                                  redshiftMaximum           , toleranceResolutionSelf, &
         &                                                                  toleranceResolutionParent , toleranceTimeEarliest
    logical                                                              :: branchIntervalStep        , ignoreNoProgress

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>mergeProbability</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The largest probability of branching allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    </inputParameter>
    <inputParameter>
      <name>accretionLimit</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The largest fractional mass change due to subresolution accretion allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d5</defaultValue>
      <description>The highest redshift to which merger trees will be built in the \cite{cole_hierarchical_2000} method.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceTimeEarliest</name>
      <source>parameters</source>
      <defaultValue>2.0d-6</defaultValue>
      <description>The fractional tolerance used to judge if a branch is at the earliest allowed time in the tree.</description>
    </inputParameter>
    <inputParameter>
      <name>branchIntervalStep</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If {\normalfont \ttfamily false} use the original \cite{cole_hierarchical_2000} method to determine whether branching occurs in a timestep. If {\normalfont \ttfamily true} draw branching intervals from a negative exponential distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceResolutionSelf</name>
      <source>parameters</source>
      <defaultValue>1.0d-6</defaultValue>
      <description>The fractional tolerance in node mass at the resolution limit below which branch mis-orderings will be ignored.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceResolutionParent</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>The fractional tolerance in parent node mass at the resolution limit below which branch mis-orderings will be ignored.</description>
    </inputParameter>
    <inputParameter>
      <name>ignoreNoProgress</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, failure to make progress on a branch will be ignored (and the branch terminated).</description>
    </inputParameter>
    <objectBuilder class="mergerTreeMassResolution"       name="mergerTreeMassResolution_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"  source="parameters"/>
    <objectBuilder class="mergerTreeBuildController"      name="mergerTreeBuildController_" source="parameters"/>
    !!]
    self   =mergerTreeBuilderCole2000Parallel(                                                                                                             &
         &                                                                                                                   mergeProbability            , &
         &                                                                                                                   accretionLimit              , &
         &                                    cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum           )), &
         &                                                                                                                   toleranceTimeEarliest       , &
         &                                                                                                                   branchIntervalStep          , &
         &                                                                                                                   toleranceResolutionSelf     , &
         &                                                                                                                   toleranceResolutionParent   , &
         &                                                                                                                   ignoreNoProgress            , &
         &                                                                                                                   mergerTreeMassResolution_   , &
         &                                                                                                                   cosmologyFunctions_         , &
         &                                                                                                                   criticalOverdensity_        , &
         &                                                                                                                   cosmologicalMassVariance_   , &
         &                                                                                                                   mergerTreeBuildController_    &
         &                                   )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeMassResolution_"      />
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="mergerTreeBuildController_"     />
    !!]
    return
  end function cole2000ParallelConstructorParameters

  function cole2000ParallelConstructorInternal(mergeProbability,accretionLimit,timeEarliest,toleranceTimeEarliest,branchIntervalStep,toleranceResolutionSelf,toleranceResolutionParent,ignoreNoProgress,mergerTreeMassResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBuildController_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeBuilderCole2000Parallel  )                        :: self
    double precision                                     , intent(in   )         :: mergeProbability          , accretionLimit         , &
         &                                                                          timeEarliest              , toleranceResolutionSelf, &
         &                                                                          toleranceResolutionParent , toleranceTimeEarliest
    logical                                              , intent(in   )         :: branchIntervalStep        , ignoreNoProgress
    class           (mergerTreeMassResolutionClass      ), intent(in   ), target :: mergerTreeMassResolution_
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (mergerTreeBuildControllerClass     ), intent(in   ), target :: mergerTreeBuildController_
    !![
    <constructorAssign variables="mergeProbability, accretionLimit, timeEarliest, toleranceTimeEarliest, branchIntervalStep, toleranceResolutionSelf, toleranceResolutionParent, ignoreNoProgress, *mergerTreeMassResolution_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBuildController_"/>
    !!]

    ! Store maximum redshift.
    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeEarliest))
    ! Initialize state.
    self%branchingIntervalDistributionInitialized=.false.
    self%timeParameterIsMassDependent            = self%criticalOverdensity_     %isMassDependent      (     ) &
         &                                        .or.                                                         &
         &                                         self%cosmologicalMassVariance_%growthIsMassDependent(     )
    self%timeNow                                 = self%cosmologyFunctions_      %cosmicTime           (1.0d0)
    ! Validate parameters.
    if (self%accretionLimit >= 1.0d0) call Error_Report('accretionLimit < 1 required'//{introspection:location})
    ! Mark workers as not yet initialized.
    self%workersInitialized=.false.
    return
  end function cole2000ParallelConstructorInternal

  subroutine cole2000ParallelBuild(self,tree)
    !!{
    Build a merger tree.
    !!}
    !$ use OMP_Lib                 , only : OMP_Get_Max_Threads          , OMP_Get_Thread_Num
    use :: Display                 , only : displayReset                 , displayMagenta
    use :: Error                   , only : Error_Report                 , Warn
    use :: Events_Hooks            , only : eventsHooksFutureThread
    use :: Galacticus_Nodes        , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: ISO_Varying_String      , only : varying_string
    use :: Kind_Numbers            , only : kind_int8
    use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
    use :: Merger_Tree_Walkers     , only : mergerTreeWalkerClass
    implicit none
    class           (mergerTreeBuilderCole2000Parallel), intent(inout), target  :: self
    type            (mergerTree                       ), intent(inout), target  :: tree
    type            (treeNode                         )               , pointer :: node                    , nodeChild
    class           (nodeComponentBasic               )               , pointer :: basic                   , basicChild
    double precision                                                            :: deltaCritical           , deltaCriticalEarliest, &
         &                                                                         massResolution          , timeNodeBase         , &
         &                                                                         rootVarianceGrowthFactor, timeOfCollapse
    integer                                                                     :: countWorkers
    
    ! Begin construction.
    self__         => self
    self%nodeIndex =  1               ! Initialize the node index counter to unity.
    node           => tree%nodeBase   ! Point to the base node.
    basic          => node%basic   () ! Get the basic component of the node.
    if (.not.self%branchingIntervalDistributionInitialized.and.self%branchIntervalStep) then
       ! Note that we use a unit rate - we will scale the results to the actual rate required.
       self%branchingIntervalDistribution           =distributionFunction1DNegativeExponential(1.0d0)
       self%branchingIntervalDistributionInitialized=.true.
    end if
    ! Get the mass resolution for this tree.
    massResolution=self%mergerTreeMassResolution_%resolution(tree)
    ! Find the critical overdensity at the earliest time to which we will build this tree. We evaluate this here to ensure that
    ! δc(t) is already evaluated over the full range of epochs that we will require - this prevents any possible retabulation
    ! during tree construction which can potentially lead to misordering of branches. This is not the ideal solution - ideally,
    ! critical overdensity classes which rely on tabulation and which have to retabulate themselves should ensure that they simply
    ! expand their range without changing any of the previous computed values (as we do for expansion factor vs. time in the
    ! cosmology function class).
    deltaCriticalEarliest   =+self%criticalOverdensity_     %value          (time              =self%timeEarliest/2.0d0,mass=basic%mass(),node=node) &
         &                   *self%cosmologicalMassVariance_%rootVariance   (time              =self%timeNow           ,mass=basic%mass()          ) &
         &                   /self%cosmologicalMassVariance_%rootVariance   (time              =self%timeEarliest/2.0d0,mass=basic%mass()          )
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this class).
    timeNodeBase            =                                                                   basic%time        ()
    rootVarianceGrowthFactor=+self%cosmologicalMassVariance_%rootVariance   (time              =      timeNodeBase     ,mass=basic%mass()          ) &
         &                   /self%cosmologicalMassVariance_%rootVariance   (time              =self %timeNow          ,mass=basic%mass()          )
    deltaCritical           =+self%criticalOverdensity_     %value          (time              =basic%time        ()   ,mass=basic%mass(),node=node) &
         &                   /rootVarianceGrowthFactor
    call basic%timeSet(deltaCritical)
    ! Ensure that the time of collapse is pre-tabulated.
    timeOfCollapse          =+self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=deltaCritical          ,mass=basic%mass(),node=node)
    ! Create copies of objects needed for evolution.
    countWorkers   =1
    !$ countWorkers=OMP_Get_Max_Threads()
    if (.not.self%workersInitialized) then
       allocate(self%workers(0:countWorkers-1))
       self%workersInitialized=.true.
       do numberWorker=0,countWorkers-1
          call eventsHooksFutureThread(numberWorker)
          allocate(self%workers(numberWorker)%cosmologyFunctions_       ,mold=self%cosmologyFunctions_       )
          allocate(self%workers(numberWorker)%criticalOverdensity_      ,mold=self%criticalOverdensity_      )
          allocate(self%workers(numberWorker)%cosmologicalMassVariance_ ,mold=self%cosmologicalMassVariance_ )
          allocate(self%workers(numberWorker)%mergerTreeBuildController_,mold=self%mergerTreeBuildController_)
          if (self%branchIntervalStep) then
             allocate(self%workers(numberWorker)%branchingIntervalDistribution)
             !![
	     <referenceConstruct object="self%workers(numberWorker)%branchingIntervalDistribution" constructor="distributionFunction1DNegativeExponential(1.0d0)"/>
             !!]
          end if
          !$omp critical(mergerTreeBuilderCole2000ParallelDeepCopy)
          !![
	  <deepCopyReset variables="self%cosmologyFunctions_ self%criticalOverdensity_ self%cosmologicalMassVariance_ self%mergerTreeBuildController_"/>
	  <deepCopy source="self%cosmologyFunctions_"        destination="self%workers(numberWorker)%cosmologyFunctions_"       />
	  <deepCopy source="self%criticalOverdensity_"       destination="self%workers(numberWorker)%criticalOverdensity_"      />
	  <deepCopy source="self%cosmologicalMassVariance_"  destination="self%workers(numberWorker)%cosmologicalMassVariance_" />
	  <deepCopy source="self%mergerTreeBuildController_" destination="self%workers(numberWorker)%mergerTreeBuildController_"/>
	  <deepCopyFinalize variables="self%workers(numberWorker)%cosmologyFunctions_ self%workers(numberWorker)%criticalOverdensity_ self%workers(numberWorker)%cosmologicalMassVariance_ self%workers(numberWorker)%mergerTreeBuildController_"/>
          !!]
          !$omp end critical(mergerTreeBuilderCole2000ParallelDeepCopy)
          call eventsHooksFutureThread()
       end do
    end if
    ! Begin parallel tree build.
    !$omp parallel
    ! Set pointer to self from all threads.
    self_        => self
    ! Determine our worker number.
    numberWorker =  OMP_Get_Thread_Num()
    ! Copy the random number generator from the tree.
    allocate(self%workers(numberWorker)%randomNumberGenerator_,mold=tree%randomNumberGenerator_)
    !$omp critical(mergerTreeConstructBuilderCole2000DeepCopyReset)
    !![
    <deepCopyReset variables="tree%randomNumberGenerator_"/>
    <deepCopy source="tree%randomNumberGenerator_" destination="self%workers(numberWorker)%randomNumberGenerator_"/>
    <deepCopyFinalize variables="self%workers(numberWorker)%randomNumberGenerator_"/>
    !!]
    !$omp end critical(mergerTreeConstructBuilderCole2000DeepCopyReset)
    call self%workers(numberWorker)%randomNumberGenerator_%seedSet(seed=tree%index+numberWorker,offset=.true.)
    !$omp barrier
    !$omp single
    ! Begin build of the main branch - other branches will be built recursively.
    !$omp taskgroup
    !$omp task untied
    call self%buildBranch(tree,massResolution,node)
    !$omp end task
    !$omp end taskgroup
    ! Convert w to time (and test for well-ordering) along the main branch. Other branches will be converted recursively.
    !$omp taskgroup
    !$omp task untied
    call convertTimeBranch(tree%nodeBase)
    !$omp end task
    !$omp end taskgroup
    ! Fix the time of the base node to that required by the tree.
    basic => tree%nodeBase%basic()
    call basic%timeSet(timeNodeBase)
    ! Check for mis-ordering of the base node and its child node(s). This can happen because we force the time of the base node to
    ! be precisely the base time, but for other nodes the time is computed by inverting the w(t)=δ_crit(t)/D(t)
    ! relation. Numerical inaccuracies in the inversion can lead to small mis-ordering in the tree times.
    if (associated(tree%nodeBase%firstChild)) then
       basicChild => tree%nodeBase%firstChild%basic()
       if (basic%time() <= basicChild%time()) then
          ! Base node is mis-ordered. Simply shift and child nodes to be slightly earlier. If this leads to mis-ordering of those
          ! child nodes it will be detected below.
          call Warn(displayMagenta()//'WARNING:'//displayReset()//' tree is not well-ordered at base node - fixing')
          nodeChild => tree%nodeBase%firstChild
          do while (associated(nodeChild))
             basicChild => nodeChild%basic()
             call basicChild%timeSet(basic%time())
             nodeChild => nodeChild%sibling
          end do
       end if
    end if
    !$omp end single
    !![
    <objectDestructor name="self%workers(numberWorker)%randomNumberGenerator_"/>
    !!]
    !$omp end parallel
    return
  
  contains

    recursive subroutine convertTimeBranch(nodeTip)
      !!{
      Convert from critical overdensity to time along a branch.
      !!}
      implicit none
      type(treeNode), intent(inout), target  :: nodeTip
      type(treeNode)               , pointer :: nodeCurrent
      
      nodeCurrent => nodeTip
      do while (associated(nodeCurrent))
         if (associated(nodeCurrent%sibling)) then
            ! Recursively convert the sibling branch.
            !$omp task untied
            call convertTimeBranch(nodeCurrent%sibling)
            !$omp end task
         end if
         ! Perform conversion and checks for this node.
         call self%convertTimeNode(                    nodeCurrent)
         call self%checkOrderNode (tree,massResolution,nodeCurrent)
         ! Move to the next child node (if available).
         if (associated(nodeCurrent)) nodeCurrent => nodeCurrent%firstChild
      end do
      return
    end subroutine convertTimeBranch

  end subroutine cole2000ParallelBuild

  recursive subroutine cole2000ParallelOnBranch(tree,massResolution,node)
    !!{
    Act on branching.
    !!}
    implicit none
    type            (mergerTree), intent(in   )          :: tree
    double precision            , intent(in   )          :: massResolution
    type            (treeNode  ), intent(inout), pointer :: node
    
    !$omp task untied
    call self__%buildBranch(tree,massResolution,node)
    !$omp end task
    return
  end subroutine cole2000ParallelOnBranch
