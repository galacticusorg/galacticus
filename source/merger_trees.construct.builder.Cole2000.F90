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
  An implementation of a merger tree builder using the algorithm of \cite{cole_hierarchical_2000}.
  !!}

  use :: Cosmological_Density_Field        , only : cosmologicalMassVarianceClass            , criticalOverdensityClass
  use :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use :: Merger_Tree_Branching             , only : mergerTreeBranchingProbabilityClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  use :: Statistics_Distributions          , only : distributionFunction1DNegativeExponential
  use :: Merger_Tree_Build_Controllers     , only : mergerTreeBuildControllerClass
  use :: Merger_Tree_Walkers               , only : mergerTreeWalkerTreeConstruction
  use :: Numerical_Random_Numbers          , only : randomNumberGeneratorClass
  use :: Kind_Numbers                      , only : kind_int8

  ! Structure used to hold worker copies of objects.
  type :: mergerTreeBuilderCole2000Worker
     class(cosmologyFunctionsClass                 ), pointer :: cosmologyFunctions_           => null()
     class(criticalOverdensityClass                ), pointer :: criticalOverdensity_          => null()
     class(cosmologicalMassVarianceClass           ), pointer :: cosmologicalMassVariance_     => null()
     class(mergerTreeBuildControllerClass          ), pointer :: mergerTreeBuildController_    => null()
     class(randomNumberGeneratorClass              ), pointer :: randomNumberGenerator_        => null()
     type(distributionFunction1DNegativeExponential), pointer :: branchingIntervalDistribution => null()
  end type mergerTreeBuilderCole2000Worker
  
  !![
  <mergerTreeBuilder name="mergerTreeBuilderCole2000">
   <description>
    A merger tree builder class which uses the algorithm described by \cite{cole_hierarchical_2000} (with minor modifications
    described below). This action of this algorithm is controlled by the following parameters:
    \begin{description}
     \item [{\normalfont \ttfamily [mergeProbability]}] The maximum probability for a binary merger allowed in a single
     timestep. This allows the probability to be kept small, such the the probability for multiple mergers within a single
     timestep is small.
     \item [{\normalfont \ttfamily [accretionLimit]}] The maximum fractional change in mass due to sub-resolution accretion
     allowed in any given timestep when building the tree.
     \item [{\normalfont \ttfamily [redshiftMaximum]}] The highest redshift to which the tree should be built. Any branch
     reaching this redshift will be terminated. Typically this should be set to a high value such that branches terminate when
     the resolution limit it reached, but specifying a maximum redshift can be useful in some situations.
     \item [{\normalfont \ttfamily [branchIntervalStep]}] If {\normalfont \ttfamily true}, instead of limiting each time step
     such that the probability of branching is less than {\normalfont \ttfamily mergerTreeBuildCole2000MergeProbability}, the
     interval to the next branching event will be drawn from a negative exponential with the appropriate rate. If this exceeds
     the maximum allowed timestep based on other considerations (e.g. the accretion limit), no branching occurs, and the
     timestep proceeds\footnote{Note that we do not have to concern ourselves in the subsequent timestep with the fact that no
     branching occurred in the previous timestep because of the memorylessness nature of the negative exponential
     distribution. That is, the distribution of branching intervals conditioned on the fact that no branching occurred in the
     previous timestep, is just the same negative exponential distribution.}. If the interval is less than the maximum allowed
     timestep, branching occurs at that point. In the regime of high branching rates (which occur when the branch being grown
     is far above the mass resolution), this approach allows for larger timesteps to be taken.
    \end{description}
    The minimum halo mass that the algorithm will follow is determined by the selection merger tree building mass resolution
    method (see \refPhysics{mergerTreeMassResolution}). Mass accretion below this scale is treated as smooth accretion
    and branches are truncated once they fall below this mass.
    
    In the original \cite{cole_hierarchical_2000}, when a branch split occurred masses, $M_2$ and $M_3$, of the two new halos
    were selected by first drawing the mass $M_2$ from the branching distribution function in the range $M_\mathrm{res}$ to
    $M_1/2$ (where $M_1$ is the mass of the parent halo, and $M_\mathrm{res}$ is the mass resolution being used for the tree),
    and then setting
    \begin{equation}
      M_3 = M_1 (1-F) - M_2
    \end{equation}
    where $F$ is the fraction of the parent halo mass gained through sub-resolution accretion in this timestep. As the
    sub-resolution accretion is removed entirely from the mass $M_3$ and not from $M_2$ this can lead to an asymmetry in
    progenitor mass functions close to $M_1/2$. Therefore, we instead set the progenitor masses by first drawing a mass
    $M_2^\prime$ from the mass branching distribution function and then setting
    \begin{eqnarray}
      M_2 &amp;=&amp; M_2^\prime (1-F), \nonumber \\
      M_3 &amp;=&amp; (M_1 - M_2^\prime) (1-F),
    \end{eqnarray}
    which ensures a symmetric treatment of subresolution accretion close to $M_1/2$.   
   </description>
   <deepCopy>
     <ignore variables="workers"/>
     <setTo    variables="workersInitialized" value=".false."/>
   </deepCopy>
  </mergerTreeBuilder>
  !!]
  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderCole2000
     !!{
     A merger tree builder class using the algorithm of \cite{cole_hierarchical_2000}.
     !!}
     private
     class           (cosmologyFunctionsClass                  ), pointer                   :: cosmologyFunctions_                      => null()
     class           (mergerTreeMassResolutionClass            ), pointer                   :: mergerTreeMassResolution_                => null()
     class           (criticalOverdensityClass                 ), pointer                   :: criticalOverdensity_                     => null()
     class           (mergerTreeBranchingProbabilityClass      ), pointer                   :: mergerTreeBranchingProbability_          => null()
     class           (cosmologicalMassVarianceClass            ), pointer                   :: cosmologicalMassVariance_                => null()
     class           (mergerTreeBuildControllerClass           ), pointer                   :: mergerTreeBuildController_               => null()
     type            (mergerTreeBuilderCole2000Worker          ), allocatable, dimension(:) :: workers
     logical                                                                                :: timeParameterIsMassDependent
     ! Node index counter.
     integer         (kind=kind_int8                           )                            :: nodeIndex
     ! Variables controlling merger tree accuracy.
     double precision                                                                       :: accretionLimit                                    , timeEarliest                    , &
          &                                                                                    mergeProbability                                  , timeNow                         , &
          &                                                                                    redshiftMaximum                                   , toleranceTimeEarliest
     logical                                                                                :: ignoreNoProgress                                  , ignoreWellOrdering
     ! Random number sequence variables
     logical                                                                                :: branchIntervalStep
     ! Interval distribution.
     logical                                                                                :: branchingIntervalDistributionInitialized          , workersInitialized       =.false.
     type            (distributionFunction1DNegativeExponential)                            :: branchingIntervalDistribution
     ! Tolerances for behavior close to the resolution limit.
     double precision                                                                       :: toleranceResolutionSelf                           , toleranceResolutionParent
   contains
     !![
     <methods>
       <method description="Build a branch of the merger tree."                    method="buildBranch"                   />
       <method description="Set the critical overdensity object."                  method="criticalOverdensityUpdate"     />
       <method description="Convert from critical overdensity to time for a node." method="convertTimeNode"               />
       <method description="Convert from time to critical overdensity for a node." method="convertCriticalOverdensityNode"/>
       <method description="Check well-ordering in time for a node."               method="checkOrderNode"                />
       <method description="Perform any required actions on branching."            method="onBranch"                      />
     </methods>
     !!]
     final             ::                                   cole2000Destructor
     procedure         :: build                          => cole2000Build
     procedure         :: timeEarliestSet                => cole2000TimeEarliestSet
     procedure, nopass :: buildBranch                    => cole2000BuildBranch
     procedure, nopass :: onBranch                       => cole2000OnBranch
     procedure, nopass :: criticalOverdensityUpdate      => cole2000CriticalOverdensityUpdate
     procedure, nopass :: convertTimeNode                => cole2000ConvertTimeNode
     procedure         :: convertCriticalOverdensityNode => cole2000ConvertCriticalOverdensityNode
     procedure, nopass :: checkOrderNode                 => cole2000CheckOrderNode
     procedure         :: stateStore                     => cole2000StateStore
     procedure         :: stateRestore                   => cole2000StateRestore
  end type mergerTreeBuilderCole2000

  interface mergerTreeBuilderCole2000
     !!{
     Constructors for the \refClass{mergerTreeBuilderCole2000} merger tree builder class.
     !!}
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerTreeBuilderCole2000

  ! Sub-module scope variables used in tree building.
  !![
  <workaround type="gfortran" PR="110547" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=110547">
    <description>
      We use a pointer to self here rather than have self be passed to various methods (which are defined as "nopass" above)
      because otherwise gfortran calls the destructor of self on exit from these functions when using OpenMP task-based
      parallelism. This may be a compiler bug.
    </description>
  </workaround>
  <workaround type="gfortran" PR="110548" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=110548">
    <description>
      We use an allocatable treeWalker here. Some subclasses will not need a treeWalker. Ideally this would therefore be passed
      to functions as an optional argument, but current optional arguments used in OpenMP tasks cause segfaults. Therefore, we
      instead use the allocation status of this variable as an indication of whether or not a treeWalker is needed.
    </description>
  </workaround>
  !!]  
  class  (mergerTreeBuilderCole2000       ), pointer     :: self_
  integer                                                :: numberWorker
  type   (mergerTreeWalkerTreeConstruction), allocatable :: treeWalker_
  !$omp threadprivate(self_,numberWorker,treeWalker_)

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} merger tree building class which reads parameters from a provided parameter list.
    !!}
    implicit none
    type            (mergerTreeBuilderCole2000          )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (mergerTreeMassResolutionClass      ), pointer       :: mergerTreeMassResolution_
    class           (criticalOverdensityClass           ), pointer       :: criticalOverdensity_
    class           (mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class           (cosmologicalMassVarianceClass      ), pointer       :: cosmologicalMassVariance_
    class           (mergerTreeBuildControllerClass     ), pointer       :: mergerTreeBuildController_
    double precision                                                     :: mergeProbability               , accretionLimit         , &
         &                                                                  redshiftMaximum                , toleranceResolutionSelf, &
         &                                                                  toleranceResolutionParent      , toleranceTimeEarliest
    logical                                                              :: branchIntervalStep             , ignoreNoProgress       , &
         &                                                                  ignoreWellOrdering

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
    <inputParameter>
      <name>ignoreWellOrdering</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, non-well-ordered tree branches are pruned away instead of causing errors..</description>
    </inputParameter>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"       name="mergerTreeMassResolution_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="mergerTreeBuildController"      name="mergerTreeBuildController_"      source="parameters"/>
    !!]
    self   =mergerTreeBuilderCole2000(                                                                                                                  &
         &                                                                                                           mergeProbability                 , &
         &                                                                                                           accretionLimit                   , &
         &                            cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum                )), &
         &                                                                                                           toleranceTimeEarliest            , &
         &                                                                                                           branchIntervalStep               , &
         &                                                                                                           toleranceResolutionSelf          , &
         &                                                                                                           toleranceResolutionParent        , &
         &                                                                                                           ignoreNoProgress                 , &
         &                                                                                                           ignoreWellOrdering               , &
         &                                                                                                           mergerTreeBranchingProbability_  , &
         &                                                                                                           mergerTreeMassResolution_        , &
         &                                                                                                           cosmologyFunctions_              , &
         &                                                                                                           criticalOverdensity_             , &
         &                                                                                                           cosmologicalMassVariance_        , &
         &                                                                                                           mergerTreeBuildController_         &
         &                           )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    <objectDestructor name="mergerTreeMassResolution_"      />
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="mergerTreeBuildController_"     />
    !!]
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(mergeProbability,accretionLimit,timeEarliest,toleranceTimeEarliest,branchIntervalStep,toleranceResolutionSelf,toleranceResolutionParent,ignoreNoProgress,ignoreWellOrdering,mergerTreeBranchingProbability_,mergerTreeMassResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBuildController_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeBuilderCole2000          )                        :: self
    double precision                                     , intent(in   )         :: mergeProbability               , accretionLimit         , &
         &                                                                          timeEarliest                   , toleranceResolutionSelf, &
         &                                                                          toleranceResolutionParent      , toleranceTimeEarliest
    logical                                              , intent(in   )         :: branchIntervalStep             , ignoreNoProgress       , &
         &                                                                          ignoreWellOrdering
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    class           (mergerTreeMassResolutionClass      ), intent(in   ), target :: mergerTreeMassResolution_
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (mergerTreeBuildControllerClass     ), intent(in   ), target :: mergerTreeBuildController_
    !![
    <constructorAssign variables="mergeProbability, accretionLimit, timeEarliest, toleranceTimeEarliest, branchIntervalStep, toleranceResolutionSelf, toleranceResolutionParent, ignoreNoProgress, ignoreWellOrdering, *mergerTreeBranchingProbability_, *mergerTreeMassResolution_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBuildController_"/>
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
    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuilderCole2000} merger tree builder class.
    !!}
    implicit none
    type   (mergerTreeBuilderCole2000), intent(inout) :: self
    integer                                           :: i
    
    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%mergerTreeMassResolution_"      />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%mergerTreeBuildController_"     />
    !!]
    if (self%workersInitialized) then
       do i=lbound(self%workers,dim=1),ubound(self%workers,dim=1)
          !![
	  <objectDestructor name="self%workers(i)%cosmologyFunctions_"          />
	  <objectDestructor name="self%workers(i)%criticalOverdensity_"         />
	  <objectDestructor name="self%workers(i)%cosmologicalMassVariance_"    />
	  <objectDestructor name="self%workers(i)%mergerTreeBuildController_"   />
	  <objectDestructor name="self%workers(i)%branchingIntervalDistribution"/>
          !!]
       end do
    end if
    return
  end subroutine cole2000Destructor

  subroutine cole2000Build(self,tree)
    !!{
    Build a merger tree.
    !!}
    use :: Display            , only : displayReset                 , displayMagenta
    use :: Error              , only : Warn
    use :: Events_Hooks       , only : eventsHooksFutureThread
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: ISO_Varying_String , only : varying_string
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class           (mergerTreeBuilderCole2000    ), intent(inout), target :: self
    type            (mergerTree                   ), intent(inout), target :: tree
    type            (treeNode                     ), pointer               :: node                 , nodeChild     , &
         &                                                                    nodeWork
    class           (nodeComponentBasic           ), pointer               :: basic                , basicChild
    type            (mergerTreeWalkerIsolatedNodes)                        :: treeWalkerIsolated
    double precision                                                       :: timeNodeBase         , massResolution, &
         &                                                                    deltaCriticalEarliest        
      
    ! Begin construction.
    self_          => self
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
    deltaCriticalEarliest=+self%criticalOverdensity_     %value       (time=self%timeEarliest/2.0d0,mass=basic%mass(),node=node) &
         &                *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow           ,mass=basic%mass()          ) &
         &                /self%cosmologicalMassVariance_%rootVariance(time=self%timeEarliest/2.0d0,mass=basic%mass()          )
    ! Convert time for all existing nodes in the tree to critical overdensity (which we use as a time coordinate in this class).
    timeNodeBase      =basic%time()
    treeWalkerIsolated=mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalkerIsolated%next(nodeWork))
       call self%convertCriticalOverdensityNode(nodeWork)
       ! Ensure that the node index is set to the largest value present in the tree so far. This ensures that node indices are not
       ! repeated.
       self%nodeIndex=max(self%nodeIndex,nodeWork%index())
    end do
    ! Determine our worker number.
    numberWorker=0
    ! Create copies of objects needed for evolution.
    if (.not.self%workersInitialized) then
       allocate(self%workers(0:numberWorker))
       self%workersInitialized=.true.
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
    end if
    ! Copy the random number generator from the tree.
    allocate(self%workers(numberWorker)%randomNumberGenerator_,mold=tree%randomNumberGenerator_)
    !$omp critical(mergerTreeBuilderCole2000ParallelDeepCopy)
    !![
    <deepCopyReset variables="tree%randomNumberGenerator_"/>
    <deepCopy source="tree%randomNumberGenerator_" destination="self%workers(numberWorker)%randomNumberGenerator_"/>
    <deepCopyFinalize variables="self%workers(numberWorker)%randomNumberGenerator_"/>
    !!]
    !$omp end critical(mergerTreeBuilderCole2000ParallelDeepCopy)
    call self%workers(numberWorker)%randomNumberGenerator_%seedSet(seed=tree%index+numberWorker,offset=.true.)
    ! Begin tree build loop.
    allocate(treeWalker_)
    treeWalker_=mergerTreeWalkerTreeConstruction(tree)
    do while (treeWalker_%next(node))
       call self%buildBranch(tree,massResolution,node)
    end do
    deallocate(treeWalker_)
    ! Walk the tree and convert w to time.
    treeWalkerIsolated=mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalkerIsolated%next(node))
       call self%convertTimeNode(node)
    end do
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
    ! Check for well-ordering in time.
    treeWalkerIsolated=mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalkerIsolated%next(node))
       call self%checkOrderNode(tree,massResolution,node)
    end do
    ! Clean up.
    !![
    <objectDestructor name="self%workers(numberWorker)%randomNumberGenerator_"/>
    !!]
    return
  end subroutine cole2000Build
  
  recursive subroutine cole2000BuildBranch(tree,massResolution,nodeTip)
    !!{
    Build a single branch of the merger tree.
    !!}
    use :: Display              , only : displayReset                 , displayMagenta               , displayGreen, displayYellow, &
         &                               displayBlue
    use :: Galacticus_Nodes     , only : mergerTree                   , nodeComponentBasic
    use :: Merger_Tree_Branching, only : mergerTreeBranchingBoundLower, mergerTreeBranchingBoundUpper
    use :: Merger_Tree_Walkers  , only : mergerTreeWalkerClass
    use :: Numerical_Comparison , only : Values_Agree
    use :: String_Handling      , only : stringXMLFormat
    implicit none
    type            (mergerTree                         ), intent(in   )           :: tree
    double precision                                     , intent(in   )           :: massResolution
    type            (treeNode                           ), intent(inout), pointer  :: nodeTip
    type            (treeNode                           )               , pointer  :: nodeNew1                              , nodeNew2                   , &
         &                                                                            nodeCurrent
    class           (nodeComponentBasic                 )               , pointer  :: basicNew1                             , basicNew2                  , &
         &                                                                            basic_
    class           (mergerTreeBranchingProbabilityClass)               , pointer  :: mergerTreeBranchingProbability_
    double precision                                     , parameter               :: toleranceDeltaCritical         =1.0d-6
    double precision                                     , parameter               :: toleranceTime                  =1.0d-6
    integer                                              , parameter               :: limitNoProgress                =10000
    integer                                                                        :: countNoProgress
    double precision                                                               :: accretionFraction                     , branchingProbability       , &
         &                                                                            deltaCritical1                        , deltaCritical2             , &
         &                                                                            deltaW                                , nodeMass1                  , &
         &                                                                            nodeMass2                             , uniformRandom              , &
         &                                                                            accretionFractionCumulative           , branchMassCurrent          , &
         &                                                                            branchDeltaCriticalCurrent            , branchingInterval          , &
         &                                                                            branchingProbabilityRate              , branchingIntervalScaleFree , &
         &                                                                            deltaWAccretionLimit                  , deltaWEarliestTime         , &
         &                                                                            collapseTimeTruncate                  , time                       , &
         &                                                                            deltaWController                      , branchDeltaCriticalPrevious, &
         &                                                                            branchMassPrevious                    , collapseTime_              , &
         &                                                                            rootVarianceGrowthFactor_             , deltaCritical_
    logical                                                                        :: doBranch                              , branchIsDone               , &
         &                                                                            snapAccretionFraction                 , snapEarliestTime           , &
         &                                                                            controlLimited                        , controlInsertNode

    nodeCurrent => nodeTip
    do while (associated(nodeCurrent))
       ! Apply control.
       if (allocated(treeWalker_)) then
          if (.not.self_%workers(numberWorker)%mergerTreeBuildController_%control(nodeCurrent,treeWalker_)) return
       else
          if (.not.self_%workers(numberWorker)%mergerTreeBuildController_%control(nodeCurrent            )) return
       end if
       ! Get the basic component of the node.
       basic_                       => nodeCurrent%basic()
       ! Initialize the state for this branch.
       accretionFractionCumulative =  0.0d0
       branchMassCurrent           =  basic_      %mass ()
       branchDeltaCriticalCurrent  =  basic_      %time ()
       ! Evolve the branch until mass falls below the resolution limit, the earliest time is reached, or the branch ends.
       branchIsDone                =  .false.
       countNoProgress             =  0
       branchDeltaCriticalPrevious =  -huge(0.0d0)
       branchMassPrevious          =  -huge(0.0d0)
       do while (.not.branchIsDone)
          ! Get the branching probability object to use for this step.
          mergerTreeBranchingProbability_ => self_%workers(numberWorker)%mergerTreeBuildController_%branchingProbabilityObject(nodeCurrent)
          ! Get the growth factor in the root variance at the mass of the current branch.
          time                     =+self_%workers(numberWorker)%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=      branchDeltaCriticalCurrent,mass=branchMassCurrent,node=nodeCurrent)
          rootVarianceGrowthFactor_=+self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance  (time               =      time                      ,mass=branchMassCurrent                 ) &
               &                    /self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance  (time               =self_%timeNow                   ,mass=branchMassCurrent                 )
          ! Check for progress.
          if (branchDeltaCriticalCurrent == branchDeltaCriticalPrevious .and. branchMassCurrent == branchMassPrevious) then
             countNoProgress            =countNoProgress+1
          else
             countNoProgress            =0
             branchDeltaCriticalPrevious=branchDeltaCriticalCurrent
             branchMassPrevious         =branchMassCurrent
          end if
          if (countNoProgress > limitNoProgress) then
             block
               type     (varying_string) :: message
               character(len=20        ) :: label
               message='branch is making no progress'
               if (time < self_%timeEarliest*1.01d0) then
                  write (label,'(e12.6)') self_%toleranceTimeEarliest
                  message=message                                                                                                                   //char(10)// &
                       &  displayGreen()//'   HELP:'//displayReset()//' branch is within 1% of the imposed earliest time'                           //char(10)// &
                       &                  '          - this may be a tolerance issue'                                                               //char(10)// &
                       &                  '          - try increasing the value of the tolerance parameter, currently (see highlighted line below):'//char(10)// &
                       & stringXMLFormat(                                                                                                                        &
                       &                  '<mergerTreeBuilder value="'//char(self_%objectType(short=.true.))//'">'                                            // &
                       &                  '**B<toleranceTimeEarliest value="'//trim(adjustl(label))//'"/>**C'                                                 // &
                       &                  '</mergerTreeBuilder>'                                                                                                 &
                       &                 )
               end if
               call Error_Report(message//{introspection:location})
             end block
          end if
          ! Process the branch.
          if     (                                                                               &
               &     branchMassCurrent <= massResolution                                         &
               &  .or.                                                                           &
               &     time              <  self_%timeEarliest*(1.0d0+self_%toleranceTimeEarliest) &
               &  .or.                                                                           &
               &   (                                                                             &
               &     countNoProgress   == limitNoProgress                                        &
               &    .and.                                                                        &
               &     self_%ignoreNoProgress                                                      &
               &   )                                                                             &
               & ) then
             ! Branch should be terminated. If we have any accumulated accretion, terminate the branch with a final node.
             if (accretionFractionCumulative > 0.0d0) then
                !$omp atomic
                self_%nodeIndex =  self_%nodeIndex+1
                nodeNew1        => treeNode(self_%nodeIndex,nodeCurrent%hostTree)
                basicNew1       => nodeNew1%basic(autoCreate=.true.)
                ! Compute new mass accounting for sub-resolution accretion.
                nodeMass1       =  basic_%mass()*(1.0d0-accretionFractionCumulative)
                call basicNew1%massSet(nodeMass1     )
                ! Compute the critical overdensity corresponding to this new node.
                deltaCritical1=self_%criticalOverdensityUpdate(branchDeltaCriticalCurrent,branchMassCurrent,nodeMass1,nodeNew1)
                call basicNew1%timeSet(deltaCritical1)
                ! Inform the build controller of this new node.
                call self_%workers(numberWorker)%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeNew1)
                ! Create links from old to new node and vice-versa.
                nodeCurrent%firstChild => nodeNew1
                nodeNew1%parent => nodeCurrent
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it).
                nodeCurrent            => nodeNew1
             end if
             ! Flag that the branch is done.
             branchIsDone=.true.
          else
             ! Find branching probability rate per unit δw.
             branchingProbabilityRate=+mergerTreeBranchingProbability_%probabilityBound     (branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,mergerTreeBranchingBoundUpper,nodeCurrent) &
                  &                   *rootVarianceGrowthFactor_
             ! Find accretion rate.
             accretionFraction       =+mergerTreeBranchingProbability_%fractionSubresolution(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution                              ,nodeCurrent) &
                  &                   *rootVarianceGrowthFactor_
             ! A negative accretion fraction indicates that the node is so close to the resolution limit that
             ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
             ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
             if (accretionFraction < 0.0d0) then
                ! Terminate the branch with a final node.
                !$omp atomic
                self_%nodeIndex    =  self_%nodeIndex+1
                nodeNew1           => treeNode      (self_%nodeIndex        ,nodeCurrent%hostTree)
                basicNew1          => nodeNew1%basic(autoCreate=.true.     )
                ! Create a node at the mass resolution.
                nodeMass1          =  massResolution
                ! Compute critical overdensity for this new node.
                deltaCritical1      =self_%criticalOverdensityUpdate(branchDeltaCriticalCurrent,branchMassCurrent,nodeMass1,nodeNew1)
                ! Ensure critical overdensity exceeds that of the current node.
                collapseTime_       =self_%workers(numberWorker)%criticalOverdensity_%timeOfCollapse(criticalOverdensity=branchDeltaCriticalCurrent,mass=branchMassCurrent,node=nodeNew1)
                collapseTimeTruncate=self_%workers(numberWorker)%criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical1            ,mass=nodeMass1        ,node=nodeNew1)
                if (collapseTimeTruncate > collapseTime_*(1.0d0+toleranceTime)) then
                   call Error_Report('truncating to resolution, but resolution node exists after parent'//{introspection:location})
                else
                   do while (collapseTimeTruncate > collapseTime_*(1.0d0-toleranceTime))
                      deltaCritical1      =deltaCritical1*(1.0d0+toleranceTime)
                      collapseTimeTruncate=self_%workers(numberWorker)%criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical1,mass=nodeMass1,node=nodeNew1)
                   end do
                end if
                ! Set properties of the new node.
                call basicNew1%massSet(nodeMass1     )
                call basicNew1%timeSet(deltaCritical1)
                ! Inform the build controller of this new node.
                call self_%workers(numberWorker)%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeNew1)
                ! Create links from old to new node and vice-versa.
                nodeCurrent%firstChild => nodeNew1
                nodeNew1   %parent     => nodeCurrent
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it), and flag that the branch is done.
                nodeCurrent            => nodeNew1
                branchIsDone           =  .true.
             else
                ! Finding maximum allowed step in w. Limit based on branching rate only if we are using the original Cole et
                ! al. (2000) algorithm.
                deltaW               =mergerTreeBranchingProbability_%stepMaximum(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution)/rootVarianceGrowthFactor_
                snapAccretionFraction=.false.
                if (accretionFraction > 0.0d0) then
                   deltaWAccretionLimit=(self_%accretionLimit-accretionFractionCumulative)/accretionFraction
                   if (deltaWAccretionLimit <= deltaW) then
                      deltaW               =deltaWAccretionLimit
                      snapAccretionFraction=.true.
                   end if
                end if
                if     (                                  &
                     &   branchingProbabilityRate > 0.0d0 &
                     &  .and.                             &
                     &   .not.self_%branchIntervalStep    &
                     & ) then
                   if (self_%mergeProbability/branchingProbabilityRate < deltaW) then
                      ! Timestep is limited by branching rate. Reduce the timestep to the allowed
                      ! size and unset the flag to snap the accretion fraction to its maximum
                      ! allowed value (since we won't reach that value with this new, reduced
                      ! timestep).
                      deltaW               =self_%mergeProbability/branchingProbabilityRate
                      snapAccretionFraction=.false.
                   end if
                end if
                ! Limit the timestep so that the maximum allowed time is not exceeded.
                deltaWEarliestTime=+self_%workers(numberWorker)%criticalOverdensity_     %value       (time=self_%timeEarliest,mass=branchMassCurrent         ,node=nodeCurrent) &
                     &             *self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance(time=self_%timeNow     ,mass=branchMassCurrent                          ) &
                     &             /self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance(time=self_%timeEarliest,mass=branchMassCurrent                          ) &
                     &             -                                                                                                branchDeltaCriticalCurrent
                if (deltaWEarliestTime < deltaW) then
                   deltaW               =deltaWEarliestTime
                   snapEarliestTime     =.true.
                   snapAccretionFraction=.false.
                else
                   snapEarliestTime     =.false.
                end if
                ! Allow the build controller to set a minimum step.
                deltaW=max(                                                                                                                                              &
                     &                                                             deltaW                                                                              , &
                     &     +self_%workers(numberWorker)%mergerTreeBuildController_%timeMinimum               (nodeCurrent,branchMassCurrent,branchDeltaCriticalCurrent)  &
                     &     -                                                       branchDeltaCriticalCurrent                                                            &
                     &    )
                ! Limit the timestep according to the build controller.
                deltaWController=+self_%workers(numberWorker)%mergerTreeBuildController_%timeMaximum               (nodeCurrent,branchMassCurrent,branchDeltaCriticalCurrent,self_%timeNow,controlInsertNode) &
                     &           -                                                       branchDeltaCriticalCurrent
                if (deltaWController < deltaW) then
                   deltaW               =deltaWController
                   controlLimited       =.true.
                   snapAccretionFraction=.false.
                   snapEarliestTime     =.false.
                else
                   controlLimited       =.false.
                end if
                ! Scale values to the determined timestep.
                branchingProbability       =0.0d0
                if (.not.self_%branchIntervalStep)                           &
                     & branchingProbability=branchingProbabilityRate*deltaW
                accretionFraction          =accretionFraction       *deltaW
                ! Accretion fraction must be less than unity. Reduce timestep (and branching
                ! probability and accretion fraction) by factors of two until this condition is
                ! satisfied.
                do while (accretionFraction+accretionFractionCumulative >= 1.0d0)
                   if (.not.self_%branchIntervalStep)                      &
                        & branchingProbability=branchingProbability*0.5d0
                   accretionFraction          =accretionFraction   *0.5d0
                   deltaW                     =deltaW              *0.5d0
                   snapAccretionFraction      =.false.
                   snapEarliestTime           =.false.
                   controlLimited             =.false.
                end do
                ! Decide if a branching occurs.
                if (self_%branchIntervalStep) then
                   ! In this case we draw intervals between branching events from a negative
                   ! exponential distribution.
                   if (branchingProbabilityRate > 1.0d-100) then
                      branchingIntervalScaleFree=0.0d0
                      do while (branchingIntervalScaleFree <= 0.0d0)
                         branchingIntervalScaleFree=self_%workers(numberWorker)%branchingIntervalDistribution%sample(randomNumberGenerator_=self_%workers(numberWorker)%randomNumberGenerator_)
                      end do
                      branchingInterval=branchingIntervalScaleFree/branchingProbabilityRate
                      ! Based on the upper bound on the rate, check if branching occurs before the maximum allowed timestep.
                      if (branchingInterval < deltaW) then
                         ! It does, so recheck using the actual branching rate.
                         branchingProbabilityRate=+mergerTreeBranchingProbability_%probability(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,nodeCurrent) &
                              &                   *rootVarianceGrowthFactor_
                         branchingInterval       =+branchingIntervalScaleFree                                                                                                &
                              &                   /branchingProbabilityRate
                         doBranch                =(branchingInterval <= deltaW)
                         if (doBranch) then
                            ! Branching occurred, adjust the accretion fraction, and timestep to their values at the branching event.
                            accretionFraction    =accretionFraction*branchingInterval/deltaW
                            deltaW               =branchingInterval
                            snapAccretionFraction=.false.
                            snapEarliestTime     =.false.
                            controlLimited       =.false.
                            ! Draw a random deviate and scale by the branching rate - this will be used to choose the branch mass.
                            uniformRandom       =self_%workers(numberWorker)%randomNumberGenerator_%uniformSample()       
                            branchingProbability=uniformRandom*branchingProbabilityRate
                         end if
                      else
                         doBranch=.false.
                      end if
                   else
                      doBranch=.false.
                   end if
                else
                   ! In this case we're using the original Cole et al. (2000) algorithm.
                   if (branchingProbability > 0.0d0) then
                      uniformRandom=self_%workers(numberWorker)%randomNumberGenerator_%uniformSample()
                      doBranch=(uniformRandom <= branchingProbability)
                      if (doBranch) then
                         branchingProbability   =+mergerTreeBranchingProbability_%probabilityBound(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,mergerTreeBranchingBoundLower,nodeCurrent) &
                              &                  *deltaW                                                                                                                                                       &
                              &                  *rootVarianceGrowthFactor_
                         if (uniformRandom <= branchingProbability) then
                            doBranch=.true.
                         else
                            branchingProbability=+mergerTreeBranchingProbability_%probability     (branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution                              ,nodeCurrent) &
                                 &               *deltaW                                                                                                                                                       &
                                 &               *rootVarianceGrowthFactor_
                            doBranch=(uniformRandom <= branchingProbability)
                         end if
                         ! First convert the realized probability back to a rate.
                         if (doBranch) branchingProbability=uniformRandom/deltaW
                      end if
                   else
                      doBranch=.false.
                   end if
                end if
                ! Determine the critical overdensity for collapse for the new halo(s).
                if (snapEarliestTime) then
                   deltaCritical_=+self_%workers(numberWorker)%criticalOverdensity_     %value       (time=self_%timeEarliest,mass=branchMassCurrent,node=nodeCurrent) &
                        &         *self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance(time=self_%timeNow     ,mass=branchMassCurrent                 ) &
                        &         /self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance(time=self_%timeEarliest,mass=branchMassCurrent                 )
                else
                   deltaCritical_=+branchDeltaCriticalCurrent                                                                                                          &
                        &         +deltaW
                end if
                if (snapAccretionFraction) then
                   accretionFractionCumulative=self_%accretionLimit
                else
                   accretionFractionCumulative=accretionFractionCumulative+accretionFraction
                end if
                ! Create new nodes.
                select case (doBranch)
                case (.true.)
                   ! Branching occurs - create two progenitors.
                   !$omp atomic
                   self_%nodeIndex =  self_%nodeIndex+1
                   nodeNew1        => treeNode(self_%nodeIndex,nodeCurrent%hostTree)
                   basicNew1       => nodeNew1%basic(autoCreate=.true.)
                   ! Compute mass of one of the new nodes.
                   nodeMass1       =  mergerTreeBranchingProbability_%massBranch(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,branchingProbability/rootVarianceGrowthFactor_,self_%workers(numberWorker)%randomNumberGenerator_,nodeCurrent)
                   nodeMass2       =  basic_%mass()-nodeMass1
                   nodeMass1       =  nodeMass1*(1.0d0-accretionFractionCumulative)
                   nodeMass2       =  nodeMass2*(1.0d0-accretionFractionCumulative)
                   ! Compute the critical overdensity of the first new node.
                   deltaCritical1 =  self_%criticalOverdensityUpdate(deltaCritical_,branchMassCurrent,nodeMass1,nodeNew1)
                   ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                   if (snapEarliestTime.and.Values_Agree(deltaCritical1,deltaCritical_,relTol=toleranceDeltaCritical)) deltaCritical1=deltaCritical_
                   call basicNew1%massSet(nodeMass1     )
                   call basicNew1%timeSet(deltaCritical1)
                   ! Create second progenitor if it would be above the mass resolution.
                   if (nodeMass2 > massResolution) then
                      !$omp atomic
                      self_%nodeIndex =  self_%nodeIndex+1
                      nodeNew2        => treeNode(self_%nodeIndex,nodeCurrent%hostTree)
                      basicNew2       => nodeNew2%basic(autoCreate=.true.)
                      ! Compute the critical overdensity of the second new node.
                      deltaCritical2=self_%criticalOverdensityUpdate(deltaCritical_,branchMassCurrent,nodeMass2,nodeNew2)
                      ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                      if (snapEarliestTime.and.Values_Agree(deltaCritical2,deltaCritical_,relTol=toleranceDeltaCritical)) deltaCritical2=deltaCritical_
                      call basicNew2%massSet(nodeMass2     )
                      call basicNew2%timeSet(deltaCritical2)
                      ! Inform the build controller of these new nodes.
                      call self_%workers(numberWorker)%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeNew1,nodeNew2,didBranch=.true.)
                      ! Create links from old to new nodes and vice-versa. (Ensure that the first child node is the more massive progenitor.)
                      if (nodeMass2 > nodeMass1) then
                         nodeCurrent%firstChild => nodeNew2
                         nodeNew2   %sibling    => nodeNew1
                         call self_%onBranch(tree,massResolution,nodeNew1)
                      else
                         nodeCurrent%firstChild => nodeNew1
                         nodeNew1   %sibling    => nodeNew2
                         call self_%onBranch(tree,massResolution,nodeNew2)
                      end if
                      nodeNew1      %parent     => nodeCurrent
                      nodeNew2      %parent     => nodeCurrent
                   else
                      ! Second branch would be subresolution - do not create it.
                      ! Inform the build controller of the new node.
                      call self_%workers(numberWorker)%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeNew1         ,didBranch=.true.)
                      ! Create links from old to new nodes and vice-versa.
                      nodeCurrent   %firstChild => nodeNew1
                      nodeNew1      %sibling    => null()
                      nodeNew1      %parent     => nodeCurrent
                   end if
                   branchIsDone                 =  .true.
                case (.false.)
                   ! No branching occurs - create one progenitor.
                   if (accretionFractionCumulative >= self_%accretionLimit .or. (controlLimited .and. controlInsertNode)) then
                      !$omp atomic
                      self_%nodeIndex =  self_%nodeIndex+1
                      nodeNew1        => treeNode(self_%nodeIndex,nodeCurrent%hostTree)
                      basicNew1       => nodeNew1%basic(autoCreate=.true.)
                      ! Compute new mass accounting for sub-resolution accretion.
                      nodeMass1       =  basic_%mass()*(1.0d0-accretionFractionCumulative)
                      ! Compute the critical overdensity corresponding to this new node.
                      deltaCritical1  =  self_%criticalOverdensityUpdate(deltaCritical_,branchMassCurrent,nodeMass1,nodeNew1)
                      call basicNew1%massSet(nodeMass1     )
                      call basicNew1%timeSet(deltaCritical1)
                      ! Inform the build controller of this new node.
                      call self_%workers(numberWorker)%mergerTreeBuildController_%nodesInserted(nodeCurrent,nodeNew1)
                      ! Create links from old to new node and vice-versa.
                      nodeCurrent %firstChild => nodeNew1
                      nodeNew1    %parent     => nodeCurrent
                      branchIsDone            =  .true.
                   else
                      ! Insufficient accretion has occurred to warrant making a new node. Simply update the mass and critical
                      ! overdensity and take another step. We update the critical overdensity by mapping to a time at the current
                      ! branch mass, then mapping back to a critical overdensity at the new branch mass. This ensures that if
                      ! critical overdensity is a function of mass we preserve correct time-ordering along the branch.
                      nodeMass1                 =basic_%mass()*(1.0d0-accretionFractionCumulative)
                      branchDeltaCriticalCurrent=self_%criticalOverdensityUpdate(deltaCritical_,branchMassCurrent,nodeMass1,nodeCurrent)
                      branchMassCurrent         =nodeMass1
                   end if
                end select
                ! If the timestep was limited by the build controller, allow the build controller to respond.
                if (controlLimited) branchIsDone=branchIsDone .or. .not.self_%workers(numberWorker)%mergerTreeBuildController_%controlTimeMaximum(nodeCurrent,branchMassCurrent,branchDeltaCriticalCurrent,self_%nodeIndex)
             end if
          end if
       end do
       if (allocated(treeWalker_)) then
          nodeCurrent => null()
       else
          nodeCurrent => nodeCurrent%firstChild
       end if
    end do
    return
  end subroutine cole2000BuildBranch

  recursive subroutine cole2000OnBranch(tree,massResolution,node)
    !!{
    Act on branching.
    !!}
    implicit none
    type            (mergerTree), intent(in   )          :: tree
    double precision            , intent(in   )          :: massResolution
    type            (treeNode  ), intent(inout), pointer :: node
    !$GLC attributes unused :: tree, massResolution, node
    
    return
  end subroutine cole2000OnBranch

  subroutine cole2000TimeEarliestSet(self,timeEarliest)
    !!{
    Set the earliest time for the tree builder.
    !!}
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout) :: self
    double precision                           , intent(in   ) :: timeEarliest

    self%timeEarliest=timeEarliest
    return
  end subroutine cole2000TimeEarliestSet

  double precision function cole2000CriticalOverdensityUpdate(deltaCritical,massCurrent,massNew,nodeNew) result(criticalOverdensity)
    !!{
    Update the critical overdensity for a new node, given that of the parent,
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    double precision          , intent(in   ) :: massCurrent  , massNew, &
         &                                       deltaCritical
    type            (treeNode), intent(inout) :: nodeNew
    double precision                          :: time
    
    if (self_%timeParameterIsMassDependent) then
       ! Critical overdensity is mass-dependent, so convert current critical overdensity to a time at the mass of the parent, and
       ! then convert that time back to a critical overdensity at the mass of the new node.
       time               =+self_%workers(numberWorker)%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=     deltaCritical,mass=massCurrent,node=nodeNew)
       criticalOverdensity=+self_%workers(numberWorker)%criticalOverdensity_     %value         (time               =     time         ,mass=massNew    ,node=nodeNew) &
            &              *self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance  (time               =self_%timeNow     ,mass=massNew                 ) &
            &              /self_%workers(numberWorker)%cosmologicalMassVariance_%rootVariance  (time               =     time         ,mass=massNew                 )
    else
       ! Critical overdensity is mass independent, so the critical overdensity for parent and child node is the same.
       criticalOverdensity=                                                                                               deltaCritical
    end if
    return
  end function cole2000CriticalOverdensityUpdate

  recursive subroutine cole2000ConvertTimeNode(nodeTip)
    !!{
    Convert from critical overdensity to time in a single node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    type            (treeNode          ), intent(inout) :: nodeTip
    class           (nodeComponentBasic), pointer       :: basic_
    double precision                                    :: timeOfCollapse_
    
    ! Get the basic component of the node.
    basic_          => nodeTip                                           %basic         (                                                                 )
    ! Compute the collapse time.
    timeOfCollapse_ =  self_  %workers(numberWorker)%criticalOverdensity_%timeOfCollapse(criticalOverdensity=basic_%time(),mass=basic_%mass(),node=nodeTip)
    call basic_%timeSet(timeOfCollapse_)
    return
  end subroutine cole2000ConvertTimeNode
      
  subroutine cole2000ConvertCriticalOverdensityNode(self,nodeTip)
    !!{
    Convert from time to critical overdensity in a single node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: nodeTip
    class           (nodeComponentBasic       ), pointer       :: basic_
    double precision                                           :: time_   , criticalOverdensity_    , &
         &                                                        timeNow_, rootVarianceGrowthFactor
    
    ! Get the basic component of the node.
    basic_                   =>  nodeTip                          %basic                   (                                             )
    ! Compute the critical overdensity.
    time_                    =   basic_                           %time                    (                                             )
    timeNow_                 =  +self                             %timeNow
    rootVarianceGrowthFactor =  +self   %cosmologicalMassVariance_%rootVariance            (time=time_   ,mass=basic_%mass()             ) &
         &                      /self   %cosmologicalMassVariance_%rootVariance            (time=timeNow_,mass=basic_%mass()             )
    criticalOverdensity_     =  +self   %criticalOverdensity_     %value                   (time=time_   ,mass=basic_%mass(),node=nodeTip) &
         &                      /                                  rootVarianceGrowthFactor
    call basic_%timeSet(criticalOverdensity_)
    return
  end subroutine cole2000ConvertCriticalOverdensityNode

  recursive subroutine cole2000CheckOrderNode(tree,massResolution,node)
    !!{
    Check well-ordering of time for the given node.
    !!}
    use :: Galacticus_Nodes              , only : nodeComponentBasic
    use :: Error                         , only : Warn
    use :: Display                       , only : displayReset                   , displayMagenta
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Unlink_Parent
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    double precision                    , intent(in   )          :: massResolution
    type            (mergerTree        ), intent(in   )          :: tree
    class           (nodeComponentBasic)               , pointer :: basic_                   , basicParent_
    logical                                                      :: closeToResolution
    logical                             , save                   :: warned           =.false.
    
    if (associated(node%parent)) then
       basic_       => node       %basic()
       basicParent_ => node%parent%basic()
       if (basicParent_%time() < basic_%time()) then
          closeToResolution= basicParent_%mass() < massResolution*(1.0d0+self_%toleranceResolutionParent) &
               &            .and.                                                                         &
               &             basic_      %mass() < massResolution*(1.0d0+self_%toleranceResolutionSelf  )
          if (self_%ignoreWellOrdering .or. closeToResolution) then
             ! Parent halo is very close to the resolution limit, or we are ignoring well-ordering errors. Simply prune away the remainder of this branch.
             if (.not.closeToResolution.and..not.warned) then
                !$omp critical(mergerTreeBuilderCole2000WellOrderingWarn)
                if (.not.warned) then
                   call Warn(displayMagenta()//'WARNING:'//displayReset()//' tree is not well-ordered - pruning branch and ignoring')
                   warned=.true.
                end if
                !$omp end critical(mergerTreeBuilderCole2000WellOrderingWarn)
             end if
             call Merger_Tree_Prune_Unlink_Parent(node,node%parent,parentWillBePruned=.false.,preservePrimaryProgenitor=.true.)
             call node%destroyBranch()
             deallocate(node)
             node => null()
          else
             ! Parent halo is not close to the resolution limit - this is an error.
             block
               character(len=20        ) :: label
               type     (varying_string) :: message

               message="branch is not well-ordered in time:"           //char(10)
               write (label,'(i20)'   ) tree       %index
               message=message//" ->      tree index = "//label        //char(10)
               write (label,'(i20)'   ) node       %index()
               message=message//" ->      node index = "//label        //char(10)
               write (label,'(i20)'   ) node%parent%index()
               message=message//" ->    parent index = "//label        //char(10)
               write (label,'(e20.14)')                                                             basic_      %time()
               message=message//" ->       node time = "//label//" Gyr"//char(10)
               write (label,'(e20.14)')                                                             basicParent_%time()
               message=message//" ->     parent time = "//label//" Gyr"//char(10)
               write (label,'(e20.14)')                                                             basic_      %mass()
               message=message//" ->       node mass = "//label//" M☉" //char(10)
               write (label,'(e20.14)')                                                             basicParent_%mass()
               message=message//" ->     parent mass = "//label//" M☉" //char(10)
               write (label,'(e20.14)') self_%workers(numberWorker)%criticalOverdensity_%value(time=basic_      %time(),mass=basic_      %mass(),node=node       )
               message=message//" ->         node δc = "//label        //char(10)
               write (label,'(e20.14)') self_%workers(numberWorker)%criticalOverdensity_%value(time=basicParent_%time(),mass=basicParent_%mass(),node=node%parent)
               message=message//" ->       parent δc = "//label        //char(10)
               basic_ => tree%nodeBase%basic()
               write (label,'(e20.14)')                                   basic_%time()
               message=message//" ->       tree time = "//label//" Gyr"//char(10)
               write (label,'(e20.14)')                                                              basic_     %mass()
               message=message//" ->       tree mass = "//label//" M☉" //char(10)
               write (label,'(e20.14)') massResolution
               message=message//" -> mass resolution = "//label//" M☉"
               call Error_Report(message//{introspection:location})
             end block
          end if
       end if
    end if
    return
  end subroutine cole2000CheckOrderNode

  subroutine cole2000StateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Custom state store function.
    !!}
    implicit none
    class  (mergerTreeBuilderCole2000), intent(inout) :: self
    integer                           , intent(in   ) :: stateFile
    type   (c_ptr                    ), intent(in   ) :: gslStateFile
    integer(c_size_t                 ), intent(in   ) :: stateOperationID
    integer                                           :: i
    
    call self%stateStore_(stateFile,gslStateFile,stateOperationID)
    if (self%workersInitialized) then
       write (stateFile) lbound(self%workers,dim=1),ubound(self%workers,dim=1)
       do i=lbound(self%workers,dim=1),ubound(self%workers,dim=1)
          call        self%workers(i)%cosmologyFunctions_          %stateStore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%criticalOverdensity_         %stateStore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%cosmologicalMassVariance_    %stateStore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%mergerTreeBuildController_   %stateStore(stateFile,gslStateFile,stateOperationID)
          if (self%branchIntervalStep) &
               & call self%workers(i)%branchingIntervalDistribution%stateStore(stateFile,gslStateFile,stateOperationID)
       end do
    end if
    return
  end subroutine cole2000StateStore

  subroutine cole2000StateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Custom state restore function.
    !!}
    use :: Events_Hooks, only : eventsHooksFutureThread
    implicit none
    class  (mergerTreeBuilderCole2000), intent(inout) :: self
    integer                           , intent(in   ) :: stateFile
    type   (c_ptr                    ), intent(in   ) :: gslStateFile
    integer(c_size_t                 ), intent(in   ) :: stateOperationID
    integer                                           :: boundLower      , boundUpper, &
         &                                               i

    call self%stateRestore_(stateFile,gslStateFile,stateOperationID)
    if (self%workersInitialized) then
       read (stateFile) boundLower,boundUpper
       allocate(self%workers(boundLower:boundUpper))
       do i=boundLower,boundUpper
          call eventsHooksFutureThread(i)
          allocate(self%workers(i)%cosmologyFunctions_       ,mold=self%cosmologyFunctions_       )
          allocate(self%workers(i)%criticalOverdensity_      ,mold=self%criticalOverdensity_      )
          allocate(self%workers(i)%cosmologicalMassVariance_ ,mold=self%cosmologicalMassVariance_ )
          allocate(self%workers(i)%mergerTreeBuildController_,mold=self%mergerTreeBuildController_)
          if (self%branchIntervalStep) then
             allocate(self%workers(i)%branchingIntervalDistribution)
             !![
	     <referenceConstruct object="self%workers(i)%branchingIntervalDistribution" constructor="distributionFunction1DNegativeExponential(1.0d0)"/>
             !!]
          end if
          !$omp critical(mergerTreeBuilderCole2000ParallelDeepCopy)
          !![
	  <deepCopyReset variables="self%cosmologyFunctions_ self%criticalOverdensity_ self%cosmologicalMassVariance_ self%mergerTreeBuildController_"/>
	  <deepCopy source="self%cosmologyFunctions_"        destination="self%workers(i)%cosmologyFunctions_"       />
	  <deepCopy source="self%criticalOverdensity_"       destination="self%workers(i)%criticalOverdensity_"      />
	  <deepCopy source="self%cosmologicalMassVariance_"  destination="self%workers(i)%cosmologicalMassVariance_" />
	  <deepCopy source="self%mergerTreeBuildController_" destination="self%workers(i)%mergerTreeBuildController_"/>
	  <deepCopyFinalize variables="self%workers(i)%cosmologyFunctions_ self%workers(i)%criticalOverdensity_ self%workers(i)%cosmologicalMassVariance_ self%workers(i)%mergerTreeBuildController_"/>
          !!]
          !$omp end critical(mergerTreeBuilderCole2000ParallelDeepCopy)
          call eventsHooksFutureThread()
          call        self%workers(i)%cosmologyFunctions_          %stateRestore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%criticalOverdensity_         %stateRestore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%cosmologicalMassVariance_    %stateRestore(stateFile,gslStateFile,stateOperationID)
          call        self%workers(i)%mergerTreeBuildController_   %stateRestore(stateFile,gslStateFile,stateOperationID)
          if (self%branchIntervalStep) &
               & call self%workers(i)%branchingIntervalDistribution%stateRestore(stateFile,gslStateFile,stateOperationID)
       end do
    end if
    return
  end subroutine cole2000StateRestore
  
