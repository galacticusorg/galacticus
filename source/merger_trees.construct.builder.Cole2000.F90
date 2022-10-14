!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  
  !![
  <mergerTreeBuilder name="mergerTreeBuilderCole2000">
   <description>
    A merger tree builder class which uses the algorithm described by \cite{cole_hierarchical_2000} (with minor modifications
    described below). This action of this algorithm is controlled by the following parameters:
    \begin{description}
     \item [{\normalfont \ttfamily [mergeProbability]}] The maximum probability for a binary merger allowed in a single
     timestep. This allows the probability to be kept small, such the the probability for multiple mergers within a single
     timestep is small.
     \item [{\normalfont \ttfamily [accretionLimit]}] The maximum fractional change in mass due to sub-esolution accretion
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
    sub-resolution accretion is removed entirely from the mass $M_3$ and not from $M_2$ this can lead to an assymetry in
    progenitor mass functions close to $M_1/2$. Therefore, we instead set the progenitor masses by first drawing a mass
    $M_2^\prime$ from the mass branching distribution function and then setting
    \begin{eqnarray}
      M_2 &amp;=&amp; M_2^\prime (1-F), \nonumber \\
      M_3 &amp;=&amp; (M_1 - M_2^\prime) (1-F),
    \end{eqnarray}
    which ensures a symmetric treatment of subresolution accretion close to $M_1/2$.
    
   </description>
  </mergerTreeBuilder>
  !!]

  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderCole2000
     !!{
     A merger tree builder class using the algorithm of \cite{cole_hierarchical_2000}.
     !!}
     private
     class           (cosmologyFunctionsClass                  ), pointer :: cosmologyFunctions_                      => null()
     class           (mergerTreeMassResolutionClass            ), pointer :: mergerTreeMassResolution_                => null()
     class           (criticalOverdensityClass                 ), pointer :: criticalOverdensity_                     => null()
     class           (mergerTreeBranchingProbabilityClass      ), pointer :: mergerTreeBranchingProbability_          => null()
     class           (cosmologicalMassVarianceClass            ), pointer :: cosmologicalMassVariance_                => null()
     class           (mergerTreeBuildControllerClass           ), pointer :: mergerTreeBuildController_               => null()
     logical                                                              :: timeParameterIsMassDependent
     ! Variables controlling merger tree accuracy.
     double precision                                                     :: accretionLimit                                    , timeEarliest             , &
          &                                                                  mergeProbability                                  , timeNow                  , &
          &                                                                  redshiftMaximum
     ! Random number sequence variables
     logical                                                              :: branchIntervalStep
     ! Interval distribution.
     logical                                                              :: branchingIntervalDistributionInitialized
     type            (distributionFunction1DNegativeExponential)          :: branchingIntervalDistribution
     ! Tolerances for behavior close to the resolution limit.
     double precision                                                     :: toleranceResolutionSelf                           , toleranceResolutionParent
   contains
     !![
     <methods>
       <method description="Set the critical overdensity object." method="criticalOverdensitySet"   />
       <method description="Set the critical overdensity object." method="criticalOverdensityUpdate"/>
     </methods>
     !!]
     final     ::                              cole2000Destructor
     procedure :: build                     => cole2000Build
     procedure :: timeEarliestSet           => cole2000TimeEarliestSet
     procedure :: criticalOverdensitySet    => cole2000CriticalOverdensitySet
     procedure :: criticalOverdensityUpdate => cole2000CriticalOverdensityUpdate
  end type mergerTreeBuilderCole2000

  interface mergerTreeBuilderCole2000
     !!{
     Constructors for the {\normalfont \ttfamily cole2000} merger tree builder class.
     !!}
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerTreeBuilderCole2000

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
         &                                                                  toleranceResolutionParent
    logical                                                              :: branchIntervalStep

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
         &                                                                                                           branchIntervalStep               , &
         &                                                                                                           toleranceResolutionSelf          , &
         &                                                                                                           toleranceResolutionParent        , &
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

  function cole2000ConstructorInternal(mergeProbability,accretionLimit,timeEarliest,branchIntervalStep,toleranceResolutionSelf,toleranceResolutionParent,mergerTreeBranchingProbability_,mergerTreeMassResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBuildController_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeBuilderCole2000          )                        :: self
    double precision                                     , intent(in   )         :: mergeProbability               , accretionLimit         , &
         &                                                                          timeEarliest                   , toleranceResolutionSelf, &
         &                                                                          toleranceResolutionParent
    logical                                              , intent(in   )         :: branchIntervalStep
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    class           (mergerTreeMassResolutionClass      ), intent(in   ), target :: mergerTreeMassResolution_
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (mergerTreeBuildControllerClass     ), intent(in   ), target :: mergerTreeBuildController_
    !![
    <constructorAssign variables="mergeProbability, accretionLimit, timeEarliest, branchIntervalStep, toleranceResolutionSelf, toleranceResolutionParent, *mergerTreeBranchingProbability_, *mergerTreeMassResolution_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBuildController_"/>
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
    Destructor for the {\normalfont \ttfamily cole2000} merger tree builder class.
    !!}
    implicit none
    type(mergerTreeBuilderCole2000), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%mergerTreeMassResolution_"      />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%mergerTreeBuildController_"     />
    !!]
    return
  end subroutine cole2000Destructor

  subroutine cole2000Build(self,tree)
    !!{
    Build a merger tree.
    !!}
    use :: Display                 , only : displayReset                 , displayMagenta
    use :: Error                   , only : Error_Report                 , Warn
    use :: Galacticus_Nodes        , only : mergerTree                   , nodeComponentBasic              , treeNode
    use :: ISO_Varying_String      , only : varying_string
    use :: Kind_Numbers            , only : kind_int8
    use :: Merger_Tree_Branching   , only : mergerTreeBranchingBoundLower, mergerTreeBranchingBoundUpper
    use :: Merger_Tree_Walkers     , only : mergerTreeWalkerIsolatedNodes, mergerTreeWalkerTreeConstruction
    use :: Numerical_Comparison    , only : Values_Agree
    implicit none
    class           (mergerTreeBuilderCole2000          ), intent(inout)         :: self
    type            (mergerTree                         ), intent(inout), target :: tree
    type            (treeNode                           ), pointer               :: nodeNew1                              , nodeNew2                   , node                      , &
         &                                                                          nodeChild
    class           (nodeComponentBasic                 ), pointer               :: basicNew1                             , basicNew2                  , basic                     , &
         &                                                                          basicParent                           , basicChild
    class           (mergerTreeBranchingProbabilityClass), pointer               :: mergerTreeBranchingProbability_
    double precision                                     , parameter             :: toleranceTimeEarliest          =2.0d-6
    double precision                                     , parameter             :: toleranceDeltaCritical         =1.0d-6
    double precision                                     , parameter             :: toleranceTime                  =1.0d-6
    type            (mergerTreeWalkerTreeConstruction   )                        :: treeWalkerConstruction
    type            (mergerTreeWalkerIsolatedNodes      )                        :: treeWalkerIsolated
    integer         (kind=kind_int8                     )                        :: nodeIndex
    double precision                                                             :: accretionFraction                     , timeNodeBase               , branchingProbability      , &
         &                                                                          collapseTime                          , deltaCritical              , deltaCritical1            , &
         &                                                                          deltaCritical2                        , deltaW                     , nodeMass1                 , &
         &                                                                          nodeMass2                             , deltaCriticalEarliest      , uniformRandom             , &
         &                                                                          massResolution                        , accretionFractionCumulative, branchMassCurrent         , &
         &                                                                          branchDeltaCriticalCurrent            , branchingInterval          , branchingIntervalScaleFree, &
         &                                                                          branchingProbabilityRate              , deltaWAccretionLimit       , deltaWEarliestTime        , &
         &                                                                          collapseTimeTruncate                  , rootVarianceGrowthFactor   , time                      , &
         &                                                                          deltaWController
    logical                                                                      :: doBranch                              , branchIsDone               , snapAccretionFraction     , &
         &                                                                          snapEarliestTime                      , controlLimited
    type            (varying_string                     )                        :: message
    character       (len=20                             )                        :: label

    ! Begin construction.
    nodeIndex =  1               ! Initialize the node index counter to unity.
    node      => tree%nodeBase   ! Point to the base node.
    basic     => node%basic   () ! Get the basic component of the node.
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
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this class).
    timeNodeBase            =                                                  basic%time        ()
    rootVarianceGrowthFactor=+self%cosmologicalMassVariance_%rootVariance(time=      timeNodeBase  ,mass=basic%mass()          ) &
         &                   /self%cosmologicalMassVariance_%rootVariance(time=self %timeNow       ,mass=basic%mass()          )
    deltaCritical           =+self%criticalOverdensity_     %value       (time=basic%time        (),mass=basic%mass(),node=node) &
         &                   /rootVarianceGrowthFactor
    call basic%timeSet(deltaCritical)
    ! Begin tree build loop.
    treeWalkerConstruction=mergerTreeWalkerTreeConstruction(tree)
    do while (treeWalkerConstruction%next(node))
       ! Apply control.
       if (.not.self%mergerTreeBuildController_%control(node,treeWalkerConstruction)) exit
       ! Get the basic component of the node.
       basic                       => node %basic()
       ! Initialize the state for this branch.
       accretionFractionCumulative =  0.0d0
       branchMassCurrent           =  basic%mass ()
       branchDeltaCriticalCurrent  =  basic%time ()
       ! Evolve the branch until mass falls below the resolution limit, the earliest time is reached, or the branch ends.
       branchIsDone                =  .false.
       do while (.not.branchIsDone)
          ! Get the branching probability object to use for this step.
          mergerTreeBranchingProbability_ => self%mergerTreeBuildController_%branchingProbabilityObject(node)
          ! Get the growth factor in the root variance at the mass of the current branch.
          time                    =+self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=     branchDeltaCriticalCurrent,mass=branchMassCurrent,node=node)
          rootVarianceGrowthFactor=+self%cosmologicalMassVariance_%rootVariance  (time               =     time                      ,mass=branchMassCurrent          ) &
               &                   /self%cosmologicalMassVariance_%rootVariance  (time               =self%timeNow                   ,mass=branchMassCurrent          )
          if     (                                                                      &
               &   branchMassCurrent <= massResolution                                  &
               &  .or.                                                                  &
               &   time              <  self%timeEarliest*(1.0d0+toleranceTimeEarliest) &
               & ) then
             ! Branch should be terminated. If we have any accumulated accretion, terminate the branch with a final node.
             if (accretionFractionCumulative > 0.0d0) then
                nodeIndex      =  nodeIndex+1
                nodeNew1       => treeNode(nodeIndex,tree)
                basicNew1      => nodeNew1%basic(autoCreate=.true.)
                ! Compute new mass accounting for sub-resolution accretion.
                nodeMass1      =  basic%mass()*(1.0d0-accretionFractionCumulative)
                call basicNew1%massSet(nodeMass1     )
                ! Compute the critical overdensity corresponding to this new node.
                deltaCritical1=self%criticalOverdensityUpdate(branchDeltaCriticalCurrent,branchMassCurrent,nodeMass1,nodeNew1)
                call basicNew1%timeSet(deltaCritical1)
                ! Create links from old to new node and vice-versa.
                node%firstChild => nodeNew1
                nodeNew1%parent => node
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it).
                node            => nodeNew1
             end if
             ! Flag that the branch is done.
             branchIsDone=.true.
          else
             ! Find branching probability rate per unit deltaW.
             branchingProbabilityRate=+mergerTreeBranchingProbability_%probabilityBound     (branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,mergerTreeBranchingBoundUpper,node) &
                  &                   *rootVarianceGrowthFactor
             ! Find accretion rate.
             accretionFraction       =+mergerTreeBranchingProbability_%fractionSubresolution(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution                              ,node) &
                  &                   *rootVarianceGrowthFactor
             ! A negative accretion fraction indicates that the node is so close to the resolution limit that
             ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
             ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
             if (accretionFraction < 0.0d0) then
                ! Terminate the branch with a final node.
                nodeIndex          =  nodeIndex+1
                nodeNew1           => treeNode      (nodeIndex        ,tree)
                basicNew1          => nodeNew1%basic(autoCreate=.true.     )
                ! Create a node at the mass resolution.
                nodeMass1          =  massResolution
                ! Compute critical overdensity for this new node.
                deltaCritical1      =self%criticalOverdensityUpdate(branchDeltaCriticalCurrent,branchMassCurrent,nodeMass1,nodeNew1)
                ! Ensure critical overdensity exceeds that of the current node.
                collapseTime        =self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=branchDeltaCriticalCurrent,mass=branchMassCurrent,node=nodeNew1)
                collapseTimeTruncate=self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical1            ,mass=nodeMass1        ,node=nodeNew1)
                if (collapseTimeTruncate > collapseTime*(1.0d0+toleranceTime)) then
                   call Error_Report('truncating to resolution, but resolution node exists after parent'//{introspection:location})
                else
                   do while (collapseTimeTruncate > collapseTime*(1.0d0-toleranceTime))
                      deltaCritical1      =deltaCritical1*(1.0d0+toleranceTime)
                      collapseTimeTruncate=self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical1,mass=nodeMass1,node=nodeNew1)
                   end do
                end if
                ! Set properties of the new node.
                call basicNew1%massSet(nodeMass1     )
                call basicNew1%timeSet(deltaCritical1)
                ! Create links from old to new node and vice-versa.
                node    %firstChild => nodeNew1
                nodeNew1%parent     => node
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it), and flag that the branch is done.
                node                => nodeNew1
                branchIsDone        =  .true.
             else
                ! Finding maximum allowed step in w. Limit based on branching rate only if we are using the original Cole et
                ! al. (2000) algorithm.
                deltaW               =mergerTreeBranchingProbability_%stepMaximum(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution)/rootVarianceGrowthFactor
                snapAccretionFraction=.false.
                if (accretionFraction > 0.0d0) then
                   deltaWAccretionLimit=(self%accretionLimit-accretionFractionCumulative)/accretionFraction
                   if (deltaWAccretionLimit <= deltaW) then
                      deltaW               =deltaWAccretionLimit
                      snapAccretionFraction=.true.
                   end if
                end if
                if     (                                                                   &
                     &   branchingProbabilityRate > 0.0d0                                  &
                     &  .and.                                                              &
                     &   .not.self%branchIntervalStep                                      &
                     & ) then
                   if (self%mergeProbability/branchingProbabilityRate < deltaW) then
                      ! Timestep is limited by branching rate. Reduce the timestep to the allowed
                      ! size and unset the flag to snap the accretion fraction to its maximum
                      ! allowed value (since we won't reach that value with this new, reduced
                      ! timestep).
                      deltaW               =self%mergeProbability/branchingProbabilityRate
                      snapAccretionFraction=.false.
                   end if
                end if
                ! Limit the timestep so that the maximum allowed time is not exceeded.
                deltaWEarliestTime=+self%criticalOverdensity_     %value       (time=self%timeEarliest,mass=branchMassCurrent         ,node=node) &
                     &             *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow     ,mass=branchMassCurrent                   ) &
                     &             /self%cosmologicalMassVariance_%rootVariance(time=self%timeEarliest,mass=branchMassCurrent                   ) &
                     &             -                                                                        branchDeltaCriticalCurrent
                if (deltaWEarliestTime < deltaW) then
                   deltaW               =deltaWEarliestTime
                   snapEarliestTime     =.true.
                   snapAccretionFraction=.false.
                else
                   snapEarliestTime     =.false.
                end if
                ! Limit the timestep according to the build controller.
                deltaWController=+self%mergerTreeBuildController_%timeMaximum(node,branchMassCurrent,branchDeltaCriticalCurrent) &
                     &           -branchDeltaCriticalCurrent
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
                if (.not.self%branchIntervalStep)                           &
                     & branchingProbability=branchingProbabilityRate*deltaW
                accretionFraction          =accretionFraction       *deltaW
                ! Accretion fraction must be less than unity. Reduce timestep (and branching
                ! probability and accretion fraction) by factors of two until this condition is
                ! satisfied.
                do while (accretionFraction+accretionFractionCumulative >= 1.0d0)
                   if (.not.self%branchIntervalStep)                      &
                        & branchingProbability=branchingProbability*0.5d0
                   accretionFraction          =accretionFraction   *0.5d0
                   deltaW                     =deltaW              *0.5d0
                   snapAccretionFraction      =.false.
                   snapEarliestTime           =.false.
                   controlLimited             =.false.
                end do
                ! Decide if a branching occurs.
                if (self%branchIntervalStep) then
                   ! In this case we draw intervals between branching events from a negative
                   ! exponential distribution.
                   if (branchingProbabilityRate > 1.0d-100) then
                      branchingIntervalScaleFree=0.0d0
                      do while (branchingIntervalScaleFree <= 0.0d0)
                         branchingIntervalScaleFree=self%branchingIntervalDistribution%sample(randomNumberGenerator_=tree%randomNumberGenerator_)
                      end do
                      branchingInterval=branchingIntervalScaleFree/branchingProbabilityRate
                      ! Based on the upper bound on the rate, check if branching occurs before the maximum allowed timestep.
                      if (branchingInterval < deltaW) then
                         ! It does, so recheck using the actual branching rate.
                         branchingProbabilityRate=+mergerTreeBranchingProbability_%probability(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,node) &
                              &                   *rootVarianceGrowthFactor
                         branchingInterval       =+branchingIntervalScaleFree                                                                                         &
                              &                   /branchingProbabilityRate
                         doBranch                =(branchingInterval <= deltaW)
                         if (doBranch) then
                            ! Branching occured, adjust the accretion fraction, and timestep to their values at the branching event.
                            accretionFraction    =accretionFraction*branchingInterval/deltaW
                            deltaW               =branchingInterval
                            snapAccretionFraction=.false.
                            snapEarliestTime     =.false.
                            controlLimited       =.false.
                            ! Draw a random deviate and scale by the branching rate - this will be used to choose the branch mass.
                            uniformRandom       =tree%randomNumberGenerator_%uniformSample()       
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
                      uniformRandom=tree%randomNumberGenerator_%uniformSample()
                      doBranch=(uniformRandom <= branchingProbability)
                      if (doBranch) then
                         branchingProbability   =+mergerTreeBranchingProbability_%probabilityBound(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,mergerTreeBranchingBoundLower,node) &
                              &                  *deltaW                                                                                                                                                &
                              &                  *rootVarianceGrowthFactor
                         if (uniformRandom <= branchingProbability) then
                            doBranch=.true.
                         else
                            branchingProbability=+mergerTreeBranchingProbability_%probability     (branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution                              ,node) &
                                 &               *deltaW                                                                                                                                                &
                                 &               *rootVarianceGrowthFactor
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
                   deltaCritical=+self%criticalOverdensity_     %value       (time=self%timeEarliest,mass=branchMassCurrent,node=node) &
                        &        *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow     ,mass=branchMassCurrent          ) &
                        &        /self%cosmologicalMassVariance_%rootVariance(time=self%timeEarliest,mass=branchMassCurrent          )
               else
                   deltaCritical=+branchDeltaCriticalCurrent &
                        &        +deltaW
                end if
                if (snapAccretionFraction) then
                   accretionFractionCumulative=self%accretionLimit
                else
                   accretionFractionCumulative=accretionFractionCumulative+accretionFraction
                end if
                ! Create new nodes.
                select case (doBranch)
                case (.true.)
                   ! Branching occurs - create two progenitors.
                   nodeIndex     =  nodeIndex+1
                   nodeNew1      => treeNode(nodeIndex,tree)
                   basicNew1     => nodeNew1%basic(autoCreate=.true.)
                   ! Compute mass of one of the new nodes.
                   nodeMass1     =  mergerTreeBranchingProbability_%massBranch(branchMassCurrent,branchDeltaCriticalCurrent,time,massResolution,branchingProbability/rootVarianceGrowthFactor,tree%randomNumberGenerator_,node)
                   nodeMass2     =  basic%mass()-nodeMass1
                   nodeMass1=nodeMass1*(1.0d0-accretionFractionCumulative)
                   nodeMass2=nodeMass2*(1.0d0-accretionFractionCumulative)
                   ! Compute the critical overdensity of the first new node.
                   deltaCritical1=self%criticalOverdensityUpdate(deltaCritical,branchMassCurrent,nodeMass1,nodeNew1)
                   ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                   if (snapEarliestTime.and.Values_Agree(deltaCritical1,deltaCritical,relTol=toleranceDeltaCritical)) deltaCritical1=deltaCritical
                   call basicNew1%massSet(nodeMass1     )
                   call basicNew1%timeSet(deltaCritical1)
                   ! Create second progenitor.
                   nodeIndex=nodeIndex+1
                   nodeNew2  => treeNode(nodeIndex,tree)
                   basicNew2 => nodeNew2%basic(autoCreate=.true.)
                   ! Compute the critical overdensity of the second new node.
                   deltaCritical2=self%criticalOverdensityUpdate(deltaCritical,branchMassCurrent,nodeMass2,nodeNew2)
                   ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                   if (snapEarliestTime.and.Values_Agree(deltaCritical2,deltaCritical,relTol=toleranceDeltaCritical)) deltaCritical2=deltaCritical
                   call basicNew2%massSet(nodeMass2     )
                   call basicNew2%timeSet(deltaCritical2)
                   !! NOTE: Here we can call our new "nodesInserted" function in the controller with the current and new nodes.
                   ! Inform the build controller of these new nodes.
                   call self%mergerTreeBuildController_%nodesInserted(node,nodeNew1,nodeNew2)
                   ! Create links from old to new nodes and vice-versa. (Ensure that the first child node is the more massive progenitor.)
                   if (nodeMass2 > nodeMass1) then
                      node    %firstChild => nodeNew2
                      nodeNew2%sibling    => nodeNew1
                   else
                      node    %firstChild => nodeNew1
                      nodeNew1%sibling    => nodeNew2
                   end if
                   nodeNew1%parent        => node
                   nodeNew2%parent        => node
                   branchIsDone           =  .true.
                case (.false.)
                   ! No branching occurs - create one progenitor.
                   if (accretionFractionCumulative >= self%accretionLimit) then
                      nodeIndex      =  nodeIndex+1
                      nodeNew1       => treeNode(nodeIndex,tree)
                      basicNew1      => nodeNew1%basic(autoCreate=.true.)
                      ! Compute new mass accounting for sub-resolution accretion.
                      nodeMass1      =  basic%mass()*(1.0d0-accretionFractionCumulative)
                      ! Compute the critical overdensity corresponding to this new node.
                      deltaCritical1=self%criticalOverdensityUpdate(deltaCritical,branchMassCurrent,nodeMass1,nodeNew1)
                      call basicNew1%massSet(nodeMass1     )
                      call basicNew1%timeSet(deltaCritical1)
                      ! Inform the build controller of this new node.
                      call self%mergerTreeBuildController_%nodesInserted(node,nodeNew1)
                      ! Create links from old to new node and vice-versa.
                      node    %firstChild => nodeNew1
                      nodeNew1%parent     => node
                      branchIsDone        =  .true.
                   else
                      ! Insufficient accretion has occured to warrant making a new node. Simply update the mass and critical
                      ! overdensity and take another step. We update the critical overdensity by mapping to a time at the current
                      ! branch mass, then mapping back to a critical overdensity at the new branch mass. This ensures that if
                      ! critical overdensity is a function of mass we preserve correct time-ordering along the branch.
                      nodeMass1                 =basic%mass()*(1.0d0-accretionFractionCumulative)
                      branchDeltaCriticalCurrent=self%criticalOverdensityUpdate(deltaCritical,branchMassCurrent,nodeMass1,node)
                      branchMassCurrent         =nodeMass1
                   end if
                end select
                ! If the timestep was limited by the build controller, allow the build controller to response.
                if (controlLimited) branchIsDone=.not.self%mergerTreeBuildController_%controlTimeMaximum(node,branchMassCurrent,branchDeltaCriticalCurrent,nodeIndex)
             end if
          end if
       end do
    end do
    ! Walk the tree and convert w to time.
    treeWalkerIsolated=mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalkerIsolated%next(node))
       ! Get the basic component of the node.
       basic        => node%basic()
       ! Compute the collapse time.
       collapseTime =  self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=basic%time(),mass=basic%mass(),node=node)
       call basic%timeSet(collapseTime)
    end do
    basic => tree%nodeBase%basic()
    call basic%timeSet(timeNodeBase)
    ! Check for mis-ordering of the base node and its child node(s). This can happen because we force the time of the base node to
    ! be precisely the base time, but for other nodes the time is computed by inverting the w(t)=delta_crit(t)/D(t)
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
       if (associated(node%parent)) then
          basic       => node       %basic()
          basicParent => node%parent%basic()
          if (basicParent%time() < basic%time()) then
             if     (                                                                            &
                  &   basicParent%mass() < massResolution*(1.0d0+self%toleranceResolutionParent) &
                  &  .and.                                                                       &
                  &   basic      %mass() < massResolution*(1.0d0+self%toleranceResolutionSelf  ) &
                  & ) then
                ! Parent halo is very close to the resolution limit. Simply prune away the remainder of this branch.
                call node%destroyBranch()
                deallocate(node)
                call treeWalkerIsolated%previous(node)
             else
                ! Parent halo is not close to the resolution limit - this is an error.
                message="branch is not well-ordered in time:"           //char(10)
                write (label,'(i20)'   ) tree       %index
                message=message//" ->      tree index = "//label        //char(10)
                write (label,'(i20)'   ) node       %index()
                message=message//" ->      node index = "//label        //char(10)
                write (label,'(i20)'   ) node%parent%index()
                message=message//" ->    parent index = "//label        //char(10)
                write (label,'(e20.14)')                                      basic      %time()
                message=message//" ->       node time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                      basicParent%time()
                message=message//" ->     parent time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                                              basic      %mass()
                message=message//" ->       node mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)')                                                              basicParent%mass()
                message=message//" ->     parent mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)') self%criticalOverdensity_%value(time=basic      %time(),mass=basic      %mass(),node=node       )
                message=message//" ->         node δc = "//label        //char(10)
                write (label,'(e20.14)') self%criticalOverdensity_%value(time=basicParent%time(),mass=basicParent%mass(),node=node%parent)
                message=message//" ->       parent δc = "//label        //char(10)
                basic => tree%nodeBase%basic()
                write (label,'(e20.14)')                                   basic%time()
                message=message//" ->       tree time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                                               basic      %mass()
                message=message//" ->       tree mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)') massResolution
                message=message//" -> mass resolution = "//label//" M☉"
                call Error_Report(message//{introspection:location})
             end if
          end if
       end if
    end do
    return
  end subroutine cole2000Build

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

  subroutine cole2000CriticalOverdensitySet(self,criticalOverdensity_)
    !!{
    Set the critical overdensity object for this tree builder.
    !!}
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout)         :: self
    class(criticalOverdensityClass ), intent(in   ), target :: criticalOverdensity_

    !![
    <objectDestructor name="self%criticalOverdensity_"/>
    !!]
    self%criticalOverdensity_         =>      criticalOverdensity_
    !![
    <referenceCountIncrement owner="self" object="criticalOverdensity_"/>
    !!]
    self%timeParameterIsMassDependent = self%criticalOverdensity_     %isMassDependent      () &
         &                             .or.                                                    &
         &                              self%cosmologicalMassVariance_%growthIsMassDependent()
    return
  end subroutine cole2000CriticalOverdensitySet

  double precision function cole2000CriticalOverdensityUpdate(self,deltaCritical,massCurrent,massNew,nodeNew)
    !!{
    Update the critical overdensity for a new node, given that of the parent,
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout) :: self
    double precision                           , intent(in   ) :: massCurrent  , massNew, &
         &                                                        deltaCritical
    type            (treeNode                 ), intent(inout) :: nodeNew
    double precision                                           :: time

    if (self%timeParameterIsMassDependent) then
       ! Critical overdensity is mass-dependent, so convert current critical overdensity to a time at the mass of the parent, and
       ! then convert that time back to a critical overdensity at the mass of the new node.
       time                             =self%criticalOverdensity_      %timeOfCollapse(criticalOverdensity     =deltaCritical,mass=massCurrent,node=nodeNew)
       cole2000CriticalOverdensityUpdate=self%criticalOverdensity_      %value         (time               =     time         ,mass=massNew    ,node=nodeNew) &
            &                            *self%cosmologicalMassVariance_%rootVariance  (time               =self%timeNow      ,mass=massNew                 ) &
            &                            /self%cosmologicalMassVariance_%rootVariance  (time               =     time         ,mass=massNew                 )
    else
       ! Critical overdensity is mass independent, so the critical overdensity for parent and child node is the same.
       cole2000CriticalOverdensityUpdate=                                                             deltaCritical
    end if
    return
  end function cole2000CriticalOverdensityUpdate
