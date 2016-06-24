!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of a merger tree builder using the algorithm of \cite{cole_hierarchical_2000}.
  use Merger_Trees_Build_Mass_Resolution

  !# <mergerTreeBuilder name="mergerTreeBuilderCole2000">
  !#  <description>Merger trees are built using the algorithm of \cite{cole_hierarchical_2000}.</description>
  !# </mergerTreeBuilder>
  use Statistics_Distributions

  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderCole2000
     !% A merger tree builder class using the algorithm of \cite{cole_hierarchical_2000}.
     private
     class           (mergerTreeMassResolutionClass), pointer :: mergerTreeMassResolution_
     ! Variables controlling merger tree accuracy.
     double precision                                         :: accretionLimit                          , timeEarliest      , &
          &                                                      mergeProbability
     ! Option controlling random number sequences.
     logical                                                  :: randomSeedsFixed
     ! Random number sequence variables
     logical                                                  :: branchIntervalStep
     ! Interval distribution.
     logical                                                  :: branchingIntervalDistributionInitialized
     type            (distributionNegativeExponential       ) :: branchingIntervalDistribution
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeBuilderCole2000</object>
     !@   <objectMethod>
     !@     <method>shouldAbort</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textless type(mergerTree)\textgreater\ tree\argin</arguments>
     !@     <description>Return true if construction of the merger tree should be aborted.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>shouldFollowBranch</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textless type(mergerTree)\textgreater\ tree\argin, \textless type(treeNode)\textgreater\ node\argin</arguments>
     !@     <description>Return true if the branch should be followed.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>validateParameters</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Validate the parameters of an instance of this class, aborting if any are invalid.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: build              => cole2000Build
     procedure :: shouldAbort        => cole2000ShouldAbort
     procedure :: shouldFollowBranch => cole2000ShouldFollowBranch
     procedure :: timeEarliestSet    => cole2000TimeEarliestSet
     procedure :: validateParameters => cole2000ValidateParameters
  end type mergerTreeBuilderCole2000

  interface mergerTreeBuilderCole2000
     !% Constructors for the {\normalfont \ttfamily cole2000} merger tree builder class.
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerTreeBuilderCole2000
  
contains

  function cole2000ConstructorParameters(parameters)
    !% Constructor for the \cite{cole_hierarchical_2000} merger tree building class which reads parameters from a provided parameter list.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    implicit none
    type            (mergerTreeBuilderCole2000)                :: cole2000ConstructorParameters    
    type            (inputParameters          ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    double precision                                           :: redshiftMaximum
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>mergeProbability</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%mergeProbability</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The largest probability of branching allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>accretionLimit</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%accretionLimit</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The largest fractional mass change due to subresolution accretion allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d5</defaultValue>
    !#   <description>The highest redshift to which merger trees will be built in the \cite{cole_hierarchical_2000} method.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomSeedsFixed</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%randomSeedsFixed</variable>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether the random number sequence should be restarted for each tree using a deterministically derived (from the tree index) seed. This allows the exact same tree to be generated even when running multiple threads.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>branchIntervalStep</name>
    !#   <source>parameters</source>
    !#   <variable>cole2000ConstructorParameters%branchIntervalStep</variable>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If {\normalfont \ttfamily false} use the original \cite{cole_hierarchical_2000} method to determine whether branching occurs in a timestep. If {\normalfont \ttfamily true} draw branching intervals from a negative exponential distribution.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="mergerTreeMassResolution" name="cole2000ConstructorParameters%mergerTreeMassResolution_" source="parameters"/>
    ! Convert maximum redshift to earliest time.
    cosmologyFunctions_ => cosmologyFunctions()
    cole2000ConstructorParameters%timeEarliest=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    ! Initialize state.
    cole2000ConstructorParameters%branchingIntervalDistributionInitialized=.false.
    ! Validate parameters.
    call cole2000ConstructorParameters%validateParameters()
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(mergeProbability,accretionLimit,timeEarliest,randomSeedsFixed,branchIntervalStep,mergerTreeMassResolution_)
    !% Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    implicit none
    type            (mergerTreeBuilderCole2000    )                        :: cole2000ConstructorInternal
    double precision                               , intent(in   )         :: mergeProbability           , accretionLimit    , &
         &                                                                    timeEarliest
    logical                                        , intent(in   )         :: randomSeedsFixed           , branchIntervalStep
    class           (mergerTreeMassResolutionClass), intent(in   ), target :: mergerTreeMassResolution_
    
    ! Store options.
    cole2000ConstructorInternal%mergeProbability          =  mergeProbability
    cole2000ConstructorInternal%accretionLimit            =  accretionLimit
    cole2000ConstructorInternal%timeEarliest              =  timeEarliest
    cole2000ConstructorInternal%randomSeedsFixed          =  randomSeedsFixed
    cole2000ConstructorInternal%branchIntervalStep        =  branchIntervalStep
    cole2000ConstructorInternal%mergerTreeMassResolution_ => mergerTreeMassResolution_
    ! Initialize state.
    cole2000ConstructorInternal%branchingIntervalDistributionInitialized=.false.
    ! Validate parameters.
    call cole2000ConstructorInternal%validateParameters()
    return
  end function cole2000ConstructorInternal

  subroutine cole2000ValidateParameters(self)
    !% Validate parameters for the {\normalfont \ttfamily Cole2000} merger tree builder class.
    use Galacticus_Error
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout) :: self

    if (self%accretionLimit >= 1.0d0) call Galacticus_Error_Report('cole2000ValidateParameters','accretionLimit < 1 required')
    return
  end subroutine cole2000ValidateParameters

  subroutine cole2000Build(self,tree)
    !% Build a merger tree.
    use Galacticus_Nodes
    use Galacticus_Error
    use Critical_Overdensities
    use Merger_Tree_Branching
    use Merger_Tree_Branching_Options
    use Pseudo_Random
    use Kind_Numbers
    use ISO_Varying_String
    use Numerical_Comparison
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout)         :: self
    type            (mergerTree               ), intent(inout), target :: tree
    type            (treeNode                 ), pointer               :: newNode1                  , newNode2                   , thisNode                  , &
         &                                                                previousNode
    class           (criticalOverdensityClass ), pointer               :: criticalOverdensity_
    class           (nodeComponentBasic       ), pointer               :: newBasic1                 , newBasic2                  , thisBasic                 , &
         &                                                                parentBasic
    double precision                           , parameter             :: toleranceResolutionParent=1.0d-3
    double precision                           , parameter             :: toleranceResolutionSelf  =1.0d-6
    integer         (kind=kind_int8           )                        :: nodeIndex
    double precision                                                   :: accretionFraction         , baseNodeTime               , branchingProbability      , &
         &                                                                collapseTime              , deltaCritical              , deltaCritical1            , &
         &                                                                deltaCritical2            , deltaW                     , nodeMass1                 , &
         &                                                                nodeMass2                 , time                       , uniformRandom             , &
         &                                                                massResolution            , accretionFractionCumulative, branchMassCurrent         , &
         &                                                                branchDeltaCriticalCurrent, branchingInterval          , branchingIntervalScaleFree, &
         &                                                                branchingProbabilityRate  , deltaWAccretionLimit       , deltaWEarliestTime
    logical                                                            :: doBranch                  , branchIsDone               , snapAccretionFraction     , &
         &                                                                snapEarliestTime
    type            (varying_string           )                        :: message
    character       (len=20                   )                        :: label
    
    ! Get default objects.
    criticalOverdensity_ => criticalOverdensity()
    ! Begin construction.
    nodeIndex =  1                   ! Initialize the node index counter to unity.
    thisNode  => tree    %baseNode   ! Point to the base node.
    thisBasic => thisNode%basic   () ! Get the basic component of the node.
    if (.not.self%branchingIntervalDistributionInitialized.and.self%branchIntervalStep) then
       ! Note that we use a unit rate - we will scale the results to the actual rate required.
       self%branchingIntervalDistribution           =distributionNegativeExponential(1.0d0)
       self%branchingIntervalDistributionInitialized=.true.
    end if
    ! Restart the random number sequence.
    call tree%randomNumberGenerator%initialize()
    uniformRandom=tree%randomNumberGenerator%sample(ompThreadOffset=.false.,incrementSeed=int(tree%index))
    ! Get the mass resolution for this tree.
    massResolution=self%mergerTreeMassResolution_%resolution(tree)
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this module).
    baseNodeTime =thisBasic%time()
    deltaCritical=criticalOverdensity_%value(time=thisBasic%time(),mass=thisBasic%mass())
    call thisBasic%timeSet(deltaCritical)
    ! Begin tree build loop.    
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasic => thisNode%basic()
       ! Initialize the state for this branch.
       accretionFractionCumulative=0.0d0
       branchMassCurrent          =thisBasic%mass()
       branchDeltaCriticalCurrent =thisBasic%time()
       ! Evolve the branch until mass falls below the resolution limit, the earliest time is reached, or the branch ends.       
       branchIsDone=.false.
       do while (.not.branchIsDone)
          if     (                                                                                                                          &
               &   branchMassCurrent                           <= massResolution                                                            &
               &  .or.                                                                                                                      &
               &   branchDeltaCriticalCurrent                  >= criticalOverdensity_%value(time=self%timeEarliest,mass=branchMassCurrent) &
               &  .or.                                                                                                                      &
               &   .not.self%shouldFollowBranch(tree,thisNode)                                                                              &
               &) then
             ! Branch should be terminated. If we have any accumulated accretion, terminate the branch with a final node.
             if (accretionFractionCumulative > 0.0d0) then
                nodeIndex      =  nodeIndex+1
                newNode1       => treeNode(nodeIndex,tree)
                newBasic1      => newNode1%basic(autoCreate=.true.)
                ! Compute new mass accounting for sub-resolution accretion.
                nodeMass1      =  thisBasic%mass()*(1.0d0-accretionFractionCumulative)
                ! Compute the time corresponding to this new node.
                time           =  criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical,mass=branchMassCurrent)
                ! Set properties of the new node.
                deltaCritical1 =  criticalOverdensity_%value         (time               =time         ,mass=nodeMass1        )
                call newBasic1%massSet(nodeMass1     )
                call newBasic1%timeSet(deltaCritical1)
                ! Create links from old to new node and vice-versa.
                thisNode%firstChild => newNode1
                newNode1%parent     => thisNode
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it).
                thisNode            => newNode1
             end if
             ! Flag that the branch is done.
             branchIsDone=.true.
          else
             ! Find branching probability rate per unit deltaW.
             branchingProbabilityRate=Tree_Branching_Probability_Bound(branchMassCurrent,branchDeltaCriticalCurrent,massResolution,boundUpper)
             ! Find accretion rate.
             accretionFraction       =Tree_Subresolution_Fraction     (branchMassCurrent,branchDeltaCriticalCurrent,massResolution           )
             ! A negative accretion fraction indicates that the node is so close to the resolution limit that
             ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
             ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
             if (accretionFraction < 0.0d0) then
                ! Terminate the branch with a final node.
                nodeIndex          =  nodeIndex+1
                newNode1           => treeNode      (nodeIndex,tree)
                newBasic1          => newNode1%basic(autoCreate=.true. )
                ! Compute new mass accounting for sub-resolution accretion.
                nodeMass1          = massResolution
                ! Compute the time corresponding to this event.
                time               = criticalOverdensity_%timeOfCollapse(criticalOverdensity=branchDeltaCriticalCurrent,mass=branchMassCurrent)
                ! Set properties of the new node.
                deltaCritical1     = criticalOverdensity_%value         (time               =time                      ,mass=nodeMass1        )
                call newBasic1%massSet(nodeMass1     )
                call newBasic1%timeSet(deltaCritical1)
                ! Create links from old to new node and vice-versa.
                thisNode%firstChild => newNode1
                newNode1%parent     => thisNode
                ! Move to the terminating node (necessary otherwise we would move to this terminating node next and continue to
                ! grow a branch from it), and flag that the branch is done.
                thisNode            => newNode1
                branchIsDone        =  .true.
             else
                ! Finding maximum allowed step in w. Limit based on branching rate only if we are using the original Cole et
                ! al. (2000) algorithm.
                deltaW               =Tree_Maximum_Step(branchMassCurrent,branchDeltaCriticalCurrent,massResolution)
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
                deltaWEarliestTime=+criticalOverdensity_%value(                        &
                     &                                         time=self%timeEarliest, &
                     &                                         mass=branchMassCurrent  &
                     &                                        )                        &
                     &             -branchDeltaCriticalCurrent
                if (deltaWEarliestTime < deltaW) then
                   deltaW               =deltaWEarliestTime
                   snapEarliestTime     =.true.
                   snapAccretionFraction=.false.
                else
                   snapEarliestTime     =.false.
                end if
                ! Scale values to the determined timestep.
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
                end do
                ! Decide if a branching occurs.
                if (self%branchIntervalStep) then
                   ! In this case we draw intervals between branching events from a negative
                   ! exponential distribution.
                   if (branchingProbabilityRate > 0.0d0) then
                      branchingIntervalScaleFree=0.0d0
                      do while (branchingIntervalScaleFree <= 0.0d0)
                         branchingIntervalScaleFree=self%branchingIntervalDistribution%sample(randomNumberGenerator=tree%randomNumberGenerator)
                      end do
                      branchingInterval=branchingIntervalScaleFree/branchingProbabilityRate
                      ! Based on the upper bound on the rate, check if branching occurs before the maximum allowed timestep.
                      if (branchingInterval < deltaW) then
                         ! It does, so recheck using the actual branching rate.
                         branchingProbabilityRate=Tree_Branching_Probability(branchMassCurrent,branchDeltaCriticalCurrent,massResolution)
                         branchingInterval       =branchingIntervalScaleFree/branchingProbabilityRate
                         doBranch                =(branchingInterval <= deltaW)
                         if (doBranch) then
                            ! Branching occured, adjust the accretion fraction, and timestep to their values at the branching event.
                            accretionFraction    =accretionFraction*branchingInterval/deltaW
                            deltaW               =branchingInterval
                            snapAccretionFraction=.false.
                            snapEarliestTime     =.false.
                            ! Draw a random deviate and scale by the branching rate - this will be used to choose the branch mass.
                            uniformRandom       =tree%randomNumberGenerator%sample()
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
                      uniformRandom=tree%randomNumberGenerator%sample()
                      doBranch=(uniformRandom <= branchingProbability)
                      if (doBranch) then
                         branchingProbability=Tree_Branching_Probability_Bound(branchMassCurrent,branchDeltaCriticalCurrent,massResolution,boundLower)*deltaW
                         if (uniformRandom <= branchingProbability) then
                            doBranch=.true.
                         else
                            branchingProbability=Tree_Branching_Probability(branchMassCurrent,branchDeltaCriticalCurrent,massResolution)*deltaW
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
                   deltaCritical                 =+criticalOverdensity_%value(                        &
                        &                                                        time=self%timeEarliest, &
                        &                                                        mass=branchMassCurrent  &
                        &                                                       )
                else
                   deltaCritical                 =+branchDeltaCriticalCurrent &
                        &                         +deltaW
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
                   newNode1      => treeNode(nodeIndex,tree)
                   newBasic1     => newNode1%basic(autoCreate=.true.)
                   ! Compute mass of one of the new nodes.
                   nodeMass1     =  Tree_Branch_Mass(branchMassCurrent,branchDeltaCriticalCurrent,massResolution,branchingProbability,tree%randomNumberGenerator)
                   nodeMass2     =  thisBasic%mass()-nodeMass1
                   nodeMass1=nodeMass1*(1.0d0-accretionFractionCumulative)
                   nodeMass2=nodeMass2*(1.0d0-accretionFractionCumulative)
                   ! Compute the time corresponding to this branching event.
                   time          =  criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical,mass=branchMassCurrent)
                   ! Set properties of first new node.
                   deltaCritical1=  criticalOverdensity_%value         (time               =time         ,mass=nodeMass1        )
                   ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                   if (snapEarliestTime.and.Values_Agree(deltaCritical1,deltaCritical,relTol=1.0d-6)) deltaCritical1=deltaCritical
                   call newBasic1%massSet(nodeMass1     )
                   call newBasic1%timeSet(deltaCritical1)
                   ! Create second progenitor.
                   nodeIndex=nodeIndex+1
                   newNode2  => treeNode(nodeIndex,tree)
                   newBasic2 => newNode2%basic(autoCreate=.true.)
                   ! Set properties of second new node.
                   deltaCritical2=criticalOverdensity_%value(time=time,mass=nodeMass2)
                   ! If we are to snap halos to the earliest time, and the computed deltaCritical is sufficiently close to that time, snap it.
                   if (snapEarliestTime.and.Values_Agree(deltaCritical2,deltaCritical,relTol=1.0d-6)) deltaCritical2=deltaCritical
                   call newBasic2%massSet(nodeMass2     )
                   call newBasic2%timeSet(deltaCritical2)
                   ! Create links from old to new nodes and vice-versa. (Ensure that child node is the more massive progenitor.)
                   if (nodeMass2 > nodeMass1) then
                      thisNode%firstChild => newNode2
                      newNode2%sibling    => newNode1
                   else
                      thisNode%firstChild => newNode1
                      newNode1%sibling    => newNode2
                   end if
                   newNode1%parent        => thisNode
                   newNode2%parent        => thisNode
                   branchIsDone           =  .true.
                case (.false.)
                   ! No branching occurs - create one progenitor.
                   if (accretionFractionCumulative >= self%accretionLimit) then
                      nodeIndex      =  nodeIndex+1
                      newNode1       => treeNode(nodeIndex,tree)
                      newBasic1      => newNode1%basic(autoCreate=.true.)
                      ! Compute new mass accounting for sub-resolution accretion.
                      nodeMass1      =  thisBasic%mass()*(1.0d0-accretionFractionCumulative)
                      ! Compute the time corresponding to this new node.
                      time           =  criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaCritical,mass=branchMassCurrent)
                      ! Set properties of the new node.
                      deltaCritical1 =  criticalOverdensity_%value         (time               =time         ,mass=nodeMass1        )
                      call newBasic1%massSet(nodeMass1     )
                      call newBasic1%timeSet(deltaCritical1)
                      ! Create links from old to new node and vice-versa.
                      thisNode%firstChild => newNode1
                      newNode1%parent     => thisNode
                      branchIsDone=.true.
                   else
                      branchMassCurrent         =thisBasic%mass()*(1.0d0-accretionFractionCumulative)
                      branchDeltaCriticalCurrent=deltaCritical
                   end if
                end select
             end if
          end if
       end do
       ! Check if tree should be aborted.
       if (self%shouldAbort(tree)) then
          thisNode => null()
       else
          ! Walk to the next node.
          thisNode => thisNode%walkTreeUnderConstruction()
       end if
    end do
    ! Walk the tree and convert w to time.
    thisNode => tree%baseNode
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasic    => thisNode%basic()
       ! Compute the collapse time.
       collapseTime =  criticalOverdensity_%timeOfCollapse(criticalOverdensity=thisBasic%time(),mass=thisBasic%mass())
       call thisBasic%timeSet(collapseTime)
       thisNode => thisNode%walkTree()
    end do
    thisBasic => tree%baseNode%basic()
    call thisBasic%timeSet(baseNodeTime)
    ! Check for well-ordering in time.
    thisNode     => tree%baseNode
    previousNode => thisNode
    do while (associated(thisNode))       
       if (associated(thisNode%parent)) then
          thisBasic   => thisNode       %basic()
          parentBasic => thisNode%parent%basic()
          if (parentBasic%time() <= thisBasic%time()) then
             if     (                                                                       &
                  &   parentBasic%mass() < massResolution*(1.0d0+toleranceResolutionParent) &
                  &  .and.                                                                  &
                  &     thisBasic%mass() < massResolution*(1.0d0+toleranceResolutionSelf  ) &
                  & ) then
                ! Parent halo is very close to the resolution limit. Simply prune away the remainder of this branch.
                call thisNode%destroyBranch()
                deallocate(thisNode)
                thisNode => previousNode
             else
                ! Parent halo is not close to the resolution limit - this is an error.
                message="branch is not well-ordered in time:"           //char(10)
                write (label,'(i20)'   ) thisNode       %index()
                message=message//" ->      node index = "//label        //char(10)
                write (label,'(i20)'   ) thisNode%parent%index()
                message=message//" ->    parent index = "//label        //char(10)
                write (label,'(e20.14)')                                   thisBasic%time()
                message=message//" ->       node time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                 parentBasic%time()
                message=message//" ->     parent time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                                           thisBasic%mass()
                message=message//" ->       node mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)')                                                         parentBasic%mass()
                message=message//" ->     parent mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)') criticalOverdensity_%value(time=  thisBasic%time(),mass=  thisBasic%mass())
                message=message//" ->         node δc = "//label        //char(10)
                write (label,'(e20.14)') criticalOverdensity_%value(time=parentBasic%time(),mass=parentBasic%mass())
                message=message//" ->       parent δc = "//label        //char(10)
                thisBasic => tree%baseNode%basic()
                write (label,'(e20.14)')                                   thisBasic%time()
                message=message//" ->       tree time = "//label//" Gyr"//char(10)
                write (label,'(e20.14)')                                                           thisBasic%mass()
                message=message//" ->       tree mass = "//label//" M☉" //char(10)
                write (label,'(e20.14)') massResolution
                message=message//" -> mass resolution = "//label//" M☉"
                call Galacticus_Error_Report('cole2000Build',message)
             end if
          end if
       end if
       previousNode => thisNode
       thisNode     => thisNode%walkTree()
    end do
    return
  end subroutine cole2000Build

  logical function cole2000ShouldAbort(self,tree)
    !% Return {\normalfont \ttfamily true} if tree construction should be aborted. In the {\normalfont \ttfamily cole2000} tree
    !% builder we never abort.
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout) :: self
    type (mergerTree               ), intent(in   ) :: tree
    !GCC$ attributes unused :: self, tree
    
    cole2000ShouldAbort=.false.
    return
  end function cole2000ShouldAbort
  
  logical function cole2000ShouldFollowBranch(self,tree,node)
    !% Return {\normalfont \ttfamily true} if tree construction should continue to follow the current branch. In the {\normalfont
    !% \ttfamily cole2000} tree builder we always continue.
    implicit none
    class(mergerTreeBuilderCole2000), intent(inout)          :: self
    type (mergerTree               ), intent(in   )          :: tree
    type (treeNode                 ), intent(inout), pointer :: node
    !GCC$ attributes unused :: self, tree, node

    cole2000ShouldFollowBranch=.true.
    return
  end function cole2000ShouldFollowBranch
  
  subroutine cole2000TimeEarliestSet(self,timeEarliest)
    !% Set the earliest time for the tree builder.
    implicit none
    class           (mergerTreeBuilderCole2000), intent(inout) :: self
    double precision                           , intent(in   ) :: timeEarliest

    self%timeEarliest=timeEarliest
    return
  end subroutine cole2000TimeEarliestSet
