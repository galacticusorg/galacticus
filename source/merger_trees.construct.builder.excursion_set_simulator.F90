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
  Implements a merger tree builder class which builds merger trees by simulating trajectories from the excursion set.
  !!}
  
  use :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  use :: Cosmological_Density_Field        , only : cosmologicalMassVarianceClass, criticalOverdensityClass
  
  !![
  <mergerTreeBuilder name="mergerTreeBuilderExcursionSetSimulator">
    <description>
      \textbf{Warning:} \emph{This is a proof-of-concept implementation---it should not be used to generate
      scientifically-reliable results.} A merger tree builder class which creates trees by simulating trajectories from the
      excursion set. As warned above, this is a proof-of-concept implementation of a merger tree builder using direct simulation
      of the excursion set. It has not been validated or calibrated to produce merger trees that accurately match the statistical
      properties of trees measured from N-body simulations---in fact, it is known to \emph{not} match the statistical properties
      of those trees. It has also not been carefully verified to be free from errors, nor has it been optimized for speed or
      memory footprint. \emph{Caveat arborist!}
    </description>
  </mergerTreeBuilder>
  !!]
  type, extends(mergerTreeBuilderClass) :: mergerTreeBuilderExcursionSetSimulator
     !!{
     A merger tree builder class which creates trees by simulating trajectories from the excursion set.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (mergerTreeMassResolutionClass), pointer :: mergerTreeMassResolution_ => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     double precision                                         :: timeEarliest                       , timeNow                 , &
          &                                                      varianceStepMinimum                , varianceStepSigmaMaximum, &
          &                                                      factorMassConsolidate              , factorTimeConsolidate   , &
          &                                                      excursionStep                      , redshiftMaximum
   contains
     final     ::          excursionSetSimulatorDestructor
     procedure :: build => excursionSetSimulatorBuild
  end type mergerTreeBuilderExcursionSetSimulator

  interface mergerTreeBuilderExcursionSetSimulator
     !!{
     Constructors for the \refClass{mergerTreeBuilderExcursionSetSimulator} merger tree constructor class.
     !!}
     module procedure excursionSetSimulatorConstructorParameters
       module procedure excursionSetSimulatorConstructorInternal
      end interface mergerTreeBuilderExcursionSetSimulator

    contains

      function excursionSetSimulatorConstructorParameters(parameters) result(self)
        !!{
        Constructor for the \refClass{mergerTreeBuilderExcursionSetSimulator} merger tree builder class which takes a parameter set as input.
        !!}
        use :: Input_Parameters, only : inputParameter, inputParameters
        implicit none
        type            (mergerTreeBuilderExcursionSetSimulator)                :: self
        type            (inputParameters                       ), intent(inout) :: parameters
        class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
        class           (mergerTreeMassResolutionClass         ), pointer       :: mergerTreeMassResolution_
        class           (cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassVariance_
        class           (criticalOverdensityClass              ), pointer       :: criticalOverdensity_
        double precision                                                        :: redshiftMaximum          , timeEarliest            , &
             &                                                                     varianceStepMinimum      , varianceStepSigmaMaximum, &
             &                                                                     factorMassConsolidate    , factorTimeConsolidate   , &
             &                                                                     excursionStep
        
        !![
        <inputParameter>
	  <name>redshiftMaximum</name>
	  <source>parameters</source>
	  <defaultValue>10.0d0</defaultValue>
	  <description>The highest redshift to which merger trees will be built.</description>
	</inputParameter>
        <inputParameter>
	  <name>varianceStepMinimum</name>
	  <source>parameters</source>
	  <defaultValue>1.0d-6</defaultValue>
	  <description>
	    The smallest step in variance to take when simulating trajectories from the excursion set. Using a smaller step will
	    result in less quantization of node masses.
	  </description>
	</inputParameter>
        <inputParameter>
	  <name>varianceStepSigmaMaximum</name>
	  <source>parameters</source>
	  <defaultValue>5.0d+0</defaultValue>
	  <description>
	    Controls the size of the steps in variance when simulating trajectories from the excursion set. Specifically, the step
	    is never smaller than $\Delta S = [ ( \delta_\mathrm{i} - \delta ) / N ]^2$, where $N=${\normalfont \ttfamily
	    [varianceStepSigmaMaximum]} is the number of standard deviations (in a standard normal distribution) that would be
	    required in a positive fluctuation to make the excursion exceed the initial excursion in this step. In this way, the
	    probability of ``missing'' a first crossing at some earlier step should be kept small (of order the fraction of a
	    standard normal distribution which lies more than $N$ standard deviations above zero).
	  </description>
	</inputParameter>
        <inputParameter>
	  <name>excursionStep</name>
	  <source>parameters</source>
	  <defaultValue>0.02d0</defaultValue>
	  <description>Sets the size of the timesteps in excursion, $\delta$, to take when constructing the tree.</description>
	</inputParameter>
        <inputParameter>
	  <name>factorMassConsolidate</name>
	  <source>parameters</source>
	  <defaultValue>0.9d0</defaultValue>
	  <description>The maximum factor ($&lt;1$) by which the mass of a halo can have changed before consolidation is no longer permitted.</description>
	</inputParameter>
        <inputParameter>
	  <name>factorTimeConsolidate</name>
	  <source>parameters</source>
	  <defaultValue>0.9d0</defaultValue>
	  <description>The maximum factor ($&lt;1$) by which the time of a halo can have changed before consolidation is no longer permitted.</description>
	</inputParameter>
 	<objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
	<objectBuilder class="mergerTreeMassResolution" name="mergerTreeMassResolution_" source="parameters"/>
	<objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
	<objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
        !!]
        timeEarliest=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
        self        =mergerTreeBuilderExcursionSetSimulator(timeEarliest,varianceStepMinimum,varianceStepSigmaMaximum,excursionStep,factorMassConsolidate,factorTimeConsolidate,cosmologyFunctions_,mergerTreeMassResolution_,criticalOverdensity_,cosmologicalMassVariance_)
        !![
        <inputParametersValidate source="parameters"/>
	<objectDestructor name="cosmologyFunctions_"      />
	<objectDestructor name="mergerTreeMassResolution_"/>
	<objectDestructor name="criticalOverdensity_"     />
	<objectDestructor name="cosmologicalMassVariance_"/>
        !!]
        return
      end function excursionSetSimulatorConstructorParameters

      function excursionSetSimulatorConstructorInternal(timeEarliest,varianceStepMinimum,varianceStepSigmaMaximum,excursionStep,factorMassConsolidate,factorTimeConsolidate,cosmologyFunctions_,mergerTreeMassResolution_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
        !!{
        Internal constructor for the \refClass{mergerTreeBuilderExcursionSetSimulator} merger tree builder class.
        !!}
        implicit none
        type            (mergerTreeBuilderExcursionSetSimulator)                        :: self
        class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
        class           (mergerTreeMassResolutionClass         ), intent(in   ), target :: mergerTreeMassResolution_
        class           (criticalOverdensityClass              ), intent(in   ), target :: criticalOverdensity_
        class           (cosmologicalMassVarianceClass         ), intent(in   ), target :: cosmologicalMassVariance_
        double precision                                        , intent(in   )         :: timeEarliest             , varianceStepMinimum  , &
             &                                                                             varianceStepSigmaMaximum , factorMassConsolidate, &
             &                                                                             factorTimeConsolidate    , excursionStep
        !![
        <constructorAssign variables="timeEarliest, varianceStepMinimum, varianceStepSigmaMaximum, excursionStep, factorMassConsolidate, factorTimeConsolidate, *cosmologyFunctions_, *mergerTreeMassResolution_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
        !!]
        
        self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeEarliest))
        self%timeNow        =self%cosmologyFunctions_%cosmicTime                 (                                                1.0d0 )
        if (factorMassConsolidate <= 0.0d0 .or. factorMassConsolidate > 1.0d0) call Error_Report('[factorMassConsolidate] ∈ [0,1)'//{introspection:location})
        if (factorTimeConsolidate <= 0.0d0 .or. factorMassConsolidate > 1.0d0) call Error_Report('[factorTimeConsolidate] ∈ [0,1)'//{introspection:location})
        return
      end function excursionSetSimulatorConstructorInternal

      subroutine excursionSetSimulatorDestructor(self)
        !!{
        Destructor for the \refClass{mergerTreeBuilderExcursionSetSimulator} merger tree constructor class.
        !!}
        implicit none
        type(mergerTreeBuilderExcursionSetSimulator), intent(inout) :: self
        
        !![
        <objectDestructor name="self%cosmologyFunctions_"      />
	<objectDestructor name="self%mergerTreeMassResolution_"/>
	<objectDestructor name="self%criticalOverdensity_"     />
	<objectDestructor name="self%cosmologicalMassVariance_"/>
        !!]
        return
      end subroutine excursionSetSimulatorDestructor

      subroutine excursionSetSimulatorBuild(self,tree)
        !!{
        Build a merger tree with a smooth mass accretion history.
        !!}
        use :: Error_Functions    , only : Error_Function
        use :: Galacticus_Nodes   , only : nodeComponentBasic              , treeNode
        use :: Kind_Numbers       , only : kind_int8
        use :: Merger_Tree_Walkers, only : mergerTreeWalkerTreeConstruction
        implicit none
        class           (mergerTreeBuilderExcursionSetSimulator), intent(inout), target  :: self
        type            (mergerTree                            ), intent(inout), target  :: tree
        type            (treeNode                              )               , pointer :: node                   , nodeNew              , &
             &                                                                              nodePrior              , nodeProgenitorPrimary
        class           (nodeComponentBasic                    )               , pointer :: basic                  , basicNew             , &
             &                                                                              basicParent
        double precision                                        , parameter              :: varianceStepTiny=1.0d-3
        integer         (kind=kind_int8                        )                         :: indexNode
        double precision                                                                 :: massNode               , timeNode             , &
             &                                                                              varianceMinimum        , varianceMaximum      , &
             &                                                                              massResolution         , alphaMinimum         , &
             &                                                                              excursion              , excursionInitial     , &
             &                                                                              variance               , rootVarianceMinimum  , &
             &                                                                              massNew                , timeNew              , &
             &                                                                              massPatchRemaining     , massProgenitorMaximum, &
             &                                                                              excursionStep          , excursionMinimum     , &
             &                                                                              varianceStep           , volumePatchRemaining , &
             &                                                                              volumePatchInitial     , excursionTarget      , &
             &                                                                              fractionAccretion
        type            (mergerTreeWalkerTreeConstruction      )                         :: treeWalker
        logical                                                                          :: acceptTrajectory
        
        ! Get the mass resolution for this tree.
        massResolution  =self%mergerTreeMassResolution_%resolution  (     tree                                 )
        ! Find the maximum variance for this tree. This corresponds to the variance at the mass resolution.
        varianceMaximum =self%cosmologicalMassVariance_%rootVariance(time=self%timeNow     ,mass=massResolution)**2
        ! Extract the base node properties.
        node      => tree%nodeBase
        indexNode =  1_kind_int8
        ! Begin tree build loop.
        treeWalker=mergerTreeWalkerTreeConstruction(tree)
        do while (treeWalker%next(node))
           ! Extract properties of the current node. Set the prior progenitor node to null (as we have no progenitors yet).
           basic     => node %basic()
           massNode  =  basic%mass ()
           timeNode  =  basic%time ()
           nodePrior =>       null ()           
           ! If this node is at or below mass resolution, or exists at or before the earliest time then it can have no
           ! progenitors, so we simply move on to the next node in the tree.
           if (massNode <= massResolution .or. timeNode <= self%timeEarliest) cycle           
           ! Find the accretion fraction for this node and timestep.
           fractionAccretion=Error_Function(                                                                                            &
                &                           +self%excursionStep                                                                         &
                &                           /sqrt(                                                                                      &
                &                                 +2.0d0                                                                                &
                &                                 *(                                                                                    &
                &                                   +                               varianceMaximum                                     &
                &                                   -self%cosmologicalMassVariance_%rootVariance   (time=self%timeNow,mass=massNode)**2 &
                &                                  )                                                                                    &
                &                                )                                                                                      &
                &                          )
           ! Find the initial excursion, δᵢ, for this node. This is the critical overdensity at the time at which the node exists,
           ! scaled to the present day following assuming linear growth.
           excursionInitial=+self%criticalOverdensity_     %value       (time=     timeNode,mass=massNode) &
                &           *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow ,mass=massNode) &
                &           /self%cosmologicalMassVariance_%rootVariance(time=     timeNode,mass=massNode)
           ! Find the target excursion for the next step.
           excursionTarget =+     excursionInitial &
                &           +self%excursionStep
           ! Set the initial remaining mass and (Lagrangian, comoving) volume of the patch of mass which we will partition into
           ! progenitors. The initial mass, Mᵢ, is just the mass of the node. The initial volume is Vᵢ = Mᵢ/ρ̅(1+δᵢ), but we
           ! actually just track vᵢ = 1/(1+δᵢ) for convenience, ignoring the constant (for this node) factors.
           volumePatchInitial  =+  1.0d0             &
                &               /(                   &
                &                 +1.0d0             &
                &                 +excursionInitial  &
                &                )
           massPatchRemaining  =+massNode            &
                &               *(                   &
                &                 +1.0d0             &
                &                 -fractionAccretion &
                &                )
           volumePatchRemaining=+volumePatchInitial  &
                &               *(                   &
                &                 +1.0d0             &
                &                 -fractionAccretion &
                &                )
           ! Begin partitioning the mass into progenitors until the remaining mass is less than the resolution limit.
           do while (massPatchRemaining > massResolution)
              ! Find the variance corresponding to the remaining mass in the patch. This is the starting variance from our
              ! trajectory.
              varianceMinimum=self%cosmologicalMassVariance_%rootVariance(time=self%timeNow,mass=massPatchRemaining)**2
              ! Determine the initial excursion, δ, for the remaining mass in the patch. Since we know the mass and volume of the
              ! remaining part of the patch the overdensity must be:
              !
              !  δ = m/v-1,
              !
              ! where m is the mass of the patch expressed as a fraction of the original mass. Note that this is at least somewhat
              ! similar to the approach of Sheth (1995; https://ui.adsabs.harvard.edu/abs/1995MNRAS.276..796S) for Poisson
              ! clustering models.
              excursionMinimum=+(                    &
                   &             +massPatchRemaining &
                   &             /massNode           &
                   &            )                    &
                   &           /volumePatchRemaining &
                   &           -1.0d0
              ! Begin simulating trajectories until we have one that is accepted.
              acceptTrajectory=.false.
              do while (.not.acceptTrajectory)
                 ! Construct the excursion. Set the initial point for the excursion.
                 variance =varianceMinimum
                 excursion=excursionMinimum
                 ! Grow the excursion in variance until a first crossing occurs. Note that the following loop condition is
                 ! unbounded - the loop will be exited once a first-crossing occurs.
                 do while (.true.)
                    ! Determine the step to take in variance. As we are interested only in first crossing events we can take a
                    ! large step if the current excursion is far below the initial excursion. Specifically, in such cases we take
                    ! a step in variance equal to ΔS = [ ( δᵢ - δ ) / N ]², where N is the number of standard deviations (in a
                    ! standard normal distribution) that would be required in a positive fluctuation to make the excursion exceed
                    ! the target excursion in this step. In this way, the probability of "missing" a first crossing at some
                    ! earlier step should be kept small (of order the fraction of a standard normal distribution which lies more
                    ! than N standard deviations above zero).
                    varianceStep=max(                                 &
                         &             self%varianceStepMinimum     , &
                         &           (                                &
                         &            +(                              &
                         &              +    excursionTarget          &
                         &              -    excursion                &
                         &             )                              &
                         &            /self%varianceStepSigmaMaximum  &
                         &           )**2                             &
                         &          )
                    ! Draw the step in excursion from a zero-mean normal distribution with variance equal to the step variance.
                    excursionStep=+tree%randomNumberGenerator_%standardNormalSample() &
                         &        *sqrt(varianceStep)
                    ! Accumulate the step.
                    variance =+ variance+ varianceStep
                    excursion=+excursion+excursionStep
                    ! Look for new first crossings (or non-crossings, which correspond to the trajectory exceeding the maximum
                    ! variance - i.e. no first-crossing above the mass resolution).
                    if (excursion >= excursionTarget) then
                       ! A new maximum excursion is reached - this corresponds to a branching event. Determine the mass at this first-crossing.
                       if (variance-varianceMinimum < varianceStepTiny) then
                          ! For very small step in variance we use a linear approximation to determine the mass of the new
                          ! node. This avoids rounding errors which could otherwise cause the mass of the new node to slightly
                          ! exceed that of the original node. The change in mass in this case is given by:
                          !
                          !  ΔM = ½ α⁻¹ (M/S) ΔS,
                          !
                          ! where α = d log σ / d log M.
                          call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(                                                          &
                               &                                                                 time                           =self%timeNow            , &
                               &                                                                 mass                           =     massPatchRemaining , &
                               &                                                                 rootVariance                   =     rootVarianceMinimum, &
                               &                                                                 rootVarianceLogarithmicGradient=     alphaMinimum         &
                               &                                                                )
                          massNew=+massPatchRemaining  &
                               &  *(                   &
                               &    +1.0d0             &
                               &    +0.5d0             &
                               &    /alphaMinimum      &
                               &    *(                 &
                               &      +variance        &
                               &      /varianceMinimum &
                               &      -1.0d0           &
                               &     )                 &
                               &   )
                       else
                          ! For larger steps in variance just compute the new mass from the variance directly.
                          massNew=self%cosmologicalMassVariance_%mass                           (                                                          &
                               &                                                                 rootVariance                   =     sqrt(variance)     , &
                               &                                                                 time                           =self%timeNow              &
                               &                                                                )
                       end if
                       ! Decide whether to accept. Draws from the excursion set essentially give us a volume-weighted sampling of
                       ! progenitor halo masses. But, we want a mass-weighted sample. To achieve this we use rejection sampling
                       ! with a rejection probability of Mₘᵢₙ/M₂, where M₂ is the mass of the potential progenitor halo.
                       acceptTrajectory=tree%randomNumberGenerator_%uniformSample() <= massResolution/massNew                       
                       if (acceptTrajectory) then
                          ! The trajectory was accepted. If the variance is below the maximum variance, then a new first crossing
                          ! occurred above the mass resolution. In such cases we add a new node to the tree. Otherwise, we have a
                          ! subresolution progenitor and no new node is added.
                          if (variance < varianceMaximum) then
                             indexNode        =  indexNode+1_kind_int8
                             timeNew          =  self    %criticalOverdensity_%timeOfCollapse(           excursionTarget,massNew)
                             nodeNew          => treeNode                                    (           indexNode      ,tree   )
                             basicNew         => nodeNew %basic                              (autoCreate=.true.                 )
                             nodeNew  %parent => node
                             call basicNew%massSet(massNew)
                             call basicNew%timeSet(timeNew)
                             ! If another progenitor already exists, link our new progenitor in as a sibling. Otherwise link our
                             ! new progenitor as the first child of the parent node.
                             if (associated(nodePrior)) then
                                nodePrior%sibling    => nodeNew
                             else
                                node     %firstChild => nodeNew
                             end if
                             ! Keep a pointer to this new node as the prior progenitor (for use in establishing sibling pointers
                             ! in any subsequent progenitors).
                             nodePrior => nodeNew
                          end if
                          ! Subtract the mass and volume of this new progenitor node from the remaining mass and volume of the
                          ! patch.
                          massPatchRemaining  =+massPatchRemaining   &
                               &               -massNew
                          volumePatchRemaining=+volumePatchRemaining &
                               &               -massNew              &
                               &               /massNode             &
                               &               /(                    &
                               &                 +1.0d0              &
                               &                 +excursionTarget    &
                               &                )
                       end if
                       ! We are done creating this progenitor, so can exit the trajectory simulation loop and move on to the
                       ! next progenitor (if any).
                       exit
                    else if (variance > varianceMaximum) then
                       ! If the trajectory is subresolution, reject it and try again.
                       acceptTrajectory=.false.
                       exit
                    end if
                 end do
              end do
           end do
           ! Ensure that the most massive progenitor is the first child.
           !! First, find the most massive progenitor.
           massProgenitorMaximum =  0.0d0
           nodeNew               => node%firstChild
           nodeProgenitorPrimary => null()
           do while (associated(nodeNew))
              basicNew => nodeNew%basic()
              if (basicNew%mass() > massProgenitorMaximum) then
                 nodeProgenitorPrimary => nodeNew
                 massProgenitorMaximum =  basicNew%mass()
              end if
              nodeNew => nodeNew%sibling
           end do
           !! If we have a most massive progenitor (we may not if all progenitors were sub-resolution), then ensure it is the
           !! first child.
           if (associated(nodeProgenitorPrimary)) then
              if (.not.associated(node%firstChild,nodeProgenitorPrimary)) then
                 nodeNew => node%firstChild
                 do while (.not.associated(nodeNew%sibling,nodeProgenitorPrimary))
                    nodeNew => nodeNew%sibling
                 end do
                 nodeNew              %sibling    => nodeProgenitorPrimary%sibling
                 nodeProgenitorPrimary%sibling    => node                 %firstChild
                 node                 %firstChild => nodeProgenitorPrimary
              end if
              ! If no branching has occurred consecutively, see if we can consolidate along the branch.        
              if     (                                               &
                   &        associated(node%parent                 ) & ! Can consolidate only if the parent node has a parent node.
                   &  .and.                                          &
                   &        associated(node%parent%firstChild ,node) & ! Can consolidate only if the parent node is the primary progenitor.
                   &  .and.                                          &
                   &   .not.associated(node%sibling                ) & ! Can consolidate only if the parent node has no sibling.
                   &  .and.                                          &
                   &   .not.associated(node%firstChild%sibling     ) & ! Can consolidate only if the parent node's primary progenitor has no sibling.
                   & ) then
                 basicNew    => node%firstChild%basic()
                 basicParent => node%parent    %basic()
                 ! Check if the change in mass and time is sufficiently small that we can consolidate.
                 if     (                                                                 &
                      &   basicNew%mass() > self%factorMassConsolidate*basicParent%mass() &
                      &  .and.                                                            &
                      &   basicNew%time() > self%factorTimeConsolidate*basicParent%time() &
                      & ) then
                    ! Consolidate the node by unlinking, and then destroying the intermediate node. Then set the current node to
                    ! be the primary progenitor.
                    nodePrior                       => node     %parent
                    node     %firstChild%parent     => node     %parent
                    node     %parent    %firstChild => node     %firstChild
                    node     %firstChild%sibling    => node     %sibling
                    call node%destroy()
                    deallocate(node)
                    node                            => nodePrior
                    call treeWalker%setNode(nodePrior)
                 end if
              end if
           end if
        end do
        return
      end subroutine excursionSetSimulatorBuild
