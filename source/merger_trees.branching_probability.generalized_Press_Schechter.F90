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
Implements a merger tree branching probability class using a generalized Press-Schechter approach.
!!}

  use :: Cosmological_Density_Field     , only : cosmologicalMassVarianceClass              , criticalOverdensityClass
  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Excursion_Sets_First_Crossings , only : excursionSetFirstCrossingClass
  use :: Merger_Tree_Branching_Modifiers, only : mergerTreeBranchingProbabilityModifierClass
  use :: Root_Finder                    , only : rootFinder

  !![
  <mergerTreeBranchingProbability name="mergerTreeBranchingProbabilityGnrlzdPrssSchchtr">
   <description>
    A merger tree branching probability class using a generalized Press-Schechter approach. Branching probabilities are
    computed from solutions to the excursion set barrier first crossing rate problem. Specifically, the branching probability
    per unit time is:
    \begin{equation}
     {\d f \over \d t} = {\mathrm{d} t \over \mathrm{d}\omega} \int_{M_\mathrm{min}}^{M/2} {M \over M^\prime} {\d f \over \d t}
     {\d S \over \d M^\prime} \left| {\d t \over \d \omega}\right| G[\omega,\sigma(M),\sigma(M^\prime)] \d M^\prime,
    \end{equation}
    where $\omega = \delta_\mathrm{c,0}/D(t)$. The rate of accretion of mass in halos below the resolution limit of the merger
    tree is
    \begin{equation}
     {\d R \over \d t} = {\mathrm{d} t \over \mathrm{d}\omega} \int_0^{M_\mathrm{min}} {\d f \over \d t} {\d S \over \d
     M^\prime} \left| {\d t \over \d \omega}\right| G[\omega,\sigma(M),\sigma(M^\prime)] \d M^\prime.
    \end{equation}
    In the above, $G[\omega,\sigma(M),\sigma(M^\prime)]$ is a modification to the merger rate as computed by the selected
    \refClass{mergerTreeBranchingProbabilityModifierClass}. If {\normalfont \ttfamily [smoothAccretion]}$=${\normalfont \ttfamily true}
    then smooth accretion (i.e. accretion of matter not in dark matter halos) is accounted for at the rate:
    \begin{equation}
     {\d R_\mathrm{s} \over \d t} = {\mathrm{d} t \over \mathrm{d}\omega} G[\omega,\sigma_\mathrm{max},\sigma(M^\prime)] {\d
     \stackrel{\sim}{f} \over \d t},
    \end{equation}
    where $\sigma_\mathrm{max}$ is the peak value of $\sigma(M)$ (for the lowest mass halos) and $\d \stackrel{\sim}{f}/\d t$
    is the rate at which excursion set trajectories \emph{fail} to cross the barrier on any mass scale.
   </description>
  </mergerTreeBranchingProbability>
  !!]
  type, extends(mergerTreeBranchingProbabilityClass) :: mergerTreeBranchingProbabilityGnrlzdPrssSchchtr
     !!{
     A merger tree branching probability class using a generalized Press-Schechter approach.
     !!}
     private
     class           (cosmologicalMassVarianceClass              ), pointer :: cosmologicalMassVariance_                  => null()
     class           (criticalOverdensityClass                   ), pointer :: criticalOverdensity_                       => null()
     class           (cosmologyFunctionsClass                    ), pointer :: cosmologyFunctions_                        => null()
     class           (excursionSetFirstCrossingClass             ), pointer :: excursionSetFirstCrossing_                 => null()
     class           (mergerTreeBranchingProbabilityModifierClass), pointer :: mergerTreeBranchingProbabilityModifier_    => null()
     type            (rootFinder                                 )          :: finder
     ! Parent halo shared variables.
     double precision                                                       :: parentDTimeDDeltaCritical                           , parentDelta           , &
          &                                                                    parentHaloMass                                      , parentSigma           , &
          &                                                                    parentSigmaSquared                                  , parentTime            , &
          &                                                                    probabilityMinimumMass                              , probabilitySeek
     ! Record of mass resolution.
     double precision                                                       :: resolutionSigma                                     , massResolutionPrevious
     ! Accuracy parameter to ensure that steps in critical overdensity do not become too large.
     double precision                                                       :: deltaStepMaximum
     ! The maximum sigma that we expect to find.
     double precision                                                       :: sigmaMaximum
     ! Record of whether we have tested the excursion set routines.
     logical                                                                :: excursionSetsTested
     ! Control for inclusion of smooth accretion rates.
     logical                                                                :: smoothAccretion
     ! Record of issued warnings.
     logical                                                                :: subresolutionFractionIntegrandFailureWarned
     ! Minimum mass to which subresolution fractions will be integrated.
     double precision                                                       :: massMinimum
     ! Current epoch.
     double precision                                                       :: timeNow
   contains
     !![
     <methods>
       <method description="Compute common factors needed for the calculations." method="computeCommonFactors"/>
       <method description="Compute common factors needed for the calculations." method="excursionSetTest"    />
     </methods>
     !!]
     final     ::                          generalizedPressSchechterDestructor
     procedure :: probability           => generalizedPressSchechterProbability
     procedure :: probabilityBound      => generalizedPressSchechterProbabilityBound
     procedure :: fractionSubresolution => generalizedPressSchechterFractionSubresolution
     procedure :: massBranch            => generalizedPressSchechterMassBranch
     procedure :: stepMaximum           => generalizedPressSchechterStepMaximum
     procedure :: computeCommonFactors  => generalizedPressSchechterComputeCommonFactors
     procedure :: excursionSetTest      => generalizedPressSchechterExcursionSetTest
  end type mergerTreeBranchingProbabilityGnrlzdPrssSchchtr

  interface mergerTreeBranchingProbabilityGnrlzdPrssSchchtr
     !!{
     Constructors for the {\normalfont \ttfamily generalizedPressSchechter} merger tree builder class.
     !!}
     module procedure generalizedPressSchechterConstructorParameters
     module procedure generalizedPressSchechterConstructorInternal
  end interface mergerTreeBranchingProbabilityGnrlzdPrssSchchtr

  ! Module-scope pointer to self used for root-finding.
  class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), pointer   :: generalizedPressSchechterSelf
  !$omp threadprivate(generalizedPressSchechterSelf)

  ! Module-scope variables used in integrands.
  type            (treeNode                                       ), pointer   :: generalizedPressSchechterNode
  !$omp threadprivate(generalizedPressSchechterNode)

  ! Branching probability integrand integration tolerance.
  double precision                                                 , parameter :: generalizedPressSchechterIntegrandToleranceRelative=1.0d-2

contains

  function generalizedPressSchechterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``generalizedPressSchechter'' merger tree branching probability class which reads parameters from a
    provided parameter list.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    class           (criticalOverdensityClass                       ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                  ), pointer       :: cosmologicalMassVariance_
    class           (cosmologyFunctionsClass                        ), pointer       :: cosmologyFunctions_
    class           (excursionSetFirstCrossingClass                 ), pointer       :: excursionSetFirstCrossing_
    class           (mergerTreeBranchingProbabilityModifierClass    ), pointer       :: mergerTreeBranchingProbabilityModifier_
    double precision                                                                 :: deltaStepMaximum                       , massMinimum
    logical                                                                          :: smoothAccretion

    !![
    <inputParameter>
      <name>deltaStepMaximum</name>
      <defaultValue>0.1d0</defaultValue>
      <description>Limits the step in $\delta_\mathrm{crit}$ when constructing merger trees using the generalized Press-Schechter branching algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d6</defaultValue>
      <description>The minimum mass to use in computing subresolution accretion rates when constructing merger trees using the generalized Press-Schechter branching algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>smoothAccretion</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to include smooth accretion in subresolution accretion rates when constructing merger trees using the generalized Press-Schechter branching algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"                    name="criticalOverdensity_"                    source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"               name="cosmologicalMassVariance_"               source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                     name="cosmologyFunctions_"                     source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing"              name="excursionSetFirstCrossing_"              source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbabilityModifier" name="mergerTreeBranchingProbabilityModifier_" source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityGnrlzdPrssSchchtr(deltaStepMaximum,massMinimum,smoothAccretion,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,excursionSetFirstCrossing_,mergerTreeBranchingProbabilityModifier_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"                   />
    <objectDestructor name="cosmologicalMassVariance_"              />
    <objectDestructor name="cosmologyFunctions_"                    />
    <objectDestructor name="excursionSetFirstCrossing_"             />
    <objectDestructor name="mergerTreeBranchingProbabilityModifier_"/>
    !!]
    return
  end function generalizedPressSchechterConstructorParameters

  function generalizedPressSchechterConstructorInternal(deltaStepMaximum,massMinimum,smoothAccretion,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,excursionSetFirstCrossing_,mergerTreeBranchingProbabilityModifier_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr)                        :: self
    class           (cosmologicalMassVarianceClass                  ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                       ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetFirstCrossingClass                 ), intent(in   ), target :: excursionSetFirstCrossing_
    class           (mergerTreeBranchingProbabilityModifierClass    ), intent(in   ), target :: mergerTreeBranchingProbabilityModifier_
    double precision                                                 , intent(in   )         :: deltaStepMaximum                              , massMinimum
    logical                                                          , intent(in   )         :: smoothAccretion
    double precision                                                 , parameter             :: toleranceAbsolute                      =0.0d+0, toleranceRelative=1.0d-9
    !![
    <constructorAssign variables="deltaStepMaximum, massMinimum, smoothAccretion, *criticalOverdensity_, *cosmologicalMassVariance_, *cosmologyFunctions_, *excursionSetFirstCrossing_, *mergerTreeBranchingProbabilityModifier_"/>
    !!]

    self%excursionSetsTested                        =.false.
    self%subresolutionFractionIntegrandFailureWarned=.false.
    self%massResolutionPrevious                     =-1.0d0
    self%timeNow                                    =self%cosmologyFunctions_%cosmicTime(1.0d0)
    self%finder                                     =rootFinder(                                                           &
         &                                                      rootFunction     =generalizedPressSchechterMassBranchRoot, &
         &                                                      toleranceAbsolute=toleranceAbsolute                      , &
         &                                                      toleranceRelative=toleranceRelative                        &
         &                                                     )
    return
  end function generalizedPressSchechterConstructorInternal

  subroutine generalizedPressSchechterDestructor(self)
    implicit none
    type(mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"                   />
    <objectDestructor name="self%cosmologicalMassVariance_"              />
    <objectDestructor name="self%cosmologyFunctions_"                    />
    <objectDestructor name="self%excursionSetFirstCrossing_"             />
    <objectDestructor name="self%mergerTreeBranchingProbabilityModifier_"/>
    !!]
    return
  end subroutine generalizedPressSchechterDestructor

  subroutine generalizedPressSchechterExcursionSetTest(self,node)
    !!{
    Make a call to excursion set routines with the maximum $\sigma$ that we will use to ensure that they can handle it.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout) :: self
    type            (treeNode                                       ), intent(inout) :: node
    double precision                                                                 :: presentTime    , testResult, &
         &                                                                              varianceMaximum

    if (.not.self%excursionSetsTested) then
       if (self%massMinimum > 0.0d0) then
          presentTime      =self%cosmologyFunctions_      %cosmicTime  (1.0d0                                                 )
          self%sigmaMaximum=self%cosmologicalMassVariance_%rootVariance(self%massMinimum                     ,presentTime     )
          varianceMaximum  =self%sigmaMaximum**2
          testResult       =self%excursionSetFirstCrossing_%rate       (0.5d0*varianceMaximum,varianceMaximum,presentTime,node)
       end if
       self%excursionSetsTested=.true.
    end if
    return
  end subroutine generalizedPressSchechterExcursionSetTest

  double precision function generalizedPressSchechterMassBranch(self,haloMass,deltaCritical,time,massResolution,probabilityFraction,randomNumberGenerator_,node)
    !!{
    Determine the mass of one of the halos to which the given halo branches, given the branching probability, {\normalfont
    \ttfamily probabilityFraction}. Typically, {\normalfont \ttfamily probabilityFraction} is found by multiplying {\normalfont \ttfamily
    Generalized\_Press\_Schechter\_Branching\_Probability()} by a random variable drawn in the interval 0--1 if a halo
    branches. This routine then finds the progenitor mass corresponding to this value.
    !!}
    use :: Display           , only : displayMessage, displayVerbosity, verbosityLevelWarn
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical                  , haloMass           , &
         &                                                                                      massResolution                 , probabilityFraction, &
         &                                                                                      time
    class           (randomNumberGeneratorClass                     ), intent(inout)         :: randomNumberGenerator_
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                 , parameter             :: smallProbabilityFraction=1.0d-3
    type            (varying_string                                 )                        :: message
    character       (len=26                                         )                        :: label
    !$GLC attributes unused :: randomNumberGenerator_

    ! Ensure excursion set calculations have sufficient range in sigma.
    call self%excursionSetTest(node)
    ! Initialize global variables.
    generalizedPressSchechterSelf => self
    self%probabilityMinimumMass   =  massResolution
    self%probabilitySeek          =  probabilityFraction
    call self%computeCommonFactors(node,haloMass,deltaCritical,time)
    ! Check that the root is bracketed.
    if     (                                                           &
         &     generalizedPressSchechterMassBranchRoot(massResolution) &
         &    *generalizedPressSchechterMassBranchRoot(0.5d0*haloMass) &
         &  >=                                                         &
         &    0.0d0                                                    &
         & ) then
       ! Warn about this situation.
       if (displayVerbosity() >= verbosityLevelWarn) then
          message="halo branching mass root is not bracketed in generalizedPressSchechterMassBranch()"
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6,a1,e12.6)') massResolution,":",generalizedPressSchechterMassBranchRoot(massResolution)
          message=" => massMinimum:rootFunction(massMinimum) = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6,a1,e12.6)') 0.5d0*haloMass,":",generalizedPressSchechterMassBranchRoot(0.5d0*haloMass)
          message=" => massMaximum:rootFunction(massMaximum) = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6)') probabilityFraction
          message=" =>                           probability = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
       end if
       ! If the root function is positive at half of the parent halo mass then we have a binary split.
       if (generalizedPressSchechterMassBranchRoot(0.5d0*haloMass) >= 0.0d0) then
          ! Check that we are sufficiently close to zero. If we're not, it might indicate a problem.
          if     (                                                                           &
               &   generalizedPressSchechterMassBranchRoot(0.5d0*haloMass)                   &
               &  >                                                                          &
               &   probabilityFraction*smallProbabilityFraction                              &
               & ) call Error_Report(                                                        &
               &                     "numerical accuracy lost in root finding calculation"// &
               &                     {introspection:location}                                &
               &                    )
          ! Return a binary split mass.
          generalizedPressSchechterMassBranch=0.5d0*haloMass
          return
       end if
    end if
    ! Find the branch mass.
    generalizedPressSchechterMassBranch=self%finder%find(rootRange=[massResolution,0.5d0*haloMass])
    return
  end function generalizedPressSchechterMassBranch

  double precision function generalizedPressSchechterMassBranchRoot(massMaximum)
    !!{
    Root function used in solving for the branch mass.
    !!}
    use :: Numerical_Integration, only : GSL_Integ_Gauss15, integrator
    implicit none
    double precision            , intent(in   ) :: massMaximum
    type            (integrator)                :: integrator_

    integrator_                            = integrator           (                                                                        &
         &                                                                           generalizedPressSchechterProbabilityIntegrand       , &
         &                                                         toleranceRelative=generalizedPressSchechterIntegrandToleranceRelative , &
         &                                                         integrationRule  =GSL_Integ_Gauss15                                     &
         &                                                        )
    generalizedPressSchechterMassBranchRoot=+                                        generalizedPressSchechterSelf%probabilitySeek         &
         &                                  -integrator_%integrate(                                                                        &
         &                                                                           generalizedPressSchechterSelf%probabilityMinimumMass, &
         &                                                                                                                    massMaximum  &
         &                                                        )
    return
  end function generalizedPressSchechterMassBranchRoot

  double precision function generalizedPressSchechterStepMaximum(self,haloMass,deltaCritical,time,massResolution)
    !!{
    Return the maximum allowed step in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily
    deltaCritical} should be allowed to take.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout) :: self
    double precision                                                 , intent(in   ) :: deltaCritical , haloMass, &
         &                                                                              massResolution, time
    !$GLC attributes unused :: deltaCritical, haloMass, massResolution, time

    generalizedPressSchechterStepMaximum=self%deltaStepMaximum
    return
  end function generalizedPressSchechterStepMaximum

  double precision function generalizedPressSchechterProbabilityBound(self,haloMass,deltaCritical,time,massResolution,bound,node)
    !!{
    Return bounds onthe probability per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily
    haloMass} at time {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass greater than
    {\normalfont \ttfamily massResolution}.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout)         :: self
    double precision                                                 , intent(in   )         :: deltaCritical , haloMass, &
         &                                                                                      massResolution, time
    type            (enumerationMergerTreeBranchingBoundType        ), intent(in   )         :: bound
    type            (treeNode                                       ), intent(inout), target :: node
    !$GLC attributes unused :: bound

    generalizedPressSchechterProbabilityBound=self%probability(haloMass,deltaCritical,time,massResolution,node)
    return
  end function generalizedPressSchechterProbabilityBound

  double precision function generalizedPressSchechterProbability(self,haloMass,deltaCritical,time,massResolution,node)
    !!{
    Return the probability per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at
    time {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass greater than {\normalfont
    \ttfamily massResolution}.
    !!}
    use :: Numerical_Integration, only : GSL_Integ_Gauss15, integrator
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical , haloMass   , &
         &                                                                                      massResolution, time
    type            (treeNode                                       ), intent(inout), target :: node
    type            (integrator                                     )                        :: integrator_
    double precision                                                                         :: massMaximum   , massMinimum

    call self%excursionSetTest(node)
    ! Get sigma and delta_critical for the parent halo.
    if (haloMass > 2.0d0*massResolution) then
       call self%computeCommonFactors(node,haloMass,deltaCritical,time)
       massMinimum                          =             massResolution
       massMaximum                          =  0.5d0*self%parentHaloMass
       generalizedPressSchechterSelf        => self
       integrator_                          =  integrator           (                                                                       &
            &                                                                          generalizedPressSchechterProbabilityIntegrand      , &
            &                                                        toleranceRelative=generalizedPressSchechterIntegrandToleranceRelative, &
            &                                                        integrationRule  =GSL_Integ_Gauss15                                    &
            &                                                       )
       generalizedPressSchechterProbability =  integrator_%integrate(                                                                       &
            &                                                                          massMinimum                                        , &
            &                                                                          massMaximum                                          &
            &                                                       )
    else
       generalizedPressSchechterProbability=0.0d0
    end if
    return
  end function generalizedPressSchechterProbability

  double precision function generalizedPressSchechterFractionSubresolution(self,haloMass,deltaCritical,time,massResolution,node)
    !!{
    Return the fraction of mass accreted in subresolution halos, i.e. those below {\normalfont \ttfamily massResolution}, per unit change in
 $\delta_\mathrm{crit}$ for a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily deltaCritical}. The integral is computed numerically.
 !!}
    use :: Display              , only : displayMagenta   , displayReset
    use :: Error                , only : Warn             , errorStatusSuccess
    use :: ISO_Varying_String   , only : varying_string
    use :: Numerical_Integration, only : GSL_Integ_Gauss15, integrator
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical                                 , haloMass   , &
         &                                                                                      massResolution                                , time
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                 , parameter             :: resolutionSigmaOverParentSigmaTolerance=1.0d-3
    double precision                                                                         :: massMaximum                                   , massMinimum, &
         &                                                                                      resolutionSigmaOverParentSigma
    type            (integrator                                     )                        :: integrator_
    integer                                                                                  :: errorStatus
    type            (varying_string                                 )                        :: message

    call self%excursionSetTest(node)
    ! Get sigma and delta_critical for the parent halo.
    call self%computeCommonFactors(node,haloMass,deltaCritical,time)
    ! If requested, compute the rate of smooth accretion.
    if (self%smoothAccretion) then
       generalizedPressSchechterFractionSubresolution=+abs(self%parentDTimeDDeltaCritical)                                                                                                            &
            &                                         *    self%excursionSetFirstCrossing_             %rateNonCrossing(              self%parentSigmaSquared                  ,self%parentTime,node) &
            &                                         *    self%mergerTreeBranchingProbabilityModifier_%rateModifier   (node,haloMass,self%parentSigma       ,self%sigmaMaximum,self%parentTime     )
    else
       generalizedPressSchechterFractionSubresolution=0.0d0
    end if

    if (massResolution /= self%massResolutionPrevious) then
       self%resolutionSigma       =self%cosmologicalMassVariance_%rootVariance(massResolution,self%parentTime)
       self%massResolutionPrevious=massResolution
    end if
    resolutionSigmaOverParentSigma=self%resolutionSigma/self%parentSigma
    if (resolutionSigmaOverParentSigma >= 1.0d0) then
       generalizedPressSchechterSelf                 =>  self
       massMinimum                                   =   self%massMinimum
       massMaximum                                   =        massResolution
       integrator_                                   =   integrator           (                                                                           &
            &                                                                                    generalizedPressSchechterFractionSubresolutionIntegrand, &
            &                                                                  toleranceRelative=1.0d-3                                                 , &
            &                                                                  integrationRule  =GSL_Integ_Gauss15                                        &
            &                                                                 )
       generalizedPressSchechterFractionSubresolution=  +generalizedPressSchechterFractionSubresolution                                                   &
            &                                           +integrator_%integrate(                                                                           &
            &                                                                                    massMinimum                                            , &
            &                                                                                    massMaximum                                            , &
            &                                                                  status           =errorStatus                                              &
            &                                                                 )
       if (errorStatus /= errorStatusSuccess) then
          if (resolutionSigmaOverParentSigma < 1.0d0+resolutionSigmaOverParentSigmaTolerance) then
             generalizedPressSchechterFractionSubresolution=-1.0d0
          else
             ! Attempt the integral again with lower tolerance. Issue a warnings if this is the first time this has happened.
             if (.not.self%subresolutionFractionIntegrandFailureWarned) then
                message=displayMagenta()//'WARNING:'                                                                                                          //displayReset(  )// &
                     &                    ' Integration of the subresolution fraction in the generalized Press-Schechter branching probability module failed.'//char        (10)// &
                     &                    'Will try again with lower tolerance. This warning will not be issued again.'                                                         // &
                     &                    {introspection:location}
                call Warn(message)
                self%subresolutionFractionIntegrandFailureWarned=.true.
             end if
             call integrator_%toleranceSet(toleranceRelative=1.0d-2)
             generalizedPressSchechterFractionSubresolution=+generalizedPressSchechterFractionSubresolution &
                  &                                         +integrator_%integrate(massMinimum,massMaximum)
          end if
       end if
    else
       generalizedPressSchechterFractionSubresolution=-1.0d0
    end if
    return
  end function generalizedPressSchechterFractionSubresolution

  double precision function generalizedPressSchechterProbabilityIntegrand(childHaloMass)
    !!{
    Integrand for the branching probability.
    !!}
    implicit none
    double precision, intent(in   ) :: childHaloMass
    double precision                :: childAlpha   , childSigma

    call generalizedPressSchechterSelf%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,generalizedPressSchechterSelf%parentTime,childSigma,childAlpha)
    generalizedPressSchechterProbabilityIntegrand=generalizedPressSchechterProgenitorMassFunction(childHaloMass,childSigma,childAlpha,generalizedPressSchechterNode)
    return
  end function generalizedPressSchechterProbabilityIntegrand

  double precision function generalizedPressSchechterFractionSubresolutionIntegrand(childHaloMass)
    !!{
    Integrand for the subresolution fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: childHaloMass
    double precision                :: childAlpha   , childSigma

    if (childHaloMass>0.0d0) then
       call generalizedPressSchechterSelf%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,generalizedPressSchechterSelf%parentTime,childSigma,childAlpha)
       generalizedPressSchechterFractionSubresolutionIntegrand=+generalizedPressSchechterProgenitorMassFunction(childHaloMass,childSigma,childAlpha,generalizedPressSchechterNode) &
            &                                                  *                              childHaloMass                                                                        &
            &                                                  /generalizedPressSchechterSelf%parentHaloMass
    else
       generalizedPressSchechterFractionSubresolutionIntegrand=0.0d0
    end if
    return
  end function generalizedPressSchechterFractionSubresolutionIntegrand

  double precision function generalizedPressSchechterProgenitorMassFunction(childHaloMass,childSigma,childAlpha,node)
    !!{
    Progenitor mass function from Press-Schechter.
    !!}
    implicit none
    double precision          , intent(in   ) :: childAlpha, childHaloMass, childSigma
    type            (treeNode), intent(inout) :: node

    generalizedPressSchechterProgenitorMassFunction=(generalizedPressSchechterSelf%parentHaloMass/childHaloMass**2)*generalizedPressSchechterMergingRate(childSigma,childAlpha,node)
    return
  end function generalizedPressSchechterProgenitorMassFunction

  double precision function generalizedPressSchechterMergingRate(childSigma,childAlpha,node)
    !!{
    Computes the merging rate of dark matter halos in the generalized Press-Schechter algorithm. This ``merging rate'' is specifically defined as
    \begin{equation}
    {\mathrm{d}^2 f \over \mathrm{d} \ln M_\mathrm{child} \mathrm{d} \delta_\mathrm{c}} = 2 \sigma^2(M_\mathrm{child}) \left.{\mathrm{d} \ln \sigma \over \mathrm{d} \ln M}\right|_{M=M_\mathrm{child}} {\mathrm{d}t\over \mathrm{d}\delta_\mathrm{c}} {\mathrm{d}f_{12}\over \mathrm{d}t},
    \end{equation}
    where $\mathrm{d} f_{12}/\mathrm{d}t$ is the excursion set barrier crossing probabilty per unit time for the effective barrier
    $B^\prime(S_\mathrm{child}|S_\mathrm{parent},t)\equiv B(S_\mathrm{child},t-\delta t)-B(S_\mathrm{parent},t)$ in the limit $\delta t
    \rightarrow 0$.
    !!}
    implicit none
    double precision          , intent(in   ) :: childAlpha       , childSigma
    type            (treeNode), intent(inout) :: node
    double precision                          :: childSigmaSquared

    childSigmaSquared                   =+childSigma**2
    generalizedPressSchechterMergingRate=-2.0d0                                                                                                                                &
         &                               *generalizedPressSchechterSelf%mergerTreeBranchingProbabilityModifier_%rateModifier(                                                  &
         &                                                                                                                                                 node              , &
         &                                                                                                                   generalizedPressSchechterSelf%parentHaloMass    , &
         &                                                                                                                   generalizedPressSchechterSelf%parentSigma       , &
         &                                                                                                                                                 childSigma        , &
         &                                                                                                                   generalizedPressSchechterSelf%parentTime          &
         &                                                                                                                  )                                                  &
         &                               *generalizedPressSchechterSelf%excursionSetFirstCrossing_             %rate        (                                                  &
         &                                                                                                                   generalizedPressSchechterSelf%parentSigmaSquared, &
         &                                                                                                                                                 childSigmaSquared , &
         &                                                                                                                   generalizedPressSchechterSelf%parentTime        , &
         &                                                                                                                                                 node                &
         &                                                                                                                  )                                                  &
         &                               *childSigmaSquared                                                                                                                    &
         &                               *abs(childAlpha)                                                                                                                      &
         &                               *generalizedPressSchechterSelf%parentDTimeDDeltaCritical
    return
  end function generalizedPressSchechterMergingRate

  subroutine generalizedPressSchechterComputeCommonFactors(self,node,haloMass,deltaCritical,time)
    !!{
    Precomputes some useful factors that are used in the generalized Press-Schechter branching integrals.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout)         :: self
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                 , intent(in   )         :: haloMass, deltaCritical, &
         &                                                                                      time

    generalizedPressSchechterNode  =>                                                                                                                                           node
    self%parentHaloMass            =                                                                                                                        haloMass
    self%parentDelta               =                                                                                                deltaCritical
    self%parentTime                =                                                                                                time
    self%parentSigma               =    self%cosmologicalMassVariance_%rootVariance                       (time               =self%parentTime   ,mass=     haloMass)
    self%parentSigmaSquared        =    self%parentSigma                                                                                                        **2
    ! "deltaCritical" here is actually δ_c(t) σ(M,t₀) / σ(M,t).
    self%parentDTimeDDeltaCritical = +1.0d0                                                                                                                                           &
         &                           /  self%cosmologicalMassVariance_%rootVariance                       (time               =self%parentTime   ,mass=self%parentHaloMass          ) &
         &                           *  self%cosmologicalMassVariance_%rootVariance                       (time               =self%timeNow      ,mass=self%parentHaloMass          ) &
         &                           /                                                                                              deltaCritical                                     &
         &                           /(                                                                                                                                               &
         &                             +self%criticalOverdensity_     %gradientTime                       (time               =self%parentTime   ,mass=self%parentHaloMass,node=node) &
         &                             /self%criticalOverdensity_     %value                              (time               =self%parentTime   ,mass=self%parentHaloMass,node=node) &
         &                             -self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(time               =self%parentTime   ,mass=self%parentHaloMass          ) &
         &                             /                                                                                       self%parentTime                                        &
         &                            )
    return
  end subroutine generalizedPressSchechterComputeCommonFactors

