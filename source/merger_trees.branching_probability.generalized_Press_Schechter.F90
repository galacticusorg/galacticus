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

!+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

!!{
Implements a merger tree branching probability class using a generalized Press-Schechter approach.
!!}

  use :: Cosmological_Density_Field     , only : cosmologicalMassVarianceClass              , criticalOverdensityClass
  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Excursion_Sets_First_Crossings , only : excursionSetFirstCrossingClass
  use :: Merger_Tree_Branching_Modifiers, only : mergerTreeBranchingProbabilityModifierClass
  use :: Root_Finder                    , only : rootFinder
  use :: Numerical_Integration          , only : integrator

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
     type            (integrator                                 )          :: integrator_                                         , integratorSubresolution_
     ! Parent halo shared variables.
     double precision                                                       :: parentDTimeDDeltaCritical                           , parentDelta                  , &
          &                                                                    parentHaloMass                                      , parentSigma                  , &
          &                                                                    parentSigmaSquared                                  , parentTime                   , &
          &                                                                    probabilityMinimumMass                              , probabilitySeek              , &
          &                                                                    normalization

     ! Record of mass resolution.
     double precision                                                       :: resolutionSigma                                     , massResolutionPrevious
     ! Record of parent time.
     double precision                                                       :: parentTimePrevious
     ! Accuracy parameter to ensure that steps in critical overdensity do not become too large.
     double precision                                                       :: deltaStepMaximum
     ! The maximum σ that we expect to find.
     double precision                                                       :: sigmaMaximum
     ! Record of whether we have tested the excursion set routines.
     logical                                                                :: excursionSetsTested
     ! Control for inclusion of smooth accretion rates.
     logical                                                                :: smoothAccretion
     ! Record of issued warnings.
     logical                                                                :: subresolutionFractionIntegrandFailureWarned
     ! Option controlling whether only lower-half of the distribution function should be used.
     logical                                                                :: distributionFunctionLowerHalfOnly                   , distributionFunctionNormalize
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
     procedure :: rate                  => generalizedPressSchechterRate
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
  class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), pointer   :: self_
  !$omp threadprivate(self_)

  ! Module-scope variables used in integrands.
  type            (treeNode                                       ), pointer   :: node_
  !$omp threadprivate(node_)

  ! Branching probability integrand integration tolerance.
  double precision                                                 , parameter :: toleranceIntegrandRelative=1.0d-2
  
contains

  function generalizedPressSchechterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily generalizedPressSchechter} merger tree branching probability class which reads parameters from a
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
    logical                                                                          :: smoothAccretion                        , distributionFunctionLowerHalfOnly, &
         &                                                                              distributionFunctionNormalize

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
    <inputParameter>
      <name>distributionFunctionLowerHalfOnly</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, only the lower half ($M &lt; M_0/2$) of the branching rate distribution function is used, as per the algorithm of \cite{cole_hierarchical_2000}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>distributionFunctionNormalize</name>
      <defaultValue>.true.</defaultValue>
      <description>
	If using the full range ($M &lt; M_0$) of the branching rate distribution, if this parameter is {\normalfont \ttfamily
	true} then divide the branching rate by 2. This is appropriate if two progenitors are to be sampled (i.e. a binary
	split). If the branching rate applies to only a single branch is it more appropriate to set this parameter to be
	{\normalfont \ttfamily true} in which case this normalization by a factor 2 is \emph{not} applied.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"                    name="criticalOverdensity_"                    source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"               name="cosmologicalMassVariance_"               source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                     name="cosmologyFunctions_"                     source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing"              name="excursionSetFirstCrossing_"              source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbabilityModifier" name="mergerTreeBranchingProbabilityModifier_" source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityGnrlzdPrssSchchtr(deltaStepMaximum,massMinimum,smoothAccretion,distributionFunctionLowerHalfOnly,distributionFunctionNormalize,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,excursionSetFirstCrossing_,mergerTreeBranchingProbabilityModifier_)
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

  function generalizedPressSchechterConstructorInternal(deltaStepMaximum,massMinimum,smoothAccretion,distributionFunctionLowerHalfOnly,distributionFunctionNormalize,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,excursionSetFirstCrossing_,mergerTreeBranchingProbabilityModifier_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} merger tree building class.
    !!}
    use :: Numerical_Integration, only : GSL_Integ_Gauss15
    implicit none
    type            (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr)                        :: self
    class           (cosmologicalMassVarianceClass                  ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                       ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetFirstCrossingClass                 ), intent(in   ), target :: excursionSetFirstCrossing_
    class           (mergerTreeBranchingProbabilityModifierClass    ), intent(in   ), target :: mergerTreeBranchingProbabilityModifier_
    double precision                                                 , intent(in   )         :: deltaStepMaximum                              , massMinimum
    logical                                                          , intent(in   )         :: smoothAccretion                               , distributionFunctionLowerHalfOnly       , &
         &                                                                                      distributionFunctionNormalize
    double precision                                                 , parameter             :: toleranceAbsolute                      =0.0d+0, toleranceRelative                =1.0d-9
    !![
    <constructorAssign variables="deltaStepMaximum, massMinimum, smoothAccretion, distributionFunctionLowerHalfOnly, distributionFunctionNormalize, *criticalOverdensity_, *cosmologicalMassVariance_, *cosmologyFunctions_, *excursionSetFirstCrossing_, *mergerTreeBranchingProbabilityModifier_"/>
    !!]

    self%excursionSetsTested                        =.false.
    self%subresolutionFractionIntegrandFailureWarned=.false.
    self%massResolutionPrevious                     =-1.0d0
    self%parentTimePrevious                         =-1.0d0
    self%timeNow                                    =self%cosmologyFunctions_      %cosmicTime  (                      1.0d0  )
    self%sigmaMaximum                               =self%cosmologicalMassVariance_%rootVariance(self%massMinimum,self%timeNow)
    self%finder                                     =rootFinder(                                                                           &
         &                                                      rootFunction     =generalizedPressSchechterMassBranchRoot                , &
         &                                                      toleranceAbsolute=toleranceAbsolute                                      , &
         &                                                      toleranceRelative=toleranceRelative                                        &
         &                                                     )
    self%integrator_                                =integrator(                                                                           &
         &                                                                        generalizedPressSchechterProbabilityIntegrand          , &
         &                                                      toleranceRelative=toleranceIntegrandRelative                             , &
         &                                                      integrationRule  =GSL_Integ_Gauss15                                        &
         &                                                     )
    self%integratorSubresolution_                   =integrator(                                                                           &
            &                                                                     generalizedPressSchechterFractionSubresolutionIntegrand, &
            &                                                   toleranceRelative=1.0d-3                                                 , &
            &                                                   integrationRule  =GSL_Integ_Gauss15                                        &
            &                                                  )
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
    double precision                                                                         :: massUpper                      , rootFunctionLower  , &
         &                                                                                      rootFunctionUpper
    !$GLC attributes unused :: randomNumberGenerator_

    ! Ensure excursion set calculations have sufficient range in σ.
    call self%excursionSetTest(node)
    ! Initialize global variables.
    self_                       => self
    self%probabilityMinimumMass =  massResolution
    self%probabilitySeek        =  probabilityFraction
    call self%computeCommonFactors(node,haloMass,deltaCritical,time)
    ! Determine the upper mass limit to use.
    if (self%distributionFunctionLowerHalfOnly) then
       massUpper            =+0.5d0*haloMass
       self%normalization   =+1.0d0
    else
       massUpper            =+      haloMass
       if (self%distributionFunctionNormalize) then
          self%normalization=+0.5d0
       else
          self%normalization=+1.0d0
       end if
    end if
    ! Check that the root is bracketed.
    rootFunctionLower=generalizedPressSchechterMassBranchRoot(massResolution)
    rootFunctionUpper=generalizedPressSchechterMassBranchRoot(massUpper     )
    if (rootFunctionLower*rootFunctionUpper >= 0.0d0) then
       ! Warn about this situation.
       if (displayVerbosity() >= verbosityLevelWarn) then
          message="halo branching mass root is not bracketed in generalizedPressSchechterMassBranch()"
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6,a1,e12.6)') massResolution,":",generalizedPressSchechterMassBranchRoot(massResolution)
          message=" => massMinimum:rootFunction(massMinimum) = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6,a1,e12.6)') massUpper     ,":",generalizedPressSchechterMassBranchRoot(massUpper     )
          message=" => massMaximum:rootFunction(massMaximum) = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
          write (label,'(e12.6)') probabilityFraction
          message=" =>                           probability = "//trim(label)
          call displayMessage(message,verbosityLevelWarn)
       end if
       ! If the root function is positive at half of the parent halo mass then we have a binary split.
       if (rootFunctionUpper >= 0.0d0) then
          ! Check that we are sufficiently close to zero. If we're not, it might indicate a problem.
          if     (                                                                           &
               &   generalizedPressSchechterMassBranchRoot(massUpper)                        &
               &  >                                                                          &
               &   probabilityFraction*smallProbabilityFraction                              &
               & ) call Error_Report(                                                        &
               &                     "numerical accuracy lost in root finding calculation"// &
               &                     {introspection:location}                                &
               &                    )
          ! Return a binary split mass.
          generalizedPressSchechterMassBranch=massUpper
          return
       end if
    end if
    ! Find the branch mass.
    generalizedPressSchechterMassBranch=self%finder%find(rootRange=[massResolution,massUpper],rootRangeValues=[rootFunctionLower,rootFunctionUpper])
    return
  end function generalizedPressSchechterMassBranch

  double precision function generalizedPressSchechterMassBranchRoot(massMaximum)
    !!{
    Root function used in solving for the branch mass.
    !!}
    implicit none
    double precision, intent(in   ) :: massMaximum
    
    generalizedPressSchechterMassBranchRoot=+                            self_%probabilitySeek         &
         &                                  -self_%normalization                                       &
         &                                  *self_%integrator_%integrate(                              &
         &                                                               self_%probabilityMinimumMass, &
         &                                                                                massMaximum  &
         &                                                              )
    return
  end function generalizedPressSchechterMassBranchRoot

  double precision function generalizedPressSchechterRate(self,mass,deltaCritical,time,massBranch,node)
    !!{
    Return the rate per unit mass and per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time
    {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass {\normalfont \ttfamily massBranch}.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical, mass     , &
         &                                                                                      massBranch   , time
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                                         :: massBranch_  , childSigma, &
         &                                                                                      childAlpha
    
    ! Always use the rate from the lower half of the mass range.
    if (self%distributionFunctionLowerHalfOnly) then
       if (massBranch > 0.5d0*mass) then
          massBranch_=+mass-massBranch
       else
          massBranch_=     +massBranch
       end if
    else
       massBranch_   =     +massBranch
    end if
    call self                          %computeCommonFactors              (node,mass,deltaCritical,time)
    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massBranch,self%parentTime,childSigma,childAlpha)
    self_                         => self
    generalizedPressSchechterRate =  generalizedPressSchechterProgenitorMassFunction(massBranch,childSigma,childAlpha,node)
    return
  end function generalizedPressSchechterRate

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
    Return bounds on the probability per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily
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
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical , haloMass   , &
         &                                                                                      massResolution, time
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                                         :: massMaximum   , massMinimum, &
         &                                                                                      normalization
    
    ! Ensure excursion set calculations have sufficient range in σ.
    call self%excursionSetTest(node)
    ! Get σ and δ_critical for the parent halo.
    call self%computeCommonFactors(node,haloMass,deltaCritical,time)
    massMinimum=massResolution
    if (self%distributionFunctionLowerHalfOnly) then
       massMaximum     =+0.5d0*self%parentHaloMass
       normalization   =+1.0d0
    else
       massMaximum     =+      self%parentHaloMass
       if (self%distributionFunctionNormalize) then
          normalization=+0.5d0
       else
          normalization=+1.0d0
       end if
    end if
    self_                                =>  self
    generalizedPressSchechterProbability =  +normalization                           &
         &                                  *self%integrator_%integrate(             &
         &                                                              massMinimum, &
         &                                                              massMaximum  &
         &                                                             )
    return
  end function generalizedPressSchechterProbability

  double precision function generalizedPressSchechterFractionSubresolution(self,haloMass,deltaCritical,time,massResolution,node)
    !!{
    Return the fraction of mass accreted in subresolution halos, i.e. those below {\normalfont \ttfamily massResolution}, per unit
    change in $\delta_\mathrm{crit}$ for a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily
    deltaCritical}. The integral is computed numerically.
    !!}
    use :: Display           , only : displayMagenta, displayReset
    use :: Error             , only : Warn          , errorStatusSuccess
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class           (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr), intent(inout), target :: self
    double precision                                                 , intent(in   )         :: deltaCritical                                 , haloMass      , &
         &                                                                                      massResolution                                , time
    type            (treeNode                                       ), intent(inout), target :: node
    double precision                                                 , parameter             :: resolutionSigmaOverParentSigmaTolerance=1.0d-3
    double precision                                                                         :: resolutionSigmaOverParentSigma                , integral
    integer                                                                                  :: errorStatus
    type            (varying_string                                 )                        :: message

    ! Ensure excursion set calculations have sufficient range in σ.
    call self%excursionSetTest(node)
    ! Get σ and δ_critical for the parent halo.
    call self%computeCommonFactors(node,haloMass,deltaCritical,time)
    ! Update the root-variances corresponding to the mass resolution and the minimum subresolution halo if the mass resolution
    ! or the time of parent halo changes.
    if (massResolution /= self%massResolutionPrevious .or. self%parentTime /= self%parentTimePrevious) then
       self%resolutionSigma       =self%cosmologicalMassVariance_%rootVariance(     massResolution,self%parentTime)
       self%massResolutionPrevious=massResolution
       if (self%parentTime /= self%parentTimePrevious) then
          self%sigmaMaximum       =self%cosmologicalMassVariance_%rootVariance(self%massMinimum   ,self%parentTime)
          self%parentTimePrevious =self%parentTime
       end if
    end if
    ! If requested, compute the rate of smooth accretion.
    if (self%smoothAccretion) then
       generalizedPressSchechterFractionSubresolution=+abs(self%parentDTimeDDeltaCritical)                                                                                                            &
            &                                         *    self%excursionSetFirstCrossing_             %rateNonCrossing(              self%parentSigmaSquared,self%massMinimum ,self%parentTime,node) &
            &                                         *    self%mergerTreeBranchingProbabilityModifier_%rateModifier   (node,haloMass,self%parentSigma       ,self%sigmaMaximum,self%parentTime     )
    else
       generalizedPressSchechterFractionSubresolution=0.0d0
    end if
    resolutionSigmaOverParentSigma=self%resolutionSigma/self%parentSigma
    if (resolutionSigmaOverParentSigma >= 1.0d0) then
       self_    => self
       integral =  self%integratorSubresolution_%integrate(                         &
            &                                                     self%massMinimum, &
            &                                                     massResolution  , &
            &                                              status=errorStatus       &
            &                                             )
       if (errorStatus /= errorStatusSuccess) then
          if (resolutionSigmaOverParentSigma > 1.0d0+resolutionSigmaOverParentSigmaTolerance) then
             ! Attempt the integral again with lower tolerance. Issue a warnings if this is the first time this has happened.
             if (.not.self%subresolutionFractionIntegrandFailureWarned) then
                message=displayMagenta()//'WARNING:'                                                                                                          //displayReset(  )// &
                     &                    ' Integration of the subresolution fraction in the generalized Press-Schechter branching probability module failed.'//char        (10)// &
                     &                    'Will try again with lower tolerance. This warning will not be issued again.'                                                         // &
                     &                    {introspection:location}
                call Warn(message)
                self%subresolutionFractionIntegrandFailureWarned=.true.
             end if
             call self%integratorSubresolution_%toleranceSet(toleranceRelative=1.0d-2)
             integral=self%integratorSubresolution_%integrate(self%massMinimum,massResolution)
          end if
       end if
       generalizedPressSchechterFractionSubresolution=+generalizedPressSchechterFractionSubresolution &
               &                                      +integral
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

    call self_%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,self_%parentTime,childSigma,childAlpha)
    generalizedPressSchechterProbabilityIntegrand=generalizedPressSchechterProgenitorMassFunction(childHaloMass,childSigma,childAlpha,node_)
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
       call self_%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,self_%parentTime,childSigma,childAlpha)
       generalizedPressSchechterFractionSubresolutionIntegrand=+generalizedPressSchechterProgenitorMassFunction(childHaloMass,childSigma,childAlpha,node_) &
            &                                                  *                                                childHaloMass                              &
            &                                                  /self_%parentHaloMass
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

    generalizedPressSchechterProgenitorMassFunction=(self_%parentHaloMass/childHaloMass**2)*generalizedPressSchechterMergingRate(childHaloMass,childSigma,childAlpha,node)
    return
  end function generalizedPressSchechterProgenitorMassFunction

  double precision function generalizedPressSchechterMergingRate(childHaloMass,childSigma,childAlpha,node)
    !!{
    Computes the merging rate of dark matter halos in the generalized Press-Schechter algorithm. This ``merging rate'' is specifically defined as
    \begin{equation}
    {\mathrm{d}^2 f \over \mathrm{d} \ln M_\mathrm{child} \mathrm{d} \delta_\mathrm{c}} = 2 \sigma^2(M_\mathrm{child}) \left.{\mathrm{d} \ln \sigma \over \mathrm{d} \ln M}\right|_{M=M_\mathrm{child}} {\mathrm{d}t\over \mathrm{d}\delta_\mathrm{c}} {\mathrm{d}f_{12}\over \mathrm{d}t},
    \end{equation}
    where $\mathrm{d} f_{12}/\mathrm{d}t$ is the excursion set barrier crossing probability per unit time for the effective barrier
    $B^\prime(S_\mathrm{child}|S_\mathrm{parent},t)\equiv B(S_\mathrm{child},t-\delta t)-B(S_\mathrm{parent},t)$ in the limit $\delta t
    \rightarrow 0$.
    !!}
    implicit none
    double precision          , intent(in   ) :: childAlpha       , childSigma        , &
         &                                       childHaloMass
    type            (treeNode), intent(inout) :: node
    double precision                          :: childSigmaSquared, childSigmaEffective

    childSigmaSquared=childSigma**2
    ! Decide what to use for the progenitor σ(M) in the modifier function.
    if (self_%distributionFunctionLowerHalfOnly) then
       ! Only the lower half of the distribution function is being used. Therefore, use the σ(M) corresponding to the mass of the
       ! child halo in this lower half.
       childSigmaEffective   =                                childSigma
    else
       ! The full range of the distribution function is being used. Use the σ(M) corresponding to the less massive of the two
       ! progenitors that will be generated.
       if (childHaloMass < 0.5d0*self_%parentHaloMass) then
          ! Actual progenitor halo is in the lower half of the mass range, use its mass directly to compute σ(M).
          childSigmaEffective=                                childSigma
       else
          ! Actual progenitor halo is in the upper half of the mass range, use the mass of the complement halo to compute σ(M).
          childSigmaEffective=self_%cosmologicalMassVariance_%rootVariance(self_%parentHaloMass-childHaloMass,self_%parentTime)
       end if
    end if
    generalizedPressSchechterMergingRate=-2.0d0                                                                                 &
         &                               *self_%mergerTreeBranchingProbabilityModifier_%rateModifier(                           &
         &                                                                                                 node               , &
         &                                                                                           self_%parentHaloMass     , &
         &                                                                                           self_%parentSigma        , &
         &                                                                                                 childSigmaEffective, &
         &                                                                                           self_%parentTime           &
         &                                                                                          )                           &
         &                               *self_%excursionSetFirstCrossing_             %rate        (                           &
         &                                                                                           self_%parentSigmaSquared , &
         &                                                                                                 childSigmaSquared  , &
         &                                                                                           self_%parentTime         , &
         &                                                                                                 node                 &
         &                                                                                          )                           &
         &                               *childSigmaSquared                                                                     &
         &                               *abs(childAlpha)                                                                       &
         &                               *self_%parentDTimeDDeltaCritical
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

    node_                          =>                                                                                                                           node
    self%parentHaloMass            =                                                                                                         haloMass
    self%parentDelta               =                                                                                 deltaCritical
    self%parentTime                =                                                                                 time
    self%parentSigma               =    self%cosmologicalMassVariance_%rootVariance                       (time=self%parentTime   ,mass=     haloMass                )
    self%parentSigmaSquared        =    self%parentSigma                                                                                                              **2
    ! "deltaCritical" here is actually δ_c(t) σ(M,t₀) / σ(M,t).
    self%parentDTimeDDeltaCritical = +1.0d0                                                                                                                               &
         &                           /  self%cosmologicalMassVariance_%rootVariance                       (time=self%parentTime   ,mass=self%parentHaloMass          )    &
         &                           *  self%cosmologicalMassVariance_%rootVariance                       (time=self%timeNow      ,mass=self%parentHaloMass          )    &
         &                           /                                                                               deltaCritical                                        &
         &                           /(                                                                                                                                   &
         &                             +self%criticalOverdensity_     %gradientTime                       (time=self%parentTime   ,mass=self%parentHaloMass,node=node)    &
         &                             /self%criticalOverdensity_     %value                              (time=self%parentTime   ,mass=self%parentHaloMass,node=node)    &
         &                             -self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(time=self%parentTime   ,mass=self%parentHaloMass          )    &
         &                             /                                                                        self%parentTime                                           &
         &                            )
    return
  end subroutine generalizedPressSchechterComputeCommonFactors

