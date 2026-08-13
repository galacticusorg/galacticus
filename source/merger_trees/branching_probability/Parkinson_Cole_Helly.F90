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

!!{RST
Implements a merger tree branching probability class using the algorithm of :cite:t:`parkinson_generating_2008`.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Numerical_Integration     , only : integrator
  use :: Numerical_Ranges          , only : rangeLattice
  use :: Tables                    , only : table1DLogarithmicLinear

  !![
  <mergerTreeBranchingProbability name="mergerTreeBranchingProbabilityParkinsonColeHelly" docformat="rst">
   <description>
   A merger tree branching probability class using the algorithm of :cite:t:`parkinson_generating_2008`. The parameters :math:`G_0`, :math:`\gamma_1` and :math:`\gamma_2` of their algorithm are specified by the input parameters ``[G0]``, ``[gamma]`` and ``[gamma2]`` respectively. Additionally, the parameter ``accuracyFirstOrder`` limits the step in :math:`\delta_\mathrm{crit}` so that it never exceeds ``accuracyFirstOrder``\ :math:`\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}`, which ensures the the first order expansion of the merging rate that is assumed is accurate. To find bounds on the branching probability, we make use of the fact that eqn. (4) of :cite:t:`parkinson_generating_2008` can be written as

   .. math::

      {\d f \over \d t} = {\mathrm{d} t \over \mathrm{d}\omega} \int_{M_\mathrm{min}}^{M/2} {M \over M^\prime} {\d f \over \d t}
      {\d S \over \d M^\prime} \left| {\d t \over \d \omega}\right| G[\omega,\sigma(M),\sigma(M^\prime)] \d M^\prime.

   By holding the :math:`M^\prime` in the denominator of the first term in the integrand fixed, we obtain an analytic solution to the integral in terms of hypergeometric functions. If we fix this :math:`M^\prime` at :math:`M_\mathrm{min}` we obtain an upper limit on the branching probability, while if we fix it to :math:`M/2` a lower limit is obtained.

   Calculation of branching probabilities involves computation of several hypergeometric functions which are numerically slow. Two parameters control the accuracy and application of these functions. First, ``[precisionHypergeometric]``\ (:math:`=10^{-6}`) specifies the fractional tolerance to which these functions should be computed. Second, if ``[hypergeometricTabulate]``\ :math:`=`\ ``true`` then these functions will be tabulated for rapid lookup (at some loss of precision).
   </description>
  </mergerTreeBranchingProbability>
  !!]
  type, extends(mergerTreeBranchingProbabilityClass) :: mergerTreeBranchingProbabilityParkinsonColeHelly
     !!{RST
     A merger tree branching probability class using the algorithm of :cite:t:`parkinson_generating_2008`.
     !!}
     private
     double precision                                         :: gamma1                                           , gamma2                                     , &
          &                                                      G0                                               , accuracyFirstOrder                         , &
          &                                                      precisionHypergeometric
     logical                                                  :: hypergeometricTabulate                           , cdmAssumptions                             , &
          &                                                      hypergeometricFailureWarned                      , tolerateRoundOffErrors
     type            (table1DLogarithmicLinear     )          :: subresolutionHypergeometric                      , upperBoundHypergeometric
     ! Lattices to which the tabulated abscissae are pinned. These are the source of truth for the extent of each tabulation. The
     ! abscissae of the upper bound tabulation are masses in units of twice the mass resolution - see
     ! `parkinsonColeHellyUpperBoundHypergeometricTabulate`.
     type            (rangeLattice                 )          :: latticeSubresolution                             , latticeUpperBound
     logical                                                  :: subresolutionHypergeometricInitialized =  .false., upperBoundHypergeometricInitialized=.false.
     double precision                                         :: massResolutionTabulated                          , factorG0Gamma2                             , &
          &                                                      branchingProbabilityPreFactor                    , sigmaParentSquared                         , &
          &                                                      sigmaParent                                      , deltaParent                                , &
          &                                                      massHaloParent                                   , probabilityMinimumMassLog                  , &
          &                                                      probabilityMaximumMassLog                        , probabilitySeek                            , &
          &                                                      probabilityGradientMinimum                       , probabilityGradientMaximum                 , &
          &                                                      probabilityMaximum                               , probabilityMinimumMass                     , &
          &                                                      haloMassPrevious                                 , deltaCriticalPrevious                      , &
          &                                                      massResolutionPrevious                           , probabilityPrevious                        , &
          &                                                      resolutionSigma                                  , resolutionAlpha                            , &
          &                                                      timeParent
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_              => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_                   => null()
     type            (integrator                   )          :: integrator_
   contains
     !![
     <methods docformat="rst">
       <method description="Compute common factors needed for the calculations."                         method="computeCommonFactors"/>
       <method description="Compute the function :math:`V(q)` from :cite:t:`parkinson_generating_2008`."          method="V"                   />
       <method description="Compute the part of the modifier term which depends on :math:`\sigma_\mathrm{s}`." method="modifier"            />
       <method description="Compute the :math:`a` parameter of the hypergeometric function."                   method="hypergeometricA"     />
     </methods>
     !!]
     final     ::                                 parkinsonColeHellyDestructor
     procedure :: V                            => parkinsonColeHellyV
     procedure :: modifier                     => parkinsonColeHellyModifier
     procedure :: progenitorMassFunctionFactor => parkinsonColeHellyProgenitorMassFunctionFactor
     procedure :: hypergeometricA              => parkinsonColeHellyHypergeometricA
     procedure :: rate                         => parkinsonColeHellyRate
     procedure :: probability                  => parkinsonColeHellyProbability
     procedure :: probabilityBound             => parkinsonColeHellyProbabilityBound
     procedure :: fractionSubresolution        => parkinsonColeHellyFractionSubresolution
     procedure :: massBranch                   => parkinsonColeHellyMassBranch
     procedure :: stepMaximum                  => parkinsonColeHellyStepMaximum
     procedure :: computeCommonFactors         => parkinsonColeHellyComputeCommonFactors
  end type mergerTreeBranchingProbabilityParkinsonColeHelly

  interface mergerTreeBranchingProbabilityParkinsonColeHelly
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeBranchingProbabilityParkinsonColeHelly` merger tree branching probability class.
     !!}
     module procedure parkinsonColeHellyConstructorParameters
     module procedure parkinsonColeHellyConstructorInternal
  end interface mergerTreeBranchingProbabilityParkinsonColeHelly

  ! Module-scope pointer to self used for root-finding.
  class           (mergerTreeBranchingProbabilityParkinsonColeHelly), pointer   :: self_
  !$omp threadprivate(self_)

  ! Branching probability integrand integration tolerance.
  double precision                                                  , parameter :: toleranceIntegralRelative=1.0d-3

  ! Limit on α for use in effective γ parameters.
  double precision                                                  , parameter :: alphaMinimum             =5.0d-3

contains

  function parkinsonColeHellyConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeBranchingProbabilityParkinsonColeHelly` merger tree branching probability class which reads parameters from a provided parameter list.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityParkinsonColeHelly)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (cosmologicalMassVarianceClass                   ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                        ), pointer       :: criticalOverdensity_
    double precision                                                                  :: gamma1                   , gamma2            , &
         &                                                                               G0                       , accuracyFirstOrder, &
         &                                                                               precisionHypergeometric
    logical                                                                           :: hypergeometricTabulate   , cdmAssumptions    , &
         &                                                                               tolerateRoundOffErrors

    ! Check and read parameters.
    !![
    <inputParameter docformat="rst">
      <name>G0</name>
      <defaultValue>0.57d0</defaultValue>
      <description>
      The parameter :math:`G_0` appearing in the modified merger rate expression of :cite:t:`parkinson_generating_2008`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gamma1</name>
      <defaultValue>0.38d0</defaultValue>
      <description>
      The parameter :math:`\gamma_1` appearing in the modified merger rate expression of :cite:t:`parkinson_generating_2008`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gamma2</name>
      <defaultValue>-0.01d0</defaultValue>
      <description>
      The parameter :math:`\gamma_2` appearing in the modified merger rate expression of :cite:t:`parkinson_generating_2008`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>accuracyFirstOrder</name>
      <defaultValue>0.1d0</defaultValue>
      <description>
      Limits the step in :math:`\delta_\mathrm{crit}` when constructing merger trees using the :cite:t:`parkinson_generating_2008` algorithm, so that it never exceeds ``accuracyFirstOrder``\ :math:`\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>precisionHypergeometric</name>
      <defaultValue>1.0d-6</defaultValue>
      <description>
      The fractional precision required in evaluates of hypergeometric functions in the modified Press-Schechter tree branching calculations.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>hypergeometricTabulate</name>
      <defaultValue>.true.</defaultValue>
      <description>
      Specifies whether hypergeometric factors should be precomputed and tabulated in modified Press-Schechter tree branching functions.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>cdmAssumptions</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, assume that :math:`\alpha(=-\mathrm{d}\log \sigma/\mathrm{d}\log M)&gt;0` and :math:`\mathrm{d}\alpha/\mathrm{d}M&gt;0` (as is true in the case of :term:`CDM`) when constructing merger trees using the :cite:t:`parkinson_generating_2008`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>tolerateRoundOffErrors</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, round-off errors in integrations of branching probability will be tolerated. This may degrade the accuracy of solutions, but can be unavoidable in models with cut-offs in their power spectra.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityParkinsonColeHelly(G0,gamma1,gamma2,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,tolerateRoundOffErrors,cosmologicalMassVariance_,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function parkinsonColeHellyConstructorParameters

  function parkinsonColeHellyConstructorInternal(G0,gamma1,gamma2,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,tolerateRoundOffErrors,cosmologicalMassVariance_,criticalOverdensity_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerTreeBranchingProbabilityParkinsonColeHelly` merger tree branching probability class.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : GSL_Integ_Gauss15
    implicit none
    type            (mergerTreeBranchingProbabilityParkinsonColeHelly)                        :: self
    double precision                                                  , intent(in   )         :: gamma1                   , gamma2            , &
         &                                                                                       G0                       , accuracyFirstOrder, &
         &                                                                                       precisionHypergeometric
    logical                                                           , intent(in   )         :: hypergeometricTabulate   , cdmAssumptions    , &
         &                                                                                       tolerateRoundOffErrors
    class           (cosmologicalMassVarianceClass                   ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                        ), intent(in   ), target :: criticalOverdensity_
    !![
    <constructorAssign variables="G0, gamma1, gamma2, accuracyFirstOrder, precisionHypergeometric, hypergeometricTabulate, cdmAssumptions, tolerateRoundOffErrors, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    !!]

    ! Validate inputs.
    if (gamma1 == 1.0d0) call Error_Report('γ₁=1 leads to divergent integrals'//{introspection:location})
    ! Initialize.
    self%subresolutionHypergeometricInitialized=.false.
    self%upperBoundHypergeometricInitialized   =.false.
    self%massResolutionTabulated               =-1.0d0
    self%haloMassPrevious                      =-1.0d0
    self%deltaCriticalPrevious                 =-1.0d0
    self%massResolutionPrevious                =-1.0d0
    self%probabilityPrevious                   =-1.0d0
    self%integrator_                           =integrator(                                                                     &
            &                                                                parkinsonColeHellyProbabilityIntegrandLogarithmic, &
            &                                              toleranceRelative=toleranceIntegralRelative                        , &
            &                                              integrationRule  =GSL_Integ_Gauss15                                  &
            &                                             )
    return
  end function parkinsonColeHellyConstructorInternal

  subroutine parkinsonColeHellyDestructor(self)
    implicit none
    type(mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self

    if (self%subresolutionHypergeometricInitialized) call self%subresolutionHypergeometric%destroy()
    if (self%upperBoundHypergeometricInitialized   ) call self%upperBoundHypergeometric   %destroy()
    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine parkinsonColeHellyDestructor

  double precision function parkinsonColeHellyMassBranch(self,haloMass,deltaCritical,time,massResolution,probabilityFraction,randomNumberGenerator_,node)
    !!{RST
    A merger tree branch split mass function.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout), target :: self
    double precision                                                  , intent(in   )         :: deltaCritical                 , haloMass              , &
         &                                                                                       massResolution                , probabilityFraction   , &
         &                                                                                       time
    class           (randomNumberGeneratorClass                      ), intent(inout)         :: randomNumberGenerator_
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                                          :: B                             , mu                    , &
         &                                                                                       beta                          , halfMassAlpha         , &
         &                                                                                       halfMassSigma                 , eta                   , &
         &                                                                                       halfMassV                     , massFractionResolution, &
         &                                                                                       halfPowerEta                  , x                     , &
         &                                                                                       massFraction                  , resolutionSigma       , &
         &                                                                                       massFractionResolutionPowerEta
    !$GLC attributes unused :: node

    ! Simply branch to the relevant function.
    if (self%cdmAssumptions) then
       parkinsonColeHellyMassBranch=massBranchCDMAssumptions()
    else
       parkinsonColeHellyMassBranch=massBranchGeneric       ()
    end if
    return

  contains

    double precision function massBranchCDMAssumptions()
      !!{RST
      A merger tree branch split mass function which assumes a :term:`CDM`-like power spectrum. With these assumptions, it can employ the mass sampling algorithm of :cite:t:`parkinson_generating_2008`. One difference with respect to the algorithm of :cite:t:`parkinson_generating_2008` is that here the normalization of their function :math:`S(q)` (eqn. A2) is irrelevant, since a branch split has already been decided to have occurred---all that remains necessary is to determine its mass. Variable and function names follow :cite:t:`parkinson_generating_2008`.
      !!}
      implicit none
      logical :: reject

      ! Get parent and half-mass σ and α.
      self%timeParent        =self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=     deltaCritical,mass=haloMass,node=node)
      self%sigmaParentSquared=self%cosmologicalMassVariance_%rootVariance  (time               =self%timeParent   ,mass=haloMass          )**2
      call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*haloMass,self%timeParent,halfMassSigma,halfMassAlpha)
      ! Compute parameters β, μ, and B.
      massFractionResolution=+massResolution                               &
           &                 /haloMass
      halfMassV             =+self%V(0.5d0,haloMass)
      beta                  =+log(                                         &
           &                      +self%V(massFractionResolution,haloMass) &
           &                      /halfMassV                               &
           &                     )                                         &
           &                 /log(                                         &
           &                      +massFractionResolution                  &
           &                      /0.5d0                                   &
           &                     )
      B                     =+halfMassV                                    &
           &                 *2.0d0    **beta
      if (self%gamma1 >= 0.0d0) then
         mu                 =-halfMassAlpha
      else
         resolutionSigma    =+self%cosmologicalMassVariance_%rootVariance(massResolution,self%timeParent)
         mu                 =-log(                        &
              &                   +resolutionSigma        &
              &                   /halfMassSigma          &
              &                  )                        &
              &              /log(                        &
              &                   +massFractionResolution &
              &                   /0.5d0                  &
              &                  )
      end if
      eta                           =+beta                        &
           &                         -1.0d0                       &
           &                         -mu                          &
           &                         *self%gamma1
      massFractionResolutionPowerEta=+massFractionResolution**eta
      halfPowerEta                  =+0.5d0                 **eta
      ! Sample from S(q), using rejection sampling on R(q) to decide whether to keep/reject the
      ! proposed q.
      reject=.true.
      do while (reject)
         ! Draw a random q from S(q).
         x           =randomNumberGenerator_%uniformSample()
         massFraction=(                                                 &
              &        +              massFractionResolutionPowerEta    &
              &        +(halfPowerEta-massFractionResolutionPowerEta)*x &
              &       )**(1.0d0/eta)
         x           =randomNumberGenerator_%uniformSample()
         reject=x > R(massFraction)
      end do
      massBranchCDMAssumptions=massFraction*haloMass
      return
    end function massBranchCDMAssumptions

    double precision function R(massFraction)
      !!{RST
      The function :math:`R(q)` from :cite:t:`parkinson_generating_2008`.
      !!}
      implicit none
      double precision, intent(in   ) :: massFraction
      double precision                :: massFractionSigma, massFractionAlpha

      call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massFraction*haloMass,self%timeParent,massFractionSigma,massFractionAlpha)
      if (massFractionSigma**2 <= self%sigmaParentSquared) then
         R      =+1.0d0
      else
         R      =+(                             &
              &    +massFractionAlpha           &
              &    /halfMassAlpha               &
              &   )                             &
              &  *self%V(massFraction,haloMass) &
              &  /B                             &
              &  /massFraction**beta            &
              &  *(                             &
              &    +(                           &
              &      +2.0d0                     &
              &      *massFraction              &
              &     )**mu                       &
              &    *massFractionSigma           &
              &    /halfMassSigma               &
              &   )**self%gamma1
      end if
      return
    end function R

    double precision function massBranchGeneric()
      !!{RST
      Determine the mass of one of the halos to which the given halo branches, given the branching probability, ``probability``. Typically, ``probabilityFraction`` is found by multiplying ``probability`` by a random variable drawn in the interval 0--1 if a halo branches. This routine then finds the progenitor mass corresponding to this value.
      !!}
      use :: Root_Finder, only : GSL_Root_fSolver_Brent, rootFinder
      implicit none
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-9
      type            (rootFinder), save      :: finder
      logical                     , save      :: finderConstructed=.false.
      !$omp threadprivate(finder,finderConstructed)
      double precision                        :: logMassMinimum           , logMassMaximum

      ! Initialize global variables.
      call self%computeCommonFactors(deltaCritical,time,haloMass,node)
      self_                            => self
      self  %probabilityMinimumMass    =            massResolution
      self  %probabilityMinimumMassLog =  log(      massResolution)
      self  %probabilityMaximumMassLog =  log(0.5d0*haloMass      )
      self  %probabilitySeek           =  probabilityFraction
      self  %probabilityMaximum        =  huge(0.0d0)
      ! Check the sign of the root function at half the halo mass.
      if (parkinsonColeHellyMassBranchRoot(self%probabilityMaximumMassLog) >= 0.0d0) then
         ! The root function is zero, or very close to it (which can happen due to rounding errors
         ! occasionally). Therefore we have an almost perfect binary split.
         massBranchGeneric=0.5d0*haloMass
      else
         ! Initialize our root finder.
         if (.not.finderConstructed) then
            finder           =rootFinder(                                                    &
                 &                       rootFunction     =parkinsonColeHellyMassBranchRoot, &
                 &                       toleranceAbsolute=toleranceAbsolute               , &
                 &                       toleranceRelative=toleranceRelative,                &
                 &                       solverType       =GSL_Root_fSolver_Brent            &
                 &                      )
            finderConstructed=.true.
         end if
         ! Split is not binary - seek the actual mass of the smaller progenitor.
         logMassMinimum                  =log(      massResolution)
         logMassMaximum                  =log(0.5d0*haloMass      )
         self_%probabilityGradientMinimum=parkinsonColeHellyMassBranchRootDerivative(logMassMinimum)
         self_%probabilityGradientMaximum=parkinsonColeHellyMassBranchRootDerivative(logMassMaximum)
         self_%probabilityMaximum        =parkinsonColeHellyMassBranchRoot          (logMassMaximum)
         massBranchGeneric=exp(finder%findWithFUpper(rootRange=[logMassMinimum,logMassMaximum],rootRangeValueHigh=self_%probabilityMaximum))
      end if
      return
    end function massBranchGeneric

  end function parkinsonColeHellyMassBranch

  double precision function parkinsonColeHellyV(self,massFraction,haloMass)
    !!{RST
    The function :math:`V(q)` from :cite:t:`parkinson_generating_2008`.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: massFraction     , haloMass
    double precision                                                                  :: childSigmaSquared, varianceDifference

    childSigmaSquared  =+self%cosmologicalMassVariance_%rootVariance(massFraction*haloMass,self%timeParent)**2
    ! (childSigma²-σ_p²)^{3/2} written as d·√d to use a square root instead of a general power.
    varianceDifference =+childSigmaSquared-self%sigmaParentSquared
    parkinsonColeHellyV=+childSigmaSquared/(varianceDifference*sqrt(varianceDifference))
    return
  end function parkinsonColeHellyV

  double precision function parkinsonColeHellyMassBranchRoot(logMassMaximum)
    !!{RST
    Used to find the mass of a merger tree branching event.
    !!}
    use :: Display        , only : displayGreen    , displayBlue, displayYellow, displayReset
    use :: Error          , only : errorStatusRound, errorStatusSuccess
    use :: String_Handling, only : stringXMLFormat
    implicit none
    double precision, intent(in   ) :: logMassMaximum
    double precision                :: integral      , massMaximum
    integer                         :: status

    if      (logMassMaximum < self_%probabilityMinimumMassLog) then
       parkinsonColeHellyMassBranchRoot=self_%probabilitySeek   +self_%probabilityGradientMinimum*(logMassMaximum-self_%probabilityMinimumMassLog)
    else if (logMassMaximum > self_%probabilityMaximumMassLog .and. self_%probabilityMaximum < huge(0.0d0)) then
       parkinsonColeHellyMassBranchRoot=self_%probabilityMaximum+self_%probabilityGradientMaximum*(logMassMaximum-self_%probabilityMaximumMassLog)
    else
       massMaximum=+exp(logMassMaximum)
       integral   =+self_%branchingProbabilityPreFactor                                                                         &
            &      *self_%integrator_                  %integrate(self_%probabilityMinimumMassLog,logMassMaximum,status=status)
       if (.not.(status == errorStatusSuccess .or. (status == errorStatusRound .and. self_%tolerateRoundOffErrors))) then
          if (status == errorStatusRound) then
             call Error_Report(                                                                                                                                                                                                                           &
                  &            'probability integral failed to converge due to round-off errors - this can happen below the cut off scale in truncated power spectra'                                                              //char(10)//           &
                  &             displayGreen()//'HELP:'//displayReset()//' set the highlighted option in your input parameter file as shown below:'                                                                                //char(10)//char(10)// &
                  &             stringXMLFormat('<mergerTreeBranchingProbability value="'//char(self_%objectType(short=.true.))//'">**B<tolerateRoundOffErrors value="true"/>**C</mergerTreeBranchingProbability>',indentInitial=3)//char(10)//char(10)// &
                  &             'to ignore round-off errors and proceed'//{introspection:location}                                                                                                                                                        &
                  &           )
          else
             call Error_Report('probability integral failed to converge'//{introspection:location})
          end if
       end if
       parkinsonColeHellyMassBranchRoot=self_%probabilitySeek-integral
    end if
    return
  end function parkinsonColeHellyMassBranchRoot

  double precision function parkinsonColeHellyMassBranchRootDerivative(logMassMaximum)
    !!{RST
    Used to find the mass of a merger tree branching event.
    !!}
    implicit none
    double precision, intent(in   ) :: logMassMaximum
    double precision                :: integral

    integral=+self_%branchingProbabilityPreFactor                                                    &
         &   *parkinsonColeHellyProbabilityIntegrandLogarithmic(                                     &
         &                                                      max(                                 &
         &                                                          logMassMaximum                 , &
         &                                                          self_%probabilityMinimumMassLog  &
         &                                                         )                                 &
         &                                                     )
    parkinsonColeHellyMassBranchRootDerivative=-integral
    return
  end function parkinsonColeHellyMassBranchRootDerivative

  double precision function parkinsonColeHellyStepMaximum(self,haloMass,deltaCritical,time,massResolution)
    !!{RST
    Return the maximum allowed step in :math:`\delta_\mathrm{crit}` that a halo of mass ``haloMass`` at time ``deltaCritical`` should be allowed to take.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: deltaCritical             , haloMass   , &
         &                                                                               massResolution            , time
    double precision                                                  , parameter     :: largeStep          =1.0d10                 !   Effectively infinitely large step in w(=delta_crit).
    double precision                                                                  :: parentHalfMassSigma       , parentSigma, &
         &                                                                               varianceResidual
    !$GLC attributes unused :: deltaCritical, time

    ! Get σ and δ_critical for the parent halo.
    if (haloMass > 2.0d0*massResolution) then
       parentSigma                  =+self%cosmologicalMassVariance_%rootVariance(      haloMass,self%timeParent)
       parentHalfMassSigma          =+self%cosmologicalMassVariance_%rootVariance(0.5d0*haloMass,self%timeParent)
       varianceResidual             =+parentHalfMassSigma**2 &
            &                        -parentSigma        **2
       if (varianceResidual > 0.0d0) then
          parkinsonColeHellyStepMaximum=+self%accuracyFirstOrder &
               &                        *sqrt(                   &
               &                              +2.0d0             &
               &                              *varianceResidual  &
               &                             )
       else
          parkinsonColeHellyStepMaximum=largeStep
       end if
    else
       parkinsonColeHellyStepMaximum=largeStep
    end if
    return
  end function parkinsonColeHellyStepMaximum

  double precision function parkinsonColeHellyRate(self,mass,deltaCritical,time,massBranch,node)
    !!{RST
    Return the rate per unit mass and per unit change in :math:`\delta_\mathrm{crit}` that a halo of mass ``haloMass`` at time ``deltaCritical`` will undergo a branching to progenitors with mass ``massBranch``.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout), target :: self
    double precision                                                  , intent(in   )         :: deltaCritical , mass, &
         &                                                                                       massBranch    , time
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                                          :: massBranch_
    
    ! Always use the rate from the lower half of the mass range.
    if (massBranch > 0.5d0*mass) then
       massBranch_=+mass-massBranch
    else
       massBranch_=     +massBranch
    end if
    call self%computeCommonFactors(deltaCritical,time,mass,node)
    self_                  =>  self
    parkinsonColeHellyRate =  +self%branchingProbabilityPreFactor                                  &
         &                    *parkinsonColeHellyProbabilityIntegrandLogarithmic(log(massBranch_)) &
         &                    /massBranch_
    return
  end function parkinsonColeHellyRate

  double precision function parkinsonColeHellyProbability(self,haloMass,deltaCritical,time,massResolution,node)
    !!{RST
    Return the probability per unit change in :math:`\delta_\mathrm{crit}` that a halo of mass ``haloMass`` at time ``deltaCritical`` will undergo a branching to progenitors with mass greater than ``massResolution``.
    !!}
    use :: Display        , only : displayGreen    , displayBlue             , displayYellow     , displayReset
    use :: Error          , only : errorStatusRound, errorStatusMaxIterations, errorStatusSuccess
    use :: String_Handling, only : stringXMLFormat
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout), target :: self
    double precision                                                  , intent(in   )         :: deltaCritical , haloMass   , &
         &                                                                                       massResolution, time
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                                          :: massMaximum   , massMinimum
    integer                                                                                   :: status
    !$GLC attributes unused :: node

    ! Recompute branching probability if necessary.
    if     (                                               &
         &   haloMass       /= self%haloMassPrevious       &
         &  .or.                                           &
         &   deltaCritical  /= self%deltaCriticalPrevious  &
         &  .or.                                           &
         &   massResolution /= self%massResolutionPrevious &
         & ) then
       self_                        => self
       self %haloMassPrevious       =  haloMass
       self %deltaCriticalPrevious  =  deltaCritical
       self %massResolutionPrevious =  massResolution
       ! Get σ and δ_critical for the parent halo.
       if (haloMass > 2.0d0*massResolution) then
          call self%computeCommonFactors(deltaCritical,time,haloMass,node)
          massMinimum             =+           massResolution
          massMaximum             =+0.5d0*self%massHaloParent
          self%probabilityPrevious=+self%branchingProbabilityPreFactor                                    &
               &                   *self%integrator_                  %integrate(                         &
               &                                                                        log(massMinimum), &
               &                                                                        log(massMaximum), &
               &                                                                 status=status            &
               &                                                                )
          if (.not.(status == errorStatusSuccess .or. ((status == errorStatusRound .or. status == errorStatusMaxIterations) .and. self_%tolerateRoundOffErrors))) then
             if    (status == errorStatusRound          ) then
                call Error_Report(                                                                                                                                                                                                                          &
                     &            'probability integral failed to converge due to round-off errors - this can happen below the cut off scale in truncated power spectra'                                                             //char(10)//           &
                     &             displayGreen()//'HELP:'//displayReset()//' set the highlighted option in your input parameter file as shown below:'                                                                               //char(10)//char(10)// &
                     &             stringXMLFormat('<mergerTreeBranchingProbability value="'//char(self%objectType(short=.true.))//'">**B<tolerateRoundOffErrors value="true"/>**C</mergerTreeBranchingProbability>',indentInitial=3)//char(10)//char(10)// &
                     &             'to ignore round-off errors and proceed'//{introspection:location}                                                                                                                                                       &
                     &           )
             else if (status == errorStatusMaxIterations) then
                call Error_Report(                                                                                                                                                                                                                          &
                     &            'probability integral failed to converge due to exceeding the maximum number of iterations - this can happen below the cut off scale in truncated power spectra'                                   //char(10)//           &
                     &             displayGreen()//'HELP:'//displayReset()//' set the highlighted option in your input parameter file as shown below:'                                                                               //char(10)//char(10)// &
                     &             stringXMLFormat('<mergerTreeBranchingProbability value="'//char(self%objectType(short=.true.))//'">**B<tolerateRoundOffErrors value="true"/>**C</mergerTreeBranchingProbability>',indentInitial=3)//char(10)//char(10)// &
                     &             'to ignore round-off errors and proceed'//{introspection:location}                                                                                                                                                       &
                     &           )
             else 
                call Error_Report('probability integral failed to converge'//{introspection:location})
             end if
          end if
       else
          self%probabilityPrevious=0.0d0
       end if
    end if
    parkinsonColeHellyProbability=self%probabilityPrevious
    return
  end function parkinsonColeHellyProbability

  double precision function parkinsonColeHellyProbabilityIntegrandLogarithmic(logChildHaloMass)
    !!{RST
    Integrand for the branching probability.
    !!}
    implicit none
    double precision, intent(in   ) :: logChildHaloMass
    double precision                :: childAlpha      , childSigma, &
         &                             childHaloMass

    childHaloMass=exp(logChildHaloMass)
    call self_%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,self_%timeParent,childSigma,childAlpha)
    parkinsonColeHellyProbabilityIntegrandLogarithmic=+parkinsonColeHellyProgenitorMassFunction(childHaloMass,childSigma,childAlpha)&
         &                                            *                                         childHaloMass
    return
  end function parkinsonColeHellyProbabilityIntegrandLogarithmic

  double precision function parkinsonColeHellyProgenitorMassFunction(childHaloMass,childSigma,childAlpha)
    !!{RST
    Progenitor mass function from Press-Schechter. The constant factor of the parent halo mass is not included here---instead it is included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
    implicit none
    double precision, intent(in   ) :: childAlpha, childHaloMass, childSigma

    ! The merging-rate × modifier product is obtained through a bindable function so that subclasses (e.g.
    ! PCHPlus) can supply an algebraically-fused form that avoids redundant power evaluations.
    parkinsonColeHellyProgenitorMassFunction=+self_%progenitorMassFunctionFactor(childSigma,childAlpha) &
         &                                   /childHaloMass**2
    return
  end function parkinsonColeHellyProgenitorMassFunction

  double precision function parkinsonColeHellyProgenitorMassFunctionFactor(self,childSigma,childAlpha)
    !!{RST
    The unnormalized progenitor mass function factor---the merging rate multiplied by the modifier. Provided
    as a bindable function so that subclasses can override it with an algebraically-fused form.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: childSigma, childAlpha

    parkinsonColeHellyProgenitorMassFunctionFactor=+parkinsonColeHellyMergingRate(childSigma,childAlpha)*self%modifier(childSigma)
    return
  end function parkinsonColeHellyProgenitorMassFunctionFactor

  double precision function parkinsonColeHellyMergingRate(childSigma,childAlpha)
    !!{RST
    Merging rate from Press-Schechter. The constant factor of :math:`\sqrt{2/\pi}` not included here---instead it is included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
    implicit none
    double precision, intent(in   ) :: childAlpha       , childSigma
    double precision                :: childSigmaSquared, varianceDifference

    childSigmaSquared=childSigma**2
    if (childSigmaSquared > self_%sigmaParentSquared .and. childAlpha < 0.0d0) then
       ! The denominator (childSigma²-σ_p²)^{3/2} is written as d·√d to use a (cheaper) square root in place
       ! of a general power.
       varianceDifference           =childSigmaSquared-self_%sigmaParentSquared
       parkinsonColeHellyMergingRate=(childSigmaSquared/(varianceDifference*sqrt(varianceDifference)))*abs(childAlpha)
    else
       parkinsonColeHellyMergingRate=0.0d0
    end if
    return
  end function parkinsonColeHellyMergingRate

  double precision function parkinsonColeHellyModifier(self,childSigma)
    !!{RST
    Empirical modification of the progenitor mass function from :cite:t:`parkinson_generating_2008`. The constant factors of :math:`G_0 (\delta_\mathrm{p}/\sigma_\mathrm{p})^{\gamma_2}` and :math:`1/\sigma_\mathrm{p}^{\gamma_1}` are not included here---instead they are included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: childSigma

    parkinsonColeHellyModifier=childSigma**self%gamma1
    return
  end function parkinsonColeHellyModifier

  function parkinsonColeHellyHypergeometricA(self,gamma) result(a)
    !!{RST
    Compute the :math:`a` parameter of the hypergeometric function.
    !!}
    implicit none
    double precision                                                  , dimension(2)  :: a
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: gamma
    
    a=[1.5d0,0.5d0-0.5d0*gamma]
    return
  end function parkinsonColeHellyHypergeometricA
  
  double precision function parkinsonColeHellyProbabilityBound(self,haloMass,deltaCritical,time,massResolution,bound,node)
    !!{RST
    Return a bound on the probability per unit change in :math:`\delta_\mathrm{crit}` that a halo of mass ``haloMass`` at time ``deltaCritical`` will undergo a branching to progenitors with mass greater than ``massResolution``.
    !!}
    use            :: Display                 , only : displayMessage    , verbosityLevelWarn, displayMagenta, displayReset
    use            :: Error                   , only : Error_Report
    use            :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use, intrinsic :: ISO_C_Binding           , only : c_int
    use            :: Interface_GSL           , only : GSL_Success
    use            :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout)         :: self
    double precision                                                  , intent(in   )         :: deltaCritical                                , haloMass                 , &
         &                                                                                       massResolution                               , time
    type            (enumerationMergerTreeBranchingBoundType         ), intent(in   )         :: bound
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                  , parameter             :: sqrtTwoOverPi                 =sqrt(2.0d0/Pi)
    double precision                                                                          :: probabilityIntegrandLower                    , probabilityIntegrandUpper, &
         &                                                                                       halfParentSigma                              , halfParentAlpha          , &
         &                                                                                       gammaEffective
    double precision                                                                          :: hyperGeometricFactorLower                    , hyperGeometricFactorUpper, &
         &                                                                                       resolutionSigmaOverParentSigma
    integer         (c_int                                           )                        :: statusLower                                  , statusUpper
    logical                                                                                   :: usingCDMAssumptions
    integer                                                                                   :: iBound
    !$GLC attributes unused :: node

    ! Get σ and δ_critical for the parent halo.
    if (haloMass <= 2.0d0*massResolution) then
       parkinsonColeHellyProbabilityBound=0.0d0
       return
    end if
    call self%computeCommonFactors(deltaCritical,time,haloMass,node)
    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,self%timeParent,self%resolutionSigma,self%resolutionAlpha)
    if (massResolution /= self%massResolutionTabulated .or. self%cosmologicalMassVariance_%growthIsMassDependent()) then
       ! Resolution changed - recompute σ and α at resolution limit. Also reset the hypergeometric factor tables since
       ! these depend on resolution.
       self%upperBoundHypergeometricInitialized=.false.
    end if
    resolutionSigmaOverParentSigma=self%resolutionSigma/self%sigmaParent
    ! Estimate probability.
    if (resolutionSigmaOverParentSigma <= 1.0d0) then
       parkinsonColeHellyProbabilityBound=-1.0d0
       return
    end if
    ! Compute relevant σ and α.
    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*self%massHaloParent,self%timeParent,halfParentSigma,halfParentAlpha)
    if (halfParentSigma <= self%sigmaParent) then
       parkinsonColeHellyProbabilityBound=-1.0d0
       return
    end if
    ! Iterative over available bounds.
    parkinsonColeHellyProbabilityBound=0.0d0
    do iBound=1,2
       ! Determine if CDM assumptions can be used. Do this only if these have been explicitly allowed, if this is our first
       ! pass through the bounds evaluation, and if both αs are sufficiently large. (This last condition is required
       ! since we raise quantities to the power of 1/α which can cause problems for very small α.)
       usingCDMAssumptions= self%cdmAssumptions                       &
            &              .and.                                      &
            &               iBound                    == 1            &
            &              .and.                                      &
            &               abs(self%resolutionAlpha) >  alphaMinimum &
            &              .and.                                      &
            &               abs(     halfParentAlpha) >  alphaMinimum
       ! Compute the effective value of γ.
       gammaEffective=self%gamma1
       if (usingCDMAssumptions) then
          select case (bound%ID)
          case (mergerTreeBranchingBoundLower%ID)
             gammaEffective=gammaEffective-1.0d0/self%resolutionAlpha
          case (mergerTreeBranchingBoundUpper%ID)
             gammaEffective=gammaEffective-1.0d0/     halfParentAlpha
          end select
       end if
       ! Compute probability factors. The logic here becomes complicated, as we use various optimizations and tabulations to
       ! speed up calculation.
       !
       ! Tabulations will only be used if self%tabulateHypergeometric is true.
       !
       ! Set status to success by default.
       statusLower=GSL_Success
       statusUpper=GSL_Success
       ! First, check if CDM assumptions are not being used and we're allowed to tabulate hypergeometric factors,
       if (.not.usingCDMAssumptions.and.self%hypergeometricTabulate) then
          ! CDM assumptions are not being used. In this case we can use the same table of hypergeometric factors as the
          ! subresolution merger fraction.
          call parkinsonColeHellySubresolutionHypergeometricTabulate(self,resolutionSigmaOverParentSigma)
          call parkinsonColeHellySubresolutionHypergeometricTabulate(self,halfParentSigma   /self%sigmaParent)
          probabilityIntegrandLower=+self%factorG0Gamma2*self%subresolutionHypergeometric%interpolate(+resolutionSigmaOverParentSigma  -1.0d0)/self%sigmaParent
          probabilityIntegrandUpper=+self%factorG0Gamma2*self%subresolutionHypergeometric%interpolate(+halfParentSigma/self%sigmaParent-1.0d0)/self%sigmaParent
       else
          ! Next, check if CDM assumptions are being used, we're allowed to tabulate hypergeometric factors, and the bound
          ! requested is the upper bound.
          if     ( usingCDMAssumptions                    &
               &  .and.                                   &
               &   self%hypergeometricTabulate            &
               &  .and.                                   &
               &   bound == mergerTreeBranchingBoundUpper &
               & ) then
             ! Use a tabulation of the hypergeometric functions for the upper bound, made using CDM assumptions. Since the
             ! tables already include the difference between the upper and lower integrand, we simply set the lower
             ! integrand to zero here.
             call parkinsonColeHellyUpperBoundHypergeometricTabulate(self,self%massHaloParent,massResolution)
             ! The upper bound tabulation is made against mass in units of twice the mass resolution - see
             ! `parkinsonColeHellyUpperBoundHypergeometricTabulate`.
             probabilityIntegrandUpper=self%factorG0Gamma2*self%upperBoundHypergeometric%interpolate(self%massHaloParent/(2.0d0*self%massResolutionTabulated))/self%resolutionSigma
             probabilityIntegrandLower=0.0d0
          else
             ! Use a direct calculation of the hypergeometric factors in this case.
             hyperGeometricFactorLower=Hypergeometric_2F1(                                                           &
                  &                                                         self%hypergeometricA(gammaEffective)   , &
                  &                                                         [      1.5d0-0.5d0*gammaEffective]     , &
                  &                                                         1.0d0/resolutionSigmaOverParentSigma**2, &
                  &                                       toleranceRelative=self%precisionHypergeometric           , &
                  &                                       status           =statusLower                              &
                  &                                      )
             if (statusLower /= GSL_Success) then
                if (usingCDMAssumptions) then
                   if (.not.self%hypergeometricFailureWarned) then
                      self%hypergeometricFailureWarned=.true.
                      call displayMessage(                                                                                                                      &
                           &              displayMagenta()//'WARNING:'//displayReset()//' hypergeometric function evaluation failed when computing'//char(10)// &
                           &              'merger tree branching probability bounds - will revert to more'                                         //char(10)// &
                           &              'robust (but less stringent) bound in this and future cases'                                                       ,  &
                           &              verbosityLevelWarn                                                                                                    &
                           &             )
                   end if
                   cycle
                else
                   parkinsonColeHellyProbabilityBound=0.0d0
                   call Error_Report('hypergeometric function evaluation failed'//{introspection:location})
                end if
             end if
             probabilityIntegrandLower=+sqrtTwoOverPi                                            &
                  &                    *(self%factorG0Gamma2/self%sigmaParent)                   &
                  &                    *(resolutionSigmaOverParentSigma**(gammaEffective-1.0d0)) &
                  &                    /(1.0d0-gammaEffective)                                   &
                  &                    *hyperGeometricFactorLower
             ! Check if we can use a table to compute the upper factor.
             hyperGeometricFactorUpper=Hypergeometric_2F1(                                                          &
                  &                                                         self%hypergeometricA(gammaEffective)  , &
                  &                                                         [      1.5d0-0.5d0*gammaEffective]    , &
                  &                                                         self%sigmaParent**2/halfParentSigma**2, &
                  &                                       toleranceRelative=self%precisionHypergeometric          , &
                  &                                       status           =statusUpper                             &
                  &                                      )
             if (statusUpper /= GSL_Success) then
                if (usingCDMAssumptions) then
                   if (.not.self%hypergeometricFailureWarned) then
                      self%hypergeometricFailureWarned=.true.
                      call displayMessage(                                                                                                                      &
                           &              displayMagenta()//'WARNING:'//displayReset()//' hypergeometric function evaluation failed when computing'//char(10)// &
                           &              'merger tree branching probability bounds - will revert to more'                                         //char(10)// &
                           &              'robust (but less stringent) bound in this and future cases'                                                       ,  &
                           &              verbosityLevelWarn                                                                                                    &
                           &             )
                   end if
                   cycle
                else
                   parkinsonColeHellyProbabilityBound=0.0d0
                   call Error_Report('hypergeometric function evaluation failed'//{introspection:location})
                end if
             end if
             probabilityIntegrandUpper=+sqrtTwoOverPi                                                &
                  &                    *(self%factorG0Gamma2/self%sigmaParent)                       &
                  &                    *((halfParentSigma/self%sigmaParent)**(gammaEffective-1.0d0)) &
                  &                    /(1.0d0-gammaEffective)                                       &
                  &                    *hyperGeometricFactorUpper
          end if
       end if
       ! Compute the bound.
       select case (bound%ID)
       case (mergerTreeBranchingBoundLower%ID)
          if (usingCDMAssumptions) then
             parkinsonColeHellyProbabilityBound=+(                               &
                  &                               +probabilityIntegrandUpper     &
                  &                               -probabilityIntegrandLower     &
                  &                              )                               &
                  &                             *self%massHaloParent             &
                  &                             /massResolution                  &
                  &                             *(                               &
                  &                               +self%resolutionSigma          &
                  &                               /self%sigmaParent              &
                  &                              )**(1.0d0/self%resolutionAlpha)
          else
             parkinsonColeHellyProbabilityBound=+(                               &
                  &                               +probabilityIntegrandUpper     &
                  &                               -probabilityIntegrandLower     &
                  &                              )                               &
                  &                             *       self%massHaloParent      &
                  &                             /(0.5d0*self%massHaloParent)
          end if
       case (mergerTreeBranchingBoundUpper%ID)
          if (usingCDMAssumptions) then
             parkinsonColeHellyProbabilityBound=+(                               &
                  &                               +probabilityIntegrandUpper     &
                  &                               -probabilityIntegrandLower     &
                  &                              )                               &
                  &                             *self%massHaloParent             &
                  &                             /massResolution                  &
                  &                             *(                               &
                  &                               +self%resolutionSigma          &
                  &                               /self%sigmaParent              &
                  &                              )**(1.0d0/halfParentAlpha)
          else
             parkinsonColeHellyProbabilityBound=+(                               &
                  &                               +probabilityIntegrandUpper     &
                  &                               -probabilityIntegrandLower     &
                  &                              )                               &
                  &                             *self%massHaloParent             &
                  &                             /massResolution
          end if
       case default
          parkinsonColeHellyProbabilityBound=-1.0d0
          call Error_Report('unknown bound type'//{introspection:location})
       end select
       if (statusUpper == GSL_Success .and. statusLower == GSL_Success) exit
    end do     
    return
  end function parkinsonColeHellyProbabilityBound

  double precision function parkinsonColeHellyFractionSubresolution(self,haloMass,deltaCritical,time,massResolution,node)
    !!{RST
    Return the fraction of mass accreted in subresolution halos, i.e. those below ``massResolution``, per unit change in :math:`\delta_\mathrm{crit}` for a halo of mass ``haloMass`` at time ``deltaCritical``. The integral is computed analytically in terms of the :math:`_2F_1` hypergeometric function.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout), target :: self
    double precision                                                  , intent(in   )         :: deltaCritical                      , haloMass                      , &
         &                                                                                       massResolution                     , time
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                  , parameter             :: sqrtTwoOverPi       =sqrt(2.0d0/Pi)
    double precision                                                                          :: hyperGeometricFactor               , resolutionSigmaOverParentSigma, &
         &                                                                                       resolutionSigma
    !$GLC attributes unused :: node

    ! Get σ and δ_critical for the parent halo.
    call self%computeCommonFactors(deltaCritical,time,haloMass,node)
    resolutionSigma               =self%cosmologicalMassVariance_%rootVariance(massResolution,self%timeParent)
    resolutionSigmaOverParentSigma=resolutionSigma/self%sigmaParent
    if (resolutionSigmaOverParentSigma > 1.0d0) then
       if (self%hypergeometricTabulate) then
          ! Use tabulation of hypergeometric factors.
          call parkinsonColeHellySubresolutionHypergeometricTabulate(self,resolutionSigmaOverParentSigma)
          parkinsonColeHellyFractionSubresolution=+self%factorG0Gamma2                                                          &
               &                                  *self%subresolutionHypergeometric%interpolate(                                &
               &                                                                                +resolutionSigmaOverParentSigma &
               &                                                                                -1.0d0                          &
               &                                                                               )                                &
               &                                  /self%sigmaParent
       else
          ! Compute hypergeometric factors directly.
          hyperGeometricFactor=Hypergeometric_2F1(                                                           &
               &                                                    self%hypergeometricA(self%gamma1)      , &
               &                                                    [      1.5d0-0.5d0*self%gamma1]        , &
               &                                                    1.0d0/resolutionSigmaOverParentSigma**2, &
               &                                  toleranceRelative=self%precisionHypergeometric             &
               &                                 )
          parkinsonColeHellyFractionSubresolution=+sqrtTwoOverPi                                             &
               &                                  *self%factorG0Gamma2                                       &
               &                                  /self%sigmaParent                                          &
               &                                  *resolutionSigmaOverParentSigma**(+self%gamma1-1.0d0)      &
               &                                  /                                (-self%gamma1+1.0d0)      &
               &                                  *hyperGeometricFactor
       end if
    else
       parkinsonColeHellyFractionSubresolution=-1.0d0
    end if
    return
  end function parkinsonColeHellyFractionSubresolution

  subroutine parkinsonColeHellyComputeCommonFactors(self,deltaParent,time,massHaloParent,node)
    !!{RST
    Precomputes some useful factors that are used in the modified Press-Schechter branching integrals.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout) :: self
    double precision                                                  , intent(in   ) :: deltaParent                 , massHaloParent, &
         &                                                                               time
    type            (treeNode                                        ), intent(inout) :: node
    double precision                                                  , parameter     :: sqrtTwoOverPi=sqrt(2.0d0/Pi)

    self%deltaParent                  =                                                                       deltaParent
    self%massHaloParent               =                                                                                             massHaloParent
    self%timeParent                   =                                                                       time
    self%sigmaParent                  =self%cosmologicalMassVariance_%rootVariance  (time               =self%timeParent ,mass=self%massHaloParent          )
    self%sigmaParentSquared           =self%sigmaParent**2
    self%factorG0Gamma2               =self%G0*((max(self%criticalOverdensity_%value(time=self%timeParent,mass=self%massHaloParent,node=node),0.0d0)/self%sigmaParent)**self%gamma2)
    self%branchingProbabilityPreFactor=sqrtTwoOverPi*self%massHaloParent*self%factorG0Gamma2/self%sigmaParent**self%gamma1
    return
  end subroutine parkinsonColeHellyComputeCommonFactors

  double precision function parkinsonColeHellyHypergeometricFactor(x,gamma,toleranceRelative) result(factor)
    !!{RST
    Evaluate the factor

    .. math::

       (1+x)^{\gamma-1} \, {}_2F_1\left(\frac{3}{2},\frac{1-\gamma}{2};\frac{3-\gamma}{2};\frac{1}{(1+x)^2}\right) \Big/ (1-\gamma)

    which appears in the subresolution merger fraction, as a function of :math:`x`---the excess of the ratio of
    :math:`\sigma(M_\mathrm{res})` to :math:`\sigma(M_\mathrm{parent})` over unity.

    Evaluating this directly is impossible for small :math:`x`. The argument of the hypergeometric function,
    :math:`z=1/(1+x)^2`, approaches its singular point :math:`z=1` as :math:`x \rightarrow 0`, and since
    :math:`1+x` cannot be distinguished from :math:`1` once :math:`x` falls below the machine epsilon, :math:`z`
    becomes *exactly* :math:`1` and the evaluation fails outright with a domain error. Well before that it loses
    precision: what the function depends on is :math:`1-z`, which is formed by a subtraction from unity and so
    carries an absolute error of order the machine epsilon---a relative error of order
    :math:`\epsilon/2x`, which reaches a part in :math:`10^8` by :math:`x=10^{-8}`.

    The cure is to evaluate in terms of :math:`w \equiv 1-z = x(2+x)/(1+x)^2`, which suffers no such cancellation,
    using the connection formula between :math:`z` and :math:`1-z` (:cite:t:`olver_nist_2010`, eq. 15.8.4). Here
    :math:`c-a-b=-1/2` exactly---never an integer, so that formula is always valid---and two simplifications
    follow: the first of its two hypergeometric functions has equal first numerator and denominator parameters, so
    it collapses to :math:`(1-w)^{(\gamma-1)/2}=(1+x)^{1-\gamma}` and cancels the prefactor entirely, leaving a
    constant; and the second reduces to :math:`{}_2F_1(-\gamma/2,1;1/2;w)`, whose series in :math:`w` has the
    simple ratio of successive terms used below. The result is

    .. math::

       \frac{-2\sqrt{\pi}\,\Gamma\left(\frac{3-\gamma}{2}\right)}{(1-\gamma)\Gamma\left(-\frac{\gamma}{2}\right)}
       + \frac{(1+x)^{\gamma-1}}{\sqrt{w}} \, {}_2F_1\left(-\frac{\gamma}{2},1;\frac{1}{2};w\right),

    which is exact, not asymptotic. The series is used below :math:`x=1/10`, where it needs some twenty terms and
    agrees with the direct evaluation to a part in :math:`10^{15}`; above that the direct evaluation is used, as
    the series converges ever more slowly while the direct evaluation is at its most accurate.
    !!}
    use :: Gamma_Functions         , only : Gamma_Function
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: x                         , gamma         , &
         &                             toleranceRelative
    ! Value of x below which the series in w is used in place of a direct evaluation.
    double precision, parameter     :: xSeriesMaximum    =0.1d+00
    ! Maximum number of terms to accumulate in that series. Twenty are needed at the switch-over, so this is never
    ! approached - it exists only to bound the loop.
    integer         , parameter     :: countTermsMaximum =1000
    ! Tolerance within which -γ/2 is judged to sit on a pole of the Γ function.
    double precision, parameter     :: toleranceGammaPole=1.0d-12
    double precision                :: w                         , term          , &
         &                             seriesSum                 , factorConstant
    integer                         :: i

    if (x > xSeriesMaximum) then
       ! Away from the singular point evaluate directly.
       factor=+(1.0d0+x)**(+gamma-1.0d0)                                       &
            & /             (-gamma+1.0d0)                                     &
            & *Hypergeometric_2F1(                                             &
            &                                       [1.5d0,0.5d0-0.5d0*gamma], &
            &                                       [      1.5d0-0.5d0*gamma], &
            &                                       1.0d0/(1.0d0+x)**2       , &
            &                     toleranceRelative=toleranceRelative          &
            &                    )
    else
       ! Close to the singular point evaluate through the connection formula, in terms of w=1-z, which is free of
       ! the cancellation that destroys z.
       w      =+       x     &
            &  *(2.0d0+x)    &
            &  /(1.0d0+x)**2
       ! Accumulate ₂F₁(-γ/2,1;1/2;w). Since the second numerator parameter is unity the factorial in the series
       ! cancels the Pochhammer symbol it would otherwise carry, leaving this ratio of successive terms.
       term     =1.0d0
       seriesSum=1.0d0
       do i=0,countTermsMaximum-1
          term     =+term                   &
               &    *(-0.5d0*gamma+dble(i)) &
               &    /(+0.5d0      +dble(i)) &
               &    *w
          seriesSum=seriesSum+term
          if (abs(term) <= toleranceRelative*abs(seriesSum)) exit
       end do
       ! The constant term. Γ(-γ/2) has poles at γ=0,2,4,…, at which this term vanishes.
       if (gamma >= -toleranceGammaPole .and. abs(gamma-2.0d0*dble(nint(gamma/2.0d0))) < toleranceGammaPole) then
          factorConstant=+0.0d0
       else
          factorConstant=-2.0d0                             &
               &         *sqrt(Pi)                          &
               &         *Gamma_Function(1.5d0-0.5d0*gamma) &
               &         /Gamma_Function(     -0.5d0*gamma) &
               &         /             (-gamma+1.0d0)
       end if
       factor=+factorConstant            &
            & +(1.0d0+x)**(+gamma-1.0d0) &
            & *seriesSum                 &
            & /sqrt(w)
    end if
    return
  end function parkinsonColeHellyHypergeometricFactor

  subroutine parkinsonColeHellySubresolutionHypergeometricTabulate(self,x,xMinimumIn,xMaximumIn)
    !!{RST
    Tabulate the hypergeometric term appearing in the subresolution merger fraction expression.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Range_Pinned          , gridSchemePerDecade
    use :: Table_Labels            , only : extrapolationTypeAbort
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout)               :: self
    double precision                                                  , intent(in   )               :: x
    double precision                                                  , intent(in   ), optional     :: xMinimumIn                    , xMaximumIn
    ! The tabulation is made at twenty points per decade. Ten points per decade leaves an interpolation error of order two parts
    ! in a thousand in the subresolution accretion rate, which is the accuracy to which that rate is checked against a direct
    ! evaluation of the hypergeometric functions - so the tabulated result sat exactly at the limit of its own tolerance, and any
    ! change in where the tabulated points fall relative to the value being interpolated could carry it over. The error falls as
    ! the square of the spacing, so this places it comfortably within that tolerance instead.
    integer                                                           , parameter                   :: xCountPerDecade=20
    double precision                                                  , parameter                   :: sqrtTwoOverPi  =sqrt(2.0d0/Pi)
    double precision                                                  , parameter                   :: xSeedMinimum   =1.0d-9        , xSeedMaximum=12.5d0
    double precision                                                                                :: xMinimum                      , xMaximum
    type            (rangeLattice                                    )                              :: latticeX
    logical                                                           , allocatable  , dimension(:) :: isComputed
    integer                                                                                         :: i
    logical                                                                                         :: tabulate

    ! Discard any previous tabulation which has been marked as invalid, so that no stale value is carried over into the new one.
    if (.not.self%subresolutionHypergeometricInitialized) then
       call self%subresolutionHypergeometric%destroy()
       self%latticeSubresolution=rangeLattice()
    end if
    ! Find the range to tabulate, pinning it to an absolute lattice so that the abscissae evaluated - and therefore every value
    ! interpolated between them - depend only on which lattice points are spanned, and not on the sequence of requests. The
    ! safety margin is one-sided - a factor of two above the request and none below it - so it is applied through the target
    ! rather than through `marginFactor`. The seed range is folded into the target, which is safe as it is a constant; the range
    ! already tabulated is unioned in through `latticeCurrent`, and never through the target, which would apply the margin to an
    ! already margined bound and so ratchet the range outward on every retabulation. The bounds are anchored to half decades
    ! rather than to whole decades because the tabulated function approaches the singular point of the hypergeometric function
    ! as its abscissa approaches zero, so the range should not be carried further into that limit than is needed.
    if (present(xMinimumIn)) then
       xMinimum=    xMinimumIn
    else
       xMinimum=min(xSeedMinimum,      (x-1.0d0))
    end if
    if (present(xMaximumIn)) then
       xMaximum=    xMaximumIn
    else
       xMaximum=max(xSeedMaximum,2.0d0*(x-1.0d0))
    end if
    latticeX=Range_Pinned(                                           &
         &                              [xMinimum,xMaximum]        , &
         &                               xCountPerDecade           , &
         &                               gridSchemePerDecade       , &
         &                marginFactor  = 1.0d0                    , &
         &                anchorEvery   = xCountPerDecade/2        , &
         &                latticeCurrent= self%latticeSubresolution  &
         &               )
    ! Decide whether the tabulation must be built or extended. The decision is taken from the pinned lattice rather than from the
    ! bounds directly, so that the decision and the range cannot disagree.
    tabulate=.not.self%subresolutionHypergeometricInitialized
    if (.not.tabulate) tabulate=.not.self%latticeSubresolution%covers(latticeX)
    if (tabulate) then
       ! Extend the tabulation onto the new lattice, preserving the values already computed - each depends only on its own
       ! abscissa and on γ₁, which is fixed, so a value carried over is precisely the value which would be computed afresh.
       call self%subresolutionHypergeometric%extend(latticeX,isComputed,tableCount=1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       do i=1,latticeX%count
          if (isComputed(i)) cycle
          call self%subresolutionHypergeometric%populate(                                                                               &
               &                                         +sqrtTwoOverPi                                                                 &
               &                                         *parkinsonColeHellyHypergeometricFactor(                                       &
               &                                                                                 self%subresolutionHypergeometric%x(i), &
               &                                                                                 self%gamma1                          , &
               &                                                                                 self%precisionHypergeometric           &
               &                                                                                )                                     , &
               &                                         i                                                                              &
               &                                        )
       end do
       self%latticeSubresolution                  =latticeX
       self%subresolutionHypergeometricInitialized=.true.
    end if
    return
  end subroutine parkinsonColeHellySubresolutionHypergeometricTabulate

  subroutine parkinsonColeHellyUpperBoundHypergeometricTabulate(self,mass,massResolution,massMinimumIn,massMaximumIn)
    !!{RST
    Tabulate the hypergeometric term appearing in the upper bound branching probability rate expression.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Range_Pinned          , gridSchemePerDecade
    use :: Table_Labels            , only : extrapolationTypeAbort
    use :: Display                 , only : displayGreen          , displayBlue, displayYellow, displayReset
    use :: Error                   , only : Error_Report
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout)               :: self
    double precision                                                  , intent(in   )               :: mass                              , massResolution
    double precision                                                  , intent(in   ), optional     :: massMinimumIn                     , massMaximumIn
    integer                                                           , parameter                   :: massCountPerDecade =30
    double precision                                                  , parameter                   :: sqrtTwoOverPi      =sqrt(2.0d0/Pi)
    double precision                                                  , parameter                   :: massSeedMaximum    =1.0d16
    double precision                                                                                :: massMinimum                       , massMaximum        , &
         &                                                                                             massScale                         , massTabulated
    type            (rangeLattice                                    )                              :: latticeMass
    logical                                                           , allocatable  , dimension(:) :: isComputed
    integer                                                                                         :: i
    logical                                                                                         :: tabulate
    double precision                                                                                :: massSigma                         , gammaEffective     , &
         &                                                                                             halfMassSigma                     , halfMassAlpha      , &
         &                                                                                             resolutionMassSigma               , resolutionMassAlpha

    ! Discard any previous tabulation which has been marked as invalid - the mass resolution on which it was built may have
    ! changed, or the growth of the mass variance may be mass dependent, in which case no value in it may be carried over.
    if (.not.self%upperBoundHypergeometricInitialized) then
       call self%upperBoundHypergeometric%destroy()
       self%latticeUpperBound      =rangeLattice()
       self%massResolutionTabulated=massResolution
    end if
    ! Find the range of masses to tabulate. Masses are tabulated *in units of twice the mass resolution*, and not as masses
    ! themselves. The lower end of this tabulation is not free: at twice the mass resolution the arguments of both hypergeometric
    ! functions below reach unity, beyond which they do not converge, so the range must begin exactly there. No lattice point
    ! coincides with that mass, but in units of it the lower bound is exactly the lattice point of index zero - which is what
    ! allows the axis to be pinned while still beginning precisely where it must. The mass resolution is fixed for the lifetime
    ! of the tabulation, being what invalidates it above, so these units do not shift under it.
    massScale  =2.0d0*self%massResolutionTabulated
    if (present(massMinimumIn)) then
       massMinimum=    massMinimumIn  /massScale
    else
       massMinimum=    1.0d0
    end if
    if (present(massMaximumIn)) then
       massMaximum=    massMaximumIn  /massScale
    else
       massMaximum=max(massSeedMaximum,2.0d0*mass)/massScale
    end if
    latticeMass=Range_Pinned(                                         &
         &                                 [massMinimum,massMaximum], &
         &                                  massCountPerDecade      , &
         &                                  gridSchemePerDecade     , &
         &                   marginFactor  = 1.0d0                  , &
         &                   limitMinimum  = 1.0d0                  , &
         &                   anchorEvery   = massCountPerDecade/2   , &
         &                   latticeCurrent= self%latticeUpperBound   &
         &                  )
    ! Decide whether the tabulation must be built or extended. The decision is taken from the pinned lattice rather than from the
    ! bounds directly, so that the decision and the range cannot disagree.
    tabulate=.not.self%upperBoundHypergeometricInitialized
    if (.not.tabulate) tabulate=.not.self%latticeUpperBound%covers(latticeMass)
    if (tabulate) then
       ! Extend the tabulation onto the new lattice, preserving the values already computed. Each depends only on its own mass
       ! and on the mass resolution - the ratios of σ which appear below are independent of epoch unless the growth of the mass
       ! variance is mass dependent, in which case the tabulation is invalidated on every use above - so a value carried over is
       ! precisely the value which would be computed afresh.
       call self%upperBoundHypergeometric%extend(latticeMass,isComputed,tableCount=1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       ! Evaluate σ and α at the mass resolution.
       call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,self%timeParent,resolutionMassSigma,resolutionMassAlpha)
       do i=1,latticeMass%count
          if (isComputed(i)) cycle
          ! Evaluate σ and α.
          massTabulated =massScale*self%upperBoundHypergeometric%x(i)
          call           self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*massTabulated,self%timeParent,halfMassSigma,halfMassAlpha)
          massSigma     =self%cosmologicalMassVariance_%rootVariance                      (      massTabulated,self%timeParent                            )
          if (abs(halfMassAlpha) <= alphaMinimum)                                                                                                                                                                                                                      &
               & call Error_Report(                                                                                                                                                                                                                                    &
               &                   'CDM assumptions were requested, but seem to be violated: unexpected |α| = |dlogσ/dlogM| ≪ 1'//char(10)                                                                                                                         // &
               &                   displayGreen()//'HELP:'//displayReset()//' set <'//displayBlue()//'cdmAssumptions'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"false"'//displayReset()//'/> '                         // &
               &                   'in class <'//displayBlue()//'mergerTreeBranchingProbability'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"parkinsonColeHelly"'//displayReset()//'/> to avoid making these assumptions'// &
               &                   {introspection:location}                                                                                                                                                                                                            &
               &                  )
          gammaEffective=self%gamma1-1.0d0/halfMassAlpha
          call self%upperBoundHypergeometric%populate(                                                                             &
               &                                      +sqrtTwoOverPi                                                               &
               &                                      *resolutionMassSigma                                                         &
               &                                      /massSigma                                                                   &
               &                                      *(                                                                           &
               &                                        +(halfMassSigma/massSigma)**(+gammaEffective-1.0d0)                        &
               &                                        /                           (-gammaEffective+1.0d0)                        &
               &                                        *Hypergeometric_2F1(                                                       &
               &                                                                             self%hypergeometricA(gammaEffective), &
               &                                                                             [      1.5d0-0.5d0*gammaEffective]  , &
               &                                                                             (massSigma/halfMassSigma)**2        , &
               &                                                           toleranceRelative=self%precisionHypergeometric          &
               &                                                          )                                                        &
               &                                        -(resolutionMassSigma/massSigma)**(+gammaEffective-1.0d0)                  &
               &                                        /                                 (-gammaEffective+1.0d0)                  &
               &                                        *Hypergeometric_2F1(                                                       &
               &                                                                             self%hypergeometricA(gammaEffective), &
               &                                                                             [      1.5d0-0.5d0*gammaEffective]  , &
               &                                                                             (massSigma/resolutionMassSigma)**2  , &
               &                                                           toleranceRelative=self%precisionHypergeometric          &
               &                                                          )                                                        &
               &                                       )                                                                         , &
               &                                      i                                                                            &
               &                                     )
       end do
       self%latticeUpperBound                  =latticeMass
       self%upperBoundHypergeometricInitialized=.true.
    end if
    return
  end subroutine parkinsonColeHellyUpperBoundHypergeometricTabulate

