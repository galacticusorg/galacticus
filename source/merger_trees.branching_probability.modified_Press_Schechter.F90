!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of branching probabilties in modified Press-Schechter theory.

module Modified_Press_Schechter_Branching
  !% Implements calculations of branching probabilties in modified Press-Schechter theory.
  use Cosmological_Mass_Variance
  use Numerical_Constants_Math
  use Tables
  implicit none
  private
  public :: Modified_Press_Schechter_Branching_Initialize   , Modified_Press_Schechter_Branching_State_Store, &
       &    Modified_Press_Schechter_Branching_State_Restore
  
  ! Parent halo shared variables.
  double precision            :: branchingProbabilityPreFactor                               , modificationG0Gamma2Factor  , &
       &                         parentDelta                                                 , parentHaloMass              , &
       &                         parentSigma                                                 , parentSigmaSquared          , &
       &                         probabilityMinimumMass                                      , probabilitySeek             , &
       &                         probabilityMaximumMassLog                                   , probabilityMinimumMassLog   , &
       &                         probabilityMaximum                                          , probabilityGradientMinimum  , &
       &                         probabilityGradientMaximum
  !$omp threadprivate(parentHaloMass,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass)
  !$omp threadprivate(modificationG0Gamma2Factor,branchingProbabilityPreFactor,probabilityMaximumMassLog,probabilityMinimumMassLog)
  !$omp threadprivate(probabilityMaximum,probabilityGradientMinimum,probabilityGradientMaximum)

  ! Parameters of the merger rate modification function.
  logical                     :: parametersRead                               =.false.
  double precision            :: modifiedPressSchechterG0                                    , modifiedPressSchechterGamma1, &
       &                         modifiedPressSchechterGamma2

  ! Accuracy parameter to ensure that merger rate function (which is correct to 1st order) is sufficiently accurate.
  double precision            :: modifiedPressSchechterFirstOrderAccuracy

  ! Tolerance parameter for hypergeometric evaluation.
  double precision            :: modifiedPressSchechterHypergeometricPrecision

  ! Tabulate hypergeometric factors?
  logical                     :: modifiedPressSchechterTabulateHypergeometricFactors
  
  ! Precomputed numerical factors.
  double precision, parameter :: sqrtTwoOverPi                                =sqrt(2.0d0/Pi)

  ! Branching probability integrand integration tolerance.
  double precision, parameter :: branchingProbabilityIntegrandToleraceRelative=1.0d-3

  ! Limit on alpha for use in effective gamma parameters.
  double precision, parameter :: alphaMinimum                                 =5.0d-3

  ! Branching assumptions.
  logical                     :: modifiedPressSchechterUseCDMAssumptions

  ! Warning records.
  logical                     :: hypergeometricFailureWarned=.false.

  ! Tables of hypergeometric factors.
  type   (table1DLogarithmicLinear) :: subresolutionHypergeometric                   , upperBoundHypergeometric
  logical                           :: subresolutionHypergeometricInitialized=.false., upperBoundHypergeometricInitialized=.false.
  double precision                  :: massResolutionTabulated               =-1.0d0
  !$omp threadprivate(subresolutionHypergeometric,upperBoundHypergeometric,subresolutionHypergeometricInitialized,upperBoundHypergeometricInitialized,massResolutionTabulated)

contains

  !# <treeBranchingMethod>
  !#  <unitName>Modified_Press_Schechter_Branching_Initialize</unitName>
  !# </treeBranchingMethod>
  subroutine Modified_Press_Schechter_Branching_Initialize(treeBranchingMethod,Tree_Branching_Probability_Bound,Tree_Branching_Probability&
       &,Tree_Subresolution_Fraction,Tree_Branch_Mass,Tree_Maximum_Step)
    !% Initialize the modified Press-Schechter branching routines.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                      ), intent(in   )          :: treeBranchingMethod
    procedure(Modified_Press_Schechter_Branching_Probability_Bound), intent(inout), pointer :: Tree_Branching_Probability_Bound
    procedure(Modified_Press_Schechter_Branching_Probability      ), intent(inout), pointer :: Tree_Branching_Probability
    procedure(Modified_Press_Schechter_Subresolution_Fraction     ), intent(inout), pointer :: Tree_Subresolution_Fraction
    procedure(Modified_Press_Schechter_Branch_Mass_Generic        ), intent(inout), pointer :: Tree_Branch_Mass
    procedure(Modified_Press_Schechter_Branching_Maximum_Step     ), intent(inout), pointer :: Tree_Maximum_Step
    
    if (treeBranchingMethod == 'modifiedPress-Schechter') then
       call Modified_Press_Schechter_Branching_Parameters()
       Tree_Branching_Probability_Bound  => Modified_Press_Schechter_Branching_Probability_Bound
       Tree_Branching_Probability        => Modified_Press_Schechter_Branching_Probability
       Tree_Subresolution_Fraction       => Modified_Press_Schechter_Subresolution_Fraction
       Tree_Maximum_Step                 => Modified_Press_Schechter_Branching_Maximum_Step
       if (modifiedPressSchechterUseCDMAssumptions) then
          Tree_Branch_Mass               => Modified_Press_Schechter_Branch_Mass_CDMAssumptions
       else
          Tree_Branch_Mass               => Modified_Press_Schechter_Branch_Mass_Generic
       end if
    end if
    return
  end subroutine Modified_Press_Schechter_Branching_Initialize

  subroutine Modified_Press_Schechter_Branching_Parameters()
    !% Read paramters for the modified Press-Schechter merger tree branching module.
    use Input_Parameters
    implicit none
    
    if (.not.parametersRead) then
       !$omp critical(Modified_Press_Schechter_Branching_Parameters)
       if (.not.parametersRead) then
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterG0</name>
          !@   <defaultValue>0.57</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter $G_0$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterG0'                ,modifiedPressSchechterG0                ,defaultValue=&
               & 0.57d0)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterGamma1</name>
          !@   <defaultValue>0.38</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter $\gamma_1$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterGamma1'            ,modifiedPressSchechterGamma1            ,defaultValue=&
               & 0.38d0)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterGamma2</name>
          !@   <defaultValue>-0.01</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter $\gamma_2$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterGamma2'            ,modifiedPressSchechterGamma2            ,defaultValue=&
               &-0.01d0)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterFirstOrderAccuracy</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Limits the step in $\delta_{\mathrm crit}$ when constructing merger trees using the \cite{parkinson_generating_2008}
          !@     algorithm, so that it never exceeds {\normalfont \ttfamily
          !@     modifiedPressSchechterFirstOrderAccuracy}$\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterFirstOrderAccuracy',modifiedPressSchechterFirstOrderAccuracy,defaultValue &
               &=0.1d0)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterHypergeometricPrecision</name>
          !@   <defaultValue>$10^{-6}$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The fractional precision required in evaluates of hypergeometric functions in the modified Press-Schechter tree branching calculations.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterHypergeometricPrecision',modifiedPressSchechterHypergeometricPrecision,defaultValue &
               &=1.0d-6)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterTabulateHypergeometricFactors</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether hypergeometric factors should be precomputed and tabulated in modified Press-Schechter tree branching functions.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterTabulateHypergeometricFactors',modifiedPressSchechterTabulateHypergeometricFactors,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>modifiedPressSchechterUseCDMAssumptions</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     If true, assume that $\alpha(=-{\mathrm d}\log \sigma/{\mathrm d}\log M)>0$ and ${\mathrm d}\alpha/{\mathrm d}M>0$ (as is true in the case of \gls{cdm}) when constructing merger trees using the \cite{parkinson_generating_2008}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('modifiedPressSchechterUseCDMAssumptions',modifiedPressSchechterUseCDMAssumptions,defaultValue=.false.)
          parametersRead=.true.
       end if
       !$omp end critical(Modified_Press_Schechter_Branching_Parameters)
    end if
    return
  end subroutine Modified_Press_Schechter_Branching_Parameters
  
  double precision function Modified_Press_Schechter_Branch_Mass_CDMAssumptions(haloMass,deltaCritical,massResolution,probability,randomNumberGenerator)
    !% A merger tree branch split mass function which assumes a \gls{cdm}-like power
    !% spectrum. With these assumptions, it can employ the mass sampling algorithm of
    !% \cite{parkinson_generating_2008}. One difference with respect to the algorithm of
    !% \cite{parkinson_generating_2008} is that here the normalization of their function $S(q)$
    !% (eqn. A2) is irrelevant, since a branch split has already been decided to have
    !% occcurred---all that remains necessary is to determine its mass. Variable and function
    !% names follow \cite{parkinson_generating_2008}.
    use Pseudo_Random
    implicit none
    double precision                               , intent(in   ) :: deltaCritical                 , haloMass       , &
         &                                                            massResolution                , probability
    type            (pseudoRandom                 ), intent(inout) :: randomNumberGenerator
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                                               :: massFractionResolution        , beta           , &
         &                                                            B                             , halfMassSigma  , &
         &                                                            halfMassAlpha                 , eta            , &
         &                                                            x                             , mu             , &
         &                                                            massFraction                  , resolutionSigma, &
         &                                                            massFractionResolutionPowerEta, halfPowerEta   , &
         &                                                            halfMassV
    logical                                                        :: reject

    ! Get parent and half-mass sigmas and alphas.
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    parentSigmaSquared=cosmologicalMassVariance_%rootVariance(haloMass)**2
    call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*haloMass,halfMassSigma,halfMassAlpha)
    ! Compute parameters beta, mu, and B.
    massFractionResolution=+massResolution                 &
         &                 /haloMass
    halfMassV             =+V(0.5d0)
    beta                  =+log(                           &
         &                      +V(massFractionResolution) &
         &                      /halfMassV                 &
         &                     )                           &
         &                 /log(                           &
         &                      +massFractionResolution    &
         &                      /0.5d0                     &
         &                     )
    B                     =+halfMassV                      &
         &                 *2.0d0    **beta
    if (modifiedPressSchechterGamma1 >= 0.0d0) then
       mu                 =-halfMassAlpha
    else
       resolutionSigma    =+cosmologicalMassVariance_%rootVariance(massResolution)
       mu                 =-log(                        &
            &                   +resolutionSigma        &
            &                   /halfMassSigma          &
            &                  )                        &
            &              /log(                        &
            &                   +massFractionResolution &
            &                   /0.5d0                  &
            &                  )
    end if
    eta                           =+beta                              &
         &                         -1.0d0                             &
         &                         -mu                                &
         &                         *modifiedPressSchechterGamma1
    massFractionResolutionPowerEta=+massFractionResolution      **eta
    halfPowerEta                  =+0.5d0                       **eta
    ! Sample from S(q), using rejection sampling on R(q) to decide whether to keep/reject the
    ! proposed q.
    reject=.true.
    do while (reject)
       ! Draw a random q from S(q).
       x           =randomNumberGenerator%sample()
       massFraction=(                                                 &
            &        +              massFractionResolutionPowerEta    &
            &        +(halfPowerEta-massFractionResolutionPowerEta)*x &
            &       )**(1.0d0/eta)
       x           =randomNumberGenerator%sample()
       reject=x > R(massFraction)
    end do
    Modified_Press_Schechter_Branch_Mass_CDMAssumptions=massFraction*haloMass
    return

  contains

    double precision function R(massFraction)
      !% The function $R(q)$ from \cite[][eqn. A3]{parkinson_generating_2008}.
      implicit none
      double precision, intent(in   ) :: massFraction
      double precision                :: massFractionSigma, massFractionAlpha

       call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massFraction*haloMass,massFractionSigma,massFractionAlpha)   
       R      =+(                               &
            &    +massFractionAlpha             &
            &    /halfMassAlpha                 &
            &   )                               &
            &  *V(massFraction)                 &
            &  /B                               &
            &  /massFraction**beta              &
            &  *(                               &
            &    +(                             &
            &      +2.0d0                       &
            &      *massFraction                &
            &     )**mu                         &
            &    *massFractionSigma             &
            &    /halfMassSigma                 &
            &   )**modifiedPressSchechterGamma1
       return
    end function R
      
    double precision function V(massFraction)
      !% The function $V(q)$ from \cite[][eqn. A4]{parkinson_generating_2008}.
      implicit none
      double precision, intent(in   ) :: massFraction
      double precision                :: childSigmaSquared

      childSigmaSquared=cosmologicalMassVariance_%rootVariance(massFraction*haloMass)**2
      V                =+   childSigmaSquared &
           &            /(                    &
           &              + childSigmaSquared &
           &              -parentSigmaSquared &
           &             )**1.5d0
      return
    end function V
    
  end function Modified_Press_Schechter_Branch_Mass_CDMAssumptions
  
  double precision function Modified_Press_Schechter_Branch_Mass_Generic(haloMass,deltaCritical,massResolution,probability,randomNumberGenerator)
    !% Determine the mass of one of the halos to which the given halo branches, given the branching probability,
    !% {\normalfont \ttfamily probability}. Typically, {\normalfont \ttfamily probabilityFraction} is found by multiplying {\tt
    !% Modified\_Press\_Schechter\_Branching\_Probability()} by a random variable drawn in the interval 0--1 if a halo
    !% branches. This routine then finds the progenitor mass corresponding to this value.
    use Pseudo_Random
    use Root_Finder
    implicit none
    double precision                               , intent(in   ) :: deltaCritical              , haloMass                , &
         &                                                            massResolution             , probability
    type            (pseudoRandom                 ), intent(inout) :: randomNumberGenerator
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                               , parameter     :: toleranceAbsolute    =0.0d0, toleranceRelative=1.0d-9
    type            (rootFinder                   ), save          :: finder
    !$omp threadprivate(finder)
    double precision                                               :: logMassMinimum             , logMassMaximum
    
    ! Get required objects.
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Initialize global variables.
    parentHaloMass           =haloMass
    parentSigma              =cosmologicalMassVariance_%rootVariance(haloMass)
    parentDelta              =deltaCritical
    probabilityMinimumMass   =massResolution
    probabilityMinimumMassLog=log(massResolution)
    probabilityMaximumMassLog=log(0.5d0*haloMass)
    probabilitySeek          =probability
    call Compute_Common_Factors
    ! Check the sign of the root function at half the halo mass.
    if (Modified_Press_Schechter_Branch_Mass_Root(probabilityMaximumMassLog) >= 0.0d0) then
       ! The root function is zero, or very close to it (which can happen due to rounding errors
       ! occasionally). Therefore we have an almost perfect binary split.
       Modified_Press_Schechter_Branch_Mass_Generic=0.5d0*haloMass
    else
       ! Initialize our root finder.
       if (.not.finder%isInitialized()) then
          call finder%rootFunctionDerivative(                                                      &
               &                             Modified_Press_Schechter_Branch_Mass_Root           , &
               &                             Modified_Press_Schechter_Branch_Mass_Root_Derivative, &
               &                             Modified_Press_Schechter_Branch_Mass_Root_Both        &
               &                            )
          call finder%tolerance             (                                                      &
               &                             toleranceAbsolute                                   , &
               &                             toleranceRelative                                     &
               &                            )
          call finder%typeDerivative        (                                                      &
               &                             FGSL_Root_fdfSolver_Steffenson                        &
               &                            )
       end if
       ! Split is not binary - seek the actual mass of the smaller progenitor.
       logMassMinimum                              =log(      massResolution)
       logMassMaximum                              =log(0.5d0*haloMass      )
       probabilityGradientMinimum                  =Modified_Press_Schechter_Branch_Mass_Root_Derivative(logMassMinimum)
       probabilityGradientMaximum                  =Modified_Press_Schechter_Branch_Mass_Root_Derivative(logMassMaximum)
       probabilityMaximum                          =Modified_Press_Schechter_Branch_Mass_Root           (logMassMaximum)
       Modified_Press_Schechter_Branch_Mass_Generic=exp(finder%find(rootRange=[logMassMinimum,logMassMaximum]))
    end if
    return
  end function Modified_Press_Schechter_Branch_Mass_Generic

  double precision function Modified_Press_Schechter_Branch_Mass_Root(logMassMaximum)
    !% Used to find the mass of a merger tree branching event.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: logMassMaximum
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    type            (c_ptr                     )                :: parameterPointer
    double precision                                            :: integral            , massMaximum

    if (logMassMaximum < probabilityMinimumMassLog) then
       Modified_Press_Schechter_Branch_Mass_Root=probabilitySeek   +probabilityGradientMinimum*(logMassMaximum-probabilityMinimumMassLog)
    else if (logMassMaximum > probabilityMaximumMassLog) then
       Modified_Press_Schechter_Branch_Mass_Root=probabilityMaximum+probabilityGradientMaximum*(logMassMaximum-probabilityMaximumMassLog)
    else
       massMaximum=+exp(logMassMaximum)
       integral   =+branchingProbabilityPreFactor                                                                        &
            &      *Integrate(                                                                                           &
            &                 probabilityMinimumMassLog                                                                , &
            &                 logMassMaximum                                                                           , &
            &                 Branching_Probability_Integrand_Logarithmic                                              , &
            &                 parameterPointer                                                                         , &
            &                 integrandFunction                                                                        , &
            &                 integrationWorkspace                                                                     , &
            &                 toleranceAbsolute                          =0.0d0                                        , &
            &                 toleranceRelative                          =branchingProbabilityIntegrandToleraceRelative, &
            &                 integrationRule                            =FGSL_Integ_Gauss15                             &
            &                )
       call Integrate_Done(integrandFunction,integrationWorkspace)
       Modified_Press_Schechter_Branch_Mass_Root=probabilitySeek-integral
    end if
    return
  end function Modified_Press_Schechter_Branch_Mass_Root

  double precision function Modified_Press_Schechter_Branch_Mass_Root_Derivative(logMassMaximum)
    !% Used to find the mass of a merger tree branching event.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Galacticus_Error
    implicit none
    double precision       , intent(in   ) :: logMassMaximum
    double precision                       :: integral        , massMaximum
    type            (c_ptr)                :: parameterPointer

    if (logMassMaximum < probabilityMinimumMassLog) then
       Modified_Press_Schechter_Branch_Mass_Root_Derivative=probabilityGradientMinimum
    else if (logMassMaximum > probabilityMaximumMassLog) then
       Modified_Press_Schechter_Branch_Mass_Root_Derivative=probabilityGradientMaximum
    else
       massMaximum=+exp(logMassMaximum)
       integral   =+branchingProbabilityPreFactor                                              &
            &      *Branching_Probability_Integrand_Logarithmic(                               &
            &                                                   max(                           &
            &                                                       log(massMaximum)         , &
            &                                                       probabilityMinimumMassLog  &
            &                                                      )                         , &
            &                                                   parameterPointer               &
            &                                                  )
       Modified_Press_Schechter_Branch_Mass_Root_Derivative=-integral
    end if
    return
  end function Modified_Press_Schechter_Branch_Mass_Root_Derivative

  subroutine Modified_Press_Schechter_Branch_Mass_Root_Both(massMaximum,massRoot,massRootDerivative)
    !% Used to find the mass of a merger tree branching event.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision, intent(in   ) :: massMaximum
    double precision, intent(  out) :: massRoot   , massRootDerivative
    
    massRoot          =Modified_Press_Schechter_Branch_Mass_Root           (massMaximum)
    massRootDerivative=Modified_Press_Schechter_Branch_Mass_Root_Derivative(massMaximum)
    return
  end subroutine Modified_Press_Schechter_Branch_Mass_Root_Both

  double precision function Modified_Press_Schechter_Branching_Maximum_Step(haloMass,deltaCritical,massResolution)
    !% Return the maximum allowed step in $\delta_{\mathrm crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time {\tt
    !% deltaCritical} should be allowed to take.
    implicit none
    double precision                               , intent(in   ) :: deltaCritical                   , haloMass, &
         &                                                            massResolution
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                               , parameter     :: largeStep                =1.0d10           !   Effectively infinitely large step in w(=delta_crit).
    double precision                                               :: parentHalfMassSigma

    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       cosmologicalMassVariance_ => cosmologicalMassVariance()
       parentSigma        =cosmologicalMassVariance_%rootVariance(      haloMass)
       parentHalfMassSigma=cosmologicalMassVariance_%rootVariance(0.5d0*haloMass)
       Modified_Press_Schechter_Branching_Maximum_Step=modifiedPressSchechterFirstOrderAccuracy*sqrt(2.0d0&
            &*(parentHalfMassSigma**2-parentSigma**2))
    else
       Modified_Press_Schechter_Branching_Maximum_Step=largeStep
    end if
    return
  end function Modified_Press_Schechter_Branching_Maximum_Step

  double precision function Modified_Press_Schechter_Branching_Probability(haloMass,deltaCritical,massResolution)
    !% Return the probability per unit change in $\delta_{\mathrm crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time {\tt
    !% deltaCritical} will undergo a branching to progenitors with mass greater than {\normalfont \ttfamily massResolution}.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                               , intent(in   ) :: deltaCritical                , haloMass                    , &
         &                                                            massResolution
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    type            (c_ptr                        )                :: parameterPointer
    type            (fgsl_function                )                :: integrandFunction
    type            (fgsl_integration_workspace   )                :: integrationWorkspace
    double precision                                               :: massMaximum                  , massMinimum
    double precision                               , save          :: haloMassPrevious      =-1.0d0, deltaCriticalPrevious=-1.0d0, &
         &                                                            massResolutionPrevious=-1.0d0, probabilityPrevious
    !$omp threadprivate(haloMassPrevious,deltaCriticalPrevious,massResolutionPrevious,probabilityPrevious)
    
    ! Recompute branching probability if necessary.
    if (haloMass /= haloMassPrevious .or. deltaCritical /= deltaCriticalPrevious .or. massResolution /= massResolutionPrevious) then
       haloMassPrevious=haloMass
       deltaCriticalPrevious=deltaCritical
       massResolutionPrevious=massResolution
       ! Get sigma and delta_critical for the parent halo.
       if (haloMass>2.0d0*massResolution) then
          cosmologicalMassVariance_ => cosmologicalMassVariance()
          parentHaloMass=haloMass
          parentSigma=cosmologicalMassVariance_%rootVariance(haloMass)
          parentDelta=deltaCritical
          call Compute_Common_Factors()
          massMinimum=      massResolution
          massMaximum=0.5d0*parentHaloMass
          probabilityPrevious=+branchingProbabilityPreFactor                                                                        &
               &              *Integrate(                                                                                           &
               &                         log(massMinimum)                                                                         , &
               &                         log(massMaximum)                                                                         , &
               &                         Branching_Probability_Integrand_Logarithmic                                              , &
               &                         parameterPointer                                                                         , &
               &                         integrandFunction                                                                        , &
               &                         integrationWorkspace                                                                     , &
               &                         toleranceAbsolute                          =0.0d0                                        , &
               &                         toleranceRelative                          =branchingProbabilityIntegrandToleraceRelative, &
               &                         integrationRule                            =FGSL_Integ_Gauss15                             &
               &                        )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       else
          probabilityPrevious=0.0d0
       end if
    end if
    Modified_Press_Schechter_Branching_Probability=probabilityPrevious
    return
  end function Modified_Press_Schechter_Branching_Probability

  double precision function Modified_Press_Schechter_Branching_Probability_Bound(haloMass,deltaCritical,massResolution,bound)
    !% Return the probability per unit change in $\delta_{\mathrm crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time {\tt
    !% deltaCritical} will undergo a branching to progenitors with mass greater than {\normalfont \ttfamily massResolution}.
    use Merger_Tree_Branching_Options
    use Galacticus_Error
    use Galacticus_Display
    use Hypergeometric_Functions
    use FGSL
    use Numerical_Comparison
    implicit none
    double precision                               , intent(in   ) :: deltaCritical                   , haloMass                 , &
         &                                                            massResolution
    integer                                        , intent(in   ) :: bound
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                               , save          :: resolutionSigma                 , resolutionAlpha
    !$omp threadprivate(resolutionSigma,resolutionAlpha)
    double precision                                               :: probabilityIntegrandLower       , probabilityIntegrandUpper, &
         &                                                            halfParentSigma                 , halfParentAlpha          , &
         &                                                            gammaEffective
    double precision                                               :: hyperGeometricFactorLower       , hyperGeometricFactorUpper, &
         &                                                            resolutionSigmaOverParentSigma
    integer         (fgsl_int                     )                :: statusLower                     , statusUpper
    logical                                                        :: usingCDMAssumptions
    integer                                                        :: iBound
    
    ! Get sigma and delta_critical for the parent halo.
    if (haloMass > 2.0d0*massResolution) then
       cosmologicalMassVariance_ => cosmologicalMassVariance()
       parentHaloMass=haloMass
       parentSigma=cosmologicalMassVariance_%rootVariance(haloMass)
       parentDelta=deltaCritical
       call Compute_Common_Factors()
       if (massResolution /= massResolutionTabulated) then
          ! Resolution changed - recompute sigma and alpha at resolution limit. Also reset the hypergeometric factor tables since
          ! these depend on resolution.
          call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,resolutionSigma,resolutionAlpha)
          upperBoundHypergeometricInitialized=.false.
       end if
       resolutionSigmaOverParentSigma=resolutionSigma/parentSigma
       ! Estimate probability.
       if (resolutionSigmaOverParentSigma > 1.0d0) then
          ! Compute relevant sigmas and alphas.
          call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*parentHaloMass,halfParentSigma,halfParentAlpha)
          ! Iterative over available bounds.
          do iBound=1,2
             ! Determine if CDM assumptions can be used. Do this only is these have been explicitly allowed, if this is our first
             ! pass through the bounds evaluation, and if both alphas are sufficiently large. (This last condition is required
             ! since we raise quantities to the power of 1/alpha which can cause problems for very small alpha.)
             usingCDMAssumptions= modifiedPressSchechterUseCDMAssumptions &
                  &              .and.                                    &
                  &               iBound               == 1               &
                  &              .and.                                    &
                  &               abs(resolutionAlpha) >  alphaMinimum    &
                  &              .and.                                    &
                  &               abs(halfParentAlpha) >  alphaMinimum
             ! Compute the effective value of gamma.             
             gammaEffective=modifiedPressSchechterGamma1
             if (usingCDMAssumptions) then
                select case (bound)
                case (boundLower)
                   gammaEffective=gammaEffective-1.0d0/resolutionAlpha
                case (boundUpper)
                   gammaEffective=gammaEffective-1.0d0/halfParentAlpha
                end select
             end if
             ! Compute probability factors. The logic here becomes complicated, as we use various optimizations and tabulations to
             ! speed up calculation.
             !
             ! Tabulations will only be used if modifiedPressSchechterTabulateHypergeometricFactors is true.
             !
             ! Set status to success by default.
             statusLower=FGSL_Success
             statusUpper=FGSL_Success
             ! First, check if CDM assumptions are not being used and we're allowed to tabulate hypergeometric factors, 
             if (.not.usingCDMAssumptions.and.modifiedPressSchechterTabulateHypergeometricFactors) then
                ! CDM assumptions are not being used. In this case we can use the same table of hypergeometric factors as the
                ! subresolution merger fraction.
                call Subresolution_Hypergeometric_Tabulate(resolutionSigmaOverParentSigma)
                call Subresolution_Hypergeometric_Tabulate(halfParentSigma   /parentSigma)
                probabilityIntegrandLower=+modificationG0Gamma2Factor*subresolutionHypergeometric%interpolate(+resolutionSigmaOverParentSigma-1.0d0)/parentSigma
                probabilityIntegrandUpper=+modificationG0Gamma2Factor*subresolutionHypergeometric%interpolate(+halfParentSigma   /parentSigma-1.0d0)/parentSigma                
             else
                ! Next, check if CDM assumptions are being used, we're allowed to tabulate hypergeometric factors, and the bound
                ! requested is the upper bound.                
                if     ( usingCDMAssumptions                                 &
                     &  .and.                                                &
                     &   modifiedPressSchechterTabulateHypergeometricFactors &
                     &  .and.                                                &
                     &   bound == boundUpper                                 &
                     & ) then
                   ! Use a tabulation of the hypergeometric functions for the upper bound, made using CDM assumptions. Since the
                   ! tables already include the difference between the upper and lower integrand, we simply set the lower
                   ! integrand to zero here.
                   call Upper_Bound_Hypergeometric_Tabulate(parentHaloMass,massResolution)
                   probabilityIntegrandUpper=modificationG0Gamma2Factor*upperBoundHypergeometric%interpolate(parentHaloMass)
                   probabilityIntegrandLower=0.0d0
                else
                   ! Use a direct calculation of the hypergeometric factors in this case.
                   hyperGeometricFactorLower=Hypergeometric_2F1(                                                                 &
                        &                                       [1.5d0,0.5d0-0.5d0*gammaEffective]                             , &
                        &                                       [      1.5d0-0.5d0*gammaEffective]                             , &
                        &                                       1.0d0/resolutionSigmaOverParentSigma**2                        , &
                        &                                       toleranceRelative=modifiedPressSchechterHypergeometricPrecision, &
                        &                                       status           =statusLower                                    &
                        &                                      )
                   if (statusLower /= FGSL_Success) then
                      if (usingCDMAssumptions) then
                         if (.not.hypergeometricFailureWarned) then
                            hypergeometricFailureWarned=.true.
                            call Galacticus_Display_Message(                                                                                &
                                 &                          'WARNING: hypergeometric function evaluation failed when computing'//char(10)// &
                                 &                          'merger tree branching probability bounds - will revert to more'   //char(10)// &
                                 &                          'robust (but less stringent) bound in this and future cases'                 ,  &
                                 &                          verbosityWarn                                                                   &
                                 &                         )
                         end if
                         cycle
                      else
                         call Galacticus_Error_Report('Modified_Press_Schechter_Branching_Probability_Bound','hypergeometric function evaluation failed')
                      end if
                   end if
                   probabilityIntegrandLower=+sqrtTwoOverPi                                              &
                        &                    *(modificationG0Gamma2Factor/parentSigma)                   &
                        &                    *(resolutionSigmaOverParentSigma**(gammaEffective-1.0d0))   &
                        &                    /(1.0d0-gammaEffective)                                     &
                        &                    *hyperGeometricFactorLower         
                   ! Check if we can use a table to compute the upper factor.
                   hyperGeometricFactorUpper=Hypergeometric_2F1(                                                                 &
                        &                                       [1.5d0,0.5d0-0.5d0*gammaEffective]                             , &
                        &                                       [      1.5d0-0.5d0*gammaEffective]                             , &
                        &                                       parentSigma**2/halfParentSigma**2                              , &
                        &                                       toleranceRelative=modifiedPressSchechterHypergeometricPrecision, &
                        &                                       status           =statusUpper                                    &
                        &                                      )
                   if (statusUpper /= FGSL_Success) then
                      if (usingCDMAssumptions) then
                         if (.not.hypergeometricFailureWarned) then
                            hypergeometricFailureWarned=.true.
                            call Galacticus_Display_Message(                                                                                &
                                 &                          'WARNING: hypergeometric function evaluation failed when computing'//char(10)// &
                                 &                          'merger tree branching probability bounds - will revert to more'   //char(10)// &
                                 &                          'robust (but less stringent) bound in this and future cases'                 ,  &
                                 &                          verbosityWarn                                                                   &
                                 &                         )
                         end if
                         cycle
                      else
                         call Galacticus_Error_Report('Modified_Press_Schechter_Branching_Probability_Bound','hypergeometric function evaluation failed')
                      end if
                   end if
                   probabilityIntegrandUpper=+sqrtTwoOverPi                                           &
                        &                    *(modificationG0Gamma2Factor/parentSigma)                &
                        &                    *((halfParentSigma/parentSigma)**(gammaEffective-1.0d0)) &
                        &                    /(1.0d0-gammaEffective)                                  &
                        &                    *hyperGeometricFactorUpper   
                end if
             end if
             ! Compute the bound.
             select case (bound)
             case (boundLower)
                if (usingCDMAssumptions) then
                   Modified_Press_Schechter_Branching_Probability_Bound=+(                           &
                        &                                                 +probabilityIntegrandUpper &
                        &                                                 -probabilityIntegrandLower &
                        &                                                 )                          &
                        &                                                *parentHaloMass             &
                        &                                                /massResolution             &
                        &                                                *(                          &
                        &                                                  +resolutionSigma          &
                        &                                                  /parentSigma              &
                        &                                                 )**(1.0d0/resolutionAlpha)
                else
                   Modified_Press_Schechter_Branching_Probability_Bound=+(                           &
                        &                                                 +probabilityIntegrandUpper &
                        &                                                 -probabilityIntegrandLower &
                        &                                                 )                          &
                        &                                                *       parentHaloMass      &
                        &                                                /(0.5d0*parentHaloMass)
                end if
             case (boundUpper)
                if (usingCDMAssumptions) then
                   Modified_Press_Schechter_Branching_Probability_Bound=+(                           &
                        &                                                 +probabilityIntegrandUpper &
                        &                                                 -probabilityIntegrandLower &
                        &                                                 )                          &
                        &                                                *parentHaloMass             &
                        &                                                /massResolution             &
                        &                                                *(                          &
                        &                                                 +resolutionSigma           &
                        &                                                 /parentSigma               &
                        &                                                )**(1.0d0/halfParentAlpha)
                else
                   Modified_Press_Schechter_Branching_Probability_Bound=+(                           &
                        &                                                 +probabilityIntegrandUpper &
                        &                                                 -probabilityIntegrandLower &
                        &                                                 )                          &
                        &                                                *parentHaloMass             &
                        &                                                /massResolution
                end if
             case default
                call Galacticus_Error_Report('Modified_Press_Schechter_Branching_Probability_Bound','unknown bound type')
             end select
             if (statusUpper == FGSL_Success .and. statusLower == FGSL_Success) exit
          end do
       else
          Modified_Press_Schechter_Branching_Probability_Bound=-1.0d0
       end if
    else
       Modified_Press_Schechter_Branching_Probability_Bound=0.0d0
    end if
    return
  end function Modified_Press_Schechter_Branching_Probability_Bound

  double precision function Modified_Press_Schechter_Subresolution_Fraction(haloMass,deltaCritical,massResolution)
    !% Return the fraction of mass accreted in subresolution halos, i.e. those below {\normalfont \ttfamily massResolution}, per unit change in
    !% $\delta_{\mathrm crit}$ for a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily deltaCritical}. The integral is computed analytically in
    !% terms of the $_2F_1$ hypergeometric function.
    use Hypergeometric_Functions
    implicit none
    double precision                               , intent(in   ) :: deltaCritical                 , haloMass                      , &
         &                                                            massResolution
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                               , save          :: massResolutionPrevious=-1.0d+0, resolutionSigma
    !$omp threadprivate(resolutionSigma,massResolutionPrevious)
    double precision                                               :: hyperGeometricFactor          , resolutionSigmaOverParentSigma

    ! Get required objects.
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Get sigma and delta_critical for the parent halo.
    parentHaloMass=haloMass
    parentSigma   =cosmologicalMassVariance_%rootVariance(haloMass)
    parentDelta   =deltaCritical
    call Compute_Common_Factors
    if (massResolution /= massResolutionPrevious) then
       resolutionSigma       =cosmologicalMassVariance_%rootVariance(massResolution)
       massResolutionPrevious=                                       massResolution
    end if
    resolutionSigmaOverParentSigma=resolutionSigma/parentSigma
    if (resolutionSigmaOverParentSigma > 1.0d0) then
       if (modifiedPressSchechterTabulateHypergeometricFactors) then
          ! Use tabulation of hypergeometric factors.
          call Subresolution_Hypergeometric_Tabulate(resolutionSigmaOverParentSigma)
          Modified_Press_Schechter_Subresolution_Fraction=+modificationG0Gamma2Factor                                              &
               &                                          *subresolutionHypergeometric%interpolate(                                &
               &                                                                                   +resolutionSigmaOverParentSigma &
               &                                                                                   -1.0d0                          &
               &                                                                                  )                                &
               &                                          /parentSigma
       else
          ! Compute hypergeometric factors directly.
          hyperGeometricFactor=Hypergeometric_2F1(                                                                 &
               &                                  [1.5d0,0.5d0-0.5d0*modifiedPressSchechterGamma1]               , &
               &                                  [      1.5d0-0.5d0*modifiedPressSchechterGamma1]               , &
               &                                  1.0d0/resolutionSigmaOverParentSigma**2                        , &
               &                                  toleranceRelative=modifiedPressSchechterHypergeometricPrecision  &
               &)
          Modified_Press_Schechter_Subresolution_Fraction=+sqrtTwoOverPi                                                         &
               &                                          *modificationG0Gamma2Factor                                            &
               &                                          /parentSigma                                                           &
               &                                          *resolutionSigmaOverParentSigma**(+modifiedPressSchechterGamma1-1.0d0) &
               &                                          /                                (-modifiedPressSchechterGamma1+1.0d0) &
               &                                          *hyperGeometricFactor
       end if
    else
       Modified_Press_Schechter_Subresolution_Fraction=-1.0d0
    end if
    return
  end function Modified_Press_Schechter_Subresolution_Fraction

  function Branching_Probability_Integrand_Logarithmic(logChildHaloMass,parameterPointer) bind(c)
    !% Integrand for the branching probability.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real (kind=c_double                )          :: Branching_Probability_Integrand_Logarithmic
    real (kind=c_double                ), value   :: logChildHaloMass
    type (     c_ptr                   ), value   :: parameterPointer
    class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
    real (kind=c_double                )          :: childAlpha                     , childSigma, childHaloMass
    
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    childHaloMass=exp(logChildHaloMass)
    call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,childSigma,childAlpha)
    Branching_Probability_Integrand_Logarithmic=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)*childHaloMass
    return
  end function Branching_Probability_Integrand_Logarithmic

  double precision function Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    !% Progenitor mass function from Press-Schechter. The constant factor of the parent halo
    !% mass is not included here---instead it is included in a multiplicative prefactor by which
    !% integrals over this function are multiplied.
    implicit none
    double precision, intent(in   ) :: childAlpha, childHaloMass, childSigma

    Progenitor_Mass_Function=Merging_Rate(childSigma,childAlpha)*Modification_Function(childSigma)/(childHaloMass**2)
    return
  end function Progenitor_Mass_Function

  double precision function Merging_Rate(childSigma,childAlpha)
    !% Merging rate from Press-Schechter. The constant factor of sqrt(2/pi) not included
    !% here---instead it is included in a multiplicative prefactor by which integrals over this
    !% function are multiplied.
    implicit none
    double precision, intent(in   ) :: childAlpha       , childSigma
    double precision                :: childSigmaSquared

    childSigmaSquared=childSigma**2
    if (childSigmaSquared > parentSigmaSquared .and. childAlpha < 0.0d0) then
       Merging_Rate=(childSigmaSquared/((childSigmaSquared-parentSigmaSquared)**1.5d0))*abs(childAlpha)
    else
       Merging_Rate=0.0d0
    end if
    return
  end function Merging_Rate

  double precision function Modification_Function(childSigma)
    !% Empirical modification of the progenitor mass function from
    !% \cite{parkinson_generating_2008}. The constant factors of $G_0 (\delta_{\rm
    !% p}/\sigma_{\mathrm p})^{\gamma_2}$ and $1/\sigma_{\mathrm p}^{\gamma_1}$ are not included
    !% here---instead they are included in a multiplicative prefactor by which integrals over
    !% this function are multiplied.
    implicit none
    double precision, intent(in   ) :: childSigma

    Modification_Function=childSigma**modifiedPressSchechterGamma1
    return
  end function Modification_Function

  subroutine Compute_Common_Factors
    !% Precomputes some useful factors that are used in the modified Press-Schechter branching integrals.
    implicit none

    parentSigmaSquared           =parentSigma**2
    modificationG0Gamma2Factor   =modifiedPressSchechterG0*((parentDelta/parentSigma)**modifiedPressSchechterGamma2)
    branchingProbabilityPreFactor=sqrtTwoOverPi*parentHaloMass*modificationG0Gamma2Factor/parentSigma**modifiedPressSchechterGamma1
    return
  end subroutine Compute_Common_Factors

  subroutine Subresolution_Hypergeometric_Tabulate(x,xMinimumIn,xMaximumIn)
    !% Tabulate the hypergeometric term appearing in the subresolution merger fraction expression.
    use Hypergeometric_Functions
    use Galacticus_Error
    use Table_Labels
    implicit none
    double precision, intent(in   )           :: x
    double precision, intent(in   ), optional :: xMinimumIn        , xMaximumIn
    integer         , parameter               :: xCountPerDecade=10
    double precision                          :: xMinimum          , xMaximum
    integer                                   :: xCount            , i
    logical                                   :: tabulate
    
    tabulate=.false.
    if (.not.subresolutionHypergeometricInitialized) then
       tabulate=.true.
       if (present(xMinimumIn)) then
          xMinimum=xMinimumIn
       else
          xMinimum=min( 1.0d-9 ,     (x-1.0d0))
       end if
       if (present(xMaximumIn)) then
          xMaximum=xMaximumIn
       else
          xMaximum=max(12.5d+0,2.0d0*(x-1.0d0))
       end if
    else
       if     (                                               &
            &   (x-1.0d0) < subresolutionHypergeometric%x(+1) &
            &  .or.                                           &
            &   (x-1.0d0) > subresolutionHypergeometric%x(-1) &
            & ) then
          tabulate=.true.
          xMinimum=min(subresolutionHypergeometric%x(+1),      (x-1.0d0))
          xMaximum=max(subresolutionHypergeometric%x(-1),2.0d0*(x-1.0d0))
       end if
    end if
    if (tabulate) then
       xCount=int(log10(xMaximum/xMinimum)*dble(xCountPerDecade))+1
       if (.not.subresolutionHypergeometricInitialized) call subresolutionHypergeometric%destroy()
       call subresolutionHypergeometric%create(xMinimum,xMaximum,xCount,1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       do i=1,xCount
          call subresolutionHypergeometric%populate(                                                                                    &
               &                                    +sqrtTwoOverPi                                                                      &
               &                                    *(subresolutionHypergeometric%x(i)+1.0d0)**(+modifiedPressSchechterGamma1-1.0d0)    &
               &                                    /                                          (-modifiedPressSchechterGamma1+1.0d0)    &
               &                                    *Hypergeometric_2F1(                                                                &
               &                                                       [1.5d0,0.5d0-0.5d0*modifiedPressSchechterGamma1]               , &
               &                                                       [      1.5d0-0.5d0*modifiedPressSchechterGamma1]               , &
               &                                                       1.0d0/(subresolutionHypergeometric%x(i)+1.0d0)**2              , &
               &                                                       toleranceRelative=modifiedPressSchechterHypergeometricPrecision  &
               &                                                      )                                                               , &
               &                                    i                                                                                   &
               &                                   )
       end do
       subresolutionHypergeometricInitialized=.true.
    end if
    return
  end subroutine Subresolution_Hypergeometric_Tabulate
  
  subroutine Upper_Bound_Hypergeometric_Tabulate(mass,massResolution,massMinimumIn,massMaximumIn)
    !% Tabulate the hypergeometric term appearing in the upper bound branching probability rate expression.
    use Hypergeometric_Functions
    use Galacticus_Error
    use Table_Labels
    implicit none
    double precision                               , intent(in   )           :: mass                        , massResolution
    double precision                               , intent(in   ), optional :: massMinimumIn               , massMaximumIn
    class           (cosmologicalMassVarianceClass), pointer                 :: cosmologicalMassVariance_
    integer                                        , parameter               :: massCountPerDecade       =30
    double precision                                                         :: massMinimum                 , massMaximum
    integer                                                                  :: massCount                   , i
    logical                                                                  :: tabulate
    double precision                                                         :: massSigma                   , gammaEffective     , &
         &                                                                      halfMassSigma               , halfMassAlpha      , &
         &                                                                      resolutionMassSigma         , resolutionMassAlpha

    tabulate=.false.
    if (.not.upperBoundHypergeometricInitialized) then
       tabulate=.true.
       if (present(massMinimumIn)) then
          massMinimum=massMinimumIn
       else
          massMinimum=           2.0d0*massResolution
       end if
       if (present(massMaximumIn)) then
          massMaximum=massMaximumIn
       else
          massMaximum=max(1.0d16,2.0d0*mass          )
       end if
    else
       if     (                                       &
            &   mass < upperBoundHypergeometric%x(+1) &
            &  .or.                                   &
            &   mass > upperBoundHypergeometric%x(-1) &
            & ) then
          tabulate=.true.
          massMinimum=                                   2.0d0*massResolution
          massMaximum=max(upperBoundHypergeometric%x(-1),2.0d0*mass          )
       end if
    end if
    if (tabulate) then
       massResolutionTabulated=massResolution
       massCount=int(log10(massMaximum/massMinimum)*dble(massCountPerDecade))+1
       if (.not.upperBoundHypergeometricInitialized) call upperBoundHypergeometric%destroy()
       call upperBoundHypergeometric%create(massMinimum,massMaximum,massCount,1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       ! Evaluate sigma and alpha at the mass resolution.
       cosmologicalMassVariance_ => cosmologicalMassVariance()
       call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,resolutionMassSigma,resolutionMassAlpha)
       do i=1,massCount
          ! Evaluate sigmas and alpha.
          call           cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*upperBoundHypergeometric%x(i),halfMassSigma,halfMassAlpha)
          massSigma     =cosmologicalMassVariance_%rootVariance                      (      upperBoundHypergeometric%x(i)                            )
          gammaEffective=modifiedPressSchechterGamma1-1.0d0/halfMassAlpha
          call upperBoundHypergeometric%populate(                                                                                         &
               &                                    +sqrtTwoOverPi                                                                        &
               &                                    /massSigma                                                                            & 
               &                                    *(                                                                                    &
               &                                      +(halfMassSigma/massSigma)**(+gammaEffective-1.0d0)                                 &
               &                                      /                           (-gammaEffective+1.0d0)                                 &
               &                                      *Hypergeometric_2F1(                                                                &
               &                                                         [1.5d0,0.5d0-0.5d0*gammaEffective]                             , &
               &                                                         [      1.5d0-0.5d0*gammaEffective]                             , &
               &                                                         (massSigma/halfMassSigma)**2                                   , &
               &                                                         toleranceRelative=modifiedPressSchechterHypergeometricPrecision  &
               &                                                        )                                                                 &
               &                                      -(resolutionMassSigma/massSigma)**(+gammaEffective-1.0d0)                           &
               &                                      /                                 (-gammaEffective+1.0d0)                           &
               &                                      *Hypergeometric_2F1(                                                                &
               &                                                         [1.5d0,0.5d0-0.5d0*gammaEffective]                             , &
               &                                                         [      1.5d0-0.5d0*gammaEffective]                             , &
               &                                                         (massSigma/resolutionMassSigma)**2                             , &
               &                                                         toleranceRelative=modifiedPressSchechterHypergeometricPrecision  &
               &                                                        )                                                                 &
               &                                     )                                                                                  , &
               &                                    i                                                                                     &
               &                                   )
       end do
       upperBoundHypergeometricInitialized=.true.
    end if
    return
  end subroutine Upper_Bound_Hypergeometric_Tabulate
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Modified_Press_Schechter_Branching_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Modified_Press_Schechter_Branching_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) subresolutionHypergeometricInitialized,upperBoundHypergeometricInitialized
    if (subresolutionHypergeometricInitialized) write (stateFile) subresolutionHypergeometric%x(1),subresolutionHypergeometric%x(-1)
    if (   upperBoundHypergeometricInitialized) write (stateFile)    upperBoundHypergeometric%x(1),   upperBoundHypergeometric%x(-1),massResolutionTabulated
    return
  end subroutine Modified_Press_Schechter_Branching_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Modified_Press_Schechter_Branching_State_Restore</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Modified_Press_Schechter_Branching_State_Restore(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer                    , intent(in   ) :: stateFile
    type            (fgsl_file), intent(in   ) :: fgslStateFile
    double precision                           :: xMinimum     , xMaximum, &
         &                                        xResolution

    read (stateFile) subresolutionHypergeometricInitialized,upperBoundHypergeometricInitialized
    if (subresolutionHypergeometricInitialized) then
       subresolutionHypergeometricInitialized=.false.
       read (stateFile) xMinimum,xMaximum
       call Modified_Press_Schechter_Branching_Parameters()
       call Subresolution_Hypergeometric_Tabulate(xMinimum            ,xMinimum,xMaximum)
    end if
    if (   upperBoundHypergeometricInitialized) then
       upperBoundHypergeometricInitialized   =.false.
       read (stateFile) xMinimum,xMaximum,xResolution
       call Modified_Press_Schechter_Branching_Parameters()
       call Upper_Bound_Hypergeometric_Tabulate  (xMinimum,xResolution,xMinimum,xMaximum)
    end if
    return
  end subroutine Modified_Press_Schechter_Branching_State_Restore

end module Modified_Press_Schechter_Branching
