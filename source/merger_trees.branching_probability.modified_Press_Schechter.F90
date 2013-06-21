!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Power_Spectra
  use Numerical_Constants_Math
  implicit none
  private
  public :: Modified_Press_Schechter_Branching_Initialize

  ! Parent halo shared variables.
  double precision            :: branchingProbabilityPreFactor                               , modificationG0Gamma2Factor  , &
       &                         parentDelta                                                 , parentHaloMass              , &
       &                         parentSigma                                                 , parentSigmaSquared          , &
       &                         probabilityMinimumMass                                      , probabilitySeek
  !$omp threadprivate(parentHaloMass,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass,modificationG0Gamma2Factor,branchingProbabilityPreFactor)
  ! Parameters of the merger rate modification function.
  double precision            :: modifiedPressSchechterG0                                    , modifiedPressSchechterGamma1, &
       &                         modifiedPressSchechterGamma2

  ! Accuracy parameter to ensure that merger rate function (which is correct to 1st order) is sufficiently accurate.
  double precision            :: modifiedPressSchechterFirstOrderAccuracy

  ! Precomputed numerical factors.
  double precision, parameter :: sqrtTwoOverPi                                =sqrt(2.0d0/Pi)

  ! Branching probability integrand integration tolerance.
  double precision, parameter :: branchingProbabilityIntegrandToleraceRelative=1.0d-3

contains

  !# <treeBranchingMethod>
  !#  <unitName>Modified_Press_Schechter_Branching_Initialize</unitName>
  !# </treeBranchingMethod>
  subroutine Modified_Press_Schechter_Branching_Initialize(treeBranchingMethod,Tree_Branching_Probability&
       &,Tree_Subresolution_Fraction,Tree_Branch_Mass,Tree_Maximum_Step)
    !% Initialize the modified Press-Schechter branching routines.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: treeBranchingMethod
    procedure(double precision), intent(inout), pointer :: Tree_Branch_Mass   , Tree_Branching_Probability , &
         &                                                 Tree_Maximum_Step  , Tree_Subresolution_Fraction

    if (treeBranchingMethod == 'modifiedPress-Schechter') then
       Tree_Branching_Probability  => Modified_Press_Schechter_Branching_Probability
       Tree_Subresolution_Fraction => Modified_Press_Schechter_Subresolution_Fraction
       Tree_Branch_Mass            => Modified_Press_Schechter_Branch_Mass
       Tree_Maximum_Step           => Modified_Press_Schechter_Branching_Maximum_Step
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
       !@     Limits the step in $\delta_{\rm crit}$ when constructing merger trees using the \cite{parkinson_generating_2008}
       !@     algorithm, so that it never exceeds {\tt
       !@     modifiedPressSchechterFirstOrderAccuracy}$\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('modifiedPressSchechterFirstOrderAccuracy',modifiedPressSchechterFirstOrderAccuracy,defaultValue &
            &=0.1d0)
    end if
    return
  end subroutine Modified_Press_Schechter_Branching_Initialize

  double precision function Modified_Press_Schechter_Branch_Mass(haloMass,deltaCritical,massResolution,probability)
    !% Determine the mass of one of the halos to which the given halo branches, given the branching probability,
    !% {\tt probability}. Typically, {\tt probabilityFraction} is found by multiplying {\tt
    !% Modified\_Press\_Schechter\_Branching\_Probability()} by a random variable drawn in the interval 0--1 if a halo
    !% branches. This routine then finds the progenitor mass corresponding to this value.
    use, intrinsic :: ISO_C_Binding
    use Root_Finder
    implicit none
    double precision            , intent(in   ) :: deltaCritical          , haloMass                , &
         &                                         massResolution         , probability
    double precision            , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    type            (rootFinder), save          :: finder
    !$omp threadprivate(finder)
    ! Initialize global variables.
    parentHaloMass        =haloMass
    parentSigma           =Cosmological_Mass_Root_Variance(haloMass)
    parentDelta           =deltaCritical
    probabilityMinimumMass=massResolution
    probabilitySeek       =probability
    call Compute_Common_Factors
    ! Check the sign of the root function at half the halo mass.
    if (Modified_Press_Schechter_Branch_Mass_Root(0.5d0*haloMass) >= 0.0d0) then
       ! The root function is zero, or very close to it (which can happen due to rounding errors
       ! occasionally). Therefore we have an almost perfect binary split.
       Modified_Press_Schechter_Branch_Mass=0.5d0*haloMass
    else
       ! Initialize our root finder.
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(Modified_Press_Schechter_Branch_Mass_Root)
          call finder%tolerance   (toleranceAbsolute,toleranceRelative      )
       end if
       ! Split is not binary - seek the actual mass of the smaller progenitor.
       Modified_Press_Schechter_Branch_Mass=finder%find(rootRange=[massResolution,0.5d0*haloMass])
    end if
    return
  end function Modified_Press_Schechter_Branch_Mass

  double precision function Modified_Press_Schechter_Branch_Mass_Root(massMaximum)
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: massMaximum
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    type            (c_ptr                     )                :: parameterPointer

    Modified_Press_Schechter_Branch_Mass_Root=probabilitySeek-branchingProbabilityPreFactor*Integrate(probabilityMinimumMass&
         &,massMaximum ,Branching_Probability_Integrand,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute&
         &=0.0d0 ,toleranceRelative=branchingProbabilityIntegrandToleraceRelative,integrationRule=FGSL_Integ_Gauss15)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Modified_Press_Schechter_Branch_Mass_Root

  double precision function Modified_Press_Schechter_Branching_Maximum_Step(haloMass,deltaCritical,massResolution)
    !% Return the maximum allowed step in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} should be allowed to take.
    implicit none
    double precision, intent(in   ) :: deltaCritical             , haloMass                                                         , massResolution
    double precision, parameter     :: largeStep          =1.0d10           !   Effectively infinitely large step in w(=delta_crit).
    double precision                :: parentHalfMassSigma

    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       parentSigma        =Cosmological_Mass_Root_Variance(      haloMass)
       parentHalfMassSigma=Cosmological_Mass_Root_Variance(0.5d0*haloMass)
       Modified_Press_Schechter_Branching_Maximum_Step=modifiedPressSchechterFirstOrderAccuracy*sqrt(2.0d0&
            &*(parentHalfMassSigma**2-parentSigma**2))
    else
       Modified_Press_Schechter_Branching_Maximum_Step=largeStep
    end if
    return
  end function Modified_Press_Schechter_Branching_Maximum_Step

  double precision function Modified_Press_Schechter_Branching_Probability(haloMass,deltaCritical,massResolution)
    !% Return the probability per unit change in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} will undergo a branching to progenitors with mass greater than {\tt massResolution}.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: deltaCritical       , haloMass   , massResolution
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: massMaximum         , massMinimum

    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       parentHaloMass=haloMass
       parentSigma=Cosmological_Mass_Root_Variance(haloMass)
       parentDelta=deltaCritical
       call Compute_Common_Factors
       massMinimum=massResolution
       massMaximum=0.5d0*parentHaloMass
       Modified_Press_Schechter_Branching_Probability=branchingProbabilityPreFactor*Integrate(massMinimum,massMaximum,Branching_Probability_Integrand &
            &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=branchingProbabilityIntegrandToleraceRelative&
            &,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       Modified_Press_Schechter_Branching_Probability=0.0d0
    end if
    return
  end function Modified_Press_Schechter_Branching_Probability

  double precision function Modified_Press_Schechter_Subresolution_Fraction(haloMass,deltaCritical,massResolution)
    !% Return the fraction of mass accreted in subresolution halos, i.e. those below {\tt massResolution}, per unit change in
    !% $\delta_{\rm crit}$ for a halo of mass {\tt haloMass} at time {\tt deltaCritical}. The integral is computed analytically in
    !% terms of the $_2F_1$ hypergeometric function.
    use Hypergeometric_Functions
    implicit none
    double precision, intent(in   ) :: deltaCritical                , haloMass                      , &
         &                             massResolution
    double precision, save          :: massResolutionPrevious=-1.0d0, resolutionSigma
    !$omp threadprivate(resolutionSigma,massResolutionPrevious)
    double precision                :: hyperGeometricFactor         , resolutionSigmaOverParentSigma

    ! Get sigma and delta_critical for the parent halo.
    parentHaloMass=haloMass
    parentSigma   =Cosmological_Mass_Root_Variance(haloMass)
    parentDelta   =deltaCritical
    call Compute_Common_Factors
    if (massResolution /= massResolutionPrevious) then
       resolutionSigma       =Cosmological_Mass_Root_Variance(massResolution)
       massResolutionPrevious=massResolution
    end if
    resolutionSigmaOverParentSigma=resolutionSigma/parentSigma
    if (resolutionSigmaOverParentSigma > 1.0d0) then
       hyperGeometricFactor=Hypergeometric_2F1([1.5d0,0.5d0-0.5d0*modifiedPressSchechterGamma1],[1.5d0-0.5d0&
            &*modifiedPressSchechterGamma1],1.0d0/resolutionSigmaOverParentSigma**2)
       Modified_Press_Schechter_Subresolution_Fraction=sqrtTwoOverPi*(modificationG0Gamma2Factor/parentSigma) &
            &*(resolutionSigmaOverParentSigma**(modifiedPressSchechterGamma1-1.0d0))/(1.0d0-modifiedPressSchechterGamma1)&
            &*hyperGeometricFactor
    else
       Modified_Press_Schechter_Subresolution_Fraction=-1.0d0
    end if
    return
  end function Modified_Press_Schechter_Subresolution_Fraction

  function Branching_Probability_Integrand(childHaloMass,parameterPointer) bind(c)
    !% Integrand for the branching probability.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: Branching_Probability_Integrand
    real(kind=c_double), value :: childHaloMass
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: childAlpha                     , childSigma

    call Cosmological_Mass_Root_Variance_Plus_Logarithmic_Derivative(childHaloMass,childSigma,childAlpha)
    Branching_Probability_Integrand=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    return
  end function Branching_Probability_Integrand

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
    Merging_Rate=(childSigmaSquared/((childSigmaSquared-parentSigmaSquared)**1.5d0))*abs(childAlpha)
    return
  end function Merging_Rate

  double precision function Modification_Function(childSigma)
    !% Empirical modification of the progenitor mass function from
    !% \cite{parkinson_generating_2008}. The constant factors of $G_0 (\delta_{\rm
    !% p}/\sigma_{\rm p})^{\gamma_2}$ and $1/\sigma_{\rm p}^{\gamma_1}$ are not included
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

end module Modified_Press_Schechter_Branching
