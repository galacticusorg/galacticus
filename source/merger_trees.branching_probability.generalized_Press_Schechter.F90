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

!% Contains a module which implements calculations of branching probabilties in generalized Press-Schechter theory.

module Generalized_Press_Schechter_Branching
  !% Implements calculations of branching probabilties in generalized Press-Schechter theory.
  use FGSL
  use Power_Spectrum
  use Numerical_Constants_Math
  implicit none
  private
  public :: Generalized_Press_Schechter_Branching_Initialize
  
  ! Parent halo shared variables.
  double precision :: parentHaloMass,parentTime,parentDTimeDDeltaCritical,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass
  !$omp threadprivate(parentHaloMass,parentTime,parentDTimeDDeltaCritical,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass)

  ! Accuracy parameter to ensure that steps in critical overdensity do not become too large.
  double precision :: generalizedPressSchechterDeltaStepMaximum

  ! Minimum mass to which subresolution fractions will be integrated.
  double precision :: generalizedPressSchechterMinimumMass

  ! Precomputed numerical factors.
  double precision, parameter  :: sqrtTwoOverPi=dsqrt(2.0d0/Pi)

  ! Branching probability integrand integration tolerance.
  double precision, parameter  :: branchingProbabilityIntegrandToleraceRelative=1.0d-2

  ! The maximum sigma that we expect to find.
  double precision             :: sigmaMaximum

  ! Record of whether we have tested the excursion set routines.
  logical                      :: excursionSetsTested=.false.

  ! Control for inclusion of smooth accretion rates.
  logical                      :: generalizedPressSchechterSmoothAccretion

contains
  
  !# <treeBranchingMethod>
  !#  <unitName>Generalized_Press_Schechter_Branching_Initialize</unitName>
  !# </treeBranchingMethod>
  subroutine Generalized_Press_Schechter_Branching_Initialize(treeBranchingMethod,Tree_Branching_Probability&
       &,Tree_Subresolution_Fraction,Tree_Branch_Mass,Tree_Maximum_Step)
    !% Initialize the generalized Press-Schechter branching routines.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: treeBranchingMethod
    procedure(double precision), pointer, intent(inout) :: Tree_Branching_Probability,Tree_Subresolution_Fraction,Tree_Branch_Mass&
         &,Tree_Maximum_Step
    
    if (treeBranchingMethod == 'generalizedPress-Schechter') then
       Tree_Branching_Probability  => Generalized_Press_Schechter_Branching_Probability
       Tree_Subresolution_Fraction => Generalized_Press_Schechter_Subresolution_Fraction
       Tree_Branch_Mass            => Generalized_Press_Schechter_Branch_Mass
       Tree_Maximum_Step           => Generalized_Press_Schechter_Branching_Maximum_Step
       !@ <inputParameter>
       !@   <name>generalizedPressSchechterDeltaStepMaximum</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Limits the step in $\delta_{\rm crit}$ when constructing merger trees using the generalized Press-Schechter branching algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('generalizedPressSchechterDeltaStepMaximum',generalizedPressSchechterDeltaStepMaximum,defaultValue &
            &=0.1d0)
       !@ <inputParameter>
       !@   <name>generalizedPressSchechterMinimumMass</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass to used in computing subresolution accretion rates when constructing merger trees using the generalized Press-Schechter branching algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('generalizedPressSchechterMinimumMass',generalizedPressSchechterMinimumMass,defaultValue &
            &=0.0d0)
       !@ <inputParameter>
       !@   <name>generalizedPressSchechterSmoothAccretion</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to include smooth accretion in subresolution accretion rates when constructing merger trees using the generalized Press-Schechter branching algorithm.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('generalizedPressSchechterSmoothAccretion',generalizedPressSchechterSmoothAccretion,defaultValue &
            &=.true.)
    end if
    return
  end subroutine Generalized_Press_Schechter_Branching_Initialize

  subroutine Excursion_Sets_Maximum_Sigma_Test()
    !% Make a call to excursion set routines with the maximum $\sigma$ that we will use to ensure that they can handle it.
    use Excursion_Sets_First_Crossings
    use Cosmology_Functions
    implicit none
    double precision :: presentTime,varianceMaximum,testResult

    !$omp critical (Excursion_Sets_Maximum_Sigma_Test)
    if (.not.excursionSetsTested) then
       presentTime    =Cosmology_Age(1.0d0)
       sigmaMaximum   =sigma_CDM(generalizedPressSchechterMinimumMass)
       varianceMaximum=sigmaMaximum**2
       testResult     =Excursion_Sets_First_Crossing_Probability(                      varianceMaximum,presentTime)
       testResult     =Excursion_Sets_First_Crossing_Rate       (0.5d0*varianceMaximum,varianceMaximum,presentTime)
       excursionSetsTested=.true.
    end if
    !$omp end critical (Excursion_Sets_Maximum_Sigma_Test)
    return
  end subroutine Excursion_Sets_Maximum_Sigma_Test

  double precision function Generalized_Press_Schechter_Branch_Mass(haloMass,deltaCritical,massResolution,probability)
    !% Determine the mass of one of the halos to which the given halo branches, given the branching probability,
    !% {\tt probability}. Typically, {\tt probabilityFraction} is found by multiplying {\tt
    !% Generalized\_Press\_Schechter\_Branching\_Probability()} by a random variable drawn in the interval 0--1 if a halo
    !% branches. This routine then finds the progenitor mass corresponding to this value.
    use, intrinsic :: ISO_C_Binding
    use ISO_Varying_String
    use Root_Finder
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    double precision,        intent(in) :: haloMass,deltaCritical,massResolution,probability
    type(fgsl_function),     save       :: rootFunction
    type(fgsl_root_fsolver), save       :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision,        parameter  :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-9
    double precision,        parameter  :: smallProbabilityFraction=1.0d-3
    type(c_ptr)                         :: parameterPointer
    type(varying_string)                :: message
    character(len=26)                   :: label

    ! Ensure excursion set calculations have sufficient range in sigma.
    call Excursion_Sets_Maximum_Sigma_Test()
    
    ! Initialize global variables.
    parentHaloMass        =haloMass
    parentSigma           =sigma_CDM(haloMass)
    parentDelta           =deltaCritical
    probabilityMinimumMass=massResolution
    probabilitySeek       =probability
    call Compute_Common_Factors

    ! Check that the root is bracketed.
    if     (                                                                                 &
         &     Generalized_Press_Schechter_Branch_Mass_Root(massResolution,parameterPointer) &
         &    *Generalized_Press_Schechter_Branch_Mass_Root(0.5d0*haloMass,parameterPointer) &
         &  >=                                                                               &
         &    0.0d0                                                                          &
         & ) then
       ! Warn about this situation.
       if (Galacticus_Verbosity_Level() >= verbosityWarn) then
          message="halo branching mass root is not bracketed in Generalized_Press_Schechter_Branch_Mass()"
          call Galacticus_Display_Message(message,verbosityWarn)
          write (label,'(e12.6,a1,e12.6)') massResolution,":",Generalized_Press_Schechter_Branch_Mass_Root(massResolution,parameterPointer)
          message=" => massMinimum:rootFunction(massMinimum) = "//trim(label)
          call Galacticus_Display_Message(message,verbosityWarn)
          write (label,'(e12.6,a1,e12.6)') 0.5d0*haloMass,":",Generalized_Press_Schechter_Branch_Mass_Root(0.5d0*haloMass,parameterPointer)
          message=" => massMaximum:rootFunction(massMaximum) = "//trim(label)
          call Galacticus_Display_Message(message,verbosityWarn)
          write (label,'(e12.6)') probability
          message=" =>                           probability = "//trim(label)
          call Galacticus_Display_Message(message,verbosityWarn)
       end if
       ! If the root function is positive at half of the parent halo mass then we have a binary split.
       if (Generalized_Press_Schechter_Branch_Mass_Root(0.5d0*haloMass,parameterPointer) >= 0.0d0) then
          ! Check that we are sufficiently close to zero. If we're not, it might indicate a problem.
          if     (                                                                                     &
               &   Generalized_Press_Schechter_Branch_Mass_Root(0.5d0*haloMass,parameterPointer)       &
               &  >                                                                                    &
               &   probability*smallProbabilityFraction                                                &
               & ) call Galacticus_Error_Report(                                                       &
               &                                "Generalized_Press_Schechter_Branch_Mass()"          , &
               &                                "numerical accuracy lost in root finding calculation"  &
               &                               )
          ! Return a binary split mass.
          Generalized_Press_Schechter_Branch_Mass=0.5d0*haloMass
          return
       end if
    end if

    Generalized_Press_Schechter_Branch_Mass=Root_Find(massResolution,0.5d0*haloMass,Generalized_Press_Schechter_Branch_Mass_Root&
         &,parameterPointer,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
    return
  end function Generalized_Press_Schechter_Branch_Mass

  function Generalized_Press_Schechter_Branch_Mass_Root(massMaximum,parameterPointer) bind(c)
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    real(c_double)                         :: Generalized_Press_Schechter_Branch_Mass_Root
    real(c_double),                  value :: massMaximum
    type(c_ptr),                     value :: parameterPointer
    type(fgsl_function)                    :: integrandFunction
    type(fgsl_integration_workspace)       :: integrationWorkspace

    Generalized_Press_Schechter_Branch_Mass_Root=probabilitySeek-Integrate(probabilityMinimumMass,massMaximum&
         &,Branching_Probability_Integrand_Generalized,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0&
         &,toleranceRelative=branchingProbabilityIntegrandToleraceRelative,integrationRule=FGSL_Integ_Gauss15)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Generalized_Press_Schechter_Branch_Mass_Root

  double precision function Generalized_Press_Schechter_Branching_Maximum_Step(haloMass,deltaCritical,massResolution)
    !% Return the maximum allowed step in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} should be allowed to take.
    use Numerical_Integration
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution

    Generalized_Press_Schechter_Branching_Maximum_Step=generalizedPressSchechterDeltaStepMaximum
    return
  end function Generalized_Press_Schechter_Branching_Maximum_Step

  double precision function Generalized_Press_Schechter_Branching_Probability(haloMass,deltaCritical,massResolution)
    !% Return the probability per unit change in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} will undergo a branching to progenitors with mass greater than {\tt massResolution}.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision,                 intent(in) :: haloMass,deltaCritical,massResolution
    type(c_ptr)                                  :: parameterPointer
    type(fgsl_function)                          :: integrandFunction
    type(fgsl_integration_workspace)             :: integrationWorkspace
    double precision                             :: massMinimum,massMaximum

    call Excursion_Sets_Maximum_Sigma_Test()
    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       parentHaloMass           =haloMass
       parentSigma              =sigma_CDM(haloMass)
       parentDelta              =deltaCritical
       call Compute_Common_Factors
       massMinimum=massResolution
       massMaximum=0.5d0*parentHaloMass
       Generalized_Press_Schechter_Branching_Probability=Integrate(massMinimum,massMaximum,Branching_Probability_Integrand_Generalized &
            &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=branchingProbabilityIntegrandToleraceRelative&
            &,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       Generalized_Press_Schechter_Branching_Probability=0.0d0
    end if
    return
  end function Generalized_Press_Schechter_Branching_Probability

  double precision function Generalized_Press_Schechter_Subresolution_Fraction(haloMass,deltaCritical,massResolution)
    !% Return the fraction of mass accreted in subresolution halos, i.e. those below {\tt massResolution}, per unit change in
    !% $\delta_{\rm crit}$ for a halo of mass {\tt haloMass} at time {\tt deltaCritical}. The integral is computed numerically.
    use, intrinsic :: ISO_C_Binding
    use Critical_Overdensity
    use Numerical_Integration
    use Excursion_Sets_First_Crossings
    use Merger_Tree_Branching_Modifiers
    implicit none
    double precision, intent(in)     :: haloMass,deltaCritical,massResolution
    double precision, save           :: resolutionSigma,massResolutionPrevious=-1.0d0
    !$omp threadprivate(resolutionSigma,massResolutionPrevious)
    double precision                 :: resolutionSigmaOverParentSigma,massMinimum,massMaximum
    type(c_ptr)                      :: parameterPointer
    type(fgsl_function)              :: integrandFunction
    type(fgsl_integration_workspace) :: integrationWorkspace

    call Excursion_Sets_Maximum_Sigma_Test()
    ! Get sigma and delta_critical for the parent halo.
    parentHaloMass           =haloMass
    parentSigma              =sigma_CDM(haloMass)
    parentDelta              =deltaCritical
    call Compute_Common_Factors
    
    ! If requested, compute the rate of smooth accretion.
    if (generalizedPressSchechterSmoothAccretion) then
       Generalized_Press_Schechter_Subresolution_Fraction=abs(parentDTimeDDeltaCritical)*Merger_Tree_Branching_Modifier(parentDelta&
            &,sigmaMaximum,parentSigma)*Excursion_Sets_Non_Crossing_Rate(parentSigmaSquared,parentTime)
    else
       Generalized_Press_Schechter_Subresolution_Fraction=0.0d0
    end if

    if (massResolution /= massResolutionPrevious) then
       resolutionSigma       =sigma_CDM(massResolution)
       massResolutionPrevious=massResolution
    end if
    resolutionSigmaOverParentSigma=resolutionSigma/parentSigma
    if (resolutionSigmaOverParentSigma >= 1.0d0) then
       massMinimum=generalizedPressSchechterMinimumMass
       massMaximum=massResolution
       Generalized_Press_Schechter_Subresolution_Fraction=Generalized_Press_Schechter_Subresolution_Fraction&
            &+Integrate(massMinimum,massMaximum,Subresolution_Fraction_Integrand_Generalized ,parameterPointer,integrandFunction&
            &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3 ,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       Generalized_Press_Schechter_Subresolution_Fraction=-1.0d0
    end if
    return
  end function Generalized_Press_Schechter_Subresolution_Fraction

  function Branching_Probability_Integrand_Generalized(childHaloMass,parameterPointer) bind(c)
    !% Integrand for the branching probability.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Branching_Probability_Integrand_Generalized
    real(c_double), value   :: childHaloMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: childSigma,childAlpha

    call sigma_CDM_Plus_Logarithmic_Derivative(childHaloMass,childSigma,childAlpha)
    Branching_Probability_Integrand_Generalized=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    return
  end function Branching_Probability_Integrand_Generalized
  
  function Subresolution_Fraction_Integrand_Generalized(childHaloMass,parameterPointer) bind(c)
    !% Integrand for the subresolution fraction.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Subresolution_Fraction_Integrand_Generalized
    real(c_double), value   :: childHaloMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: childSigma,childAlpha

    if (childHaloMass>0.0d0) then
       call sigma_CDM_Plus_Logarithmic_Derivative(childHaloMass,childSigma,childAlpha)
       Subresolution_Fraction_Integrand_Generalized=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)*(childHaloMass&
            &/parentHaloMass)
    else
       Subresolution_Fraction_Integrand_Generalized=0.0d0
    end if
    return
  end function Subresolution_Fraction_Integrand_Generalized
  
  double precision function Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    !% Progenitor mass function from Press-Schechter.
    implicit none
    double precision, intent(in) :: childHaloMass,childSigma,childAlpha

    Progenitor_Mass_Function=(parentHaloMass/childHaloMass**2)*Merging_Rate(childSigma,childAlpha)
    return
  end function Progenitor_Mass_Function
  
  double precision function Merging_Rate(childSigma,childAlpha)
    !% Computes the merging rate of dark matter halos in the generalized Press-Schechter algorithm. This ``merging rate'' is specifically defined as
    !% \begin{equation}
    !% {{\rm d}^2 f \over {\rm d} \ln M_{\rm child} {\rm d} \delta_{\rm c}} = 2 \sigma^2(M_{\rm child}) \left.{{\rm d} \ln \sigma \over {\rm d} \ln M}\right|_{M=M_{\rm child}} {{\rm d}t\over {\rm d}\delta_{\rm c}} {{\rm d}f_{12}\over {\rm d}t},
    !% \end{equation}
    !% where ${\rm d} f_{12}/{\rm d}t$ is the excursion set barrier crossing probabilty per unit time for the effective barrier
    !% $B^\prime(S_{\rm child}|S_{\rm parent},t)\equiv B(S_{\rm child},t-\delta t)-B(S_{\rm parent},t)$ in the limit $\delta t
    !% \rightarrow 0$.
    use Merger_Tree_Branching_Modifiers
    use Excursion_Sets_First_Crossings
    implicit none
    double precision, intent(in) :: childSigma,childAlpha
    double precision             :: childSigmaSquared

    childSigmaSquared=childSigma**2
    Merging_Rate=-2.0d0*Excursion_Sets_First_Crossing_Rate(parentSigmaSquared,childSigmaSquared,parentTime)*childSigmaSquared&
         &*abs(childAlpha)*parentDTimeDDeltaCritical*Merger_Tree_Branching_Modifier(parentDelta,childSigma,parentSigma)
   return
  end function Merging_Rate
    
  subroutine Compute_Common_Factors
    !% Precomputes some useful factors that are used in the generalized Press-Schechter branching integrals.
    use Critical_Overdensity
    implicit none

    parentSigmaSquared       =parentSigma**2
    parentTime               =Time_of_Collapse(parentDelta,parentHaloMass)
    parentDTimeDDeltaCritical=1.0d0/Critical_Overdensity_for_Collapse_Time_Gradient(parentTime,mass=parentHaloMass)
    return
  end subroutine Compute_Common_Factors

end module Generalized_Press_Schechter_Branching
