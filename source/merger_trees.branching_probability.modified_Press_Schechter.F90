!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of branching probabilties in modified Press-Schechter theory.

module Modified_Press_Schechter_Branching
  !% Implements calculations of branching probabilties in modified Press-Schechter theory.
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  use CDM_Power_Spectrum
  use Numerical_Constants_Math
  private
  public :: Modified_Press_Schechter_Branching_Initialize
  
  ! Parent halo shared variables.
  double precision :: parentHaloMass,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass,modificationG0Gamma2Factor
  !$omp threadprivate(parentHaloMass,parentSigma,parentSigmaSquared,parentDelta,probabilitySeek,probabilityMinimumMass,modificationG0Gamma2Factor)

  ! Parameters of the merger rate modification function.
  double precision :: modifiedPressSchechterG0,modifiedPressSchechterGamma1,modifiedPressSchechterGamma2

  ! Accuracy parameter to ensure that merger rate function (which is correct to 1st order) is sufficiently accurate.
  double precision :: modifiedPressSchechterFirstOrderAccuracy

  ! Precomputed numerical factors.
  double precision, parameter  :: sqrtTwoOverPi=dsqrt(2.0d0/Pi)

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
    type(varying_string),          intent(in)    :: treeBranchingMethod
    procedure(),          pointer, intent(inout) :: Tree_Branching_Probability,Tree_Subresolution_Fraction,Tree_Branch_Mass&
         &,Tree_Maximum_Step
    
    if (treeBranchingMethod.eq.'modified Press-Schechter') then
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
    use Root_Finder
    implicit none
    double precision,        intent(in) :: haloMass,deltaCritical,massResolution,probability
    type(fgsl_function),     save       :: rootFunction
    type(fgsl_root_fsolver), save       :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision,        parameter  :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-9
    type(c_ptr)                         :: parameterPointer

    probabilityMinimumMass=massResolution
    probabilitySeek       =probability
    Modified_Press_Schechter_Branch_Mass=Root_Find(massResolution,0.5d0*haloMass,Modified_Press_Schechter_Branch_Mass_Root&
         &,parameterPointer,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
    return
  end function Modified_Press_Schechter_Branch_Mass

  function Modified_Press_Schechter_Branch_Mass_Root(massMaximum,parameterPointer) bind(c)
    use Numerical_Integration
    implicit none
    real(c_double)                         :: Modified_Press_Schechter_Branch_Mass_Root
    real(c_double),                  value :: massMaximum
    type(c_ptr),                     value :: parameterPointer
    type(fgsl_function)                    :: integrandFunction
    type(fgsl_integration_workspace)       :: integrationWorkspace

    Modified_Press_Schechter_Branch_Mass_Root=probabilitySeek-Integrate(probabilityMinimumMass,massMaximum&
         &,Branching_Probability_Integrand,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0&
         &,toleranceRelative=1.0d-6)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Modified_Press_Schechter_Branch_Mass_Root

  double precision function Modified_Press_Schechter_Branching_Maximum_Step(haloMass,deltaCritical,massResolution)
    !% Return the maximum allowed step in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} should be allowed to take.
    use Numerical_Integration
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution
    double precision, parameter  :: largeStep=1.0d10 ! Effectively infinitely large step in w(=delta_crit).
    double precision             :: parentHalfMassSigma

    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       parentSigma        =sigma_CDM(      haloMass)
       parentHalfMassSigma=sigma_CDM(0.5d0*haloMass)
       Modified_Press_Schechter_Branching_Maximum_Step=modifiedPressSchechterFirstOrderAccuracy*dsqrt(2.0d0&
            &*(parentHalfMassSigma**2-parentSigma**2))
    else
       Modified_Press_Schechter_Branching_Maximum_Step=largeStep
    end if
    return
  end function Modified_Press_Schechter_Branching_Maximum_Step

  double precision function Modified_Press_Schechter_Branching_Probability(haloMass,deltaCritical,massResolution)
    !% Return the probability per unit change in $\delta_{\rm crit}$ that a halo of mass {\tt haloMass} at time {\tt
    !% deltaCritical} will undergo a branching to progenitors with mass greater than {\tt massResolution}.
    use Numerical_Integration
    implicit none
    double precision,                 intent(in) :: haloMass,deltaCritical,massResolution
    type(c_ptr)                                  :: parameterPointer
    type(fgsl_function)                          :: integrandFunction
    type(fgsl_integration_workspace)             :: integrationWorkspace
    double precision                             :: massMinimum,massMaximum

    ! Get sigma and delta_critical for the parent halo.
    if (haloMass>2.0d0*massResolution) then
       parentHaloMass=haloMass
       parentSigma=sigma_CDM(haloMass)
       parentDelta=deltaCritical
       call Compute_Common_Factors
       massMinimum=massResolution
       massMaximum=0.5d0*parentHaloMass
       Modified_Press_Schechter_Branching_Probability=Integrate(massMinimum,massMaximum,Branching_Probability_Integrand&
            &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
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
    double precision, intent(in) :: haloMass,deltaCritical,massResolution
    double precision, save       :: resolutionSigma,massResolutionPrevious=-1.0d0
    !$omp threadprivate(resolutionSigma,massResolutionPrevious)
    double precision             :: resolutionSigmaOverParentSigma,hyperGeometricFactor

    ! Get sigma and delta_critical for the parent halo.
    parentHaloMass=haloMass
    parentSigma   =sigma_CDM(haloMass)
    parentDelta   =deltaCritical
    call Compute_Common_Factors
    if (massResolution /= massResolutionPrevious) then
       resolutionSigma       =sigma_CDM(massResolution)
       massResolutionPrevious=massResolution
    end if
    resolutionSigmaOverParentSigma=resolutionSigma/parentSigma
    hyperGeometricFactor=Hypergeometric_2F1([1.5d0,0.5d0-0.5d0*modifiedPressSchechterGamma1],[1.5d0-0.5d0&
         &*modifiedPressSchechterGamma1],1.0d0/resolutionSigmaOverParentSigma**2)
    Modified_Press_Schechter_Subresolution_Fraction=sqrtTwoOverPi*(modificationG0Gamma2Factor/parentSigma) &
         &*(resolutionSigmaOverParentSigma**(modifiedPressSchechterGamma1-1.0d0))/(1.0d0-modifiedPressSchechterGamma1)&
         &*hyperGeometricFactor

    return
  end function Modified_Press_Schechter_Subresolution_Fraction

  function Branching_Probability_Integrand(childHaloMass,parameterPointer) bind(c)
    !% Integrand for the branching probability.
    implicit none
    real(c_double)          :: Branching_Probability_Integrand
    real(c_double), value   :: childHaloMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: childSigma,childAlpha

    call sigma_CDM_Plus_Logarithmic_Derivative(childHaloMass,childSigma,childAlpha)
    Branching_Probability_Integrand=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    return
  end function Branching_Probability_Integrand
  
  function Subresolution_Fraction_Integrand(childHaloMass,parameterPointer) bind(c)
    !% Integrand for the subresolution fraction.
    implicit none
    real(c_double)          :: Subresolution_Fraction_Integrand
    real(c_double), value   :: childHaloMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: childSigma,childAlpha

    if (childHaloMass>0.0d0) then
       call sigma_CDM_Plus_Logarithmic_Derivative(childHaloMass,childSigma,childAlpha)
       Subresolution_Fraction_Integrand=Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)*(childHaloMass&
            &/parentHaloMass)
    else
       Subresolution_Fraction_Integrand=0.0d0
    end if
    return
  end function Subresolution_Fraction_Integrand
  
  double precision function Progenitor_Mass_Function(childHaloMass,childSigma,childAlpha)
    !% Progenitor mass function from Press-Schechter.
    implicit none
    double precision, intent(in) :: childHaloMass,childSigma,childAlpha

    Progenitor_Mass_Function=(parentHaloMass/childHaloMass**2)*Merging_Rate(childSigma,childAlpha)&
         &*Modification_Function(childSigma)
    return
  end function Progenitor_Mass_Function
  
  double precision function Merging_Rate(childSigma,childAlpha)
    !% Merging rate from Press-Schechter.
    implicit none
    double precision, intent(in) :: childSigma,childAlpha
    double precision             :: childSigmaSquared

    childSigmaSquared=childSigma**2
    Merging_Rate=sqrtTwoOverPi*(childSigmaSquared/((childSigmaSquared-parentSigmaSquared)**1.5d0))*dabs(childAlpha)
    return
  end function Merging_Rate
  
  double precision function Modification_Function(childSigma)
    !% Empirical modification of the progenitor mass function from \cite{parkinson_generating_2008}.
    implicit none
    double precision, intent(in) :: childSigma

    Modification_Function=modificationG0Gamma2Factor*((childSigma/parentSigma)**modifiedPressSchechterGamma1)
    return
  end function Modification_Function
  
  subroutine Compute_Common_Factors
    !% Precomputes some useful factors that are used in the modified Press-Schechter branching integrals.
    implicit none

    parentSigmaSquared        =parentSigma**2
    modificationG0Gamma2Factor=modifiedPressSchechterG0*((parentDelta/parentSigma)**modifiedPressSchechterGamma2)
    return
  end subroutine Compute_Common_Factors

end module Modified_Press_Schechter_Branching
