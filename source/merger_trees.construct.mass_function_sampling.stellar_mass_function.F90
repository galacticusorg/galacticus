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

!% Contains a module which implements halo mass function sampling optimized to minimize variance in the model stellar mass function.

module Merger_Trees_Mass_Function_Sampling_Stellar_MF
  !% Implements halo mass function sampling optimized to minimize variance in the model stellar mass function.
  private
  public :: Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize

  ! Global halo mass used in integrand function.
  double precision :: massHalo

contains

  !# <haloMassFunctionSamplingMethod>
  !#  <unitName>Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize</unitName>
  !# </haloMassFunctionSamplingMethod>
  subroutine Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize(haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get)
    !% Initializes the ``stellarMassFunction'' halo mass function sampling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: haloMassFunctionSamplingMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Construct_Mass_Function_Sampling_Get
    
    if (haloMassFunctionSamplingMethod == 'stellarMassFunction') Merger_Tree_Construct_Mass_Function_Sampling_Get => Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF
    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF(mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a power-law distribution.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Cosmology_Functions
    use Halo_Mass_Function
    use Galacticus_Meta_Compute_Times
    use Galacticus_Display
    use Numerical_Integration
    implicit none
    double precision,                 intent(in) :: mass,time,massMinimum,massMaximum
    double precision,                 parameter  :: logStellarMassMinimum=8.10d0,logStellarMassMaximum=11.88d0
    double precision,                 parameter  :: toleranceAbsolute=1.0d-3,toleranceRelative=1.0d-2
    double precision                             :: xi,xiIntegral,haloMassFunction,treeComputeTime
    type(fgsl_function)                          :: integrandFunction
    type(fgsl_integration_workspace)             :: integrationWorkspace
    type(c_ptr)                                  :: parameterPointer

    ! Check that the time given is consistent with the redshift of this constraint.
    if (Expansion_Factor(time) < 0.877d0) call Galacticus_Display_Message('WARNING - z≅0.07 stellar mass function constraint used&
         & define halo mass sampling density function, but halos are rooted at a significantly earlier time',verbosity&
         &=verbosityWarn)

    ! Get the halo mass function, defined per logarithmic interval in halo mass.
    haloMassFunction=mass*Halo_Mass_Function_Differential(time,mass)

    ! Compute the integral that appears in the "xi" function.
    massHalo        =mass
    xiIntegral      =Integrate(logStellarMassMinimum,logStellarMassMaximum,Xi_Integrand,parameterPointer,integrandFunction &
         &,integrationWorkspace,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative)
    call Integrate_Done(integrandFunction,integrationWorkspace)

    ! Compute the "xi" function.
    xi              =haloMassFunction**2*xiIntegral

    ! Get the time taken to compute a tree of this mass.
    treeComputeTime =Galacticus_Time_Per_Tree(mass)

    ! Compute the optimal weighting for trees of this mass.
    Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF=sqrt(xi/treeComputeTime)
    return
  end function Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF

  function Xi_Integrand(logStellarMass,parameterPointer) bind(c)
    !% The integrand appearing in the $\xi$ function.
    use, intrinsic :: ISO_C_Binding
    use Conditional_Stellar_Mass_Functions
    implicit none
    real(c_double)              :: Xi_Integrand
    real(c_double),   value     :: logStellarMass
    type(c_ptr),      value     :: parameterPointer
    double precision, parameter :: observedErrorPhi0   = 1.000d-3
    double precision, parameter :: observedErrorAlpha  =-0.300d0
    double precision, parameter :: observedErrorMstar  = 4.500d10
    double precision, parameter :: observedErrorMinimum= 1.000d-7
    double precision, parameter :: deltaLogStellarMass = 0.097d0 ! The width of bins in log of stellar mass used by Li & White (2009).
    double precision            :: stellarMass,stellarMassMinimum,stellarMassMaximum,conditionalMassFunctionVariance&
         &,stellarMassFunctionObservedError

    ! Compute the stellar mass and range corresponding to data bins.
    stellarMass       =10.0d0** logStellarMass
    stellarMassMinimum=10.0d0**(logStellarMass-0.5d0*deltaLogStellarMass)
    stellarMassMaximum=10.0d0**(logStellarMass+0.5d0*deltaLogStellarMass)

    ! Compute the variance in the model conditional stellar mass function.
    conditionalMassFunctionVariance=Cumulative_Conditional_Stellar_Mass_Function_Variance(massHalo,stellarMassMinimum,stellarMassMaximum)

    ! Compute the error in the observed stellar mass. We use a simple Schechter function (plus minimum error) fit to the
    ! observational results of Li & White (2009).
    stellarMassFunctionObservedError= observedErrorPhi0                                        &
         &                           *exp(-stellarMass/observedErrorMstar)                     &
         &                           *   ( stellarMass/observedErrorMstar)**observedErrorAlpha &
         &                           +observedErrorMinimum

    ! Compute the integrand for the xi function integral.
    Xi_Integrand=conditionalMassFunctionVariance/stellarMassFunctionObservedError**2
    return
  end function Xi_Integrand

end module Merger_Trees_Mass_Function_Sampling_Stellar_MF
