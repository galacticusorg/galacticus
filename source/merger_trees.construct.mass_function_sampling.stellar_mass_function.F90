!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements halo mass function sampling optimized to minimize variance in the model stellar mass function.

module Merger_Trees_Mass_Function_Sampling_Stellar_MF
  !% Implements halo mass function sampling optimized to minimize variance in the model stellar mass function.
  private
  public :: Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize

  ! Global halo mass used in integrand function.
  double precision :: massHalo

  ! Parameters of the mass function error model.
  double precision :: haloMassFunctionSamplingStellarMassFunctionErrorAlpha      , haloMassFunctionSamplingStellarMassFunctionErrorBeta       , &
       &              haloMassFunctionSamplingStellarMassFunctionErrorConstant   , haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth, &
       &              haloMassFunctionSamplingStellarMassFunctionErrorMassMaximum, haloMassFunctionSamplingStellarMassFunctionErrorMassMinimum, &
       &              haloMassFunctionSamplingStellarMassFunctionErrorMstar      , haloMassFunctionSamplingStellarMassFunctionErrorPhi0

contains

  !# <haloMassFunctionSamplingMethod>
  !#  <unitName>Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize</unitName>
  !# </haloMassFunctionSamplingMethod>
  subroutine Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize(haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get)
    !% Initializes the ``stellarMassFunction'' halo mass function sampling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: haloMassFunctionSamplingMethod
    procedure(Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF), intent(inout), pointer :: Merger_Tree_Construct_Mass_Function_Sampling_Get

    if (haloMassFunctionSamplingMethod == 'stellarMassFunction') then
       ! Set the pointer to point to our function.
       Merger_Tree_Construct_Mass_Function_Sampling_Get => Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF
       ! Read parameters of the assumed error model.
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorPhi0</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The value $\phi_0$ in a Schechter function, $\sigma(M) = \phi_0 (M/M_\star)^\alpha \exp(-[M/M_\star]^\beta)$, describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorPhi0',haloMassFunctionSamplingStellarMassFunctionErrorPhi0)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorAlpha</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The value $\alpha$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorAlpha',haloMassFunctionSamplingStellarMassFunctionErrorAlpha)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorBeta</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The value $\beta$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorBeta',haloMassFunctionSamplingStellarMassFunctionErrorBeta)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorMstar</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The value $M_\star$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorMstar',haloMassFunctionSamplingStellarMassFunctionErrorMstar)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorConstant</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The constant error contribution to the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorConstant',haloMassFunctionSamplingStellarMassFunctionErrorConstant)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The logarithmic width of bins in the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth',haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorMassMinimum</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum stellar mass to consider when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorMassMinimum',haloMassFunctionSamplingStellarMassFunctionErrorMassMinimum)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingStellarMassFunctionErrorMassMaximum</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum stellar mass to consider when computing the optimal sampling density function for tree masses.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingStellarMassFunctionErrorMassMaximum',haloMassFunctionSamplingStellarMassFunctionErrorMassMaximum)
    end if
    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Stellar_MF_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling_Stellar_MF(mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a power-law distribution.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Halo_Mass_Function
    use Galacticus_Meta_Compute_Times
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: mass                        , massMaximum                 , &
         &                                                         massMinimum                 , time
    double precision                            , parameter     :: toleranceAbsolute    =1.0d-3, toleranceRelative    =1.0d-2
    double precision                                            :: haloMassFunction            , logStellarMassMaximum       , &
         &                                                         logStellarMassMinimum       , treeComputeTime             , &
         &                                                         xi                          , xiIntegral
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    type            (c_ptr                     )                :: parameterPointer

    ! Get the halo mass function, defined per logarithmic interval in halo mass.
    haloMassFunction=mass*Halo_Mass_Function_Differential(time,mass)

    ! Compute the integral that appears in the "xi" function.
    massHalo             =mass
    logStellarMassMinimum=log10(haloMassFunctionSamplingStellarMassFunctionErrorMassMinimum)
    logStellarMassMaximum=log10(haloMassFunctionSamplingStellarMassFunctionErrorMassMaximum)
    xiIntegral           =Integrate(logStellarMassMinimum,logStellarMassMaximum,Xi_Integrand,parameterPointer,integrandFunction &
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
    use Conditional_Mass_Functions
    implicit none
    real            (kind=c_double               )          :: Xi_Integrand
    real            (kind=c_double               ), value   :: logStellarMass
    type            (c_ptr                       ), value   :: parameterPointer
    class           (conditionalMassFunctionClass), pointer :: conditionalMassFunction_
    double precision                                        :: conditionalMassFunctionVariance , stellarMass       , &
         &                                                     stellarMassFunctionObservedError, stellarMassMaximum, &
         &                                                     stellarMassMinimum

    ! Compute the stellar mass and range corresponding to data bins.
    stellarMass       =10.0d0** logStellarMass
    stellarMassMinimum=10.0d0**(logStellarMass-0.5d0*haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth)
    stellarMassMaximum=10.0d0**(logStellarMass+0.5d0*haloMassFunctionSamplingStellarMassFunctionErrorLogBinWidth)

    ! Compute the variance in the model conditional stellar mass function.
    conditionalMassFunction_        => conditionalMassFunction()
    conditionalMassFunctionVariance =  conditionalMassFunction_%massFunctionVariance(massHalo,stellarMassMinimum,stellarMassMaximum)

    ! Compute the error in the observed stellar mass. We use a simple Schechter function (plus minimum error) fit.
    stellarMassFunctionObservedError= haloMassFunctionSamplingStellarMassFunctionErrorPhi0                                                                             &
         &                           *exp(-(stellarMass/haloMassFunctionSamplingStellarMassFunctionErrorMstar)**haloMassFunctionSamplingStellarMassFunctionErrorBeta ) &
         &                           *     (stellarMass/haloMassFunctionSamplingStellarMassFunctionErrorMstar)**haloMassFunctionSamplingStellarMassFunctionErrorAlpha  &
         &                           +haloMassFunctionSamplingStellarMassFunctionErrorConstant

    ! Compute the integrand for the xi function integral.
    Xi_Integrand=conditionalMassFunctionVariance/stellarMassFunctionObservedError**2
    return
  end function Xi_Integrand

end module Merger_Trees_Mass_Function_Sampling_Stellar_MF
