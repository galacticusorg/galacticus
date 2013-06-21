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

!% Contains a module which implements top-hat window function for power spectrum variance computation.

module Power_Spectrum_Window_Functions_TH_KSS_Hybrid
  !% Implements top-hat window function for power spectrum variance computation.
  implicit none
  private
  public :: Power_Spectrum_Window_Functions_TH_KSS_Hybrid_Initialize,Power_Spectrum_Window_Function_TH_KSS_Hybrid

  ! Parameter controlling the normalization between mass and cut-off wavenumber.
  double precision :: cutOffNormalization                                    
  
  ! Parameter controlling the ratio of radii in k-space sharp and top-hat window functions.
  double precision :: powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio 
  
contains

  !# <powerSpectrumWindowFunctionMethod>
  !#  <unitName>Power_Spectrum_Window_Functions_TH_KSS_Hybrid_Initialize</unitName>
  !# </powerSpectrumWindowFunctionMethod>
  subroutine Power_Spectrum_Window_Functions_TH_KSS_Hybrid_Initialize(powerSpectrumWindowFunctionMethod,Power_Spectrum_Window_Function_Get,Power_Spectrum_Window_Function_Wavenumber_Maximum_Get)
    !% Initializes the ``topHatKSpaceSharpHybrid'' power spectrum variance window function module.
    use Numerical_Constants_Math
    use Cosmological_Parameters
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type            (varying_string                                                 ), intent(in   )          :: powerSpectrumWindowFunctionMethod                       
    procedure       (Power_Spectrum_Window_Function_TH_KSS_Hybrid                   ), intent(inout), pointer :: Power_Spectrum_Window_Function_Get                      
    procedure       (Power_Spectrum_Window_Function_Wavenumber_Maximum_TH_KSS_Hybrid), intent(inout), pointer :: Power_Spectrum_Window_Function_Wavenumber_Maximum_Get   
    character       (len=32                                                         )                         :: powerSpectrumWindowFunctionSharpKSpaceNormalizationText 
    double precision                                                                                          :: powerSpectrumWindowFunctionSharpKSpaceNormalization     
    
    if (powerSpectrumWindowFunctionMethod == 'topHatKSpaceSharpHybrid') then
       ! Set a pointer to our function.
       Power_Spectrum_Window_Function_Get                    => Power_Spectrum_Window_Function_TH_KSS_Hybrid
       Power_Spectrum_Window_Function_Wavenumber_Maximum_Get => Power_Spectrum_Window_Function_Wavenumber_Maximum_TH_KSS_Hybrid
       ! Get parameters. 
       !@ <inputParameter>
       !@   <name>powerSpectrumWindowFunctionSharpKSpaceNormalization</name>
       !@   <defaultValue>natural</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $a$ in the relation $k_{\rm s} = a/r_{\rm s}$, where $k_{\rm s}$ is the cut-off wavenumber for
       !@     the sharp $k$-space window function and $r_{\rm s}$ is the radius of a sphere (in real-space) enclosing the
       !@     requested smoothing mass. Alternatively, a value of {\tt natural} will be supplied in which case the normalization
       !@     is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass
       !@     of $M=6 \pi^2 \bar{\rho} k_{\rm s}^{-3}$.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>       
       call Get_Input_Parameter('powerSpectrumWindowFunctionSharpKSpaceNormalization'&
            &,powerSpectrumWindowFunctionSharpKSpaceNormalizationText,defaultValue="natural")
       if (powerSpectrumWindowFunctionSharpKSpaceNormalizationText == "natural") then
          cutOffNormalization=(6.0d0*Pi**2*Omega_Matter()*Critical_Density())**(1.0d0/3.0d0)
       else
          read (powerSpectrumWindowFunctionSharpKSpaceNormalizationText,*) powerSpectrumWindowFunctionSharpKSpaceNormalization
          cutOffNormalization=powerSpectrumWindowFunctionSharpKSpaceNormalization*(4.0d0*Pi*Omega_Matter()*Critical_Density()&
               &/3.0d0)**(1.0d0/3.0d0)
       end if
       !@ <inputParameter>
       !@   <name>powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta$ in the relation $r_{\rm s}=\beta r_{\rm th}$ between $k$-space sharp and top-hat window function radii in the hybrid window function used for computing the variance in the power spectrum.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>       
       call Get_Input_Parameter('powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio'&
            &,powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio,defaultValue=1.0d0)
    end if
    return
  end subroutine Power_Spectrum_Window_Functions_TH_KSS_Hybrid_Initialize
  
  double precision function Power_Spectrum_Window_Function_TH_KSS_Hybrid(wavenumber,smoothingMass)
    !% Computes a window function for calculations of the variance in the power spectrum. Specifically, uses a convolution of
    !% top-hat real-space and sharp $k$-space window functions. The top-hat radius is $r_{\rm th}$, while the $k$-space cut-off
    !% wavenumber is $k_{\rm s}=a/r_{\rm s}$, where $a=${\tt [powerSpectrumWindowFunctionSharpKSpaceNormalization]}. The two radii
    !% are chosen such that $r_{\rm th}^2 + r_{\rm s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$ and $r_{\rm s}=\beta r_{\rm th}$ where
    !% $\beta=${\tt [powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio]}.
    use Numerical_Constants_Math
    use Cosmological_Parameters
    implicit none
    double precision, intent(in   ) :: smoothingMass           , wavenumber                                
    double precision, parameter     :: xSeriesMaximum   =1.0d-3                                            
    double precision                :: kSpaceSharpRadius       , topHatRadius    , topHatWindowFunction, & 
         &                             totalRadius             , wavenumberCutOff, x                   , & 
         &                             xSquared                                                            
    
    ! Find the radius enclosing this mass.
    totalRadius=((3.0d0/4.0d0/Pi)*smoothingMass/Omega_Matter()/Critical_Density())**(1.0d0/3.0d0)

    ! Find the top-hat and sharp k-space radii, and the k-space wavenumber.
    topHatRadius     =totalRadius/sqrt(1.0d0+powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio**2)
    kSpaceSharpRadius=powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio*topHatRadius
    wavenumberCutOff =cutOffNormalization/kSpaceSharpRadius

    ! Compute the top-hat window function.
    x=wavenumber*topHatRadius
    if      (x <= 0.0d0) then
       topHatWindowFunction=0.0d0
    else if (x <= xSeriesMaximum) then 
       ! Use a series expansion of the window function for small x.
       xSquared=x**2
       topHatWindowFunction= 1.0d0                      &
            &               +xSquared*(-1.0d0/   10.0d0 &
            &               +xSquared*(+1.0d0/  280.0d0 &
            &               +xSquared*(-1.0d0/15120.0d0 &
            &               )))
    else
       ! For larger x, use the full expression.
       topHatWindowFunction=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if

    ! Compute k-space sharp window function.
    if      (wavenumber <=            0.0d0) then
       wavenumberCutOff=0.0d0
    else if (wavenumber <= wavenumberCutOff) then 
       wavenumberCutOff=1.0d0
    else
       wavenumberCutOff=0.0d0
    end if

    ! Compute the convolution (which is just the multiplication in k-space).
    Power_Spectrum_Window_Function_TH_KSS_Hybrid=wavenumberCutOff*topHatWindowFunction
    return
  end function Power_Spectrum_Window_Function_TH_KSS_Hybrid
  
  double precision function Power_Spectrum_Window_Function_Wavenumber_Maximum_TH_KSS_Hybrid(smoothingMass)
    !% Computes the maximum wavenumber at which the window function for calculations of the variance in the power spectrum is
    !% non-zero. Specifically, uses a convolution of top-hat real-space and sharp $k$-space window functions. The top-hat radius
    !% is $r_{\rm th}$, while the $k$-space cut-off wavenumber is $k_{\rm s}=a/r_{\rm s}$, where $a=${\tt
    !% [powerSpectrumWindowFunctionSharpKSpaceNormalization]}. The two radii are chosen such that $r_{\rm th}^2 + r_{\rm s}^2 = (3
    !% M / 4 \pi \bar{rho})^{1/3}$ and $r_{\rm s}=\beta r_{\rm th}$ where $\beta=${\tt
    !% [powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio]}.
    implicit none
    double precision, intent(in   ) :: smoothingMass                                  
    double precision, parameter     :: wavenumberLarge=1.0d30 !   Effective infinity. 
    
    Power_Spectrum_Window_Function_Wavenumber_Maximum_TH_KSS_Hybrid=wavenumberLarge
    return
  end function Power_Spectrum_Window_Function_Wavenumber_Maximum_TH_KSS_Hybrid
  
end module Power_Spectrum_Window_Functions_TH_KSS_Hybrid

