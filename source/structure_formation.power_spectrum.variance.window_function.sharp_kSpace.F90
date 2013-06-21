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

!% Contains a module which implements sharp in $k$-space window function for power spectrum variance computation.

module Power_Spectrum_Window_Functions_Sharp_kSpace
  !% Implements a sharp in $k$-space window function for power spectrum variance computation.
  implicit none
  private
  public :: Power_Spectrum_Window_Functions_Sharp_kSpace_Initialize

  ! Parameter controlling the normalization between mass and cut-off wavenumber.
  double precision :: cutOffNormalization  
                                        
contains

  !# <powerSpectrumWindowFunctionMethod>
  !#  <unitName>Power_Spectrum_Window_Functions_Sharp_kSpace_Initialize</unitName>
  !# </powerSpectrumWindowFunctionMethod>
  subroutine Power_Spectrum_Window_Functions_Sharp_kSpace_Initialize(powerSpectrumWindowFunctionMethod&
       &,Power_Spectrum_Window_Function_Get,Power_Spectrum_Window_Function_Wavenumber_Maximum_Get)
    !% Initializes the ``kSpaceSharp'' power spectrum variance window function module.
    use Numerical_Constants_Math
    use Cosmological_Parameters
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type            (varying_string                                                ), intent(in   )          :: powerSpectrumWindowFunctionMethod                        
    procedure       (Power_Spectrum_Window_Function_Sharp_kSpace                   ), intent(inout), pointer :: Power_Spectrum_Window_Function_Get                       
    procedure       (Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace), intent(inout), pointer :: Power_Spectrum_Window_Function_Wavenumber_Maximum_Get    
    character       (len=32                                                        )                         :: powerSpectrumWindowFunctionSharpKSpaceNormalizationText  
    double precision                                                                                         :: powerSpectrumWindowFunctionSharpKSpaceNormalization      
                                                                                                                                                                      
    if (powerSpectrumWindowFunctionMethod == 'kSpaceSharp') then
       ! Assign function pointer.
       Power_Spectrum_Window_Function_Get                    => Power_Spectrum_Window_Function_Sharp_kSpace
       Power_Spectrum_Window_Function_Wavenumber_Maximum_Get => Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace
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
    end if
    return
  end subroutine Power_Spectrum_Window_Functions_Sharp_kSpace_Initialize
  
  double precision function Power_Spectrum_Window_Function_Sharp_kSpace(wavenumber,smoothingMass)
    !% Top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the power
    !% spectrum. The normalization of the filter is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass
    !% of $M=6 \pi^2 \bar{\rho} k_{\rm s}^{-3}$ if $k_{\rm s}$ is the cut-off wavelength for the filter.
    implicit none
    double precision, intent(in   ) :: smoothingMass   , wavenumber  
    double precision                :: wavenumberCutOff              
                                                                  
    wavenumberCutOff=Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace(smoothingMass)
    if      (wavenumber <=            0.0d0) then
       Power_Spectrum_Window_Function_Sharp_kSpace=0.0d0
    else if (wavenumber <= wavenumberCutOff) then 
       Power_Spectrum_Window_Function_Sharp_kSpace=1.0d0
    else
       Power_Spectrum_Window_Function_Sharp_kSpace=0.0d0
    end if
    return
  end function Power_Spectrum_Window_Function_Sharp_kSpace
  
  double precision function Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace(smoothingMass)
    !% Top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the power
    !% spectrum. The normalization of the filter is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass
    !% of $M=6 \pi^2 \bar{\rho} k_{\rm s}^{-3}$ if $k_{\rm s}$ is the cut-off wavelength for the filter.
    implicit none
    double precision, intent(in   ) :: smoothingMass  
                                                   
    Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace=cutOffNormalization/smoothingMass**(1.0d0/3.0d0)
    return
  end function Power_Spectrum_Window_Function_Wavenumber_Maximum_Sharp_kSpace
  
end module Power_Spectrum_Window_Functions_Sharp_kSpace

