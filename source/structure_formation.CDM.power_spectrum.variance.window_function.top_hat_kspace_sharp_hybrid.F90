!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
  subroutine Power_Spectrum_Window_Functions_TH_KSS_Hybrid_Initialize(powerSpectrumWindowFunctionMethod,Power_Spectrum_Window_Function_Get)
    !% Initializes the ``topHatKSpaceSharpHybrid'' power spectrum variance window function module.
    use Numerical_Constants_Math
    use Cosmological_Parameters
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ),          intent(in   ) :: powerSpectrumWindowFunctionMethod
    procedure(double precision), pointer, intent(inout) :: Power_Spectrum_Window_Function_Get
    character(len=32          )                         :: powerSpectrumWindowFunctionSharpKSpaceNormalizationText
    double precision                                    :: powerSpectrumWindowFunctionSharpKSpaceNormalization

    if (powerSpectrumWindowFunctionMethod == 'topHatKSpaceSharpHybrid') then
       ! Set a pointer to our function.
       Power_Spectrum_Window_Function_Get => Power_Spectrum_Window_Function_TH_KSS_Hybrid
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
    double precision, intent(in) :: wavenumber,smoothingMass
    double precision, parameter  :: xSeriesMaximum=1.0d-3
    double precision             :: totalRadius,topHatRadius,kSpaceSharpRadius,wavenumberCutOff,x,xSquared,topHatWindowFunction&
         &,sharpKWindowFunction

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
       topHatWindowFunction=3.0d0*(dsin(x)-x*dcos(x))/(x**3)
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
  
end module Power_Spectrum_Window_Functions_TH_KSS_Hybrid

