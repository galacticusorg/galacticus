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

module Power_Spectrum_Window_Functions_Top_Hat
  !% Implements top-hat window function for power spectrum variance computation.
  implicit none
  private
  public :: Power_Spectrum_Window_Functions_Top_Hat_Initialize,Power_Spectrum_Window_Function_Top_Hat

contains

  !# <powerSpectrumWindowFunctionMethod>
  !#  <unitName>Power_Spectrum_Window_Functions_Top_Hat_Initialize</unitName>
  !# </powerSpectrumWindowFunctionMethod>
  subroutine Power_Spectrum_Window_Functions_Top_Hat_Initialize(powerSpectrumWindowFunctionMethod,Power_Spectrum_Window_Function_Get)
    !% Initializes the ``topHat'' power spectrum variance window function module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ),          intent(in   ) :: powerSpectrumWindowFunctionMethod
    procedure(double precision), pointer, intent(inout) :: Power_Spectrum_Window_Function_Get
    
    if (powerSpectrumWindowFunctionMethod == 'topHat') Power_Spectrum_Window_Function_Get => Power_Spectrum_Window_Function_Top_Hat
    return
  end subroutine Power_Spectrum_Window_Functions_Top_Hat_Initialize
  
  double precision function Power_Spectrum_Window_Function_Top_Hat(wavenumber,smoothingMass)
    !% Top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the power
    !% spectrum.
    use Numerical_Constants_Math
    use Cosmological_Parameters
    implicit none
    double precision, intent(in) :: wavenumber,smoothingMass
    double precision, parameter  :: xSeriesMaximum=1.0d-3
    double precision             :: topHatRadius,x,xSquared

    topHatRadius=((3.0d0/4.0d0/Pi)*smoothingMass/Omega_Matter()/Critical_Density())**(1.0d0/3.0d0)
    x=wavenumber*topHatRadius
    if      (x <= 0.0d0) then
       Power_Spectrum_Window_Function_Top_Hat=0.0d0
    else if (x <= xSeriesMaximum) then 
       ! Use a series expansion of the window function for small x.
       xSquared=x**2
       Power_Spectrum_Window_Function_Top_Hat=1.0d0 &
            & +xSquared*(-1.0d0/   10.0d0           &
            & +xSquared*(+1.0d0/  280.0d0           &
            & +xSquared*(-1.0d0/15120.0d0           &
            &           )))
    else
       ! For larger x, use the full expression.
       Power_Spectrum_Window_Function_Top_Hat=3.0d0*(dsin(x)-x*dcos(x))/(x**3)
    end if
    return
  end function Power_Spectrum_Window_Function_Top_Hat
  
end module Power_Spectrum_Window_Functions_Top_Hat

