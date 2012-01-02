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


!% Contains a module which implements the \cite{madau_radiative_1995} calculation of the attenuation of spectra by the intergalactic medium.

module Stellar_Population_Spectra_Postprocess_Madau1995
  !% Implements the \cite{madau_radiative_1995} calculation of the attenuation of spectra by the intergalactic medium.
  use ISO_Varying_String

  public :: Stellar_Population_Spectra_Postprocess_Madau1995_Initialize

contains
  
  !# <stellarPopulationSpectraPostprocessMethod>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Madau1995_Initialize</unitName>
  !# </stellarPopulationSpectraPostprocessMethod>
  subroutine Stellar_Population_Spectra_Postprocess_Madau1995_Initialize(stellarPopulationSpectraPostprocessMethod,Stellar_Population_Spectra_Postprocess_Get)
    !% Initializes the ``Madau1995'' stellar spectrum postprocessing module.
    implicit none
    type(varying_string),                 intent(in)    :: stellarPopulationSpectraPostprocessMethod
    procedure(double precision), pointer, intent(inout) :: Stellar_Population_Spectra_Postprocess_Get
    
    if (stellarPopulationSpectraPostprocessMethod == 'Madau1995') Stellar_Population_Spectra_Postprocess_Get =>&
         & Stellar_Population_Spectra_Postprocess_Madau1995_Get
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Madau1995_Initialize

  double precision function Stellar_Population_Spectra_Postprocess_Madau1995_Get(wavelength,redshift)
    !% Computes the factor by which the spectrum of a galaxy at given {\tt redshift} is attenuated at the given {\tt wavelength}
    !% by the intervening intergalactic medium according to \cite{madau_radiative_1995}.
    use Numerical_Constants_Atomic
    use Factorials
    use Gamma_Functions
    implicit none
    double precision, intent(in)              :: wavelength,redshift
    double precision, dimension(9), parameter :: opticalDepthLymanLinesCoefficients=[0.00360d0,0.00170d0,0.00120d0,0.00093d0&
         &,0.00093d0,0.00093d0,0.00093d0,0.00093d0,0.00093d0]
    double precision, dimension(9)            :: opticalDepthLymanLines
    integer                                   :: iLine
    double precision                          :: opticalDepth,emissionFactor,continuumFactor,wavelengthObservedLymanContinuum
    
    ! Check if this is a zero redshift case.
    if (redshift <= 0.0d0) then
       ! It is, so return no attenuation.
       Stellar_Population_Spectra_Postprocess_Madau1995_Get=1.0d0
    else
       ! Compute the observed wavelength in units of the Lyman-continuum wavelength.
       wavelengthObservedLymanContinuum=wavelength*(1.0d0+redshift)/ionizationWavelengthHydrogen

       ! Compute contribution to optical depth from Lyman-series lines.
       forall(iLine=1:9)
          opticalDepthLymanLines(iLine)=opticalDepthLymanLinesCoefficients(iLine)*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0&
               &/dble((iLine+1)**2)))**3.46d0
       end forall

       ! Cumulate optical depth from lines.
       opticalDepth=0.0d0
       do iLine=1,9
          if (wavelengthObservedLymanContinuum < (1.0d0+redshift)/(1.0d0-1.0d0/dble((iLine+1)**2))) opticalDepth=opticalDepth&
               &+opticalDepthLymanLines(iLine)
       end do
       
       ! Add in continuum optical depth.
       if (wavelengthObservedLymanContinuum < 1.0d0+redshift) then
          emissionFactor =1.0d0+redshift
          continuumFactor=wavelengthObservedLymanContinuum
          opticalDepth=opticalDepth+0.25d0*(continuumFactor**3)*(emissionFactor**0.46d0-continuumFactor**0.46d0)+9.4d0&
               &*(continuumFactor**1.5d0)*(emissionFactor**0.18d0-continuumFactor**0.18d0)-0.7d0*(continuumFactor**3)*(1.0d0&
               &/continuumFactor**1.32d0-1.0d0/emissionFactor**1.32d0)-0.023d0*(emissionFactor**1.68d0-continuumFactor**1.68d0)
       end if
       
       ! Compute attenuation from optical depth.
       Stellar_Population_Spectra_Postprocess_Madau1995_Get=dexp(-opticalDepth) 
    end if
    return
  end function Stellar_Population_Spectra_Postprocess_Madau1995_Get

end module Stellar_Population_Spectra_Postprocess_Madau1995
