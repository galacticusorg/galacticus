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


!% Contains a module which implements the \cite{meiksin_colour_2006} calculation of the attenuation of spectra by the intergalactic medium.

module Stellar_Population_Spectra_Postprocess_Meiksin2006
  !% Implements the \cite{meiksin_colour_2006} calculation of the attenuation of spectra by the intergalactic medium.
  use ISO_Varying_String

  public :: Stellar_Population_Spectra_Postprocess_Meiksin2006_Initialize

contains
  
  !# <stellarPopulationSpectraPostprocessMethod>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Meiksin2006_Initialize</unitName>
  !# </stellarPopulationSpectraPostprocessMethod>
  subroutine Stellar_Population_Spectra_Postprocess_Meiksin2006_Initialize(stellarPopulationSpectraPostprocessMethod,Stellar_Population_Spectra_Postprocess_Get)
    !% Initializes the ``Meiksin2006'' stellar spectrum postprocessing module.
    implicit none
    type(varying_string),                 intent(in)    :: stellarPopulationSpectraPostprocessMethod
    procedure(double precision), pointer, intent(inout) :: Stellar_Population_Spectra_Postprocess_Get
    
    if (stellarPopulationSpectraPostprocessMethod == 'Meiksin2006') Stellar_Population_Spectra_Postprocess_Get =>&
         & Stellar_Population_Spectra_Postprocess_Meiksin2006_Get
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Meiksin2006_Initialize

  double precision function Stellar_Population_Spectra_Postprocess_Meiksin2006_Get(wavelength,redshift)
    !% Computes the factor by which the spectrum of a galaxy at given {\tt redshift} is attenuated at the given {\tt wavelength}
    !% by the intervening intergalactic medium according to \cite{meiksin_colour_2006}.
    use Numerical_Constants_Atomic
    use Factorials
    use Gamma_Functions
    implicit none
    double precision, intent(in)    :: wavelength,redshift
    ! Parameters of the Lyman-limit system distribution.
    double precision, parameter     :: N0   =0.25d0
    double precision, parameter     :: beta =1.50d0
    double precision, parameter     :: gamma=1.50d0
    double precision, dimension(31) :: opticalDepthLymanLines,redshiftLymanLines
    integer                         :: iLine
    double precision                :: seriesSolutionTermA,seriesSolutionTermB,wavelengthObservedLymanContinuum,nFactorial,opticalDepth

    ! Check if this is a zero redshift case.
    if (redshift <= 0.0d0) then
       ! It is, so return no attenuation.
       Stellar_Population_Spectra_Postprocess_Meiksin2006_Get=1.0d0
    else
       ! Compute the observed wavelength in units of the Lyman-continuum wavelength.
       wavelengthObservedLymanContinuum=wavelength*(1.0d0+redshift)/ionizationWavelengthHydrogen
       
       ! Evaluate redshifts of various Lyman-series lines.
       forall (iLine=3:9)
          redshiftLymanLines(iLine)=wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2))-1.0d0
       end forall
       
       ! Evaluate optical depths relative to Lyman-alpha.
       opticalDepthLymanLines(2)=1.0d0 ! By definition.
       if (redshiftLymanLines(3) < 3.0d0) then
          opticalDepthLymanLines(3)=0.348d0*(0.25d0*(1.0+redshiftLymanLines(3)))**0.3333d0
       else
          opticalDepthLymanLines(3)=0.348d0*(0.25d0*(1.0+redshiftLymanLines(3)))**0.1667d0
       end if
       if (redshiftLymanLines(4) < 3.0d0) then
          opticalDepthLymanLines(4)=0.179d0*(0.25d0*(1.0d0+redshiftLymanLines(4)))**0.3333d0
       else
          opticalDepthLymanLines(4)=0.179d0*(0.25d0*(1.0d0+redshiftLymanLines(4)))**0.1667d0
       end if
       if (redshiftLymanLines(5) < 3.0d0) then
          opticalDepthLymanLines(5)=0.109d0*(0.25d0*(1.0d0+redshiftLymanLines(5)))**0.3333d0
       else
          opticalDepthLymanLines(5)=0.109d0*(0.25d0*(1.0d0+redshiftLymanLines(5)))**0.1667d0
       end if
       opticalDepthLymanLines(6)=0.0722d0*(0.25d0*(1.0d0+redshiftLymanLines(6)))**0.3333d0
       opticalDepthLymanLines(7)=0.0508d0*(0.25d0*(1.0d0+redshiftLymanLines(7)))**0.3333d0
       opticalDepthLymanLines(8)=0.0373d0*(0.25d0*(1.0d0+redshiftLymanLines(8)))**0.3333d0
       opticalDepthLymanLines(9)=0.0283d0*(0.25d0*(1.0d0+redshiftLymanLines(9)))**0.3333d0
       forall (iLine=10:31)
          opticalDepthLymanLines(iLine)=opticalDepthLymanLines(9)*720.0d0/dble(iLine)/dble(iLine**2-1)
       end forall
       
       ! Scale optical depths by Lyman-alpha optical depth.
       if (redshift <= 4.0d0) then
          forall (iLine=2:31)
             opticalDepthLymanLines(iLine)=opticalDepthLymanLines(iLine)*0.00211d0*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2)))**3.70d0
          end forall
       else
          forall (iLine=2:31)
             opticalDepthLymanLines(iLine)=opticalDepthLymanLines(iLine)*0.00058d0*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2)))**4.50d0
          end forall
       end if
       
       ! Accumulate optical depths if line falls within the required redshift range.
       opticalDepth=0.0d0
       do iLine=2,31
          if (wavelengthObservedLymanContinuum < (1.0d0+redshift)/(1.0d0-1.0d0/(dble(iLine)**2))) opticalDepth=opticalDepth+opticalDepthLymanLines(iLine)
       end do
       
       if (wavelengthObservedLymanContinuum < (1.0d0+redshift)) then
          ! Add in photoelectric absorption contributions.
          seriesSolutionTermA=0.0d0
          seriesSolutionTermB=0.0d0
          do iLine=0,9
             nFactorial=Factorial(iLine)
             seriesSolutionTermA=seriesSolutionTermA+dble(-1**iLine)*(beta-1.0d0)/(dble(iLine)+1.0d0-beta)/nFactorial
             seriesSolutionTermB=seriesSolutionTermB+dble(-1**iLine)*(beta-1.0d0)*(((1.0d0+redshift)**(gamma+1.0d0-3.0d0 &
                  &*dble(iLine))*(wavelengthObservedLymanContinuum**(3.0d0*dble(iLine))) -(wavelengthObservedLymanContinuum**(gamma &
                  &+1.0d0))))/(dble(iLine)+1.0d0-beta)/(3.0d0*dble(iLine)-gamma-1.0d0)/nFactorial
          end do
          ! Add contribution due to Lyman-limit systems.
          opticalDepth=opticalDepth+N0*(dexp(Gamma_Function_Logarithmic(2.0d0-beta))-dexp(-1.0d0)-seriesSolutionTermA)*(((1.0d0+redshift)**(&
               &-3.0d0*(beta-1.0d0)+gamma+1.0d0))*(wavelengthObservedLymanContinuum**(3.0d0*(beta-1.0d0)))&
               &-(wavelengthObservedLymanContinuum**(gamma+1.0d0)))/(4.0d0+gamma-3.0d0*beta)-N0*seriesSolutionTermB
          ! Add contribution due to optically thin systems.
          opticalDepth=opticalDepth+0.805d0*(wavelengthObservedLymanContinuum**3)*(1.0d0/wavelengthObservedLymanContinuum-1.0d0 &
               &/(1.0d0+redshift))
       end if
       
       ! Compute attenuation from optical depth.
       Stellar_Population_Spectra_Postprocess_Meiksin2006_Get=dexp(-opticalDepth)
    end if
    return
  end function Stellar_Population_Spectra_Postprocess_Meiksin2006_Get

end module Stellar_Population_Spectra_Postprocess_Meiksin2006
