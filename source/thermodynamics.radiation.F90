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


module Thermodynamics_Radiation
  !% Implements calculations of thermal radiation.
  private
  public :: Blackbody_Emission

  ! Define labels that specify if the computed radiance is per unit wavelength or per unit frequency.
  integer, public, parameter :: radianceTypeWavelength=0
  integer, public, parameter :: radianceTypeFrequency =1

contains

  double precision function Blackbody_Emission(wavelength,temperature,radianceType)
    !% Compute the Planck blackbody spectral radiance (defined per unit wavelength, in units of J s$^{-1}$ m$^{-2}$ sr$^{-1}$
    !% \AA$^{-1}$) or J s$^{-1}$ m$^{-2}$ sr$^{-1}$ Hz$^{-1}$ depending on the optional {\tt radianceType argument}). Input {\tt
    !% wavelength} is in Angstroms, input temperature is in Kelvin.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision, intent(in)           :: wavelength,temperature
    integer,          intent(in), optional :: radianceType
    double precision, parameter            :: exponentialArgumentMaximum=100.0d0
    double precision                       :: exponentialArgument
    integer                                :: radianceTypeActual

    ! Return zero emission for non-positive temperatures.
    if (temperature <= 0.0d0) then
       Blackbody_Emission=0.0d0
       return
    end if

    ! Determine what type of radiance to use.
    if (present(radianceType)) then
       radianceTypeActual=radianceType
    else
       radianceTypeActual=radianceTypeWavelength
    end if

    ! Compute the argument appearing in the exponential term.
    exponentialArgument=plancksConstant*speedLight*angstromsPerMeter/wavelength/boltzmannsConstant/temperature
    ! If the exponential argument is not too large, then compute the spectrum, otherwise return zero.
    if (exponentialArgument < exponentialArgumentMaximum) then
       select case (radianceTypeActual)
       case (radianceTypeWavelength)
          Blackbody_Emission=2.0d0*plancksConstant*speedLight**2*angstromsPerMeter**4/wavelength**5/(dexp(exponentialArgument)-1.0d0)
       case (radianceTypeFrequency )
          Blackbody_Emission=2.0d0*plancksConstant*speedLight   *angstromsPerMeter**3/wavelength**3/(dexp(exponentialArgument)-1.0d0)
       end select
    else
       Blackbody_Emission=0.0d0
    end if
    return
  end function Blackbody_Emission

end module Thermodynamics_Radiation
