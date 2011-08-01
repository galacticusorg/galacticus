!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which computes the contribution to the cooling function due to Compton cooling off of the cosmic microwave
!% background.

module Cooling_Functions_CMB_Compton
  !% Computes the contribution to the cooling function due to Compton cooling off of the cosmic microwave background.
  use ISO_Varying_String
  implicit none
  private
  public :: Cooling_Function_CMB_Compton_Initialize, Cooling_Function_CMB_Compton,&
       & Cooling_Function_Density_Slope_CMB_Compton, Cooling_Function_Temperature_Slope_CMB_Compton
  
  ! Flag indicating whether or not this cooling function is selected.
  logical                     :: functionSelected=.false.

contains
  
  !# <coolingFunctionMethods>
  !#  <unitName>Cooling_Function_CMB_Compton_Initialize</unitName>
  !# </coolingFunctionMethods>
  subroutine Cooling_Function_CMB_Compton_Initialize(coolingFunctionMethods)
    !% Initializes the ``atomic CIE cooling function from {\sc Cloudy}'' module.
    implicit none
    type(varying_string), intent(in) :: coolingFunctionMethods(:)
 
    ! Check if this cooling function has been selected.
    if (any(coolingFunctionMethods == 'CMB_Compton')) functionSelected=.true.

    return
  end subroutine Cooling_Function_CMB_Compton_Initialize

  !# <coolingFunctionCompute>
  !#   <unitName>Cooling_Function_CMB_Compton</unitName>
  !# </coolingFunctionCompute>
  subroutine Cooling_Function_CMB_Compton(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunction
    double precision,                  parameter   :: comptonRateNormalization=4.0d0*thomsonCrossSection*radiationConstant&
         &*boltzmannsConstant/electronMass/speedLight/ergs
    double precision                               :: electronDensity

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Get the electron density.
       electronDensity=Electron_Density(temperature,numberDensityHydrogen,abundances,radiation)

       ! Compute the Compton cooling rate.
       ! <gfortran 4.6> Would like to use the radiation%temperature() type-bound procedure form here, but type-bound procedures with optional arguments are buggy under gfortran 4.4
       coolingFunction=comptonRateNormalization*electronDensity*(Radiation_Temperature(radiation,[radiationTypeCMB])**4)&
            &*(temperature-Radiation_Temperature(radiation,[radiationTypeCMB]))

    else

       ! Not selected, return zero.
       coolingFunction=0.0d0

    end if

    return
  end subroutine Cooling_Function_CMB_Compton
  
  !# <coolingFunctionDensitySlopeCompute>
  !#   <unitName>Cooling_Function_Density_Slope_CMB_Compton</unitName>
  !# </coolingFunctionDensitySlopeCompute>
  subroutine Cooling_Function_Density_Slope_CMB_Compton(coolingFunctionDensitySlope,temperature,numberDensityHydrogen,abundances&
       &,chemicalDensities ,radiation)
    !% Return the gradient with respect to density of cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionDensitySlope
    double precision                               :: coolingFunction,electronDensityDensityLogSlope

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Get the cooling function.
       call Cooling_Function_CMB_Compton(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
       
       ! Get logarithmic slope of electron density with density.
       electronDensityDensityLogSlope=Electron_Density_Density_Log_Slope(temperature,numberDensityHydrogen,abundances&
            &,radiation)

       ! Depends only on the behavior of electron density with density.
       coolingFunctionDensitySlope=electronDensityDensityLogSlope*coolingFunction/numberDensityHydrogen
       
    else
       
       ! Not selected, return zero.
       coolingFunctionDensitySlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Density_Slope_CMB_Compton
  
  !# <coolingFunctionTemperatureSlopeCompute>
  !#   <unitName>Cooling_Function_Temperature_Slope_CMB_Compton</unitName>
  !# </coolingFunctionTemperatureSlopeCompute>
  subroutine Cooling_Function_Temperature_Slope_CMB_Compton(coolingFunctionTemperatureSlope,temperature,numberDensityHydrogen&
       &,abundances,chemicalDensities,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionTemperatureSlope
    double precision                               :: coolingFunction,electronDensityTemperatureLogSlope

    ! Check if this cooling function has been selected.
    if (functionSelected) then
              
       ! Get the cooling function.
       call Cooling_Function_CMB_Compton(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
       
       ! Get logarithmic slope of electron density with temperature.
       electronDensityTemperatureLogSlope=Electron_Density_Temperature_Log_Slope(temperature,numberDensityHydrogen,abundances&
            &,radiation)

       ! Compute the partial derivative of the cooling rate with respect to temperature.
       coolingFunctionTemperatureSlope=coolingFunction*(electronDensityTemperatureLogSlope/temperature+1.0d0/(temperature&
            &-Radiation_Temperature(radiation,[radiationTypeCMB])))
       
    else
       
       ! Not selected, return zero.
       coolingFunctionTemperatureSlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Temperature_Slope_CMB_Compton
  
end module Cooling_Functions_CMB_Compton
