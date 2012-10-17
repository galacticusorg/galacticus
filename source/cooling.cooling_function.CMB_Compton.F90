!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  subroutine Cooling_Function_CMB_Compton_Initialize(coolingFunctionMethods,coolingFunctionsMatched)
    !% Initializes the ``atomic CIE cooling function from {\sc Cloudy}'' module.
    implicit none
    type(varying_string), intent(in   ) :: coolingFunctionMethods(:)
    integer,              intent(inout) :: coolingFunctionsMatched

    ! Check if this cooling function has been selected.
    if (any(coolingFunctionMethods == 'CMBCompton')) then
       functionSelected=.true.
       coolingFunctionsMatched=coolingFunctionsMatched+1
    end if

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
       coolingFunction=comptonRateNormalization*electronDensity*(radiation%temperature([radiationTypeCMB])**4)&
            &*(temperature-radiation%temperature([radiationTypeCMB]))

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
