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

!% Contains a module which implements a simple cooling time calculation (based on the ratio of the thermal energy density to the
!% volume cooling rate).

module Cooling_Times_Simple
  !% Implements a simple cooling time calculation (based on the ratio of the thermal energy density to the volume cooling rate).
  implicit none
  private
  public :: Cooling_Time_Simple_Initialize

  ! Number of degrees of freedom assumed for the cooling time calculation.
  double precision :: coolingTimeSimpleDegreesOfFreedom

contains

  !# <coolingTimeMethod>
  !#  <unitName>Cooling_Time_Simple_Initialize</unitName>
  !# </coolingTimeMethod>
  subroutine Cooling_Time_Simple_Initialize(coolingTimeMethod,Cooling_Time_Get,Cooling_Time_Density_Log_Slope_Get,Cooling_Time_Temperature_Log_Slope_Get)
    !% Initializes the ``simple'' cooling time module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: coolingTimeMethod
    procedure(Cooling_Time_Simple), intent(inout), pointer :: Cooling_Time_Density_Log_Slope_Get
    procedure(Cooling_Time_Density_Log_Slope_Simple), intent(inout), pointer :: Cooling_Time_Get
    procedure(Cooling_Time_Temperature_Log_Slope_Simple), intent(inout), pointer :: Cooling_Time_Temperature_Log_Slope_Get

    if (coolingTimeMethod == 'simple') then
       Cooling_Time_Get => Cooling_Time_Simple
       Cooling_Time_Density_Log_Slope_Get => Cooling_Time_Density_Log_Slope_Simple
       Cooling_Time_Temperature_Log_Slope_Get => Cooling_Time_Temperature_Log_Slope_Simple
       !@ <inputParameter>
       !@   <name>coolingTimeSimpleDegreesOfFreedom</name>
       !@   <defaultValue>3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Number of degrees of freedom to assume when computing the energy density of cooling gas in the ``simple'' cooling time module.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingTimeSimpleDegreesOfFreedom',coolingTimeSimpleDegreesOfFreedom,defaultValue=3.0d0)
    end if
    return
  end subroutine Cooling_Time_Simple_Initialize

  double precision function Cooling_Time_Simple(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Compute the cooling time (in Gyr) for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt gasAbundances} and experiencing a radiation field as described by {\tt radiation}.
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    use Chemical_States
    implicit none
    double precision                    , intent(in   ) :: density                       , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation
    ! Effectively infinite time (for arbitrarily long cooling times).
    double precision                    , parameter     :: largeTime              =1.0d10
    double precision                                    :: coolingFunction               , energyDensityThermal , &
         &                                                 numberDensityAllSpecies       , numberDensityHydrogen

    ! Compute number density of hydrogen (in cm^-3).
    numberDensityHydrogen=density*gasAbundances%hydrogenMassFraction()*massSolar/massHydrogenAtom/(hecto*megaParsec)**3

    ! Get the number density of all species, including electrons.
    numberDensityAllSpecies= numberDensityHydrogen/gasAbundances%hydrogenNumberFraction() &
         &                  +Electron_Density(temperature,numberDensityHydrogen,gasAbundances,radiation)

    ! Get the cooling function (in ergs cm^-3 s^-1).
    coolingFunction=Cooling_Function(temperature,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)

    ! Determine the thermal energy density of the gas (in ergs cm^-3).
    energyDensityThermal=(coolingTimeSimpleDegreesOfFreedom/2.0d0)*boltzmannsConstant*temperature*numberDensityAllSpecies/ergs

    ! Compute the cooling time.
    if (coolingFunction > 0.0d0) then
       Cooling_Time_Simple=energyDensityThermal/coolingFunction/gigaYear
    else
       Cooling_Time_Simple=largeTime
    end if
    return
  end function Cooling_Time_Simple

  double precision function Cooling_Time_Density_Log_Slope_Simple(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Return $\d\ln t_{\rm cool}/\d\ln \rho$ for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt gasAbundances} and experiencing a radiation field as described by {\tt radiation}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    implicit none
    double precision                    , intent(in   ) :: density          , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation

    Cooling_Time_Density_Log_Slope_Simple=1.0d0-Cooling_Function_Density_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)
    return
  end function Cooling_Time_Density_Log_Slope_Simple

  double precision function Cooling_Time_Temperature_Log_Slope_Simple(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Return $\d\ln t_{\rm cool}/\d\ln T$ for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt gasAbundances} and experiencing a radiation field as described by {\tt radiation}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    implicit none
    double precision                    , intent(in   ) :: density          , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation

    Cooling_Time_Temperature_Log_Slope_Simple=-Cooling_Function_Temperature_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)
    return
  end function Cooling_Time_Temperature_Log_Slope_Simple

end module Cooling_Times_Simple
