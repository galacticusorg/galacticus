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


!% Contains a module which implements a simple cooling time calculation (based on the ratio of the thermal energy density to the
!% volume cooling rate).

module Cooling_Times_Simple
  !% Implements a simple cooling time calculation (based on the ratio of the thermal energy density to the volume cooling rate).
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
    type(varying_string),          intent(in)    :: coolingTimeMethod
    procedure(),          pointer, intent(inout) :: Cooling_Time_Get,Cooling_Time_Density_Log_Slope_Get&
         &,Cooling_Time_Temperature_Log_Slope_Get
    
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
       !@ </inputParameter>
       call Get_Input_Parameter('coolingTimeSimpleDegreesOfFreedom',coolingTimeSimpleDegreesOfFreedom,defaultValue=3.0d0)
    end if
    return
  end subroutine Cooling_Time_Simple_Initialize

  double precision function Cooling_Time_Simple(temperature,density,abundances,radiation)
    !% Compute the cooling time (in Gyr) for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt abundances} and experiencing a radiation field as described by {\tt radiation}.
    use Numerical_Constants_Atomic
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Units
    use Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    implicit none
    double precision,          intent(in) :: temperature,density
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    ! Effectively infinite time (for arbitrarily long cooling times).
    double precision,          parameter  :: largeTime=1.0d10
    double precision                      :: numberDensityHydrogen,numberDensityAllSpecies,coolingFunction,energyDensityThermal

    ! Compute number density of hydrogen (in cm^-3).
    numberDensityHydrogen=density*abundances%hydrogenMassFraction()*massSolar/massHydrogenAtom/(hecto*megaParsec)**3

    ! Get the number density of all species.
    numberDensityAllSpecies=numberDensityHydrogen/abundances%hydrogenNumberFraction()

    ! Get the cooling function (in ergs cm^-3 s^-1).
    coolingFunction=Cooling_Function(temperature,numberDensityHydrogen,abundances,radiation)

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
  
  double precision function Cooling_Time_Density_Log_Slope_Simple(temperature,density,abundances,radiation)
    !% Return $\d\ln t_{\rm cool}/\d\ln \rho$ for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt abundances} and experiencing a radiation field as described by {\tt radiation}.
    use Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    implicit none
    double precision,          intent(in) :: temperature,density
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    Cooling_Time_Density_Log_Slope_Simple=1.0d0-Cooling_Function_Density_Log_Slope(temperature,density,abundances,radiation)
    return
  end function Cooling_Time_Density_Log_Slope_Simple
  
  double precision function Cooling_Time_Temperature_Log_Slope_Simple(temperature,density,abundances,radiation)
    !% Return $\d\ln t_{\rm cool}/\d\ln T$ for gas at the given {\tt temperature} (in Kelvin), {\tt density} (in $M_\odot$
    !% Mpc$^{-3}$), composition specified by {\tt abundances} and experiencing a radiation field as described by {\tt radiation}.
    use Abundances_Structure
    use Radiation_Structure
    use Cooling_Functions
    implicit none
    double precision,          intent(in) :: temperature,density
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    Cooling_Time_Temperature_Log_Slope_Simple=-Cooling_Function_Temperature_Log_Slope(temperature,density,abundances,radiation)
    return
  end function Cooling_Time_Temperature_Log_Slope_Simple
  
end module Cooling_Times_Simple
