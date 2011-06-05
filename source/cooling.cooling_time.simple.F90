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
