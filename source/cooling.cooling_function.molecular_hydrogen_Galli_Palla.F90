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


!% Contains a module which computes the contribution to the cooling function from molecular hydrogen using the cooling function of
!% \cite{galli_chemistry_1998}.

module Cooling_Functions_Molecular_Hydrogen_Galli_Palla
  !% Computes the contribution to the cooling function from molecular hydrogen using the cooling function of
  !\cite{galli_chemistry_1998}.
  use ISO_Varying_String
  private
  public :: Cooling_Function_Molecular_Hydrogen_GP_Initialize, Cooling_Function_Molecular_Hydrogen_GP,&
       & Cooling_Function_Density_Slope_Molecular_Hydrogen_GP, Cooling_Function_Temperature_Slope_Molecular_Hydrogen_GP
  
  ! Flag indicating whether or not this cooling function is selected.
  logical :: functionSelected=.false.

  ! Indices of "molecules" (includes atoms, atomic ions and electrons also) used in this cooling function.
  integer :: electronMoleculeIndex,atomicHydrogenMoleculeIndex,molecularHydrogenCationMoleculeIndex,molecularHydrogenMoleculeIndex
  
  ! Parameters for Hollenbach & McKee cooling function fits.
  double precision, parameter :: coolingFunctionRotationalLambda1      = 9.50d-22
  double precision, parameter :: coolingFunctionRotationalLambda2      = 3.00d-24
  double precision, parameter :: coolingFunctionRotationalTemperature1 = 0.13d0
  double precision, parameter :: coolingFunctionRotationalTemperature2 = 0.51d0
  double precision, parameter :: coolingFunctionRotationalExponent1    = 3.76d0
  double precision, parameter :: coolingFunctionRotationalExponent2    = 2.10d0
  double precision, parameter :: coolingFunctionRotationalCoefficient1 = 0.12d0
  double precision, parameter :: coolingFunctionVibrationalLambda1     = 6.70d-19
  double precision, parameter :: coolingFunctionVibrationalLambda2     = 1.60d-18
  double precision, parameter :: coolingFunctionVibrationalTemperature1= 5.86d0
  double precision, parameter :: coolingFunctionVibrationalTemperature2=11.70d0

  ! Parameters for low-density limit cooling function.
  double precision, parameter, dimension(0:4) :: coolingFunctionLowDensityLimitCoefficient=[ &
       & -103.0000d0,+97.5900d0,-48.0500d0,+10.8000d0,-0.9032d0                              &
       &                                                                                   ]

  ! Parameters for H_2^+ - e^- cooling function.
  double precision, parameter, dimension(0:2) :: coolingFunctionH2PlusElectronCoefficient=[ &
       & -33.3299d0,+5.56465d0,-4.67461d-1                                                  &
       &                                                                                   ]

  ! Parameters for H_2^+ - H cooling function.
  double precision, parameter, dimension(0:2) :: coolingFunctionH2PlusHCoefficient        =[ &
       & -35.2804d0,+5.86234d0,-5.12276d-1                                                   &
       &                                                                                   ]

contains
  
  !# <coolingFunctionMethods>
  !#  <unitName>Cooling_Function_Molecular_Hydrogen_GP_Initialize</unitName>
  !# </coolingFunctionMethods>
  subroutine Cooling_Function_Molecular_Hydrogen_GP_Initialize(coolingFunctionMethods)
    !% Initializes the ``molecular hydrogen Galli \& Palla'' cooling function module.
    use Molecular_Abundances_Structure
    implicit none
    type(varying_string), intent(in) :: coolingFunctionMethods(:)
 
    ! Check if this cooling function has been selected.
    if (any(coolingFunctionMethods == 'molecularHydrogenGalliPalla')) then
       ! Flag that this cooling function has been selected.
       functionSelected=.true.

       ! Get the indices of molecules that will be used.
       electronMoleculeIndex               =Molecules_Index("Electron"               )
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       molecularHydrogenMoleculeIndex      =Molecules_Index("MolecularHydrogen"      )
    end if

    return
  end subroutine Cooling_Function_Molecular_Hydrogen_GP_Initialize

  !# <coolingFunctionCompute>
  !#   <unitName>Cooling_Function_Molecular_Hydrogen_GP</unitName>
  !# </coolingFunctionCompute>
  subroutine Cooling_Function_Molecular_Hydrogen_GP(coolingFunction,temperature,numberDensityHydrogen,abundances&
       &,molecularDensities,radiation)
    !% Return the cooling function due to molecular hydrogen using the cooling function of \cite{galli_chemistry_1998} (which
    !% refers to the local thermodynamic equilibrium cooling function of \cite{hollenbach_molecule_1979}). Cooling functions
    !% involving H$_2^+$ are computed using polynomial fits to the results of \cite{suchkov_cooling_1978} found by Andrew
    !% Benson by measuring curves from the original paper.
    use Ionization_States
    use Abundances_Structure
    use Molecular_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision,                   intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),          intent(in)  :: abundances
    type(molecularAbundancesStructure), intent(in)  :: molecularDensities
    type(radiationStructure),           intent(in)  :: radiation
    double precision,                   intent(out) :: coolingFunction
    double precision                                :: logarithmic10Temperature,coolingFunctionLowDensityLimit &
         &,coolingFunctionLocalThermodynamicEquilibrium ,molecularHydrogenDensity,atomicHydrogenDensity,electronDensity&
         &,molecularHydrogenCationDensity ,numberDensityCriticalOverNumberDensityHydrogen&
         &,coolingFunctionAtomicHydrogenMolecularHydrogenCation ,coolingFunctionElectronMolecularHydrogenCation

    ! Check if this cooling function has been selected.
    if (functionSelected) then

       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,molecularDensities)

       ! H_2^+ - e^- cooling function.
       coolingFunction=coolingFunction+Cooling_Function_GP_H2Plus_Electron(temperature,molecularDensities)

       ! H - H_2^+ cooling function.
       coolingFunction=coolingFunction+Cooling_Function_GP_H_H2Plus       (temperature,molecularDensities)

    else

       ! Not selected, return zero.
       coolingFunction=0.0d0

    end if

    return
  end subroutine Cooling_Function_Molecular_Hydrogen_GP
  
  !# <coolingFunctionDensitySlopeCompute>
  !#   <unitName>Cooling_Function_Density_Slope_Molecular_Hydrogen_GP</unitName>
  !# </coolingFunctionDensitySlopeCompute>
  subroutine Cooling_Function_Density_Slope_Molecular_Hydrogen_GP(coolingFunctionDensitySlope,temperature,numberDensityHydrogen&
       &,abundances,molecularDensities ,radiation)
    !% Return the gradient with respect to density of the cooling function due to molecular hydrogen using the cooling function
    !% of \cite{galli_chemistry_1998}.
    use Ionization_States
    use Abundances_Structure
    use Molecular_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),          intent(in)  :: abundances
    type(molecularAbundancesStructure), intent(in)  :: molecularDensities
    type(radiationStructure),           intent(in)  :: radiation
    double precision,                   intent(out) :: coolingFunctionDensitySlope
    double precision                                :: coolingFunction,numberDensityCriticalOverNumberDensityHydrogen&
         &,coolingFunctionLocalThermodynamicEquilibrium ,coolingFunctionLowDensityLimit

    ! Check if this cooling function has been selected.
    if (functionSelected) then

       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,molecularDensities)
       call Number_Density_Critical_Over_Number_Density_Hydrogen(numberDensityHydrogen,temperature &
            &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
            &,coolingFunctionLowDensityLimit)
       if (coolingFunctionLowDensityLimit > 0.0d0) then
          coolingFunctionDensitySlope=(coolingFunction/numberDensityHydrogen)*(1.0d0&
               &+numberDensityCriticalOverNumberDensityHydrogen/(1.0d0+numberDensityCriticalOverNumberDensityHydrogen))
       else
          coolingFunctionDensitySlope=(coolingFunction/numberDensityHydrogen)*2.0d0
       end if

       ! H_2^+ - e^- cooling function.
       coolingFunction=Cooling_Function_GP_H2Plus_Electron(temperature,molecularDensities)
       coolingFunctionDensitySlope=coolingFunctionDensitySlope+(coolingFunction/numberDensityHydrogen)*2.0d0

       ! H - H_2^+ cooling function.
       coolingFunction=Cooling_Function_GP_H_H2Plus(temperature,molecularDensities)
       coolingFunctionDensitySlope=coolingFunctionDensitySlope+(coolingFunction/numberDensityHydrogen)*2.0d0

    else
       
       ! Not selected, return zero.
       coolingFunctionDensitySlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Density_Slope_Molecular_Hydrogen_GP
  
  !# <coolingFunctionTemperatureSlopeCompute>
  !#   <unitName>Cooling_Function_Temperature_Slope_Molecular_Hydrogen_GP</unitName>
  !# </coolingFunctionTemperatureSlopeCompute>
  subroutine Cooling_Function_Temperature_Slope_Molecular_Hydrogen_GP(coolingFunctionTemperatureSlope,temperature &
       &,numberDensityHydrogen,abundances,molecularDensities,radiation)
    !% Return the gradient with respect to temperature of the cooling function due to molecular hydrogen using the cooling
    !% function of \cite{galli_chemistry_1998}.
    use Ionization_States
    use Abundances_Structure
    use Molecular_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Prefixes
    implicit none
    double precision,                   intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),          intent(in)  :: abundances
    type(molecularAbundancesStructure), intent(in)  :: molecularDensities
    type(radiationStructure),           intent(in)  :: radiation
    double precision,                   intent(out) :: coolingFunctionTemperatureSlope
    double precision,                   save        :: temperaturePrevious1=-1,temperaturePrevious2=-1,temperaturePrevious3=-1&
         &,temperaturePrevious4=-1,coolingFunctionLowDensityLimitTemperatureLogGradient&
         &,coolingFunctionH2PlusElectronTemperatureLogGradient,coolingFunctionH2PlusHTemperatureLogGradient
    !$omp threadprivate(temperaturePrevious1,temperaturePrevious2,temperaturePrevious3,temperaturePrevious4)
    !$omp threadprivate(coolingFunctionLowDensityLimitTemperatureLogGradient)
    !$omp threadprivate(coolingFunctionH2PlusElectronTemperatureLogGradient,coolingFunctionH2PlusHTemperatureLogGradient)
    double precision                                :: coolingFunction,numberDensityCriticalOverNumberDensityHydrogen&
         &,coolingFunctionLocalThermodynamicEquilibrium ,coolingFunctionLowDensityLimit&
         &,coolingFunctionLowDensityLimitTemperatureGradient&
         &,temperatureThousand,coolingFunctionRotationalTemperatureGradient,coolingFunctionVibrationalTemperatureGradient&
         &,coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient,logarithmic10Temperature

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,molecularDensities)
       call Number_Density_Critical_Over_Number_Density_Hydrogen(numberDensityHydrogen,temperature &
            &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
            &,coolingFunctionLowDensityLimit)
       ! Check if we need to recompute the low density limit cooling function gradient.
       if (temperature /= temperaturePrevious1) then
          logarithmic10Temperature=dlog10(temperature)
          coolingFunctionLowDensityLimitTemperatureLogGradient=                                                 &
               &                    +                          coolingFunctionLowDensityLimitCoefficient(1)     &
               &                    +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(2)     &
               &                    +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(3)     &
               &                    +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(4))))
          temperaturePrevious1=temperature
       end if
       if (coolingFunctionLowDensityLimit > 0.0d0) then
          ! Get the temperature gradient of the low density limit cooling function.
          coolingFunctionLowDensityLimitTemperatureGradient=(coolingFunctionLowDensityLimit/temperature)&
               &*coolingFunctionLowDensityLimitTemperatureLogGradient
          ! Check if we need to compute gradients of rotational and vibrational cooling functions.
          if (temperature /= temperaturePrevious2) then
             ! Get the temperature gradient of the LTE cooling function.
             temperatureThousand=temperature*milli ! Convert to units of 1,000K.
             coolingFunctionRotationalTemperatureGradient =                                                                                                  &
                  & +(  coolingFunctionRotationalLambda1*(temperatureThousand**(coolingFunctionRotationalExponent1-1.0d0))                                   &
                  &    /(1.0d0+coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2))                              &
                  &  )                                                                                                                                       &
                  &   *dexp(   -(coolingFunctionRotationalTemperature1/temperatureThousand)**3)*(                                                            &
                  &      +3.0d0*(coolingFunctionRotationalTemperature1/temperatureThousand)**3                                                               &
                  &      +coolingFunctionRotationalExponent1                                                                                                 &
                  &      -coolingFunctionRotationalExponent2*coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2) &
                  &         /(1.0d0+coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2))                         &
                  &                                                                             )                                                            &
                  & +(coolingFunctionRotationalLambda2/temperatureThousand)                                                                                  &
                  &   *    ( coolingFunctionRotationalTemperature2/temperatureThousand)                                                                      &
                  &   *dexp(-coolingFunctionRotationalTemperature2/temperatureThousand)
             coolingFunctionVibrationalTemperatureGradient=(                                                                                                 &
                  &                                     coolingFunctionVibrationalLambda1*(coolingFunctionVibrationalTemperature1/temperatureThousand)       &
                  &                                       *dexp(-coolingFunctionVibrationalTemperature1/temperatureThousand)                                 &
                  &                                    +coolingFunctionVibrationalLambda2*(coolingFunctionVibrationalTemperature2/temperatureThousand)       &
                  &                                       *dexp(-coolingFunctionVibrationalTemperature2/temperatureThousand)                                 &
                  &                                        )/temperatureThousand
             temperaturePrevious2=temperature
          end if
          coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient=milli*( coolingFunctionRotationalTemperatureGradient  &
               &                                                                 +coolingFunctionVibrationalTemperatureGradient &
               &                                                                )/numberDensityHydrogen
          ! Get the temperature gradient of the cooling function.
          coolingFunctionTemperatureSlope=coolingFunction*(                                                                   &
               &+coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient/coolingFunctionLocalThermodynamicEquilibrium &
               &-(1.0d0/coolingFunctionLowDensityLimit/(1.0d0+numberDensityCriticalOverNumberDensityHydrogen))*(              &
               &      +coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient                                        &
               &      -numberDensityCriticalOverNumberDensityHydrogen*coolingFunctionLowDensityLimitTemperatureGradient       &
               &                                                                                               )              &
               &                                          )
       else
          coolingFunctionTemperatureSlope=coolingFunction*coolingFunctionLowDensityLimitTemperatureLogGradient/temperature
       end if

       ! H_2^+ - e^- cooling function.
       if (temperature > 1.0d3 .and. temperature < 1.0d4) then
          coolingFunction=Cooling_Function_GP_H2Plus_Electron(temperature,molecularDensities)
          if (temperature /= temperaturePrevious3) then
             logarithmic10Temperature=dlog10(temperature)
             coolingFunctionH2PlusElectronTemperatureLogGradient=                         &
                  & +coolingFunctionH2PlusElectronCoefficient(1)                          &
                  & +coolingFunctionH2PlusElectronCoefficient(2)*logarithmic10Temperature
             temperaturePrevious3=temperature
          end if
          coolingFunctionTemperatureSlope= coolingFunctionTemperatureSlope &
               &                          +(coolingFunction/temperature)*coolingFunctionH2PlusElectronTemperatureLogGradient
       end if

       ! H - H_2^+ cooling function.
       if (temperature > 1.0d3 .and. temperature < 1.0d4) then
          coolingFunction=Cooling_Function_GP_H_H2Plus(temperature,molecularDensities)
          if (temperature /= temperaturePrevious4) then
             logarithmic10Temperature=dlog10(temperature)
             coolingFunctionH2PlusHTemperatureLogGradient=                         &
                  & +coolingFunctionH2PlusHCoefficient(1)                          &
                  & +coolingFunctionH2PlusHCoefficient(2)*logarithmic10Temperature
             temperaturePrevious4=temperature
          end if
          coolingFunctionTemperatureSlope= coolingFunctionTemperatureSlope &
               &                          +(coolingFunction/temperature)*coolingFunctionH2PlusHTemperatureLogGradient
       end if

    else
       
       ! Not selected, return zero.
       coolingFunctionTemperatureSlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Temperature_Slope_Molecular_Hydrogen_GP

  double precision function Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,molecularDensities)
    !% Compute the cooling function due to H--H$_2$ interactions.
    use Molecular_Abundances_Structure
    implicit none
    double precision,                   intent(in) :: numberDensityHydrogen,temperature
    type(molecularAbundancesStructure), intent(in) :: molecularDensities
    double precision                               :: atomicHydrogenDensity,molecularHydrogenDensity&
         &,numberDensityCriticalOverNumberDensityHydrogen ,coolingFunctionLocalThermodynamicEquilibrium &
         &,coolingFunctionLowDensityLimit

    ! Get the relevant densities.
    atomicHydrogenDensity   =molecularDensities%abundance(atomicHydrogenMoleculeIndex   )
    molecularHydrogenDensity=molecularDensities%abundance(molecularHydrogenMoleculeIndex)

    ! Compute the cooling function.
    call Number_Density_Critical_Over_Number_Density_Hydrogen(numberDensityHydrogen,temperature &
         &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
         &,coolingFunctionLowDensityLimit)
    if (coolingFunctionLowDensityLimit > 0.0d0) then
       Cooling_Function_GP_H_H2=coolingFunctionLocalThermodynamicEquilibrium/(1.0d0+numberDensityCriticalOverNumberDensityHydrogen)
    else
       Cooling_Function_GP_H_H2=coolingFunctionLowDensityLimit
    end if

    ! Scale to the density of hydrogen and molecular hydrogen.
    Cooling_Function_GP_H_H2=Cooling_Function_GP_H_H2*molecularHydrogenDensity*atomicHydrogenDensity
    return
  end function Cooling_Function_GP_H_H2

  subroutine Number_Density_Critical_Over_Number_Density_Hydrogen(numberDensityHydrogen,temperature&
       &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium&
       &,coolingFunctionLowDensityLimit)
    !% Compute the ratio of critical number density to the hydrogen number density for use in molecular hydrogen cooling functions.
    use Numerical_Constants_Prefixes
    implicit none
    double precision, intent(in)  :: numberDensityHydrogen,temperature
    double precision, intent(out) :: numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
         &,coolingFunctionLowDensityLimit
    double precision, save        :: temperaturePrevious=-1.0d0,coolingFunctionRotationalTemperaturePart&
         &,coolingFunctionVibrationalTemperaturePart,coolingFunctionLowDensityLimitStored
    !$omp threadprivate(temperaturePrevious,coolingFunctionRotationalTemperaturePart,coolingFunctionVibrationalTemperaturePart)
    !$omp threadprivate(coolingFunctionLowDensityLimitStored)
    double precision              :: temperatureThousand,coolingFunctionRotational,coolingFunctionVibrational,logarithmic10Temperature

    if (temperature /= temperaturePrevious) then
       ! The expression from Galli & Palla (1998), assumes an equilibrium (1:3) ratio of para:ortho.
       logarithmic10Temperature=dlog10(temperature)
       coolingFunctionLowDensityLimitStored=10.0d0**(coolingFunctionLowDensityLimitCoefficient(0)     &
            &             +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(1)     &
            &             +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(2)     &
            &             +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(3)     &
            &             +logarithmic10Temperature*(coolingFunctionLowDensityLimitCoefficient(4))))) &
            &                                 )
       ! Rotational and vibrational cooling functions from Hollenbach & McKee (1979; their equations 6.37 and 6.38).
       temperatureThousand=temperature*milli ! Convert to units of 1,000K.
       coolingFunctionRotationalTemperaturePart =((coolingFunctionRotationalLambda1*(temperatureThousand**coolingFunctionRotationalExponent1)) &
            & /(1.0d0+coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2)))                        &
            & *dexp(-(coolingFunctionRotationalTemperature1/temperatureThousand)**3)                                                           &
            &                                     +coolingFunctionRotationalLambda2                                                            &
            & *dexp(-coolingFunctionRotationalTemperature2/temperatureThousand)
       coolingFunctionVibrationalTemperaturePart= coolingFunctionVibrationalLambda1                                    &
            &                                       *dexp(-coolingFunctionVibrationalTemperature1/temperatureThousand) &
            &                                    +coolingFunctionVibrationalLambda2                                    &
            &                                       *dexp(-coolingFunctionVibrationalTemperature2/temperatureThousand)
       
       ! Record the temperature used.
       temperaturePrevious=temperature
    end if

    ! Get the low density limit cooling function.
    coolingFunctionLowDensityLimit=coolingFunctionLowDensityLimitStored

    ! Compute the LTE cooling function.
    coolingFunctionLocalThermodynamicEquilibrium=(coolingFunctionRotationalTemperaturePart&
         &+coolingFunctionVibrationalTemperaturePart)/numberDensityHydrogen

    ! Compute the number density ratio.
    if (coolingFunctionLowDensityLimit > 0.0d0) then
       numberDensityCriticalOverNumberDensityHydrogen=coolingFunctionLocalThermodynamicEquilibrium&
            &/coolingFunctionLowDensityLimit
    else
       ! Unphysical value, should not be used.
       numberDensityCriticalOverNumberDensityHydrogen=-1.0d0
    end if
    return
  end subroutine Number_Density_Critical_Over_Number_Density_Hydrogen

  double precision function Cooling_Function_GP_H2Plus_Electron(temperature,molecularDensities)
    !% Compute the cooling function due to H$_2^+$--e$^-$ interactions.
    use Molecular_Abundances_Structure
    implicit none
    double precision,                   intent(in) :: temperature
    type(molecularAbundancesStructure), intent(in) :: molecularDensities
    double precision,                   save       :: temperaturePrevious,coolingFunctionElectronMolecularHydrogenCation
    !$omp threadprivate(temperaturePrevious,coolingFunctionElectronMolecularHydrogenCation)
    double precision                               :: electronDensity,molecularHydrogenCationDensity ,logarithmic10Temperature
    
    ! Get the relevant densities.
    electronDensity               =molecularDensities%abundance(electronMoleculeIndex               )
    molecularHydrogenCationDensity=molecularDensities%abundance(molecularHydrogenCationMoleculeIndex)
    
    ! H_2^+ - e^- cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= temperaturePrevious) then
          logarithmic10Temperature=dlog10(temperature)
          coolingFunctionElectronMolecularHydrogenCation=10.0d0**(                        &
               &  coolingFunctionH2PlusElectronCoefficient(0)                             &
               & +coolingFunctionH2PlusElectronCoefficient(1)*logarithmic10Temperature    &
               & +coolingFunctionH2PlusElectronCoefficient(2)*logarithmic10Temperature**2 &
               &                                                 )
          temperaturePrevious=temperature
       end if
       Cooling_Function_GP_H2Plus_Electron=coolingFunctionElectronMolecularHydrogenCation*molecularHydrogenCationDensity &
            &*electronDensity
    else
       Cooling_Function_GP_H2Plus_Electron=0.0d0
    end if
    
    return
  end function Cooling_Function_GP_H2Plus_Electron

  double precision function Cooling_Function_GP_H_H2Plus(temperature,molecularDensities)
    !% Compute the cooling function due to H--H$_2^+$ interactions.
    use Molecular_Abundances_Structure
    implicit none
    double precision,                   intent(in) :: temperature
    type(molecularAbundancesStructure), intent(in) :: molecularDensities
    double precision,                   save       :: temperaturePrevious,coolingFunctionAtomicHydrogenMolecularHydrogenCation
    !$omp threadprivate(temperaturePrevious,coolingFunctionAtomicHydrogenMolecularHydrogenCation)
    double precision                               :: atomicHydrogenDensity,molecularHydrogenCationDensity,logarithmic10Temperature
    
    ! Get the relevant densities.
    atomicHydrogenDensity         =molecularDensities%abundance(atomicHydrogenMoleculeIndex         )
    molecularHydrogenCationDensity=molecularDensities%abundance(molecularHydrogenCationMoleculeIndex)
    
    ! H - H_2^+ cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= temperaturePrevious) then
          logarithmic10Temperature=dlog10(temperature)
          coolingFunctionAtomicHydrogenMolecularHydrogenCation=10.0d0**(           &
               &  coolingFunctionH2PlusHCoefficient(0)                             &
               & +coolingFunctionH2PlusHCoefficient(1)*logarithmic10Temperature    &
               & +coolingFunctionH2PlusHCoefficient(2)*logarithmic10Temperature**2 &
               &                                                       )
          temperaturePrevious=temperature
       end if
       Cooling_Function_GP_H_H2Plus=coolingFunctionAtomicHydrogenMolecularHydrogenCation*atomicHydrogenDensity&
            &*molecularHydrogenCationDensity
    else
       Cooling_Function_GP_H_H2Plus=0.0d0
    end if
    
    return
  end function Cooling_Function_GP_H_H2Plus

end module Cooling_Functions_Molecular_Hydrogen_Galli_Palla
