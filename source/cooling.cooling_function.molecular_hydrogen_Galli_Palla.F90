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

!% Contains a module which computes the contribution to the cooling function from molecular hydrogen using the cooling function of
!% \cite{galli_chemistry_1998}.

module Cooling_Functions_Molecular_Hydrogen_Galli_Palla
  !% Computes the contribution to the cooling function from molecular hydrogen using the cooling function of
  !\cite{galli_chemistry_1998}.
  use ISO_Varying_String
  implicit none
  private
  public :: Cooling_Function_Molecular_Hydrogen_GP_Initialize, Cooling_Function_Molecular_Hydrogen_GP,&
       & Cooling_Function_Density_Slope_Molecular_Hydrogen_GP, Cooling_Function_Temperature_Slope_Molecular_Hydrogen_GP
  
  ! Flag indicating whether or not this cooling function is selected.
  logical                                     :: functionSelected                         =.false.                                                                                     
  
  ! Indices of "chemical species" (includes atoms, atomic ions and electrons also) used in this cooling function.
  integer                                     :: atomicHydrogenChemicalIndex                                                                       , electronChemicalIndex         , & 
       &                                         molecularHydrogenCationChemicalIndex                                                              , molecularHydrogenChemicalIndex    
  
  ! Parameters for Hollenbach & McKee cooling function fits.
  double precision                , parameter :: coolingFunctionRotationalLambda1         =9.50d-22                                                                                    
  double precision                , parameter :: coolingFunctionRotationalLambda2         =3.00d-24                                                                                    
  double precision                , parameter :: coolingFunctionRotationalTemperature1    =0.13d0                                                                                      
  double precision                , parameter :: coolingFunctionRotationalTemperature2    =0.51d0                                                                                      
  double precision                , parameter :: coolingFunctionRotationalExponent1       =3.76d0                                                                                      
  double precision                , parameter :: coolingFunctionRotationalExponent2       =2.10d0                                                                                      
  double precision                , parameter :: coolingFunctionRotationalCoefficient1    =0.12d0                                                                                      
  double precision                , parameter :: coolingFunctionVibrationalLambda1        =6.70d-19                                                                                    
  double precision                , parameter :: coolingFunctionVibrationalLambda2        =1.60d-18                                                                                    
  double precision                , parameter :: coolingFunctionVibrationalTemperature1   =5.86d0                                                                                      
  double precision                , parameter :: coolingFunctionVibrationalTemperature2   =11.70d0                                                                                     
  
  ! Parameters for low-density limit cooling function.
  double precision, dimension(0:4), parameter :: coolingFunctionLowDensityLimitCoefficient=[-103.0000d0,+97.5900d0,-48.0500d0,+10.8000d0,-0.9032d0]                                    
  
  ! Parameters for H_2^+ - e^- cooling function.
  double precision, dimension(0:2), parameter :: coolingFunctionH2PlusElectronCoefficient =[-33.3299d0,+5.56465d0,-4.67461d-1]                                                         
  
  ! Parameters for H_2^+ - H cooling function.
  double precision, dimension(0:2), parameter :: coolingFunctionH2PlusHCoefficient        =[-35.2804d0,+5.86234d0,-5.12276d-1]                                                         
  
contains
  
  !# <coolingFunctionMethods>
  !#  <unitName>Cooling_Function_Molecular_Hydrogen_GP_Initialize</unitName>
  !#  <methodName>molecularHydrogenGalliPalla</methodName>
  !# </coolingFunctionMethods>
  subroutine Cooling_Function_Molecular_Hydrogen_GP_Initialize(coolingFunctionMethods,coolingFunctionsMatched)
    !% Initializes the ``molecular hydrogen Galli \& Palla'' cooling function module.
    use Chemical_Abundances_Structure
    implicit none
    type   (varying_string), intent(in   ) :: coolingFunctionMethods (:) 
    integer                , intent(inout) :: coolingFunctionsMatched    
    
    ! Check if this cooling function has been selected.
    if (any(coolingFunctionMethods == 'molecularHydrogenGalliPalla')) then
       ! Flag that this cooling function has been selected.
       functionSelected=.true.
       coolingFunctionsMatched=coolingFunctionsMatched+1
       
       ! Get the indices of chemicals that will be used.
       electronChemicalIndex               =Chemicals_Index("Electron"               )
       atomicHydrogenChemicalIndex         =Chemicals_Index("AtomicHydrogen"         )
       molecularHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
       molecularHydrogenChemicalIndex      =Chemicals_Index("MolecularHydrogen"      )
    end if

    return
  end subroutine Cooling_Function_Molecular_Hydrogen_GP_Initialize

  !# <coolingFunctionCompute>
  !#   <unitName>Cooling_Function_Molecular_Hydrogen_GP</unitName>
  !# </coolingFunctionCompute>
  subroutine Cooling_Function_Molecular_Hydrogen_GP(coolingFunction,temperature,numberDensityHydrogen,gasAbundances&
       &,chemicalDensities,radiation)
    !% Return the cooling function due to molecular hydrogen using the cooling function of \cite{galli_chemistry_1998} (which
    !% refers to the local thermodynamic equilibrium cooling function of \cite{hollenbach_molecule_1979}). Cooling functions
    !% involving H$_2^+$ are computed using polynomial fits to the results of \cite{suchkov_cooling_1978} found by Andrew
    !% Benson by measuring curves from the original paper.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature 
    type            (abundances        ), intent(in   ) :: gasAbundances                      
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                  
    type            (radiationStructure), intent(in   ) :: radiation                          
    double precision                    , intent(  out) :: coolingFunction                    
    
    ! Check if this cooling function has been selected and the hydrogen density is positive.
    if (functionSelected.and.numberDensityHydrogen > 0.0d0) then

       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,chemicalDensities)

       ! H_2^+ - e^- cooling function.
       coolingFunction=coolingFunction+Cooling_Function_GP_H2Plus_Electron(temperature,chemicalDensities)

       ! H - H_2^+ cooling function.
       coolingFunction=coolingFunction+Cooling_Function_GP_H_H2Plus       (temperature,chemicalDensities)

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
       &,gasAbundances,chemicalDensities ,radiation)
    !% Return the gradient with respect to density of the cooling function due to molecular hydrogen using the cooling function
    !% of \cite{galli_chemistry_1998}.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen         , temperature                                       
    type            (abundances        ), intent(in   ) :: gasAbundances                                                                     
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                                                                 
    type            (radiationStructure), intent(in   ) :: radiation                                                                         
    double precision                    , intent(  out) :: coolingFunctionDensitySlope                                                       
    double precision                                    :: coolingFunction               , coolingFunctionLocalThermodynamicEquilibrium  , & 
         &                                                 coolingFunctionLowDensityLimit, numberDensityCriticalOverNumberDensityHydrogen    
    
    ! Check if this cooling function has been selected and the hydrogen density is positive.
    if (functionSelected.and.numberDensityHydrogen > 0.0d0) then

       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,chemicalDensities)
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
       coolingFunction=Cooling_Function_GP_H2Plus_Electron(temperature,chemicalDensities)
       coolingFunctionDensitySlope=coolingFunctionDensitySlope+(coolingFunction/numberDensityHydrogen)*2.0d0

       ! H - H_2^+ cooling function.
       coolingFunction=Cooling_Function_GP_H_H2Plus(temperature,chemicalDensities)
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
       &,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)
    !% Return the gradient with respect to temperature of the cooling function due to molecular hydrogen using the cooling
    !% function of \cite{galli_chemistry_1998}.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Prefixes
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen                                             , temperature                                        
    type            (abundances        ), intent(in   ) :: gasAbundances                                                                                                          
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                                                                                                      
    type            (radiationStructure), intent(in   ) :: radiation                                                                                                              
    double precision                    , intent(  out) :: coolingFunctionTemperatureSlope                                                                                        
    double precision                    , save          :: coolingFunctionH2PlusElectronTemperatureLogGradient               , coolingFunctionH2PlusHTemperatureLogGradient   , & 
         &                                                 coolingFunctionLowDensityLimitTemperatureLogGradient              , temperaturePrevious1                        =-1, & 
         &                                                 temperaturePrevious2                                           =-1, temperaturePrevious3                        =-1, & 
         &                                                 temperaturePrevious4                                           =-1                                                     
    !$omp threadprivate(temperaturePrevious1,temperaturePrevious2,temperaturePrevious3,temperaturePrevious4)
    !$omp threadprivate(coolingFunctionLowDensityLimitTemperatureLogGradient)
    !$omp threadprivate(coolingFunctionH2PlusElectronTemperatureLogGradient,coolingFunctionH2PlusHTemperatureLogGradient)
    double precision                                    :: coolingFunction                                                   , coolingFunctionLocalThermodynamicEquilibrium   , & 
         &                                                 coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient   , coolingFunctionLowDensityLimit                 , & 
         &                                                 coolingFunctionLowDensityLimitTemperatureGradient                 , coolingFunctionRotationalTemperatureGradient   , & 
         &                                                 coolingFunctionVibrationalTemperatureGradient                     , logarithmic10Temperature                       , & 
         &                                                 numberDensityCriticalOverNumberDensityHydrogen                    , temperatureThousand                                
    
    ! Check if this cooling function has been selected and the hydrogen density is positive.
    if (functionSelected.and.numberDensityHydrogen > 0.0d0) then
       ! H - H_2 cooling function.
       coolingFunction=Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,chemicalDensities)
       call Number_Density_Critical_Over_Number_Density_Hydrogen(numberDensityHydrogen,temperature &
            &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
            &,coolingFunctionLowDensityLimit)
       ! Check if we need to recompute the low density limit cooling function gradient.
       if (temperature /= temperaturePrevious1) then
          logarithmic10Temperature=log10(temperature)
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
                  &   *exp(   -(coolingFunctionRotationalTemperature1/temperatureThousand)**3)*(                                                            &
                  &      +3.0d0*(coolingFunctionRotationalTemperature1/temperatureThousand)**3                                                               &
                  &      +coolingFunctionRotationalExponent1                                                                                                 &
                  &      -coolingFunctionRotationalExponent2*coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2) &
                  &         /(1.0d0+coolingFunctionRotationalCoefficient1*(temperatureThousand**coolingFunctionRotationalExponent2))                         &
                  &                                                                             )                                                            &
                  & +(coolingFunctionRotationalLambda2/temperatureThousand)                                                                                  &
                  &   *    ( coolingFunctionRotationalTemperature2/temperatureThousand)                                                                      &
                  &   *exp(-coolingFunctionRotationalTemperature2/temperatureThousand)
             coolingFunctionVibrationalTemperatureGradient=(                                                                                                 &
                  &                                     coolingFunctionVibrationalLambda1*(coolingFunctionVibrationalTemperature1/temperatureThousand)       &
                  &                                       *exp(-coolingFunctionVibrationalTemperature1/temperatureThousand)                                 &
                  &                                    +coolingFunctionVibrationalLambda2*(coolingFunctionVibrationalTemperature2/temperatureThousand)       &
                  &                                       *exp(-coolingFunctionVibrationalTemperature2/temperatureThousand)                                 &
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
          coolingFunction=Cooling_Function_GP_H2Plus_Electron(temperature,chemicalDensities)
          if (temperature /= temperaturePrevious3) then
             logarithmic10Temperature=log10(temperature)
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
          coolingFunction=Cooling_Function_GP_H_H2Plus(temperature,chemicalDensities)
          if (temperature /= temperaturePrevious4) then
             logarithmic10Temperature=log10(temperature)
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

  double precision function Cooling_Function_GP_H_H2(numberDensityHydrogen,temperature,chemicalDensities)
    !% Compute the cooling function due to H--H$_2$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen                         , temperature                                     
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                                                                               
    double precision                                    :: atomicHydrogenDensity                         , coolingFunctionLocalThermodynamicEquilibrium, & 
         &                                                 coolingFunctionLowDensityLimit                , molecularHydrogenDensity                    , & 
         &                                                 numberDensityCriticalOverNumberDensityHydrogen                                                  
    
    ! Get the relevant densities.
    atomicHydrogenDensity   =chemicalDensities%abundance(atomicHydrogenChemicalIndex   )
    molecularHydrogenDensity=chemicalDensities%abundance(molecularHydrogenChemicalIndex)

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
    double precision, intent(in   ) :: numberDensityHydrogen                         , temperature                                        
    double precision, intent(  out) :: coolingFunctionLocalThermodynamicEquilibrium  , coolingFunctionLowDensityLimit                 , & 
         &                             numberDensityCriticalOverNumberDensityHydrogen                                                     
    double precision, save          :: coolingFunctionLowDensityLimitStored          , coolingFunctionRotationalTemperaturePart       , & 
         &                             coolingFunctionVibrationalTemperaturePart     , temperaturePrevious                     =-1.0d0    
    !$omp threadprivate(temperaturePrevious,coolingFunctionRotationalTemperaturePart,coolingFunctionVibrationalTemperaturePart)
    !$omp threadprivate(coolingFunctionLowDensityLimitStored)
    double precision                :: logarithmic10Temperature                      , temperatureThousand                                
    
    if (temperature /= temperaturePrevious) then
       ! The expression from Galli & Palla (1998), assumes an equilibrium (1:3) ratio of para:ortho.
       logarithmic10Temperature=log10(temperature)
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
            & *exp(-(coolingFunctionRotationalTemperature1/temperatureThousand)**3)                                                           &
            &                                     +coolingFunctionRotationalLambda2                                                            &
            & *exp(-coolingFunctionRotationalTemperature2/temperatureThousand)
       coolingFunctionVibrationalTemperaturePart= coolingFunctionVibrationalLambda1                                    &
            &                                       *exp(-coolingFunctionVibrationalTemperature1/temperatureThousand) &
            &                                    +coolingFunctionVibrationalLambda2                                    &
            &                                       *exp(-coolingFunctionVibrationalTemperature2/temperatureThousand)
       
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

  double precision function Cooling_Function_GP_H2Plus_Electron(temperature,chemicalDensities)
    !% Compute the cooling function due to H$_2^+$--e$^-$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    double precision                    , intent(in   ) :: temperature                                                                 
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                                                           
    double precision                    , save          :: coolingFunctionElectronMolecularHydrogenCation, temperaturePrevious         
    !$omp threadprivate(temperaturePrevious,coolingFunctionElectronMolecularHydrogenCation)
    double precision                                    :: electronDensity                               , logarithmic10Temperature, & 
         &                                                 molecularHydrogenCationDensity                                              
    
    ! Get the relevant densities.
    electronDensity               =chemicalDensities%abundance(electronChemicalIndex               )
    molecularHydrogenCationDensity=chemicalDensities%abundance(molecularHydrogenCationChemicalIndex)
    
    ! H_2^+ - e^- cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= temperaturePrevious) then
          logarithmic10Temperature=log10(temperature)
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

  double precision function Cooling_Function_GP_H_H2Plus(temperature,chemicalDensities)
    !% Compute the cooling function due to H--H$_2^+$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    double precision                    , intent(in   ) :: temperature                                                                       
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities                                                                 
    double precision                    , save          :: coolingFunctionAtomicHydrogenMolecularHydrogenCation, temperaturePrevious         
    !$omp threadprivate(temperaturePrevious,coolingFunctionAtomicHydrogenMolecularHydrogenCation)
    double precision                                    :: atomicHydrogenDensity                               , logarithmic10Temperature, & 
         &                                                 molecularHydrogenCationDensity                                                    
    
    ! Get the relevant densities.
    atomicHydrogenDensity         =chemicalDensities%abundance(atomicHydrogenChemicalIndex         )
    molecularHydrogenCationDensity=chemicalDensities%abundance(molecularHydrogenCationChemicalIndex)
    
    ! H - H_2^+ cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= temperaturePrevious) then
          logarithmic10Temperature=log10(temperature)
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
