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


!% Contains a module which implements calculations of molecular reaction rates for hydrogen using the fits from
!% \cite{abel_modeling_1997} and \cite{tegmark_small_1997}.

module Molecular_Hydrogen_Rates
  !% Implements calculations of molecular reaction rates for hydrogen using the fits from \cite{abel_modeling_1997} and
  !% \cite{tegmark_small_1997}.
  use Molecular_Structures
  use Molecular_Abundances_Structure
  private
  public :: Molecular_Hydrogen_Rates_Initialize, Molecular_Hydrogen_Rates_Compute

  ! Flag indicating if these rates have been selected.
  logical          :: ratesSelected=.false.

  ! Flag indicating whether fast rate calculations should be used.
  logical          :: hydrogenNetworkFast

  ! Flag indicating whether to use CMB only when computing rates for some radiative processes.
  logical          :: hydrogenNetworkCMBOnly

  ! Indices of molecules.
  integer          :: atomicHydrogenIndex,atomicHydrogenCationIndex,electronIndex,atomicHydrogenAnionIndex

  ! Density of atomic hydrogen anion to use in all rate calculations.
  double precision :: densityAtomicHydrogenAnion
  !$omp threadprivate(densityAtomicHydrogenAnion)

contains

  !# <molecularReactionRates>
  !#  <unitName>Molecular_Hydrogen_Rates_Initialize</unitName>
  !# </molecularReactionRates>
  subroutine Molecular_Hydrogen_Rates_Initialize(molecularReactionRatesMethods)
    !% Initializes the molecular hydrogen reaction network module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string), intent(in) :: molecularReactionRatesMethods(:)
 
    ! Check if this cooling function has been selected.
    if (any(molecularReactionRatesMethods == 'hydrogenNetwork')) then
       ! Flag that these rates have been selected.
       ratesSelected=.true.

       ! Get parameters controlling implementation.
       !@ <inputParameter>
       !@   <name>hydrogenNetworkFast</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not to use simplifying assumptions to speed the hydrogen network calculation. If true, H$^-$
       !@    is assumed to be at equilibrium abundance, H$_2^+$ reactions are ignored and other slow reactions are ignored (see
       !@    \citealt{abel_modeling_1997}).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hydrogenNetworkFast',hydrogenNetworkFast,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>hydrogenNetworkCMBOnly</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not to use the cosmic microwave background only when computed certain radiative rates.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hydrogenNetworkCMBOnly',hydrogenNetworkCMBOnly,defaultValue=.true.)

       ! Get indices for molecules as necessary.
       if (hydrogenNetworkFast) then
          ! Get indices of molecules needed for equilbrium hydrogen anion calculation.
          atomicHydrogenIndex      =Molecules_Index("AtomicHydrogen"      )
          atomicHydrogenCationIndex=Molecules_Index("AtomicHydrogenCation")
          electronIndex            =Molecules_Index("Electron"            )
          if (atomicHydrogenIndex       <= 0) call Galacticus_Error_Report('Molecular_Hydrogen_Rates_Initialize','atomic&
               & hydrogen must be included for fast hydrogen network calculation')
          if (atomicHydrogenCationIndex <= 0) call Galacticus_Error_Report('Molecular_Hydrogen_Rates_Initialize'&
               &,'hydrogen cation must be included for fast hydrogen network calculation')
          if (electronIndex             <= 0) call Galacticus_Error_Report('Molecular_Hydrogen_Rates_Initialize'&
               &,'electrons must be included for fast hydrogen network calculation'      )
       else
          ! Get actual hydrogen anion index.
          atomicHydrogenAnionIndex =Molecules_Index("AtomicHydrogenAnion")
       end if
       
    end if

    return
  end subroutine Molecular_Hydrogen_Rates_Initialize

  !# <molecularRatesCompute>
  !#  <unitName>Molecular_Hydrogen_Rates_Compute</unitName>
  !# </molecularRatesCompute>
  subroutine Molecular_Hydrogen_Rates_Compute(temperature,moleculeDensity,radiation,moleculeRates)
    !% Compute rates of change of molecular abundances due to reactions involving molecular hydrogen species.
    use Abundances_Structure
    use Radiation_Structure
    use Galacticus_Error
    implicit none
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    double precision                                  :: creationTerm,destructionTerm

    ! Return if not selected.
    if (.not.ratesSelected) return

    ! Determine the atomic hydrogen anion density to use.
    if (hydrogenNetworkFast) then
       ! For the fast network, assume the hydrogen anion is always at equilibrium density. (Following eqn. 24 of Abel et al. 1997.)
       creationTerm   = H_Electron_to_Hminus_Photon_Rate_Coefficient   (temperature)*moleculeDensity%abundance(atomicHydrogenIndex      ) &
            &                                                                       *moleculeDensity%abundance(electronIndex            )
       destructionTerm= H_Hminus_to_H2_Electron_Rate_Coefficient       (temperature)*moleculeDensity%abundance(atomicHydrogenIndex      ) &
            &          +Hminus_Hplus_to_2H_Rate_Coefficient            (temperature)*moleculeDensity%abundance(atomicHydrogenCationIndex) &
            &          +Hminus_Electron_to_H_2Electron_Rate_Coefficient(temperature)*moleculeDensity%abundance(electronIndex            )
       if (destructionTerm /= 0.0d0) then
          densityAtomicHydrogenAnion=creationTerm/destructionTerm
       else
          if (creationTerm > 0.0d0) call Galacticus_Error_Report('Molecular_Hydrogen_Rates_Compute','hydrogen anion equilibrium density is infinite')
          densityAtomicHydrogenAnion=0.0d0
       end if
    else if (atomicHydrogenAnionIndex > 0) then
       ! For the slow network, if we have the hydrogen anion then use its density directly.
       densityAtomicHydrogenAnion=moleculeDensity%abundance(atomicHydrogenAnionIndex)
    else
       ! Otherwise, we have no way to compute the hydrogen anion density, so set it to zero.
       densityAtomicHydrogenAnion=0.0d0
    end if

    ! Compute rates. References after each call refer to the rate coefficient in Tegmark et al. (1997) and the equation number in
    ! Abel et al. (1997) respectively.
    call Molecular_Hydrogen_Rate_H_Electron_to_Hplus_2Electron  (temperature,radiation,moleculeDensity,moleculeRates) !    ;  1
    call Molecular_Hydrogen_Rate_Hplus_Electron_to_H_Photon     (temperature,radiation,moleculeDensity,moleculeRates) ! k_1;  2
    call Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon    (temperature,radiation,moleculeDensity,moleculeRates) ! k_2;  7
    call Molecular_Hydrogen_Rate_H_Hminus_to_H2_Electron        (temperature,radiation,moleculeDensity,moleculeRates) ! k_3;  8
    call Molecular_Hydrogen_Rate_H_Hplus_to_H2plus_Photon       (temperature,radiation,moleculeDensity,moleculeRates) ! k_5;  9
    call Molecular_Hydrogen_Rate_H2plus_H_to_H2_Hplus           (temperature,radiation,moleculeDensity,moleculeRates) ! k_6; 10
    call Molecular_Hydrogen_Rate_H2_Hplus_to_H2plus_H           (temperature,radiation,moleculeDensity,moleculeRates) !    ; 11
    call Molecular_Hydrogen_Rate_H2_Electron_to_2H_Electron     (temperature,radiation,moleculeDensity,moleculeRates) !    ; 12
    call Molecular_Hydrogen_Rate_H2_H_to_3H                     (temperature,radiation,moleculeDensity,moleculeRates) !    ; 13
    call Molecular_Hydrogen_Rate_Hminus_Electron_to_H_2Electron (temperature,radiation,moleculeDensity,moleculeRates) !    ; 14
    call Molecular_Hydrogen_Rate_Hminus_H_to_2H_Electron        (temperature,radiation,moleculeDensity,moleculeRates) !    ; 15
    call Molecular_Hydrogen_Rate_Hminus_Hplus_to_2H             (temperature,radiation,moleculeDensity,moleculeRates) !    ; 16
    call Molecular_Hydrogen_Rate_Hminus_Hplus_to_H2plus_Electron(temperature,radiation,moleculeDensity,moleculeRates) !    ; 17
    call Molecular_Hydrogen_Rate_H2plus_Electron_to_2H          (temperature,radiation,moleculeDensity,moleculeRates) !    ; 18
    call Molecular_Hydrogen_Rate_H2plus_Hminus_to_H2_H          (temperature,radiation,moleculeDensity,moleculeRates) !    ; 19
    call Molecular_Hydrogen_Rate_H_Gamma_to_Hplus_Electron      (temperature,radiation,moleculeDensity,moleculeRates) !    ; 20
    call Molecular_Hydrogen_Rate_Hminus_Gamma_to_H_Electron     (temperature,radiation,moleculeDensity,moleculeRates) ! k_4; 23
    call Molecular_Hydrogen_Rate_H2_Gamma_to_H2plus_Electron    (temperature,radiation,moleculeDensity,moleculeRates) !    ; 24
    call Molecular_Hydrogen_Rate_H2plus_Gamma_to_H_Hplus        (temperature,radiation,moleculeDensity,moleculeRates) ! k_7; 25
    call Molecular_Hydrogen_Rate_H2plus_Gamma_to_2Hplus_Electron(temperature,radiation,moleculeDensity,moleculeRates) !    ; 26
    call Molecular_Hydrogen_Rate_H2_Gamma_to_H2star_to_2H       (temperature,radiation,moleculeDensity,moleculeRates) !    ; 27
    call Molecular_Hydrogen_Rate_H2_Gamma_to_2H                 (temperature,radiation,moleculeDensity,moleculeRates) !    ; 28
    return
  end subroutine Molecular_Hydrogen_Rates_Compute

  subroutine Molecular_Hydrogen_Rate_H_Electron_to_Hplus_2Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^+ + 2\hbox{e}^-$.
    use Radiation_Structure
    use Atomic_Rates_Ionization_Collisional
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hplus_2Electron_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex      =Molecules_Index("AtomicHydrogen"      )
       atomicHydrogenCationMoleculeIndex=Molecules_Index("AtomicHydrogenCation")
       electronMoleculeIndex            =Molecules_Index("Electron"            )
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0 .and. electronMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hplus_2Electron_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=Atomic_Rate_Ionization_Collisional(1,1,temperature)

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(electronMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex      ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex            , &
            & moleculeRates%abundance   (electronMoleculeIndex            ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H_Electron_to_Hplus_2Electron

  subroutine Molecular_Hydrogen_Rate_Hplus_Electron_to_H_Photon(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^+ + \hbox{e}^- \rightarrow \hbox{H} + \gamma$.
    use Radiation_Structure
    use Atomic_Rates_Recombination_Radiative
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_Hplus_Electron_to_H_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex      =Molecules_Index("AtomicHydrogen"      )
       atomicHydrogenCationMoleculeIndex=Molecules_Index("AtomicHydrogenCation")
       electronMoleculeIndex            =Molecules_Index("Electron"            )
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0 .and. electronMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_Hplus_Electron_to_H_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=Atomic_Rate_Recombination_Radiative(1,1,temperature)

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(electronMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex            , &
            & moleculeRates%abundance   (electronMoleculeIndex            ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex      ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hplus_Electron_to_H_Photon

  subroutine Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^- + \gamma$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenAnionMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex     =Molecules_Index("AtomicHydrogen"     )
       atomicHydrogenAnionMoleculeIndex=Molecules_Index("AtomicHydrogenAnion")
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0 .and. electronMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Get the rate coefficient.
       rateCoefficient=H_Electron_to_Hminus_Photon_Rate_Coefficient(temperature)

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(electronMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex     , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex     ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex           , &
            & moleculeRates%abundance   (electronMoleculeIndex           ) &
            & -rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon

  double precision function H_Electron_to_Hminus_Photon_Rate_Coefficient(temperature)
    !% Computes the rate coefficient (in units of cm$^3$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^- + \gamma$.
    implicit none
    double precision, intent(in) :: temperature
    double precision, save       :: temperaturePrevious,rateCoefficientStored
    !$omp threadprivate(temperaturePrevious,rateCoefficientStored)
    double precision             :: log10Temperature

    ! Determine if we need to recompute the rate coefficient.
    if (temperature /= temperaturePrevious) then
       ! Store the new temperature.
       temperaturePrevious=temperature

       ! Compute base 10 logarithm of temperature.
       log10Temperature=dlog10(temperature)
       
       ! Compute rate coefficient.
       if      (temperature <=    1.0d0) then
          rateCoefficientStored=1.429d-18
       else if (temperature <= 6000.0d0) then
          rateCoefficientStored=1.429d-18*(temperature**0.7620d0)*(temperature**(0.1523d0*log10Temperature))*(temperature**(-3.274d-2&
               &*(log10Temperature**2)))
       else
          rateCoefficientStored=3.802d-17*(temperature**(0.1998d0*log10Temperature))*(10.0d0**((4.0415d-5*(log10Temperature**2)-5.447d-3)*(log10Temperature**4)))
       end if
    end if
    
    ! Return the store rate coefficient.
    H_Electron_to_Hminus_Photon_Rate_Coefficient=rateCoefficientStored
    return
  end function H_Electron_to_Hminus_Photon_Rate_Coefficient

  subroutine Molecular_Hydrogen_Rate_H_Hminus_to_H2_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{e}^-$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,molecularHydrogenMoleculeIndex&
         &,electronMoleculeIndex,atomicHydrogenAnionMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Hminus_to_H2_Electron_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex     =Molecules_Index("AtomicHydrogen"     )
       atomicHydrogenAnionMoleculeIndex=Molecules_Index("AtomicHydrogenAnion")
       molecularHydrogenMoleculeIndex  =Molecules_Index("MolecularHydrogen"  )
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if all species were found.
       reactionActive=       atomicHydrogenMoleculeIndex      > 0                           &
            &         .and.  molecularHydrogenMoleculeIndex   > 0                           & 
            &         .and.  electronMoleculeIndex            > 0                           &
            &         .and. (atomicHydrogenAnionMoleculeIndex > 0 .or. hydrogenNetworkFast)
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Hminus_to_H2_Electron_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Get the rate coefficient.
       rateCoefficient=H_Hminus_to_H2_Electron_Rate_Coefficient(temperature)

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)*densityAtomicHydrogenAnion

       ! Record rate.
       if (.not.hydrogenNetworkFast) &
    &  call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex     , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex     ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex  , &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex  ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex           , &
            & moleculeRates%abundance   (electronMoleculeIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H_Hminus_to_H2_Electron

  double precision function H_Hminus_to_H2_Electron_Rate_Coefficient(temperature)
    !% Computes the rate coefficient (in units of c$^3$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{e}^-$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision, intent(in) :: temperature
    double precision             :: temperatureElectronVolts,logNaturalTemperatureElectronVolts
    
    ! Compute the temperature in electron volts.
    temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
    
    ! Compute the rate coefficient.
    if (temperatureElectronVolts >= 0.1d0) then
       ! Get the natural logarithm of the temperature in electron volts.
       logNaturalTemperatureElectronVolts=dlog(temperatureElectronVolts)

       H_Hminus_to_H2_Electron_Rate_Coefficient=dexp(               & 
            &                                      -20.069138970d0  &
            & +logNaturalTemperatureElectronVolts*(+ 0.228980000d0  &
            & +logNaturalTemperatureElectronVolts*(+ 3.599837700d-2 &
            & +logNaturalTemperatureElectronVolts*(- 4.555120000d-3 &
            & +logNaturalTemperatureElectronVolts*(- 3.105115440d-4 &
            & +logNaturalTemperatureElectronVolts*(+ 1.073294000d-4 &
            & +logNaturalTemperatureElectronVolts*(- 8.366719600d-6 &
            & +logNaturalTemperatureElectronVolts*(+ 2.238306230d-7 &
            &                                     )))))))           &
            &                                      )
    else
       H_Hminus_to_H2_Electron_Rate_Coefficient=1.428d-9
    end if
    return
  end function H_Hminus_to_H2_Electron_Rate_Coefficient

  subroutine Molecular_Hydrogen_Rate_H_Hplus_to_H2plus_Photon(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \gamma$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex&
         &,molecularHydrogenCationMoleculeIndex
    double precision                                  :: temperatureElectronVolts,rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0 .and. molecularHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt

       ! Compute the rate coefficient.
       if (temperatureElectronVolts < 0.577d0) then
          rateCoefficient=3.833d-16*(temperatureElectronVolts**1.8d0)
       else
          rateCoefficient=5.810d-16*((0.20651d0*temperatureElectronVolts)**(-0.2891d0*dlog(0.20651d0*temperatureElectronVolts)))
       end if

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H_Hplus_to_H2plus_Photon

  subroutine Molecular_Hydrogen_Rate_H2plus_H_to_H2_Hplus(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex&
         &,molecularHydrogenMoleculeIndex,molecularHydrogenCationMoleculeIndex
    double precision,                   parameter     :: rateCoefficient=6.4d-10
    double precision                                  :: rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2plus_H_to_H2_Hplus_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       molecularHydrogenMoleculeIndex      =Molecules_Index("MolecularHydrogen"      )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0 .and.&
            & molecularHydrogenMoleculeIndex > 0 .and. molecularHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2plus_H_to_H2_Hplus_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute the rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenCationMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex      ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2plus_H_to_H2_Hplus

  subroutine Molecular_Hydrogen_Rate_H2_Hplus_to_H2plus_H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{H}$.
    use Radiation_Structure
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: molecularHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex&
         &,molecularHydrogenCationMoleculeIndex,atomicHydrogenMoleculeIndex
    double precision                                  :: temperatureElectronVolts,logNaturalTemperatureElectronVolts &
         &,rateCoefficient,rate

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       molecularHydrogenMoleculeIndex      =Molecules_Index("MolecularHydrogen"      )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0 .and.&
            & molecularHydrogenMoleculeIndex > 0 .and. molecularHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       logNaturalTemperatureElectronVolts=dlog(temperatureElectronVolts)

       ! Compute rate coefficient.
       if (temperature < 1000.0d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=dexp(                                          &
               & -24.24914687d0                                          &
               & +3.400824440d0 * logNaturalTemperatureElectronVolts     &
               & -3.898003960d0 *(logNaturalTemperatureElectronVolts**2) &
               & +2.045558782d0 *(logNaturalTemperatureElectronVolts**3) &
               & -0.541618285d0 *(logNaturalTemperatureElectronVolts**4) &
               & +8.410775030d-2*(logNaturalTemperatureElectronVolts**5) &
               & -7.879026150d-3*(logNaturalTemperatureElectronVolts**6) &
               & +4.138398420d-4*(logNaturalTemperatureElectronVolts**7) &
               & -9.363458880d-6*(logNaturalTemperatureElectronVolts**8) &
               &                     )
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex      ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_Hplus_to_H2plus_H

  subroutine Molecular_Hydrogen_Rate_H2_Electron_to_2H_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{e}^- \rightarrow 2\hbox{H} + \hbox{e}^-$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,molecularHydrogenMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex   =Molecules_Index("AtomicHydrogen"   )
       molecularHydrogenMoleculeIndex=Molecules_Index("MolecularHydrogen")
       electronMoleculeIndex         =Molecules_Index("Electron"         )
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. molecularHydrogenMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       if (temperature < 1000.0d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=5.6d-11*dsqrt(temperature)*dexp(-102124.0d0/temperature)
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)*moleculeDensity%abundance(electronMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex   ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_Electron_to_2H_Electron

  subroutine Molecular_Hydrogen_Rate_H2_H_to_3H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H} \rightarrow 3\hbox{H}$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,molecularHydrogenMoleculeIndex
    double precision                                  :: temperatureElectronVolts,log10Temperature,rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex   =Molecules_Index("AtomicHydrogen"   )
       molecularHydrogenMoleculeIndex=Molecules_Index("MolecularHydrogen")
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. molecularHydrogenMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt

       ! Compute base 10 logarithm of temperature.
       log10Temperature=dlog10(temperature)

       ! Compute the rate coefficient.
       if (log10Temperature < 3.0d0 .or. log10Temperature > 5.4d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=1.067d-10*(temperatureElectronVolts**2.012d0)*dexp(-(4.463d0/temperatureElectronVolts)*((1.0d0+0.2472d0&
               &*temperatureElectronVolts)**3.512d0))
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex   ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_H_to_3H

  subroutine Molecular_Hydrogen_Rate_Hminus_Electron_to_H_2Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenAnionMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex     =Molecules_Index("AtomicHydrogen"     )
       atomicHydrogenAnionMoleculeIndex=Molecules_Index("AtomicHydrogenAnion")
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Get rate coefficient.
       rateCoefficient=Hminus_Electron_to_H_2Electron_Rate_Coefficient(temperature)

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)*moleculeDensity%abundance(electronMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex     , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex     ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex           , &
            & moleculeRates%abundance   (electronMoleculeIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hminus_Electron_to_H_2Electron

  double precision function Hminus_Electron_to_H_2Electron_Rate_Coefficient(temperature)
    !% Computes the rate coefficient (in units of cm$^3$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    double precision, intent(in) :: temperature
    double precision             :: temperatureElectronVolts,logNaturalTemperatureElectronVolts

    ! Compute the temperature in electron volts.
    temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
    logNaturalTemperatureElectronVolts=dlog(temperatureElectronVolts)
    
    ! Compute rate coefficient.
    Hminus_Electron_to_H_2Electron_Rate_Coefficient=dexp(       &
         &                                      -18.01849334d0  &
         & +logNaturalTemperatureElectronVolts*(+ 2.36085220d0  &
         & +logNaturalTemperatureElectronVolts*(- 0.28274430d0  &
         & +logNaturalTemperatureElectronVolts*(+ 1.62331664d-2 &
         & +logNaturalTemperatureElectronVolts*(- 3.36501203d-2 &
         & +logNaturalTemperatureElectronVolts*(+ 1.17832978d-2 &
         & +logNaturalTemperatureElectronVolts*(- 1.65619470d-3 &
         & +logNaturalTemperatureElectronVolts*(+ 1.06827520d-4 &
         & +logNaturalTemperatureElectronVolts*(- 2.63128581d-6 &
         &                                     ))))))))         &
         &                                              )
    return
  end function Hminus_Electron_to_H_2Electron_Rate_Coefficient

  subroutine Molecular_Hydrogen_Rate_Hminus_H_to_2H_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H} \rightarrow 2 \hbox{H} + \hbox{e}^-$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenAnionMoleculeIndex,electronMoleculeIndex
    double precision                                  :: temperatureElectronVolts,logNaturalTemperatureElectronVolts&
         &,rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex     =Molecules_Index("AtomicHydrogen"     )
       atomicHydrogenAnionMoleculeIndex=Molecules_Index("AtomicHydrogenAnion")
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Electron_to_Hminus_Photon_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       logNaturalTemperatureElectronVolts=dlog(temperatureElectronVolts)

       ! Compute rate coefficient.
       if (temperatureElectronVolts >= 0.1d0) then
          rateCoefficient=dexp(                                          &
               & -20.37260896d0                                          &
               & + 1.13944330d0 * logNaturalTemperatureElectronVolts     &
               & - 0.14210136d0 *(logNaturalTemperatureElectronVolts**2) &
               & + 8.46445540d-3*(logNaturalTemperatureElectronVolts**3) &
               & - 1.43276410d-3*(logNaturalTemperatureElectronVolts**4) &
               & + 2.01225030d-4*(logNaturalTemperatureElectronVolts**5) &
               & + 8.66396320d-5*(logNaturalTemperatureElectronVolts**6) &
               & - 2.58500970d-5*(logNaturalTemperatureElectronVolts**7) &
               & + 2.45550120d-6*(logNaturalTemperatureElectronVolts**8) &
               & - 8.06838250d-8*(logNaturalTemperatureElectronVolts**9) &
               &                     )
       else
          rateCoefficient=2.5634d-9*(temperatureElectronVolts**1.78186d0)
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex     , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex     ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex           , &
            & moleculeRates%abundance   (electronMoleculeIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hminus_H_to_2H_Electron

  subroutine Molecular_Hydrogen_Rate_Hminus_Hplus_to_2H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenAnionMoleculeIndex&
         &,atomicHydrogenCationMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_Hminus_Hplus_to_2H_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenCationMoleculeIndex=Molecules_Index("AtomicHydrogenCation")
       atomicHydrogenAnionMoleculeIndex =Molecules_Index("AtomicHydrogenAnion" )
       atomicHydrogenMoleculeIndex      =Molecules_Index("AtomicHydrogen"      )
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_Hminus_Hplus_to_2H_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Get the rate coefficient.
       rateCoefficient=Hminus_Hplus_to_2H_Rate_Coefficient(temperature)

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenCationMoleculeIndex)*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex , &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex      ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hminus_Hplus_to_2H

  double precision function Hminus_Hplus_to_2H_Rate_Coefficient(temperature)
    !% Compute the rate coefficient (in units of c$^3$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    implicit none
    double precision, intent(in) :: temperature
    
    Hminus_Hplus_to_2H_Rate_Coefficient=7.0d-8/dsqrt(temperature/100.0d0)
    return
  end function Hminus_Hplus_to_2H_Rate_Coefficient
  
  subroutine Molecular_Hydrogen_Rate_Hminus_Hplus_to_H2plus_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{e}^-$.
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenCationMoleculeIndex,atomicHydrogenAnionMoleculeIndex&
         &,molecularHydrogenCationMoleculeIndex,electronMoleculeIndex
    double precision                                  :: temperatureElectronVolts,rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_Hminus_Hplus_to_H2plus_Electron_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       atomicHydrogenAnionMoleculeIndex    =Molecules_Index("AtomicHydrogenAnion"    )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenCationMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0 .and.&
            & molecularHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_Hminus_Hplus_to_H2plus_Electron_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then
       
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       
       ! Compute rate coefficient.
       if (temperatureElectronVolts < 1.719d0) then
          rateCoefficient=2.2910d-10/(temperatureElectronVolts**0.4d0)
       else
          rateCoefficient=8.4258d-10/(temperatureElectronVolts**1.4d0)*dexp(-1.301d0/temperatureElectronVolts)
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)&
            &*moleculeDensity%abundance(atomicHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex    , &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex    ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex               , &
            & moleculeRates%abundance   (electronMoleculeIndex               ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hminus_Hplus_to_H2plus_Electron

  subroutine Molecular_Hydrogen_Rate_H2plus_Electron_to_2H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{e}^- \rightarrow 2\hbox{H}$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenMoleculeIndex,molecularHydrogenCationMoleculeIndex,electronMoleculeIndex
    double precision                                  :: rateCoefficient,rate

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2plus_Electron_to_2H_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       electronMoleculeIndex               =Molecules_Index("Electron"               )
       ! This reaction is active if both species were found.
       reactionActive=atomicHydrogenMoleculeIndex > 0 .and. molecularHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2plus_Electron_to_2H_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       if (temperature < 617.0d0) then
          rateCoefficient=1.0d-8
       else
          rateCoefficient=1.32d-6/(temperature**0.76d0)
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenCationMoleculeIndex)*moleculeDensity%abundance(electronMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & +rate*2.0d0               )
       call   moleculeRates%abundanceSet(electronMoleculeIndex               , &
            & moleculeRates%abundance   (electronMoleculeIndex               ) &
            & -rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2plus_Electron_to_2H

  subroutine Molecular_Hydrogen_Rate_H2plus_Hminus_to_H2_H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{H}$.
    use Radiation_Structure
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: molecularHydrogenCationMoleculeIndex,atomicHydrogenAnionMoleculeIndex&
         &,molecularHydrogenMoleculeIndex,atomicHydrogenMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (hydrogenNetworkFast) return

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2plus_Hminus_to_H2_H_Init)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       atomicHydrogenAnionMoleculeIndex    =Molecules_Index("AtomicHydrogenAnion"    )
       molecularHydrogenMoleculeIndex      =Molecules_Index("MolecularHydrogen"      )
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       ! This reaction is active if both species were found.
       reactionActive=molecularHydrogenCationMoleculeIndex > 0 .and. atomicHydrogenAnionMoleculeIndex > 0 .and.&
            & molecularHydrogenMoleculeIndex > 0 .and. atomicHydrogenMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2plus_Hminus_to_H2_H_Init)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=5.0d-7*dsqrt(100.0d0/temperature)

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenCationMoleculeIndex)&
            &*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex    , &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex    ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex      ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & +rate                     )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2plus_Hminus_to_H2_H

  subroutine Molecular_Hydrogen_Rate_Hminus_Gamma_to_H_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \gamma \rightarrow \hbox{H} + \hbox{e}^-$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    ! Energy range for the cross-section.
    double precision,                   parameter     :: crossSectionEnergyLow     =0.755d0
    ! Wavelength range for the cross-section.
    double precision,                   parameter     :: crossSectionWavelengthHigh=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyLow
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: atomicHydrogenAnionMoleculeIndex,electronMoleculeIndex&
         &,atomicHydrogenMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_Hminus_Gamma_to_H_Electron)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenAnionMoleculeIndex=Molecules_Index("AtomicHydrogenAnion")
       atomicHydrogenMoleculeIndex     =Molecules_Index("AtomicHydrogen"     )
       electronMoleculeIndex           =Molecules_Index("Electron"           )
       ! This reaction is active if all species were found.
       reactionActive=atomicHydrogenAnionMoleculeIndex > 0 .and. atomicHydrogenMoleculeIndex > 0 .and. electronMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_Hminus_Gamma_to_H_Electron)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       if (hydrogenNetworkCMBOnly) then
          ! <gfortran 4.6> Would like to use the radiation%temperature() type-bound procedure form here, but type-bound procedures
          ! with optional arguments are buggy under gfortran 4.4
          rateCoefficient=0.144d0*(Radiation_Temperature(radiation,[radiationTypeCMB])**2.13d0)*dexp(-8650.0d0&
               &/Radiation_Temperature(radiation,[radiationTypeCMB]))
       else
          rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_Hminus_Gamma_to_H_Electron,[0.0d0&
               &,crossSectionWavelengthHigh])
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenAnionMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenAnionMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenAnionMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex     , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex     ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex           , &
            & moleculeRates%abundance   (electronMoleculeIndex           ) &
            & +rate                     )

    end if
    return
  end subroutine Molecular_Hydrogen_Rate_Hminus_Gamma_to_H_Electron

  double precision function Cross_Section_Hminus_Gamma_to_H_Electron(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}^- + \gamma \rightarrow \hbox{H} + \hbox{e}^-$
    !% using the fitting function given by \cite{shapiro_hydrogen_1987}, renormalized\footnote{It seems unclear what units were
    !% used in \protect\cite{shapiro_hydrogen_1987}, hence the recalibration.} to match the results of
    !% \cite{nascimento_photodetachment_1977}.
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: wavelength
    double precision, parameter  :: energyThreshold=0.755d0
    double precision             :: energy

    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*angstromsPerMeter/electronVolt/wavelength

    ! Evaluate the fitting function for the cross-section.
    if (energy >=  energyThreshold) then
       Cross_Section_Hminus_Gamma_to_H_Electron=2.085d-16*(energy-energyThreshold)**1.5d0/energy**3
    else
       Cross_Section_Hminus_Gamma_to_H_Electron= 0.0d0
    end if
    return
  end function Cross_Section_Hminus_Gamma_to_H_Electron

  subroutine Molecular_Hydrogen_Rate_H2plus_Gamma_to_H_Hplus(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow \hbox{H} + \hbox{H}^+$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    ! Energy range for the cross-section.
    double precision,                   parameter     :: crossSectionEnergyLow     = 2.65d0
    double precision,                   parameter     :: crossSectionEnergyHigh    =21.00d0
    ! Wavelength range for the cross-section.
    double precision,                   parameter     :: crossSectionWavelengthLow =plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyHigh
    double precision,                   parameter     :: crossSectionWavelengthHigh=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyLow
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: molecularHydrogenCationMoleculeIndex,atomicHydrogenCationMoleculeIndex&
         &,atomicHydrogenMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2plus_Gamma_to_H_Hplus)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       atomicHydrogenMoleculeIndex         =Molecules_Index("AtomicHydrogen"         )
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       ! This reaction is active if all species were found.
       reactionActive=molecularHydrogenCationMoleculeIndex > 0 .and. atomicHydrogenMoleculeIndex > 0 .and. atomicHydrogenCationMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2plus_Gamma_to_H_Hplus)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       if (hydrogenNetworkCMBOnly) then
          ! <gfortran 4.6> Would like to use the radiation%temperature() type-bound procedure form here, but type-bound procedures
          ! with optional arguments are buggy under gfortran 4.4
          rateCoefficient=6.36d5*dexp(-71600.0d0/Radiation_Temperature(radiation,[radiationTypeCMB]))
       else
          rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_H2plus_Gamma_to_H_Hplus,[crossSectionWavelengthLow&
               &,crossSectionWavelengthHigh])
       end if

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex         , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex         ) &
            & +rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & +rate                     )

    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2plus_Gamma_to_H_Hplus

  double precision function Cross_Section_H2plus_Gamma_to_H_Hplus(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow \hbox{H} + \hbox{H}^+$
    !% as given by \cite{shapiro_hydrogen_1987}.
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: wavelength
    double precision             :: energy

    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*angstromsPerMeter/electronVolt/wavelength

    ! Evaluate the fitting function for the cross-section.
    if      (energy >=  2.65d0 .and. energy < 11.27d0) then
       Cross_Section_H2plus_Gamma_to_H_Hplus=10.0d0**(-40.97d0+energy*(+6.030d+0 &
            &                                                 +energy*(-0.504d+0 &
            &                                                 +energy*(+1.387d-2 &
            &                                                         )))        &
            &                                        )
    else if (energy >= 11.27d0 .and. energy < 21.00d0) then
       Cross_Section_H2plus_Gamma_to_H_Hplus=10.0d0**(-30.26d0+energy*(+2.790d+0 &
            &                                                 +energy*(-0.184d+0 &
            &                                                 +energy*(+3.535d-3 &
            &                                                         )))        &
            &                                        )
    else
       Cross_Section_H2plus_Gamma_to_H_Hplus= 0.0d0
    end if
    return
  end function Cross_Section_H2plus_Gamma_to_H_Hplus

  subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_H2star_to_2H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow H_2^* \rightarrow
    !% 2\hbox{H}$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    integer,                            save          :: molecularHydrogenMoleculeIndex,atomicHydrogenMoleculeIndex
    ! Median energy of the Lyman band in molecular hydrogen (in eV).
    double precision,                   parameter     :: energyLymanBand=12.87d0
    ! Corresponding median wavelength of the Lyman band in molecular hydrogen (in Angstroms).
    double precision,                   parameter     :: wavelengthLymanBand=angstromsPerMeter*plancksConstant*speedLight&
         &/(energyLymanBand*electronVolt)
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2_Gamma_to_H2star_to_2H)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       molecularHydrogenMoleculeIndex=Molecules_Index("MolecularHydrogen")
       atomicHydrogenMoleculeIndex   =Molecules_Index("AtomicHydrogen"   )
       ! This reaction is active if all species were found.
       reactionActive=molecularHydrogenMoleculeIndex > 0 .and. atomicHydrogenMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2_Gamma_to_H2star_to_2H)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=1.1d8*(4.0d0*Pi*Radiation_Flux(radiation,wavelengthLymanBand))

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex) &
            & -      rate               )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex   ) &
            & +2.0d0*rate               )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_H2star_to_2H

  subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_H2plus_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow \hbox{H}_2^+ + \hbox{e}^-$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    ! Energy of the edge in the cross-section.
    double precision,                   parameter     :: crossSectionEdgeEnergy    =15.42d0
    ! Wavelength of the edge in the cross-section.
    double precision,                   parameter     :: crossSectionEdgeWavelength=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEdgeEnergy
    integer,                            save          :: molecularHydrogenMoleculeIndex,molecularHydrogenCationMoleculeIndex &
         &,electronMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2_Gamma_to_H2plus_Electron)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       molecularHydrogenMoleculeIndex      =Molecules_Index("MolecularHydrogen"      )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       electronMoleculeIndex               =Molecules_Index("Electron"               )
       ! This reaction is active if all species were found.
       reactionActive=       molecularHydrogenMoleculeIndex       > 0 &
            &          .and. molecularHydrogenCationMoleculeIndex > 0 &
            &          .and. electronMoleculeIndex                > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2_Gamma_to_H2plus_Electron)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_H2_Gamma_to_H2plus_Electron,[0.0d0&
            &,crossSectionEdgeWavelength])

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex      ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex               , &
            & moleculeRates%abundance   (electronMoleculeIndex               ) &
            & +rate                     )

    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_H2plus_Electron

  double precision function Cross_Section_H2_Gamma_to_H2plus_Electron(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow \hbox{H}_2^+ +
    !% \hbox{e}^-$ as given by\footnote{\protect\cite{abel_modeling_1997} cite ``O'Neil \& Reinhardt (1978)'' as the source for
    !% this fit, but it is not listed in their bibliography, and I have not been able to locate by any other means.}
    !% \cite{abel_modeling_1997}.
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: wavelength
    double precision             :: energy

    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*angstromsPerMeter/electronVolt/wavelength

    ! Evaluate the fitting function for the cross-section.
    if      (energy < 15.42d0) then
       Cross_Section_H2_Gamma_to_H2plus_Electron=0.0d0
    else if (energy < 16.50d0) then
       Cross_Section_H2_Gamma_to_H2plus_Electron=6.2d-18*energy-9.40d-17
    else if (energy < 17.70d0) then
       Cross_Section_H2_Gamma_to_H2plus_Electron=1.4d-18*energy-1.48d-17
    else
       Cross_Section_H2_Gamma_to_H2plus_Electron=2.5d-14/(energy**2.71d0)
    end if
    return
  end function Cross_Section_H2_Gamma_to_H2plus_Electron

  subroutine Molecular_Hydrogen_Rate_H2plus_Gamma_to_2Hplus_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ +
    !% \hbox{e}^-$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    ! Energy range for the cross-section.
    double precision,                   parameter     :: crossSectionEnergyLow     =30.0d0
    double precision,                   parameter     :: crossSectionEnergyHigh    =90.0d0
    ! Wavelength range for the cross-section.
    double precision,                   parameter     :: crossSectionWavelengthLow =plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyHigh
    double precision,                   parameter     :: crossSectionWavelengthHigh=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyLow
    integer,                            save          :: atomicHydrogenCationMoleculeIndex,molecularHydrogenCationMoleculeIndex &
         &,electronMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2plus_Gamma_to_2Hplus_Electron)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenCationMoleculeIndex   =Molecules_Index("AtomicHydrogenCation"   )
       molecularHydrogenCationMoleculeIndex=Molecules_Index("MolecularHydrogenCation")
       electronMoleculeIndex               =Molecules_Index("Electron"               )
       ! This reaction is active if all species were found.
       reactionActive=       atomicHydrogenCationMoleculeIndex    > 0 &
            &          .and. molecularHydrogenCationMoleculeIndex > 0 &
            &          .and. electronMoleculeIndex                > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2plus_Gamma_to_2Hplus_Electron)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_H2plus_Gamma_to_2Hplus_Electron&
            &,[crossSectionWavelengthLow,crossSectionWavelengthHigh])

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenCationMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenCationMoleculeIndex) &
            & -      rate               )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex   ) &
            & +2.0d0*rate               )
       call   moleculeRates%abundanceSet(electronMoleculeIndex               , &
            & moleculeRates%abundance   (electronMoleculeIndex               ) &
            & +     rate                )
    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2plus_Gamma_to_2Hplus_Electron

  double precision function Cross_Section_H2plus_Gamma_to_2Hplus_Electron(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ +
    !% \hbox{e}^-$ as given by \cite{shapiro_hydrogen_1987}.
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: wavelength
    double precision             :: energy

    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*angstromsPerMeter/electronVolt/wavelength

    ! Evaluate the fitting function for the cross-section.
    if (energy >= 30.0d0 .and. energy <= 90.0d0) then
       Cross_Section_H2plus_Gamma_to_2Hplus_Electron=10.0d0**(         -16.926d+0 &
            &                                                 +energy*(- 4.528d-2 &
            &                                                 +energy*(  2.238d-4 &
            &                                                 +energy*(  4.245d-7 &
            &                                                         )))         &
            &                                                )
    else
       Cross_Section_H2plus_Gamma_to_2Hplus_Electron=0.0d0
    end if
    return
  end function Cross_Section_H2plus_Gamma_to_2Hplus_Electron

  subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_2H(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ +
    !% \hbox{e}^-$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    ! Energy range for the cross-section.
    double precision,                   parameter     :: crossSectionEnergyLow     =14.159d0
    double precision,                   parameter     :: crossSectionEnergyHigh    =17.700d0
    ! Wavelength range for the cross-section.
    double precision,                   parameter     :: crossSectionWavelengthLow =plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyHigh
    double precision,                   parameter     :: crossSectionWavelengthHigh=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyLow
    integer,                            save          :: atomicHydrogenMoleculeIndex,molecularHydrogenMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H2_Gamma_to_2H)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex   =Molecules_Index("AtomicHydrogen"   )
       molecularHydrogenMoleculeIndex=Molecules_Index("MolecularHydrogen")
       ! This reaction is active if all species were found.
       reactionActive=       atomicHydrogenMoleculeIndex    > 0 &
            &          .and. molecularHydrogenMoleculeIndex > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H2_Gamma_to_2H)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_H2_Gamma_to_2H&
            &,[crossSectionWavelengthLow,crossSectionWavelengthHigh])

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(molecularHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(molecularHydrogenMoleculeIndex, &
            & moleculeRates%abundance   (molecularHydrogenMoleculeIndex) &
            & -      rate               )
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex   , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex   ) &
            & +2.0d0*rate               )

    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H2_Gamma_to_2H

  double precision function Cross_Section_H2_Gamma_to_2H(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$ as given by
    !% \cite{abel_modeling_1997}.
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: wavelength
    double precision, parameter  :: ratioOrthoToPara=0.0d0 ! Assume all H_2 is in the para- configuration.
    double precision             :: energy,crossSectionLymanPara,crossSectionWernerPara,crossSectionLymanOrtho&
         &,crossSectionWernerOrtho

    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*angstromsPerMeter/electronVolt/wavelength
    
    ! Evaluate the Lyman and Wener band cross sections for para- and ortho- configurations.
    if      (energy > 14.675d0 .and. energy <= 16.820d0) then
       crossSectionLymanPara  =10.0d0**(-18.0d0+15.1289d0-1.0513900000d+0*energy                       )
    else if (                        energy <= 17.600d0) then
       crossSectionLymanPara  =10.0d0**(-18.0d0-31.4100d0+1.8042000000d-2*energy**3-4.2339d-5*energy**5)
    else
       crossSectionLymanPara  = 0.0d0
    end if
    if      (energy > 14.675d0 .and. energy <= 17.700d0) then
       crossSectionWernerPara =10.0d0**(-18.0d0+13.5311d0-0.9182618000d0*energy                        )
    else
       crossSectionWernerPara = 0.0d0
    end if
    if      (energy > 14.159d0 .and. energy <= 15.302d0) then
       crossSectionLymanOrtho =10.0d0**(-18.0d0+12.0218406d0-0.8194290d0*energy                        )
    else if (                        energy <= 17.200d0) then
       crossSectionLymanOrtho =10.0d0**(-18.0d0+16.0464400d0-1.0824380d0*energy                        )
    else
       crossSectionLymanOrtho = 0.0d0
    end if
    if      (energy > 14.159d0 .and. energy <= 17.200d0) then
       crossSectionWernerOrtho=10.0d0**(-18.0d0+12.8736700d0-0.85088597d0*energy                       )
    else
       crossSectionWernerOrtho= 0.0d0
    end if

    ! Construct the combined cross-section weighted by the appropriate ortho- to para- ratio.
    Cross_Section_H2_Gamma_to_2H= (      1.0d0/(ratioOrthoToPara+1.0d0))*(crossSectionLymanPara +crossSectionWernerPara ) &
         &                       +(1.0d0-1.0d0/(ratioOrthoToPara+1.0d0))*(crossSectionLymanOrtho+crossSectionWernerOrtho)

    return
  end function Cross_Section_H2_Gamma_to_2H

  subroutine Molecular_Hydrogen_Rate_H_Gamma_to_Hplus_Electron(temperature,radiation,moleculeDensity,moleculeRates)
    !% Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \gamma \rightarrow \hbox{H}^+ +
    !% \hbox{e}^-$.
    use Radiation_Structure
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    double precision,                   intent(in)    :: temperature
    type(radiationStructure),           intent(in)    :: radiation
    type(molecularAbundancesStructure), intent(in)    :: moleculeDensity
    type(molecularAbundancesStructure), intent(inout) :: moleculeRates
    logical,                            save          :: reactionInitialized=.false.,reactionActive=.false.
    ! Energy range for the cross-section (in eV).
    double precision,                   parameter     :: crossSectionEnergyLow     =13.60d0
    ! Wavelength range for the cross-section (in Angstroms).
    double precision,                   parameter     :: crossSectionWavelengthHigh=plancksConstant*speedLight*angstromsPerMeter&
         &/electronVolt/crossSectionEnergyLow
    integer,                            save          :: atomicHydrogenMoleculeIndex,atomicHydrogenCationMoleculeIndex&
         &,electronMoleculeIndex
    double precision                                  :: rate,rateCoefficient

    ! Check if this reaction needs initializing.
    !$omp critical(Molecular_Hydrogen_Rate_H_Gamma_to_H_Electron)
    if (.not.reactionInitialized) then
       ! Find the molecules in this reaction.
       atomicHydrogenMoleculeIndex      =Molecules_Index("AtomicHydrogen"      )
       atomicHydrogenCationMoleculeIndex=Molecules_Index("AtomicHydrogenCation")
       electronMoleculeIndex            =Molecules_Index("Electron"            )
       ! This reaction is active if all species were found.
       reactionActive=       atomicHydrogenMoleculeIndex       > 0 &
            &          .and. atomicHydrogenCationMoleculeIndex > 0 &
            &          .and. electronMoleculeIndex             > 0
       ! Flag that the reaction is now initialized.
       reactionInitialized=.true.
    end if
    !$omp end critical(Molecular_Hydrogen_Rate_H_Gamma_to_H_Electron)

    ! Do calculation if this reaction is active.
    if (reactionActive) then

       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection(Cross_Section_H_Gamma_to_Hplus_Electron&
            &,[0.0d0,crossSectionWavelengthHigh])

       ! Compute rate.
       rate=rateCoefficient*moleculeDensity%abundance(atomicHydrogenMoleculeIndex)

       ! Record rate.
       call   moleculeRates%abundanceSet(atomicHydrogenMoleculeIndex      , &
            & moleculeRates%abundance   (atomicHydrogenMoleculeIndex      ) &
            & -rate                     )
       call   moleculeRates%abundanceSet(atomicHydrogenCationMoleculeIndex, &
            & moleculeRates%abundance   (atomicHydrogenCationMoleculeIndex) &
            & +rate                     )
       call   moleculeRates%abundanceSet(electronMoleculeIndex            , &
            & moleculeRates%abundance   (electronMoleculeIndex            ) &
            & +rate                     )

    end if
    return
  end subroutine Molecular_Hydrogen_Rate_H_Gamma_to_Hplus_Electron

  double precision function Cross_Section_H_Gamma_to_Hplus_Electron(wavelength)
    !% Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$ as given by
    !% \cite{abel_modeling_1997}.
    use Atomic_Cross_Sections_Ionization_Photo
    implicit none
    double precision, intent(in) :: wavelength

    ! Use the hydrogen photoionization cross section method.
    Cross_Section_H_Gamma_to_Hplus_Electron=Atomic_Cross_Section_Ionization_Photo(1,1,1,wavelength)

    return
  end function Cross_Section_H_Gamma_to_Hplus_Electron

end module Molecular_Hydrogen_Rates
