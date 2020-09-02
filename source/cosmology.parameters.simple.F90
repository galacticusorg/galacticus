!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% A simple implementation of the cosmological parameters class.

  !# <cosmologyParameters name="cosmologyParametersSimple">
  !#  <description>Provides basic cosmological parameters: $(H_0,\Omega_\mathrm{M},\Omega_\Lambda,\Omega_\mathrm{b},T_\mathrm{CMB})$. Also provides derived quantities $(\Omega_\mathrm{K},\Omega_\mathrm{r},\rho_\mathrm{crit})$.</description>
  !# </cosmologyParameters>
  type, extends(cosmologyParametersClass) :: cosmologyParametersSimple
     !% A simple cosmological parameters class.
     private
     double precision :: HubbleConstantValue, OmegaBaryonValue   , OmegaDarkEnergyValue, &
          &              OmegaMatterValue   , temperatureCMBValue
   contains
     procedure :: OmegaMatter     => simpleOmegaMatter
     procedure :: OmegaDarkEnergy => simpleOmegaDarkEnergy
     procedure :: OmegaBaryon     => simpleOmegaBaryon
     procedure :: OmegaRadiation  => simpleOmegaRadiation
     procedure :: OmegaCurvature  => simpleOmegaCurvature
     procedure :: HubbleConstant  => simpleHubbleConstant
     procedure :: temperatureCMB  => simpleTemperatureCMB
     procedure :: densityCritical => simpleDensityCritical
  end type cosmologyParametersSimple

  interface cosmologyParametersSimple
     !% Constructors for the simple cosmological parameters class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface cosmologyParametersSimple

contains

  function simpleConstructorParameters(parameters)
    !% Constructor for the simple cosmological parameters class which takes a parameter set as input.
    use :: Galacticus_Error, only : Galacticus_Warn
    implicit none
    type(cosmologyParametersSimple)                :: simpleConstructorParameters
    type(inputParameters          ), intent(inout) :: parameters

    !# <inputParameter>
    !#   <name>OmegaMatter</name>
    !#   <source>parameters</source>
    !#   <variable>simpleConstructorParameters%OmegaMatterValue</variable>
    !#   <defaultValue>0.3153d0</defaultValue>
    !#   <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
    !#   <description>The density of matter in the Universe in units of the critical density.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>OmegaBaryon</name>
    !#   <source>parameters</source>
    !#   <variable>simpleConstructorParameters%OmegaBaryonValue</variable>
    !#   <defaultValue>0.04930d0</defaultValue>
    !#   <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
    !#   <description>The density of baryons in the Universe in units of the critical density.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>OmegaDarkEnergy</name>
    !#   <source>parameters</source>
    !#   <variable>simpleConstructorParameters%OmegaDarkEnergyValue</variable>
    !#   <defaultValue>0.6847d0</defaultValue>
    !#   <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
    !#   <description>The density of dark energy in the Universe in units of the critical density.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>temperatureCMB</name>
    !#   <source>parameters</source>
    !#   <variable>simpleConstructorParameters%temperatureCMBValue</variable>
    !#   <defaultValue>2.72548d0</defaultValue>
    !#   <defaultSource>\citep{fixsen_temperature_2009}</defaultSource>
    !#   <description>The present day temperature of the \gls{cmb} in units of Kelvin.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>HubbleConstant</name>
    !#   <source>parameters</source>
    !#   <variable>simpleConstructorParameters%HubbleConstantValue</variable>
    !#   <defaultValue>67.36d0</defaultValue>
    !#   <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
    !#   <description>The present day value of the Hubble parameter in units of km/s/Mpc.</description>
    !# </inputParameter>
    ! Validate the input.
    if (simpleConstructorParameters%HubbleConstantValue <= 0.0d0)                                                                                    &
         & call Galacticus_Warn("WARNING [cosmologyParametersSimple::simpleConstructorParameters]: H_0 ≤ 0 - are you sure this is what you wanted? "//{introspection:location})
    !# <inputParametersValidate source="parameters"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(OmegaMatter,OmegaBaryon,OmegaDarkEnergy,temperatureCMB,HubbleConstant)
    !% Internal constructor for the simple cosmological parameters class.
    implicit none
    type            (cosmologyParametersSimple)                :: simpleConstructorInternal
    double precision                           , intent(in   ) :: HubbleConstant   , OmegaBaryon, &
         &                                                        OmegaDarkEnergy  , OmegaMatter, &
         &                                                        temperatureCMB

    simpleConstructorInternal%OmegaMatterValue    =OmegaMatter
    simpleConstructorInternal%OmegaBaryonValue    =OmegaBaryon
    simpleConstructorInternal%OmegaDarkEnergyValue=OmegaDarkEnergy
    simpleConstructorInternal%temperatureCMBValue =temperatureCMB
    simpleConstructorInternal%HubbleConstantValue =HubbleConstant
    return
  end function simpleConstructorInternal

  double precision function simpleOmegaMatter(self)
    !% Return the cosmological matter density in units of the critical density at the present day.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleOmegaMatter=self%OmegaMatterValue
    return
  end function simpleOmegaMatter

  double precision function simpleOmegaBaryon(self)
    !% Return the cosmological baryon density in units of the critical density at the present day.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleOmegaBaryon=self%OmegaBaryonValue
    return
  end function simpleOmegaBaryon

  double precision function simpleOmegaDarkEnergy(self)
    !% Return the cosmological dark energy density in units of the critical density at the present day.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleOmegaDarkEnergy=self%OmegaDarkEnergyValue
    return
  end function simpleOmegaDarkEnergy

  double precision function simpleOmegaRadiation(self)
    !% Return the cosmological radiation density in units of the critical density at the present day.
    use :: Numerical_Constants_Astronomical, only : massSolar        , megaParsec
    use :: Numerical_Constants_Physical    , only : radiationConstant, speedLight
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleOmegaRadiation=+radiationConstant         &
         &               *megaParsec            **3 &
         &               /massSolar                 &
         &               /speedLight            **2 &
         &               *self%temperatureCMB ()**4 &
         &               /self%densityCritical()
    return
  end function simpleOmegaRadiation

  double precision function simpleOmegaCurvature(self)
    !% Return the cosmological curvature density in units of the critical density at the present day.
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleOmegaCurvature=+1.0d0                  &
         &               -self%OmegaMatter    () &
         &               -self%OmegaDarkEnergy()
    return
  end function simpleOmegaCurvature

  double precision function simpleTemperatureCMB(self)
    !% Return the present day temperature of the \gls{cmb}.
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleTemperatureCMB=self%temperatureCMBValue
    return
  end function simpleTemperatureCMB

  double precision function simpleHubbleConstant(self,units)
    !% Return the present day value of the Hubble constant.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : gigaYear               , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (cosmologyParametersSimple), intent(inout)           :: self
    integer                                    , intent(in   ), optional :: units
    double precision                           , parameter               :: HubbleConstantNormalization=100.0d0

    !# <optionalArgument name="units" defaultsTo="hubbleUnitsStandard" />
    select case (units_)
    case (hubbleUnitsStandard)
       simpleHubbleConstant=+self%HubbleConstantValue
    case (hubbleUnitsTime    )
       simpleHubbleConstant=+self%HubbleConstantValue         &
            &                    *gigaYear                    &
            &                    *kilo                        &
            &                    /megaParsec
    case (hubbleUnitsLittleH )
       simpleHubbleConstant=+self%HubbleConstantValue         &
            &                    /HubbleConstantNormalization
    case default
       simpleHubbleConstant=0.0d0
       call Galacticus_Error_Report('unknown units for Hubble parameter'//{introspection:location})
    end select
    return
  end function simpleHubbleConstant

  double precision function simpleDensityCritical(self)
    !% Return the present day critical density of the Universe in units of $M_\odot$/Mpc$^3$.
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    simpleDensityCritical=+3.0d0                                &
         &                *self%HubbleConstant            ()**2 &
         &                /8.0d0                                &
         &                /Pi                                   &
         &                /gravitationalConstantGalacticus
    return
  end function simpleDensityCritical
