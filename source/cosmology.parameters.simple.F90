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

!% A simple implementation of the cosmological parameters class.

  !# <cosmologyParameters name="cosmologyParametersSimple" />

  type, extends(cosmologyParametersClass) :: cosmologyParametersSimple
     !% A simple cosmological parameters class.
     private
     double precision :: HubbleConstantValue, OmegaBaryonValue   , OmegaDarkEnergyValue, &
          &              OmegaMatterValue   , temperatureCMBValue
   contains
     final            :: cosmologyParametersSimpleDestructor
     procedure :: OmegaMatter    =>OmegaMatterSimple
     procedure :: OmegaDarkEnergy=>OmegaDarkEnergySimple
     procedure :: OmegaBaryon    =>OmegaBaryonSimple
     procedure :: OmegaRadiation =>OmegaRadiationSimple
     procedure :: OmegaCurvature =>OmegaCurvatureSimple
     procedure :: HubbleConstant =>HubbleConstantSimple
     procedure :: temperatureCMB =>temperatureCMBSimple
     procedure :: densityCritical=>densityCriticalSimple
  end type cosmologyParametersSimple

  interface cosmologyParametersSimple
     !% Constructors for the simple cosmological parameters class.
     module procedure cosmologyParametersSimpleDefaultConstructor
     module procedure cosmologyParametersSimpleConstructor
  end interface cosmologyParametersSimple

contains

  function cosmologyParametersSimpleDefaultConstructor()
    !% Default constructor for the simple cosmological parameters class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(cosmologyParametersSimple) :: cosmologyParametersSimpleDefaultConstructor

    !@ <inputParameter>
    !@   <name>Omega_Matter</name>
    !@   <defaultValue>0.2812 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The density of matter in the Universe in units of the critical density.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('Omega_Matter',cosmologyParametersSimpleDefaultConstructor%OmegaMatterValue,defaultValue=0.2812d0  )
    !@ <inputParameter>
    !@   <name>Omega_b</name>
    !@   <defaultValue>0.04611 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The density of baryons in the Universe in units of the critical density.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('Omega_b',cosmologyParametersSimpleDefaultConstructor%OmegaBaryonValue     ,defaultValue=0.04611d0 )
    !@ <inputParameter>
    !@   <name>Omega_DE</name>
    !@   <defaultValue>0.7188 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The density of dark energy in the Universe in units of the critical density.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('Omega_DE',cosmologyParametersSimpleDefaultConstructor%OmegaDarkEnergyValue,defaultValue=0.7188d0  )
    !@ <inputParameter>
    !@   <name>T_CMB</name>
    !@   <defaultValue>2.72548 \citep{fixsen_temperature_2009}</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The present day temperature of the \gls{cmb} in units of Kelvin.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('T_CMB',cosmologyParametersSimpleDefaultConstructor%temperatureCMBValue     ,defaultValue=2.72548d0)
    !@ <inputParameter>
    !@   <name>H_0</name>
    !@   <defaultValue>69.7 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The present day value of the Hubble parameter in units of km/s/Mpc.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('H_0',cosmologyParametersSimpleDefaultConstructor%HubbleConstantValue      ,defaultValue=69.7d0     )
    ! Validate the input.
    if (cosmologyParametersSimpleDefaultConstructor%HubbleConstantValue <= 0.0d0) &
         & call Galacticus_Display_Message("WARNING: H_0<=0 - are you sure this is what you wanted?",verbosityWarn)
    return
  end function cosmologyParametersSimpleDefaultConstructor

  function cosmologyParametersSimpleConstructor(OmegaMatter,OmegaBaryon,OmegaDarkEnergy,temperatureCMB,HubbleConstant)
    !% User-defined constructor for the simple cosmological parameters class.
    implicit none
    type            (cosmologyParametersSimple)                :: cosmologyParametersSimpleConstructor
    double precision                           , intent(in   ) :: HubbleConstant                      , OmegaBaryon, &
         &                                                        OmegaDarkEnergy                     , OmegaMatter, &
         &                                                        temperatureCMB

    cosmologyParametersSimpleConstructor%OmegaMatterValue    =OmegaMatter
    cosmologyParametersSimpleConstructor%OmegaBaryonValue    =OmegaBaryon
    cosmologyParametersSimpleConstructor%OmegaDarkEnergyValue=OmegaDarkEnergy
    cosmologyParametersSimpleConstructor%temperatureCMBValue =temperatureCMB
    cosmologyParametersSimpleConstructor%HubbleConstantValue =HubbleConstant
    return
  end function cosmologyParametersSimpleConstructor

  elemental subroutine cosmologyParametersSimpleDestructor(self)
    !% Default constructor for the simple cosmological parameters class.
    implicit none
    type(cosmologyParametersSimple), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine cosmologyParametersSimpleDestructor

  double precision function OmegaMatterSimple(self)
    !% Return the cosmological matter density in units of the critical density at the present day.
    use Galacticus_Error
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    OmegaMatterSimple=self%OmegaMatterValue
    return
  end function OmegaMatterSimple

  double precision function OmegaBaryonSimple(self)
    !% Return the cosmological baryon density in units of the critical density at the present day.
    use Galacticus_Error
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    OmegaBaryonSimple=self%OmegaBaryonValue
    return
  end function OmegaBaryonSimple

  double precision function OmegaDarkEnergySimple(self)
    !% Return the cosmological dark energy density in units of the critical density at the present day.
    use Galacticus_Error
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    OmegaDarkEnergySimple=self%OmegaDarkEnergyValue
    return
  end function OmegaDarkEnergySimple

  double precision function OmegaRadiationSimple(self)
    !% Return the cosmological radiation density in units of the critical density at the present day.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    OmegaRadiationSimple=radiationConstant*(self%temperatureCMB()**4)*megaParsec**3/massSolar/speedLight**2/self%densityCritical()
    return
  end function OmegaRadiationSimple

  double precision function OmegaCurvatureSimple(self)
    !% Return the cosmological curvature density in units of the critical density at the present day.
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    OmegaCurvatureSimple=1.0d0-self%OmegaMatter()-self%OmegaDarkEnergy()
    return
  end function OmegaCurvatureSimple

  double precision function temperatureCMBSimple(self)
    !% Return the present day temperature of the \gls{cmb}.
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    temperatureCMBSimple=self%temperatureCMBValue
    return
  end function TemperatureCMBSimple

  double precision function HubbleConstantSimple(self,units)
    !% Return the present day value of the Hubble constant.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    class           (cosmologyParametersSimple), intent(inout)           :: self
    integer                                    , intent(in   ), optional :: units
    integer                                                              :: unitsActual
    double precision                           , parameter               :: HubbleConstantNormalization=100.0d0

    unitsActual=unitsStandard
    if (present(units)) unitsActual=units
    select case (unitsActual)
    case (unitsStandard)
       HubbleConstantSimple=self%HubbleConstantValue
    case (unitsTime    )
       HubbleConstantSimple=self%HubbleConstantValue*gigaYear*kilo/megaParsec
    case (unitsLittleH )
       HubbleConstantSimple=self%HubbleConstantValue/HubbleConstantNormalization
    end select
    return
  end function HubbleConstantSimple

  double precision function densityCriticalSimple(self)
    !% Return the present day critical density of the Universe in units of $M_\odot$/Mpc$^3$.
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    class(cosmologyParametersSimple), intent(inout) :: self

    densityCriticalSimple=3.0d0*(self%HubbleConstant()**2)/8.0d0/Pi/gravitationalConstantGalacticus
    return
  end function densityCriticalSimple
