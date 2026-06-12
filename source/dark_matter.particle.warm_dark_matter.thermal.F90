!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{RST
Implements a thermal warm dark matter particle class.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <darkMatterParticle name="darkMatterParticleWDMThermal" docformat="rst">
   <description>
   Provides a thermal relic warm dark matter particle whose free-streaming length suppresses small-scale structure formation. The particle mass in keV is set by ``[mass]``, and the effective number of relativistic degrees of freedom at freeze-out by ``[degreesOfFreedomEffective]``.
   </description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleClass) :: darkMatterParticleWDMThermal
     !!{RST
     A thermal warm dark matter particle class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: massValue                     , degreesOfFreedomEffectiveValue
   contains
     !![
     <methods>
       <method description="Return the effective number of degrees of freedom of the thermal wark dark matter particle." method="degreesOfFreedomEffective" />
       <method description="Return the effective number of relativisitc degrees of freedom in the universe at the time at which the thermal wark dark matter particle decoupled." method="degreesOfFreedomEffectiveDecoupling" />
     </methods>
     !!]
     final     ::                                        wdmThermalDestructor
     procedure :: mass                                => wdmThermalMass
     procedure :: degreesOfFreedomEffective           => wdmThermalDegreesOfFreedomEffective
     procedure :: degreesOfFreedomEffectiveDecoupling => wdmThermalDegreesOfFreedomEffectiveDecoupling
  end type darkMatterParticleWDMThermal

  interface darkMatterParticleWDMThermal
     !!{RST
     Constructors for the :galacticus-class:`darkMatterParticleWDMThermal` dark matter particle class.
     !!}
     module procedure wdmThermalConstructorParameters
     module procedure wdmThermalConstructorInternal
  end interface darkMatterParticleWDMThermal

contains

  function wdmThermalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterParticleWDMThermal` dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleWDMThermal)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    double precision                                              :: massValue           , degreesOfFreedomEffectiveValue

    !![
    <inputParameter docformat="rst">
      <name>mass</name>
      <source>parameters</source>
      <variable>massValue</variable>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The mass (in keV) of the thermal warm dark matter particle, which determines the free-streaming length and the suppression of small-scale structure relative to cold dark matter; smaller masses give stronger suppression.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>degreesOfFreedomEffective</name>
      <source>parameters</source>
      <variable>degreesOfFreedomEffectiveValue</variable>
      <defaultValue>1.5d0</defaultValue>
      <description>
      The effective number of degrees of freedom for the thermal warm dark matter particle.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=darkMatterParticleWDMThermal(massValue,degreesOfFreedomEffectiveValue,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function wdmThermalConstructorParameters

  function wdmThermalConstructorInternal(mass,degreesOfFreedomEffective,cosmologyParameters_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterParticleWDMThermal` dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleWDMThermal)                        :: self
    double precision                              , intent(in   )         :: mass                , degreesOfFreedomEffective
    class           (cosmologyParametersClass    ), intent(inout), target :: cosmologyParameters_
    !![
    <constructorAssign variables="*cosmologyParameters_"/>
    !!]

    self%massValue                     =mass
    self%degreesOfFreedomEffectiveValue=degreesOfFreedomEffective
    return
  end function wdmThermalConstructorInternal

  subroutine wdmThermalDestructor(self)
    !!{RST
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(darkMatterParticleWDMThermal), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine wdmThermalDestructor

  double precision function wdmThermalMass(self)
    !!{RST
    Return the mass, in units of keV, of a thermal warm dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleWDMThermal), intent(inout) :: self

    wdmThermalMass=self%massValue
    return
  end function wdmThermalMass

  double precision function wdmThermalDegreesOfFreedomEffective(self)
    !!{RST
    Return the effective number of degrees of freedom of a thermal warm dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleWDMThermal), intent(inout) :: self

    wdmThermalDegreesOfFreedomEffective=self%degreesOfFreedomEffectiveValue
    return
  end function wdmThermalDegreesOfFreedomEffective

  double precision function wdmThermalDegreesOfFreedomEffectiveDecoupling(self)
    !!{RST
    Return the effective number of relativistic degrees of freedom at the time of decoupling of a thermal warm dark matter particle. The effective number of relativistic degrees of freedom at the time of decoupling is determined from the requirement that the warm dark matter particle have the correct relic density to provide the entire mass of dark matter :cite:p:`hogan_warm_1999`.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    implicit none
    class(darkMatterParticleWDMThermal), intent(inout) :: self

    wdmThermalDegreesOfFreedomEffectiveDecoupling=+78.3d0                                                                      &
         &                                        *  self                     %degreesOfFreedomEffective(                  )   &
         &                                        *  self                     %mass                     (                  )   &
         &                                        /(                                                                           &
         &                                          +self%cosmologyParameters_%OmegaMatter              (                  )   &
         &                                          -self%cosmologyParameters_%OmegaBaryon              (                  )   &
         &                                        )                                                                            &
         &                                        /  self%cosmologyParameters_%HubbleConstant           (hubbleUnitsLittleH)**2
    return
  end function wdmThermalDegreesOfFreedomEffectiveDecoupling

