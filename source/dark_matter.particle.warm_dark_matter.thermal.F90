!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a thermal warm dark matter particle class.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <darkMatterParticle name="darkMatterParticleWDMThermal">
   <description>Provides a thermal warm dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleClass) :: darkMatterParticleWDMThermal
     !!{
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
     !!{
     Constructors for the \normalfont \ttfamily WDMThermal} dark matter particle class.
     !!}
     module procedure wdmThermalConstructorParameters
     module procedure wdmThermalConstructorInternal
  end interface darkMatterParticleWDMThermal

contains

  function wdmThermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \normalfont \ttfamily WDMThermal} dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleWDMThermal)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    double precision                                              :: massValue           , degreesOfFreedomEffectiveValue

    !![
    <inputParameter>
      <name>mass</name>
      <source>parameters</source>
      <variable>massValue</variable>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass (in keV) of the thermal warm dark matter particle.</description>
    </inputParameter>
    <inputParameter>
      <name>degreesOfFreedomEffective</name>
      <source>parameters</source>
      <variable>degreesOfFreedomEffectiveValue</variable>
      <defaultValue>1.5d0</defaultValue>
      <description>The effective number of degrees of freedom for the thermal warm dark matter particle.</description>
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
    !!{
    Internal constructor for the \normalfont \ttfamily WDMThermal} dark matter particle class.
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
    !!{
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
    !!{
    Return the mass, in units of keV, of a thermal warm dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleWDMThermal), intent(inout) :: self

    wdmThermalMass=self%massValue
    return
  end function wdmThermalMass

  double precision function wdmThermalDegreesOfFreedomEffective(self)
    !!{
    Return the effective number of degrees of freedom of a thermal warm dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleWDMThermal), intent(inout) :: self

    wdmThermalDegreesOfFreedomEffective=self%degreesOfFreedomEffectiveValue
    return
  end function wdmThermalDegreesOfFreedomEffective

  double precision function wdmThermalDegreesOfFreedomEffectiveDecoupling(self)
    !!{
    Return the effective number of relativistic degrees of freedom at the time of decoupling of a thermal warm dark matter
    particle. The effective number of relativistic degrees of freedom at the time of decoupling is determined from the
    requirement that the warm dark matter particle have the correct relic density to provide the entire mass of dark matter
    \citep[][eqn.~17]{hogan_warm_1999}.
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

