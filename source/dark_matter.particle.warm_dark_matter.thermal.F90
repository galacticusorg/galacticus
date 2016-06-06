!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a thermal warm dark matter particle class.

  use Cosmology_Parameters
  
  !# <darkMatterParticle name="darkMatterParticleWDMThermal">
  !#  <description>Provides a thermal warm dark matter particle.</description>
  !# </darkMatterParticle>
  type, extends(darkMatterParticleClass) :: darkMatterParticleWDMThermal
     !% A thermal warm dark matter particle class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     double precision                                    :: massValue           , degreesOfFreedomEffectiveValue
   contains
     procedure :: mass                                => wdmThermalMass
     procedure :: degreesOfFreedomEffective           => wdmThermalDegreesOfFreedomEffective
     procedure :: degreesOfFreedomEffectiveDecoupling => wdmThermalDegreesOfFreedomEffectiveDecoupling
  end type darkMatterParticleWDMThermal
  
  interface darkMatterParticleWDMThermal
     !% Constructors for the ``{\normalfont \ttfamily WDMThermal}'' dark matter particle class.
     module procedure wdmThermalConstructorParameters
     module procedure wdmThermalConstructorInternal
  end interface darkMatterParticleWDMThermal
  
contains

  function wdmThermalConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily WDMThermal}'' dark matter particle class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (darkMatterParticleWDMThermal)                :: wdmThermalConstructorParameters
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_    
    double precision                                              :: mass                           , degreesOfFreedomEffective
    !# <inputParameterList label="allowedParameterNames" />

    !# <inputParameter>
    !#   <name>mass</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The mass (in keV) of the theral warm dark matter particle.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>degreesOfFreedomEffective</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.5d0</defaultValue>
    !#   <description>The effective number of degrees of freedom for the thermal warm dark matter particle.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    wdmThermalConstructorParameters=darkMatterParticleWDMThermal(mass,degreesOfFreedomEffective,cosmologyParameters_)
    return
  end function wdmThermalConstructorParameters

  function wdmThermalConstructorInternal(mass,degreesOfFreedomEffective,cosmologyParameters_)
    !% Internal constructor for the ``{\normalfont \ttfamily WDMThermal}'' dark matter particle class.
    use Input_Parameters2
    implicit none
    type            (darkMatterParticleWDMThermal)                        :: wdmThermalConstructorInternal
    double precision                              , intent(in   )         :: mass                         , degreesOfFreedomEffective
    class           (cosmologyParametersClass    ), intent(inout), target :: cosmologyParameters_    
    
    wdmThermalConstructorInternal%massValue                      =  mass
    wdmThermalConstructorInternal%degreesOfFreedomEffectiveValue =  degreesOfFreedomEffective
    wdmThermalConstructorInternal%cosmologyParameters_           => cosmologyParameters_
    return
  end function wdmThermalConstructorInternal
  
  subroutine wdmThermalDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (darkMatterParticleWDMThermal), intent(inout) :: self
    type     (inputParameters             ), intent(inout) :: descriptor
    type     (inputParameters             )                :: subParameters
    character(len=10                      )                :: parameterLabel

    call descriptor%addParameter("darkMatterParticleMethod","WDMThermal")
    subParameters=descriptor%subparameters("criticalOverdensityMethod")
    write (parameterLabel,'(f10.6)') self%massValue
    call subParameters%addParameter("mass"                     ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%degreesOfFreedomEffectiveValue
    call subParameters%addParameter("degreesOfFreedomEffective",trim(adjustl(parameterLabel)))
    return
  end subroutine wdmThermalDescriptor

  double precision function wdmThermalMass(self)
    !% Return the mass, in units of keV, of a thermal warm dark matter particle.
    implicit none
    class(darkMatterParticleWDMThermal), intent(in   ) :: self

    wdmThermalMass=self%massValue
    return
  end function wdmThermalMass
  
  double precision function wdmThermalDegreesOfFreedomEffective(self)
    !% Return the effective number of degrees of freedom of a thermal warm dark matter particle.
    implicit none
    class(darkMatterParticleWDMThermal), intent(in   ) :: self

    wdmThermalDegreesOfFreedomEffective=self%degreesOfFreedomEffectiveValue
    return
  end function wdmThermalDegreesOfFreedomEffective
  
  double precision function wdmThermalDegreesOfFreedomEffectiveDecoupling(self)
    !% Return the effective number of relativistic degrees of freedom at the time of decoupling of a thermal warm dark matter
    !% particle. The effective number of relativistic degrees of freedom at the time of decoupling is determined from the
    !% requirement that the warm dark matter particle have the correct relic density to provide the entire mass of dark matter
    !% \citep[][eqn.~17]{hogan_warm_1999}.
    implicit none
    class(darkMatterParticleWDMThermal), intent(in   ) :: self

    wdmThermalDegreesOfFreedomEffectiveDecoupling=+78.3d0                                                                      &
         &                                        *  self                     %degreesOfFreedomEffective(                  )   &
         &                                        *  self                     %mass                     (                  )   &
         &                                        /(                                                                           &
         &                                          +self%cosmologyParameters_%OmegaMatter             (                  )    &
         &                                          -self%cosmologyParameters_%OmegaBaryon             (                  )    &
         &                                        )                                                                            &
         &                                        /  self%cosmologyParameters_%HubbleConstant          (hubbleUnitsLittleH)**2
    return
  end function wdmThermalDegreesOfFreedomEffectiveDecoupling
  
