!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a dark matter halo mass function class which modifies another mass function using a simple model for systematics.

  !# <haloMassFunction name="haloMassFunctionSimpleSystematic">
  !#  <description>The halo mass function is computed by modifying another halo mass function using a simple model for systematic errors.</description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionSimpleSystematic
     !% A halo mass function class which modifies another mass function using a simple model for systematics.
     private
     double precision                                 :: alpha                , beta
     class           (haloMassFunctionClass), pointer :: referenceMassFunction
    contains
     final     ::                 simpleSystematicDestructor
     procedure :: differential => simpleSystematicDifferential
  end type haloMassFunctionSimpleSystematic

  interface haloMassFunctionSimpleSystematic
     !% Constructors for the {\normalfont \ttfamily simpleSystematic} halo mass function class.
     module procedure simpleSystematicConstructorParameters
     module procedure simpleSystematicConstructorInternal
  end interface haloMassFunctionSimpleSystematic

contains

  function simpleSystematicConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily simpleSystematic} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionSimpleSystematic)                :: simpleSystematicConstructorParameters
    type(inputParameters                 ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <source>parameters</source>
    !#   <variable>simpleSystematicConstructorParameters%alpha</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Parameter $\alpha$ appearing in model for simple systematic shift in the halo mass function.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <source>parameters</source>
    !#   <variable>simpleSystematicConstructorParameters%beta</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Parameter $\beta$ appearing in model for simple systematic shift in the halo mass function.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="simpleSystematicConstructorParameters%cosmologyParameters_"  source="parameters"/>
    !# <objectBuilder class="haloMassFunction"    name="simpleSystematicConstructorParameters%referenceMassFunction" source="parameters"/>
   return
  end function simpleSystematicConstructorParameters

  function simpleSystematicConstructorInternal(cosmologyParameters_,referenceMassFunction)
    !% Internal constructor for the {\normalfont \ttfamily simpleSystematic} halo mass function class.
    implicit none
    type (haloMassFunctionSimpleSystematic)                        :: simpleSystematicConstructorInternal
    class(cosmologyParametersClass        ), target, intent(in   ) :: cosmologyParameters_    
    class(haloMassFunctionClass           ), target, intent(in   ) :: referenceMassFunction

    simpleSystematicConstructorInternal%cosmologyParameters_  => cosmologyParameters_
    simpleSystematicConstructorInternal%referenceMassFunction => referenceMassFunction
    return
  end function simpleSystematicConstructorInternal
  
  subroutine simpleSystematicDestructor(self)
    !% Destructor for the {\normalfont \ttfamily simpleSystematic} halo mass function class.
    implicit none
    type(haloMassFunctionSimpleSystematic), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%referenceMassFunction" />
    return
  end subroutine simpleSystematicDestructor

  double precision function simpleSystematicDifferential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    implicit none
    class           (haloMassFunctionSimpleSystematic), intent(inout) :: self
    double precision                                  , intent(in   ) :: time                , mass
    double precision                                  , parameter     :: massZeroPoint=1.0d12
    
    ! Compute the systematic shift in the halo mass function.
    simpleSystematicDifferential=+self%referenceMassFunction%differential(time,mass) &
         &                       *(                                                  &
         &                         +1.0d0                                            &
         &                         +self%alpha                                       &
         &                         +self%beta                                        &
         &                         *log10(                                           &
         &                                 mass                                      &
         &                                /massZeroPoint                             &
         &                               )                                           &
         &                        )
    return
  end function simpleSystematicDifferential
