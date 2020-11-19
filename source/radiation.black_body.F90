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

  !% Implements a class for blackbody radiation fields.

  !# <radiationField name="radiationFieldBlackBody">
  !#  <description>A radiation field class for blackbody fields.</description>
  !# </radiationField>
  type, extends(radiationFieldClass) :: radiationFieldBlackBody
     !% A radiation field class for blackbody fields.
     private
     double precision :: temperature_
   contains
     !# <methods>
     !#   <method description="Return the temperature of te black-body radiation field." method="temperature" />
     !# </methods>
     procedure :: flux        => blackBodyFlux
     procedure :: temperature => blackBodyTemperature
  end type radiationFieldBlackBody

  interface radiationFieldBlackBody
     !% Constructors for the {\normalfont \ttfamily blackBody} radiation field class.
     module procedure blackBodyConstructorParameters
     module procedure blackBodyConstructorInternal
  end interface radiationFieldBlackBody

contains

  function blackBodyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily blackBody} radiation field class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (radiationFieldBlackBody)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    double precision                                         :: temperature_

    !# <inputParameter>
    !#   <name>temperature</name>
    !#   <variable>temperature_</variable>
    !#   <description>The temperature of the black body radiation field.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=radiationFieldBlackBody(temperature_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function blackBodyConstructorParameters

  function blackBodyConstructorInternal(temperature_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily blackBody} radiation field class.
    implicit none
    type            (radiationFieldBlackBody)                :: self
    double precision                         , intent(in   ) :: temperature_
    !# <constructorAssign variables="temperature_"/>

    return
  end function blackBodyConstructorInternal

  double precision function blackBodyFlux(self,wavelength,node)
    !% Return the flux of a blackBody radiation field.
    use :: Numerical_Constants_Prefixes, only : centi
    use :: Numerical_Constants_Units   , only : ergs
    use :: Thermodynamics_Radiation    , only : Blackbody_Emission, radianceTypeFrequency
    implicit none
    class           (radiationFieldBlackBody), intent(inout) :: self
    double precision                         , intent(in   ) :: wavelength
    type            (treeNode               ), intent(inout) :: node
    !$GLC attributes unused :: node

    blackBodyFlux=+centi**2                                                               &
         &        /ergs                                                                   &
         &        *Blackbody_Emission(wavelength,self%temperature_,radianceTypeFrequency)
    return
  end function blackBodyFlux

  double precision function blackBodyTemperature(self)
    !% Return the temperature of a black body radiation field.
    implicit none
    class(radiationFieldBlackBody), intent(inout) :: self

    blackBodyTemperature=self%temperature_
    return
  end function blackBodyTemperature
