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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a dust attenuation class in which no attenuation occurs.
  !!}

  !![
  <dustAttenuation name="dustAttenuationZero" docformat="rst">
   <description>
   A dust attenuation class in which nothing is attenuated: the transmission is unity at all wavelengths, for all
   components, and at all ages. This is the default, and is useful as a control against which the effect of a dust
   model can be measured.

   Because it depends on nothing, it requests no decomposition at all---not even by component---so that wrapping a
   luminosity in this attenuator costs a single parcel per output element.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationZero
     !!{RST
     A dust attenuation class in which no attenuation occurs.
     !!}
     private
   contains
     procedure :: transmission => zeroTransmission
     procedure :: request      => zeroRequest
  end type dustAttenuationZero

  interface dustAttenuationZero
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationZero` dust attenuation class.
     !!}
     module procedure zeroConstructorParameters
  end interface dustAttenuationZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationZero` dust attenuation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(dustAttenuationZero)                :: self
    type(inputParameters    ), intent(inout) :: parameters

    self=dustAttenuationZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  function zeroTransmission(self,node,descriptors) result(transmission)
    !!{RST
    Return unit transmission---no attenuation---for every parcel.
    !!}
    implicit none
    class           (dustAttenuationZero), intent(inout)                               :: self
    type            (treeNode           ), intent(inout), target                       :: node
    type            (emissionDescriptor ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                    , dimension(size(descriptors)) :: transmission
    !$GLC attributes unused :: self, node, descriptors

    transmission=1.0d0
    return
  end function zeroTransmission

  function zeroRequest(self) result(request)
    !!{RST
    Return a decomposition request specifying no resolution at all, since the transmission depends on nothing.
    !!}
    implicit none
    type (decompositionRequest)                :: request
    class(dustAttenuationZero ), intent(inout) :: self
    !$GLC attributes unused :: self

    request%resolveComponents =.false.
    request%resolveMetallicity=.false.
    request%resolveRadius     =.false.
    return
  end function zeroRequest
