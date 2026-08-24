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
  Implements a dust attenuation class for dust mixed uniformly with the emitting stars.
  !!}

  !![
  <dustAttenuation name="dustAttenuationMixedSlab" docformat="rst">
   <description>
   Attenuation by dust mixed uniformly with the emitting stars, rather than lying in a screen in front of them. For a
   slab in which emitters and absorbers are uniformly interspersed the emergent fraction is

   .. math::

      T = \frac{1-\mathrm{e}^{-\tau}}{\tau},

   which tends to unity as :math:`\tau \rightarrow 0` and, crucially, falls only as :math:`1/\tau` at large optical
   depth rather than exponentially: stars on the near side of the slab always escape, so a mixed geometry can never
   extinguish a component completely. This is the qualitative difference from a screen, and it matters whenever the
   inferred optical depths are large.

   This class is a decorator: it takes the transmission :math:`T_\mathrm{screen}` computed by any other attenuator,
   recovers the equivalent optical depth :math:`\tau = -\ln T_\mathrm{screen}`, and re-applies it in the mixed
   geometry. Any prescription for how much dust a galaxy has can therefore be used in either geometry without being
   reimplemented.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationMixedSlab
     !!{RST
     Attenuation by dust mixed uniformly with the emitting stars.
     !!}
     private
     class(dustAttenuationClass), pointer :: dustAttenuation_ => null()
   contains
     final     ::                      mixedSlabDestructor
     procedure :: transmission      => mixedSlabTransmission
     procedure :: request           => mixedSlabRequest
     procedure :: supportsComponent => mixedSlabSupportsComponent
  end type dustAttenuationMixedSlab

  interface dustAttenuationMixedSlab
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationMixedSlab` dust attenuation class.
     !!}
     module procedure mixedSlabConstructorParameters
     module procedure mixedSlabConstructorInternal
  end interface dustAttenuationMixedSlab

contains

  function mixedSlabConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationMixedSlab` dust attenuation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (dustAttenuationMixedSlab)                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(dustAttenuationClass    ), pointer       :: dustAttenuation_

    !![
    <objectBuilder class="dustAttenuation" name="dustAttenuation_" source="parameters"/>
    !!]
    self=dustAttenuationMixedSlab(dustAttenuation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustAttenuation_"/>
    !!]
    return
  end function mixedSlabConstructorParameters

  function mixedSlabConstructorInternal(dustAttenuation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationMixedSlab` dust attenuation class.
    !!}
    implicit none
    type (dustAttenuationMixedSlab)                        :: self
    class(dustAttenuationClass    ), intent(in   ), target :: dustAttenuation_
    !![
    <constructorAssign variables="*dustAttenuation_"/>
    !!]

    return
  end function mixedSlabConstructorInternal

  subroutine mixedSlabDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationMixedSlab` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationMixedSlab), intent(inout) :: self

    !![
    <objectDestructor name="self%dustAttenuation_"/>
    !!]
    return
  end subroutine mixedSlabDestructor

  function mixedSlabTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through dust mixed uniformly with the emitting stars.
    !!}
    implicit none
    class           (dustAttenuationMixedSlab), intent(inout)                               :: self
    type            (treeNode                ), intent(inout), target                       :: node
    type            (emissionDescriptor      ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                    , intent(in   ), optional           :: inclination
    double precision                                         , dimension(size(descriptors)) :: transmission
    ! Below this optical depth the series expansion of [1-exp(-τ)]/τ is used: the direct expression suffers
    ! cancellation as τ →, where both numerator and denominator vanish.
    double precision                                         , parameter                    :: depthOpticalSmall=1.0d-6
    integer                                                                                 :: i
    double precision                                                                        :: depthOptical
    !$GLC attributes unused :: inclination

    transmission=self%dustAttenuation_%transmission(node,descriptors)
    do i=1,size(transmission)
       if (transmission(i) >= 1.0d0) then
          ! No attenuation at all.
          transmission(i)=1.0d0
       else if (transmission(i) <= 0.0d0) then
          ! The screen is completely opaque, so the equivalent optical depth is infinite and the mixed transmission
          ! vanishes -- albeit only as 1/τ.
          transmission(i)=0.0d0
       else
          depthOptical=-log(transmission(i))
          if (depthOptical < depthOpticalSmall) then
             transmission(i)=+1.0d0                 &
                  &          -depthOptical   /2.0d0 &
                  &          +depthOptical**2/6.0d0
          else
             transmission(i)=+(                    &
                  &            +1.0d0              &
                  &            -exp(-depthOptical) &
                  &           )                    &
                  &          /depthOptical
          end if
       end if
    end do
    return
  end function mixedSlabTransmission

  function mixedSlabRequest(self) result(request)
    !!{RST
    Return the decomposition request of the wrapped attenuator: changing the geometry does not change which
    properties the attenuation depends upon.
    !!}
    implicit none
    type (decompositionRequest    )                :: request
    class(dustAttenuationMixedSlab), intent(inout) :: self

    request=self%dustAttenuation_%request()
    return
  end function mixedSlabRequest

  logical function mixedSlabSupportsComponent(self,componentType) result(supportsComponent)
    !!{RST
    Return whether the wrapped attenuator supports the given component.
    !!}
    implicit none
    class(dustAttenuationMixedSlab    ), intent(inout) :: self
    type (enumerationComponentTypeType), intent(in   ) :: componentType

    supportsComponent=self%dustAttenuation_%supportsComponent(componentType)
    return
  end function mixedSlabSupportsComponent
