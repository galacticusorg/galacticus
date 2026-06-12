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

  !![
  <radiativeTransferPhotonPacket name="radiativeTransferPhotonPacketSimple" docformat="rst">
   <description>
   A simple photon packet class for Monte Carlo radiative transfer that tracks the minimum required properties: wavelength, luminosity, three-dimensional position, and propagation direction. The initial wavelength and luminosity of each packet are set by ``[wavelength]`` and ``[luminosity]``, with wavelength bounds given by ``[wavelengthMinimum]`` and ``[wavelengthMaximum]``.
   </description>
  </radiativeTransferPhotonPacket>
  !!]
  type, extends(radiativeTransferPhotonPacketClass) :: radiativeTransferPhotonPacketSimple
     !!{RST
     Implementation of a simple photon packet class which tracks only wavelength, position, and direction.
     !!}
     private
     double precision               :: wavelengthMinimum_, wavelengthMaximum_, &
          &                            wavelength_       , luminosity_
     double precision, dimension(3) :: position_         , direction_
     integer                        :: sourceType_
   contains
     procedure :: wavelength           => simpleWavelength
     procedure :: wavelengthSet        => simpleWavelengthSet
     procedure :: wavelengthMinimum    => simpleWavelengthMinimum
     procedure :: wavelengthMinimumSet => simpleWavelengthMinimumSet
     procedure :: wavelengthMaximum    => simpleWavelengthMaximum
     procedure :: wavelengthMaximumSet => simpleWavelengthMaximumSet
     procedure :: luminosity           => simpleLuminosity
     procedure :: luminositySet        => simpleLuminositySet
     procedure :: position             => simplePosition
     procedure :: positionSet          => simplePositionSet
     procedure :: direction            => simpleDirection
     procedure :: directionSet         => simpleDirectionSet
     procedure :: sourceType           => simpleSourceType
     procedure :: sourceTypeSet        => simpleSourceTypeSet
  end type radiativeTransferPhotonPacketSimple

  interface radiativeTransferPhotonPacketSimple
     !!{RST
     Constructors for the :galacticus-class:`radiativeTransferPhotonPacketSimple` radiative transfer photon packet class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface radiativeTransferPhotonPacketSimple
  
contains

  function simpleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`radiativeTransferPhotonPacketSimple` radiative transfer photon packet class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferPhotonPacketSimple)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: wavelengthMinimum, wavelengthMaximum, &
         &                                                                  wavelength       , luminosity

    !![
    <inputParameter docformat="rst">
      <name>wavelength</name>
      <defaultValue>1.0d4</defaultValue>
      <description>
      The wavelength of the photon packet (in \AA).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthMinimum</name>
      <defaultValue>0.5d4</defaultValue>
      <description>
      The minimum wavelength of the photon packet (in \AA).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthMaximum</name>
      <defaultValue>2.0d4</defaultValue>
      <description>
      The maximum wavelength of the photon packet (in \AA).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>luminosity</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The luminosity of the photon packet (in :math:`L_\odot`).
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=radiativeTransferPhotonPacketSimple(wavelength,wavelengthMinimum,wavelengthMaximum,luminosity)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(wavelength,wavelengthMinimum,wavelengthMaximum,luminosity) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`radiativeTransferPhotonPacketSimple` radiative transfer photon packet class.
    !!}
    implicit none
    type            (radiativeTransferPhotonPacketSimple)                :: self
    double precision                                     , intent(in   ) :: wavelengthMinimum, wavelengthMaximum, &
         &                                                                  wavelength       , luminosity

    self%wavelength_       =wavelength
    self%wavelengthMinimum_=wavelengthMinimum
    self%wavelengthMaximum_=wavelengthMaximum
    self%luminosity_       =luminosity
    self%sourceType_       =0
    return
  end function simpleConstructorInternal

  subroutine simpleWavelengthSet(self,wavelength)
    !!{RST
    Set the wavelength of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self
    double precision                                     , intent(in   ) :: wavelength

    self%wavelength_=wavelength
    return
  end subroutine simpleWavelengthSet

  double precision function simpleWavelength(self)
    !!{RST
    Return the wavelength of the photon packet.
    !!}
    implicit none
    class(radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleWavelength=self%wavelength_
    return
  end function simpleWavelength

  subroutine simpleWavelengthMinimumSet(self,wavelength)
    !!{RST
    Set the minimum wavelength of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self
    double precision                                     , intent(in   ) :: wavelength

    self%wavelengthMinimum_=wavelength
    return
  end subroutine simpleWavelengthMinimumSet

  double precision function simpleWavelengthMinimum(self)
    !!{RST
    Return the minimum wavelength of the photon packet.
    !!}
    implicit none
    class(radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleWavelengthMinimum=self%wavelengthMinimum_
    return
  end function simpleWavelengthMinimum

  subroutine simpleWavelengthMaximumSet(self,wavelength)
    !!{RST
    Set the maximum wavelength of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self
    double precision                                     , intent(in   ) :: wavelength

    self%wavelengthMaximum_=wavelength
    return
  end subroutine simpleWavelengthMaximumSet

  double precision function simpleWavelengthMaximum(self)
    !!{RST
    Return the maximum wavelength of the photon packet.
    !!}
    implicit none
    class(radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleWavelengthMaximum=self%wavelengthMaximum_
    return
  end function simpleWavelengthMaximum

  subroutine simpleLuminositySet(self,luminosity)
    !!{RST
    Set the luminosity of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self
    double precision                                     , intent(in   ) :: luminosity

    self%luminosity_=luminosity
    return
  end subroutine simpleLuminositySet

  double precision function simpleLuminosity(self)
    !!{RST
    Return the luminosity of the photon packet.
    !!}
    implicit none
    class(radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleLuminosity=self%luminosity_
    return
  end function simpleLuminosity
  
  subroutine simplePositionSet(self,position)
    !!{RST
    Set the position of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout)               :: self
    double precision                                     , intent(in   ), dimension(3) :: position

    self%position_=position
    return
  end subroutine simplePositionSet

  function simplePosition(self)
    !!{RST
    Return the position of the photon packet.
    !!}
    implicit none
    double precision                                     , dimension(3)  :: simplePosition
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simplePosition=self%position_
    return
  end function simplePosition

  subroutine simpleDirectionSet(self,direction)
    !!{RST
    Set the direction of the photon packet.
    !!}
    implicit none
    class           (radiativeTransferPhotonPacketSimple), intent(inout)               :: self
    double precision                                     , intent(in   ), dimension(3) :: direction

    self%direction_=direction
    return
  end subroutine simpleDirectionSet

  function simpleDirection(self)
    !!{RST
    Return the direction of the photon packet.
    !!}
    implicit none
    double precision                                     , dimension(3)  :: simpleDirection
    class           (radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleDirection=self%direction_
    return
  end function simpleDirection

  integer function simpleSourceType(self)
    !!{RST
    Return the source type for this photon packet.
    !!}
    implicit none
    class(radiativeTransferPhotonPacketSimple), intent(inout) :: self

    simpleSourceType=self%sourceType_
    return
  end function simpleSourceType

  subroutine simpleSourceTypeSet(self,sourceType)
    !!{RST
    Set the source type for this photon packet.
    !!}
    implicit none
    class  (radiativeTransferPhotonPacketSimple), intent(inout) :: self
    integer                                     , intent(in   ) :: sourceType

    self%sourceType_=sourceType
    return
  end subroutine simpleSourceTypeSet
