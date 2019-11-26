!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Radiative_Transfer_Spectra, only : radiativeTransferSpectrumClass
  
  !# <radiativeTransferSource name="radiativeTransferSourcePoint">
  !#  <description>A photon source class for point sources.</description>
  !# </radiativeTransferSource>
  type, extends(radiativeTransferSourceClass) :: radiativeTransferSourcePoint
     !% Implementation of a point photon packet class which tracks only wavelength, position, and direction.
     private
     class           (radiativeTransferSpectrumClass), pointer      :: radiativeTransferSpectrum_ => null()
     class           (randomNumberGeneratorClass    ), pointer      :: randomNumberGenerator_     => null()
     double precision                                , dimension(3) :: position
   contains
     final     ::                           pointDestructor
     procedure :: initializePhotonPacket => pointInitializePhotonPacket
     procedure :: spectrum               => pointSpectrum
  end type radiativeTransferSourcePoint

  interface radiativeTransferSourcePoint
     !% Constructors for the {\normalfont \ttfamily point} radiative transfer photon packet class.
     module procedure pointConstructorParameters
     module procedure pointConstructorInternal
  end interface radiativeTransferSourcePoint
  
contains

  function pointConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily point} radiative transfer photon packet class which takes a parameter set as
    !% input.
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferSourcePoint  )                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                , dimension(3)  :: position
    class           (radiativeTransferSpectrumClass), pointer       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass    ), pointer       :: randomNumberGenerator_

    !# <inputParameter>
    !#   <name>position</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>[0.0d0,0.0d0,0.0d0]</defaultValue>
    !#   <description>The position of the point source.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="radiativeTransferSpectrum" name="radiativeTransferSpectrum_" source="parameters"/>
    !# <objectBuilder class="randomNumberGenerator"     name="randomNumberGenerator_"     source="parameters"/>
    self=radiativeTransferSourcePoint(position,radiativeTransferSpectrum_,randomNumberGenerator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="radiativeTransferSpectrum_"/>
    !# <objectDestructor name="randomNumberGenerator_"    />
    return
  end function pointConstructorParameters

  function pointConstructorInternal(position,radiativeTransferSpectrum_,randomNumberGenerator_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily point} radiative transfer source class.
    implicit none
    type            (radiativeTransferSourcePoint  )                              :: self
    double precision                                , intent(in   ), dimension(3) :: position
    class           (radiativeTransferSpectrumClass), intent(in   ), target       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass    ), intent(in   ), target       :: randomNumberGenerator_
    !# <constructorAssign variables="position, *radiativeTransferSpectrum_, *randomNumberGenerator_"/>
    
    return
  end function pointConstructorInternal

  subroutine pointDestructor(self)
    !% Destructor for the {\normalfont \ttfamily point} radiative transfer source class.
    implicit none
    type(radiativeTransferSourcePoint), intent(inout) :: self

    !# <objectDestructor name="self%radiativeTransferSpectrum_"/>
    !# <objectDestructor name="self%randomNumberGenerator_"    />
    return
  end subroutine pointDestructor

  subroutine pointInitializePhotonPacket(self,photonPacket)
    !% Set the wavelength of the photon packet.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (radiativeTransferSourcePoint      ), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    double precision                                                    :: theta       , phi

    ! Choose a propagation direction isotropically at random.
    phi  =+     2.0d0*Pi*self%randomNumberGenerator_%uniformSample()
    theta=+acos(2.0d0   *self%randomNumberGenerator_%uniformSample()-1.0d0)
    ! Set photon packet properties.
    call photonPacket%positionSet  (self%position)
    call photonPacket%directionSet (                      &
         &                          [                     &
         &                           sin(theta)*cos(phi), &
         &                           sin(theta)*sin(phi), &
         &                           cos(theta)           &
         &                          ]                     &
         &                         )
    call photonPacket%luminositySet(self%radiativeTransferSpectrum_%luminosity(photonPacket%wavelengthMinimum(),photonPacket%wavelengthMaximum()))
    return
  end subroutine pointInitializePhotonPacket

  double precision function pointSpectrum(self,wavelength)
    !% Return the spectrum of the point source.
    implicit none
    class           (radiativeTransferSourcePoint), intent(inout) :: self
    double precision                              , intent(in   ) :: wavelength

    pointSpectrum=self%radiativeTransferSpectrum_%spectrum(wavelength)
    return
  end function pointSpectrum
