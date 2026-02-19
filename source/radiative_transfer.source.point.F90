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

  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Radiative_Transfer_Spectra, only : radiativeTransferSpectrumClass
  
  !![
  <radiativeTransferSource name="radiativeTransferSourcePoint">
   <description>A photon source class for point sources.</description>
  </radiativeTransferSource>
  !!]
  type, extends(radiativeTransferSourceClass) :: radiativeTransferSourcePoint
     !!{
     Implementation of a point source class for radiative transfer calculations.
     !!}
     private
     class           (radiativeTransferSpectrumClass), pointer      :: radiativeTransferSpectrum_ => null()
     class           (randomNumberGeneratorClass    ), pointer      :: randomNumberGenerator_     => null()
     double precision                                , dimension(3) :: position
     type            (varying_string                )               :: label
   contains
     final     ::                           pointDestructor
     procedure :: initializePhotonPacket => pointInitializePhotonPacket
     procedure :: luminosity             => pointLuminosity
     procedure :: spectrum               => pointSpectrum
     procedure :: sourceTypeCount        => pointSourceTypeCount
     procedure :: sourceTypeName         => pointSourceTypeName
  end type radiativeTransferSourcePoint

  interface radiativeTransferSourcePoint
     !!{
     Constructors for the \refClass{radiativeTransferSourcePoint} radiative transfer source class.
     !!}
     module procedure pointConstructorParameters
     module procedure pointConstructorInternal
  end interface radiativeTransferSourcePoint
  
contains

  function pointConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferSourcePoint} radiative transfer source class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (radiativeTransferSourcePoint  )                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                , dimension(3)  :: position
    class           (radiativeTransferSpectrumClass), pointer       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass    ), pointer       :: randomNumberGenerator_
    type            (varying_string                )                :: label

    !![
    <inputParameter>
      <name>position</name>
      <defaultValue>[0.0d0,0.0d0,0.0d0]</defaultValue>
      <description>The position of the point source.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>A descriptive label for the source.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="radiativeTransferSpectrum" name="radiativeTransferSpectrum_" source="parameters"/>
    <objectBuilder class="randomNumberGenerator"     name="randomNumberGenerator_"     source="parameters"/>
    !!]
    self=radiativeTransferSourcePoint(position,label,radiativeTransferSpectrum_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="radiativeTransferSpectrum_"/>
    <objectDestructor name="randomNumberGenerator_"    />
    !!]
    return
  end function pointConstructorParameters

  function pointConstructorInternal(position,label,radiativeTransferSpectrum_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferSourcePoint} radiative transfer source class.
    !!}
    implicit none
    type            (radiativeTransferSourcePoint  )                              :: self
    double precision                                , intent(in   ), dimension(3) :: position
    type            (varying_string                ), intent(in   )               :: label
    class           (radiativeTransferSpectrumClass), intent(in   ), target       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass    ), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="position, label, *radiativeTransferSpectrum_, *randomNumberGenerator_"/>
    !!]
    
    return
  end function pointConstructorInternal

  subroutine pointDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferSourcePoint} radiative transfer source class.
    !!}
    implicit none
    type(radiativeTransferSourcePoint), intent(inout) :: self

    !![
    <objectDestructor name="self%radiativeTransferSpectrum_"/>
    <objectDestructor name="self%randomNumberGenerator_"    />
    !!]
    return
  end subroutine pointDestructor

  subroutine pointInitializePhotonPacket(self,photonPacket)
    !!{
    Initialize the wavelength of the photon packet.
    !!}
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
    call photonPacket%luminositySet(self%luminosity(photonPacket%wavelengthMinimum(),photonPacket%wavelengthMaximum()))
    ! Set a source type.
    call photonPacket%sourceTypeSet(1)
    return
  end subroutine pointInitializePhotonPacket

  double precision function pointSpectrum(self,wavelength,sourceType)
    !!{
    Return the spectrum of the point source.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiativeTransferSourcePoint), intent(inout)           :: self
    double precision                              , intent(in   )           :: wavelength
    integer                                       , intent(in   ), optional :: sourceType

    if (present(sourceType)) then
       if (sourceType /= 0 .and. sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    end if
    pointSpectrum=self%radiativeTransferSpectrum_%spectrum(wavelength)
    return
  end function pointSpectrum

  double precision function pointLuminosity(self,wavelengthMinimum,wavelengthMaximum,sourceType)
    !!{
    Return the luminosity of the point source.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiativeTransferSourcePoint), intent(inout)           :: self
    double precision                              , intent(in   )           :: wavelengthMinimum, wavelengthMaximum
    integer                                       , intent(in   ), optional :: sourceType

    if (present(sourceType)) then
       if (sourceType /= 0 .and. sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    end if
    pointLuminosity=self%radiativeTransferSpectrum_%luminosity(wavelengthMinimum,wavelengthMaximum)
    return
  end function pointLuminosity

  integer function pointSourceTypeCount(self)
    !!{
    Return the number of source types provided.
    !!}
    implicit none
    class(radiativeTransferSourcePoint), intent(inout) :: self
    !$GLC attributes unused :: self
    
    pointSourceTypeCount=1
    return
  end function pointSourceTypeCount

  function pointSourceTypeName(self,sourceType)
    !!{
    Return the name of the source type.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (varying_string              )                :: pointSourceTypeName
    class  (radiativeTransferSourcePoint), intent(inout) :: self
    integer                              , intent(in   ) :: sourceType

    if (sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    pointSourceTypeName=self%label
    return
  end function pointSourceTypeName
