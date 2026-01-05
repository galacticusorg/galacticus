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

  use :: Mass_Distributions        , only : massDistributionClass
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Radiative_Transfer_Spectra, only : radiativeTransferSpectrumClass
  
  !![
  <radiativeTransferSource name="radiativeTransferSourceDistributed">
   <description>A photon source class for distributed sources.</description>
  </radiativeTransferSource>
  !!]
  type, extends(radiativeTransferSourceClass) :: radiativeTransferSourceDistributed
     !!{
     Implementation of a distributed source class for radiative transfer calculations.
     !!}
     private
     class           (massDistributionClass         ), pointer      :: massDistribution_          => null()
     class           (radiativeTransferSpectrumClass), pointer      :: radiativeTransferSpectrum_ => null()
     class           (randomNumberGeneratorClass    ), pointer      :: randomNumberGenerator_     => null()
     double precision                                , dimension(3) :: position
     type            (varying_string                )               :: label
   contains
     final     ::                           distributedDestructor
     procedure :: initializePhotonPacket => distributedInitializePhotonPacket
     procedure :: luminosity             => distributedLuminosity
     procedure :: spectrum               => distributedSpectrum
     procedure :: sourceTypeCount        => distributedSourceTypeCount
     procedure :: sourceTypeName         => distributedSourceTypeName
  end type radiativeTransferSourceDistributed

  interface radiativeTransferSourceDistributed
     !!{
     Constructors for the \refClass{radiativeTransferSourceDistributed} radiative transfer source class.
     !!}
     module procedure distributedConstructorParameters
     module procedure distributedConstructorInternal
  end interface radiativeTransferSourceDistributed
  
contains

  function distributedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferSourceDistributed} radiative transfer source class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (radiativeTransferSourceDistributed)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                    , dimension(3)  :: position
    class           (massDistributionClass             ), pointer       :: massDistribution_
    class           (radiativeTransferSpectrumClass    ), pointer       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass        ), pointer       :: randomNumberGenerator_
    type            (varying_string                    )                :: label

    !![
    <inputParameter>
      <name>position</name>
      <defaultValue>[0.0d0,0.0d0,0.0d0]</defaultValue>
      <description>The position of the distributed source.</description>
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
    <objectBuilder class="massDistribution"          name="massDistribution_"          source="parameters"/>
    !!]
    self=radiativeTransferSourceDistributed(position,label,massDistribution_,radiativeTransferSpectrum_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="radiativeTransferSpectrum_"/>
    <objectDestructor name="randomNumberGenerator_"    />
    <objectDestructor name="massDistribution_"         />
    !!]
    return
  end function distributedConstructorParameters

  function distributedConstructorInternal(position,label,massDistribution_,radiativeTransferSpectrum_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferSourceDistributed} radiative transfer source class.
    !!}
    implicit none
    type            (radiativeTransferSourceDistributed)                              :: self
    double precision                                    , intent(in   ), dimension(3) :: position
    type            (varying_string                    ), intent(in   )               :: label
    class           (massDistributionClass             ), intent(in   ), target       :: massDistribution_
    class           (radiativeTransferSpectrumClass    ), intent(in   ), target       :: radiativeTransferSpectrum_
    class           (randomNumberGeneratorClass        ), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="position, label, *massDistribution_, *radiativeTransferSpectrum_, *randomNumberGenerator_"/>
    !!]
 
    return
  end function distributedConstructorInternal

  subroutine distributedDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferSourceDistributed} radiative transfer source class.
    !!}
    implicit none
    type(radiativeTransferSourceDistributed), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"         />
    <objectDestructor name="self%radiativeTransferSpectrum_"/>
    <objectDestructor name="self%randomNumberGenerator_"    />
    !!]
    return
  end subroutine distributedDestructor

  subroutine distributedInitializePhotonPacket(self,photonPacket)
    !!{
    Initialize the wavelength of the photon packet.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (radiativeTransferSourceDistributed), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    double precision                                                    :: theta       , phi

    ! Choose a position in the distributed source.
    call photonPacket%positionSet(                                                                    &
         &                        +self%massDistribution_%positionSample(self%randomNumberGenerator_) &
         &                        +self                  %position                                    &
         &                       )
    ! Choose a propagation direction isotropically at random.
    phi  =+     2.0d0*Pi*self%randomNumberGenerator_%uniformSample()
    theta=+acos(2.0d0   *self%randomNumberGenerator_%uniformSample()-1.0d0)
    call photonPacket%directionSet(                      &
         &                         [                     &
         &                          sin(theta)*cos(phi), &
         &                          sin(theta)*sin(phi), &
         &                          cos(theta)           &
         &                         ]                     &
         &                        )
    call photonPacket%luminositySet(self%luminosity(photonPacket%wavelengthMinimum(),photonPacket%wavelengthMaximum()))
    ! Set a source type.
    call photonPacket%sourceTypeSet(1)
    return
  end subroutine distributedInitializePhotonPacket

  double precision function distributedSpectrum(self,wavelength,sourceType)
    !!{
    Return the spectrum of the distributed source.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiativeTransferSourceDistributed), intent(inout)           :: self
    double precision                                    , intent(in   )           :: wavelength
    integer                                             , intent(in   ), optional :: sourceType

    if (present(sourceType)) then
       if (sourceType /= 0 .and. sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    end if
    distributedSpectrum=self%radiativeTransferSpectrum_%spectrum(wavelength)
    return
  end function distributedSpectrum

  double precision function distributedLuminosity(self,wavelengthMinimum,wavelengthMaximum,sourceType)
    !!{
    Return the luminosity of the distributed source.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiativeTransferSourceDistributed), intent(inout)           :: self
    double precision                                    , intent(in   )           :: wavelengthMinimum, wavelengthMaximum
    integer                                             , intent(in   ), optional :: sourceType

    if (present(sourceType)) then
       if (sourceType /= 0 .and. sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    end if
    distributedLuminosity=self%radiativeTransferSpectrum_%luminosity(wavelengthMinimum,wavelengthMaximum)
    return
  end function distributedLuminosity

  integer function distributedSourceTypeCount(self)
    !!{
    Return the number of source types provided.
    !!}
    implicit none
    class(radiativeTransferSourceDistributed), intent(inout) :: self
    !$GLC attributes unused :: self
    
    distributedSourceTypeCount=1
    return
  end function distributedSourceTypeCount

  function distributedSourceTypeName(self,sourceType)
    !!{
    Return the name of the source type.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (varying_string                    )                :: distributedSourceTypeName
    class  (radiativeTransferSourceDistributed), intent(inout) :: self
    integer                                    , intent(in   ) :: sourceType

    if (sourceType /= 1) call Error_Report('sourceType is out of range'//{introspection:location})
    distributedSourceTypeName=self%label
    return
  end function distributedSourceTypeName
