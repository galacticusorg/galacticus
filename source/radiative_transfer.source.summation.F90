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

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  type, public :: radiativeTransferSourceList
     class(radiativeTransferSourceClass), pointer :: radiativeTransferSource => null()
     type (radiativeTransferSourceList ), pointer :: next                    => null()
  end type radiativeTransferSourceList

  !![
  <radiativeTransferSource name="radiativeTransferSourceSummation">
   <description>A photon source class for summation sources.</description>
   <linkedList type="radiativeTransferSourceList" variable="radiativeTransferSources" next="next" object="radiativeTransferSource" objectType="radiativeTransferSourceClass"/>
  </radiativeTransferSource>
  !!]
  type, extends(radiativeTransferSourceClass) :: radiativeTransferSourceSummation
     !!{
     Implementation of a summation source class for radiative transfer calculations.
     !!}
     private
     type   (radiativeTransferSourceList), pointer                   :: radiativeTransferSources => null()
     class  (randomNumberGeneratorClass ), pointer                   :: randomNumberGenerator_   => null()
     integer                                                         :: sourceTypeCount_
     type   (varying_string             ), allocatable, dimension(:) :: sourceTypeName_
   contains
     final     ::                           summationDestructor
     procedure :: initializePhotonPacket => summationInitializePhotonPacket
     procedure :: spectrum               => summationSpectrum
     procedure :: luminosity             => summationLuminosity
     procedure :: sourceTypeCount        => summationSourceTypeCount
     procedure :: sourceTypeName         => summationSourceTypeName
  end type radiativeTransferSourceSummation

  interface radiativeTransferSourceSummation
     !!{
     Constructors for the \refClass{radiativeTransferSourceSummation} radiative transfer source class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface radiativeTransferSourceSummation
  
contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferSourceSummation} radiative transfer source class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (radiativeTransferSourceSummation), target        :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    type   (radiativeTransferSourceList     ), pointer       :: radiativeTransferSource_
    integer                                                  :: i

    radiativeTransferSource_ => null()
    self%sourceTypeCount_    =  parameters%copiesCount('radiativeTransferSource',zeroIfNotPresent=.true.)
    allocate(self%sourceTypeName_(self%sourceTypeCount_))
    do i=1,parameters%copiesCount('radiativeTransferSource',zeroIfNotPresent=.true.)
       if (associated(radiativeTransferSource_)) then
          allocate(radiativeTransferSource_%next)
          radiativeTransferSource_ => radiativeTransferSource_%next
       else
          allocate(self%radiativeTransferSources)
          radiativeTransferSource_ => self%radiativeTransferSources
       end if
       !![
       <objectBuilder class="radiativeTransferSource" name="radiativeTransferSource_%radiativeTransferSource" source="parameters" copy="i" />
       !!]
       if (radiativeTransferSource_%radiativeTransferSource%sourceTypeCount() /= 1) &
            & call Error_Report('more than 1 source type is not presently supported'//{introspection:location})
       self%sourceTypeName_(i)=radiativeTransferSource_%radiativeTransferSource%sourceTypeName(1)
    end do
    !![
    <objectBuilder class="randomNumberGenerator" name="self%randomNumberGenerator_" source="parameters"/>
    <inputParametersValidate source="parameters" multiParameters="radiativeTransferSource"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(radiativeTransferSources,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferSourceSummation} radiative transfer source class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (radiativeTransferSourceSummation)                         :: self
    type (radiativeTransferSourceList     ), target , intent(in   ) :: radiativeTransferSources
    type (radiativeTransferSourceList     ), pointer                :: radiativeTransferSource_
    class(randomNumberGeneratorClass      ), target , intent(in   ) :: randomNumberGenerator_
    !![
    <constructorAssign variables="*randomNumberGenerator_"/>
    !!]

    self                    %radiativeTransferSources => radiativeTransferSources
    radiativeTransferSource_                          => radiativeTransferSources
    self%sourceTypeCount_                             =  0
    do while (associated(radiativeTransferSource_))
       !![
       <referenceCountIncrement owner="radiativeTransferSource_" object="radiativeTransferSource"/>
       !!]
       self%sourceTypeCount_    =  self%sourceTypeCount_        +1
       radiativeTransferSource_ => radiativeTransferSource_%next
    end do
    allocate(self%sourceTypeName_(self%sourceTypeCount_))
    radiativeTransferSource_ => radiativeTransferSources
    self%sourceTypeCount_    =  0
    do while (associated(radiativeTransferSource_))
       self%sourceTypeCount_    =  self%sourceTypeCount_        +1
       if (radiativeTransferSource_%radiativeTransferSource%sourceTypeCount() /= 1) &
            & call Error_Report('more than 1 source type is not presently supported'//{introspection:location})
       self%sourceTypeName_(self%sourceTypeCount_)=radiativeTransferSource_%radiativeTransferSource%sourceTypeName(1)
       radiativeTransferSource_ => radiativeTransferSource_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferSourceSummation} radiative transfer source class.
    !!}
    implicit none
    type(radiativeTransferSourceSummation), intent(inout) :: self
    type(radiativeTransferSourceList     ), pointer       :: radiativeTransferSource, radiativeTransferSourceNext

    radiativeTransferSource => self%radiativeTransferSources
    do while (associated(radiativeTransferSource))
       radiativeTransferSourceNext => radiativeTransferSource%next
       !![
       <objectDestructor name="radiativeTransferSource%radiativeTransferSource"/>
       !!]
       deallocate(radiativeTransferSource)
       radiativeTransferSource => radiativeTransferSourceNext
    end do
    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine summationDestructor
  
  subroutine summationInitializePhotonPacket(self,photonPacket)
    !!{
    Initialize the photon packet.
    !!}
    implicit none
    class           (radiativeTransferSourceSummation  ), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    type            (radiativeTransferSourceList       ), pointer       :: radiativeTransferSource
    double precision                                                    :: luminosityTotal        , luminosity
    integer                                                             :: sourceType

    ! Find the total luminosity in this wavelength range across all sources.
    luminosityTotal=self%luminosity(photonPacket%wavelengthMinimum(),photonPacket%wavelengthMaximum())
    ! Sample a luminosity at random.
    luminosity=luminosityTotal*self%randomNumberGenerator_%uniformSample()
    ! Find which source contributes this luminosity
    sourceType              =  0
    radiativeTransferSource => self%radiativeTransferSources
    do while (associated(radiativeTransferSource))
       luminosity=+                                                luminosity                                                                    &
            &     -radiativeTransferSource%radiativeTransferSource%luminosity(photonPacket%wavelengthMinimum(),photonPacket%wavelengthMaximum())
       sourceType=+sourceType                                                                                                                    &
            &     +1
       if (luminosity <= 0.0d0) exit
       radiativeTransferSource => radiativeTransferSource%next
    end do
    ! Set the properties from that source.
    call radiativeTransferSource%radiativeTransferSource%initializePhotonPacket(photonPacket)
    call photonPacket                                   %sourceTypeSet         (sourceType  )
    ! The luminosity is the sum of all sources - when randomly sampled this ensures that we get the correct total luminosity.
    call photonPacket%luminositySet(luminosityTotal)
    return
  end subroutine summationInitializePhotonPacket

  double precision function summationLuminosity(self,wavelengthMinimum,wavelengthMaximum,sourceType)
    !!{
    Return the spectrum of the summation source.
    !!}
    implicit none
    class           (radiativeTransferSourceSummation), intent(inout)           :: self
    double precision                                  , intent(in   )           :: wavelengthMinimum      , wavelengthMaximum
    integer                                           , intent(in   ), optional :: sourceType
    type            (radiativeTransferSourceList     ), pointer                 :: radiativeTransferSource
    integer                                                                     :: sourceIndex
    logical                                                                     :: accumulate

    summationLuminosity     =  0.0d0
    sourceIndex             =  0
    radiativeTransferSource => self%radiativeTransferSources
    do while (associated(radiativeTransferSource))
       sourceIndex                =  +sourceIndex                                                                                     &
            &                        +1
       accumulate                 =.not.present(sourceType)
       if (.not.accumulate) accumulate=sourceIndex == sourceType .or. sourceType == 0
       if (accumulate)                                                                                                                &
            & summationLuminosity =  +summationLuminosity                                                                             &
            &                        +radiativeTransferSource%radiativeTransferSource%luminosity(wavelengthMinimum,wavelengthMaximum)
       radiativeTransferSource    =>  radiativeTransferSource%next
    end do
    return
  end function summationLuminosity

  double precision function summationSpectrum(self,wavelength,sourceType)
    !!{
    Return the spectrum of the summation source.
    !!}
    implicit none
    class           (radiativeTransferSourceSummation), intent(inout)           :: self
    double precision                                  , intent(in   )           :: wavelength
    integer                                           , intent(in   ), optional :: sourceType
    type            (radiativeTransferSourceList     ), pointer                 :: radiativeTransferSource
    integer                                                                     :: sourceIndex
    logical                                                                     :: accumulate

    summationSpectrum       =  0.0d0
    sourceIndex             =  0
    radiativeTransferSource => self%radiativeTransferSources
    do while (associated(radiativeTransferSource))
       sourceIndex              =  +sourceIndex                                                          &
            &                      +1
       accumulate               =.not.present(sourceType)
       if (.not.accumulate) accumulate=sourceIndex == sourceType .or. sourceType == 0
       if (accumulate)                                                                                   &
            & summationSpectrum =  +summationSpectrum                                                    &
            &                      +radiativeTransferSource%radiativeTransferSource%spectrum(wavelength)
       radiativeTransferSource  =>  radiativeTransferSource%next
    end do
    return
  end function summationSpectrum

  integer function summationSourceTypeCount(self)
    !!{
    Return the number of source types provided.
    !!}
    implicit none
    class(radiativeTransferSourceSummation), intent(inout) :: self
    
    summationSourceTypeCount=self%sourceTypeCount_
    return
  end function summationSourceTypeCount

  function summationSourceTypeName(self,sourceType)
    !!{
    Return the name of the source type.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (varying_string                  )                :: summationSourceTypeName
    class  (radiativeTransferSourceSummation), intent(inout) :: self
    integer                                  , intent(in   ) :: sourceType

    if (sourceType < 1 .or. sourceType > self%sourceTypeCount_) call Error_Report('sourceType is out of range'//{introspection:location})
    summationSourceTypeName=self%sourceTypeName_(sourceType)
    return
  end function summationSourceTypeName
