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
  <radiativeTransferOutputter name="radiativeTransferOutputterLymanContinuumRate">
   <description>A radiative transfer outputter class which outputs the Lyman continuum photon emission rate.</description>
  </radiativeTransferOutputter>
  !!]
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterLymanContinuumRate
     !!{
     Implementation of a radiative transfer outputter class which outputs the Lyman continuum photon emission rate.
     !!}
     private
     double precision                                            :: lymanContinuumRateEscaping
     double precision                , allocatable, dimension(:) :: lymanContinuumRateEscapingTagged
     type            (varying_string), allocatable, dimension(:) :: sourceTypeName
   contains
     procedure :: reset               => lymanContinuumRateReset
     procedure :: sourceProperties    => lymanContinuumRateSourceProperties
     procedure :: photonPacketEscapes => lymanContinuumRatePhotonPacketEscapes
     procedure :: finalize            => lymanContinuumRateFinalize
     procedure :: output              => lymanContinuumRateOutput
  end type radiativeTransferOutputterLymanContinuumRate

  interface radiativeTransferOutputterLymanContinuumRate
     !!{
     Constructors for the \refClass{radiativeTransferOutputterLymanContinuumRate} radiative transfer outputter packet class.
     !!}
     module procedure lymanContinuumRateConstructorParameters
     module procedure lymanContinuumRateConstructorInternal
  end interface radiativeTransferOutputterLymanContinuumRate
  
contains

  function lymanContinuumRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferOutputterLymanContinuumRate} radiative transfer outputter class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(radiativeTransferOutputterLymanContinuumRate)                :: self
    type(inputParameters                             ), intent(inout) :: parameters
    
    self=radiativeTransferOutputterLymanContinuumRate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function lymanContinuumRateConstructorParameters
  
  function lymanContinuumRateConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferOutputterLymanContinuumRate} radiative transfer outputter class.
    !!}
    implicit none
    type(radiativeTransferOutputterLymanContinuumRate) :: self

    self%lymanContinuumRateEscaping=0.0d0
    return
  end function lymanContinuumRateConstructorInternal
  
  subroutine lymanContinuumRateReset(self)
    !!{
    Reset the accumulated Lyman continuum photon escape rate.
    !!}
    implicit none
    class(radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self

    self%lymanContinuumRateEscaping=0.0d0
    if (allocated(self%lymanContinuumRateEscapingTagged)) self%lymanContinuumRateEscapingTagged=0.0d0
    return
  end subroutine lymanContinuumRateReset

  subroutine lymanContinuumRateSourceProperties(self,radiativeTransferSource_,outputGroup)
    !!{
    Compute and output the Lyman continuum photon emission rate.
    !!}
    use :: HDF5_Access               , only : hdf5Access
    use :: ISO_Varying_String        , only : var_str                                  , operator(//)
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    use :: Numerical_Integration     , only : integrator
    use :: MPI_Utilities             , only : mpiSelf
    use :: String_Handling           , only : String_Upper_Case_First
    implicit none
    class           (radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    class           (radiativeTransferSourceClass                ), intent(inout) :: radiativeTransferSource_
    type            (hdf5Object                                  ), intent(inout) :: outputGroup
    type            (integrator                                  )                :: integrator_
    double precision                                                              :: rateLymanContinuum
    integer                                                                       :: sourceIndex
    type            (varying_string                              )                :: label

    if (.not.allocated(self%lymanContinuumRateEscapingTagged)) then
       allocate(self%lymanContinuumRateEscapingTagged(radiativeTransferSource_%sourceTypeCount()))
       allocate(self%sourceTypeName                  (radiativeTransferSource_%sourceTypeCount()))
       self%lymanContinuumRateEscapingTagged=0.0d0
       do sourceIndex=1,size(self%lymanContinuumRateEscapingTagged)
          self%sourceTypeName(sourceIndex)=radiativeTransferSource_%sourceTypeName(sourceIndex)
       end do
    end if
    if (mpiSelf%isMaster()) then
       integrator_       = integrator           (                                                             &
            &                                                      integrand                                , &
            &                                    toleranceRelative=1.0d-2                                     &
            &                                   )
       do sourceIndex=0,size(self%lymanContinuumRateEscapingTagged)
          rateLymanContinuum=+integrator_%integrate(                                                             &
               &                                                      1.0d-6*lymanSeriesLimitWavelengthHydrogen_atomic, &
               &                                                             lymanSeriesLimitWavelengthHydrogen_atomic  &
               &                                   )
          label=var_str('rateLymanContinuumEmitted')
          if (sourceIndex > 0) label=label//String_Upper_Case_First(char(radiativeTransferSource_%sourceTypeName(sourceIndex)))
          !$ call hdf5Access%set  ()
          call outputGroup%writeAttribute(rateLymanContinuum,char(label))
          !$ call hdf5Access%unset()
       end do
    end if
    return

  contains

    double precision function integrand(wavelength)
      !!{
      Integrand over the source spectrum.
      !!}
      use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
      use :: Numerical_Constants_Units       , only : metersToAngstroms
      use :: Numerical_Constants_Astronomical, only : luminositySolar
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: energyPhoton
      
      energyPhoton=+plancksConstant                                           &
           &       *speedLight                                                &
           &       *metersToAngstroms                                         &
           &       /wavelength
      integrand   =+radiativeTransferSource_%spectrum(wavelength,sourceIndex) &
           &       *luminositySolar                                           &
           &       /energyPhoton
      return
    end function integrand

  end subroutine lymanContinuumRateSourceProperties

  subroutine lymanContinuumRatePhotonPacketEscapes(self,photonPacket)
    !!{
    Process an escaping photon packet.
    !!}
    use :: Numerical_Constants_Atomic      , only : lymanSeriesLimitWavelengthHydrogen_atomic
    use :: Numerical_Constants_Physical    , only : plancksConstant                   , speedLight
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    class           (radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass          ), intent(inout) :: photonPacket
    double precision                                                              :: energyPhoton, rateEscape
    integer                                                                       :: sourceIndex
    
    if (photonPacket%wavelength() < lymanSeriesLimitWavelengthHydrogen_atomic) then
       energyPhoton                                      =+plancksConstant                                            &
            &                                             *speedLight                                                 &
            &                                             *metersToAngstroms                                          &
            &                                             /photonPacket%wavelength                      (           )
       rateEscape                                        =+photonPacket%luminosity                      (           ) &
            &                                             *luminositySolar                                            &
            &                                             /energyPhoton
       self%lymanContinuumRateEscaping                   =+self        %lymanContinuumRateEscaping                    &
            &                                             +rateEscape
       sourceIndex                                       =+photonPacket%sourceType()
       self%lymanContinuumRateEscapingTagged(sourceIndex)=+self        %lymanContinuumRateEscapingTagged(sourceIndex) &
            &                                             +rateEscape
    end if
    return
  end subroutine lymanContinuumRatePhotonPacketEscapes

  subroutine lymanContinuumRateFinalize(self)
    !!{
    Finalize the Lyman continuum photon escape rate.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class(radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self

    ! Sum the Lyc rate across all MPI processes.
    self%lymanContinuumRateEscaping      =mpiSelf%sum(self%lymanContinuumRateEscaping      )
    self%lymanContinuumRateEscapingTagged=mpiSelf%sum(self%lymanContinuumRateEscapingTagged)
    return
  end subroutine lymanContinuumRateFinalize

  subroutine lymanContinuumRateOutput(self,outputGroup)
    !!{
    Output the Lyman continuum photon escape rate.
    !!}
    use :: HDF5_Access    , only : hdf5Access
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    class  (radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    type   (hdf5Object                                  ), intent(inout) :: outputGroup
    integer                                                              :: sourceIndex

    !$ call hdf5Access%set  ()
    call outputGroup%writeAttribute(self%lymanContinuumRateEscaping,'rateLymanContinuumEscaping')
    !$ call hdf5Access%unset()
    do sourceIndex=1,size(self%lymanContinuumRateEscapingTagged)
       !$ call hdf5Access%set  ()
       call outputGroup%writeAttribute(self%lymanContinuumRateEscapingTagged(sourceIndex),'rateLymanContinuumEscaping'//String_Upper_Case_First(char(self%sourceTypeName(sourceIndex))))
       !$ call hdf5Access%unset()
    end do
    return
  end subroutine lymanContinuumRateOutput
