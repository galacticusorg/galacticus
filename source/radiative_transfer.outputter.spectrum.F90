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

  use, intrinsic :: ISO_C_Binding, only : c_size_t
  
  !![
  <radiativeTransferOutputter name="radiativeTransferOutputterSpectrum">
   <description>A radiative transfer outputter class which outputs a binned spectrum of both total emitted and emergent photons.</description>
  </radiativeTransferOutputter>
  !!]
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterSpectrum
     !!{
     Implementation of a radiative transfer outputter class which outputs a binned spectrum of both total emitted and emergent photons.
     !!}
     private
     integer         (c_size_t     )                                 :: wavelengthCountPerDecade, countWavelengths  , &
          &                                                             countThetas
     double precision                                                :: wavelengthMinimum       , wavelengthMaximum , &
          &                                                             thetaMinimum            , thetaMaximum
     double precision                , allocatable, dimension(:    ) :: wavelengths             , thetas            , &
          &                                                             wavelengthsMinimum      , wavelengthsMaximum, &
          &                                                             thetasMinimum           , thetasMaximum     , &
          &                                                             solidAngles             , energies
     double precision                , allocatable, dimension(:,:  ) :: spectrumEmergent
     double precision                , allocatable, dimension(:,:,:) :: spectrumEmergentTagged
     type            (varying_string), allocatable, dimension(:    ) :: sourceTypeName
   contains
     procedure :: reset               => spectrumReset
     procedure :: sourceProperties    => spectrumSourceProperties
     procedure :: photonPacketEscapes => spectrumPhotonPacketEscapes
     procedure :: finalize            => spectrumFinalize
     procedure :: output              => spectrumOutput
  end type radiativeTransferOutputterSpectrum

  interface radiativeTransferOutputterSpectrum
     !!{
     Constructors for the \refClass{radiativeTransferOutputterSpectrum} radiative transfer outputter class.
     !!}
     module procedure spectrumConstructorParameters
     module procedure spectrumConstructorInternal
  end interface radiativeTransferOutputterSpectrum
  
contains

  function spectrumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferOutputterSpectrum} radiative transfer outputter class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters        , only : inputParameters
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (radiativeTransferOutputterSpectrum)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    integer         (c_size_t                          )                :: wavelengthCountPerDecade, countThetas
    double precision                                                    :: wavelengthMinimum       , wavelengthMaximum, &
         &                                                                 thetaMinimum            , thetaMaximum

    !![
    <inputParameter>
      <name>wavelengthMinimum</name>
      <defaultValue>0.3d4</defaultValue>
      <description>The minimum wavelength at which to compute spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <defaultValue>10.0d4</defaultValue>
      <description>The maximum wavelength at which to compute spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthCountPerDecade</name>
      <defaultValue>10_c_size_t</defaultValue>
      <description>The number of wavelengths per decade at which to compute spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>thetaMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum angle $\theta$ at which to bin the emergent spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>thetaMaximum</name>
      <defaultValue>+Pi</defaultValue>
      <description>The maximum angle $\theta$ at which to bin the emergent spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countThetas</name>
      <defaultValue>1_c_size_t</defaultValue>
      <description>The number of bins in angle $\theta$ at which to compute spectra.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=radiativeTransferOutputterSpectrum(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,thetaMinimum,thetaMaximum,countThetas)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function spectrumConstructorParameters

  function spectrumConstructorInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,thetaMinimum,thetaMaximum,countThetas) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferOutputterSpectrum} radiative transfer outputter class.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : plancksConstant  , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms, electronVolt
    use :: Numerical_Ranges            , only : Make_Range       , rangeTypeLogarithmic, rangeTypeLinear
    implicit none
    type            (radiativeTransferOutputterSpectrum)                :: self
    integer         (c_size_t                          ), intent(in   ) :: wavelengthCountPerDecade, countThetas
    double precision                                    , intent(in   ) :: wavelengthMinimum       , wavelengthMaximum, &
         &                                                                 thetaMinimum            , thetaMaximum
    !![
    <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, wavelengthCountPerDecade, thetaMinimum, thetaMaximum, countThetas"/>
    !!]

    self%countWavelengths=int(log10(wavelengthMaximum/wavelengthMinimum)*dble(wavelengthCountPerDecade),c_size_t)+1_c_size_t
    allocate(self%wavelengths       (self%countWavelengths                 ))
    allocate(self%energies          (self%countWavelengths                 ))
    allocate(self%wavelengthsMinimum(self%countWavelengths                 ))
    allocate(self%wavelengthsMaximum(self%countWavelengths                 ))
    allocate(self%solidAngles       (                      self%countThetas))
    allocate(self%thetas            (                      self%countThetas))
    allocate(self%thetasMinimum     (                      self%countThetas))
    allocate(self%thetasMaximum     (                      self%countThetas))
    allocate(self%spectrumEmergent  (self%countWavelengths,self%countThetas))
    self%wavelengths       =Make_Range(wavelengthMinimum,wavelengthMaximum,int(self%countWavelengths),rangeTypeLogarithmic,rangeBinned=.true.)
    self%thetas            =Make_Range(thetaMinimum     ,thetaMaximum     ,int(self%countThetas     ),rangeTypeLinear     ,rangeBinned=.true.)
    self%wavelengthsMinimum=self%wavelengths/sqrt(self%wavelengths(2)/self%wavelengths(1))
    self%wavelengthsMaximum=self%wavelengths*sqrt(self%wavelengths(2)/self%wavelengths(1))
    self%thetasMinimum     =self%thetas     -0.5d0*(thetaMaximum-thetaMinimum)/dble(countThetas)
    self%thetasMaximum     =self%thetas     +0.5d0*(thetaMaximum-thetaMinimum)/dble(countThetas)
    self%spectrumEmergent  =0.0d0
    ! Compute energies in electronvolts.
    self%energies   =+plancksConstant           &
         &           *speedLight                &
         &           *metersToAngstroms         &
         &           /self%wavelengths          &
         &           /electronVolt
    ! Compute the solid angle associated with each theta bin.
    self%solidAngles=+2.0d0                     &
         &           *Pi                        &
         &           *(                         &
         &             -cos(self%thetasMaximum) &
         &             +cos(self%thetasMinimum) &
         &            )
    return
  end function spectrumConstructorInternal
  
  subroutine spectrumReset(self)
    !!{
    Reset the accumulated spectrum.
    !!}
    implicit none
    class(radiativeTransferOutputterSpectrum), intent(inout) :: self

    self%spectrumEmergent=0.0d0
    if (allocated(self%spectrumEmergentTagged)) self%spectrumEmergentTagged=0.0d0
    return
  end subroutine spectrumReset

  subroutine spectrumSourceProperties(self,radiativeTransferSource_,outputGroup)
    !!{
    Compute and output the emission spectrum.
    !!}
    use :: HDF5_Access                     , only : hdf5Access
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    use :: Numerical_Integration           , only : integrator
    use :: MPI_Utilities                   , only : mpiSelf
    use :: String_Handling                 , only : String_Upper_Case_First
    implicit none
    class           (radiativeTransferOutputterSpectrum), intent(inout)                    :: self
    class           (radiativeTransferSourceClass      ), intent(inout)                    :: radiativeTransferSource_
    type            (hdf5Object                        ), intent(inout)                    :: outputGroup
    type            (integrator                        )                                   :: integrator_
    double precision                                    , dimension(self%countWavelengths) :: spectrumEmitted
    integer         (c_size_t                          )                                   :: i
    integer                                                                                :: sourceIndex
    type            (varying_string                    )                                   :: label                   , description

    if (.not.allocated(self%spectrumEmergentTagged)) then
       allocate(self%spectrumEmergentTagged(self%countWavelengths,self%countThetas,radiativeTransferSource_%sourceTypeCount()))
       allocate(self%sourceTypeName        (                                       radiativeTransferSource_%sourceTypeCount()))
       self%spectrumEmergentTagged=0.0d0
       do sourceIndex=1,radiativeTransferSource_%sourceTypeCount()
          self%sourceTypeName(sourceIndex)=radiativeTransferSource_%sourceTypeName(sourceIndex)
       end do
    end if
    if (mpiSelf%isMaster()) then
       !$ call hdf5Access%set  ()
       call outputGroup%writeDataset(self%wavelengths,'spectrumWavelength','Central wavelengths of spectral bins' )
       call outputGroup%writeDataset(self%energies   ,'spectrumEnergies'  ,'Central energies of spectral bins'    )
       call outputGroup%writeDataset(self%thetas     ,'spectrumTheta'     ,'Central angles theta of spectral bins')
       call outputGroup%writeDataset(self%solidAngles,'spectrumSolidAngle','Solid angles of spectral bins'        )
       !$ call hdf5Access%unset()
       integrator_=integrator(integrand,toleranceRelative=1.0d-2)
       do sourceIndex=0,radiativeTransferSource_%sourceTypeCount()
          do i=1_c_size_t,self%countWavelengths
             spectrumEmitted(i)=+integrator_%integrate(                            &
                  &                                    self%wavelengthsMinimum(i), &
                  &                                    self%wavelengthsMaximum(i)  &
                  &                                   )                            &
                  &             *luminositySolar                                   &
                  &             /(                                                 &
                  &               +                    self%wavelengthsMaximum(i)  &
                  &               -                    self%wavelengthsMinimum(i)  &
                  &              )
          end do
          label      =var_str('spectrumEmitted'       )
          if (sourceIndex > 0) label      =label      //      String_Upper_Case_First(char(radiativeTransferSource_%sourceTypeName(sourceIndex)))
          description=var_str('Total emitted spectrum')
          if (sourceIndex > 0) description=description//" ["//String_Upper_Case_First(char(radiativeTransferSource_%sourceTypeName(sourceIndex)))//"]"
          !$ call hdf5Access%set  ()
          call outputGroup%writeDataset(spectrumEmitted,char(label),'Total emitted spectrum'               )
          !$ call hdf5Access%unset()
       end do
    end if
    return
    
  contains

    double precision function integrand(wavelength)
      !!{
      Integrand over the source spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength
      
      integrand=radiativeTransferSource_%spectrum(wavelength,sourceIndex)
      return
    end function integrand

  end subroutine spectrumSourceProperties

  subroutine spectrumPhotonPacketEscapes(self,photonPacket)
    !!{
    Process an escaping photon packet.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (radiativeTransferOutputterSpectrum), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    integer         (c_size_t                          )                :: iWavelength , iTheta    , &
         &                                                                 sourceIndex
    double precision                                                    :: photonTheta

    photonTheta=acos(Dot_Product(photonPacket%direction(),[0.0d0,0.0d0,1.0d0]))
    if     (                                                                                                              &
         &   photonPacket%wavelength() >= self%wavelengthMinimum .and. photonPacket%wavelength() < self%wavelengthMaximum &
         &  .and.                                                                                                         &
         &   photonTheta               >= self%thetaMinimum      .and. photonTheta               < self%thetaMaximum      &
         & ) then
       if (photonPacket%wavelength() > self%wavelengthsMinimum(self%countWavelengths)) then
          iWavelength=self%countWavelengths
       else
          iWavelength                              =searchArray(self%wavelengthsMinimum,photonPacket%wavelength())
       end if
       if (photonTheta               > self%thetasMinimum      (self%countThetas    )) then
          iTheta     =self%countThetas
       else
          iTheta                                   =searchArray(self%thetasMinimum     ,photonTheta              )
       end if
       sourceIndex                                                =+photonPacket%sourceType            (                              )
       self%spectrumEmergent      (iWavelength,iTheta            )=+self        %spectrumEmergent      (iWavelength,iTheta            ) &
            &                                                      +photonPacket%luminosity            (                              )
       self%spectrumEmergentTagged(iWavelength,iTheta,sourceIndex)=+self        %spectrumEmergentTagged(iWavelength,iTheta,sourceIndex) &
            &                                                      +photonPacket%luminosity            (                              )
    end if
    return
  end subroutine spectrumPhotonPacketEscapes

  subroutine spectrumFinalize(self)
    !!{
    Finalize the spectrum.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class(radiativeTransferOutputterSpectrum), intent(inout) :: self

    ! Sum the emergent spectrum across all MPI processes.
    self%spectrumEmergent      =mpiSelf%sum(self%spectrumEmergent      )
    self%spectrumEmergentTagged=mpiSelf%sum(self%spectrumEmergentTagged)
    return
  end subroutine spectrumFinalize

  subroutine spectrumOutput(self,outputGroup)
    !!{
    Output the spectrum.
    !!}
    !$ use :: HDF5_Access                     , only : hdf5Access
    use    :: Numerical_Constants_Astronomical, only : luminositySolar
    use    :: String_Handling                 , only : String_Upper_Case_First
    implicit none
    class  (radiativeTransferOutputterSpectrum), intent(inout) :: self
    type   (hdf5Object                        ), intent(inout) :: outputGroup
    integer(c_size_t                          )                :: i          , sourceIndex

    ! Apply unit conversions and make the spectrum differential with respect to wavelength.
    do i=1_c_size_t,self%countWavelengths
       self%spectrumEmergent      (i,:  )=+  self%spectrumEmergent      (i,:  ) &
            &                             *       luminositySolar               &
            &                             /(                                    &
            &                               +self%wavelengthsMaximum    (i    ) &
            &                               -self%wavelengthsMinimum    (i    ) &
            &                              )
       self%spectrumEmergentTagged(i,:,:)=+  self%spectrumEmergentTagged(i,:,:) &
            &                             *       luminositySolar               &
            &                             /(                                    &
            &                               +self%wavelengthsMaximum    (i    ) &
            &                               -self%wavelengthsMinimum    (i    ) &
            &                              )
    end do
    !$ call hdf5Access%set  ()
    call outputGroup%writeDataset(self%spectrumEmergent,'spectrumEmergent','Total emergent spectrum')
    do sourceIndex=1,size(self%spectrumEmergentTagged,dim=3)
       call outputGroup%writeDataset(self%spectrumEmergentTagged(:,:,sourceIndex),'spectrumEmergent'//String_Upper_Case_First(char(self%sourceTypeName(sourceIndex))),'Total emergent spectrum ['//String_Upper_Case_First(char(self%sourceTypeName(sourceIndex)))//']')
    end do
    !$ call hdf5Access%unset()
    return
  end subroutine spectrumOutput
