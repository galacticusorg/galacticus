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

  !% An implementation of the accretion disk spectra class for tabulated spectra read from file.

  !# <accretionDiskSpectra name="accretionDiskSpectraFile">
  !#  <description>Accretion disk spectra are interpolated from tables read from file.</description>
  !# </accretionDiskSpectra>

  use FGSL, only : fgsl_interp_accel

  type, extends(accretionDiskSpectraClass) :: accretionDiskSpectraFile
     !% An accretion disk spectra class which interpolates in spectra read from file.
     private
     type            (varying_string   )                              :: fileName
     double precision                   , allocatable, dimension(:  ) :: luminosity                        , wavelength
     double precision                   , allocatable, dimension(:,:) :: SED
     type            (fgsl_interp_accel)                              :: interpolationAcceleratorLuminosity, interpolationAcceleratorWavelength
     logical                                                          :: resetLuminosity                   , resetWavelength
   contains
     !@ <objectMethods>
     !@   <object>accretionDiskSpectraFile</object>
     !@   <objectMethod>
     !@     <method>loadFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Load a file of AGN spectra.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::               fileDestructor
     procedure :: spectrum   => fileSpectrum
     procedure :: loadFile   => fileLoadFile
  end type accretionDiskSpectraFile

  interface accretionDiskSpectraFile
     !% Constructors for the {\normalfont \ttfamily file} accretion disk spectra class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface accretionDiskSpectraFile

  ! Supported file format.
  integer :: fileFormatCurrent=2

contains

  function fileConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily file} accretion disk spectra class which takes a parameter set as input.
    implicit none
    type(accretionDiskSpectraFile)                :: fileConstructorParameters
    type(inputParameters         ), intent(inout) :: parameters
    type(varying_string          )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of a file from which to read tabulated spectra of accretion disks.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    fileConstructorParameters=fileConstructorInternal(char(fileName))
    !# <inputParametersValidate source="parameters"/>
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName)
    !% Internal constructor for the {\normalfont \ttfamily file} accretion disk spectra class.
    use Galacticus_Error, only : Galacticus_Error_Report  , Galacticus_Component_List
    use Array_Utilities
    use Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none
    type     (accretionDiskSpectraFile), target        :: fileConstructorInternal
    character(len=*                   ), intent(in   ) :: fileName
 
    ! Ensure that the required methods are supported.
    if     (                                                                                                                           &
            &  .not.(                                                                                                                  &
            &         defaultBlackHoleComponent%radiativeEfficiencyIsGettable()                                                        &
            &        .and.                                                                                                             &
            &         defaultBlackHoleComponent%      accretionRateIsGettable()                                                        &
            &       )                                                                                                                  &
            & ) call Galacticus_Error_Report                                                                                           &
            & (                                                                                                                        &
            &  'This method requires that the "radiativeEfficiency", and "accretionRate" properties of the black hole are gettable.'// &
            &  Galacticus_Component_List(                                                                                              &
            &                            'blackHole'                                                                                ,  &
            &                             defaultBlackHoleComponent%radiativeEfficiencyAttributeMatch(requireGettable=.true.)          &
            &                            .intersection.                                                                                &
            &                             defaultBlackHoleComponent%      accretionRateAttributeMatch(requireGettable=.true.)          &
            &                           )                                                                                           // &
            &  {introspection:location}                                                                                                &
            & )
    ! Load the file.
    fileConstructorInternal%fileName=fileName
    call fileConstructorInternal%loadFile(fileName)
    ! Initialize interpolators.
    fileConstructorInternal%resetLuminosity=.true.
    fileConstructorInternal%resetWavelength=.true.
    return
  end function fileConstructorInternal

  subroutine fileDestructor(self)
    !% Default destructor for the {\normalfont \ttfamily file} accretion disk spectra class.
    use Numerical_Interpolation
    use Memory_Management
    implicit none
    type(accretionDiskSpectraFile), intent(inout) :: self

    if (allocated(self%wavelength)) then
       call deallocateArray(self%wavelength)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorWavelength,reset=self%resetWavelength)
    end if
    if (allocated(self%luminosity)) then
       call deallocateArray(self%luminosity)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorLuminosity,reset=self%resetLuminosity)
    end if
    if (allocated(self%SED       )) call deallocateArray(self%SED)
    return
  end subroutine fileDestructor

  subroutine fileLoadFile(self,fileName)
    !% Load a file of AGN spectra.
    use IO_HDF5
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (accretionDiskSpectraFile), intent(inout) :: self
    character(len=*                   ), intent(in   ) :: fileName
    integer                                            :: fileFormatFile
    type     (hdf5Object              )                :: spectraFile

    ! Open the file.
    !$ call hdf5Access%set()
    call spectraFile%openFile(fileName,readOnly=.true.)
    ! Check file format.
    call spectraFile%readAttribute('fileFormat',fileFormatFile,allowPseudoScalar=.true.)
    if (fileFormatFile /= fileFormatCurrent) call Galacticus_Error_Report('file format mismatch'//{introspection:location})
    ! Read datasets.
    call spectraFile%readDataset('wavelength'          ,self%wavelength)
    call spectraFile%readDataset('bolometricLuminosity',self%luminosity)
    call spectraFile%readDataset('SED'                 ,self%SED       )
    ! Close the file.
    call spectraFile%close()
    !$ call hdf5Access%unset()
    ! Convert luminosities to logarithmic form for interpolation.
    self%luminosity=log(self%luminosity)
    return
  end subroutine fileLoadFile

  double precision function fileSpectrum(self,node,wavelength)
    !% Return the accretion disk spectrum for tabulated spectra.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Galacticus_Nodes                , only : nodeComponentBlackHole
    implicit none
    class           (accretionDiskSpectraFile), intent(inout)  :: self
    type            (treeNode                ), intent(inout)  :: node
    double precision                          , intent(in   )  :: wavelength
    class           (nodeComponentBlackHole  ), pointer        :: blackHole
    double precision                          , dimension(0:1) :: hLuminosity         , hWavelength
    integer         (c_size_t                )                 :: iLuminosity         , iWavelength, &
         &                                                        jLuminosity         , jWavelength
    double precision                                           :: luminosityBolometric

    ! Initialize to zero spectrum.
    fileSpectrum=0.0d0
    ! Get the bolometric luminosity.
    blackHole => node%blackHole()
    luminosityBolometric=blackHole%radiativeEfficiency()*blackHole%accretionRate()*massSolar*speedLight**2/gigaYear/luminositySolar
    ! Return on non-positive luminosities.
    if (luminosityBolometric <= 0.0d0) return
    ! Assume zero flux outside of the tabulated wavelength range.
    if     (                                                     &
         &   wavelength < self%wavelength(                   1 ) &
         & .or.                                                  &
         &   wavelength > self%wavelength(size(self%wavelength)) &
         & ) return
    ! Get the interpolating factors.
    iLuminosity=Interpolate_Locate                 (self%luminosity,self%interpolationAcceleratorLuminosity,log(luminosityBolometric),self%resetLuminosity)
    iWavelength=Interpolate_Locate                 (self%wavelength,self%interpolationAcceleratorWavelength,    wavelength           ,self%resetWavelength)
    hLuminosity=Interpolate_Linear_Generate_Factors(self%luminosity,iLuminosity                            ,log(luminosityBolometric)                     )
    hWavelength=Interpolate_Linear_Generate_Factors(self%wavelength,iWavelength                            ,    wavelength                                )
    ! Do the interpolation.
    do jLuminosity=0,1
       do jWavelength=0,1
          fileSpectrum=                                           &
               &       +fileSpectrum                              &
               &       +self        %SED(                         &
               &                         iWavelength+jWavelength, &
               &                         iLuminosity+jLuminosity  &
               &                        )                         &
               &       *hLuminosity     (            jLuminosity) &
               &       *hWavelength     (            jWavelength)
       end do
    end do
    ! Prevent interpolation from returning negative fluxes.
    fileSpectrum=max(fileSpectrum,0.0d0)
    return
  end function fileSpectrum
  
