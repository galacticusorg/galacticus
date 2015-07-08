!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  use FGSL

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
     procedure :: descriptor => fileDescriptor
  end type accretionDiskSpectraFile

  interface accretionDiskSpectraFile
     !% Constructors for the {\normalfont \ttfamily file} accretion disk spectra class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface accretionDiskSpectraFile

  ! Supported file format.
  integer :: fileFormatCurrent=1

contains

  function fileConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily file} accretion disk spectra class which takes a parameter set as input.
    implicit none
    type(accretionDiskSpectraFile)                :: fileConstructorParameters
    type(inputParameters         ), intent(in   ) :: parameters
    type(varying_string          )                :: fileName
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of a file from which to read tabulated spectra of accretion disks.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    fileConstructorParameters=fileConstructorInternal(char(fileName))
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName)
    !% Internal constructor for the {\normalfont \ttfamily file} accretion disk spectra class.
    use Galacticus_Error
    use Array_Utilities
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
            &  'fileConstructorInternal'                                                                                             , &
            &  'This method requires that the "radiativeEfficiency", and "accretionRate" properties of the black hole are gettable.'// &
            &  Galacticus_Component_List(                                                                                              &
            &                            'blackHole'                                                                                ,  &
            &                             defaultBlackHoleComponent%radiativeEfficiencyAttributeMatch(requireGettable=.true.)          &
            &                            .intersection.                                                                                &
            &                             defaultBlackHoleComponent%      accretionRateAttributeMatch(requireGettable=.true.)          &
            &                           )                                                                                              &
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

    if (allocated(self%wavelength)) call Dealloc_Array(self%wavelength)
    if (allocated(self%luminosity)) call Dealloc_Array(self%luminosity)
    if (allocated(self%SED       )) call Dealloc_Array(self%SED       )
    call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorWavelength,reset=self%resetWavelength)
    call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorLuminosity,reset=self%resetLuminosity)
    return
  end subroutine fileDestructor

  subroutine fileLoadFile(self,fileName)
    !% Load a file of AGN spectra.
    use IO_HDF5
    use Galacticus_Error
    implicit none
    class    (accretionDiskSpectraFile), intent(inout) :: self
    character(len=*                   ), intent(in   ) :: fileName
    integer                                            :: fileFormatFile
    type     (hdf5Object              )                :: spectraFile

    ! Open the file.
    !$omp critical(HDF5_Access)
    call spectraFile%openFile(fileName,readOnly=.true.)
    ! Check file format.
    call spectraFile%readAttribute('fileFormat',fileFormatFile,allowPseudoScalar=.true.)
    if (fileFormatFile /= fileFormatCurrent) call Galacticus_Error_Report('fileLoadFile','file format mismatch')
    ! Read datasets.
    call spectraFile%readDataset('wavelength'          ,self%wavelength)
    call spectraFile%readDataset('bolometricLuminosity',self%luminosity)
    call spectraFile%readDataset('SED'                 ,self%SED       )
    ! Close the file.
    call spectraFile%close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine fileLoadFile

  double precision function fileSpectrum(self,node,wavelength)
    !% Return the accretion disk spectrum for tabulated spectra.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
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
    iLuminosity=Interpolate_Locate                 (self%luminosity,self%interpolationAcceleratorLuminosity,log10(luminosityBolometric),self%resetLuminosity)
    iWavelength=Interpolate_Locate                 (self%wavelength,self%interpolationAcceleratorWavelength,      wavelength           ,self%resetWavelength)
    hLuminosity=Interpolate_Linear_Generate_Factors(self%luminosity,iLuminosity                            ,log10(luminosityBolometric)                     )
    hWavelength=Interpolate_Linear_Generate_Factors(self%wavelength,iWavelength                            ,      wavelength                                )
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
  
  subroutine fileDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(accretionDiskSpectraFile), intent(inout) :: self
    type (inputParameters         ), intent(inout) :: descriptor
    type (node                    ), pointer       :: parameterNode
    type (inputParameters         )                :: subParameters

    call descriptor%addParameter("accretionDiskSpectraMethod","file")
    parameterNode => descriptor%node("accretionDiskSpectraMethod")
    subParameters=inputParameters(parameterNode)
    call subParameters%addParameter("fileName",char(self%fileName))
    return
  end subroutine fileDescriptor
