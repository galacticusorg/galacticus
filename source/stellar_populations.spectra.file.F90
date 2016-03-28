!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Implements a file-based stellar population spectra class.

  use FGSL
  
  type spectralTable
     !% Structure to hold tabulated stellar population  data.
     ! The spectra tables.
     integer                                                            :: agesCount                     , metallicityCount              , &
          &                                                                wavelengthsCount
     double precision                   , allocatable, dimension(:    ) :: ages                          , metallicities                 , &
          &                                                                wavelengths
     double precision                   , allocatable, dimension(:,:,:) :: table
     ! Interpolation structures.
     logical                                                            :: resetAge               =.true., resetMetallicity       =.true., &
          &                                                                resetWavelength        =.true.
     type            (fgsl_interp_accel)                                :: acceleratorAge                , acceleratorMetallicity        , &
          &                                                                acceleratorWavelength
   contains
     final :: spectralTableDestructorScalar, spectralTableDestructor1D
  end type spectralTable

  !# <stellarPopulationSpectra name="stellarPopulationSpectraFile">
  !#  <description>
  !#   Provides stellar population spectra via interpolation in a tabulation read from file. File names are specified via the parameters {\normalfont \ttfamily [fileNameForXXXXIMF]} parameter where {\normalfont \ttfamily XXXX} is the name of the \gls{imf}. This should specify an HDF5 file with the following structure:
  !#   \begin{verbatim}
  !#    ages                     Dataset {ageCount}
  !#    metallicities            Dataset {metallicityCount}
  !#    spectra                  Dataset {metallicityCount, ageCount, metallicityCount}
  !#    wavelengths              Dataset {wavelengthCount}
  !#   \end{verbatim}
  !#   where the datasets contain the tabulated ages (in Gyr), metallicities (logarithmic, relative to Solar), wavelengths (in \AA) and spectra (in $L_\odot$ Hz$^{-1}$).
  !#
  !#   Currently, the following pre-computed stellar spectra files are available as a separate download from \href{http://users.obs.carnegiescience.edu/abenson/galacticus/data/Galacticus_SSP_Data.tar.bz2}{\normalfont \ttfamily http://users.obs.carnegiescience.edu/abenson/galacticus/data/Galacticus\_SSP\_Data.tar.bz2}:
  !#   \begin{description}
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.0\_imfSalpeter.hdf5}] Corresponds to a Salpeter IMF computed using v2.0 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.1\_imfSalpeter.hdf5}]  Corresponds to a Salpeter IMF computed using v2.1 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.1\_imfChabrier.hdf5}]  Corresponds to a Chabrier IMF computed using v2.1 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.2\_imfChabrier.hdf5}]  Corresponds to a Chabrier IMF computed using v2.2 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.2\_imfKennicutt.hdf5}]  Corresponds to a Kennicutt IMF computed using v2.2 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Conroy-et-al\_v2.2\_imfBaugh2005TopHeavy.hdf5}]  Corresponds to the top-heavy IMF of \cite{baugh_can_2005} computed using v2.2 of the {\normalfont \ttfamily FSPS} code;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Maraston\_hbMorphologyRed\_imfKroupa.hdf5}] The spectra from \cite{maraston_evolutionary_2005} for a Kroupa IMF and a red horizontal branch morphology;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Maraston\_hbMorphologyRed\_imfSalpeter.hdf5}] The spectra from \cite{maraston_evolutionary_2005} for a Salpeter IMF and a red horizontal branch morphology; 
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_BC2003\_highResolution\_imfChabrier.hdf5}] The (high resolution) spectra from \cite{bruzual_stellar_2003} for a Chabrier IMF, using Padova 1994 tracks;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_BC2003\_highResolution\_imfSalpeter.hdf5}] The (high resolution) spectra from \cite{bruzual_stellar_2003} for a Salpeter IMF, using Padova 1994 tracks;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_BC2003\_lowResolution\_imfChabrier.hdf5}] The (low resolution) spectra from \cite{bruzual_stellar_2003} for a Chabrier IMF, using Padova 1994 tracks;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_BC2003\_lowResolution\_imfSalpeter.hdf5}] The (low resolution) spectra from \cite{bruzual_stellar_2003} for a Salpeter IMF, using Padova 1994 tracks;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Grasil\_gkn15rd\_ken.hdf5}] The spectra used by {\normalfont \scshape Grasil} for a Kennicutt IMF;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Grasil\_gkn1rd\_ken.hdf5}] The spectra used by {\normalfont \scshape Grasil} for a Kennicutt IMF;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Grasil\_gsrdk0b\_sal.hdf5}] The spectra used by {\normalfont \scshape Grasil} for a Salpeter IMF;
  !#    \item [{\normalfont \ttfamily data/stellarPopulations/SSP\_Spectra\_Grasil\_imf27\_kro.hdf5}] Spectra used by {\normalfont \scshape Grasil}.
  !#   \end{description}
  !#   Note that the high resolution spectra from \cite{bruzual_stellar_2003} may require you to adjust the {\normalfont \ttfamily [stellarPopulationLuminosityIntegrationToleranceRelative]} parameter to a larger value\footnote{Or, alternatively, set {\normalfont \ttfamily [stellarPopulationLuminosityIntegrationToleranceDegrade]}$=${\normalfont \ttfamily true}. This will cause \glc\ to increase the tolerance as necessary to get the integrals to converge---issuing warnings each time the tolerance is increased.}. The sharp features in these high resolution spectra can be difficult to integrate. Scripts to convert the data provided by \cite{maraston_evolutionary_2005} and \cite{bruzual_stellar_2003} into \glc's format are provided in the {\normalfont \ttfamily scripts/ssps} folder. Spectra for other initial mass functions will be computed automatically when using the \cite{conroy_propagation_2009} population synthesis models.
  !#  </description>
  !# </stellarPopulationSpectra>
  type, extends(stellarPopulationSpectraClass) :: stellarPopulationSpectraFile
     !% A stellar population spectra class which interpolates spectra given in a file.
     private
     type   (spectralTable ), allocatable, dimension(:) :: spectra              ! Spectral data tables.
     integer                , allocatable, dimension(:) :: imfLookup            ! Look-up array to cross-reference IMF indices to our internal data structure.
     logical                                            :: forceZeroMetallicity ! Interpolation option.
     type   (varying_string), allocatable, dimension(:) :: fileName             ! File names per IMF.
   contains
     !@ <objectMethods>
     !@   <object>stellarPopulationSpectraFile</object>
     !@   <objectMethod>
     !@     <method>readFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless char(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the named stellar population spectra file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>imfInitialize</method>
     !@     <type>void</type>
     !@     <arguments>\intzero\ imfIndex\argin</arguments>
     !@     <description>Ensure that spectra are available for the specified \gls{imf} index.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                  fileDestructor
     procedure :: readFile      => fileReadFile
     procedure :: luminosity    => fileLuminosity
     procedure :: tabulation    => fileTabulation
     procedure :: imfInitialize => fileIMFInitialize
     procedure :: descriptor    => fileDescriptor
  end type stellarPopulationSpectraFile

  interface stellarPopulationSpectraFile
     !% Constructors for the file stellar spectra class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface stellarPopulationSpectraFile

  ! The current file format version.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters)
    !% Constructor for the file stellar spectra class which takes a parameter set as input.
    use Star_Formation_IMF
    use Galacticus_Input_Paths
    implicit none
    type   (stellarPopulationSpectraFile)                             :: fileConstructorParameters
    type   (inputParameters             ), intent(inout)              :: parameters
    integer                              , dimension(:) , allocatable :: imfIndices
    type   (varying_string              ), dimension(:) , allocatable :: fileNames
    logical                                                           :: forceZeroMetallicity
    integer                                                           :: i                        , imfCount
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>forceZeroMetallicity</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Force the use of zero metallicity (or lowest metallicity available) for all stellar populations.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Allocate storage for IMF filenames.
    imfCount=IMF_Available_Count()
    allocate(imfIndices(imfCount))
    allocate(fileNames (imfCount))
    forall(i=1:imfCount)
       imfIndices(i)=i
    end forall
    ! Get all file names
    !# <inputParameter>
    !#   <iterator>fileNameFor(#imfRegisterName->name)IMF</iterator>
    !#   <source>parameters</source>
    !#   <variable>fileNames(IMF_Index("$1"))</variable>
    !#   <defaultValue>Galacticus_Input_Path()//"data/SSP_Spectra_imf$1.hdf5"</defaultValue>
    !#   <description>The name of the file of stellar populations to use for the named \gls{imf}.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Construct the instance.    
    fileConstructorParameters=fileConstructorInternal(forceZeroMetallicity,imfIndices,fileNames)
    return
  end function fileConstructorParameters
  
  function fileConstructorInternal(forceZeroMetallicity,imfIndices,fileNames)
    !% Internal constructor for the file stellar spectra class.
    use Star_Formation_IMF
    use Galacticus_Error
    implicit none
    type   (stellarPopulationSpectraFile)                              :: fileConstructorInternal
    logical                              , intent(in   )               :: forceZeroMetallicity
    integer                              , intent(in   ), dimension(:) :: imfIndices
    type   (varying_string              ), intent(in   ), dimension(:) :: fileNames
    integer                                                            :: imfCount               , i
    
    ! Store options.
    fileConstructorInternal%forceZeroMetallicity=forceZeroMetallicity
    ! Allocate storage for IMF file names.
    imfCount=IMF_Available_Count()
    allocate(fileConstructorInternal%fileName(imfCount))
    fileConstructorInternal%fileName="?"
    do i=1,size(imfIndices)
       if (imfIndices(i) < 1 .or. imfIndices(i) > imfCount) then
          call Galacticus_Error_Report('fileConstructorInternal','IMF index is out of range')
       else
          fileConstructorInternal%fileName(imfIndices(i))=fileNames(i)
       end if
    end do
    return
  end function fileConstructorInternal
  
  subroutine fileDestructor(self)
    !% Destructor for the file stellar spectra class.
    implicit none
    type(stellarPopulationSpectraFile), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine fileDestructor

  subroutine spectralTableDestructorScalar(self)
    !% Destructor for the stellar spectra table class.
    use Numerical_Interpolation
    implicit none
    type(spectralTable), intent(inout) :: self

    ! Free all FGSL objects.
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorAge        , &
         &                reset                   =self%      resetAge          &
         &               )
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorWavelength , &
         &                reset                   =self%      resetWavelength   &
         &               )
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorMetallicity, &
         &                reset                   =self%      resetMetallicity  &
         &               )
    return
  end subroutine spectralTableDestructorScalar

  subroutine spectralTableDestructor1D(self)
    !% Destructor for the stellar spectra table class.
    use Numerical_Interpolation
    implicit none
    type   (spectralTable), intent(inout), dimension(:) :: self
    integer                                             :: i

    do i=1,size(self)
       call spectralTableDestructorScalar(self(i))
    end do
    return
  end subroutine spectralTableDestructor1D
  
  double precision function fileLuminosity(self,abundancesStellar,age,wavelength,imfIndex,status)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\normalfont \ttfamily abundances}, of the
    !% given {\normalfont \ttfamily age} (in Gyr) and the specified {\normalfont \ttfamily wavelength} (in Angstroms). This is found by interpolating in tabulated
    !% spectra.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Numerical_Interpolation
    use               Galacticus_Error
    implicit none
    class           (stellarPopulationSpectraFile), intent(inout)            :: self
    type            (abundances                  ), intent(in   )            :: abundancesStellar
    double precision                              , intent(in   )            :: age                        , wavelength
    integer                                       , intent(in   )            :: imfIndex
    integer                                       , intent(  out) , optional :: status
    double precision                              , dimension(0:1)           :: hAge                       , hMetallicity  , &
         &                                                                      hWavelength
    double precision                              , parameter                :: metallicityTolerance=0.01d0
    integer         (c_size_t                    )                           :: iAge                       , iMetallicity  , &
         &                                                                      iWavelength
    integer                                                                  :: jWavelength                , imfLookupIndex, &
         &                                                                      jAge                       , jMetallicity
    double precision                                                         :: metallicity
    type            (varying_string              )                           :: message
    character       (len=12                      )                           :: metallicityLabel           , label

    ! Ensure that this IMF is initialized.
    call self%imfInitialize(imfIndex)
    ! Set status code if present.
    if (present(status)) status=errorStatusSuccess
    ! Find the internal lookup index for this IMF.
    imfLookupIndex=self%imfLookup(imfIndex)
    ! Check for out of range conditions.
    if (age > self%spectra(imfLookupIndex)%ages(self%spectra(imfLookupIndex)%agesCount)) then
       write (label,'(e12.4)') age 
       message='age ['//trim(label)//'] exceeds the maximum tabulated ['
       write (label,'(e12.4)') self%spectra(imfLookupIndex)%ages(self%spectra(imfLookupIndex)%agesCount)
       message=message//trim(label)//']'
       call Galacticus_Error_Report('fileLuminosity',message)
    end if
    if (self%forceZeroMetallicity) then
       metallicity=logMetallicityZero
    else
       metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    end if
    if (metallicity > self%spectra(imfLookupIndex)%metallicities(self%spectra(imfLookupIndex)%metallicityCount)+metallicityTolerance) then
       if (present(status)) then
          metallicity=self%spectra(imfLookupIndex)%metallicities(self%spectra(imfLookupIndex)%metallicityCount)
          status     =errorStatusInputDomain
       else
          write (metallicityLabel,'(f12.6)') metallicity
          message='metallicity ['//trim(adjustl(metallicityLabel))//'] exceeds the maximum tabulated ['
          write (metallicityLabel,'(f12.6)') self%spectra(imfLookupIndex)%metallicities(self%spectra(imfLookupIndex)%metallicityCount)
          message=message//trim(adjustl(metallicityLabel))//']'
          call Galacticus_Error_Report('fileLuminosity',message)
       end if
    end if
    ! Assume zero flux outside of the tabulated wavelength range.
    if     (                                                                                                      &
         &   wavelength < self%spectra(imfLookupIndex)%wavelengths(                                            1) &
         & .or.                                                                                                   &
         &   wavelength > self%spectra(imfLookupIndex)%wavelengths(self%spectra(imfLookupIndex)%wavelengthsCount) &
         & ) then
       fileLuminosity=0.0d0
       return
    end if
    ! Get the interpolating factors.
    iAge       =Interpolate_Locate                 (                                                     &
         &                                          self%spectra(imfLookupIndex)%           ages       , &
         &                                          self%spectra(imfLookupIndex)%acceleratorAge        , &
         &                                          age                                                , &
         &                                          self%spectra(imfLookupIndex)%      resetAge          &
         &                                         )
    iWavelength=Interpolate_Locate                 (                                                     &
         &                                          self%spectra(imfLookupIndex)%           wavelengths, &
         &                                          self%spectra(imfLookupIndex)%acceleratorWavelength , &
         &                                                                                  wavelength , &
         &                                          self%spectra(imfLookupIndex)%      resetWavelength   &
         &                                         )
    hAge       =Interpolate_Linear_Generate_Factors(                                                     &
         &                                          self%spectra(imfLookupIndex)%           ages       , &
         &                                                                                 iAge        , &
         &                                                                                  age          &
         &                                         )
    hWavelength=Interpolate_Linear_Generate_Factors(                                                     &
         &                                          self%spectra(imfLookupIndex)%           wavelengths, &
         &                                                                                 iWavelength , &
         &                                                                                  wavelength   &
         &                                         )
    if      (                                                                                                         &
         &    metallicity == logMetallicityZero                                                                       &
         &   .or.                                                                                                     &
         &    metallicity <  self%spectra(imfLookupIndex)%metallicities(                                           1) &
         &  ) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (                                                                                                         &
         &   metallicity >  self%spectra(imfLookupIndex)%metallicities(self%spectra(imfLookupIndex)%metallicityCount) &
         &  ) then
       iMetallicity=self%spectra(imfLookupIndex)%metallicityCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate                 (                                                       &
            &                                           self%spectra(imfLookupIndex)%           metallicities, &
            &                                           self%spectra(imfLookupIndex)%acceleratorMetallicity  , &
            &                                                                                   metallicity  , &
            &                                           self%spectra(imfLookupIndex)%      resetMetallicity    &
            &                                          )
       hMetallicity=Interpolate_Linear_Generate_Factors(                                                       &
            &                                           self%spectra(imfLookupIndex)%           metallicities, &
            &                                                                                  iMetallicity  , &
            &                                                                                   metallicity    &
            &                                          )
    end if
    ! Do the interpolation.
    fileLuminosity=0.0d0
    do jAge=0,1
       do jWavelength=0,1
          do jMetallicity=0,1
             fileLuminosity=+fileLuminosity                                                &
                  &         +self                                                          &
                  &              %spectra(imfLookupIndex)                                  &
                  &                                      %table(                           &
                  &                                             iWavelength +jWavelength , &
                  &                                             iAge        +jAge        , &
                  &                                             iMetallicity+jMetallicity  &
                  &                                            )                           &
                  &         *hAge                              (             jAge        ) &
                  &         *hWavelength                       (             jWavelength ) &
                  &         *hMetallicity                      (             jMetallicity)
           end do
        end do
     end do
     ! Prevent interpolation from returning negative fluxes.
     fileLuminosity=max(fileLuminosity,0.0d0)
    return
  end function fileLuminosity

  subroutine fileImfInitialize(self,imfIndex)
    !% Ensure that data is loaded for the requested IMF.
    use Galacticus_Input_Paths
    implicit none
    class  (stellarPopulationSpectraFile), intent(inout) :: self
    integer                              , intent(in   ) :: imfIndex
    logical                                              :: readFile
    type   (varying_string              )                :: defaultFile  , imfName                         , &
         &                                                  parameterName, stellarPopulationSpectraFileName

    ! Decide if we need to read the file.
    readFile=.not.allocated(self%imfLookup)
    if (.not.readFile) then
       if (size(self%imfLookup) < imfIndex) then
          readFile=.true.
       else
          readFile=(self%imfLookup(imfIndex)==0)
       end if
    end if
    ! Read the file if necessary.
    if (readFile) call self%readFile(imfIndex)
    return
  end subroutine fileImfInitialize

  subroutine fileReadFile(self,imfIndex)
    !% Read a file of simple stellar population spectra.
    use Galacticus_Error
    use Memory_Management
    use IO_HDF5
    use Star_Formation_IMF
    implicit none
    class  (stellarPopulationSpectraFile), intent(inout)             :: self
    integer                              , intent(in   )             :: imfIndex
    integer                              , allocatable, dimension(:) :: imfLookupTemporary
    type   (spectralTable               ), allocatable, dimension(:) :: spectraTemporary
    integer                                                          :: fileFormatVersion , imfLookupIndex
    type   (hdf5Object                  )                            :: spectraFile

    ! Ensure that array for IMF index mappings is sufficiently large.
    if (allocated(self%imfLookup)) then
       if (size(self%imfLookup) < imfIndex) then
          call Move_Alloc(self%imfLookup,imfLookupTemporary)
          call Alloc_Array(self%imfLookup,[imfIndex])
          self%imfLookup(                         1:size(     imfLookupTemporary))=imfLookupTemporary
          self%imfLookup(size(imfLookupTemporary)+1:size(self%imfLookup         ))=0
          call Dealloc_Array(imfLookupTemporary)
       end if
    else
       call Alloc_Array(self%imfLookup,[imfIndex])
       self%imfLookup=0
    end if
    ! If this IMF has not already been read, then assign it a lookup index and expand the spectra array appropriately.
    if (self%imfLookup(imfIndex) == 0) then
       if (allocated(self%spectra)) then
          self%imfLookup(imfIndex)=size(self%spectra)+1
          call Move_Alloc(self%spectra,spectraTemporary)
          allocate(self%spectra(self%imfLookup(imfIndex)))
          self%spectra(1:size(spectraTemporary))=spectraTemporary
          deallocate(spectraTemporary)
          call Memory_Usage_Record(sizeof(self%spectra(1)),blockCount=0)
       else
          self%imfLookup(imfIndex)=1
          allocate(self%spectra(1))
          call Memory_Usage_Record(sizeof(self%spectra))
       end if
       imfLookupIndex=self%imfLookup(imfIndex)
       ! Check that we have a file name for this IMF.
       if (self%fileName(imfIndex) == "?") call Galacticus_Error_Report('fileReadFile',"no file name specified for '"//IMF_Name(imfIndex)//"' IMF")
       !$omp critical(HDF5_Access)
       ! Open the HDF5 file.
       call spectraFile%openFile(char(self%fileName(imfIndex)),readOnly=.true.)
       ! Check that this file has the correct format.
       call spectraFile%readAttribute('fileFormat',fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('fileReadFile','format of stellar tracks file is out of date')
       ! Read the wavelengths array.
       call spectraFile%readDataset('wavelengths'                 ,self%spectra(imfLookupIndex)%wavelengths  )
       self%spectra(imfLookupIndex)%wavelengthsCount=size(self%spectra(imfLookupIndex)%wavelengths  )
       ! Read the ages array.
       call spectraFile%readDataset('ages'                        ,self%spectra(imfLookupIndex)%ages         )
       self%spectra(imfLookupIndex)%agesCount       =size(self%spectra(imfLookupIndex)%ages         )
       ! Read the metallicity array.
       call spectraFile%readDataset('metallicities'               ,self%spectra(imfLookupIndex)%metallicities)
       self%spectra(imfLookupIndex)%metallicityCount=size(self%spectra(imfLookupIndex)%metallicities)
       ! Read the spectra.
       call spectraFile%readDataset('spectra'                     ,self%spectra(imfLookupIndex)%Table        )
       ! Close the HDF5 file.
       call spectraFile%close()
       !$omp end critical(HDF5_Access)
       ! Force interpolation accelerators to be reset.
       self%spectra(imfLookupIndex)%resetAge        =.true.
       self%spectra(imfLookupIndex)%resetWavelength =.true.
       self%spectra(imfLookupIndex)%resetMetallicity=.true.
    end if
    return
  end subroutine fileReadFile

  subroutine fileTabulation(self,imfIndex,agesCount,metallicitiesCount,ages,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified \gls{imf} should be tabulated.
    use Memory_Management
    use Numerical_Constants_Astronomical
    implicit none
    class           (stellarPopulationSpectraFile)                           , intent(inout) :: self
    integer                                                                  , intent(in   ) :: imfIndex
    integer                                                                  , intent(  out) :: agesCount     , metallicitiesCount
    double precision                              , allocatable, dimension(:), intent(  out) :: ages          , metallicity
    integer                                                                                  :: imfLookupIndex

    ! Ensure that this IMF is initialized.
    call self%imfInitialize(imfIndex)
    ! Return the relevant data.
    imfLookupIndex    =self%imfLookup(imfIndex)
    agesCount         =self%spectra(imfLookupIndex)%agesCount
    metallicitiesCount=self%spectra(imfLookupIndex)%metallicityCount
    call Alloc_Array(ages       ,[agesCount         ])
    call Alloc_Array(metallicity,[metallicitiesCount])
    ages              =         self%spectra(imfLookupIndex)%ages
    metallicity       =(10.0d0**self%spectra(imfLookupIndex)%metallicities)*metallicitySolar
    return
  end subroutine fileTabulation

  subroutine fileDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    use Star_Formation_IMF
    implicit none
    class    (stellarPopulationSpectraFile), intent(inout) :: self
    type     (inputParameters             ), intent(inout) :: descriptor
    type     (inputParameters             )                :: subParameters
    character(len=10                      )                :: parameterLabel
    integer                                                :: i
    
    call descriptor%addParameter("stellarPopulationSpectraMethod","file")
    subParameters=descriptor%subparameters("stellarPopulationSpectraMethod")
    if (self%forceZeroMetallicity) then
       parameterLabel="true"
    else
       parameterLabel="false"
    end if
    call subParameters%addParameter("forceZeroMetallicity",trim(adjustl(parameterLabel)))
    do i=1,size(self%fileName)
       call subParameters%addParameter("fileNameFor"//char(IMF_Name(i))//"IMF",char(self%fileName(i)))
    end do
    return
  end subroutine fileDescriptor
