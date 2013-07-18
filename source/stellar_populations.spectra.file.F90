!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which reads and interpolates a file of stellar population spectra.

module Stellar_Population_Spectra_File
  !% Reads and interpolates a file of stellar population spectra.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Stellar_Population_Spectra_File_Initialize, Stellar_Population_Spectra_File_Read,&
       & Stellar_Population_Spectra_File_Interpolate, Stellar_Population_Spectra_File_Tabulation,&
       & Stellar_Population_Spectra_File_Format_Current

  type spectralTable
     !% Structure to hold spectral data.
     ! The spectra tables.
     integer                                                            :: stellarPopulationSpectraAgesNumberPoints              , stellarPopulationSpectraMetallicityNumberPoints       , &
          &                                                                stellarPopulationSpectraWavelengthsNumberPoints
     double precision                   , allocatable, dimension(:)     :: stellarPopulationSpectraAges                          , stellarPopulationSpectraMetallicities                 , &
          &                                                                stellarPopulationSpectraWavelengths
     double precision                   , allocatable, dimension(:,:,:) :: stellarPopulationSpectraTable

     ! Interpolation structures.
     logical                                                            :: resetAge                                       =.true., resetMetallicity                               =.true., &
          &                                                                resetWavelength                                =.true.
     type            (fgsl_interp_accel)                                :: interpolationAcceleratorAge                           , interpolationAcceleratorMetallicity                   , &
          &                                                                interpolationAcceleratorWavelength
  end type spectralTable

  ! Array of spectral tables.
  type   (spectralTable), allocatable, dimension(:) :: spectra

  ! IMF data.
  integer                                           :: imfCount                =0 !   Number of IMFs currently held.
  integer               , allocatable, dimension(:) :: imfLookup                  !   Look-up array to cross-reference IMF indices to our internal data structure.

  ! The current file format version.
  integer               , parameter                 :: fileFormatVersionCurrent=1

contains

  integer function Stellar_Population_Spectra_File_Format_Current()
    !% Return the current file format version for stellar spectra files.
    implicit none

    Stellar_Population_Spectra_File_Format_Current=fileFormatVersionCurrent
    return
  end function Stellar_Population_Spectra_File_Format_Current

  !# <stellarPopulationSpectraMethod>
  !#  <unitName>Stellar_Population_Spectra_File_Initialize</unitName>
  !# </stellarPopulationSpectraMethod>
  subroutine Stellar_Population_Spectra_File_Initialize(stellarPopulationSpectraMethod,Stellar_Population_Spectra_Get&
       &,Stellar_Population_Spectrum_Tabulation_Get)
    !% Initializes the ``stellar population spectra from file'' module.
    implicit none
    type     (varying_string                            ), intent(in   )          :: stellarPopulationSpectraMethod
    procedure(Stellar_Population_Spectra_File_Get       ), intent(inout), pointer :: Stellar_Population_Spectra_Get
    procedure(Stellar_Population_Spectra_File_Tabulation), intent(inout), pointer :: Stellar_Population_Spectrum_Tabulation_Get

    if (stellarPopulationSpectraMethod == 'file') then
       Stellar_Population_Spectra_Get             => Stellar_Population_Spectra_File_Get
       Stellar_Population_Spectrum_Tabulation_Get => Stellar_Population_Spectra_File_Tabulation
    end if
    return
  end subroutine Stellar_Population_Spectra_File_Initialize

  double precision function Stellar_Population_Spectra_File_Get(abundancesStellar,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms). This is found by interpolating in tabulated
    !% spectra.
    use Abundances_Structure
    implicit none
    type            (abundances), intent(in   ) :: abundancesStellar
    double precision            , intent(in   ) :: age              , wavelength
    integer                     , intent(in   ) :: imfIndex

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)

    ! Call routine to interpolate in the tabulated function.
    Stellar_Population_Spectra_File_Get=Stellar_Population_Spectra_File_Interpolate(abundancesStellar,age,wavelength,imfIndex)

    return
  end function Stellar_Population_Spectra_File_Get

  subroutine Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)
    !% Ensure that data is loaded for the requested IMF.
    use Input_Parameters
    use Star_Formation_IMF
    use Galacticus_Input_Paths
    implicit none
    integer                , intent(in   ) :: imfIndex
    logical                                :: readFile
    type   (varying_string)                :: defaultFile  , imfName                     , &
         &                                    parameterName, stellarPopulationSpectraFile

    ! Decide if we need to read the file.
    readFile=.not.allocated(imfLookup) ! If out lookup array is not allocated, we must read the file.
    if (.not.readFile) then
       if (size(imfLookup) < imfIndex) then
          readFile=.true.
       else
          readFile=(imfLookup(imfIndex)==0)
       end if
    end if

    ! If we must read the file, find the file name.
    if (readFile) then

       ! Get the name of this IMF.
       imfName=IMF_Name(imfIndex)

       ! Name of the parameter to be used for this IMF.
       parameterName='stellarPopulationSpectraFor'//imfName//'IMF'

       ! Default file name for this IMF.
       defaultFile=char(Galacticus_Input_Path())//'data/SSP_Spectra_imf'//imfName//'.hdf5'

       ! Get the file name.
       call Get_Input_Parameter(char(parameterName),stellarPopulationSpectraFile,defaultValue=char(defaultFile))

       ! Call routine to read in the tabulated data.
       call Stellar_Population_Spectra_File_Read(imfIndex,stellarPopulationSpectraFile)

    end if
    return
  end subroutine Stellar_Population_Spectra_File_Initialize_IMF

  double precision function Stellar_Population_Spectra_File_Interpolate(abundancesStellar,age,wavelength,imfIndex)
    !% Compute the stellar spectrum by interpolation in the tabulated data.
    use Abundances_Structure
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    type            (abundances    ), intent(in   )  :: abundancesStellar
    double precision                , intent(in   )  :: age                        , wavelength
    integer                         , intent(in   )  :: imfIndex
    double precision                , dimension(0:1) :: hAge                       , hMetallicity  , &
         &                                              hWavelength
    double precision                , parameter      :: metallicityTolerance=0.01d0
    integer                                          :: iAge                       , iMetallicity  , &
         &                                              iWavelength                , imfLookupIndex, &
         &                                              jAge                       , jMetallicity  , &
         &                                              jWavelength
    double precision                                 :: metallicity
    type            (varying_string)                 :: message
    character       (len=12        )                 :: metallicityLabel

    ! Find the internal lookup index for this IMF.
    imfLookupIndex=imfLookup(imfIndex)

    ! Check for out of range conditions.
    if (age > spectra(imfLookupIndex)%stellarPopulationSpectraAges(spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints))&
         & call Galacticus_Error_Report('Stellar_Population_Spectra_File_Interpolate','age exceeds the maximum tabulated')
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=logarithmicByMassSolar)
    if (metallicity > spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints)+metallicityTolerance) then
       write (metallicityLabel,'(f12.6)') metallicity
       message='metallicity ['//trim(adjustl(metallicityLabel))//'] exceeds the maximum tabulated ['
       write (metallicityLabel,'(f12.6)') spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints)
       message=message//trim(adjustl(metallicityLabel))//']'
       call Galacticus_Error_Report('Stellar_Population_Spectra_File_Interpolate',message)
    end if
    
    ! Assume zero flux outside of the tabulated wavelength range.
    if     (                                                                                                                                                   &
         &   wavelength < spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths(                                                                      1) &
         & .or.                                                                                                                                                &
         &   wavelength > spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints) &
         & ) then
       Stellar_Population_Spectra_File_Interpolate=0.0d0
       return
    end if
    
    ! Get the interpolations.
    iAge=Interpolate_Locate(spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints&
         &,spectra(imfLookupIndex)%stellarPopulationSpectraAges,spectra(imfLookupIndex)%interpolationAcceleratorAge,age&
         &,spectra(imfLookupIndex)%resetAge)
    hAge=Interpolate_Linear_Generate_Factors(spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints &
         &,spectra(imfLookupIndex)%stellarPopulationSpectraAges,iAge,age)
    iWavelength=Interpolate_Locate(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints &
         &,spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths &
         &,spectra(imfLookupIndex)%interpolationAcceleratorWavelength,wavelength,spectra(imfLookupIndex)%resetWavelength)
    hWavelength=Interpolate_Linear_Generate_Factors(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints &
         &,spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths,iWavelength,wavelength)

    if (metallicity == logMetallicityZero .or. metallicity < spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints)) then
       iMetallicity=spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints &
            &,spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities &
            &,spectra(imfLookupIndex)%interpolationAcceleratorMetallicity,metallicity,spectra(imfLookupIndex)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints &
            &,spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities,iMetallicity,metallicity)
    end if

     ! Do the interpolation.
     Stellar_Population_Spectra_File_Interpolate=0.0d0
     do jAge=0,1
        do jWavelength=0,1
           do jMetallicity=0,1
              Stellar_Population_Spectra_File_Interpolate=                                                            &
                   &                                       Stellar_Population_Spectra_File_Interpolate                &
                   &                                      +spectra(imfLookupIndex)                                    &
                   &                                        %stellarPopulationSpectraTable(                           &
                   &                                                                       iWavelength +jWavelength , &
                   &                                                                       iAge        +jAge        , &
                   &                                                                       iMetallicity+jMetallicity  &
                   &                                                                      )                           &
                   &                                      *hAge                           (             jAge        ) &
                   &                                      *hWavelength                    (             jWavelength ) &
                   &                                      *hMetallicity                   (             jMetallicity)
           end do
        end do
     end do
     ! Prevent interpolation from returning negative fluxes.
     Stellar_Population_Spectra_File_Interpolate=max(Stellar_Population_Spectra_File_Interpolate,0.0d0)

     return
  end function Stellar_Population_Spectra_File_Interpolate

  subroutine Stellar_Population_Spectra_File_Read(imfIndex,stellarPopulationSpectraFileToRead)
    !% Read a file of simple stellar population spectra.
    use Galacticus_Error
    use Memory_Management
    use IO_HDF5
    implicit none
    integer                                           , intent(in   ) :: imfIndex
    type   (varying_string)                           , intent(in   ) :: stellarPopulationSpectraFileToRead
    integer                , allocatable, dimension(:)                :: imfLookupTemporary
    type   (spectralTable ), allocatable, dimension(:)                :: spectraTemporary
    integer                                                           :: fileFormatVersion                 , imfLookupIndex
    type   (hdf5Object    )                                           :: stellarPopulationSpectraFile

    ! Ensure that array for IMF index mappings is sufficiently large.
    if (allocated(imfLookup)) then
       if (size(imfLookup) < imfIndex) then
          call Move_Alloc(imfLookup,imfLookupTemporary)
          call Alloc_Array(imfLookup,[imfIndex])
          imfLookup(1:size(imfLookupTemporary))=imfLookupTemporary
          imfLookup(size(imfLookupTemporary)+1:size(imfLookup))=0
          call Dealloc_Array(imfLookupTemporary)
       end if
    else
       call Alloc_Array(imfLookup,[imfIndex])
       imfLookup=0
    end if

    ! If this IMF has not already been read, then assign it a lookup index and expand the spectra array appropriately.
    if (imfLookup(imfIndex) == 0) then
       if (allocated(spectra)) then
          imfLookup(imfIndex)=size(spectra)+1
          call Move_Alloc(spectra,spectraTemporary)
          allocate(spectra(imfLookup(imfIndex)))
          spectra(1:size(spectraTemporary))=spectraTemporary
          deallocate(spectraTemporary)
          call Memory_Usage_Record(sizeof(spectra(1)),blockCount=0)
       else
          imfLookup(imfIndex)=1
          allocate(spectra(1))
          call Memory_Usage_Record(sizeof(spectra))
       end if
       imfLookupIndex=imfLookup(imfIndex)

       !$omp critical(HDF5_Access)
       ! Open the HDF5 file.
       call stellarPopulationSpectraFile%openFile(char(stellarPopulationSpectraFileToRead),readOnly=.true.)

       ! Check that this file has the correct format.
       call stellarPopulationSpectraFile%readAttribute('fileFormat',fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Stellar_Population_Spectra_File_Read','format of stellar tracks file is out of date')

       ! Read the wavelengths array.
       call stellarPopulationSpectraFile%readDataset('wavelengths'                 ,spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths  )
       spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints=size(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths  )

       ! Read the ages array.
       call stellarPopulationSpectraFile%readDataset('ages'                        ,spectra(imfLookupIndex)%stellarPopulationSpectraAges         )
       spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints       =size(spectra(imfLookupIndex)%stellarPopulationSpectraAges         )

       ! Read the metallicity array.
       call stellarPopulationSpectraFile%readDataset('metallicities'               ,spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities)
       spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints=size(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities)

       ! Read the spectra.
       call stellarPopulationSpectraFile%readDataset('spectra'                     ,spectra(imfLookupIndex)%stellarPopulationSpectraTable        )

       ! Close the HDF5 file.
       call stellarPopulationSpectraFile%close()
       !$omp end critical(HDF5_Access)

       ! Force interpolation accelerators to be reset.
       spectra(imfLookupIndex)%resetAge        =.true.
       spectra(imfLookupIndex)%resetWavelength =.true.
       spectra(imfLookupIndex)%resetMetallicity=.true.

    end if

    return
  end subroutine Stellar_Population_Spectra_File_Read

  subroutine Stellar_Population_Spectra_File_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    use Memory_Management
    use Numerical_Constants_Astronomical
    implicit none
    integer                                    , intent(in   ) :: imfIndex
    integer                                    , intent(  out) :: agesCount     , metallicitiesCount
    double precision, allocatable, dimension(:), intent(  out) :: age           , metallicity
    integer                                                    :: imfLookupIndex

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)

    ! Return the relevant data.
    imfLookupIndex=imfLookup(imfIndex)
    agesCount         =spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints
    metallicitiesCount=spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints
    call Alloc_Array(age        ,[agesCount         ])
    call Alloc_Array(metallicity,[metallicitiesCount])
    age               =         spectra(imfLookupIndex)%stellarPopulationSpectraAges
    metallicity       =(10.0d0**spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities)*metallicitySolar

    return
  end subroutine Stellar_Population_Spectra_File_Tabulation

end module Stellar_Population_Spectra_File
