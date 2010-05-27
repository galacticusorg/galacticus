!% Contains a module which reads and interpolates a file of stellar population spectra.

module Stellar_Population_Spectra_File
  !% Reads and interpolates a file of stellar population spectra.
  use ISO_Varying_String
  use FGSL
  use Stellar_Population_Spectra_Table
  private
  public :: Stellar_Population_Spectra_File_Initialize, Stellar_Population_Spectra_File_Read,&
       & Stellar_Population_Spectra_File_Interpolate, Stellar_Population_Spectra_File_Tabulation
  
  ! Array of spectral tables.
  type(spectralTable), allocatable, dimension(:) :: spectra

  ! IMF data.
  integer                            :: imfCount=0 ! Number of IMFs currently held.
  integer, allocatable, dimension(:) :: imfLookup ! Look-up array to cross-reference IMF indices to our internal data structure.

contains
  
  !# <stellarPopulationSpectraMethod>
  !#  <unitName>Stellar_Population_Spectra_File_Initialize</unitName>
  !# </stellarPopulationSpectraMethod>
  subroutine Stellar_Population_Spectra_File_Initialize(stellarPopulationSpectraMethod,Stellar_Population_Spectra_Get&
       &,Stellar_Population_Spectrum_Tabulation_Get)
    !% Initializes the ``stellar population spectra from file'' module.
    implicit none
    type(varying_string),          intent(in)    :: stellarPopulationSpectraMethod
    procedure(),          pointer, intent(inout) :: Stellar_Population_Spectra_Get,Stellar_Population_Spectrum_Tabulation_Get
    
    if (stellarPopulationSpectraMethod == 'file') then
       Stellar_Population_Spectra_Get             => Stellar_Population_Spectra_File_Get
       Stellar_Population_Spectrum_Tabulation_Get => Stellar_Population_Spectra_File_Tabulation
    end if
    return
  end subroutine Stellar_Population_Spectra_File_Initialize

  double precision function Stellar_Population_Spectra_File_Get(abundances,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms). This is found by interpolating in tabulated
    !% spectra.
    use Abundances_Structure
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    double precision,          intent(in) :: age,wavelength
    integer,                   intent(in) :: imfIndex

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)

    ! Call routine to interpolate in the tabulated function.
    Stellar_Population_Spectra_File_Get=Stellar_Population_Spectra_File_Interpolate(abundances,age,wavelength,imfIndex)

    return
  end function Stellar_Population_Spectra_File_Get

  subroutine Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)
    !% Ensure that data is loaded for the requested IMF.
    use Input_Parameters
    use Star_Formation_IMF
    use ISO_Varying_String
    implicit none
    integer, intent(in)   :: imfIndex
    logical               :: readFile
    type(varying_string)  :: parameterName,defaultFile,stellarPopulationSpectraFile,imfName

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
       defaultFile='data/SSP_Spectra_imf'//imfName//'.hdf5'

       ! Get the file name.
       call Get_Input_Parameter(char(parameterName),stellarPopulationSpectraFile,defaultValue=char(defaultFile))

       ! Call routine to read in the tabulated data.
       call Stellar_Population_Spectra_File_Read(imfIndex,stellarPopulationSpectraFile)
       
    end if
    return
  end subroutine Stellar_Population_Spectra_File_Initialize_IMF
    
  double precision function Stellar_Population_Spectra_File_Interpolate(abundances,age,wavelength,imfIndex)
    !% Compute the stellar spectrum by interpolation in the tabulated data.
    use Abundances_Structure
    use Numerical_Constants_Astronomical
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    type(abundancesStructure), intent(in)     :: abundances
    double precision,          intent(in)     :: age,wavelength
    integer,                   intent(in)     :: imfIndex
    integer                                   :: imfLookupIndex,iAge,iWavelength,iMetallicity,jAge,jWavelength,jMetallicity
    double precision                          :: metallicity
    double precision,          dimension(0:1) :: hAge,hMetallicity,hWavelength
 
    ! Find the internal lookup index for this IMF.
    imfLookupIndex=imfLookup(imfIndex)

    ! Check for out of range conditions.
    if (age > spectra(imfLookupIndex)%stellarPopulationSpectraAges(spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints))&
         & call Galacticus_Error_Report('Stellar_Population_Spectra_File_Interpolate','age exceeds the maximum tabulated')
    if (wavelength < spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths(1) .or.&
         & wavelength > spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints))&
         & call Galacticus_Error_Report('Stellar_Population_Spectra_File_Interpolate','wavelength is out of range')
    metallicity=Abundances_Get_Metallicity(abundances,metallicityType=logarithmicByMassSolar)
    if (metallicity > spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints))&
         & call Galacticus_Error_Report('Stellar_Population_Spectra_File_Interpolate','metallicity exceeds the maximum tabulated')

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
              Stellar_Population_Spectra_File_Interpolate=Stellar_Population_Spectra_File_Interpolate&
                   &+spectra(imfLookupIndex)%stellarPopulationSpectraTable(iMetallicity+jMetallicity,iAge+jAge,iWavelength&
                   &+jWavelength)*hAge(jAge)*hWavelength(jWavelength)*hMetallicity(jMetallicity)
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
    use HDF5
    use H5Lt
    implicit none
    integer,              intent(in)                :: imfIndex
    type(varying_string), intent(in)                :: stellarPopulationSpectraFileToRead
    integer,              allocatable, dimension(:) :: imfLookupTemporary
    type(spectralTable),  allocatable, dimension(:) :: spectraTemporary
    integer                                         :: errorCode,typeClass,imfLookupIndex
    integer(HSIZE_T)                                :: dimensions(3)
    integer(HID_T)                                  :: fileIndex
    integer(SIZE_T)                                 :: typeSize

    ! Ensure that array for IMF index mappings is sufficiently large.
    if (allocated(imfLookup)) then
       if (size(imfLookup) < imfIndex) then
          call Move_Alloc(imfLookup,imfLookupTemporary)
          call Alloc_Array(imfLookup,imfIndex,'imfLookup')
          imfLookup(1:size(imfLookupTemporary))=imfLookupTemporary
          imfLookup(size(imfLookupTemporary)+1:size(imfLookup))=0
          call Dealloc_Array(imfLookupTemporary)
       end if
    else
       call Alloc_Array(imfLookup,imfIndex,'imfLookup')
       imfLookup=0
    end if

    ! If this IMF has not already been read, then assign it a lookup index and expand the spectra array appropriately.
    if (imfLookup(imfIndex) == 0) then
       if (allocated(spectra)) then
          imfLookup(imfIndex)=size(spectra)+1
          call Move_Alloc(spectra,spectraTemporary)
          call Alloc_Array(spectra,imfLookup(imfIndex),'spectra')
          spectra(1:size(spectraTemporary))=spectraTemporary
          call Dealloc_Array(spectraTemporary)
       else
          imfLookup(imfIndex)=1
          allocate(spectra(1))
       end if
       imfLookupIndex=imfLookup(imfIndex)
     
       ! Initialize the HDF5 system.
       call IO_HDF5_Initialize
       
       ! Open the HDF5 file.
       call h5fopen_f(char(stellarPopulationSpectraFileToRead),H5F_ACC_RDONLY_F,fileIndex,errorCode)
       
       ! Read the wavelengths array.
       call h5ltget_dataset_info_f(fileIndex,'wavelengths',dimensions,typeClass,typeSize,errorCode)
       spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints=dimensions(1)
       call Alloc_Array(spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths&
            &,spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints &
            &,'spectra()%stellarPopulationSpectraWavelengths')
       call h5ltread_dataset_double_f(fileIndex,'wavelengths',spectra(imfLookupIndex)%stellarPopulationSpectraWavelengths&
            &,dimensions,errorCode)
       
       ! Read the ages array.
       call h5ltget_dataset_info_f(fileIndex,'ages',dimensions,typeClass,typeSize,errorCode)
       spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints=dimensions(1)
       call Alloc_Array(spectra(imfLookupIndex)%stellarPopulationSpectraAges&
            &,spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints ,'spectra()%stellarPopulationSpectraAges')
       call h5ltread_dataset_double_f(fileIndex,'ages',spectra(imfLookupIndex)%stellarPopulationSpectraAges,dimensions,errorCode)
       
       ! Read the metallicity array.
       call h5ltget_dataset_info_f(fileIndex,'metallicities',dimensions,typeClass,typeSize,errorCode)
       spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints=dimensions(1)
       call Alloc_Array(spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities&
            &,spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints &
            &,'spectra()%stellarPopulationSpectraMetallicities')
       call h5ltread_dataset_double_f(fileIndex,'metallicities',spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities&
            &,dimensions ,errorCode)
       
       ! Allocate space for the spectra.
       call Alloc_Array(spectra(imfLookupIndex)%stellarPopulationSpectraTable&
            &,spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints &
            &,spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints&
            &,spectra(imfLookupIndex)%stellarPopulationSpectraWavelengthsNumberPoints ,'spectra%stellarPopulationSpectraTable')
       
       ! Read the spectra.
       call h5ltread_dataset_double_f(fileIndex,'spectra',spectra(imfLookupIndex)%stellarPopulationSpectraTable,dimensions&
            &,errorCode)
       
       ! Close the HDF5 file. 
       call h5fclose_f(fileIndex,errorCode)       

       ! Force interpolation accelerators to be reset.
       spectra(imfLookupIndex)%resetAge        =.true.
       spectra(imfLookupIndex)%resetWavelength =.true.
       spectra(imfLookupIndex)%resetMetallicity=.true.
       
       ! Uninitialize the HDF5 system.
       call IO_HDF5_Uninitialize
       
    end if

    return
  end subroutine Stellar_Population_Spectra_File_Read

  subroutine Stellar_Population_Spectra_File_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    use Memory_Management
    use Numerical_Constants_Astronomical
    implicit none
    integer,          intent(in)                             :: imfIndex
    integer,          intent(out)                            :: agesCount,metallicitiesCount
    double precision, intent(out), allocatable, dimension(:) :: age,metallicity
    integer                                                  :: imfLookupIndex

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_File_Initialize_IMF(imfIndex)

    ! Return the relevant data.
    imfLookupIndex=imfLookup(imfIndex)
    agesCount         =spectra(imfLookupIndex)%stellarPopulationSpectraAgesNumberPoints
    metallicitiesCount=spectra(imfLookupIndex)%stellarPopulationSpectraMetallicityNumberPoints
    call Alloc_Array(age        ,agesCount         ,'age'        )
    call Alloc_Array(metallicity,metallicitiesCount,'metallicity')
    age               =         spectra(imfLookupIndex)%stellarPopulationSpectraAges
    metallicity       =(10.0d0**spectra(imfLookupIndex)%stellarPopulationSpectraMetallicities)*metallicitySolar

    return
  end subroutine Stellar_Population_Spectra_File_Tabulation
  
end module Stellar_Population_Spectra_File
