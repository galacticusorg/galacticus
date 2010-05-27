!% Contains a module which handles stellar spectra using the \cite{conroy_propagation_2009} package.

module Stellar_Population_Spectra_Conroy
  !% Handles stellar spectra using the \cite{conroy_propagation_2009} package.
  use ISO_Varying_String
  private
  public :: Stellar_Population_Spectra_Conroy_Initialize
  
  ! List of which IMFs have been generated and read so far.
  logical, allocatable, dimension(:) :: imfRead

contains
  
  !# <stellarPopulationSpectraMethod>
  !#  <unitName>Stellar_Population_Spectra_Conroy_Initialize</unitName>
  !# </stellarPopulationSpectraMethod>
  subroutine Stellar_Population_Spectra_Conroy_Initialize(stellarPopulationSpectraMethod,Stellar_Population_Spectra_Get&
       &,Stellar_Population_Spectrum_Tabulation_Get)
    !% Initializes the ``Conroy, White \& Gunn'' module.
    implicit none
    type(varying_string),          intent(in)    :: stellarPopulationSpectraMethod
    procedure(),          pointer, intent(inout) :: Stellar_Population_Spectra_Get,Stellar_Population_Spectrum_Tabulation_Get
    
    if (stellarPopulationSpectraMethod == 'Conroy, White & Gunn') then
       Stellar_Population_Spectra_Get             => Stellar_Population_Spectra_Conroy_Get
       Stellar_Population_Spectrum_Tabulation_Get => Stellar_Population_Spectrum_Tabulation_Conroy
    end if
    return
  end subroutine Stellar_Population_Spectra_Conroy_Initialize

  subroutine Stellar_Population_Spectra_Conroy_Initialize_IMF(imfIndex)
    !% Ensure that the requested IMF has been generated and loaded.
    use Abundances_Structure
    use ISO_Varying_String
    use Input_Parameters
    use Star_Formation_IMF
    use File_Utilities
    use Memory_Management
    use System_Command
    use Stellar_Population_Spectra_File
    implicit none
    integer,              intent(in)                :: imfIndex
    logical,              allocatable, dimension(:) :: imfReadTemporary
    double precision,     allocatable, dimension(:) :: imfMass,imfPhi
    integer                                         :: imfUnit,iIMF
    type(varying_string)                            :: stellarPopulationSpectraFile,imfName,command

    ! Ensure that array for IMF index mappings is sufficiently large.
    if (allocated(imfRead)) then
       if (size(imfRead) < imfIndex) then
          call Move_Alloc(imfRead,imfReadTemporary)
          call Alloc_Array(imfRead,imfIndex,'imfRead')
          imfRead(1:size(imfReadTemporary))=imfReadTemporary
          imfRead(size(imfReadTemporary)+1:size(imfRead))=.false.
          call Dealloc_Array(imfReadTemporary)
       end if
    else
       call Alloc_Array(imfRead,imfIndex,'imfRead')
       imfRead=.false.
    end if

    ! If we must read the file, find the file name.
    if (.not.imfRead(imfIndex)) then

       ! Get the name of this IMF.
       imfName=IMF_Name(imfIndex)

       ! Name of the parameter to be used for this IMF.
       stellarPopulationSpectraFile='data/SSP_Spectra_Conroy-et-al_v2.0_imf'//imfName//'.hdf5'

       ! Generate the IMF tabulation.
       if (allocated(imfMass)) call Dealloc_Array(imfMass)
       if (allocated(imfPhi )) call Dealloc_Array(imfPhi )
       call IMF_Tabulate(imfIndex,imfMass,imfPhi)
       imfUnit=File_Units_Get()
       open(unit=imfUnit,file="galacticus.imf",status="unknown",form="formatted")
       do iIMF=1,size(imfMass)
          write (imfUnit,'(2(1x,e12.6))') imfMass(iIMF),imfPhi(iIMF)
       end do
       close(imfUnit)
       call Dealloc_Array(imfMass)
       call Dealloc_Array(imfPhi )

       ! Call the driver script to generate this file.
       command='./scripts/aux/Conroy_SPS_Driver.pl '//imfName//' '//stellarPopulationSpectraFile
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       call Stellar_Population_Spectra_File_Read(imfIndex,stellarPopulationSpectraFile)

       ! Flag that this IMF has been read.
       imfRead(imfIndex)=.true.
       
    end if

    return
  end subroutine Stellar_Population_Spectra_Conroy_Initialize_IMF

  double precision function Stellar_Population_Spectra_Conroy_Get(abundances,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms). This is computed using the
    !% \cite{conroy_propagation_2009} package.
    use Stellar_Population_Spectra_File
    use Abundances_Structure
    implicit none
    type(abundancesStructure), intent(in)                :: abundances
    double precision,          intent(in)                :: age,wavelength
    integer,                   intent(in)                :: imfIndex
 
    ! Ensure that IMF has been initialized. 
    call Stellar_Population_Spectra_Conroy_Initialize_IMF(imfIndex)
    
    ! Call routine to interpolate in the tabulated function.
    Stellar_Population_Spectra_Conroy_Get=Stellar_Population_Spectra_File_Interpolate(abundances,age,wavelength,imfIndex)

    return
  end function Stellar_Population_Spectra_Conroy_Get
    
  subroutine Stellar_Population_Spectrum_Tabulation_Conroy(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    use Memory_Management
    use Stellar_Population_Spectra_File
    implicit none
    integer,          intent(in)                             :: imfIndex
    integer,          intent(out)                            :: agesCount,metallicitiesCount
    double precision, intent(out), allocatable, dimension(:) :: age,metallicity

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_Conroy_Initialize_IMF(imfIndex)

    ! Call the file tabulation routine to get the required data.
    call Stellar_Population_Spectra_File_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)

    return
  end subroutine Stellar_Population_Spectrum_Tabulation_Conroy
  
end module Stellar_Population_Spectra_Conroy
