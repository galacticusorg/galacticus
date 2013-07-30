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

!% Contains a module which handles stellar spectra using the \cite{conroy_propagation_2009} package.

module Stellar_Population_Spectra_Conroy
  !% Handles stellar spectra using the \cite{conroy_propagation_2009} package.
  use ISO_Varying_String
  implicit none
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
    !% Initializes the ``Conroy-White-Gunn2009'' module.
    implicit none
    type     (varying_string                               ), intent(in   )          :: stellarPopulationSpectraMethod
    procedure(Stellar_Population_Spectra_Conroy_Get        ), intent(inout), pointer :: Stellar_Population_Spectra_Get
    procedure(Stellar_Population_Spectrum_Tabulation_Conroy), intent(inout), pointer :: Stellar_Population_Spectrum_Tabulation_Get

    if (stellarPopulationSpectraMethod == 'Conroy-White-Gunn2009') then
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
    use IO_HDF5
    use Star_Formation_IMF
    use Memory_Management
    use File_Utilities
    use System_Command
    use Stellar_Population_Spectra_File
    use Galacticus_Input_Paths
    use String_Handling
    use Tables
    implicit none
    integer                         , intent(in   )               :: imfIndex
    logical                         , allocatable  , dimension(:) :: imfReadTemporary
    class           (table1D       ), allocatable                 :: imf
    integer                                                       :: iIMF                            , imfUnit, &
         &                                                           fileFormatVersion
    logical                                                       :: remakeFile
    type            (varying_string)                              :: command                         , imfName, &
         &                                                           stellarPopulationSpectraFileName
    type            (hdf5Object    )                              :: stellarPopulationSpectraFile

    ! Ensure that array for IMF index mappings is sufficiently large.
    if (allocated(imfRead)) then
       if (size(imfRead) < imfIndex) then
          call Move_Alloc(imfRead,imfReadTemporary)
          call Alloc_Array(imfRead,[imfIndex])
          imfRead(1:size(imfReadTemporary))=imfReadTemporary
          imfRead(size(imfReadTemporary)+1:size(imfRead))=.false.
          call Dealloc_Array(imfReadTemporary)
       end if
    else
       call Alloc_Array(imfRead,[imfIndex])
       imfRead=.false.
    end if

    ! If we must read the file, find the file name.
    if (.not.imfRead(imfIndex)) then

       ! Get the name of this IMF.
       imfName=IMF_Descriptor(imfIndex)

       ! Name of the parameter to be used for this IMF.
       stellarPopulationSpectraFileName=char(Galacticus_Input_Path())//'data/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.4_imf'//imfName//'.hdf5'

       ! Check if the file exists and has the correct version.
       remakeFile=.false.
       if (File_Exists(char(stellarPopulationSpectraFileName))) then
          !$omp critical(HDF5_Access)
          call stellarPopulationSpectraFile%openFile     (char(stellarPopulationSpectraFileName),readOnly=.true.)
          call stellarPopulationSpectraFile%readAttribute('fileFormat',fileFormatVersion                        )
          if (fileFormatVersion /= Stellar_Population_Spectra_File_Format_Current()) remakeFile=.true.
          !$omp end critical(HDF5_Access)
       else
          remakeFile=.true.
       end if
       ! Generate the file if necessary.
       if (remakeFile) then
          ! Generate the IMF tabulation.
          call IMF_Tabulate(imfIndex,imf)
          open(newunit=imfUnit,file="galacticus.imf",status="unknown",form="formatted")
          do iIMF=1,imf%size()
             write (imfUnit,'(2(1x,e12.6))') imf%x(iIMF),imf%y(iIMF)
          end do
          close(imfUnit)
          call imf%destroy()
          ! Call the driver script to generate this file.
          command=char(Galacticus_Input_Path())//'scripts/aux/Conroy_SPS_Driver.pl '//imfName//' '//stellarPopulationSpectraFileName
          command=command//' '//Stellar_Population_Spectra_File_Format_Current()
          call System_Command_Do(command)
       end if
       ! Call routine to read in the tabulated data.
       call Stellar_Population_Spectra_File_Read(imfIndex,stellarPopulationSpectraFileName)
       ! Flag that this IMF has been read.
       imfRead(imfIndex)=.true.
    end if
    return
  end subroutine Stellar_Population_Spectra_Conroy_Initialize_IMF

  double precision function Stellar_Population_Spectra_Conroy_Get(abundancesStellar,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms). This is computed using the
    !% \cite{conroy_propagation_2009} package.
    use Stellar_Population_Spectra_File
    use Abundances_Structure
    implicit none
    type            (abundances), intent(in   ) :: abundancesStellar
    double precision            , intent(in   ) :: age              , wavelength
    integer                     , intent(in   ) :: imfIndex

    ! Ensure that IMF has been initialized.
    call Stellar_Population_Spectra_Conroy_Initialize_IMF(imfIndex)

    ! Call routine to interpolate in the tabulated function.
    Stellar_Population_Spectra_Conroy_Get=Stellar_Population_Spectra_File_Interpolate(abundancesStellar,age,wavelength,imfIndex)

    return
  end function Stellar_Population_Spectra_Conroy_Get

  subroutine Stellar_Population_Spectrum_Tabulation_Conroy(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    use Memory_Management
    use Stellar_Population_Spectra_File
    implicit none
    integer                                    , intent(in   ) :: imfIndex
    integer                                    , intent(  out) :: agesCount, metallicitiesCount
    double precision, allocatable, dimension(:), intent(  out) :: age      , metallicity

    ! Ensure that this IMF is initialized.
    call Stellar_Population_Spectra_Conroy_Initialize_IMF(imfIndex)

    ! Call the file tabulation routine to get the required data.
    call Stellar_Population_Spectra_File_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)

    return
  end subroutine Stellar_Population_Spectrum_Tabulation_Conroy

end module Stellar_Population_Spectra_Conroy
