!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
       stellarPopulationSpectraFile='data/SSP_Spectra_Conroy-et-al_v2.1_imf'//imfName//'.hdf5'

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
