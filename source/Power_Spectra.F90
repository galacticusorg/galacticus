!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a code which computes power spectra and associated data.

program Power_Spectra
  !% Computes power spectra and associated data.
  use Galacticus_Error
  use ISO_Varying_String
  use Input_Parameters
  use Memory_Management
  use Power_Spectrum_Tasks
  implicit none
  integer,                             parameter :: fileNameLengthMaximum=1024
  character(len=fileNameLengthMaximum)           :: fileCharacter
  type(varying_string)                           :: parameterFile,powerSpectrumOutputFileName

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Power_Spectra.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Power_Spectra.exe <parameterFile> <outputFile>")
  parameterFile=fileCharacter
  call Get_Command_Argument(2,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Power_Spectra.exe <parameterFile> <outputFile>")
  powerSpectrumOutputFileName=fileCharacter

  ! Open the output file.
  call Power_Spectrum_Open_File(powerSpectrumOutputFileName)
  
  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile,powerSpectrumOutputFile)

  ! Compute the power spectrum.
  call Power_Spectrum_Compute

  ! Output the power spectrum.
  call Power_Spectrum_Output

  ! Close the parameter file.
  call Input_Parameters_File_Close
  
  ! Close the output file.
  call Power_Spectrum_Close_File

  ! All done, finish.
end program Power_Spectra
