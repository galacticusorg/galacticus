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

!% {\sc Galacticus} is a semi-analytic model of galaxy formation written by Andrew Benson \href{mailto:abenson@caltech.edu}{{\tt
!% <abenson@caltech.edu>}}.

program Galacticus
  !% The main {\sc Galacticus} program.
  use Galacticus_Banner
  use Galacticus_Error
  use Galacticus_Tasks
  use Galacticus_HDF5
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  implicit none
  integer,                             parameter :: fileNameLengthMaximum=1024
  character(len=fileNameLengthMaximum)           :: parameterFileCharacter
  type(varying_string)                           :: parameterFile

  ! Show the Galacticus banner.
  call Galacticus_Banner_Show

  ! Register error handlers.
  call Galacticus_Error_Handler_Register

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Galacticus.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Galacticus.exe <parameterFile>")
  parameterFile=parameterFileCharacter

  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile,galacticusOutputFile,allowedParametersFile='Galacticus.parameters.xml')

  ! Perform tasks, until all tasks are done.
  call Galacticus_Task_Do

  ! Close the parameter file.
  call Input_Parameters_File_Close
  
  ! All done, finish.
end program Galacticus
