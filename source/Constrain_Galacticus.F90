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

!% Contains a program that determines constraints on \glc\ models.

program Galacticus_Constrain
  !% Determines constraints on \glc\ models.
  use Memory_Management
  use Input_Parameters
  use ISO_Varying_String
  use Galacticus_Error
  use MPI_Utilities
  use Constraints_Constrain
  use Galacticus_Display
  implicit none
  integer                             , parameter :: fileNameLengthMaximum =1024
  character(len=fileNameLengthMaximum)            :: parameterFileCharacter     , configFileCharacter
  type     (varying_string           )            :: parameterFile              , configFile
  integer                                         :: verbosityLevel

  ! Initialize MPI.
  call mpiInitialize()
  ! Read in basic code memory usage.
  call Code_Memory_Usage('Constrain_Galacticus.size')
  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Constrain_Galacticus.exe <parameterFile> <simulationConfig>")
  parameterFile=parameterFileCharacter
  call Get_Command_Argument(2,configFileCharacter)
  if (len_trim(   configFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Constrain_Galacticus.exe <parameterFile> <simulationConfig>")
  configFile=configFileCharacter
  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile,allowedParametersFile='Constrain_Galacticus.parameters.xml')
  ! Set the verbosity level.
  !@ <inputParameter>
  !@   <name>verbosityLevel</name>
  !@   <defaultValue>1</defaultValue>
  !@   <attachedTo>module</attachedTo>
  !@   <description>
  !@     The level of verbosity for the {\tt Constrain\_Galacticus} code (higher values give more verbosity).
  !@   </description>
  !@   <type>integer</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('verbosityLevel',verbosityLevel,1)
  call Galacticus_Verbosity_Level_Set(verbosityLevel)
  ! Run the simulation
  call Constrain(configFile)
  ! Close the parameter file.
  call Input_Parameters_File_Close()
  ! Finalize MPI.
  call mpiFinalize()
end program Galacticus_Constrain
