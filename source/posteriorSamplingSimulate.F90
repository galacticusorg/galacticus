!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a program that determines constraints on \glc\ models.

program Posterior_Sampling_Simulate
  !% Determines constraints on \glc\ models.
  use Memory_Management
  use Input_Parameters
  use ISO_Varying_String
  use Galacticus_Error
  use MPI_Utilities
  use Galacticus_Display
  use Galacticus_Nodes
  use Node_Components
  use Posterior_Sampling_Simulation
  implicit none
  integer                                  , parameter :: fileNameLengthMaximum     =1024
  character(len=fileNameLengthMaximum     )            :: parameterFileCharacter
  type     (varying_string                )            :: parameterFile
  integer                                              :: verbosityLevel
  type     (inputParameters               )            :: parameters
  class    (posteriorSampleSimulationClass), pointer   :: posteriorSampleSimulation_

  ! Initialize MPI.
  call mpiInitialize()
  ! Establish error handlers.
  call Galacticus_Error_Handler_Register() 
  ! Read in basic code memory usage.
  call Code_Memory_Usage('posteriorSamplingSimulate.size')
  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: posteriorSamplingSimulate.exe <parameterFile>")
  parameterFile=parameterFileCharacter
  ! Open the parameter file.
  parameters=inputParameters(parameterFile,allowedParametersFile='posteriorSamplingSimulate.parameters.xml')
  call parameters%markGlobal()
  ! Set the verbosity level.
  !# <inputParameter>
  !#   <name>verbosityLevel</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1</defaultValue>
  !#   <description>The level of verbosity for the {\normalfont \ttfamily Constrain\_Galacticus} code (higher values give more verbosity).</description>
  !#   <source>globalParameters</source>
  !#   <type>integer</type>
  !# </inputParameter>
  call Galacticus_Verbosity_Level_Set(verbosityLevel)
  ! Ensure the nodes objects are initialized.
  call nodeClassHierarchyInitialize()
  call Node_Components_Initialize  ()
  ! Get the simulator.
  posteriorSampleSimulation_ => posteriorSampleSimulation()
  ! Perform the simulation.
  if (mpiSelf%isMaster()) call Galacticus_Display_Indent  ('Begin simulation')
  call posteriorSampleSimulation_%simulate()
  if (mpiSelf%isMaster()) call Galacticus_Display_Unindent('Simulation done' )
  ! Finalize MPI.
  call mpiFinalize()
end program Posterior_Sampling_Simulate
