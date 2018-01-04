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

!% Contains a code which performs analysis on N-body simulation halos.

program nbodyAnalyze
  !% Performs analysis on N-body simulation halos.
  use Galacticus_Error
  use IO_HDF5
  use Memory_Management
  use Input_Parameters
  use NBody_Simulation_Data
  use NBody_Importers
  use NBody_Operators
  implicit none
  integer                             , parameter :: fileNameLengthMaximum=1024
  class    (nBodyImporterClass       ), pointer   :: nBodyImporter_
  class    (nBodyOperatorClass       ), pointer   :: nBodyOperator_
  type     (inputParameters          )            :: parameters
  type     (nBodyData                )            :: simulation
  character(len=fileNameLengthMaximum)            :: parameterFileName         , nbodyFileName  

  ! Read in basic code memory usage.
  call Code_Memory_Usage('nbodyAnalyze.size')
  ! Read arguments.
  if (Command_Argument_Count() /= 2) call Galacticus_Error_Report(message="Usage: nbodyAnalyze.exe <parameterFile> <nbodyFile>")
  call Get_Command_Argument(1,parameterFileName)
  call Get_Command_Argument(2,    nbodyFileName)
  ! Open the parameter file.
  parameters=inputParameters(parameterFileName)
  call parameters%markGlobal()
  ! Load the N-body data.
  nBodyImporter_ => nBodyImporter        (             )
  simulation     =  nBodyImporter_%import(nbodyFileName)
  ! Operate on the N-body data.
  nBodyOperator_ => nBodyOperator()
  call nBodyOperator_%operate(simulation)
  ! Close the analysis group.
  call simulation%analysis%close()
end program nbodyAnalyze
