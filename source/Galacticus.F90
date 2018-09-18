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

!% {\normalfont \scshape Galacticus} is a semi-analytic model of galaxy formation written by Andrew Benson
!% \href{mailto:abenson@carnegiescience.edu}{{\normalfont \ttfamily <abenson@carnegiescience.edu>}}.

program Galacticus
  !% The main {\normalfont \scshape Galacticus} program.
  !$ use OMP_Lib
  use Galacticus_Banner
  use Galacticus_Error
  use Galacticus_Output_Open
  use Galacticus_Display_Verbosity
  use Galacticus_HDF5
  use Tasks
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Functions_Global_Utilities
  use Galacticus_Error_Wait
  use System_Limits
#ifdef USEMPI
  use MPI
  use MPI_Utilities
#endif
  implicit none
  integer                             , parameter :: fileNameLengthMaximum =1024
  class    (taskClass                ), pointer   :: task_
  character(len=fileNameLengthMaximum)            :: parameterFileCharacter
  type     (varying_string           )            :: parameterFile
  type     (inputParameters          )            :: parameters
#ifdef USEMPI
  integer                                         :: parentCommunicator         , status
#endif
  
  ! Initialize MPI.
#ifdef USEMPI
  !$ if (OMP_Get_Max_Threads() > 1) then
  !$  call mpiInitialize(MPI_Thread_Multiple)
  !$ else
      call mpiInitialize(MPI_Thread_Single  )
  !$ end if
#endif
  ! Register error handlers.
  call Galacticus_Error_Handler_Register()
  ! Read in basic code memory usage.
  call Code_Memory_Usage('Galacticus.size')
  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Galacticus.exe <parameterFile>")
  parameterFile=parameterFileCharacter
  ! Open the parameter file.
  parameters=inputParameters(parameterFile,allowedParametersFile='Galacticus.parameters.xml')
  call parameters%markGlobal()
  ! Tell OpenMP that nested parallelism is allowed.
  !$ call OMP_Set_Nested(.true.)
  ! Establish global functions.
  call Functions_Global_Set()
  ! Set verbosity.
  call Galacticus_Verbosity_Set_From_Parameters()
  ! Set error wait times.
  call Galacticus_Error_Wait_Set_From_Parameters()
  ! Set resource limits.
  Call System_Limits_Set()
  ! Show the Galacticus banner.
  call Galacticus_Banner_Show()
  ! Perform task.
  task_ => task()
  if (task_%requiresOutputFile()) call Galacticus_Output_Open_File ()
  call task_%perform()
  if (task_%requiresOutputFile()) call Galacticus_Output_Close_File()
  ! Finalize MPI.
#ifdef USEMPI
  call MPI_Comm_Get_Parent(parentCommunicator,status)
  if (parentCommunicator /= MPI_Comm_Null) call MPI_Barrier(parentCommunicator,status)
  call mpiFinalize()
#endif
  ! All done, finish.
end program Galacticus
