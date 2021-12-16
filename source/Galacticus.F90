!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!!{
{\normalfont \scshape Galacticus} is a semi-analytic model of galaxy formation written by Andrew Benson
\href{mailto:abenson@carnegiescience.edu}{{\normalfont \ttfamily <abenson@carnegiescience.edu>}}.
!!}

program Galacticus
  !!{
  The main {\normalfont \scshape Galacticus} program.
  !!}
  use    :: Display_Verbosity         , only : displayVerbositySetFromParameters
  use    :: Events_Hooks              , only : eventsHooksInitialize
  use    :: Functions_Global_Utilities, only : Functions_Global_Set
  use    :: Galacticus_Banner         , only : Galacticus_Banner_Show
  use    :: Galacticus_Error          , only : Galacticus_Error_Handler_Register        , Galacticus_Error_Report
  use    :: Galacticus_Error_Wait     , only : Galacticus_Error_Wait_Set_From_Parameters
  use    :: Galacticus_Output_Open    , only : Galacticus_Output_Close_File             , Galacticus_Output_Open_File
  use    :: ISO_Varying_String        , only : assignment(=)                            , varying_string
  use    :: Input_Parameters          , only : inputParameter                           , inputParameters
#ifdef USEMPI
  use    :: MPI                       , only : MPI_Comm_World                           , MPI_Thread_Multiple        , MPI_Thread_Single
#endif
#ifdef USEMPI
  use    :: MPI_Utilities             , only : mpiFinalize                              , mpiInitialize
#endif
  !$ use :: OMP_Lib                   , only : OMP_Get_Max_Threads                      , OMP_Set_Nested
  use    :: System_Limits             , only : System_Limits_Set
  use    :: Tasks                     , only : task                                     , taskClass
  implicit none
#ifdef USEMPI
  integer                                         :: status
#endif
  integer                             , parameter :: fileNameLengthMaximum =1024
  class    (taskClass                ), pointer   :: task_
  character(len=fileNameLengthMaximum)            :: parameterFileCharacter
  type     (varying_string           )            :: parameterFile
  type     (inputParameters          )            :: parameters

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
  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Galacticus.exe <parameterFile>")
  parameterFile=parameterFileCharacter
  ! Open the parameter file.
  parameters=inputParameters(parameterFile)
  ! Tell OpenMP that nested parallelism is allowed.
  !$ call OMP_Set_Nested(.true.)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Establish global functions.
  call Functions_Global_Set()
  ! Set verbosity.
  call displayVerbositySetFromParameters        (parameters)
  ! Set error wait times.
  call Galacticus_Error_Wait_Set_From_Parameters(parameters)
  ! Set resource limits.
  call System_Limits_Set                        (parameters)
  ! Show the Galacticus banner.
  call Galacticus_Banner_Show()
  ! Validate parameter file.
  call parameters%checkParameters()
  ! Perform task.
  !![
  <objectBuilder class="task" name="task_" source="parameters"/>
  !!]
  if (task_%requiresOutputFile()) call Galacticus_Output_Open_File (parameters)
  call task_     %perform()
  call parameters%destroy()
  if (task_%requiresOutputFile()) call Galacticus_Output_Close_File(          )
  !![
  <objectDestructor name="task_"/>
  !!]
  ! Finalize MPI.
#ifdef USEMPI
  call MPI_Barrier(MPI_Comm_World,status)
  if (status /= 0) call Galacticus_Error_Report('MPI barrier failed'//{introspection:location})
  call mpiFinalize()
#endif
  ! All done, finish.
end program Galacticus
