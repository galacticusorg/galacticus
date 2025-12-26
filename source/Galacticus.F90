!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  use    :: Display                   , only : displayMessage                   , displayMagenta           , displayGreen                 , displayReset
  use    :: Display_Verbosity         , only : displayVerbositySetFromParameters
  use    :: Events_Hooks              , only : eventsHooksInitialize
  use    :: Functions_Global_Utilities, only : Functions_Global_Set
  use    :: Display_Banner            , only : Display_Banner_Show
  use    :: Error                     , only : Error_Handler_Register           , Error_Report             , errorStatusSuccess
  use    :: Error_Utilities           , only : Error_Wait_Set_From_Parameters
  use    :: Output_HDF5_Open          , only : Output_HDF5_Close_File           , Output_HDF5_Open_File    , Output_HDF5_Completion_Status
  use    :: ISO_Varying_String        , only : assignment(=)                    , varying_string           , var_str                        , operator(//)
  use    :: Input_Parameters          , only : inputParameter                   , inputParameters
  use    :: Input_Paths               , only : inputPath                        , pathTypeDataStatic
#ifdef USEMPI
  use    :: MPI_F08                   , only : MPI_Comm_World                   , MPI_Thread_Multiple      , MPI_Thread_Single
#endif
#ifdef USEMPI
  use    :: MPI_Utilities             , only : mpiFinalize                      , mpiInitialize
#endif
  !$ use :: OMP_Lib                   , only : OMP_Get_Max_Threads              , OMP_Set_Max_Active_Levels, OMP_Get_Supported_Active_Levels
  use    :: System_Limits             , only : System_Limits_Set
  use    :: Tasks                     , only : task                             , taskClass
  implicit none
  integer                             , parameter                 :: fileNameLengthMaximum =1024
  class    (taskClass                ), pointer                   :: task_
  integer                                                         :: status
  character(len=fileNameLengthMaximum)                            :: parameterFileCharacter    , option
  type     (varying_string           )                            :: parameterFile
  type     (varying_string           ), allocatable, dimension(:) :: changeFiles
  type     (inputParameters          )                            :: parameters
  logical                                                         :: outputFileIsRequired      , dryRun          , &
       &                                                             outputParameters
  integer                                                         :: i                         , countChangeFiles

  ! Initialize MPI.
#ifdef USEMPI
  !$ if (OMP_Get_Max_Threads() > 1) then
  !$  call mpiInitialize(MPI_Thread_Multiple)
  !$ else
      call mpiInitialize(MPI_Thread_Single  )
  !$ end if
#endif
  ! Register error handlers.
  call Error_Handler_Register()
  ! Show the Galacticus banner.
  call Display_Banner_Show()
  ! Check that we have at least one command line argument.
  if (Command_Argument_Count() < 1) call usageError()
  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  parameterFile=parameterFileCharacter
  ! Look for any parameter change files.
  countChangeFiles=0
  if (Command_Argument_Count() > 1) then
     do i=2,Command_Argument_Count()
        call Get_Command_Argument(i,option)
        if (index(option,"--") == 1) exit
        countChangeFiles=countChangeFiles+1
     end do
     allocate(changeFiles(countChangeFiles))
     countChangeFiles=0
     do i=2,Command_Argument_Count()
        call Get_Command_Argument(i,option)
        if (index(option,"--") == 1) exit
        countChangeFiles=countChangeFiles+1
        changeFiles(countChangeFiles)=option
     end do
  end if
  ! Parse any options.
  dryRun          =.false.
  outputParameters=.false.
  if (Command_Argument_Count() > 1) then
     i=1+countChangeFiles
     do while (i < Command_Argument_Count())
        i=i+1
        call Get_Command_Argument(i,option)
        select case (trim(option))
        case ('--dry-run'                    )
           dryRun          =.true.
        case ('--output-processed-parameters')
           outputParameters=.true.
           i=i+1
           if (i > Command_Argument_Count()) call usageError()
           call Get_Command_Argument(i,parameterFileCharacter)
           if (parameterFileCharacter(1:1) == "-") call usageError()
        case default
           call usageError(option=trim(option))
        end select
     end do
  end if
  ! Report on datasets location.
  call displayMessage(displayGreen()//"NOTE: "//displayReset()//"datasets are being read from "//inputPath(pathTypeDataStatic))
  ! Open the parameter file.
  parameters=inputParameters(parameterFile,changeFiles=changeFiles)
  ! Output the processed parameter file.
  if (outputParameters) call parameters%serializeToXML(var_str(trim(parameterFileCharacter)))
  ! Tell OpenMP that nested parallelism is allowed.
  !$ call OMP_Set_Max_Active_Levels(OMP_Get_Supported_Active_Levels())
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Establish global functions.
  call Functions_Global_Set()
  ! Set verbosity.
  call displayVerbositySetFromParameters(parameters)
  ! Set error wait times.
  call Error_Wait_Set_From_Parameters   (parameters)
  ! Set resource limits.
  call System_Limits_Set                (parameters)
  ! Validate parameter file.
  call parameters%checkParameters()
  ! Perform task.
  !![
  <objectBuilder class="task" name="task_" source="parameters"/>
  !!]
  if (.not.dryRun) then
     outputFileIsRequired=task_%requiresOutputFile()
  else
     outputFileIsRequired=.true.
  end if
  if (outputFileIsRequired) call Output_HDF5_Open_File (parameters)
  if (.not.dryRun) then
     call task_%perform(status)
  else
     status=errorStatusSuccess
  end if
  call Output_HDF5_Completion_Status(status)
  if (status /= errorStatusSuccess) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' task failed')
  call parameters%reset  ()
  call parameters%destroy()
  !![
  <objectDestructor name="task_"/>
  !!]
  if (outputFileIsRequired) call Output_HDF5_Close_File()
  ! Finalize MPI.
#ifdef USEMPI
  call MPI_Barrier(MPI_Comm_World,status)
  if (status /= 0) call Error_Report('MPI barrier failed'//{introspection:location})
  call mpiFinalize()
#endif
  ! All done, finish.

contains

  subroutine usageError(option)
    !!{
    Report a usage error.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    character(len=*         ), intent(in   ), optional :: option
    type     (varying_string)                          :: message

    message=""
    if (present(option)) message="unknown option `"//option//"`"//char(10)
    message=message//"Usage: Galacticus.exe <parameterFile> [<changeFile1>...<changeFileN> options...]"                            //char(10)//char(10)
    message=message//"  --dry-run                                  do not perform any task, just parse the parameter file and exit"//char(10)
    message=message//"  --output-processed-parameters <filename>   output the processed parameters to the named file"
    call Error_Report(message)
    return
  end subroutine usageError
  
end program Galacticus
