!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which provides various interfaces to the \gls{recfast} code.

module Interfaces_RecFast
  !% Provides various interfaces to the \gls{recfast} code.
  private
  public :: Interface_RecFast_Initialize

contains

  subroutine Interface_RecFast_Initialize(recfastPath,recfastVersion,static)
    !% Initialize the interface with RecFast, including downloading and compiling RecFast if necessary.
    use :: File_Utilities    , only : Directory_Make            , File_Exists        , File_Lock   , File_Unlock , &
         &                            lockDescriptor
    use :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: Galacticus_Paths  , only : galacticusPath            , pathTypeDataDynamic, pathTypeExec
    use :: ISO_Varying_String, only : varying_string            , assignment(=)      , char        , operator(//)
    use :: System_Command    , only : System_Command_Do
    implicit none
    type     (varying_string), intent(  out)           :: recfastPath, recfastVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status     , recFastUnit
    character(len=32        )                          :: line       , versionLabel
    type     (varying_string)                          :: command
    type     (lockDescriptor)                          :: fileLock
    !# <optionalArgument name="static" defaultsTo=".false." />

    ! Set path.
    recfastPath=galacticusPath(pathTypeDataDynamic)//"RecFast/"
    ! Build the code if the executable does not exist.
    if (.not.File_Exists(recfastPath//"recfast.exe")) then
       call Directory_Make(recfastPath)
       call File_Lock(char(recfastPath//"recfast.exe"),fileLock,lockIsShared=.false.)
       ! Patch the code if not already patched.
       if (.not.File_Exists(recfastPath//"patched")) then
          ! Download the code if not already downloaded.
          if (.not.File_Exists(recfastPath//"recfast.for")) then
             call Galacticus_Display_Message("downloading RecFast code....",verbosityWorking)
             call System_Command_Do("wget --no-check-certificate https://www.astro.ubc.ca/people/scott/recfast.for -O "//recfastPath//"recfast.for")
             if (.not.File_Exists(recfastPath//"recfast.for")) &
                  & call Galacticus_Error_Report("failed to download RecFast code"//{introspection:location})
          end if
          call Galacticus_Display_Message("patching RecFast code....",verbosityWorking)
          call System_Command_Do("cp "//galacticusPath(pathTypeExec)//"aux/RecFast_Galacticus_Modifications/recfast.for.patch "//recfastPath//"; cd "//recfastPath//"; patch < recfast.for.patch",status)
          if (status /= 0) call Galacticus_Error_Report("failed to patch RecFast file 'recfast.for'"//{introspection:location})
          call System_Command_Do("touch "//recfastPath//"patched")
       end if
       call Galacticus_Display_Message("compiling RecFast code....",verbosityWorking)
       command="cd "//recfastPath//"; gfortran recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none"
       if (static_) command=command//" -static"
       call System_Command_Do(char(command))
       if (.not.File_Exists(recfastPath//"recfast.exe")) &
            & call Galacticus_Error_Report("failed to build RecFast code"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    ! Determine the version.
    call File_Lock(char(recfastPath//"currentVersion"),fileLock,lockIsShared=.false.)
    if (.not.File_Exists(recfastPath//"currentVersion")) then
       recFastVersion="unknown"
       open(newUnit=recFastUnit,file=char(recfastPath)//"recfast.for",status='old',form='formatted',ioStat=status)
       do while (status == 0)
          read (recFastUnit,'(a)',ioStat=status) line
          if (line(1:2) == "CV" .and. line(4:11) == "Version:") then
             read (line(13:),'(a)') versionLabel
             recFastVersion=versionLabel
          end if
       end do
       close(recFastUnit)
       open(newUnit=recFastUnit,file=char(recfastPath)//"currentVersion",status='new',form='formatted')
       write (recfastUnit,'(a)') char(recFastVersion)
       close(recFastUnit)
    else
       open(newUnit=recFastUnit,file=char(recfastPath)//"currentVersion",status='old',form='formatted')
       read (recfastUnit,'(a)') versionLabel
       close(recFastUnit)
       recFastVersion=versionLabel
    end if
    call File_Unlock(fileLock)
    return
  end subroutine Interface_RecFast_Initialize

end module Interfaces_RecFast
