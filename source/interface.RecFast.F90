!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module which provides various interfaces to the \gls{recfast} code.
!!}

module Interfaces_RecFast
  !!{
  Provides various interfaces to the \gls{recfast} code.
  !!}
  private
  public :: Interface_RecFast_Initialize

contains

  subroutine Interface_RecFast_Initialize(recfastPath,recfastVersion,static)
    !!{
    Initialize the interface with RecFast, including downloading and compiling RecFast if necessary.
    !!}
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : Directory_Make   , File_Exists          , File_Lock         , File_Unlock   , &
          &                           lockDescriptor
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic  , pathTypeDataStatic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//)      , varying_string
    use :: System_Command    , only : System_Command_Do
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , languageFortran
    implicit none
    type     (varying_string), intent(  out)           :: recfastPath, recfastVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status     , recFastUnit
    character(len=32        )                          :: line       , versionLabel
    type     (varying_string)                          :: command
    type     (lockDescriptor)                          :: fileLock
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Set path.
    recfastPath=inputPath(pathTypeDataDynamic)//"RecFast/"
    ! Build the code if the executable does not exist.
    if (.not.File_Exists(recfastPath//"recfast.exe")) then
       call Directory_Make(     recfastPath                                              )
       call File_Lock     (char(recfastPath//"recfast.exe"),fileLock,lockIsShared=.false.)
       ! Patch the code if not already patched.
       if (.not.File_Exists(recfastPath//"patched")) then
          ! Download the code if not already downloaded.
          if (.not.File_Exists(recfastPath//"recfast.for")) then
             call displayMessage("downloading RecFast code....",verbosityLevelWorking)
             call download("https://www.astro.ubc.ca/people/scott/recfast.for",char(recfastPath)//"recfast.for",retries=5,retryWait=60)
             if (.not.File_Exists(recfastPath//"recfast.for")) &
                  & call Error_Report("failed to download RecFast code"//{introspection:location})
          end if
          call displayMessage("patching RecFast code....",verbosityLevelWorking)
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/RecFast/recfast.for.patch "//recfastPath//"; cd "//recfastPath//"; patch < recfast.for.patch",status)
          if (status /= 0) call Error_Report("failed to patch RecFast file 'recfast.for'"//{introspection:location})
          call System_Command_Do("touch "//recfastPath//"patched")
       end if
       call displayMessage("compiling RecFast code....",verbosityLevelWorking)
       command="cd "//recfastPath//"; "//compiler(languageFortran)//" recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none"
       if (static_) command=command//" -static"
       call System_Command_Do(char(command))
       if (.not.File_Exists(recfastPath//"recfast.exe")) &
            & call Error_Report("failed to build RecFast code"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    ! Determine the version.
    call File_Lock(char(inputPath(pathTypeDataDynamic))//'RecFast.currentVersion',fileLock,lockIsShared=.false.)
    if (.not.File_Exists(recfastPath//"currentVersion")) then
       recFastVersion="unknown"
       open(newUnit=recFastUnit,file=char(recfastPath)//"recfast.for",status='old',form='formatted',ioStat=status)
       do while (status == 0)
          read (recFastUnit,'(a)',ioStat=status) line
          if (line(1:2) == "CV" .and. line(4:11) == "Version:") then
             read (line(13:),'(a)') versionLabel
             recFastVersion=trim(versionLabel)
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
       recFastVersion=trim(versionLabel)
    end if
    call File_Unlock(fileLock)
    return
  end subroutine Interface_RecFast_Initialize

end module Interfaces_RecFast
