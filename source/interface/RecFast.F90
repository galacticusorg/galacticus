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

!!{RST
Contains a module which provides various interfaces to the :term:`RecFast` code.
!!}

module Interfaces_RecFast
  !!{RST
  Provides various interfaces to the :term:`RecFast` code.
  !!}
  private
  public :: Interface_RecFast_Initialize

contains

  subroutine Interface_RecFast_Initialize(recfastPath,recfastVersion,static)
    !!{RST
    Initialize the interface with RecFast, including downloading and compiling RecFast if necessary.
    !!}
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : Directory_Make   , File_Exists          , File_Lock         , File_Unlock   , &
          &                           lockDescriptor
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeTools        , pathTypeDataStatic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//)      , varying_string, &
         &                            var_str
    use :: System_Command    , only : System_Command_Do, shellEscape
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , languageFortran      , compilerValidate
    implicit none
    type     (varying_string), intent(  out)           :: recfastPath, recfastVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status     , recFastUnit
    character(len=32        )                          :: line       , versionLabel
    type     (varying_string)                          :: command           , pathExe          , &
         &                                                pathPatched       , pathFor          , &
         &                                                pathVersion       , escapedRecfastPath, &
         &                                                escapedPatchSource, escapedPatched
    type     (lockDescriptor)                          :: fileLock
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Set path.
    recfastPath=inputPath(pathTypeTools)//"RecFast/"
    ! Build the code if the executable does not exist.
    pathExe=recfastPath//"recfast.exe"
    if (.not.File_Exists(pathExe)) then
       call compilerValidate(languageFortran,'RecFast')
       call Directory_Make(recfastPath                              )
       call File_Lock     (pathExe    ,fileLock,lockIsShared=.false.)
       ! Patch the code if not already patched.
       pathPatched=recfastPath//"patched"
       if (.not.File_Exists(pathPatched)) then
          ! Download the code if not already downloaded.
          pathFor=recfastPath//"recfast.for"
          if (.not.File_Exists(pathFor)) then
             block
               type(varying_string), dimension(2) :: urls

               call displayMessage("downloading RecFast code....",verbosityLevelWorking)
               urls(1)=var_str("https://www.astro.ubc.ca/people/scott/recfast.for"                                              )
               urls(2)=var_str("https://web.archive.org/web/20250818121544im_/https://www.astro.ubc.ca/people/scott/recfast.for")
               call download(urls,pathFor,retries=5,retryWait=10)
               if (.not.File_Exists(pathFor)) &
                  & call Error_Report("failed to download RecFast code"//{introspection:location})
             end block
          end if
          call displayMessage("patching RecFast code....",verbosityLevelWorking)
          escapedRecfastPath=shellEscape(recfastPath       )
          escapedPatchSource=inputPath  (pathTypeDataStatic)//"patches/RecFast/recfast.for.patch"
          escapedPatchSource=shellEscape(escapedPatchSource)
          command="cp "//escapedPatchSource//" "//escapedRecfastPath//"; cd "//escapedRecfastPath//"; patch < recfast.for.patch"
          call System_Command_Do(command,status)
          if (status /= 0) call Error_Report("failed to patch RecFast file 'recfast.for'"//{introspection:location})
          escapedPatched=shellEscape(recfastPath//"patched")
          command="touch "//escapedPatched
          call System_Command_Do(command)
       end if
       call displayMessage("compiling RecFast code....",verbosityLevelWorking)
       escapedRecfastPath=shellEscape(recfastPath)
       command="cd "//escapedRecfastPath//"; "//compiler(languageFortran)//" recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none"
       if (static_) command=command//" -static"
       call System_Command_Do(char(command))
       if (.not.File_Exists(pathExe)) &
            & call Error_Report("failed to build RecFast code"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    ! Determine the version.
    pathVersion=recfastPath//'currentVersion'
    call File_Lock(pathVersion,fileLock,lockIsShared=.false.)
    if (.not.File_Exists(pathVersion)) then
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
       open(newUnit=recFastUnit,file=char(pathVersion),status='new',form='formatted')
       write (recfastUnit,'(a)') char(recFastVersion)
       close(recFastUnit)
    else
       open(newUnit=recFastUnit,file=char(pathVersion),status='old',form='formatted')
       read (recfastUnit,'(a)') versionLabel
       close(recFastUnit)
       recFastVersion=trim(versionLabel)
    end if
    call File_Unlock(fileLock)
    return
  end subroutine Interface_RecFast_Initialize

end module Interfaces_RecFast
