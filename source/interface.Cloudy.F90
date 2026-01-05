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
Contains a module which provides various interfaces to the \gls{cloudy} code.
!!}

module Interfaces_Cloudy
  !!{
  Provides various interfaces to the \gls{cloudy} code.
  !!}
  private
  public :: Interface_Cloudy_Initialize

contains

  subroutine Interface_Cloudy_Initialize(cloudyPath,cloudyVersion,static)
    !!{
    Initialize the interface with Cloudy, including downloading and compiling Cloudy if necessary.
    !!}
    use :: Dependencies      , only : dependencyVersion
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : File_Exists
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//)     , varying_string
    use :: String_Handling   , only : stringSubstitute
    use :: System_Command    , only : System_Command_Do
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageCPlusPlus
    implicit none
    type     (varying_string), intent(  out)           :: cloudyPath   , cloudyVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status       , statusEnvironment, &
         &                                                statusPath
    character(len=   3      )                          :: staticDefault
    character(len=1024      )                          :: compilerPath
    type     (varying_string)                          :: command      , cloudyVersionMajor
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Specify Cloudy version.
    cloudyVersion     ="c"//dependencyVersion("cloudy"                 )
    cloudyVersionMajor=     dependencyVersion("cloudy",majorOnly=.true.)
    ! Specify Cloudy path.
    cloudyPath   =inputPath(pathTypeDataDynamic)//cloudyVersion
    ! Check for existence of executable - build if necessary.
    if (.not.File_Exists(cloudyPath//"/source/cloudy.exe")) then
       ! Check for existence of source code - unpack and patch if necessary.
       if (.not.File_Exists(cloudyPath)) then
          ! Check for existence of tarball - download the Cloudy code if necessary.
          if (.not.File_Exists(cloudyPath//".tar.gz")) then
             call displayMessage("downloading Cloudy code....",verbosityLevelWorking)
             call download('"http://data.nublado.org/cloudy_releases/c'//char(cloudyVersionMajor)//'/'//char(cloudyVersion)//'.tar.gz"',char(cloudyPath)//'.tar.gz',status=status,retries=5,retryWait=60)
             if (status /= 0) then
                ! Try the "old/" subdirectory. Sometimes the source is moved to this path prior to a new release.
                call download('"http://data.nublado.org/cloudy_releases/c'//char(cloudyVersionMajor)//'/old/'//char(cloudyVersion)//'.tar.gz"',char(cloudyPath)//'.tar.gz',status=status,retries=5,retryWait=60)
                if (status /= 0) call Error_Report("failed to download Cloudy code"//{introspection:location})
             end if
          end if
          ! Unpack and patch the code.
          call displayMessage("unpacking and patching Cloudy code....",verbosityLevelWorking)
          call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//" -f "//cloudyPath//".tar.gz",status)
          if (status /= 0 .or. .not.File_Exists(cloudyPath)) call Error_Report("failed to unpack Cloudy code"//{introspection:location})
          call System_Command_Do('sed -i~ -E s/"^#\!\/bin\/sh"/"#\!\/usr\/bin\/env bash"/ '//cloudyPath//'/source/configure.sh',status)
          if (status /= 0                                  ) call Error_Report("failed to patch Cloudy code"//{introspection:location})
          call System_Command_Do('sed -i~ -E s/"^is_repo=(.*)"/"is_repo=\"false\""/ '//cloudyPath//'/source/gitversion.sh',status)
          if (status /= 0                                  ) call Error_Report("failed to patch Cloudy code"//{introspection:location})
          call System_Command_Do('sed -i~ -E s/"\\\$res[[:space:]]+\.=[[:space:]]+\"native \""/"print \"skip march=native as it breaks the build\\n\""/ '//cloudyPath//'/source/capabilities.pl',status)
          if (status /= 0                                  ) call Error_Report("failed to patch Cloudy code"//{introspection:location})
          call System_Command_Do('sed -i~ -E s/"which[[:space:]]+g\+\+"/"which '//compiler(languageCPlusPlus)//'"/ '//cloudyPath//'/source/Makefile',status)
          if (status /= 0                                  ) call Error_Report("failed to patch Cloudy code"//{introspection:location})
       end if
       ! Build the code.
       call displayMessage("compiling Cloudy code....",verbosityLevelWorking)
       if (.not.present(static)) then
          call Get_Environment_Variable('CLOUDY_STATIC_BUILD',staticDefault,status=statusEnvironment)
          if (statusEnvironment <= 0) static_=staticDefault == "yes"
       end if
       if (static_) then
          call System_Command_Do("cd "//cloudyPath//"/source; sed -i~ -E s/'^EXTRA[[:space:]]*=.*'/'EXTRA = -static "//stringSubstitute(stringSubstitute(compilerOptions(languageCPlusPlus),"/","\/"),"-","\-")//"'/g Makefile")
       else
          call System_Command_Do("cd "//cloudyPath//"/source; sed -i~ -E s/'^EXTRA[[:space:]]*=.*'/'EXTRA = "        //stringSubstitute(stringSubstitute(compilerOptions(languageCPlusPlus),"/","\/"),"-","\-")//"'/g Makefile")
       end if
       command="cd "//cloudyPath//"/source; chmod u=wrx configure.sh capabilities.pl;"
       call Get_Environment_Variable('CLOUDY_COMPILER_PATH',compilerPath,status=statusPath)
       if      (statusPath == -1) then
          call Error_Report("can not read Cloudy compiler path environment variable"//{introspection:location})
       else if (statusPath ==  0) then
          command=command//"export PATH="//trim(compilerPath)//":$PATH; "
       end if
       command=command//" make"
       call System_Command_Do(char(command),status)
       if (status /= 0 .or. .not.File_Exists(cloudyPath//"/source/cloudy.exe")) call Error_Report("failed to build Cloudy code"//{introspection:location})
    end if
    ! Append backslash to path before returning.
    cloudyPath=cloudyPath//"/"
    return
  end subroutine Interface_Cloudy_Initialize

end module Interfaces_Cloudy
