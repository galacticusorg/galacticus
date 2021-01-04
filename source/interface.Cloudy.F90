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

!% Contains a module which provides various interfaces to the \gls{cloudy} code.

module Interfaces_Cloudy
  !% Provides various interfaces to the \gls{cloudy} code.
  private
  public :: Interface_Cloudy_Initialize

contains

  subroutine Interface_Cloudy_Initialize(cloudyPath,cloudyVersion,static)
    !% Initialize the interface with Cloudy, including downloading and compiling Cloudy if necessary.
    use :: File_Utilities    , only : File_Exists
    use :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: Galacticus_Paths  , only : galacticusPath            , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string            , operator(//)       , assignment(=), char
    use :: System_Command    , only : System_Command_Do
    implicit none
    type     (varying_string), intent(  out)           :: cloudyPath   , cloudyVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status       , statusEnvironment, &
         &                                                statusPath
    character(len=3         )                          :: staticDefault
    character(len=1024      )                          :: compilerPath
    type     (varying_string)                          :: command
    !# <optionalArgument name="static" defaultsTo=".false." />

    ! Specify Cloudy version.
    cloudyVersion="c17.02"
    ! Specify Cloudy path.
    cloudyPath   =galacticusPath(pathTypeDataDynamic)//cloudyVersion
    ! Check for existance of executable - build if necessary.
    if (.not.File_Exists(cloudyPath//"/source/cloudy.exe")) then
       ! Check for existance of source code - unpack and patch if necessary.
       if (.not.File_Exists(cloudyPath)) then
          ! Check for existance of tarball - download the Cloudy code if necessary.
          if (.not.File_Exists(cloudyPath//".tar.gz")) then
             call Galacticus_Display_Message("downloading Cloudy code....",verbosityWorking)
             call System_Command_Do('wget "http://data.nublado.org/cloudy_releases/c17/'//cloudyVersion//'.tar.gz" -O '//cloudyPath//'.tar.gz',status)
             if (status /= 0) call Galacticus_Error_Report("failed to download Cloudy code"//{introspection:location})
          end if
          ! Unpack and patch the code.
          call Galacticus_Display_Message("unpacking and patching Cloudy code....",verbosityWorking)
          call System_Command_Do("tar -x -v -z -C "//galacticusPath(pathTypeDataDynamic)//" -f "//cloudyPath//".tar.gz",status)
          if (status /= 0 .or. .not.File_Exists(cloudyPath)) call Galacticus_Error_Report("failed to unpack Cloudy code"//{introspection:location})
          call System_Command_Do('sed -i~ -r s/"\\\$res\s+\.=\s+\"native \""/"print \"skip march=native as it breaks the build\\n\""/ '//cloudyPath//'/source/capabilities.pl',status)
          if (status /= 0                                  ) call Galacticus_Error_Report("failed to patch Cloudy code"//{introspection:location})
       end if
       ! Build the code.
       call Galacticus_Display_Message("compiling Cloudy code....",verbosityWorking)
       if (.not.present(static)) then
          call Get_Environment_Variable('CLOUDY_STATIC_BUILD',staticDefault,status=statusEnvironment)
          if (statusEnvironment <= 0) static_=staticDefault == "yes"
       end if
       if (static_) then
          call System_Command_Do("cd "//cloudyPath//"/source; sed -i~ -r s/'^EXTRA\s*=.*'/'EXTRA = -static'/g Makefile")
       else
          call System_Command_Do("cd "//cloudyPath//"/source; sed -i~ -r s/'^EXTRA\s*=.*'/'EXTRA ='/g Makefile")
       end if
       command="cd "//cloudyPath//"/source; chmod u=wrx configure.sh capabilities.pl;"
       call Get_Environment_Variable('CLOUDY_COMPILER_PATH',compilerPath,status=statusPath)
       if      (statusPath == -1) then
          call Galacticus_Error_Report("can not read Cloudy compiler path environment variable"//{introspection:location})
       else if (statusPath ==  0) then
          command=command//"export PATH="//trim(compilerPath)//":$PATH; "
       end if
       command=command//" make"
       call System_Command_Do(char(command),status)
       if (status /= 0 .or. .not.File_Exists(cloudyPath//"/source/cloudy.exe")) call Galacticus_Error_Report("failed to build Cloudy code"//{introspection:location})
    end if
    ! Append backslash to path before returning.
    cloudyPath=cloudyPath//"/"
    return
  end subroutine Interface_Cloudy_Initialize

end module Interfaces_Cloudy
