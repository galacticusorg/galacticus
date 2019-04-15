!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  subroutine Interface_Cloudy_Initialize(cloudyPath,cloudyVersion)
    !% Initialize the interface with Cloudy, including downloading and compiling Cloudy if necessary.
    use ISO_Varying_String
    use Galacticus_Paths
    use File_Utilities
    use System_Command
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    type   (varying_string), intent(  out) :: cloudyPath, cloudyVersion
    integer                                :: status

    ! Specify Cloudy version.
    cloudyVersion="c17.01"
    ! Specify Cloudy path.
    cloudyPath   =galacticusPath(pathTypeDataDynamic)//cloudyVersion
    ! Download the Cloudy code.
    if (.not.File_Exists(cloudyPath//".tar.gz")) then
       call Galacticus_Display_Message("downloading Cloudy code....",verbosityWorking)
       call System_Command_Do('wget "http://data.nublado.org/cloudy_releases/c17/'//cloudyVersion//'.tar.gz" -O '//cloudyPath//'.tar.gz',status)
       if (status /= 0) call Galacticus_Error_Report("failed to download Cloudy code"//{introspection:location})
    end if
    ! Unpack and patch the code.
    if (.not.File_Exists(cloudyPath)) then
       call Galacticus_Display_Message("unpacking and patching Cloudy code....",verbosityWorking)
       call System_Command_Do("tar -x -v -z -C "//galacticusPath(pathTypeDataDynamic)//" -f "//cloudyPath//".tar.gz",status)
       if (status /= 0 .or. .not.File_Exists(cloudyPath)) call Galacticus_Error_Report("failed to unpack Cloudy code"//{introspection:location})
       call System_Command_Do('sed -i~ -r s/"\\\$res\s+\.=\s+\"native \""/"print \"skip march=native as it breaks the build\\n\""/ '//cloudyPath//'/source/capabilities.pl',status)
       if (status /= 0                                  ) call Galacticus_Error_Report("failed to patch Cloudy code"//{introspection:location})
    end if
    ! Build the code.
    if (.not.File_Exists(cloudyPath//"/source/cloudy.exe")) then
       call Galacticus_Display_Message("compiling Cloudy code....",verbosityWorking)
       call System_Command_Do("cd "//cloudyPath//"/source; chmod u=wrx configure.sh capabilities.pl; make",status)
       if (status /= 0 .or. .not.File_Exists(cloudyPath//"/source/cloudy.exe")) call Galacticus_Error_Report("failed to build Cloudy code"//{introspection:location})
    end if
    ! Append backslash to path before returning.
    cloudyPath=cloudyPath//"/"
    return
  end subroutine Interface_Cloudy_Initialize
  
end module Interfaces_Cloudy
