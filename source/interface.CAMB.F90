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

!% Contains a module which provides various interfaces to the \gls{camb} code.

module Interfaces_CAMB
  !% Provides various interfaces to the \gls{camb} code.
  private
  public :: Interface_CAMB_Initialize
  
contains

  subroutine Interface_CAMB_Initialize(cambPath,cambVersion)
    !% Initialize the interface with CAMB, including downloading and compiling CAMB if necessary.
    use ISO_Varying_String
    use Galacticus_Paths
    use File_Utilities
    use System_Command
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    type   (varying_string), intent(  out) :: cambPath, cambVersion
    integer                                :: status

    ! Set path and version
    cambPath   =galacticusPath(pathTypeDataDynamic)//"CAMB/"
    cambVersion="?"
    ! Download CAMB if necessary.
    if (.not.File_Exists(galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz")) then
       call Galacticus_Display_Message("downloading CAMB code....",verbosityWorking)
       call System_Command_Do("wget http://camb.info/CAMB.tar.gz -O "//galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz",status)
       if (status /= 0 .or. .not.File_Exists(galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz")) call Galacticus_Error_Report("unable to download CAMB"//{introspection:location})
    end if
    ! Unpack the code.
    if (.not.File_Exists(cambPath)) then
       call Galacticus_Display_Message("unpacking CAMB code....",verbosityWorking)
       call System_Command_Do("tar -x -v -z -C "//galacticusPath(pathTypeDataDynamic)//" -f "//galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz");
       if (status /= 0 .or. .not.File_Exists(cambPath)) call Galacticus_Error_Report('failed to unpack CAMB code'//{introspection:location})
    end if
    ! Build the CAMB code.
    if (.not.File_Exists(cambPath//"camb")) then
       call Galacticus_Display_Message("compiling CAMB code",verbosityWorking)
       call System_Command_Do('cd '//cambPath//'; sed -r -i~ s/"ifortErr\s*=.*"/"ifortErr = 1"/ Makefile; sed -r -i~ s/"gfortErr\s*=.*"/"gfortErr = 0"/ Makefile; sed -r -i~ s/"^FFLAGS\s*\+=\s*\-march=native"/"FFLAGS+="/ Makefile; sed -r -i~ s/"^FFLAGS\s*=\s*.*"/"FFLAGS = -Ofast -fopenmp"/ Makefile; find . -name "*.f90" | xargs sed -r -i~ s/"error stop"/"error stop "/; make -j1',status);
       if (status /= 0 .or. .not.File_Exists(cambPath//"camb")) call Galacticus_Error_Report("failed to build CAMB code"//{introspection:location})
    end if
    return
  end subroutine Interface_CAMB_Initialize
  
end module Interfaces_CAMB
