!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements writing of \glc\ build information to the \glc\ output file.

! Specify an explicit dependence on the gsl_Version.o object file.
!: ./work/build/gsl_Version.o

module Galacticus_Build
  !% Implements writing of \glc\ build information to the \glc\ output file.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Galacticus_Build_Output

  interface
     function GSL_Get_Version() bind(c,name='GSL_Get_Version')
       !% Template for a C function that returns the GSL version string.
       import
       type(c_ptr) :: GSL_Get_Version
     end function GSL_Get_Version
  end interface

contains

  !# <outputFileOpenTask>
  !#  <unitName>Galacticus_Build_Output</unitName>
  !# </outputFileOpenTask>
  subroutine Galacticus_Build_Output
    !% Output build information to the main output file.
    use Galacticus_HDF5
    use Galacticus_Input_Paths
    use IO_HDF5
    use FoX_Common
    use HDF5
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Error
    use FGSL
    use File_Utilities
    implicit none
    character(kind=c_char,len=1), dimension(:), pointer :: charVersionString
    type(varying_string),         dimension(1)          :: changeSet
    integer,                      parameter             :: versionStringlengthMaximum=10
    type(hdf5Object)                                    :: buildGroup
    integer                                             :: hdfVersionMajor,hdfVersionMinor,hdfVersionRelease,hdfError,iChr
    type(c_ptr)                                         :: charVersionPointer
    type(varying_string)                                :: versionString,PREPROCESSOR,FCCOMPILER,CCOMPILER,CPPCOMPILER,MODULETYPE,FCFLAGS &
      &  ,FCFLAGS_NOOPT,CFLAGS,CPPFLAGS,LIBS,FCCOMPILER_VERSION,CCOMPILER_VERSION,CPPCOMPILER_VERSION

    ! Include build environment definitions.
    include 'galacticus.output.build.environment.inc' ! NO_USES

    ! Create a group for build information.
    buildGroup=galacticusOutputFile%openGroup('Build','Build information for this model.')

    ! Write FGSL library version string.
    call buildGroup%writeAttribute(FGSL_Version,'FGSL_library_version')

    ! Write GSL library version string.
    charVersionPointer=GSL_Get_Version()
    call c_f_pointer(charVersionPointer,charVersionString,[versionStringlengthMaximum])
    versionString=""
    iChr=0
    do while (iChr < versionStringlengthMaximum)
      iChr=iChr+1
      if (charVersionString(iChr) == c_null_char) exit
      versionString=versionString//charVersionString(iChr)
    end do
    call buildGroup%writeAttribute(versionString,'GSL_library_version' )

    ! Write FoX library version string.
    call buildGroup%writeAttribute(Fox_Version,'FoX_library_version' )

    ! Write HDF5 library version string.
    call h5get_libversion_f(hdfVersionMajor,hdfVersionMinor,hdfVersionRelease,hdfError)
    if (hdfError /= 0) call Galacticus_Error_Report('Galacticus_Build_Output','unable to get HDF5 library version number')
    versionString=''
    versionString=versionString//hdfVersionMajor//"."
    versionString=versionString//hdfVersionMinor//"."
    versionString=versionString//hdfVersionRelease
    call buildGroup%writeAttribute(versionString,'HDF5_library_version' )

    ! Write Make environment variables.
    call buildGroup%writeAttribute(FCCOMPILER         ,'make_FCCOMPILER         ')
    call buildGroup%writeAttribute(PREPROCESSOR       ,'make_PREPROCESSOR       ')
    call buildGroup%writeAttribute(CCOMPILER          ,'make_CCOMPILER          ')
    call buildGroup%writeAttribute(CPPCOMPILER        ,'make_CPPCOMPILER        ')
    call buildGroup%writeAttribute(MODULETYPE         ,'make_MODULETYPE         ')
    call buildGroup%writeAttribute(FCFLAGS            ,'make_FCFLAGS            ')
    call buildGroup%writeAttribute(FCFLAGS_NOOPT      ,'make_FCFLAGS_NOOPT      ')
    call buildGroup%writeAttribute(CFLAGS             ,'make_CFLAGS             ')
    call buildGroup%writeAttribute(CPPFLAGS           ,'make_CPPFLAGS           ')
    call buildGroup%writeAttribute(FCCOMPILER_VERSION ,'make_FCCOMPILER_VERSION ')
    call buildGroup%writeAttribute(CCOMPILER_VERSION  ,'make_CCOMPILER_VERSION  ')
    call buildGroup%writeAttribute(CPPCOMPILER_VERSION,'make_CPPCOMPILER_VERSION')

    ! Add Bazaar changeset information.
    if (File_Exists(Galacticus_Input_Path()//"work/build/galacticus.hg.patch")) then
       call changeSet(1)%loadFromFile(char(Galacticus_Input_Path()//'work/build/galacticus.hg.patch'))
       if (changeSet(1) /= "" ) call buildGroup%writeDataset(changeSet,'sourceChangeSetDiff','Output of "hg diff" - gives the uncommitted source changeset')
       call changeSet(1)%destroy()
    end if
    if (File_Exists(Galacticus_Input_Path()//"work/build/galacticus.hg.bundle")) then
       call changeSet(1)%loadFromFile(char(Galacticus_Input_Path()//'work/build/galacticus.hg.bundle'))
       if (changeSet(1) /= "" ) call buildGroup%writeDataset(changeSet,'sourceChangeSetBundle','Output of "hg bundle -t none" - gives the committed source changeset')
       call changeSet(1)%destroy()
    end if

    ! Close the build group.
    call buildGroup%close()
    return
  end subroutine Galacticus_Build_Output

end module Galacticus_Build
