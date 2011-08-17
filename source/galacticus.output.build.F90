!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
    use IO_HDF5
    use FoX_Common
    use HDF5
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Error
    use FGSL
    implicit none
    character(kind=c_char,len=1), dimension(:), pointer :: charVersionString
    integer,                      parameter             :: versionStringlengthMaximum=10
    type(hdf5Object)                                    :: buildGroup
    integer                                             :: hdfVersionMajor,hdfVersionMinor,hdfVersionRelease,hdfError,iChr
    type(varying_string)                                :: versionString
    type(c_ptr)                                         :: charVersionPointer
    type(varying_string)                                :: PREPROCESSOR,F03COMPILER,CCOMPILER,CPPCOMPILER,MODULETYPE,F03FLAGS &
      &  ,F03FLAGS_NOOPT,CFLAGS,CPPFLAGS,LIBS,F03COMPILER_VERSION,CCOMPILER_VERSION,CPPCOMPILER_VERSION

    ! Include build environment definitions.
    include 'galacticus.output.build.environment.inc' ! NO_USES

    ! Create a group for build information.
    buildGroup=IO_HDF5_Open_Group(galacticusOutputFile,'Build','Build information for this model.')

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
    call buildGroup%writeAttribute(F03COMPILER        ,'make_F03COMPILER        ')
    call buildGroup%writeAttribute(PREPROCESSOR       ,'make_PREPROCESSOR       ')
    call buildGroup%writeAttribute(CCOMPILER          ,'make_CCOMPILER          ')
    call buildGroup%writeAttribute(CPPCOMPILER        ,'make_CPPCOMPILER        ')
    call buildGroup%writeAttribute(MODULETYPE         ,'make_MODULETYPE         ')
    call buildGroup%writeAttribute(F03FLAGS           ,'make_F03FLAGS           ')
    call buildGroup%writeAttribute(F03FLAGS_NOOPT     ,'make_F03FLAGS_NOOPT     ')
    call buildGroup%writeAttribute(CFLAGS             ,'make_CFLAGS             ')
    call buildGroup%writeAttribute(CPPFLAGS           ,'make_CPPFLAGS           ')
    call buildGroup%writeAttribute(F03COMPILER_VERSION,'make_F03COMPILER_VERSION')
    call buildGroup%writeAttribute(CCOMPILER_VERSION  ,'make_CCOMPILER_VERSION  ')
    call buildGroup%writeAttribute(CPPCOMPILER_VERSION,'make_CPPCOMPILER_VERSION')

    ! Close the version group.
    call buildGroup%close()
    return
  end subroutine Galacticus_Build_Output

end module Galacticus_Build
