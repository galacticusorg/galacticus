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
Contains a module which implements writing of \glc\ build information to the \glc\ output file.
!!}

! Specify an explicit dependence on the gsl_Version.o object file.
!: $(BUILDPATH)/gsl_Version.o

module Output_Build
  !!{
  Implements writing of \glc\ build information to the \glc\ output file.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_ptr
  implicit none
  private
  public :: Output_Build_String

  interface
     function GSL_Get_Version() bind(c,name='GSL_Get_Version')
       !!{
       Template for a C function that returns the GSL version string.
       !!}
       import
       type(c_ptr) :: GSL_Get_Version
     end function GSL_Get_Version
  end interface

contains

  function Output_Build_String()
    !!{
    Returns a string describing the build environment of \glc.
    !!}
    use            :: FoX_Common        , only : Fox_Version
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : h5get_libversion_f
    use, intrinsic :: ISO_C_Binding     , only : c_char            , c_f_pointer , c_null_char
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//), varying_string
    use            :: String_Handling   , only : operator(//)
    implicit none
    type     (varying_string   )                        :: Output_Build_String
    character(kind=c_char,len=1), dimension(:), pointer :: charVersionString
    integer                     , parameter             :: versionStringLengthMaximum=10
    integer                                             :: hdfError                     , hdfVersionMajor   , &
         &                                                 hdfVersionMinor              , hdfVersionRelease , &
         &                                                 iChr
    type     (c_ptr            )                        :: charVersionPointer
    type     (varying_string   )                        :: CCOMPILER                    , CCOMPILER_VERSION , &
         &                                                 CFLAGS                       , CPPCOMPILER       , &
         &                                                 CPPCOMPILER_VERSION          , CPPFLAGS          , &
         &                                                 FCCOMPILER                   , FCCOMPILER_VERSION, &
         &                                                 FCFLAGS                      , FCFLAGS_NOOPT     , &
         &                                                 LIBS                         , PREPROCESSOR      , &
         &                                                 versionString

    ! Include build environment definitions.
    include 'output.build.environment.inc' ! NO_USES
    ! Initialize an empty string.
    Output_Build_String=""
    ! Write GSL library version string.
    charVersionPointer=GSL_Get_Version()
    call c_f_pointer(charVersionPointer,charVersionString,[versionStringLengthMaximum])
    versionString=""
    iChr=0
    do while (iChr < versionStringLengthMaximum)
      iChr=iChr+1
      if (charVersionString(iChr) == c_null_char) exit
      versionString=versionString//charVersionString(iChr)
    end do
    Output_Build_String=Output_Build_String//":GSL_version["//versionString//"]"
    ! Write FoX library version string.
    Output_Build_String=Output_Build_String//":FoX_version["//Fox_Version//"]"
    ! Write HDF5 library version string.
    call h5get_libversion_f(hdfVersionMajor,hdfVersionMinor,hdfVersionRelease,hdfError)
    if (hdfError /= 0) call Error_Report('unable to get HDF5 library version number'//{introspection:location})
    versionString=''
    versionString=versionString//hdfVersionMajor//"."
    versionString=versionString//hdfVersionMinor//"."
    versionString=versionString//hdfVersionRelease
    Output_Build_String=Output_Build_String//":HDF5_version["//versionString//"]"
    ! Write Make environment variables.
    Output_Build_String=Output_Build_String//":FCCOMPILER["         //FCCOMPILER         //"]"
    Output_Build_String=Output_Build_String//":PREPROCESSOR["       //PREPROCESSOR       //"]"
    Output_Build_String=Output_Build_String//":CCOMPILER["          //CCOMPILER          //"]"
    Output_Build_String=Output_Build_String//":CPPCOMPILER["        //CPPCOMPILER        //"]"
    Output_Build_String=Output_Build_String//":FCFLAGS["            //FCFLAGS            //"]"
    Output_Build_String=Output_Build_String//":FCFLAGS_NOOPT["      //FCFLAGS_NOOPT      //"]"
    Output_Build_String=Output_Build_String//":CFLAGS["             //CFLAGS             //"]"
    Output_Build_String=Output_Build_String//":CPPFLAGS["           //CPPFLAGS           //"]"
    Output_Build_String=Output_Build_String//":FCCOMPILER_VERSION[" //FCCOMPILER_VERSION //"]"
    Output_Build_String=Output_Build_String//":CCOMPILER_VERSION["  //CCOMPILER_VERSION  //"]"
    Output_Build_String=Output_Build_String//":CPPCOMPILER_VERSION["//CPPCOMPILER_VERSION//"]"
    return
  end function Output_Build_String
  
  !![
  <outputFileOpen function="Output_Build_Output"/>
  !!]
  subroutine Output_Build_Output()
    !!{
    Output build information to the main output file.
    !!}
    use            :: File_Utilities    , only : File_Exists
    use            :: FoX_Common        , only : Fox_Version
    use            :: Error             , only : Error_Report
    use            :: Output_HDF5       , only : outputFile
    use            :: Input_Paths       , only : inputPath         , pathTypeExec
    use            :: HDF5              , only : h5get_libversion_f
    use            :: HDF5_Access       , only : hdf5Access
    use            :: IO_HDF5           , only : hdf5Object
    use, intrinsic :: ISO_C_Binding     , only : c_char            , c_f_pointer , c_null_char
    use            :: ISO_Varying_String, only : assignment(=)     , char        , operator(//), operator(/=), &
          &                                      varying_string
    use            :: String_Handling   , only : operator(//)
    implicit none
    character(kind=c_char,len=1), dimension(:), pointer :: charVersionString
    type     (varying_string   ), dimension(1)          :: changeSet
    integer                     , parameter             :: versionStringLengthMaximum=10
    type     (hdf5Object       )                        :: buildGroup
    integer                                             :: hdfError                     , hdfVersionMajor   , &
         &                                                 hdfVersionMinor              , hdfVersionRelease , &
         &                                                 iChr
    type     (c_ptr            )                        :: charVersionPointer
    type     (varying_string   )                        :: CCOMPILER                    , CCOMPILER_VERSION , &
         &                                                 CFLAGS                       , CPPCOMPILER       , &
         &                                                 CPPCOMPILER_VERSION          , CPPFLAGS          , &
         &                                                 FCCOMPILER                   , FCCOMPILER_VERSION, &
         &                                                 FCFLAGS                      , FCFLAGS_NOOPT     , &
         &                                                 LIBS                         , PREPROCESSOR      , &
         &                                                 versionString

    ! Include build environment definitions.
    include 'output.build.environment.inc' ! NO_USES

    ! Create a group for build information.
    !$ call hdf5Access%set()
    buildGroup=outputFile%openGroup('Build','Build information for this model.')

    ! Write GSL library version string.
    charVersionPointer=GSL_Get_Version()
    call c_f_pointer(charVersionPointer,charVersionString,[versionStringLengthMaximum])
    versionString=""
    iChr=0
    do while (iChr < versionStringLengthMaximum)
      iChr=iChr+1
      if (charVersionString(iChr) == c_null_char) exit
      versionString=versionString//charVersionString(iChr)
    end do
    call buildGroup%writeAttribute(versionString,'GSL_library_version' )

    ! Write FoX library version string.
    call buildGroup%writeAttribute(Fox_Version,'FoX_library_version' )

    ! Write HDF5 library version string.
    call h5get_libversion_f(hdfVersionMajor,hdfVersionMinor,hdfVersionRelease,hdfError)
    if (hdfError /= 0) call Error_Report('unable to get HDF5 library version number'//{introspection:location})
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
    call buildGroup%writeAttribute(FCFLAGS            ,'make_FCFLAGS            ')
    call buildGroup%writeAttribute(FCFLAGS_NOOPT      ,'make_FCFLAGS_NOOPT      ')
    call buildGroup%writeAttribute(CFLAGS             ,'make_CFLAGS             ')
    call buildGroup%writeAttribute(CPPFLAGS           ,'make_CPPFLAGS           ')
    call buildGroup%writeAttribute(FCCOMPILER_VERSION ,'make_FCCOMPILER_VERSION ')
    call buildGroup%writeAttribute(CCOMPILER_VERSION  ,'make_CCOMPILER_VERSION  ')
    call buildGroup%writeAttribute(CPPCOMPILER_VERSION,'make_CPPCOMPILER_VERSION')

    ! Add Mercurial changeset information.
    if (File_Exists(inputPath(pathTypeExec)//BUILDPATH//"/galacticus.git.patch")) then
       call changeSet(1)%loadFromFile(char(inputPath(pathTypeExec)//BUILDPATH//'/galacticus.git.patch'))
       if (changeSet(1) /= "" ) call buildGroup%writeDataset(changeSet,'sourceChangeSetDiff','Output of "git diff" - gives the uncommitted source changeset')
       call changeSet(1)%destroy()
    end if
    if (File_Exists(inputPath(pathTypeExec)//BUILDPATH//"/galacticus.git.bundle")) then
       call changeSet(1)%loadFromFile(char(inputPath(pathTypeExec)//BUILDPATH//'/galacticus.git.bundle'))
       if (changeSet(1) /= "" ) call buildGroup%writeDataset(changeSet,'sourceChangeSetBundle','Output of "git bundle HEAD ^origin" - gives changesets not in the remote repo')
       call changeSet(1)%destroy()
    end if

    ! Close the build group.
    call    buildGroup%close()
    !$ call hdf5Access%unset()
   return
  end subroutine Output_Build_Output

end module Output_Build
