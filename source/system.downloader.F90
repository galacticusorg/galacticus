!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which downloads content from a supplied URL.
!!}

module System_Download
  !!{
  Downloads content from a supplied URL.
  !!}
  implicit none
  private
  public :: download

  interface download
     module procedure downloadChar
     module procedure downloadVarStr
  end interface download

  ! Available downloaders.
  logical :: downloadUsingWget  =.false., downloadUsingCurl=.false., &
       &     downloadInitialized=.false.
  
contains

  subroutine downloadInitialize()
    !!{
    Determine which downloaders are available.
    !!}
    use :: Galacticus_Error, only : errorStatusSuccess
    use :: System_Command  , only : System_Command_Do
    implicit none
    integer :: status

    if (.not.downloadInitialized) then
       !$omp critical(downloadInitialize)
       if (.not.downloadInitialized) then
          call System_Command_Do("which wget > /dev/null 2>&1",status)
          downloadUsingWget=status == errorStatusSuccess
          call System_Command_Do("which curl > /dev/null 2>&1",status)
          downloadUsingCurl=status == errorStatusSuccess
          downloadInitialized=.true.
       end if
       !$omp end critical(downloadInitialize)
    end if    
    return
  end subroutine downloadInitialize
  
  subroutine downloadVarStr(url,outputFileName,status)
    !!{
    Download content from the given {\normalfont url} to the given {\normalfont \ttfamily outputFileName}.
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type   (varying_string), intent(in   )           :: url   , outputFileName
    integer                , intent(  out), optional :: status
    
    call download(char(url),char(outputFileName),status)
    return
  end subroutine downloadVarStr

  subroutine downloadChar(url,outputFileName,status)
    !!{
    Download content from the given {\normalfont url} to the given {\normalfont \ttfamily outputFileName}.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report, errorStatusSuccess, errorStatusFail
    use :: System_Command  , only : System_Command_Do
    implicit none
    character(len=*), intent(in   )           :: url    , outputFileName
    integer         , intent(  out), optional :: status
    integer                                   :: status_
    
    call downloadInitialize()
    status_=errorStatusFail
    if      (downloadUsingWget) then
       call System_Command_Do('wget --no-check-certificate "'//trim(url)//'" -O '      //trim(outputFileName),status_)
    else if (downloadUsingCurl) then
       call System_Command_Do('curl --insecure --location "' //trim(url)//'" --output '//trim(outputFileName),status_)
    else if (.not.present(status)) then
       call Galacticus_Error_Report('no downloader available'//{introspection:location})
    end if
    if (present(status)) then
       status=status_
    else if (status_ /= 0) then
       call Galacticus_Error_Report('failed to download "'//trim(url)//'"'//{introspection:location})
    end if
    return
  end subroutine downloadChar

end module System_Download
