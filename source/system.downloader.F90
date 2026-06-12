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
Contains a module which downloads content from a supplied URL (or list of URLs).
!!}

module System_Download
  !!{
  Downloads content from a supplied URL. A scalar URL may be provided, or a 1D array of URLs to be tried in turn---if the
  download from one URL fails, the next is used as a fallback.
  !!}
  implicit none
  private
  public :: download

  interface download
     module procedure downloadCharChar
     module procedure downloadVarStrVarStr
     module procedure downloadVarStrChar
     module procedure downloadCharVarStr
     module procedure downloadCharArrayChar
     module procedure downloadCharArrayVarStr
     module procedure downloadVarStrArrayChar
     module procedure downloadVarStrArrayVarStr
  end interface download

  ! Available downloaders.
  logical :: downloadUsingWget  =.false., downloadUsingCurl=.false., &
       &     downloadInitialized=.false.

contains

  subroutine downloadInitialize()
    !!{
    Determine which downloaders are available.
    !!}
    use :: Error         , only : errorStatusSuccess
    use :: System_Command, only : System_Command_Do
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

  subroutine downloadVarStrVarStr(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the given {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   )           :: url    , outputFileName
    integer                , intent(in   ), optional :: retries, retryWait     , &
         &                                              timeout
    integer                , intent(  out), optional :: status
    type   (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status)
    return
  end subroutine downloadVarStrVarStr

  subroutine downloadVarStrChar(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the given {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type     (varying_string), intent(in   )           :: url
    character(len=*         ), intent(in   )           :: outputFileName
    integer                  , intent(in   ), optional :: retries       , retryWait, &
         &                                                timeout
    integer                  , intent(  out), optional :: status
    type   (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status)
    return
  end subroutine downloadVarStrChar

  subroutine downloadCharVarStr(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the given {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : char, var_str, varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   )           :: url
    type     (varying_string), intent(in   )           :: outputFileName
    integer                  , intent(in   ), optional :: retries       , retryWait, &
         &                                                timeout
    integer                  , intent(  out), optional :: status
    type     (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status)
    return
  end subroutine downloadCharVarStr

  subroutine downloadCharChar(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the given {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   )           :: url    , outputFileName
    integer                  , intent(in   ), optional :: retries, retryWait     , &
         &                                                timeout
    integer                  , intent(  out), optional :: status
    type     (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status)
    return
  end subroutine downloadCharChar

  subroutine downloadVarStrArrayVarStr(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the first available URL in {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   ), dimension(:) :: url
    type   (varying_string), intent(in   )               :: outputFileName
    integer                , intent(in   ), optional     :: retries       , retryWait, &
         &                                                  timeout
    integer                , intent(  out), optional     :: status

    call downloadMultiple(url,char(outputFileName),retries,retryWait,timeout,status)
    return
  end subroutine downloadVarStrArrayVarStr

  subroutine downloadVarStrArrayChar(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the first available URL in {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type     (varying_string), intent(in   ), dimension(:) :: url
    character(len=*         ), intent(in   )               :: outputFileName
    integer                  , intent(in   ), optional     :: retries       , retryWait, &
         &                                                    timeout
    integer                  , intent(  out), optional     :: status

    call downloadMultiple(url,outputFileName,retries,retryWait,timeout,status)
    return
  end subroutine downloadVarStrArrayChar

  subroutine downloadCharArrayVarStr(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the first available URL in {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : char, varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   ), dimension(:        ) :: url
    type     (varying_string), intent(in   )                       :: outputFileName
    integer                  , intent(in   ), optional             :: retries       , retryWait, &
         &                                                            timeout
    integer                  , intent(  out), optional             :: status
    type     (varying_string)               , dimension(size(url)) :: urls
    integer                                                        :: i

    do i=1,size(url)
       urls(i)=url(i)
    end do
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status)
    return
  end subroutine downloadCharArrayVarStr

  subroutine downloadCharArrayChar(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content from the first available URL in {\normalfont url} to the given \mono{outputFileName}.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=)
    implicit none
    character(len=*       ), intent(in   ), dimension(:        ) :: url
    character(len=*       ), intent(in   )                       :: outputFileName
    integer                , intent(in   ), optional             :: retries       , retryWait, &
         &                                                          timeout
    integer                , intent(  out), optional             :: status
    type   (varying_string)               , dimension(size(url)) :: urls
    integer                                                      :: i

    do i=1,size(url)
       urls(i)=url(i)
    end do
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status)
    return
  end subroutine downloadCharArrayChar

  subroutine downloadMultiple(url,outputFileName,retries,retryWait,timeout,status)
    !!{
    Download content to the given \mono{outputFileName}, trying each URL in {\normalfont url} in turn. If the download from one
    URL fails (even after any retries), the next URL is used as a fallback. The download is considered successful as soon as any
    URL succeeds.
    !!}
    use :: Error             , only : Error_Report     , errorStatusFail, errorStatusSuccess
    use :: File_Utilities    , only : File_Exists      , File_Remove
    use :: ISO_Varying_String, only : varying_string   , char           , operator(//)      , assignment(=)
    use :: System_Command    , only : System_Command_Do
    implicit none
    type     (varying_string), intent(in   ), dimension(:) :: url
    character(len=*         ), intent(in   )               :: outputFileName
    integer                  , intent(in   ), optional     :: retries       , retryWait, &
         &                                                    timeout
    integer                  , intent(  out), optional     :: status
    integer                                                :: status_       , tries    , i
    type     (varying_string)                              :: urlList
    character(len=12        )                              :: timeoutLabel
    !![
    <optionalArgument name="retries"   defaultsTo="0"  />
    <optionalArgument name="retryWait" defaultsTo="60" />
    <optionalArgument name="timeout"   defaultsTo="300"/>
    !!]

    call downloadInitialize()
    ! Build a string representation of the per-attempt timeout (in seconds) for use in the downloader commands below.
    write (timeoutLabel,'(i0)') timeout_
    if (present(status)) status=0
    status_=errorStatusFail
    do i=1,size(url)
       tries=0
       do while (tries <= retries_)
          status_=errorStatusFail
          if      (downloadUsingWget) then
             ! Force `wget` to make only a single attempt (its default is `--tries=20`). We handle retries ourselves via the loop
             ! here, so allowing `wget` to also retry internally results in a multiplicative number of attempts (and can cause the
             ! download to far exceed any time limit when each internal attempt hangs until its read-timeout). The `--timeout`
             ! option bounds the time spent on DNS lookup, connection, and reads for the single attempt.
             call System_Command_Do('wget --no-check-certificate --tries=1 --timeout='//trim(timeoutLabel)//' "'//char(url(i))//'" -O '      //trim(outputFileName),status_)
          else if (downloadUsingCurl) then
             ! Force `curl` to make only a single attempt (i.e. disable its own retrying) so that retries are handled solely by the
             ! loop here, consistent with the behavior of `wget` above. The `--max-time` option bounds the total time allowed for
             ! the single attempt.
             call System_Command_Do('curl --insecure --location --retry 0 --max-time '//trim(timeoutLabel)//' "' //char(url(i))//'" --output '//trim(outputFileName),status_)
          else if (.not.present(status)) then
             call Error_Report('no downloader available'//{introspection:location})
          end if
          if (status_ == errorStatusSuccess) then
             if (present(status)) status=status_
             return
          end if
          tries=tries+1
          if (File_Exists(outputFileName)) call File_Remove(outputFileName)
          ! Wait before the next attempt, unless this was the final attempt of the final URL.
          if (tries <= retries_ .or. i < size(url)) call sleep(retryWait_)
       end do
    end do
    if (present(status)) status=status_
    if (.not.present(status) .and. status_ /= errorStatusSuccess) then
       urlList=''
       do i=1,size(url)
          if (i > 1) urlList=urlList//', '
          urlList=urlList//char(url(i))
       end do
       call Error_Report('failed to download from "'//char(urlList)//'"'//{introspection:location})
    end if
    return
  end subroutine downloadMultiple

end module System_Download
