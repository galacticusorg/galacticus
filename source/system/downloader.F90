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

!+    The fallback from `wget` to `curl` within a single download attempt, and the use of `--ciphers=DEFAULT` to avoid hosts which
!+    reject `wget` on the basis of its TLS handshake, were diagnosed and drafted with assistance from Claude, and reviewed and
!+    verified by Andrew Benson.

!!{RST
Contains a module which downloads content from a supplied URL (or list of URLs).
!!}

module System_Download
  !!{RST
  Downloads content from a supplied URL. A scalar URL may be provided, or a 1D array of URLs to be tried in turn---if the download from one URL fails, the next is used as a fallback.

  Downloads verify the TLS certificate of the server unless the optional ``allowInsecure`` argument is set to true---see
  ``downloadMultiple`` for the conditions under which that is permissible.
  !!}
  implicit none
  private
  public :: download, downloadCommand, downloaderWget, downloaderCurl

  ! Identifiers for the supported downloader tools.
  integer, parameter :: downloaderWget=1, downloaderCurl=2

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
    !!{RST
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

  subroutine downloadVarStrVarStr(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the given url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   )           :: url          , outputFileName
    integer                , intent(in   ), optional :: retries      , retryWait     , &
         &                                              timeout
    logical                , intent(in   ), optional :: allowInsecure
    integer                , intent(  out), optional :: status
    type   (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadVarStrVarStr

  subroutine downloadVarStrChar(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the given url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type     (varying_string), intent(in   )           :: url
    character(len=*         ), intent(in   )           :: outputFileName
    integer                  , intent(in   ), optional :: retries       , retryWait, &
         &                                                timeout
    logical                  , intent(in   ), optional :: allowInsecure
    integer                  , intent(  out), optional :: status
    type     (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadVarStrChar

  subroutine downloadCharVarStr(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the given url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : char, var_str, varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   )           :: url
    type     (varying_string), intent(in   )           :: outputFileName
    integer                  , intent(in   ), optional :: retries       , retryWait, &
         &                                                timeout
    logical                  , intent(in   ), optional :: allowInsecure
    integer                  , intent(  out), optional :: status
    type     (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadCharVarStr

  subroutine downloadCharChar(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the given url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   )           :: url          , outputFileName
    integer                  , intent(in   ), optional :: retries      , retryWait     , &
         &                                                timeout
    logical                  , intent(in   ), optional :: allowInsecure
    integer                  , intent(  out), optional :: status
    type     (varying_string), dimension(1)            :: urls

    urls(1)=url
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadCharChar

  subroutine downloadVarStrArrayVarStr(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the first available URL in url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   ), dimension(:) :: url
    type   (varying_string), intent(in   )               :: outputFileName
    integer                , intent(in   ), optional     :: retries       , retryWait, &
         &                                                  timeout
    logical                , intent(in   ), optional     :: allowInsecure
    integer                , intent(  out), optional     :: status

    call downloadMultiple(url,char(outputFileName),retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadVarStrArrayVarStr

  subroutine downloadVarStrArrayChar(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the first available URL in url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type     (varying_string), intent(in   ), dimension(:) :: url
    character(len=*         ), intent(in   )               :: outputFileName
    integer                  , intent(in   ), optional     :: retries       , retryWait, &
         &                                                    timeout
    logical                  , intent(in   ), optional     :: allowInsecure
    integer                  , intent(  out), optional     :: status

    call downloadMultiple(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadVarStrArrayChar

  subroutine downloadCharArrayVarStr(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the first available URL in url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : char, varying_string, assignment(=)
    implicit none
    character(len=*         ), intent(in   ), dimension(:        ) :: url
    type     (varying_string), intent(in   )                       :: outputFileName
    integer                  , intent(in   ), optional             :: retries       , retryWait, &
         &                                                            timeout
    logical                  , intent(in   ), optional             :: allowInsecure
    integer                  , intent(  out), optional             :: status
    type     (varying_string)               , dimension(size(url)) :: urls
    integer                                                        :: i

    do i=1,size(url)
       urls(i)=url(i)
    end do
    call downloadMultiple(urls,char(outputFileName),retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadCharArrayVarStr

  subroutine downloadCharArrayChar(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content from the first available URL in url to the given ``outputFileName``.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=)
    implicit none
    character(len=*       ), intent(in   ), dimension(:        ) :: url
    character(len=*       ), intent(in   )                       :: outputFileName
    integer                , intent(in   ), optional             :: retries       , retryWait, &
         &                                                          timeout
    logical                , intent(in   ), optional             :: allowInsecure
    integer                , intent(  out), optional             :: status
    type   (varying_string)               , dimension(size(url)) :: urls
    integer                                                      :: i

    do i=1,size(url)
       urls(i)=url(i)
    end do
    call downloadMultiple(urls,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    return
  end subroutine downloadCharArrayChar

  subroutine downloadMultiple(url,outputFileName,retries,retryWait,timeout,status,allowInsecure)
    !!{RST
    Download content to the given ``outputFileName``, trying each URL in url in turn. If the download from one URL fails (even
    after any retries), the next URL is used as a fallback. The download is considered successful as soon as any URL succeeds.

    Each attempt uses ``wget``, falling back to ``curl`` if ``wget`` is unavailable or fails. Both are tried because a failure of
    one does not establish that the content is unavailable---a host which fingerprints and rejects a client at the TLS layer will
    reject it for every URL and on every retry, while accepting the other client.

    TLS certificates are verified by default. Setting ``allowInsecure`` to true disables that verification, which leaves the
    download open to interception and modification by anyone able to intercept the connection---content fetched this way must be
    treated as untrusted. It exists only for hosts which serve an unusable certificate and for which no alternative source is
    available, and must never be used for content which is subsequently compiled or executed unless its integrity is established
    by some other means.
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
    logical                  , intent(in   ), optional     :: allowInsecure
    integer                  , intent(  out), optional     :: status
    integer                                                :: status_       , tries    , &
         &                                                    i
    type     (varying_string)                              :: urlList       , command
    !![
    <optionalArgument name="retries"       defaultsTo="0"      />
    <optionalArgument name="retryWait"     defaultsTo="60"     />
    <optionalArgument name="timeout"       defaultsTo="300"    />
    <optionalArgument name="allowInsecure" defaultsTo=".false."/>
    !!]

    call downloadInitialize()
    if (present(status)) status=0
    status_=errorStatusFail
    do i=1,size(url)
       tries=0
       do while (tries <= retries_)
          status_=errorStatusFail
          if (downloadUsingWget) then
             command=downloadCommand(url(i),outputFileName,timeout_,allowInsecure_,downloaderWget)
             call System_Command_Do(char(command),status_)
          end if
          ! If `wget` is unavailable, or failed, try `curl`. A failure of one downloader does not imply that the content is
          ! unavailable: hosts which fingerprint and reject a client at the TLS layer do so irrespective of the URL requested, so
          ! the other downloader may well succeed where this one could never do so no matter how many times it is retried.
          if (status_ /= errorStatusSuccess .and. downloadUsingCurl) then
             ! Remove any partial file left behind by a failed `wget` attempt, so that no part of it can survive into the file
             ! written by `curl`.
             if (File_Exists(outputFileName)) call File_Remove(outputFileName)
             command=downloadCommand(url(i),outputFileName,timeout_,allowInsecure_,downloaderCurl)
             call System_Command_Do(char(command),status_)
          end if
          if (.not.downloadUsingWget .and. .not.downloadUsingCurl .and. .not.present(status)) &
               & call Error_Report('no downloader available'//{introspection:location})
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

  function downloadCommand(url,outputFileName,timeout,allowInsecure,downloader) result(command)
    !!{RST
    Build the shell command used to download ``url`` to ``outputFileName`` using the tool identified by ``downloader`` (one of
    ``downloaderWget`` or ``downloaderCurl``), allowing at most ``timeout`` seconds for the attempt.

    Both the URL and the output file name are shell-escaped, so that each is passed to the tool as a single, literal token: URLs
    and paths can contain characters which are special to the shell (and, for URLs built from parameter or dataset content, are
    not necessarily under our control), so neither may be interpolated unescaped. The escaped URL is preceded by ``--``, marking
    the end of options, so that a URL beginning with ``-`` can not be interpreted as an option by the tool.

    The tool is instructed to verify the server's TLS certificate unless ``allowInsecure`` is true---see ``downloadMultiple``.

    This is separated from ``downloadMultiple`` so that the constructed command can be tested directly.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : varying_string, char, operator(//), assignment(=)
    use :: System_Command    , only : shellEscape
    implicit none
    type     (varying_string)                :: command
    type     (varying_string), intent(in   ) :: url
    character(len=*         ), intent(in   ) :: outputFileName
    integer                  , intent(in   ) :: timeout          , downloader
    logical                  , intent(in   ) :: allowInsecure
    type     (varying_string)                :: escapedOutputFile, escapedURL
    character(len=12        )                :: timeoutLabel

    ! Build a string representation of the per-attempt timeout (in seconds).
    write (timeoutLabel,'(i0)') timeout
    ! Escape both arguments. Note that the result of `shellEscape` is assigned to a local variable before use, as required.
    escapedURL       =url
    escapedURL       =shellEscape(escapedURL          )
    escapedOutputFile=shellEscape(trim(outputFileName))
    select case (downloader)
    case (downloaderWget)
       ! Force `wget` to make only a single attempt (its default is `--tries=20`). Retries are handled by the calling loop in
       ! `downloadMultiple`, so allowing `wget` to also retry internally results in a multiplicative number of attempts (and can
       ! cause the download to far exceed any time limit when each internal attempt hangs until its read-timeout). The `--timeout`
       ! option bounds the time spent on DNS lookup, connection, and reads for the single attempt.
       !
       ! `--ciphers=DEFAULT` replaces `wget`'s own, more restrictive, default cipher list with that of the underlying TLS library.
       ! `wget`'s default list yields a TLS handshake which some content delivery networks (bitbucket.org, which hosts the
       ! NGenHalofit source, among them) fingerprint as a bot and reject---returning a "404 Not Found" for every URL on the host,
       ! including its front page, so that the failure masquerades as missing content rather than a refused client. This does not
       ! weaken certificate verification, which remains controlled by `--no-check-certificate` below. A `wget` too old to accept
       ! `--ciphers` will reject the option and exit without contacting the host---the fallback to `curl` in `downloadMultiple`
       ! covers that case, so downloads still succeed wherever `curl` is available.
       command='wget --tries=1 --ciphers=DEFAULT --timeout='//trim(timeoutLabel)//' '
       if (allowInsecure) command=command//'--no-check-certificate '
       command=command//'-O '//escapedOutputFile//' -- '//escapedURL
    case (downloaderCurl)
       ! Force `curl` to make only a single attempt (i.e. disable its own retrying) so that retries are handled solely by the
       ! calling loop, consistent with the behavior of `wget` above. The `--max-time` option bounds the total time allowed for the
       ! single attempt.
       !
       ! `--fail` is essential: without it `curl` treats an HTTP error response as a successful transfer, exiting with zero status
       ! after writing the server's error page to the output file. The download would then be reported as having succeeded, and
       ! the error page used as though it were the requested content.
       command='curl --location --retry 0 --fail --max-time '//trim(timeoutLabel)//' '
       if (allowInsecure) command=command//'--insecure '
       command=command//'--output '//escapedOutputFile//' -- '//escapedURL
    case default
       command=''
       call Error_Report('unknown downloader'//{introspection:location})
    end select
    return
  end function downloadCommand

end module System_Download
