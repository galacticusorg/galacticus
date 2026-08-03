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

!!{RST
Contains a program to test construction of the commands used to download content.
!!}

program Test_System_Download
  !!{RST
  Tests the construction of the commands used to download content. Two properties are checked:

  #. TLS certificates are verified unless verification is explicitly disabled by the caller; and
  #. the URL is passed to the downloader as a single, literal token, so that a URL containing characters which are special to
     the shell can not result in execution of arbitrary commands.

  The second property is checked both by inspecting the constructed command, and by actually executing it with URLs carrying
  command-injection payloads and verifying that the injected command did not run. The URLs used for that check use the reserved
  ``.invalid`` top-level domain (RFC 6761), which is guaranteed not to resolve, so the download itself always fails quickly and
  without network access. If the downloader is not installed the command simply fails, which can not cause a spurious failure of
  the test.
  !!}
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: File_Utilities    , only : File_Exists        , File_Remove
  use :: ISO_Varying_String, only : varying_string     , char                  , assignment(=)       , operator(//)
  use :: System_Command    , only : System_Command_Do
  use :: System_Download   , only : downloadCommand    , downloaderWget        , downloaderCurl
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type     (varying_string)                          :: url           , command
  type     (varying_string), dimension(2)            :: payload
  character(len=29        ), parameter               :: sentinel      ='downloadInjectionSentinel.tmp'
  character(len=16        ), parameter               :: outputFile    ='downloadTest.tmp'
  integer                  , parameter, dimension(2) :: downloader    =[downloaderWget,downloaderCurl]
  character(len=4         ), parameter, dimension(2) :: downloaderName=['wget','curl']
  character(len=20        ), parameter, dimension(2) :: payloadName   =['command substitution','quote break-out     ']
  integer                                            :: i             , j     , &
       &                                                status

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("System download")

  ! Certificate verification must be enabled unless explicitly disabled by the caller.
  url='https://example.org/data.txt'
  command=downloadCommand(url,trim(outputFile),300,.false.,downloaderWget)
  call Assert('wget verifies certificates by default'  ,index(char(command),'--no-check-certificate') == 0,.true.)
  command=downloadCommand(url,trim(outputFile),300,.false.,downloaderCurl)
  call Assert('curl verifies certificates by default'  ,index(char(command),'--insecure'            ) == 0,.true.)
  ! ... and must be disabled when, and only when, the caller asks for that.
  command=downloadCommand(url,trim(outputFile),300,.true. ,downloaderWget)
  call Assert('wget verification disabled on request'  ,index(char(command),'--no-check-certificate') >  0,.true.)
  command=downloadCommand(url,trim(outputFile),300,.true. ,downloaderCurl)
  call Assert('curl verification disabled on request'  ,index(char(command),'--insecure'            ) >  0,.true.)

  ! The URL must be quoted, and preceded by an end-of-options marker so that a URL beginning with "-" is not read as an option.
  command=downloadCommand(url,trim(outputFile),300,.false.,downloaderWget)
  call Assert('URL is quoted and follows "--"'         ,index(char(command),"-- 'https://example.org/data.txt'") > 0,.true.)

  ! A URL containing shell metacharacters must not result in execution of the injected command. The first payload uses command
  ! substitution, which is active inside double quotes (as used by an earlier implementation of this module); the second closes a
  ! single quote before injecting a command, and so targets the single-quoting used now.
  payload(1)='https://invalid.invalid/x$(touch '//trim(sentinel)//')'
  payload(2)="https://invalid.invalid/x';touch "//trim(sentinel)//";'"
  do i=1,size(payload)
     do j=1,size(downloader)
        if (File_Exists(trim(sentinel))) call File_Remove(trim(sentinel))
        url=payload(i)
        command=downloadCommand(url,trim(outputFile),5,.false.,downloader(j))
        call System_Command_Do(char(command),status)
        call Assert(                                                                                            &
             &      'no command injection via URL ['//trim(downloaderName(j))//', '//trim(payloadName(i))//']', &
             &      File_Exists(trim(sentinel))                                                               , &
             &      .false.                                                                                     &
             &     )
        if (File_Exists(trim(sentinel  ))) call File_Remove(trim(sentinel  ))
        if (File_Exists(trim(outputFile))) call File_Remove(trim(outputFile))
     end do
  end do

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_System_Download
