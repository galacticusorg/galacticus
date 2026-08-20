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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test that registered signal handlers are called by ``Error_Report``.
!!}

program Test_Error_Report_Handlers
  !!{RST
  Tests that handlers registered with ``signalHandlerRegister`` are called when a fatal error is
  reported, and not only when a signal is caught.

  Handlers exist so that a failure can dump context which would otherwise be lost---the
  posterior-sampling likelihood, for example, registers one which writes the parameter set being
  evaluated to a file. A deliberate fatal error is far more common than a crash, so a handler which
  ran only on a signal discarded that context in precisely the usual case.

  The handler is passed ``signalNone`` when called other than in response to a signal. This is
  checked too: a handler which switches on the signal number must be able to distinguish the two,
  and one which decodes it must not attempt to decode a value which is not a signal.

  The program necessarily aborts, since that is what ``Error_Report`` does, so the handler records
  what it saw in a file which ``test-error-report-handlers.py`` then checks. Reporting the result
  from the program itself is not possible: nothing it writes after ``Error_Report`` is reached.
  !!}
  use :: Display, only : displayVerbositySet   , verbosityLevelStandard
  use :: Error  , only : Error_Report          , signalHandlerInterface, signalHandlerRegister, signalNone
  implicit none
  procedure(signalHandlerInterface), pointer :: handler
  character(len=*                 ), parameter :: fileNameRecord='testSuite/outputs/errorReportHandlers.record'

  call displayVerbositySet   (verbosityLevelStandard)
  handler                    => testHandler
  call signalHandlerRegister (handler               )
  ! Abort. The handler must run before the process exits, and must be passed `signalNone`.
  call Error_Report          ('deliberate error, to test that registered handlers are called'//{introspection:location})

contains

  subroutine testHandler(signal)
    !!{RST
    Record the signal with which this handler was called.
    !!}
    implicit none
    integer, intent(in   ) :: signal
    integer                :: fileRecord

    open (newunit=fileRecord,file=fileNameRecord,status='unknown',form='formatted')
    write (fileRecord,'(a,i0)') 'signal=',signal
    write (fileRecord,'(a,l1)') 'isSignalNone=',signal == signalNone
    close(fileRecord)
    return
  end subroutine testHandler

end program Test_Error_Report_Handlers
