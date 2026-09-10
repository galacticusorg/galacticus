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

  That the handler is passed ``signalNone`` is checked too: a handler which interprets the signal
  must be able to distinguish "no signal" from a real one, and must not attempt to decode a value
  which is not a member of the signal enumeration.

  The handler exits the process, which is what makes this testable. ``Error_Report`` aborts, so
  nothing this program writes after calling it is ever reached, and the only channel the test
  harness observes is the exit status and the output. Exiting from the handler puts "the handler
  ran" into that exit status: reaching the end of ``Error_Report`` instead means it did not run, and
  the process exits non-zero, failing the test.
  !!}
  use, intrinsic :: ISO_Fortran_Env, only : output_unit
  use            :: Display        , only : displayVerbositySet   , verbosityLevelStandard
  use            :: Error          , only : Error_Report          , signalHandlerInterface, signalHandlerRegister, signalNone
  implicit none
  procedure(signalHandlerInterface), pointer :: handler

  call displayVerbositySet   (verbosityLevelStandard)
  handler                    => testHandler
  call signalHandlerRegister (handler               )
  ! Report a fatal error. The handler must run, and must exit the process; if it does not run,
  ! `Error_Report` exits non-zero and the test fails.
  call Error_Report          ('deliberate error, to test that registered handlers are called'//{introspection:location})

contains

  subroutine testHandler(signal)
    !!{RST
    Report whether this handler was called as expected, and exit.
    !!}
    implicit none
    integer, intent(in   ) :: signal

    if (signal == signalNone) then
       write (output_unit,'(a)'    ) 'SUCCESS: the registered handler was called by `Error_Report`, and passed `signalNone`'
    else
       write (output_unit,'(a,i0)' ) 'FAILED: the registered handler was called with a signal number instead of `signalNone`: ',signal
    end if
    call Flush(output_unit)
    call Exit (0          )
    return
  end subroutine testHandler

end program Test_Error_Report_Handlers
