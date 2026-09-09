.. _manual-sec-errorHandling:

Error Handling
==============

Galacticus distinguishes two kinds of failure, and the distinction determines
almost everything else on this page:

* A **fatal error** is a condition from which the run cannot continue: a
  parameter file that asks for something impossible, a data file that is not
  where it should be, an internal inconsistency that indicates a bug. These are
  reported with ``Error_Report``, which prints the message and aborts.

* A **recoverable failure** is one that a caller further up the stack is
  expected to handle. The canonical case is the failure to evolve a single
  merger tree: with ``<tolerateFailures value="true"/>`` that tree is abandoned
  and the remaining trees are evolved as usual. These are reported by setting an
  optional ``status`` argument, and returning.

Most code needs only ``Error_Report``. Adding a ``status`` argument is
worthwhile only where a caller genuinely exists that can do something other than
abort.

.. _manual-sec-errorHandlingFatal:

Reporting a fatal error
-----------------------

``Error_Report`` takes the message to report, to which the location at which the
error was raised **must** be appended:

.. code-block:: fortran

   use :: Error, only : Error_Report
   ...
   if (radius <= 0.0d0) call Error_Report('non-positive radius'//{introspection:location})

``{introspection:location}`` is expanded by the preprocessor into the name of
the enclosing function or subroutine, the source file, and the line number.
Without it the user is left with only the message and a backtrace into the
*preprocessed* sources under ``work/build/``, whose line numbers do not
correspond to those of the file being edited.

This is checked mechanically. Check 5 of ``scripts/aux/staticAnalyzer.py``
reports any ``Error_Report`` call which builds its message from a string literal
and does not append the location; it runs both in the *Fortran-Static-Analysis*
CI job and in the pre-commit hook. Only literal-built messages are checked,
because a message passed as a bare variable may have had the location appended
where that variable was assigned:

.. code-block:: fortran

   ! Checked - the message is a literal, so the location must be appended here.
   call Error_Report('unrecognized option'//{introspection:location})

   ! Not checked - the location may already be part of `message`.
   call Error_Report(message)

Writing a good message
~~~~~~~~~~~~~~~~~~~~~~

The message should say what went wrong in terms the user can act on. Two
conventions help:

**Name the values involved.** An error that reports only that something was
wrong, without saying what was expected or what was found, forces the reader
into the source to learn anything at all. Where an object is not of the expected
class, name both classes: every ``functionClass`` has an ``objectType`` method
returning the name of its concrete class, and it is declared on the
``functionClass`` base type, so the name can be obtained through any
``class(functionClass)`` reference:

.. code-block:: fortran

   select type (self)
   type is (galacticStructureSolverLinear)
      call self%solve(node)
   class is (functionClass)
      call Error_Report('object is not of [galacticStructureSolverLinear] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
   class default
      call Error_Report('object is not of [galacticStructureSolverLinear] class'//{introspection:location})
   end select

The ``class default`` branch is kept as the genuine catch-all, so an object
which is not a ``functionClass`` at all still raises an error rather than
falling through.

**Say how to fix it.** Where a remedy exists, append a ``HELP:`` hint. The
convention is the literal text ``HELP:`` in green, followed by the advice:

.. code-block:: fortran

   use :: Display, only : displayGreen, displayReset
   ...
   call Error_Report(                                                                         &
        &            'Unable to find data file "'//char(fileName)//'"'//char(10)           // &
        &            displayGreen()//'HELP:'//displayReset()                               // &
        &            ' this file is provided by the Galacticus datasets repository. If the'// &
        &            ' path above begins with "./" then the `GALACTICUS_DATA_PATH`'        // &
        &            ' environment variable is not set'                                    // &
        &            {introspection:location}                                                 &
        &           )

A hint is most valuable on the errors users meet most often --- a parameter that
is absent, a value out of range, a data file that cannot be found --- and least
valuable on internal consistency checks, which a user cannot act on in any case.
Where the error corresponds to a section of the
:doc:`troubleshooting guide <../user-guide/troubleshooting/run-time-errors>`, phrase the
message so that the two agree: a user searching for the message text should find
the section.

.. _manual-sec-errorHandlingStatus:

Reporting a recoverable failure
-------------------------------

A procedure which can fail recoverably takes an ``optional``, ``intent(out)``
``status`` argument. The convention is that the procedure behaves as it always
did when ``status`` is absent --- that is, it reports a fatal error --- and
returns a status code instead when it is present. The canonical form is:

.. code-block:: fortran

   if (present(status)) then
      call displayMessage('ODE integration failed')
      status=errorStatusUnderflow
      return
   else
      call Error_Report('ODE integration failed '//{introspection:location})
   end if

Two properties of this idiom matter. Existing callers are unaffected, since they
do not pass ``status``; and a failure is never silent, because the caller which
*does* pass ``status`` has explicitly taken responsibility for it. A procedure
with a ``status`` argument should set it to ``errorStatusSuccess`` on entry, so
that every path which returns normally reports success without having to say so.

The status codes are declared in ``source/error/_module.F90``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Code
     - Meaning
   * - ``errorStatusSuccess``
     - No error.
   * - ``errorStatusFail``
     - Generic failure, where no more specific code applies.
   * - ``errorStatusInputDomain``
     - An argument was outside the domain on which the procedure is defined.
   * - ``errorStatusRound``
     - Failure due to rounding error.
   * - ``errorStatusOutOfRange``
     - The result was outside the representable range.
   * - ``errorStatusDivideByZero``
     - Division by zero.
   * - ``errorStatusUnderflow``
     - Floating point underflow.
   * - ``errorStatusMaxIterations``
     - An iterative method did not converge within its iteration limit.
   * - ``errorStatusXCPU``
     - A CPU or wall time limit was exceeded.
   * - ``errorStatusNotExist``
     - The entity requested does not exist.

Those which have a natural GSL counterpart are defined *as* the corresponding
GSL code, so that a status returned by a GSL routine can be propagated without
translation. The remainder take values above 1024, chosen so that they cannot
collide with a GSL code.

Passing ``status`` onward
~~~~~~~~~~~~~~~~~~~~~~~~~

A procedure which itself calls a procedure with an optional ``status`` argument
must pass its own ``status`` on only when it has one. Passing an absent optional
argument onward is legal, but the two-branch form is often clearer, and where the
call has several such arguments the combinations multiply. The
:ref:`conditionalCall directive <manual-sec-conditionalCall>` generates the
combinations from a single written form, and exists largely for this case.

.. _manual-sec-errorHandlingWarnings:

Warnings
--------

``Warn`` reports something the user should know about but which does not stop the
run:

.. code-block:: fortran

   use :: Error, only : Warn
   ...
   call Warn('tabulated data required extrapolation beyond its range')

A warning behaves differently depending on the verbosity in force. At
``verbosityLevel="warn"`` or above it is displayed as it is issued. Below that
level it is *accumulated* instead, and the accumulated warnings are replayed if
the run subsequently ends in a fatal error --- so that a user who sees only the
error still sees the warnings that preceded it.

That replay is a bounded resource, and this is the constraint that governs what
should be a warning at all. Warnings are deduplicated on their message text, and
at most ``warningsUniqueMaximum`` (100) distinct messages are retained; once the
limit is reached, further *distinct* warnings are counted but not recorded.

.. warning::

   Do not use ``Warn`` to report something which happens as a matter of course.
   Because the retained set is filled in the order warnings are issued, a routine
   condition arising during start-up will exhaust the limit before evolution
   begins, and genuine warnings issued later can then never appear in the replay
   at all.

   Reporting the use of a default parameter value was once such a case: a small
   model defaults of order a hundred parameters, which filled the replay
   completely. That reporting was removed, and the information is instead
   recorded in the output file, where each defaulted parameter carries a
   ``{defaulted}`` marker alongside its value in the ``Parameters`` group.

The test for something that belongs in a warning is whether it indicates that
the *model* may not be what the user intended. Anything which is merely a record
of what the code did belongs in the output file, or at
``verbosityLevelInfo``, not in the warning channel.

.. _manual-sec-errorHandlingAllocation:

Allocation failures
-------------------

An ``allocate`` which fails without ``stat=`` aborts the process with no
indication of what was being allocated. For the large allocations --- those
whose size is set by a count read at run time, rather than by a fixed dimension
--- pass ``stat=`` and report:

.. code-block:: fortran

   integer :: allocErr
   ...
   allocate(self,stat=allocErr)
   if (allocErr /= 0) call Error_Report('unable to allocate node'//{introspection:location})

This costs a scalar integer and one branch. Note what it deliberately does *not*
do: it does not build a message reporting the size requested, because doing so
requires a ``type(varying_string)`` local, which is then constructed and
destroyed on every call rather than only on the failure that never normally
happens. Where the size is worth reporting, it belongs in a ``block`` or a helper
called only on the failing branch.

.. _manual-sec-errorHandlingSignals:

Failure context
---------------

Handlers registered with ``signalHandlerRegister`` are called when the run
fails, and exist so that context which would otherwise be lost can be dumped ---
the posterior-sampling likelihood, for example, registers a handler which writes
the parameter set being evaluated to a file, so that the failure can be
reproduced.

They are called in two circumstances:

* when the process receives ``SIGSEGV``, ``SIGFPE``, ``SIGBUS``, ``SIGILL`` or
  ``SIGXCPU``, in which case the handler is passed the signal number; and
* from ``Error_Report``, in which case it is passed ``signalNone``.

The second matters more often than the first. A crash is rare; a deliberate
fatal error is not, and a handler called only on a signal would discard its
context in precisely the usual case.

.. warning::

   A handler which interprets the signal it is passed must account for
   ``signalNone``, which is **not** a member of the ``signal`` enumeration.
   Decoding it --- with ``enumerationSignalDecode``, say, to name a file after
   the cause --- raises an error of its own, and since that error is raised from
   *inside* the handler it replaces the error actually being reported. Test for
   ``signalNone`` first.

Handlers are ``threadprivate``, so a failure dumps the context of the thread
which failed, and that thread only --- which is what is wanted, since that is
the thread whose state is relevant. A handler must otherwise be careful to do as
little as possible: it runs in a process which has already failed, and anything
it does which itself fails will replace the original diagnosis with its own.
