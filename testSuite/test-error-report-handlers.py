#!/usr/bin/env python3
"""Test that registered signal handlers are called by `Error_Report` (issue #1353).

Handlers registered with `signalHandlerRegister` exist so that a failure can dump context which
would otherwise be lost - the posterior-sampling likelihood registers one which writes the parameter
set being evaluated to a file. They were called only when a signal was caught, so a *deliberate*
fatal error, which is far the more common case, discarded that context.

`tests.error_report_handlers.exe` registers a handler and then reports a fatal error. It necessarily
aborts, so the handler records what it saw in a file, which is checked here: the handler must have
run at all, and must have been passed `signalNone` rather than a signal number, since a handler
which decodes the signal cannot decode a value which is not one.

Andrew Benson
"""

import os
import subprocess
import sys

executable = "tests.error_report_handlers.exe"
record     = "outputs/errorReportHandlers.record"

subprocess.run("mkdir -p outputs",shell=True)
# Remove any record from an earlier run, so that its presence below is meaningful.
if os.path.exists(record):
    os.remove(record)

failed = False

# The program aborts by design, so a non-zero return code is expected; it is the handler's record,
# not the return code, which is under test.
process = subprocess.run(
    f"cd ..; ./{executable}",
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
)
output = process.stdout or ""
with open("outputs/test-error-report-handlers.log","w") as logFile:
    logFile.write(output)

if process.returncode == 0:
    print("FAILED: the test program was expected to abort, but exited successfully")
    failed = True

if not failed and not os.path.exists(record):
    print("FAILED: the registered handler was not called by `Error_Report` - it wrote no record")
    failed = True

if not failed:
    values = {}
    with open(record) as recordFile:
        for line in recordFile:
            key,_,value = line.strip().partition("=")
            values[key] = value
    if values.get("isSignalNone") != "T":
        print(f"FAILED: the handler was called with signal={values.get('signal')}, but should have been passed `signalNone` when called other than in response to a signal")
        failed = True
    else:
        print("SUCCESS: the registered handler was called by `Error_Report`, and passed `signalNone`")

# The error being reported must still reach the user: the handler runs alongside it, not instead of
# it.
if not failed and "deliberate error, to test that registered handlers are called" not in output:
    print("FAILED: the error message was not reported")
    failed = True

# Failure is signaled by the text above, not by the exit status.
if failed:
    print("FAILED: error report handlers")
else:
    print("SUCCESS: error report handlers")
sys.exit(0)
