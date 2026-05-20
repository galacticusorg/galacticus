#!/bin/sh

# Wrapper script called by Make when generating timing and memory profile reports of builds.
# Andrew Benson (21-June-2021)

shift  # get rid of the '-c' supplied by make.

# Capture the start time, run the supplied command, and capture the stop time. Where possible the
# peak memory usage of the command is also measured.
start=`date --rfc-3339=seconds`
# Peak resident set size of the task, in kilobytes. A value of -1 indicates that memory usage could
# not be measured (for example, because GNU `time` is unavailable).
maxRSS=-1
if [ -x /usr/bin/time ]; then
    # Use GNU `time` to measure the peak resident set size of the task, writing the result to a
    # temporary file so that it is not mixed into the build log.
    memFile=`mktemp`
    /usr/bin/time -f '%M' -o "$memFile" sh -c "$@"
    if [ -s "$memFile" ]; then
        read maxRSS < "$memFile"
    fi
    rm -f "$memFile"
    # Guard against any non-numeric output (e.g. an error message from `time`).
    case "$maxRSS" in
        ''|*[!0-9]*) maxRSS=-1 ;;
    esac
else
    # GNU `time` is unavailable, so run the command without measuring memory usage.
    eval "$@"
fi
stop=`date --rfc-3339=seconds`

# Write our report.
echo "++Task: {$start|$stop|$maxRSS} '$@'"
