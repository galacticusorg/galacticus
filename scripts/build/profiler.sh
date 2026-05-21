#!/bin/sh

# Wrapper script called by Make when generating timing and memory profile reports of builds.
# Andrew Benson (21-June-2021)

shift  # get rid of the '-c' supplied by make.

# Determine whether GNU `time` is available. The presence of an executable `/usr/bin/time` is not
# sufficient: on some systems (e.g. macOS) this is BSD `time`, which does not support the `-f` and
# `-o` options used below. We therefore probe explicitly for support of the GNU syntax.
gnuTime=
if [ -x /usr/bin/time ] && /usr/bin/time -f '%M' -o /dev/null true >/dev/null 2>&1; then
    gnuTime=/usr/bin/time
fi

# Ensure the temporary file (if any) is removed when the script exits, including if the task is
# interrupted by a signal.
memFile=
trap 'test -n "$memFile" && rm -f "$memFile"'           EXIT
trap 'test -n "$memFile" && rm -f "$memFile"; exit 1'   HUP INT TERM

# Capture the start time, run the supplied command, and capture the stop time. Where possible the
# peak memory usage of the command is also measured.
start=`date --rfc-3339=seconds`
# Peak resident set size of the task, in kilobytes. A value of -1 indicates that memory usage could
# not be measured (for example, because GNU `time` is unavailable).
maxRSS=-1
if [ -n "$gnuTime" ] && memFile=`mktemp 2>/dev/null`; then
    # Use GNU `time` to measure the peak resident set size of the task, writing the result to the
    # temporary file so that it is not mixed into the build log.
    "$gnuTime" -f '%M' -o "$memFile" sh -c "$@"
    if [ -s "$memFile" ]; then
        read maxRSS < "$memFile"
    fi
    # Guard against any non-numeric output (e.g. an error message from `time`).
    case "$maxRSS" in
        ''|*[!0-9]*) maxRSS=-1 ;;
    esac
else
    # GNU `time` is unavailable (or a temporary file could not be created), so run the command
    # without measuring memory usage.
    memFile=
    eval "$@"
fi
stop=`date --rfc-3339=seconds`

# Write our report.
echo "++Task: {$start|$stop|$maxRSS} '$@'"
