#!/bin/sh

# Wrapper script call by Make when generating timing profile reports of builds.
# Andrew Benson (21-June-2021)

shift  # get rid of the '-c' supplied by make.

# Capture start and stop times around the supplied command.
start=`date --rfc-3339=seconds`
eval "$@"
stop=`date --rfc-3339=seconds`

# Write our report.
echo "++Task: {$start|$stop} '$@'"
