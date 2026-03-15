#!/usr/bin/env python3
import sys
import os
import stat
import subprocess
import threading

# A simple wrapper script which launches multiple threads to process the input list of scripts.
# Andrew Benson (31-May-2016)

tasks = list(sys.argv[1:])

# Validate that an even number of arguments (script/log pairs) has been provided.
if len(tasks) % 2 != 0:
    prog = os.path.basename(sys.argv[0]) if sys.argv else "localLaunchWrapper.py"
    usage = (
        f"Usage: {prog} SCRIPT_FILE_1 LOG_FILE_1 "
        "[SCRIPT_FILE_2 LOG_FILE_2 ...]\n"
    )
    sys.stderr.write(usage)
    sys.exit(1)

threads = []


def run_script(script_file, log_file):
    # Ensure the script is executable (chmod u+x).
    current_mode = os.stat(script_file).st_mode
    os.chmod(script_file, current_mode | stat.S_IXUSR)
    with open(log_file, 'w') as log:
        # Execute the script directly without invoking a shell to avoid shell metacharacter issues.
        subprocess.run([script_file], stdout=log, stderr=subprocess.STDOUT)


while tasks:
    script_file = tasks.pop(0)
    log_file    = tasks.pop(0)
    t = threading.Thread(target=run_script, args=(script_file, log_file))
    threads.append(t)
    t.start()

for t in threads:
    t.join()
