#!/usr/bin/env bash

# Run the convergence validations on a cluster, and publish their results.
#
# The convergence validations (`testSuite/validate-mergerTreeConvergence.py` and `testSuite/validate-concentrationConvergence.py`)
# build merger trees at very high resolution. They need far more memory and time than a CI runner has (tens of gigabytes per
# tree, and around two weeks), so are run on a cluster instead, where they submit their models as SLURM jobs. This script does
# the whole of such a run, unattended:
#
#   1. brings the Galacticus checkout that it is in up to date with `origin/master`;
#   2. builds `Galacticus.exe`, as a SLURM job;
#   3. runs the validations, which submit their models as further SLURM jobs and wait for them;
#   4. publishes the `results.json` of each validation that succeeded to the `gh-pages` branch, where they appear on the
#      validation dashboard (see `scripts/aux/publishValidationResults.sh`).
#
# It is intended to be run each month from `scrontab`, from a checkout kept for this purpose alone (the script updates that
# checkout, and refuses to do so if it has uncommitted changes). For example:
#
#   #SCRON --job-name=galacticusConvergence
#   #SCRON --partition=obs
#   #SCRON --cpus-per-task=1
#   #SCRON --mem=32G
#   #SCRON --output=/home/me/convergence/logs/driver-%j.log
#   0 3 1 * * /home/me/convergence/galacticus/scripts/aux/convergenceValidation.sh --work /home/me/convergence --environment /home/me/convergence/environment.sh --partition obs
#
# The job that runs this script does little itself - it waits for the jobs it submits, then reads their output and makes plots -
# so needs only one core, but it must be allowed to run for as long as the validations take.
#
# Usage:
#
#   convergenceValidation.sh --work <directory> [options]
#
#   --work <directory>         where each run's output is kept, under `runs/<run name>/` (required)
#   --environment <file>       a file to source first, to set up the environment (e.g. containing `module load galacticus`)
#   --partition <partition>    the SLURM partition to submit jobs to (default: the cluster's default)
#   --job-maximum <n>          the most model jobs that each validation may have active at once (default: 10)
#   --build-cpus <n>           the number of cores for the build job (default: 24)
#   --build-memory <memory>    the memory for the build job, as SLURM writes it (default: 96G)
#   --validations "<names>"    the validations to run (default: "mergerTreeConvergence concentrationConvergence")
#   --run-name <name>          the name of this run (default: the year and month, e.g. 2026-11)
#   --branch <branch>          the branch to validate (default: master)
#   --skip-update              do not update the checkout
#   --skip-build               do not build `Galacticus.exe`
#   --skip-publish             do not publish the results
#   --dry-run-publish          do everything for publishing the results except the final push
#
# A run can be repeated under the same name: the validations do not re-run a model whose output already exists, so a run that
# failed part way through carries on from where it stopped. To continue with the same build, add `--skip-update --skip-build`.
#
# The environment must provide what a Galacticus build needs, the datasets (`GALACTICUS_DATA_PATH`), and the Python modules the
# validations use. If the checkout has no `galacticusConfig.xml`, one is written which selects SLURM on any host.
#
# The exit status is zero only if every validation succeeded and (unless skipped) its results were published. The outcome of
# each validation is also reported on a line beginning `SUCCESS:` or `FAILED:`.

set -uo pipefail

usage() {
    sed -n '3,54p' "${BASH_SOURCE[0]}" | sed -e 's/^# \{0,1\}//'
}

fail() {
    echo "FAILED: $*" >&2
    exit 1
}

# Report whether the run that holds a lock is still active. The lock records where that run is: a SLURM job, or a process.
lockIsActive() {
    local holder
    holder=$(cat "$1/holder" 2> /dev/null) || return 1
    case "${holder}" in
        slurm:*)
            [[ -n "$(squeue --noheader --jobs "${holder#slurm:}" 2> /dev/null)" ]]
            ;;
        process:*)
            local host="${holder#process:}"
            host="${host%%:*}"
            # A process on another host can not be checked from here - assume that it is still running.
            [[ "${host}" != "$(hostname)" ]] || kill -0 "${holder##*:}" 2> /dev/null
            ;;
        *)
            return 1
            ;;
    esac
}

main() {
    local pathWork="" fileEnvironment="" partition="" jobMaximum="10" buildCPUs="24" buildMemory="96G"
    local validations="mergerTreeConvergence concentrationConvergence" runName="" branch="master"
    local skipUpdate="no" skipBuild="no" skipPublish="no" dryRunPublish="no"
    local arguments=("$@")
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --work|--environment|--partition|--job-maximum|--build-cpus|--build-memory|--validations|--run-name|--branch)
                [[ $# -ge 2 ]] || fail "'$1' needs an argument"
                case "$1" in
                    --work)         pathWork="$2"        ;;
                    --environment)  fileEnvironment="$2" ;;
                    --partition)    partition="$2"       ;;
                    --job-maximum)  jobMaximum="$2"      ;;
                    --build-cpus)   buildCPUs="$2"       ;;
                    --build-memory) buildMemory="$2"     ;;
                    --validations)  validations="$2"     ;;
                    --run-name)     runName="$2"         ;;
                    --branch)       branch="$2"          ;;
                esac
                shift 2
                ;;
            --skip-update)     skipUpdate="yes"    ; shift ;;
            --skip-build)      skipBuild="yes"     ; shift ;;
            --skip-publish)    skipPublish="yes"   ; shift ;;
            --dry-run-publish) dryRunPublish="yes" ; shift ;;
            -h|--help)         usage               ; exit 0 ;;
            *)                 usage >&2           ; fail "unrecognized option '$1'" ;;
        esac
    done
    [[ -n "${pathWork}" ]] || { usage >&2; fail "'--work' is required"; }
    [[ "${jobMaximum}" =~ ^[0-9]+$ && "${buildCPUs}" =~ ^[0-9]+$ ]] || fail "'--job-maximum' and '--build-cpus' must be integers"
    [[ -n "${runName}" ]] || runName=$(date +%Y-%m)
    [[ "${runName}" =~ ^[A-Za-z0-9_.-]+$ ]] || fail "'${runName}' is not a valid run name"
    local validation
    for validation in ${validations}; do
        [[ "${validation}" =~ ^[A-Za-z0-9]+$ ]] || fail "'${validation}' is not a valid validation name"
    done

    # Set up the environment.
    if [[ -n "${fileEnvironment}" ]]; then
        [[ -r "${fileEnvironment}" ]] || fail "can not read the environment file '${fileEnvironment}'"
        # Environment files commonly refer to variables which are not set, which `set -u` would treat as an error.
        set +u
        # shellcheck disable=SC1090
        source "${fileEnvironment}" || fail "the environment file '${fileEnvironment}' failed"
        set -u
    fi
    local pathRepository pathStart="${PWD}"
    pathRepository=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
    # SLURM runs a script given directly to `sbatch` from a copy, which is not in the checkout - and so can not find it.
    [[ -e "${pathRepository}/.git" && -f "${pathRepository}/testSuite/validate-mergerTreeConvergence.py" ]] \
        || fail "this script is not in a Galacticus checkout - run it by its path there (from scrontab, or with 'sbatch --wrap')"
    export GALACTICUS_EXEC_PATH="${pathRepository}"
    export PYTHONPATH="${pathRepository}/python${PYTHONPATH:+:${PYTHONPATH}}"
    # The queue manager identifies the host from this, which a job started by `scrontab` need not have in its environment.
    export HOSTNAME="${HOSTNAME:-$(hostname)}"
    # Plots are written to files - there is no display.
    export MPLBACKEND="Agg"
    [[ -n "${GALACTICUS_DATA_PATH:-}" ]] || fail "GALACTICUS_DATA_PATH is not set - set it in the environment file"
    local command
    for command in git make python3 sbatch squeue sinfo; do
        command -v "${command}" > /dev/null 2>&1 || fail "'${command}' was not found - is the environment set up?"
    done

    mkdir -p "${pathWork}/runs" || fail "can not create '${pathWork}/runs'"
    pathWork=$(cd "${pathWork}" && pwd)
    local pathRun="${pathWork}/runs/${runName}"

    # Allow only one run at a time: a new run would update and rebuild the checkout under the feet of one still in progress.
    local pathLock="${pathWork}/lock"
    if [[ "${CONVERGENCE_VALIDATION_LOCKED:-}" != "${pathLock}" ]]; then
        if ! mkdir "${pathLock}" 2> /dev/null; then
            if lockIsActive "${pathLock}"; then
                echo "A previous run ($(cat "${pathLock}/holder")) is still active - not starting another."
                exit 0
            fi
            echo "Removing the lock of a previous run ($(cat "${pathLock}/holder" 2> /dev/null)) which is no longer active."
            rm -rf "${pathLock}"
            mkdir "${pathLock}" || fail "can not take the lock '${pathLock}'"
        fi
        if [[ -n "${SLURM_JOB_ID:-}" ]]; then
            echo "slurm:${SLURM_JOB_ID}" > "${pathLock}/holder"
        else
            echo "process:$(hostname):$$" > "${pathLock}/holder"
        fi
        export CONVERGENCE_VALIDATION_LOCKED="${pathLock}"
    fi
    # shellcheck disable=SC2064
    trap "rm -rf '${pathLock}'" EXIT

    echo "Convergence validation run '${runName}' started at $(date)"
    cd "${pathRepository}" || fail "can not enter '${pathRepository}'"

    # Bring the checkout up to date.
    if [[ "${skipUpdate}" == "no" ]]; then
        [[ -z "$(git status --porcelain --untracked-files=no)" ]] || fail "the checkout '${pathRepository}' has uncommitted changes - not updating it"
        git fetch --quiet origin "${branch}" || fail "could not fetch '${branch}'"
        local commitBefore commitAfter
        commitBefore=$(git rev-parse HEAD)
        git checkout --quiet --detach FETCH_HEAD || fail "could not check out '${branch}'"
        commitAfter=$(git rev-parse HEAD)
        if [[ "${commitBefore}" != "${commitAfter}" ]]; then
            # The update may have changed this script, which bash reads as it runs - so start again, from the new version.
            echo "Updated the checkout from ${commitBefore} to ${commitAfter} - restarting"
            trap - EXIT
            # Restart from where this run started, so that any relative paths among the arguments mean what they did.
            cd "${pathStart}" || fail "can not return to '${pathStart}'"
            exec "${pathRepository}/scripts/aux/convergenceValidation.sh" "${arguments[@]}" --skip-update
        fi
    fi
    local commit
    commit=$(git rev-parse HEAD)
    echo "Validating commit ${commit}"

    # Check that Python has the modules needed, now, and not after a build and weeks of models.
    python3 -c "import h5py, numpy, lxml, matplotlib, scipy, git, num2tex" || fail "Python is missing modules that the validations need"

    # The queue manager needs a configuration file, which is not kept in the repository. Write one, selecting SLURM on any host,
    # if there is none - the partition and job limit are given to the validations on the command line.
    if [[ ! -e galacticusConfig.xml ]]; then
        cat > galacticusConfig.xml <<'CONFIG'
<config>
  <queueManager>
    <host>
      <name>.*</name>
      <manager>slurm</manager>
    </host>
  </queueManager>
  <slurm>
    <host>
      <name>.*</name>
    </host>
  </slurm>
</config>
CONFIG
        echo "Wrote a galacticusConfig.xml selecting SLURM on any host"
    fi

    mkdir -p "${pathRun}/logs" || fail "can not create '${pathRun}/logs'"
    local optionsPartition=()
    [[ -z "${partition}" ]] || optionsPartition=(--partition "${partition}")

    # Build, as a job of its own - this job has only the one core that waiting needs.
    if [[ "${skipBuild}" == "no" ]]; then
        echo "Building Galacticus.exe (log: ${pathRun}/logs/build.log)"
        sbatch --wait --quiet --job-name=galacticusConvergenceBuild ${optionsPartition[@]+"${optionsPartition[@]}"} \
               --nodes=1 --ntasks=1 --cpus-per-task="${buildCPUs}" --mem="${buildMemory}"                         \
               --output="${pathRun}/logs/build.log" --error="${pathRun}/logs/build.log"                           \
               --wrap="make -j${buildCPUs} Galacticus.exe"                                                        \
            || fail "the build failed - see ${pathRun}/logs/build.log"
    fi
    [[ -x Galacticus.exe ]] || fail "there is no Galacticus.exe in '${pathRepository}'"

    # Run the validations, all at once: each spends its time waiting for the jobs that it submits.
    local pids=() pid
    for validation in ${validations}; do
        [[ -f "testSuite/validate-${validation}.py" ]] || fail "there is no validation 'testSuite/validate-${validation}.py'"
        echo "Running the ${validation} validation (log: ${pathRun}/logs/${validation}.log)"
        python3 "testSuite/validate-${validation}.py" --outputPath "${pathRun}/${validation}" ${optionsPartition[@]+"${optionsPartition[@]}"} \
                                                      --jobMaximum "${jobMaximum}" --waitOnSubmit 1 --waitOnActive 60                         \
                                                      > "${pathRun}/logs/${validation}.log" 2>&1 &
        pids+=("$!")
    done
    local status=0 statusValidation index=0 published=()
    for validation in ${validations}; do
        pid="${pids[${index}]}"
        index=$((index+1))
        wait "${pid}"
        statusValidation=$?
        # The validations follow the test suite's convention: they exit with status zero even on failure, which is reported by
        # a line beginning "FAILED". So a validation succeeded only if it said so, and left its results.
        if    [[ ${statusValidation} -eq 0 ]]                                                \
           && grep -q '^SUCCESS' "${pathRun}/logs/${validation}.log"                         \
           && ! grep -q -e 'FAIL' -e 'Traceback' "${pathRun}/logs/${validation}.log"         \
           && [[ -s "${pathRun}/${validation}/results.json" ]]; then
            echo "SUCCESS: the ${validation} validation completed"
            published+=(--results "${validation}=${pathRun}/${validation}/results.json")
        else
            echo "FAILED: the ${validation} validation failed (exit status ${statusValidation}) - see ${pathRun}/logs/${validation}.log"
            status=1
        fi
    done

    # Publish the results of the validations which succeeded.
    if [[ "${skipPublish}" == "no" && ${#published[@]} -gt 0 ]]; then
        local optionsPublish=()
        [[ "${dryRunPublish}" == "no" ]] || optionsPublish=(--dry-run)
        ./scripts/aux/publishValidationResults.sh "${published[@]}" ${optionsPublish[@]+"${optionsPublish[@]}"} || status=1
    fi

    echo "Convergence validation run '${runName}' of commit ${commit} finished at $(date)"
    exit ${status}
}

main "$@"
