#!/usr/bin/env bash

# Publish validation results to the `gh-pages` branch, and regenerate the status dashboard there.
#
# The CI/CD workflow does this for the validations that it runs (see the `Validate` job in `.github/workflows/cicd.yml`). This
# script does the same for validations run elsewhere - such as the convergence validations, which are too large for a CI runner
# and are run on a cluster (see `scripts/aux/convergenceValidation.sh`).
#
# Usage:
#
#   publishValidationResults.sh --results <name>=<file> [--results <name>=<file> ...] [--remote <url>] [--dry-run]
#
# where `<name>` is the name of a validation in `scripts/build/ghPagesMetricsManifest.yml` (it is also the name of its directory
# under `dev/valid/` on the `gh-pages` branch), and `<file>` is the `results.json` file which that validation wrote. For example:
#
#   publishValidationResults.sh --results mergerTreeConvergence=/scratch/mergerTreeConvergence/results.json
#
# With `--dry-run` everything is done except the push, and a summary of what would have been pushed is shown.
#
# Pushing needs write access to the repository. If the GitHub CLI (`gh`) is installed and logged in, its credentials are used;
# otherwise `git` is left to find credentials as it normally does.
#
# The exit status is zero if the results were published (or were already up to date), and non-zero otherwise.

set -euo pipefail

usage() {
    sed -n '3,23p' "${BASH_SOURCE[0]}" | sed -e 's/^# \{0,1\}//'
}

main() {
    local remote="https://github.com/galacticusorg/galacticus.git"
    local branch="gh-pages"
    local dryRun="no"
    local results=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --results)
                [[ $# -ge 2 ]] || { echo "FAILED: '--results' needs an argument" >&2; exit 1; }
                results+=("$2")
                shift 2
                ;;
            --remote)
                [[ $# -ge 2 ]] || { echo "FAILED: '--remote' needs an argument" >&2; exit 1; }
                remote="$2"
                shift 2
                ;;
            --dry-run)
                dryRun="yes"
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                echo "FAILED: unrecognized option '$1'" >&2
                usage >&2
                exit 1
                ;;
        esac
    done
    [[ ${#results[@]} -gt 0 ]] || { echo "FAILED: no results were given" >&2; usage >&2; exit 1; }

    # The scripts and manifest used to regenerate the dashboard are those of the checkout this script is in.
    local pathRepository
    pathRepository=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
    local manifest="${pathRepository}/scripts/build/ghPagesMetricsManifest.yml"
    local thresholds="${pathRepository}/scripts/build/ghPagesAnalysisThresholds.yml"
    local indexer="${pathRepository}/scripts/build/buildGhPagesIndex.py"

    # Check each result before touching anything.
    local result name file names=()
    for result in "${results[@]}"; do
        name="${result%%=*}"
        file="${result#*=}"
        [[ "${result}" == *=* && -n "${name}" && -n "${file}" ]] || { echo "FAILED: '${result}' is not of the form <name>=<file>" >&2; exit 1; }
        [[ "${name}" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "FAILED: '${name}' is not a valid validation name" >&2; exit 1; }
        # A name which the manifest does not list would be published, but would never appear on the dashboard.
        grep -q "^  ${name}:" "${manifest}" || { echo "FAILED: validation '${name}' is not listed in ${manifest}" >&2; exit 1; }
        [[ -s "${file}" ]] || { echo "FAILED: results file '${file}' does not exist, or is empty" >&2; exit 1; }
        # Results files are JavaScript assigning the results to a variable, which is what the pages on `gh-pages` load.
        head -c 200 "${file}" | grep -q -E '^window\.[A-Z_]+ = ' || { echo "FAILED: '${file}' does not look like a validation results file" >&2; exit 1; }
        names+=("${name}")
    done

    # Use the GitHub CLI's credentials if it is installed and logged in.
    local gitAuthentication=()
    if command -v gh > /dev/null 2>&1 && gh auth status > /dev/null 2>&1; then
        gitAuthentication=(-c credential.helper= -c "credential.helper=!gh auth git-credential")
    fi

    local pathWork
    pathWork=$(mktemp -d)
    # shellcheck disable=SC2064
    trap "rm -rf '${pathWork}'" EXIT

    # The CI/CD workflow also pushes to `gh-pages`, so a push here can be rejected because the branch moved after it was cloned.
    # Each attempt therefore starts afresh from the branch as it then is.
    local attempt pathPages path paths
    for attempt in 1 2 3; do
        pathPages="${pathWork}/pages${attempt}"
        git ${gitAuthentication[@]+"${gitAuthentication[@]}"} clone --quiet --depth 1 --single-branch --branch "${branch}" "${remote}" "${pathPages}"
        for result in "${results[@]}"; do
            name="${result%%=*}"
            file="${result#*=}"
            mkdir -p "${pathPages}/dev/valid/${name}"
            cp "${file}" "${pathPages}/dev/valid/${name}/results.json"
        done
        # If the branch already has these results there is nothing to do. (Checked before regenerating the dashboard, which
        # records the time at which it was generated and so always differs.)
        if [[ -z "$(git -C "${pathPages}" status --porcelain -- dev/valid)" ]]; then
            echo "SUCCESS: results are already up to date on ${branch} (${names[*]})"
            exit 0
        fi
        # Regenerate the dashboard, so that it shows these results at once. If that is not possible here the results are still
        # published, and the dashboard catches up when the CI/CD workflow next regenerates it.
        if python3 -c "import yaml" > /dev/null 2>&1; then
            python3 "${indexer}" --manifest "${manifest}" --thresholds "${thresholds}" "${pathPages}"
        else
            echo "Python has no 'yaml' module here - not regenerating the dashboard (the CI/CD workflow will do so when it next runs)"
        fi
        paths=()
        for path in _config.yml Gemfile index.md dashboard.md _data groups dev/valid; do
            [[ -e "${pathPages}/${path}" ]] && paths+=("${path}")
        done
        git -C "${pathPages}" add -A -- "${paths[@]}"
        # Commit as whoever is running this, if git knows who that is.
        git -C "${pathPages}" config user.name  > /dev/null 2>&1 || git -C "${pathPages}" config user.name  "galacticus-validation"
        git -C "${pathPages}" config user.email > /dev/null 2>&1 || git -C "${pathPages}" config user.email "github@users.noreply.github.com"
        git -C "${pathPages}" commit --quiet -m "Update validation results (${names[*]})"
        if [[ "${dryRun}" == "yes" ]]; then
            git -C "${pathPages}" show --stat --format='Would push to %d:%n' HEAD
            echo "SUCCESS: dry run - nothing was pushed (${names[*]})"
            exit 0
        fi
        if git -C "${pathPages}" ${gitAuthentication[@]+"${gitAuthentication[@]}"} push --quiet origin "${branch}"; then
            echo "SUCCESS: published validation results to ${branch} (${names[*]})"
            exit 0
        fi
        echo "Push was rejected (attempt ${attempt} of 3) - trying again from the current ${branch}"
        sleep $((attempt*20))
    done
    echo "FAILED: could not push validation results to ${branch} (${names[*]})" >&2
    exit 1
}

main "$@"
