#!/usr/bin/env python3
"""
Launch sets of Galacticus models, iterating over sets of parameters and
performing analysis on the results.  Supports launching on a variety of
platforms via launch-method hooks.

Python port of ./scripts/aux/launch.py
Andrew Benson (11-June-2010; major re-write 01-February-2014; Python port 2026)
"""

import sys
import os
import re
import copy
import hashlib
import json
import subprocess
import threading
import time
from xml.etree import ElementTree as ET

# ---------------------------------------------------------------------------
# XML utilities
# ---------------------------------------------------------------------------

# Element names that Perl XML::Simple ForceArray would always wrap in a list.
_FORCE_ARRAY = frozenset({"modify", "value", "parameter", "parameters", "requirement"})


def _xml_to_dict(element, force_array=_FORCE_ARRAY):
    """Recursively convert an XML element to nested Python dicts/lists/strings.

    Mirrors Perl XML::Simple with KeyAttr="" and the given force_array set:
    - Elements with only text content  => plain string
    - Elements with attributes/children => dict
    - Multiple same-tag siblings       => list
    - Tags in force_array              => always list (even if only one element)
    """
    result = dict(element.attrib)

    children_by_tag: dict = {}
    for child in element:
        children_by_tag.setdefault(child.tag, []).append(child)

    for tag, children in children_by_tag.items():
        converted = [_xml_to_dict(child, force_array) for child in children]
        result[tag] = converted if (tag in force_array or len(converted) > 1) else converted[0]

    text = (element.text or "").strip()
    if text:
        if result:
            result["content"] = text
        else:
            return text  # pure-text element

    return result


def _dict_to_xml_elem(tag: str, data, parent=None, _nested: bool = False):
    """Recursively build an XML element tree from a Python dict/list/string.

    Mirrors XML::Simple XMLout default behaviour:
    - Root level (_nested=False): string values  => child text elements.
    - Nested level (_nested=True): string values => XML attributes.
    - dict values always => child elements (recursed with _nested=True).
    - list values        => multiple sibling elements with the same tag.
    """
    if isinstance(data, list):
        for item in data:
            _dict_to_xml_elem(tag, item, parent, _nested)
        return None

    elem = ET.SubElement(parent, tag) if parent is not None else ET.Element(tag)

    if isinstance(data, str):
        elem.text = data
    elif isinstance(data, dict):
        for key, val in data.items():
            if key == "content":
                elem.text = str(val)
            elif isinstance(val, str):
                if _nested:
                    elem.set(key, val)
                else:
                    child = ET.SubElement(elem, key)
                    child.text = val
            elif isinstance(val, dict):
                _dict_to_xml_elem(key, val, elem, _nested=True)
            elif isinstance(val, list):
                for item in val:
                    _dict_to_xml_elem(key, item, elem, _nested=True)

    return elem


def _xml_out(data: dict, root_name: str = "parameters") -> str:
    """Serialise *data* as pretty-printed XML with *root_name* as the root tag."""
    root = _dict_to_xml_elem(root_name, data, _nested=False)
    _add_indent(root)
    return '<?xml version="1.0" standalone="yes"?>\n' + ET.tostring(root, encoding="unicode")


def _add_indent(elem, level: int = 0):
    """Add pretty-print indentation in-place (pure stdlib)."""
    pad = "\n" + "  " * level
    if len(elem):
        if not (elem.text or "").strip():
            elem.text = pad + "  "
        for child in elem:
            _add_indent(child, level + 1)
        if not (child.tail or "").strip():
            child.tail = pad
    if level and not (elem.tail or "").strip():
        elem.tail = pad


# ---------------------------------------------------------------------------
# Options / Config
# ---------------------------------------------------------------------------

def _parse_options(argv: list) -> tuple:
    """Parse ``--key value`` arguments; return (positional_args, options_dict)."""
    positional: list = []
    options: dict = {"instance": "1:1"}
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            key = argv[i][2:]
            if i + 1 >= len(argv):
                sys.exit(f"launch.py: missing value after '--{key}'")
            if argv[i + 1].startswith("--"):
                sys.exit(f"launch.py: missing value after '--{key}'")
            val = argv[i + 1]
            i += 2
            if val.startswith('"'):
                while not val.endswith('"'):
                    if i >= len(argv):
                        sys.exit(f"launch.py: unterminated quoted value for '--{key}'")
                    val += " " + argv[i]
                    i += 1
                val = val[1:-1]
            if key in options:
                existing = options[key]
                options[key] = (existing + [val] if isinstance(existing, list) else val)
            else:
                options[key] = val
        else:
            positional.append(argv[i])
            i += 1
    return positional, options


def _load_config() -> dict:
    """Load galacticusConfig.xml (mirrors Galacticus::Options::LoadConfig)."""
    exec_path = os.environ.get("GALACTICUS_EXEC_PATH", "")
    for path in (
        os.path.join(exec_path, "galacticusConfig.xml"),
        os.path.join(os.path.expanduser("~"), ".galacticusConfig.xml"),
    ):
        if os.path.isfile(path):
            return _xml_to_dict(ET.parse(path).getroot())
    return {}


# ---------------------------------------------------------------------------
# Launch script parsing
# ---------------------------------------------------------------------------

def _parse_launch_script(filename: str) -> dict:
    """Parse the launch XML file and apply defaults."""
    script = _xml_to_dict(ET.parse(filename).getroot(), force_array=_FORCE_ARRAY)
    defaults: dict = {
        "verbosity":          0,
        "md5Names":           "no",
        "useStateFile":       "no",
        "compressModels":     "no",
        "launchMethod":       "local",
        "modelRootDirectory": "./models",
        "baseParameters":     "",
        "doAnalysis":         "no",
        "splitModels":        1,
    }
    for key, val in defaults.items():
        script.setdefault(key, val)
    script["verbosity"]   = int(script["verbosity"])
    script["splitModels"] = int(script["splitModels"])
    return script


# ---------------------------------------------------------------------------
# Parameter unfolding
# ---------------------------------------------------------------------------

def _unfold_parameters(parameter_set: dict) -> list:
    """Expand array-valued parameters into a flat list of single-value dicts.

    Python port of Perl's unfoldParameters().  Two-phase algorithm:
      Phase 1 - expand arrays (iterable elements) into separate parameter sets.
      Phase 2 - move parameters with a 'parameterLevel' attribute up the hierarchy.
    """

    # -- Phase 1 : expand iterables -------------------------------------------
    #
    # Each item in params_in is {"parameter": <data>, "parent": <parent_dict_or_None>,
    #                             "name": <key_in_parent_or_None>}
    params_in: list = [{"parameter": parameter_set, "parent": None, "name": None}]
    params_intermediate: list = []

    while params_in:
        params = params_in.pop(0)
        stack  = [params]
        cloned = False

        while stack:
            node = stack.pop(0)
            p    = node["parameter"]

            if isinstance(p, list):
                # Array-valued parameter: unfold into one clone per iterable element.
                p_snapshot   = copy.deepcopy(p)
                clones_added = False
                for i, elem in enumerate(p_snapshot):
                    if isinstance(elem, dict) and elem.get("iterable") == "no":
                        # Non-iterable: push original element for child traversal.
                        stack.append({
                            "parameter": p[i],
                            "parent":    node["parent"],
                            "name":      node["name"],
                        })
                    else:
                        # Iterable: mutate parent slot then deep-clone the whole set.
                        node["parent"][node["name"]] = elem
                        params_in.append(copy.deepcopy(params))
                        clones_added = True

                if clones_added:
                    cloned = True
                    break

            elif isinstance(p, dict):
                for key, val in p.items():
                    stack.append({"parameter": val, "parent": p, "name": key})

        if not cloned:
            params_intermediate.append(params["parameter"])

    # -- Phase 2 : handle parameterLevel --------------------------------------
    params_out: list = []

    while params_intermediate:
        params   = params_intermediate.pop(0)
        modified = False

        # Build BFS-ordered list of dict nodes.
        bfs: list     = [{"parameter": params, "parent": None, "name": None}]
        ordered: list = []
        while bfs:
            node = bfs.pop(0)
            if isinstance(node["parameter"], dict):
                ordered.append(node)
                for key, val in node["parameter"].items():
                    bfs.append({"parameter": val, "parent": node, "name": key})

        # Process deepest-first (reverse BFS approximates post-order).
        for node in reversed(ordered):
            p = node["parameter"]
            if not isinstance(p, dict) or "parameterLevel" not in p:
                continue

            level_raw = p["parameterLevel"]
            try:
                level_num = float(level_raw)
                is_numeric = True
            except (TypeError, ValueError):
                is_numeric = False

            if is_numeric and level_num < 0:
                # Move parameter up the hierarchy.
                # Perl loop: while level <= 0 -> move up, increment level.
                cur_level = level_num
                target    = node
                while cur_level <= 0:
                    cur_level += 1
                    target = target["parent"]
                del p["parameterLevel"]
                del node["parent"]["parameter"][node["name"]]
                target["parameter"].pop(node["name"], None)
                target["parameter"][node["name"]] = p

            elif level_raw == "top":
                del p["parameterLevel"]
                params[node["name"]] = p
                del node["parent"]["parameter"][node["name"]]

            else:
                sys.exit(f"launch.py: unknown parameterLevel '{level_raw}'")

            modified = True
            break

        if modified:
            params_intermediate.append(params)
        else:
            params_out.append(params)

    return params_out


# ---------------------------------------------------------------------------
# Post-process helpers
# ---------------------------------------------------------------------------

def _post_analyze(job: dict, launch_script: dict):
    """Run analysis for a completed model (mirrors PostProcess::Analyze)."""
    mg = str(job["mergeGroup"])
    launch_script.setdefault("mergeGroups", {}).setdefault(mg, []).append(job["directory"])
    if len(launch_script["mergeGroups"][mg]) < launch_script["splitModels"]:
        return
    if launch_script["splitModels"] > 1:
        hdf5_files = " ".join(
            d + "/galacticus.hdf5" for d in launch_script["mergeGroups"][mg]
        )
        subprocess.run(
            f"{os.environ.get('GALACTICUS_EXEC_PATH', '')}/scripts/aux/Merge_Models.pl"
            f" {hdf5_files} {launch_script['mergeGroups'][mg][0]}/galacticusMerged.hdf5",
            shell=True,
        )
    if job.get("analysis"):
        sh = os.path.join(job["directory"], "analysis.sh")
        with open(sh, "w") as f:
            f.write(job["analysis"])
        subprocess.run(
            f"chmod u=wrx {sh}; {sh}",
            shell=True,
            stdout=open(os.path.join(job["directory"], "analysis.out"), "w"),
            stderr=subprocess.STDOUT,
        )


def _post_failed(job: dict, launch_script: dict):
    """Report a failed Galacticus run (mirrors PostProcess::Failed)."""
    import shutil as _sh
    for fname in os.listdir("."):
        if re.match(r"core\.\d+$", fname):
            _sh.move(fname, os.path.join(job["directory"], "core"))
    msg = "FAILED: A Galacticus model failed to finish:\n\n"
    if "HOSTNAME" in os.environ:
        msg += f"  Host:\t{os.environ['HOSTNAME']}\n"
    if "USER" in os.environ:
        msg += f"  User:\t{os.environ['USER']}\n\n"
    msg += f"Model output is in: {job['directory']}\n\n"
    print(msg)
    log = os.path.join(job["directory"], "galacticus.log")
    if os.path.isfile(log):
        with open(log) as f:
            print("Log follows:\n" + f.read())


def _post_cleanup(job: dict, launch_script: dict):
    """Compress model output if requested (mirrors PostProcess::CleanUp)."""
    if launch_script.get("compressModels") == "yes":
        subprocess.run(
            f"{os.environ.get('GALACTICUS_EXEC_PATH', '')}/scripts/aux/Compress_Directory.pl"
            f" {job['directory']}",
            shell=True,
        )


# ---------------------------------------------------------------------------
# Launch method hooks - local
# ---------------------------------------------------------------------------

def _local_validate(launch_script: dict):
    local = launch_script.setdefault("local", {})
    local.setdefault("threadCount", 1)
    local.setdefault("ompThreads",  "maximum")
    local.setdefault("executable",  "Galacticus.exe")
    if str(local["threadCount"]) == "maximum":
        local["threadCount"] = os.cpu_count() or 1
    if str(local["ompThreads"]) == "maximum":
        local["ompThreads"] = os.cpu_count() or 1


def _local_output_file(output_file: str, _launch_script: dict) -> str:
    return output_file


def _local_launch(jobs: list, launch_script: dict, options: dict):
    thread_count = int(launch_script["local"]["threadCount"])
    if "threadMaximum" in options:
        thread_count = int(options["threadMaximum"])

    if thread_count <= 1:
        _local_run_models(0, 1, jobs, launch_script, options)
    else:
        threads = [
            threading.Thread(
                target=_local_run_models,
                args=(i, thread_count, jobs, launch_script, options),
                daemon=True,
            )
            for i in range(thread_count)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()


def _local_run_models(i_thread: int, thread_count: int,
                      jobs: list, launch_script: dict, options: dict):
    omp       = options.get("ompThreads", launch_script["local"]["ompThreads"])
    verbose   = int(options.get("verbosity", launch_script["verbosity"]))
    exe       = options.get("executable",   launch_script["local"]["executable"])
    exec_path = os.environ.get("GALACTICUS_EXEC_PATH", "")

    for i, job in enumerate(jobs):
        if i % thread_count != i_thread:
            continue
        if verbose > 0:
            print(f" -> thread {i_thread} running job: {job['label']}")
        log = os.path.join(job["directory"], "galacticus.log")
        rc = subprocess.run(
            f"ulimit -t unlimited; ulimit -c unlimited; "
            f"export GFORTRAN_ERROR_DUMPCORE=YES; "
            f"export OMP_NUM_THREADS={omp}; "
            f"{exec_path}/{exe} {job['directory']}/parameters.xml",
            shell=True,
            stdout=open(log, "w"), stderr=subprocess.STDOUT,
        ).returncode
        if rc == 0:
            _post_analyze(job, launch_script)
        else:
            _post_failed(job, launch_script)
        _post_cleanup(job, launch_script)


# ---------------------------------------------------------------------------
# Launch method hooks - PBS
# ---------------------------------------------------------------------------

def _pbs_validate(launch_script: dict):
    pbs = launch_script.setdefault("pbs", {})
    for key, val in {
        "mpiLaunch":               "yes",
        "mpiRun":                  ("mpirun --map-by node --mca mpi_preconnect_mpi 1"
                                    " -hostfile $PBS_NODEFILE"),
        "mpiProcesses":            1,
        "maxJobsInQueue":          -1,
        "postSubmitSleepDuration": 10,
        "jobWaitSleepDuration":    60,
        "analyze":                 "yes",
    }.items():
        pbs.setdefault(key, val)
    if pbs.get("analyze") == "yes" and launch_script["splitModels"] > 1:
        sys.exit("PBS::Validate: cannot analyze models on PBS when splitting models")


def _pbs_output_file(output_file: str, launch_script: dict) -> str:
    pbs = launch_script.get("pbs", {})
    if "scratchPath" in pbs:
        return (f"{pbs['scratchPath']}/model_{launch_script['modelCounter']}"
                f"_{os.getpid()}/galacticus.hdf5")
    return output_file


def _pbs_launch(jobs: list, launch_script: dict, options: dict):
    pbs = launch_script.get("pbs", {})
    model_queue: list = []
    for job in jobs:
        script = os.path.join(job["directory"], "pbsLaunch.sh")
        _write_pbs_script(script, job, launch_script, pbs)
        model_queue.append({"script": script, "job": job})

    job_max  = int(options.get("pbsJobMaximum",       pbs.get("maxJobsInQueue",          -1)))
    post_slp = int(options.get("submitSleepDuration", pbs.get("postSubmitSleepDuration", 10)))
    wait_slp = int(options.get("waitSleepDuration",   pbs.get("jobWaitSleepDuration",    60)))

    if launch_script["verbosity"] > 0:
        print(" -> waiting for PBS jobs to finish...")

    pbs_jobs: dict = {}
    while pbs_jobs or model_queue:
        running  = _qstat_running()
        finished = False
        for jid in list(pbs_jobs):
            if jid not in running:
                print(f" -> PBS job {jid} has finished; post-processing....")
                finished = True
                _post_analyze(pbs_jobs[jid]["job"], launch_script)
                _post_cleanup(pbs_jobs[jid]["job"], launch_script)
                del pbs_jobs[jid]

        if model_queue and (len(pbs_jobs) < job_max or job_max < 0):
            item = model_queue.pop(0)
            print(f" -> submitting script: {item['script']}")
            jid = _qsub(item["script"])
            if jid:
                pbs_jobs[jid] = item
            time.sleep(post_slp)
        elif not finished:
            time.sleep(wait_slp)


def _write_pbs_script(path: str, job: dict, launch_script: dict, pbs: dict):
    exec_path = os.environ.get("GALACTICUS_EXEC_PATH", "")
    with open(path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#PBS -N Galacticus_{job['label']}\n")
        email = _config_email(launch_script)
        if email and launch_script.get("emailReport") == "yes":
            f.write(f"#PBS -M {email}\n#PBS -m bea\n")
        if "wallTime" in pbs:
            f.write(f"#PBS -l walltime={pbs['wallTime']}\n")
        if "memory" in pbs:
            f.write(f"#PBS -l mem={pbs['memory']}\n")
        if pbs.get("mpiLaunch") == "yes":
            if "mpiNodes" in pbs:
                nodes = int(pbs["mpiNodes"])
                ppn   = int(pbs["mpiProcesses"]) * int(pbs.get("ompThreads", 1)) // nodes
            else:
                nodes = 1
                ppn   = int(pbs.get("mpiProcesses", 1)) * int(pbs.get("ompThreads", 1))
        else:
            nodes, ppn = 1, int(pbs.get("ompThreads", 1))
        f.write(f"#PBS -l nodes={nodes}:ppn={ppn}\n")
        f.write("#PBS -j oe\n")
        pwd = "" if job["directory"].startswith("/") else os.getcwd() + "/"
        f.write(f"#PBS -o {pwd}{job['directory']}/galacticus.log\n")
        if "queue" in pbs:
            f.write(f"#PBS -q {pbs['queue']}\n")
        f.write("#PBS -V\n")
        f.write("if [ ! -z ${PBS_O_WORKDIR+x} ]; then\n  cd $PBS_O_WORKDIR\n"
                "elif [ ! -z ${SLURM_SUBMIT_DIR+x} ]; then\n  cd $SLURM_SUBMIT_DIR\nfi\n")
        if pbs.get("coreDump") == "yes":
            f.write("ulimit -c unlimited\nexport GFORTRAN_ERROR_DUMPCORE=YES\n")
        else:
            f.write("ulimit -c 0\nexport GFORTRAN_ERROR_DUMPCORE=NO\n")
        f.write("ulimit -t unlimited\n")
        if "ompThreads" in pbs:
            f.write(f"export OMP_NUM_THREADS={pbs['ompThreads']}\n")
        if "scratchPath" in pbs:
            sp = f"{pbs['scratchPath']}/model_{job['modelCounter']}_{os.getpid()}/"
            f.write(f"mkdir -p {sp}\n")
        for cmd in _as_list(pbs.get("preCommand")):
            f.write(cmd + "\n")
        executable = pbs.get("executable") or "Galacticus.exe"
        if pbs.get("mpiLaunch") == "yes":
            f.write(f"{pbs['mpiRun']} -np {pbs.get('mpiProcesses', 1)} ")
        f.write(f"{exec_path}/{executable} {job['directory']}/parameters.xml\n")
        if "scratchPath" in pbs:
            sb = f"{pbs['scratchPath']}/model_{job['modelCounter']}_{os.getpid()}"
            f.write(f"mv {sb}/galacticus.hdf5 {job['directory']}/galacticus.hdf5\n")
            if launch_script.get("useStateFile") == "yes":
                f.write(f"mv {sb}/galacticus_{job['modelCounter']}_{os.getpid()}.state*"
                        f" {job['directory']}/\n")
        f.write('if compgen -G "core.*" > /dev/null; then\n'
                f'  mv core* {job["directory"]}/\nfi\n')
        if job.get("analysis") and pbs.get("analyze") == "yes":
            f.write(job["analysis"])


def _qstat_running() -> set:
    try:
        out = subprocess.run(["qstat", "-f"], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return set()
    running: set = set()
    jid_cur = None
    for line in out.splitlines():
        m = re.match(r"^Job\s+Id:\s+(\S+)", line)
        if m:
            jid_cur = m.group(1)
        if re.match(r"^\s*job_state\s*=\s*[^C]\s*$", line) and jid_cur:
            running.add(jid_cur)
    return running


def _qsub(script: str) -> str:
    try:
        out = subprocess.run(["qsub", script], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return ""
    for line in out.splitlines():
        m = re.match(r"^(\d+\S+)", line)
        if m:
            return m.group(1)
    return ""


# ---------------------------------------------------------------------------
# Launch method hooks - SLURM
# ---------------------------------------------------------------------------

def _slurm_validate(launch_script: dict):
    slurm = launch_script.setdefault("slurm", {})
    for key, val in {
        "mpiLaunch":               "yes",
        "mpiRun":                  "mpirun --bynode",
        "maxJobsInQueue":          -1,
        "postSubmitSleepDuration": 10,
        "jobWaitSleepDuration":    60,
        "analyze":                 "yes",
    }.items():
        slurm.setdefault(key, val)
    if slurm.get("analyze") == "yes" and launch_script["splitModels"] > 1:
        sys.exit("SLURM::Validate: cannot analyze models on SLURM when splitting models")


def _slurm_output_file(output_file: str, launch_script: dict) -> str:
    slurm = launch_script.get("slurm", {})
    if "scratchPath" in slurm:
        return (f"{slurm['scratchPath']}/model_{launch_script['modelCounter']}"
                f"_{os.getpid()}/galacticus.hdf5")
    return output_file


def _slurm_launch(jobs: list, launch_script: dict, options: dict):
    slurm = launch_script.get("slurm", {})
    model_queue: list = []
    for job in jobs:
        script = os.path.join(job["directory"], "slurmLaunch.sh")
        _write_slurm_script(script, job, launch_script, slurm)
        model_queue.append({"script": script, "job": job})

    max_q    = int(slurm.get("maxJobsInQueue",          -1))
    post_slp = int(slurm.get("postSubmitSleepDuration", 10))
    wait_slp = int(slurm.get("jobWaitSleepDuration",    60))

    if launch_script["verbosity"] > 0:
        print(" -> waiting for SLURM jobs to finish...")

    slurm_jobs: dict = {}
    while slurm_jobs or model_queue:
        running = _squeue_running()
        for jid in list(slurm_jobs):
            if jid not in running:
                print(f" -> SLURM job {jid} has finished; post-processing....")
                _post_analyze(slurm_jobs[jid]["job"], launch_script)
                _post_cleanup(slurm_jobs[jid]["job"], launch_script)
                del slurm_jobs[jid]

        if model_queue and (len(slurm_jobs) < max_q or max_q < 0):
            item = model_queue.pop(0)
            print(f" -> submitting script: {item['script']}")
            jid = _sbatch(item["script"])
            if jid:
                print(f"    -> job ID: {jid}")
                slurm_jobs[jid] = item
            time.sleep(post_slp)
        else:
            time.sleep(wait_slp)


def _write_slurm_script(path: str, job: dict, launch_script: dict, slurm: dict):
    exec_path = os.environ.get("GALACTICUS_EXEC_PATH", "")
    with open(path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH -J Galacticus_{job['label']}\n")
        f.write(f"#SBATCH --chdir={os.getcwd()}\n")
        email = _config_email(launch_script)
        if email and launch_script.get("emailReport") == "yes":
            f.write(f"#SBATCH --mail-user={email}\n#SBATCH --mail-type=end\n")
        if "wallTime" in slurm:
            f.write(f"#SBATCH --time={slurm['wallTime']}\n")
        if "memory" in slurm:
            f.write(f"#SBATCH -l mem={slurm['memory']}\n")
        f.write("#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --cpus-per-task={slurm.get('ompThreads', 1)}\n")
        pwd = "" if job["directory"].startswith("/") else os.getcwd() + "/"
        f.write(f"#SBATCH -o {pwd}{job['directory']}/galacticus.log\n")
        if "account" in slurm:
            f.write(f"#SBATCH -A {slurm['account']}\n")
        for env in _as_list(slurm.get("environment")):
            f.write(f"export {env}\n")
        for mod in _as_list(slurm.get("module")):
            f.write(f"module load {mod}\n")
        if slurm.get("coreDump") == "yes":
            f.write("ulimit -c unlimited\nexport GFORTRAN_ERROR_DUMPCORE=YES\n")
        else:
            f.write("ulimit -c 0\nexport GFORTRAN_ERROR_DUMPCORE=NO\n")
        f.write("ulimit -t unlimited\npwd\n")
        if "ompThreads" in slurm:
            f.write(f"export OMP_NUM_THREADS={slurm['ompThreads']}\n")
        if "scratchPath" in slurm:
            sp = f"{slurm['scratchPath']}/model_{job['modelCounter']}_{os.getpid()}/"
            f.write(f"mkdir -p {sp}\n")
        for cmd in _as_list(slurm.get("preCommand")):
            f.write(cmd + "\n")
        if slurm.get("mpiLaunch") == "yes":
            f.write(f"{slurm['mpiRun']} -np 1 ")
        f.write(f"{exec_path}/Galacticus.exe {job['directory']}/parameters.xml\n")
        if "scratchPath" in slurm:
            sb = f"{slurm['scratchPath']}/model_{job['modelCounter']}_{os.getpid()}"
            f.write(f"mv {sb}/galacticus.hdf5 {job['directory']}/galacticus.hdf5\n")
            if launch_script.get("useStateFile") == "yes":
                f.write(f"mv {sb}/galacticus_{job['modelCounter']}_{os.getpid()}.state*"
                        f" {job['directory']}/\n")
        f.write(f"mv core* {job['directory']}/\n")
        if job.get("analysis") and slurm.get("analyze") == "yes":
            f.write(job["analysis"])


def _squeue_running() -> set:
    try:
        out = subprocess.run(["squeue"], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return set()
    running: set = set()
    for line in out.splitlines():
        m = re.match(r"^\s+(\S+)", line)
        if m:
            running.add(m.group(1))
    return running


def _sbatch(script: str) -> str:
    try:
        out = subprocess.run(["sbatch", script], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return ""
    for line in out.splitlines():
        m = re.match(r"^Submitted batch job (\d+)", line)
        if m:
            return m.group(1)
    return ""


# ---------------------------------------------------------------------------
# Launch method hooks - MonolithicPBS
# ---------------------------------------------------------------------------

def _mono_pbs_validate(launch_script: dict):
    mono = launch_script.setdefault("monolithicPBS", {})
    ncpu = os.cpu_count() or 1
    for key, val in {
        "mpiRun":               "mpirun",
        "mpiOptions":           "",
        "nodes":                1,
        "threadsPerNode":       ncpu,
        "ompThreads":           ncpu,
        "shell":                "bash",
        "analyze":              "yes",
        "jobWaitSleepDuration": 60,
    }.items():
        mono.setdefault(key, val)
    nodes = int(mono["nodes"])
    tpn   = int(mono["threadsPerNode"])
    omp   = int(mono["ompThreads"])
    if tpn % omp != 0:
        sys.exit("MonolithicPBS::Validate: ompThreads must be a factor of threadsPerNode")
    mono["mpiThreads"] = nodes * tpn // omp

    model_root = launch_script["modelRootDirectory"]
    os.makedirs(model_root, exist_ok=True)
    c_src = os.path.join(model_root, "mpiRank.c")
    c_exe = os.path.join(model_root, "mpiRank")
    with open(c_src, "w") as f:
        f.write("#include <stdio.h>\n#include <mpi.h>\n"
                "int main(int argc, char *argv[]) {\n"
                "  int rank, size; char command[1024];\n"
                "  if (argc != 2) { printf(\"usage: %s <executable>\\n\", argv[0]); return 1; }\n"
                "  MPI_Init(&argc, &argv);\n"
                "  MPI_Comm_rank(MPI_COMM_WORLD, &rank);\n"
                "  MPI_Comm_size(MPI_COMM_WORLD, &size);\n"
                "  sprintf(command, \"%s %d %d\", argv[1], rank, size);\n"
                "  printf(\"mpiRank running: %s\\n\", command);\n"
                "  system(command);\n"
                "  MPI_Finalize();\n"
                "  return 0;\n}\n")
    cmd = f"mpicc {c_src} -o {c_exe}"
    for p in _as_list(mono.get("includePath")):
        cmd += f" -I{p}"
    for p in _as_list(mono.get("libraryPath")):
        cmd += f" -L{p}"
    if subprocess.run(cmd, shell=True).returncode != 0:
        sys.exit("MonolithicPBS::Validate: failed to compile mpiRank helper")


def _mono_pbs_output_file(output_file: str, launch_script: dict) -> str:
    mono = launch_script.get("monolithicPBS", {})
    if "scratchPath" in mono:
        return (f"{mono['scratchPath']}/model_{launch_script['modelCounter']}"
                f"_{os.getpid()}/galacticus.hdf5")
    return output_file


def _mono_pbs_launch(jobs: list, launch_script: dict, _options: dict):
    mono        = launch_script.get("monolithicPBS", {})
    model_root  = launch_script["modelRootDirectory"]
    shell       = mono.get("shell", "bash")
    exec_path   = os.environ.get("GALACTICUS_EXEC_PATH", "")
    mpi_threads = int(mono["mpiThreads"])

    single_rank = -1
    single_lines = [f"#!/usr/bin/env {shell}\n", "MPI_RANK=$1\n"]

    for job in jobs:
        launch_sh = os.path.join(job["directory"], "launch.sh")
        with open(launch_sh, "w") as f:
            f.write(f"#!/usr/bin/env {shell}\n")
            f.write(_mono_set_env(mono, shell))
            f.write("export GFORTRAN_ERROR_DUMPCORE=YES\n"
                    "ulimit -t unlimited\nulimit -c unlimited\n")
            f.write(f"export OMP_NUM_THREADS={mono['ompThreads']}\n")
            if "scratchPath" in mono:
                sp = (f"{mono['scratchPath']}/model_{launch_script['modelCounter']}"
                      f"_{os.getpid()}/")
                f.write(f"mkdir -p {sp}\n")
            f.write(f"{exec_path}/Galacticus.exe {job['directory']}/parameters.xml\n")
            if "scratchPath" in mono:
                f.write(f"mv {mono['scratchPath']}galacticus.hdf5"
                        f" {job['directory']}/galacticus.hdf5\n")
                if launch_script.get("useStateFile") == "yes":
                    mc = launch_script['modelCounter']
                    f.write(f"mv {mono['scratchPath']}galacticus_{mc}_{os.getpid()}.state*"
                            f" {job['directory']}/\n")
            f.write(f"mv core* {job['directory']}/\n")
            if job.get("analysis") and mono.get("analyze") == "yes":
                f.write(job["analysis"])
        single_rank = (single_rank + 1) % mpi_threads
        single_lines += [f"if [ $MPI_RANK -eq {single_rank} ]; then\n",
                         f"  source {launch_sh}\n", "fi\n"]
    single_lines.append("exit\n")

    launch_galacticus = os.path.join(model_root, "launchGalacticus.sh")
    with open(launch_galacticus, "w") as f:
        f.writelines(single_lines)

    res_req = mono.get("resourceRequest",
                       f"#PBS -l nodes={mono['nodes']}:ppn={mono['threadsPerNode']}")
    res_req = res_req.replace("%%NODES%%",   str(mono["nodes"]))
    res_req = res_req.replace("%%THREADS%%", str(mono["threadsPerNode"]))

    pbs_lines = [f"#!/usr/bin/env {shell}\n", "#PBS -N Galacticus\n"]
    email = _config_email(launch_script)
    if email and launch_script.get("emailReport") == "yes":
        pbs_lines += [f"#PBS -M {email}\n", "#PBS -m bea\n"]
    if "wallTime" in mono:
        pbs_lines.append(f"#PBS -l walltime={mono['wallTime']}\n")
    if "memory" in mono:
        pbs_lines.append(f"#PBS -l mem={mono['memory']}\n")
    pbs_lines += [res_req + "\n", "#PBS -j oe\n",
                  f"#PBS -o {model_root}/galacticus.log\n"]
    if "queue" in mono:
        pbs_lines.append(f"#PBS -q {mono['queue']}\n")
    pbs_lines += ["#PBS -V\n", "cd $PBS_O_WORKDIR\n",
                  f"chmod u=wrx {launch_galacticus}\n"]
    for cmd in _as_list(mono.get("pbsCommand")):
        pbs_lines.append(cmd + "\n")
    pbs_lines += ["setenv MPI_DSM_DISTRIBUTE 0\n",
                  "setenv KMP_AFFINITY disabled\n",
                  f"setenv OMP_NUM_THREADS {mono['ompThreads']}\n",
                  f"{mono['mpiRun']} -np {mpi_threads} {mono['mpiOptions']}"
                  f" {model_root}/mpiRank {launch_galacticus}\n"]

    pbs_file = os.path.join(model_root, "pbsLaunch.sh")
    with open(pbs_file, "w") as f:
        f.writelines(pbs_lines)

    print(f" -> submitting script: {pbs_file}")
    job_id = _qsub(pbs_file)
    if not job_id:
        sys.exit("MonolithicPBS: failed to capture job ID")

    if launch_script["verbosity"] > 0:
        print(" -> waiting for PBS jobs to finish...")

    wait_slp = int(mono.get("jobWaitSleepDuration", 60))
    while True:
        if job_id not in _qstat_running():
            break
        time.sleep(wait_slp)

    for job in jobs:
        if mono.get("analyze") != "yes":
            _post_analyze(job, launch_script)
        _post_cleanup(job, launch_script)


def _mono_set_env(mono: dict, shell: str) -> str:
    envs = _as_list(mono.get("environment"))
    if not envs:
        return ""
    if shell == "bash":
        return "".join(f"export {e}\n" for e in envs)
    if shell == "csh":
        parts = []
        for e in envs:
            k, _, v = e.partition("=")
            parts.append(f"setenv {k} {v}\n" if v else f"setenv {e}\n")
        return "".join(parts)
    return ""


# ---------------------------------------------------------------------------
# Hook registry
# ---------------------------------------------------------------------------

_MODULE_HOOKS: dict = {
    "local": {
        "validate":       _local_validate,
        "outputFileName": _local_output_file,
        "launch":         _local_launch,
    },
    "pbs": {
        "validate":       _pbs_validate,
        "outputFileName": _pbs_output_file,
        "launch":         _pbs_launch,
    },
    "slurm": {
        "validate":       _slurm_validate,
        "outputFileName": _slurm_output_file,
        "launch":         _slurm_launch,
    },
    "monolithicPBS": {
        "validate":       _mono_pbs_validate,
        "outputFileName": _mono_pbs_output_file,
        "launch":         _mono_pbs_launch,
    },
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _as_list(val) -> list:
    """Return *val* wrapped in a list if it is not already one."""
    if val is None:
        return []
    return val if isinstance(val, list) else [val]


def _config_email(launch_script: dict) -> str:
    try:
        return launch_script["config"]["contact"]["email"]
    except (KeyError, TypeError):
        return ""


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def _construct_models(launch_script: dict) -> list:
    """Build model directories and parameter files; return list of job dicts."""
    random_seed = 219
    launch_script["modelCounter"] = -1
    jobs: list = []
    method = launch_script["launchMethod"]

    parameters_list = launch_script.get("parameters", [])
    if isinstance(parameters_list, dict):
        parameters_list = [parameters_list]

    for i_model_set, parameter_set in enumerate(parameters_list):
        # Base name for this parameter set (from optional <label> element).
        raw_label  = parameter_set.get("label")
        model_base = (raw_label if isinstance(raw_label, str) and raw_label
                      else "galacticus")

        unfolded    = _unfold_parameters(parameter_set)
        i_model     = 0
        merge_group = 0

        for parameter_data in unfolded:
            merge_group += 1

            for worker in range(1, launch_script["splitModels"] + 1):
                i_model += 1
                launch_script["modelCounter"] += 1
                random_seed += 1

                run_on = (launch_script["modelCounter"] % launch_script["instanceCount"]) + 1

                if launch_script["splitModels"] > 1:
                    parameter_data["treeEvolveWorkerNumber"] = {"value": str(worker)}
                    parameter_data["treeEvolveWorkerCount"]  = {
                        "value": str(launch_script["splitModels"])
                    }

                # Determine model label.
                lbl = parameter_data.get("label")
                if isinstance(lbl, dict):
                    model_label = lbl.get("value", f"{i_model_set}:{i_model}")
                elif isinstance(lbl, str) and lbl:
                    model_label = lbl
                else:
                    model_label = f"{i_model_set}:{i_model}"

                output_dir = (f"{launch_script['modelRootDirectory']}"
                              f"/{model_base}_{model_label}")

                descriptor = None
                if launch_script["md5Names"] == "yes":
                    descriptor = json.dumps(parameter_data, sort_keys=True, default=str)
                    md5_hash   = hashlib.md5(descriptor.encode()).hexdigest()
                    output_dir = (f"{launch_script['modelRootDirectory']}"
                                  f"/{model_base}_{md5_hash}")

                # Skip if directory already exists or this is the wrong instance.
                if os.path.exists(output_dir) or run_on != launch_script["thisInstance"]:
                    continue

                os.makedirs(output_dir, exist_ok=True)

                if descriptor is not None:
                    with open(os.path.join(output_dir, "md5Descriptor.txt"), "w") as fh:
                        fh.write(descriptor + "\n")

                output_hdf5 = os.path.join(output_dir, "galacticus.hdf5")

                # Load optional base parameters.
                parameters: dict = {}
                base_file = launch_script.get("baseParameters", "")
                if base_file:
                    parameters = _xml_to_dict(
                        ET.parse(base_file).getroot(),
                        force_array=frozenset(),
                    )

                # Set output file name via launch-method hook.
                parameters["outputFileName"] = {
                    "value": _MODULE_HOOKS[method]["outputFileName"](output_hdf5, launch_script)
                }

                # Set random seed if not already in base parameters.
                if "randomNumberGenerator" not in parameters:
                    parameters["randomNumberGenerator"] = {
                        "value": "GSL",
                        "seed":  {"value": str(random_seed)},
                    }

                # State restore file.
                if launch_script.get("useStateFile") == "yes":
                    state_root = re.sub(r"\.hdf5$", "",
                                        parameters["outputFileName"]["value"])
                    parameters["stateFileRoot"] = {"value": state_root}

                # Merge model-specific parameters (overrides base).
                for key, val in parameter_data.items():
                    parameters[key] = val

                # Write parameters.xml.
                with open(os.path.join(output_dir, "parameters.xml"), "w") as fh:
                    fh.write(_xml_out(parameters, root_name="parameters"))

                # Analysis command.
                analysis = None
                if launch_script.get("doAnalysis") == "yes":
                    analysis = f"{launch_script['analysisScript']} {output_dir}\n"

                jobs.append({
                    "label":        model_label,
                    "directory":    output_dir,
                    "analysis":     analysis,
                    "mergeGroup":   merge_group,
                    "modelCounter": launch_script["modelCounter"],
                })

    return jobs


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: launch.py <runFile>")

    positional, arguments = _parse_options(sys.argv[1:])
    if not positional:
        sys.exit("Usage: launch.py <runFile>")

    launch_script = _parse_launch_script(positional[0])
    launch_script["config"] = _load_config()

    if "launchMethod" in arguments:
        launch_script["launchMethod"] = arguments["launchMethod"]

    m = re.match(r"(\d+):(\d+)", arguments.get("instance", "1:1"))
    if not m:
        sys.exit("launch.py: 'instance' argument syntax error")
    launch_script["thisInstance"]  = int(m.group(1))
    launch_script["instanceCount"] = int(m.group(2))
    if launch_script["verbosity"] > 0:
        print(f" -> launching instance {launch_script['thisInstance']}"
              f" of {launch_script['instanceCount']}")

    method = launch_script["launchMethod"]
    if method not in _MODULE_HOOKS:
        sys.exit(f"launch.py: unrecognized launch method '{method}'")

    _MODULE_HOOKS[method]["validate"](launch_script)
    jobs = _construct_models(launch_script)
    _MODULE_HOOKS[method]["launch"](jobs, launch_script, arguments)


if __name__ == "__main__":
    main()
