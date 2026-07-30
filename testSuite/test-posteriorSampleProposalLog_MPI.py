#!/usr/bin/env python3
"""Test the proposed-state log and the step labelling of sampled model vectors.

Differential evolution needs at least two chains, hence MPI.

Two behaviours are checked, both concerning whether output recorded *during* a
likelihood evaluation can be joined back to the parameters that produced it:

1. ``[logProposals]`` writes every evaluated proposal, accepted or not.  The
   chain log records only the state retained at each step, so this is the only
   place a rejected proposal's parameters appear.
2. Model vectors written via ``[pathSamples]`` are tagged with the simulation
   step, which must match the step under which the chain log records that same
   evaluation.  These previously differed by one, silently mispairing every
   record.

The checks are separated from the model run so they can be exercised against
synthetic logs without an MPI build.

Andrew Benson
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np

CHAINS = "outputs/posteriorSampleProposalLogChains"
SAMPLES = "outputs/posteriorSampleProposalLogSamples"


def read_table(path, boolean_columns=()):
    """Read a whitespace-delimited log, mapping T/F columns to 1/0."""
    rows = []
    with open(path) as fh:
        for line in fh:
            stripped = line.lstrip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            for column in boolean_columns:
                tokens[column] = "1" if tokens[column] in ("T", ".true.") else "0"
            rows.append([float(t) for t in tokens])
    return np.array(rows) if rows else np.empty((0, 0))


def check(processes, chains=CHAINS, samples=SAMPLES):
    """Return a list of failure messages; empty means the test passed."""
    failures = []
    warnings = []

    chain_files = sorted(glob.glob(f"{chains}_[0-9][0-9][0-9][0-9].log"))
    proposal_files = sorted(glob.glob(f"{chains}Proposals_[0-9][0-9][0-9][0-9].log"))
    if len(chain_files) != processes:
        failures.append(f"expected {processes} chain logs, found {len(chain_files)}")
    if len(proposal_files) != processes:
        failures.append(
            f"expected {processes} proposal logs, found {len(proposal_files)} "
            f"- is [logProposals] honored?"
        )
    if failures:
        return failures, warnings

    for rank in range(processes):
        chain = read_table(f"{chains}_{rank:04d}.log", boolean_columns=[3])
        proposals = read_table(
            f"{chains}Proposals_{rank:04d}.log", boolean_columns=[2]
        )
        if chain.size == 0 or proposals.size == 0:
            failures.append(f"rank {rank}: chain or proposal log is empty")
            continue

        # chain     = step, chain, time, converged, logPost, logLike, state...
        # proposals = step, chain, accepted, logPost, logLike, state...
        chain_step = chain[:, 0].astype(int)
        chain_state = chain[:, 6:]
        proposal_step = proposals[:, 0].astype(int)
        proposal_accepted = proposals[:, 2].astype(bool)
        proposal_state = proposals[:, 5:]

        if proposal_state.shape[1] != chain_state.shape[1]:
            failures.append(
                f"rank {rank}: proposal log has {proposal_state.shape[1]} parameter "
                f"columns, chain log has {chain_state.shape[1]}"
            )
            continue
        if not np.all(np.isin(proposal_step, chain_step)):
            failures.append(
                f"rank {rank}: proposal log contains steps absent from the chain log"
            )
            continue
        if len(set(proposal_step.tolist())) != proposal_step.size:
            failures.append(f"rank {rank}: proposal log contains duplicate steps")

        # A step is an acceptance exactly when the chain state changed there.  The
        # first row has no predecessor and so cannot be classified.
        changed = np.zeros(chain_step.size, bool)
        changed[1:] = np.any(chain_state[1:] != chain_state[:-1], axis=1)
        accepted_steps = set(chain_step[changed].tolist())
        index = {s: i for i, s in enumerate(chain_step)}

        classifiable = proposal_step != chain_step[0]
        n_rejected = 0
        rejected_differs = False
        for step, accepted, state in zip(
            proposal_step[classifiable],
            proposal_accepted[classifiable],
            proposal_state[classifiable],
        ):
            if accepted != (step in accepted_steps):
                failures.append(
                    f"rank {rank}, step {step}: proposal log says accepted="
                    f"{accepted}, chain log says {step in accepted_steps}"
                )
                break
            retained = chain_state[index[step]]
            if accepted:
                # An accepted proposal is the state the chain then retained.
                if not np.allclose(state, retained, rtol=1e-12):
                    failures.append(
                        f"rank {rank}, step {step}: accepted proposal does not "
                        f"match the state recorded in the chain log"
                    )
                    break
            else:
                n_rejected += 1
                if not np.allclose(state, retained, rtol=1e-12):
                    rejected_differs = True

        # Rejected proposals must differ from the retained state; that is the
        # point of logging them, and it confirms the proposed rather than the
        # retained state is written.
        if n_rejected == 0:
            warnings.append(
                f"rank {rank} had no rejected proposals; acceptance test is weak"
            )
        elif not rejected_differs:
            failures.append(
                f"rank {rank}: every rejected proposal equals the retained state, "
                f"suggesting the retained rather than proposed state is logged"
            )

        # Sampled model vectors must carry the same step as the chain log.
        sample_files = sorted(glob.glob(f"{samples}/*_{rank:04d}.txt"))
        if not sample_files:
            failures.append(f"rank {rank}: no sampled model vectors found in {samples}")
            continue
        sample_step = read_table(sample_files[0])[:, 0].astype(int)

        if not np.all(np.isin(sample_step, chain_step)):
            failures.append(
                f"rank {rank}: sampled model vector steps are absent from the chain "
                f"log - the step labelling of [pathSamples] output is offset"
            )
            continue
        # Every acceptance was necessarily evaluated, so must have been recorded.
        missing = accepted_steps - set(sample_step.tolist())
        if missing:
            failures.append(
                f"rank {rank}: {len(missing)} accepted steps have no sampled model "
                f"vector (e.g. {sorted(missing)[:5]}) - step labelling is offset"
            )
        # The proposal log and the sample file describe the same evaluations.
        if set(sample_step.tolist()) != set(proposal_step.tolist()):
            failures.append(
                f"rank {rank}: sampled model vectors and proposals disagree on "
                f"which steps were evaluated"
            )

    return failures, warnings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processesPerNode", type=int, default=1)
    parser.add_argument("--allowRunAsRoot", type=str, default="no")
    args, _ = parser.parse_known_args()

    # Differential evolution requires at least two chains.
    if args.processesPerNode < 2:
        print("SKIPPED: at least 2 processes per node are required for this test")
        return
    allow_run_as_root = " --allow-run-as-root" if args.allowRunAsRoot == "yes" else ""
    processes = min(args.processesPerNode, 4)

    # Clear previous output so stale files cannot mask a failure.
    subprocess.run(f"rm -rf {CHAINS}* {SAMPLES}", shell=True)
    os.makedirs("outputs", exist_ok=True)

    status = subprocess.run(
        f"export OMP_NUM_THREADS=1; cd ..; mpirun --oversubscribe -np {processes}"
        f"{allow_run_as_root} Galacticus.exe "
        f"testSuite/parameters/posteriorSampleProposalLog.xml",
        shell=True,
    )
    if status.returncode != 0:
        print("FAILED: Galacticus model failed to run")
        return

    failures, warnings = check(processes)
    for warning in warnings:
        print(f"WARNING: {warning}")
    if failures:
        for failure in failures:
            print(f"FAILED: {failure}")
        return
    print("SUCCESS: proposed-state log and sampled model vector step labelling")


if __name__ == "__main__":
    main()
    sys.exit(0)
