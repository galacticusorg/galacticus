---
name: galacticus-review
description: Review a GitHub pull request, or a branch in a local checkout, for physics/numerical correctness and for bugs, then post the review as a comment on the PR or write REVIEW.md into the checkout. Also reviews everything merged to Galacticus master since the previous such review, recording the results in a review directory. In "self" mode, the session that made the changes starts a reviewer agent and iterates with it (fixing or rebutting findings) until the reviewer accepts. Use when the user asks to review a PR, a branch, "the changes here", or recent merges to master, especially in a Galacticus repository.
argument-hint: "[PR number or URL | branch | path to checkout | master [--max N] | self [branch | path] [--rounds N]] [--draft]"
---

# Reviewing a PR or a local branch

The reader of this review is the author of the change (a human, or an agent
working in that checkout). Every finding must be something they can act on, and
every claim must have been checked. A short review of real problems beats a
long one padded with guesses.

## 1. Work out the target

| Argument | Target | Output |
|---|---|---|
| a number, `#N`, or a `github.com/.../pull/N` URL | that PR | a review comment on the PR |
| a path to a checkout | that checkout's current branch | `REVIEW.md` at its root |
| a branch name | that branch, in the current repository | `REVIEW.md` at the root of the checkout that has it checked out, or in the current one |
| nothing | the current checkout's branch | `REVIEW.md` at its root |
| `master` (or "recent merges", "since the last review") | every PR merged to `origin/master` since the previous such review | a report in `$GALACTICUS_REVIEW_DIR/master/` (ask the user for the directory if that variable is unset): **follow `master-review.md` in this skill's directory**, which uses sections 2-6 below for each PR |
| `self` [branch or path] [`--rounds N`] | the current checkout (or the one named), reviewed by a separate agent while this session responds | `REVIEW.md`, revised each round until the reviewer accepts: **follow `self-review.md` in this skill's directory**. Needs an agent that can start, and send messages to, a separate sub-agent (Claude Code can) |

`--draft` (or the user saying "don't post") means: write the PR review to the
scratchpad and show it, without posting.

**PR:**
```bash
gh pr view N --json number,title,body,author,baseRefName,headRefName,headRefOid,files,commits,url
gh pr view N --comments            # earlier review rounds: don't repeat settled points
git fetch origin "pull/N/head:refs/review/pr-N" "+refs/heads/<base>:refs/remotes/origin/<base>"
```
Read the code at the PR head without touching the user's checkout: use
`git show refs/review/pr-N:<path>`, or, if you need to build or run anything,
`git worktree add <scratch-or-work-dir>/pr-N refs/review/pr-N` (remove it when
done). Never `gh pr checkout` in the user's working tree.

**Local branch:** `git fetch origin` first (the base may have moved). The base is
`origin/master` unless the branch tracks something else or the user says
otherwise.

**The diff to review** is `git diff $(git merge-base origin/<base> <head>) <head>`
(`gh pr diff N` for a PR). For a local checkout, review uncommitted changes too
(`git diff $(git merge-base origin/<base> HEAD)` with no second ref, plus
`git status` for untracked files) and say in the output which you included. If a
branch has merged the base in, the merge-base form above already excludes the
base's own changes. Do not use `<base>..<head>` (two dots) on the code, which
picks up base commits too.

Record the exact head commit you reviewed; the output must name it.

## 2. Gather context before judging

- Read the PR description (or the branch's commit messages) and any linked
  issue. Know what the change claims to do before deciding whether it does it.
- In a Galacticus repository, read the repository's `docs/agents/fortran.md`
  (all of it) and, in `docs/agents/testing.md`, *Test coverage traps* and the
  guidance on independent references it links to. Treat their lists as a
  checklist for the changed code.
- Read each changed file in full around the hunks, not just the diff: the
  surrounding routine, the class's other methods (constructor, destructor,
  `deepCopy`, `stateStore`), and the sibling implementations of the same
  `functionClass`. Many bugs are an inconsistency between a changed line and
  an unchanged one.
- For each changed public interface, `grep` its callers.
- Skim bulk mechanical changes (parameter-file migrations, `lastModified`
  stamps, regenerated data, dependency bumps), and say that you skimmed them.

## 3. Physics and numerical correctness

For every changed equation, rate, or physical quantity:

- **Against the source.** Find the reference it cites (`!!{ !!}` docs,
  `:cite:` keys, `docs/Galacticus.bib`). Fetch the paper (arXiv/ADS) when
  reachable, and compare term by term: coefficients, exponents, signs, which
  mass/radius/density definition, log bases, h-scaling, fitting-function
  validity ranges. If the paper can't be retrieved, say so and check
  internal consistency only. Don't report a departure that a comment documents.
- **Units.** Galacticus internal units are Msun, Mpc, Gyr, km/s; constants in
  `source/numerical/constants/` are mostly SI. Check every conversion
  (`kilo`, `gigaYear`, `megaParsec`, `speedLight` in m/s, the ready-made
  `gravitationalConstant_internal` and `MpcPerKmPerSToGyr`), comoving vs.
  physical, and per-unit-mass vs. total. Do a dimensional check on each new
  expression and write it out in the finding when it fails.
- **Limits and edge cases.** r → 0 and r → ∞, zero mass/gas/metals, empty
  components, z → 0 and high z, a node with no progenitor or no parent,
  satellites vs. centrals. Look for division by zero, `log(0)`, `sqrt` of a
  negative, `exp` overflow, and silent clamping that hides a real error.
- **Numerics.** Catastrophic cancellation (differences of nearly equal
  quantities, `1-cos`, `1-exp(-x)` for small x: use `expm1`/series), loss of
  precision from single-precision literals (`1.0e-3` vs `1.0d-3`, `0.1` in a
  `double` expression), tolerances that are absolute where they should be
  relative, root finders without a valid bracket, interpolation beyond the
  tabulated range, integrals with singular or oscillatory integrands, ODE
  derivatives that are discontinuous or that drive quantities negative.
- **Conservation and bookkeeping.** Mass, metals, energy and angular momentum
  moved between components must leave one and arrive in another, at the same
  rate. Check that rates added in one place are subtracted in the other, and
  that both the rate and any "rate-limited" fraction use the same quantity.
- **Consistency.** The same physical quantity computed elsewhere in the code
  (a sibling class, a helper, an output property) should agree with this one.
  Parameter defaults should match the cited paper, or the docs should say why
  not.

Where a short calculation settles a question, do it (Python, or a computer
algebra tool if available) and quote the numbers. Where only running the code
settles it, see section 5.

## 4. Bugs

- **functionClass plumbing:** a new component missing from the constructor,
  destructor (`<objectDestructor>`), `deepCopy`, `stateStore`/`stateRestore`,
  or the parameter-driven constructor's `<inputParametersValidate>`;
  `<objectBuilder>` vs. bare pointer assignment (see the reference-count notes
  in `docs/agents/fortran.md`); state that must reset between nodes or trees
  (`autoHook`, `calculationReset`, memoized `...Previous` keys that don't
  include every argument).
- **OpenMP:** module variables used from parallel code without
  `!$omp threadprivate`; wrong `shared`/`private` clauses; derived types with
  reference-counted components in `private` clauses; unprotected updates to
  shared state; HDF5 access outside `hdf5Access` locks.
- **Fortran:** uninitialized variables (especially on one branch of an `if`),
  optional arguments used without `present()`, array bounds and off-by-one
  errors, integer overflow from default-kind arithmetic (`c_size_t`,
  `kind_int8` needed), `intent(out)` on something that must keep its value,
  `varying_string` temporaries passed as arguments (leak), `allocate(source=)`
  bypassing reference counting.
- **Error handling:** paths that silently continue, return a fallback value,
  or `print` where they should `call Error_Report(... //{introspection:location})`.
- **Copy-paste:** a block copied from a sibling class that still names the
  sibling's quantity, component, or parameter.
- **Tests, docs, and parameters** (only what CI does not already catch): new
  physics with no test, or a test that cannot fail (see *Test coverage traps*);
  a test script that exits non-zero or lacks `SUCCESS:`/`FAILED:` markers; a
  renamed or restructured parameter with no migration entry; docs or default
  values that contradict the code. Don't report formatting, spelling, or
  contributor markers that the hooks and CI already enforce.
- **Python and build scripts:** the same care applies; check them against the
  minimal CI Python environment described in `docs/agents/testing.md`.

## 5. Verify every finding before it goes in

Every candidate finding gets a deliberate attempt to refute it:

- Re-read the code path, including callers. Is the "bad" input reachable from
  a real parameter file (`parameters/`, `testSuite/`) or a real caller? Is it
  guarded elsewhere? Does a documented convention explain it?
- **Measure, don't theorize.** Where a runtime claim is cheap to check (a
  standalone calculation, a unit test, a small model), check it. Builds and
  model runs go in a separate worktree, never the user's checkout, with a `-j`
  that fits the machine's memory (see *Building* in `AGENTS.md`), run in the
  background. Ask before a build that would take more than ~30 minutes.
- Keep only findings that survive. A finding that you could neither confirm
  nor refute goes under "Unconfirmed", stated as a question, not a claim.
- Limit sub-agent use: do the review yourself; at most one sub-agent with a
  tight brief (e.g. "try to refute finding 3"). Never several in parallel.
- Before writing a PR comment, re-check each finding against the PR's
  *current* head (the author may have pushed since you fetched).

## 6. Write the review

Severity: **high** (wrong results or a crash in normal use), **medium** (wrong
in a reachable non-default configuration, or a misleading failure), **low**
(latent, cosmetic, documentation). Confidence: **confirmed** (demonstrated by a
calculation, a run, or an unambiguous code path) or **likely** (strong
evidence, not demonstrated).

Each finding: a short title; `file:line` at the reviewed head; what is wrong;
the evidence (derivation, dimensional check, numbers, the path from a real
input); a concrete suggested fix where one is clear.

Template (used for both outputs):

```markdown
# Review: <PR #N title | branch name>

Reviewed `<head sha>` against `<base>` @ `<merge-base sha>` on YYYY-MM-DD.
Scope: <files reviewed; what was skimmed; whether uncommitted changes were included>.

## Summary
<two or three sentences: what the change does, and the overall verdict>

## Findings
### 1. [high, confirmed] <title>
- **Where:** `path/to/file.F90:123`
- **Problem:** ...
- **Evidence:** ...
- **Suggested fix:** ...

## Unconfirmed / questions for the author
- ...

## Checked and found no problems
<short list of the areas you checked closely, so the reader knows the
coverage, e.g. "unit conversions in X", "deepCopy/stateStore of Y">
```

Order findings by severity. If there are no findings, say so plainly, and keep
the "Checked" section, since a clean result means little without it.

### Local checkout → `REVIEW.md`

- Write it at the checkout's root (`git rev-parse --show-toplevel`).
- If a `REVIEW.md` is already there, read it first. Replace it, and add a
  section *Since the previous review* listing each earlier finding as fixed,
  still open, or withdrawn.
- Keep it out of commits: if `git check-ignore -q REVIEW.md` fails, append
  `REVIEW.md` to `$(git rev-parse --git-common-dir)/info/exclude` (local only;
  don't edit `.gitignore`).
- Don't change any other file in the checkout.

### PR → review comment

- Post a single review with event **COMMENT**, never APPROVE or
  REQUEST_CHANGES. Start the body with a line saying it is an AI-assisted
  review by <the agent, e.g. Claude Code>, run at the request of <the user>.
- Put each finding whose line is part of the diff as an inline comment, and
  the rest in the review body:
  ```bash
  # review.json: {"commit_id": "<head sha>", "event": "COMMENT", "body": "...",
  #               "comments": [{"path": "...", "line": 123, "side": "RIGHT", "body": "..."}]}
  gh api -X POST repos/<owner>/<repo>/pulls/N/reviews --input review.json
  ```
  Build the JSON with Python (`json.dump`), not shell string pasting. A line
  outside the diff makes the whole call fail with HTTP 422; if that happens,
  move those comments into the body and retry. As a last resort,
  `gh pr review N --comment --body-file review.md`.
- Keep a copy of what was posted in the scratchpad, and give the user the
  review's URL.

## 7. Report back

Tell the user, briefly: what was reviewed (head sha), the number of findings by
severity, the top one or two, and where the review is (the `REVIEW.md` path or
the PR review URL).
