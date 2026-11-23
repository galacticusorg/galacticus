# Committing, migrating parameters, and pushing

*Read this before committing, updating a branch from master, changing a parameter, or pushing.* The short
rules are in the root `AGENTS.md`; this file has the detail and the evidence behind them. These entries
record traps found in practice; each names when it was last verified, so check that a named file, flag or
command still exists before relying on it.

## Merge master in; never rebase a branch carrying a migrations entry
*Added 2026-09-16; verified 2026-10-30*

`scripts/aux/migrations.xml` registers each migration against a **commit hash**, and `parametersMigrate.py`
matches it by exact hash against the ancestry of the parameter file's last-modified revision. Rebasing
rewrites that hash, so the entry becomes unreachable and the migration silently never fires: no error, no
warning.

The same applies to the `lastModified` stamps a migration writes: they record HEAD at run time, so amending
or rebasing that commit afterwards orphans every stamp (712 files in one checkout, 2026-10-14). Finish all
history edits *before* migrating, then check with `scripts/aux/auditRevisions.py --check` (it scans tracked
files only).

When master and the branch have both bulk-migrated, every bundled file conflicts on its stamp (716 files,
2026-10-19). Resolve by keeping *both* sides' content and either stamp. Do not re-run the migration instead:
`git rev-list --ancestry-path` never reaches a migration on the other side of the merge, so its element
would silently be lost.

A **new** parameter file added on a branch has the same blind spot (2026-10-27, PR #1690): its stamp is not
a descendant of migrations that master gains from other branches, so after updating from master
`parametersMigrate.py` applies nothing. The model then runs with an "unrecognized parameter" warning, and
the CI Test-Preprocessor validation fails. Re-migrate each new file with
`--lastModifiedRevision <earliest new migration>^ --validate no` (the structural check rejects most test
files for unrelated reasons), and grep run logs for `unrecognized parameter`, not just for failure markers.

**Why:** "update the branch against master" ordinarily means rebase, and here that would quietly disable a
migration.

**How to apply:** merge `origin/master` *into* the branch. This matches how the project takes PRs anyway
(merge commits), so branch hashes survive into master. Git's default merge message fails the `commit-msg`
hook, which enforces Conventional Commits, so give the merge a message such as
`chore(merge): bring master into <branch>`. Do not use `--no-verify`. Afterwards, check the hash is still
reachable with `git merge-base --is-ancestor <hash> HEAD`, and re-run a migration to confirm it fires.

## A parameter change needs a two-commit migration
*Added 2026-09-21; verified 2026-09-21*

`scripts/aux/migrations.xml` keys each entry to the hash of the commit that changed the parameters, so that
entry can only go in a *follow-up* commit. If you use the optional Claude review hook, it blocks the first
commit for the unregistered migration: override that one commit with `SKIP_CLAUDE=1` (master's
`d8a194faa0`/`8cd8171137` follow the same pattern), and make the second with hooks active.

Parameters *removed* by a migration must also be excepted in the test fixture. CI's Test-Preprocessor job
runs `scripts/build/parameterValidate.py` with `--exceptions testSuite/parameterValidationExceptions.json`,
which otherwise rejects the deliberately outdated `testSuite/parameters/parameterMigration.xml`. Add one
entry per removed parameter (file, kind `unknown`, name, reason). To reproduce locally, build the catalog
first with `scripts/build/parameterCatalog.py`.

**Why:** the hook rejection looks like a fault in the commit, and the CI failure appears in a job that does
not name the migration.

**How to apply:** plan two commits whenever a commit renames or removes a parameter.

## Never write a full commit hash by hand
*Added 2026-10-21; verified 2026-10-23*

Expanding a known 10-character prefix to 40 characters invents the rest of the hash. This applies wherever
a full hash is written: `<lastModified revision="...">` stamps, `<migration commit="...">` entries in
`scripts/aux/migrations.xml`, and test parameter files. The CI *Validate-Revisions* job
(`scripts/aux/auditRevisions.py --check`) fails on a stamp that is not a real ancestor commit, but no local
hook catches it. A bad `migrations.xml` hash raises no error at all: the migration just never fires.

**Why:** PR #1658 went red for exactly this (2026-10-21). On 2026-10-23 a hand-expanded `56b5f82d59…` in a
test file made `parametersMigrate.py` fail with "Invalid revision range".

**How to apply:** substitute the hash directly, e.g. `revision="$(git rev-parse <base>)"` or
`H=$(git rev-parse HEAD)`, and never retype one. Run `python3 scripts/aux/auditRevisions.py --check` before
pushing a new parameter file. After adding a migration, run it on an old-style file to confirm it fires.

## Bulk migration with `migrateAllParameterFiles.py`
*Added 2026-10-14; verified 2026-10-21*

The tool (in `galacticusDevTools/maintenance/`; see `AGENTS.md`) migrates each file on its own. So a
parameter added to two `parameters/reference/` files that models include together ends up duplicated, with
the warning "multiple copies of parameter ... only the first will be utilized" (fixed for
`subhaloOrbits.xml`/`evolutionGalaxyFormation.xml` in PR #1662, 2026-10-21).

**How to apply:** after a bulk run, start a model that includes both files and grep its log for
`multiple copies`. Validate with `parameterValidate.py` (catalog and exceptions, as in `cicd.yml`) against a
`git stash` baseline, since the regression and constraint files carry pre-existing findings.

## The git hooks are slow, and the hooks repository has its own copies of some checks
*Added 2026-07-02; verified 2026-09-13*

- `git commit` runs `pre-commit.d/00-galacticus` (static analysis, formatting, spelling, schema checks and
  more) and, if enabled, the optional Claude review hook. A commit can take several minutes.
- The optional Claude review hook sometimes hits its turn limit and prints a JSON error, but still allows
  the commit. That is not a finding.
- `00-galacticus` has its *own* copy of the embedded-fragment checks, separate from
  `scripts/aux/embeddedAnalyzer.py`, which CI runs. A change to one usually needs the same change to the
  other, in the [`gitHooks`](https://github.com/galacticusorg/gitHooks) repository.

**Why:** a foreground commit can hit an agent's tool timeout and leave changes staged but uncommitted, which
looks like a hook rejection.

**How to apply:** run `git commit` and `git push` in the background, logging to a file, and check
`git log` before retrying.

## Checks CI enforces that the hooks miss
*Added 2026-09-17 to 2026-11-10; verified 2026-11-10*

Each of these passed every local hook and then failed a PR in CI:

- **A new `:cite:` key must be in `docs/Galacticus.bib`.** A key with no entry fails the
  *Validate-Docstrings-RST* job with `unknown citation key: <key>`, while the hook's "Validate bibliography"
  check passes (2026-09-16). Run `python3 scripts/aux/docBlockAudit.py --check` (it is not executable, so
  `./` fails with "Permission denied"); it audits every source file and the manuals in seconds. Add missing
  entries in alphabetical order by key. Entries are written `@Article{<tab>  key,`, so a hand check with
  `grep "{key,"` gives a false "missing".
- **Every `<methods>` (and similar) directive in a new file needs `docformat="rst"`.**
  *Validate-Docstrings-RST* also runs `scripts/doc/convertDocstringsToRST.py --check source`, which fails
  ("contain old-style (LaTeX) docstrings") on any directive with descriptions but no `docformat="rst"`, e.g.
  a bare `<methods>` block inside a type. The hook's "Residual LaTeX in docstrings" check and
  `docBlockAudit.py --check` both pass such a file (PR #1776). Copy the new file into an empty `source/`
  directory somewhere temporary and run `python3 scripts/doc/convertDocstringsToRST.py --check` on that
  directory. Without `--check` it rewrites files, so never point it at the checkout.
- **Non-ASCII characters in XML files are rejected, even in comments.** `00-galacticus` fails a commit with
  "XML file has non-ascii character" for any `Å`, `≤` or similar character in a staged `.xml` file, although
  `xmllint` accepts it. Write "Angstroms" and "at most" in parameter-file comments, or check first with
  `grep -nP '[^\x00-\x7F]' file.xml`.
- **Warnings fail CI.** CI fails on unsuppressed compiler warnings that a local build only prints. Grep build
  logs for `Warning:` as well as errors (gfortran quotes names with curly quotes: `‘name’`).

**How to apply:** run the "Before you push" checklist in `AGENTS.md`.
