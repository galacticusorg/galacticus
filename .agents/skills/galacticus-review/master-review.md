# Reviewing everything merged to master since the last review

This mode is meant to soak up spare usage before a limit reset, so it may be
stopped at any moment. Work oldest PR first, write each PR's results to disk
as soon as they're done, and advance the state after each PR. A run that is cut
off part-way loses at most the PR it was on, and the next run carries on from
there.

## Files

All under `$GALACTICUS_REVIEW_DIR/master/`. If `GALACTICUS_REVIEW_DIR` is unset,
ask the user where reviews should go; don't guess. The directory may be shared
across machines (through Dropbox or a sync tool, say), so re-read a file
immediately before each write.

- `STATE.md`: the resume point. Its first line is
  `last-reviewed: <sha>`, the newest first-parent commit of `origin/master`
  whose changes have been fully reviewed. Below that, a short log, one line per
  run: date, host, range reviewed, report file, how it ended.
- `YYYY-MM-DD.md`: one report per run, named for the day the run started.
- `README.md`: what this directory is. Don't edit it during a run.

## 1. Start or resume

1. `git fetch origin` in any Galacticus checkout (the current one is fine; you
   only read `origin/master` refs, and never change its working tree).
2. Read `STATE.md`. If it does not exist, or has no `last-reviewed` line, **ask
   the user** for the starting point (a sha, a PR number, or a date) and create
   the file. Don't guess.
3. If the newest report has `Status: in progress`, check its log line in
   `STATE.md` for the host. If the report was modified in the last two hours
   (`stat`), another session may be running it: ask the user before touching it.
   Otherwise resume into that report rather than starting a new one.
4. List the work, oldest first:
   ```bash
   git log --first-parent --reverse --format='%H %s' <last-reviewed>..origin/master
   ```
   Each merge commit is one PR (number from the subject, `Merge pull request #N`).
   A non-merge first-parent commit is a direct push to master; review it as a
   unit of its own. Get each PR's title, body and author with
   `gh pr view N --json number,title,body,author,url`.
5. If the user passed `--max N` (or "just the next N"), stop after N units.
6. Write the report header with the plan (range, number of units, list) and
   `Status: in progress` *before* reviewing anything.

## 2. Each PR

- **The diff is `git diff <merge>^1 <merge>`**: exactly what that merge brought
  into master. Don't use `merge-base..merge`: branches that merged master in
  before landing make that include hundreds of files from other PRs.
- **Triage first.** Review by hand, briefly, without the full procedure: data
  or bibliography refreshes, dependabot and other dependency bumps, CI-only
  changes, and bulk parameter-file migrations. Give each a one-line verdict.
- Otherwise apply sections 2-5 of `SKILL.md` to the diff: context, physics
  and numerics, bugs, and verification. Read code at the merge commit with
  `git show <merge>:<path>`. Any build or run happens in a separate worktree
  (e.g. `git worktree add ~/galacticus-wt/master-review origin/master`), never
  in the user's checkout. Remove it at the end of the run.
- **Check each finding against current `origin/master`.** If a later commit
  already fixed it, drop it with a one-line note naming the fixing commit. If
  the line has moved, give the current `file:line` as well as the one at the
  merge.
- Review only. Don't commit fixes, open issues, or comment on the merged PRs
  unless the user asks.

After each unit, in this order:
1. Append its section to the report (format below).
2. Set `last-reviewed:` in `STATE.md` to that unit's commit sha.

## 3. Report format

```markdown
# Review of Galacticus master, <first sha>..<last sha>

**Status: in progress | done | stopped after #N (<reason>)**

Started YYYY-MM-DD on <host> in <checkout>. Range `<last-reviewed>..<origin/master sha>`:
<k> PRs and <m> direct commits.

## Summary
<filled in at the end: findings by severity, and the most important ones with links to their sections>

## #N <PR title> (<author>, merged YYYY-MM-DD, `<merge sha>`)
<one-line description of the change; "triaged: <verdict>" for a triaged PR>

### 1. [high, confirmed] <title>
- **Where:** `file.F90:123` at the merge (`file.F90:130` on current master)
- **Problem / Evidence / Suggested fix:** as in `SKILL.md` section 6

**Unconfirmed:** ...
**Checked, no problems:** ...
```

## 4. Finish

- Fill in the Summary, set `Status: done` (or `stopped after #N` and why, if
  you are stopping early), and add the run's log line to `STATE.md`.
- Tell the user: the range covered, the number of findings by severity, the top
  few, the report's path, and, if stopped early, where the next run will resume.
  Offer to turn high-severity findings into fix branches or issues.
