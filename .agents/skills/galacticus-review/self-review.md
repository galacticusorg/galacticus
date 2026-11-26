# Self review: iterate with a reviewer agent until it accepts

In this mode the session that made the changes (the **author**) does not review
them itself. It starts one **reviewer** agent, which reviews the checkout with
sections 2-6 of `SKILL.md` and writes `REVIEW.md`. The author then fixes or
rebuts each finding, and the reviewer checks the response. They repeat until the
reviewer accepts, or until a stopping rule hands the decision to the user.

The value of this mode is a review that is not anchored on the author's
reasoning. Everything below protects that: the reviewer gets the code, never the
author's conclusions, and the author's rebuttals must carry evidence the
reviewer can check for itself.

This mode needs an agent that can start a separate sub-agent with a fresh
context and keep sending it messages. The tool names below are Claude Code's
(`Agent`, `SendMessage`). In an agent without those, review the branch from a
fresh session instead (the default mode), and respond to its `REVIEW.md` by hand.

## Arguments

`self [path | branch] [--rounds N]`. The target is resolved as in section 1 of
`SKILL.md` (default: the current checkout, including uncommitted changes).
`--rounds` caps the number of review rounds (default 4). A PR target is not
supported in this mode: the author is a local session, so review the branch.

## 1. Before starting

- Make sure the work is in a reviewable state: it builds, and the tests you
  meant to run have run. A reviewer finding "this doesn't compile" wastes a
  round.
- Record the base and head (`git rev-parse HEAD`, the merge-base with
  `origin/<base>`), and whether uncommitted changes are included.
- If a `REVIEW.md` already exists at the root, leave it: the reviewer reads it
  as the previous review (section 6 of `SKILL.md`).
- Don't commit or push as part of this mode unless the user asked for that.

## 2. Start the reviewer

Start exactly **one** agent (in Claude Code, `subagent_type: general-purpose`,
not `fork`: a fork inherits the author's context, which defeats the point).
Never start a second reviewer in parallel, and never block waiting on its
output; wait for the completion notification. The brief, filled in:

> You are reviewing changes made by another agent session, which will
> respond to your findings. Read `.agents/skills/galacticus-review/SKILL.md`
> in the checkout at `<checkout path>` and follow sections 1-6 for a local
> checkout, with these differences:
>
> - Target: `<checkout path>`, branch `<branch>`, head `<sha>`, base
>   `origin/<base>` @ `<merge-base sha>`, uncommitted changes
>   <included / not included>. What the change is for: <one or two sentences
>   from the user's request or the commit messages, describing the goal, not
>   how it was done or why it is correct>.
> - Write `REVIEW.md` at the checkout root as section 6 describes. Change no
>   other file in the checkout. You are read-only there apart from `REVIEW.md`.
> - Do not run `make` or a model in that checkout or anywhere else: the author
>   may be building, and two builds at once can exhaust memory. If a finding
>   needs a build, a run, or a test to settle, put it in a section
>   **Requested checks** with the exact command and what result would confirm
>   or refute the finding. The author will run it and send you the raw output.
>   Standalone calculations (Python, a few-line script in your scratchpad) are
>   fine to run yourself.
> - End your reply to me with a short summary (number of open findings by
>   severity, the requested checks) and a final line that is exactly one of
>   `VERDICT: ACCEPT` or `VERDICT: CHANGES REQUESTED`.
> - Accept only when no high or medium finding is open, every requested check
>   has come back and been judged, and every question under "Unconfirmed" has
>   been answered. A low finding may stay open if the author has deferred it
>   with a reason you consider sound; list it in `REVIEW.md` as deferred.
> - In later rounds you will receive the author's response. Do not take it on
>   trust. For a fix, re-read the changed code (the diff may have grown) and
>   check that it fixes the problem without introducing a new one. For a
>   rebuttal, check its evidence yourself; withdraw the finding only if the
>   evidence holds, and say why. Review any new changes in the diff with the
>   same care as the first round.

## 3. Respond to a round

When the reviewer finishes, read `REVIEW.md` (the file, not just the summary).
If the verdict is `ACCEPT`, go to section 5. Otherwise, for **each** open
finding and question, decide on one of:

- **Fixed:** make the change, then verify it the same way you would verify any
  other change (rebuild, rerun the relevant tests, measure). State where the
  fix is and what verified it.
- **Disputed:** the finding is wrong. Give evidence the reviewer can check: the
  caller or guard it missed (`file:line`), a calculation with numbers, a run's
  output. "This is intended" is not evidence; a comment or a reference that
  documents the intention is.
- **Deferred** (low findings only): a sound reason it is out of scope here,
  for example a pre-existing problem that the change does not touch. Offer to
  open an issue for it, but don't open one without the user's say-so.
- **For the user:** the finding turns on a physics or design choice that is the
  user's to make (resolution, which physics to include, an API change). Don't
  guess; collect these for section 4.

Run each **requested check** exactly as asked, and send its raw output (the
relevant lines, not a paraphrase). If a check is too expensive (a build or run
of more than ~30 minutes), ask the user first.

Then continue the **same** reviewer (in Claude Code, with `SendMessage`, a
deferred tool: load it with `ToolSearch` `select:SendMessage` first), so it
keeps its context. The message: the new head sha and whether uncommitted
changes are included, then one entry per finding, numbered as in `REVIEW.md`,
with its disposition and evidence, then the requested checks' outputs. Ask it to
update `REVIEW.md` and give a new verdict.

If the reviewer agent has gone (it errored, or was killed by a rate limit),
start a fresh one with the brief from section 2, adding: "`REVIEW.md` holds the
previous round's review; the author's response to it follows", then the
response.

## 4. Stopping rules

Stop iterating and go to the user when any of these happens:

- the round cap is reached;
- the same finding has been disputed twice and the reviewer still stands by it
  (a genuine disagreement: show the user both sides, briefly, and let them
  decide);
- a finding was marked "for the user";
- a round produced no change: no fixes, and the reviewer withdrew nothing.

When the user decides a disputed point, tell the reviewer so it records the
decision in `REVIEW.md`, and carry on if rounds remain.

## 5. Report back

Tell the user, briefly:

- the final verdict and the number of rounds;
- the head sha the reviewer last reviewed, and whether uncommitted changes were
  included;
- what was fixed, with the one or two most significant fixes named;
- what was disputed and withdrawn, and what was deferred;
- anything left for the user to decide;
- where the record is (`REVIEW.md`, which carries the per-round history in its
  *Since the previous review* sections).

The reviewer accepting the changes is not the human verification that
`CONTRIBUTING.md` requires for AI-generated work. Say so if the user is about to
commit or open a PR on the strength of it.
