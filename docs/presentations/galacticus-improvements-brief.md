# Design Brief — "Galacticus: Recent Improvements" (research-group talk)

**Instructions for Claude Design.** Build a slide presentation from this brief. Everything
below is factual and specific to the Galacticus codebase; use it as the source of truth for
slide content. Where a slide calls for a live demo, a link-out, or a plot the presenter will
supply, that is flagged explicitly — leave a clearly-styled placeholder, do not invent data.

---

## 1. Deck-level parameters

- **Audience:** The presenter's own research group — students and postdocs who already use
  Galacticus. Assume fluency with the code, merger trees, and running models. You can be
  technical: show real file paths, real flags, real XML. Don't over-explain what Galacticus is.
- **Length / pacing:** ~45 minutes. Budget roughly **30–38 slides**. Six content sections plus
  title, agenda, and a closing summary. Two sections (parameter schema, PyPI) include a live
  demo — leave room and pacing for those.
- **Tone:** Technical, upbeat, "here's what got better and why it matters to your workflow."
  Every section should answer *"what changed"* and *"why you should care."*
- **Visual language:** Galacticus is an astrophysics / galaxy-formation code. Lean into a
  clean, modern dark-friendly scientific aesthetic — deep space-navy or near-black backgrounds,
  a single bright accent (electric cyan or violet), crisp monospaced type for code/paths, and a
  clean sans for prose. Keep it restrained and legible from the back of a seminar room. Use
  before/after columns, arrows, and small diagrams rather than walls of text. Code and file
  paths should always be in a monospaced, syntax-tinted block.
- **Recurring motif:** Consider a small consistent "before → after" visual pattern (left = old,
  right = new, arrow between) reused across sections, so each improvement reads the same way.
- **Timeline anchor:** Nearly all of this landed in **June 2026**. A single small timeline strip
  (early-June → early-July 2026) reused in a corner or on section dividers reinforces "all of
  this is recent, all at once."

**Suggested running order** (roughly easy/user-facing → deep/technical, ending on the two
crowd-pleasers to demo live):

1. Documentation → Read the Docs
2. `source/` hierarchical reorganization
3. Performance optimizations (with presenter's profiling plot)
4. Parameter validation via generated schema (**live VSCode demo**)
5. The parameter `changes`-file resolver
6. PyPI install (**live/quasi-live demo**)

---

## 2. Title + agenda

**Slide 1 — Title.**
- Title: **"Galacticus: What's New"** (or "Recent Improvements to Galacticus").
- Subtitle: "Docs, structure, speed, and tooling — a June 2026 roundup."
- Presenter name + date. Tasteful galaxy/large-scale-structure background image if available.

**Slide 2 — Agenda.**
- Six numbered items matching the running order above. Use one-line teasers:
  1. Docs are now on Read the Docs (and much bigger)
  2. `source/` is now a real directory tree
  3. The model runs faster — four targeted optimizations
  4. Parameter files now validate as you type
  5. Compose & migrate parameter files with change-files
  6. `pip install galacticus`

---

## 3. Section 1 — Documentation moved to Read the Docs

**Headline:** Documentation is now a single, searchable, versioned site — and it's much bigger.

**Key facts (accurate):**
- Live at **https://galacticus.readthedocs.io/** (RTD project slug `galacticus`). README carries the
  RTD badge.
- Built with **Sphinx** (furo theme), source in **reStructuredText**. Config: `.readthedocs.yaml`
  + `docs/conf.py`. Extensions: MathJax, bibtex (from `docs/Galacticus.bib`), and **Mermaid**
  flowchart diagrams.
- **Replaces the old LaTeX/PDF manual pipeline** (the retired `GALACTICUS_BUILD_DOCS` build and
  the "Documentation PDFs" release asset are gone). One HTML site instead of downloadable PDFs.
- **Wiki content was absorbed.** Installation and troubleshooting guides "existed only on the
  GitHub wiki" and were migrated in; tutorials, data-format specs, and physics overviews
  followed. RTD is now the single home; wiki pages were reduced to pointers.
- **Auto-generated reference.** Physics-class reference, glossary, and bibliography pages are
  generated **at build time** from RST docstrings embedded in the Fortran source under `source/`,
  via `scripts/doc/extractDocsRST.py`. So the class reference can't drift from the code.
- **Scope:** ~70 hand-written `.rst` pages plus the generated reference. Structure:
  - **User Guide** — getting started; installation (binary Linux/macOS, Docker, from source,
    **pip**); running; analysis; python interface; a full **troubleshooting** set (compile-time,
    run-time, debugging, tree-deadlocks, profiling); data-format specs; **18 tutorials** (halo
    mass functions incl. WDM, N-body & DMO merger trees, subhalo evolution, power spectra,
    stellar spectra/luminosities, reionization, lightcone mock catalogs, MCMC constraints, star
    formation histories, excursion-set, …).
  - **Developer Guide** — coding conventions, methods, creating a new class, traversing a merger
    tree, building Docker images, CI, versions & releases.
  - **Physics Guide** — components, definitions, and **11 Mermaid "how it works" flowcharts**
    (black holes, CGM & cooling, galactic structure, merger-tree building, outflows, star
    formation, structure formation, subhalo evolution, …).
  - **Reference (generated)** — physics classes, modules, enumerations, constants, glossary,
    references, workarounds.

**Slide plan (3–4 slides):**
- **S1.1 — "Docs have a new home."** Big screenshot/mock of the RTD landing page. One line:
  *from scattered LaTeX PDFs + a wiki → one Sphinx site.* **[LINK-OUT: open
  https://galacticus.readthedocs.io/ live in the browser here.]**
- **S1.2 — "What's in it now."** The three-guide structure (User / Developer / Physics) as three
  columns, with the tutorial count (18) and physics-flowchart count (11) called out as stat
  badges. Mention the wiki was folded in.
- **S1.3 — "It can't go stale."** Explain the build-time generation: docstrings in the Fortran
  source → `extractDocsRST.py` → class reference. Small diagram: `source/*.F90` (docstring) →
  arrow → generated RST → arrow → HTML page. This is the "why it's better than the old manual"
  slide.
- **S1.4 (optional) — Mermaid physics overview.** Drop in one of the actual Mermaid flowcharts
  (e.g. star formation or merger-tree building) as a visual payoff. **[Presenter: grab a
  screenshot from the live site.]**

---

## 4. Section 2 — `source/` hierarchical reorganization

**Headline:** `source/` went from ~1,965 files in one flat directory to a real physics-organized
tree.

**Key facts (accurate):**
- **Before:** essentially flat — **~1,965 files directly in `source/`**, with only 6
  subdirectories (all vendored libraries). Files used **dot-separated pseudo-namespace names**,
  e.g. `accretion.halo.total.simple.F90`, `utility.units.F90`.
- **After:** **62 top-level subdirectories**, only 4 loose files left at the top
  (`Galacticus.F90`, `libgalacticus.F90`, `state.F90`, `libraryClasses.xml`). The dot segments
  became **real path separators**: `accretion.halo.simple.F90` → `accretion/halo/simple.F90`.
- **Organizing principle: by physics component / functional domain.** Directories named for the
  science and infrastructure: `accretion`, `black_holes`, `cooling`, `cosmology`,
  `dark_matter_halos`, `dark_matter_profiles`, `merger_trees`, `satellites`, `star_formation`,
  `stellar_populations`, `structure_formation`, `mass_distributions`, … plus infrastructure:
  `external` (vendored third-party), `system` (OS/C glue), `data`, `numerical`, `utility`,
  `objects`, `nodes`, `tests`.
- **New convention:** within each class directory the abstract `functionClass` base is now
  `_class.F90`, with concrete implementations alongside (e.g. `accretion/halo/_class.F90`,
  `accretion/halo/simple.F90`, `accretion/halo/cold_mode.F90`).
- **Scale of the move:** the reorg commit renamed **~1,970 files** in one shot (PR #1185, June
  2026). Largest directories today: `nodes/` (273 files), `merger_trees/` (151), `tests/` (136),
  `output/` (135), `structure_formation/` (132).
- **Build system kept up.** The `Makefile` now discovers sources with a recursive-wildcard
  helper (`rwildcard`) and enumerates every subdirectory (`rsubdirs`) so `-I` include paths and
  `vpath` resolve across the tree — deliberately using only Make built-ins (no `$(shell)`) to
  stay compatible with the profiler's shell.

**Slide plan (2–3 slides):**
- **S2.1 — "Before → after."** The signature visual: left panel, a long scrolling list of
  `accretion.halo.simple.F90`-style flat filenames; right panel, a tidy nested tree. Big stat:
  **~1,965 flat files → 62 organized directories.** This slide should land purely visually.
- **S2.2 — "Organized by physics."** Show the directory tree grouped: physics domains vs
  infrastructure (`external`/`system`/`data`). Call out the `_class.F90` convention with a tiny
  example directory (`dark_matter_profiles/` or `accretion/halo/`).
- **S2.3 — "Why it matters."** Bullets: easier to navigate and onboard; find-in-directory now
  meaningful; `_class.F90` makes the abstract/concrete split obvious; the doc generator and build
  both walk the tree automatically. Note it was a mechanical, low-risk rename (git tracked all
  ~1,970 as renames) so history is preserved.

---

## 5. Section 3 — Performance optimizations

**Headline:** Four targeted optimizations, each aimed at a real profiling hotspot in the
lightcone/N-body benchmark.

> **[PRESENTER-SUPPLIED PLOT — this is the anchor of the section.]** The presenter has profiling
> data and will produce a run-time plot showing the cumulative effect of these optimizations.
> Reserve a **full-slide hero placeholder** for it (styled frame, title "Run-time impact", axis
> stubs) — Claude Design should NOT fabricate numbers. Suggested framing: before → after run time,
> ideally broken down by contribution. Put this slide either at the top of the section (as the
> "here's the payoff, now here's how" hook) or at the end (as the "and here's the total win"
> capstone) — recommend **end-of-section capstone.**

Do a short intro slide, then one slide per optimization, then the plot.

**S3.0 — Intro.** "We profiled the lightcone benchmark and attacked the top hotspots. Four
changes, all June 2026, all verified bit-for-bit identical outputs." Set expectation that these
are surgical, not a rewrite.

**S3.1 — Link-Time Optimization (LTO).**
- *What:* Turned on `-flto` across Fortran/C/C++ flags; LTO defers optimization to link time so
  the compiler can **inline and optimize across separately-compiled files** (whole-program
  optimization) instead of one file at a time. Enabled by default (`LTO=disabled` opts out);
  compiler is gfortran.
- *The interesting engineering detail (good story beat):* naive `-flto=auto` made **every**
  `gfortran` process detect **all** CPUs on a shared HPC node (e.g. 640) regardless of the
  requested `make -jN`; with many executables linking at once the parallel link workers
  exploded (N links × N jobs → resource blow-up). Fix: switch to **`-flto=jobserver`** so all
  link-time recompilation shares GNU make's single token pool, capped at `-jN`. Needed a `+`
  prefix on the executable link recipe so make exposes its jobserver.
- *Slide:* one-liner "optimize across the whole program, not file-by-file," plus a small
  before/after of the jobserver fix (640 phantom CPUs → bounded by `make -j`). Nice "gotcha we
  had to solve" beat.

**S3.2 — Object pooling for `massDistributions`.**
- *What:* NFW mass-distribution objects were being allocated and destroyed on **every** call
  (per node). Now they're drawn from a **reusable, reference-counted object pool** (`source/objects/pool.F90`,
  first used by `source/dark_matter_profiles_DMO/NFW.F90`); an object is "free" when only the pool
  still references it.
- *Why faster:* avoids constant heap allocate/deallocate churn. Measured: NFW object reuse stays
  ~100%; instructions-retired within +0.22% of the pooled baseline, versus **+7.5%** if pooling
  were lost.
- *Slide:* small diagram — repeated `allocate → use → deallocate` (red, churning) vs
  `acquire from pool → use → return` (green, steady). Stat: avoiding **~7.5%** instruction
  overhead.

**S3.3 — Meta-property pointer-returning getter (avoid big copies).**
- *What:* Added a pointer-returning sibling to the by-value meta-property getter for rank≥1
  properties. The old getter **allocated a fresh array and copied** the stored data on every
  call; the new one returns a **pointer** to the stored array — no allocation, no copy. Purely
  additive and read-only by contract.
- *Why it mattered:* the per-node, per-timestep read of the lightcone
  `positionInterpolatedCoefficients` meta-property was the **single hottest symbol** in the
  profile. Returning a pointer eliminated a heap allocation + full-array copy on every read.
- *Slide:* before/after code sketch — `coefficients = getCopy(...)` (allocates + copies) vs
  `coefficients => getReference(...)` (pointer). Note the Fortran detail that made it work
  (`self` declared `target`; name abbreviated to fit the 63-char identifier limit) as a small
  aside for the Fortran crowd.

**S3.4 — In-place star-formation-history (SFH) accumulation.**
- *What:* Three composable fixes to how per-ODE-substep SFH rates are accumulated (all in
  `source/objects/history.F90`):
  1. **Zero-skip** — skip history bins whose contribution is identically zero (the `fixedAges`
     rate writes exactly one non-zero bin out of ~51, so ~50/51 of the per-bin binary searches
     vanish).
  2. **Aligned-grid fast path** — when both histories share a time grid (the common case), the
     overlap integration collapses to a single element-wise add, eliminating the array search
     entirely.
  3. **`incrementSerialized`** — accumulate the serialized ODE rate-vector **in place** (rank-remap
     a slice to a 2D view) instead of making three full array copies per operator per substep.
- *Why faster:* this was a top benchmark hotspot. Measured (perf, 8 threads): `history_increment`
  self-time **~7.4% → ~1%**; the array-search cost went from significant to ~0.01%; memory-move
  traffic dropped too. **Verified bit-for-bit identical** across all 123 per-node datasets.
- *Slide:* the three fixes as three stacked mini-cards, each with its measured win. Emphasize
  "same answer, less work" — bit-identical is a selling point to a science audience.

**S3.5 — The profiling plot (presenter-supplied).** See the placeholder note above.

---

## 6. Section 4 — Parameter validation via a generated schema (LIVE VSCode DEMO)

**Headline:** Galacticus now generates an XSD schema from its own source, so parameter files
validate and autocomplete **as you type** in an editor.

**Key facts (accurate):**
- **Generator:** `scripts/build/parameterSchema.py` introspects the Fortran source tree — it
  discovers every `functionClass` base and its implementations, harvests each implementation's
  `<inputParameter>`/`<objectBuilder>` directive comments, and infers a type per parameter, reusing
  the existing source-tree parser. Run via `make parameters-schema`.
- **Output:** a committed **XSD 1.0** schema at `schema/parameters.xsd` (~5,000 lines).
  Deliberately XSD 1.0 so any XSD-aware editor can consume it.
- **What it enforces (smartly):** each class-selector element's `value` is restricted (via
  enumeration) to that class's valid implementation labels; enumeration-valued parameters are
  restricted to their allowed labels. Everything else is lax (unknown/global/meta parameters and
  arbitrary nesting still allowed) — it's built for **editor assistance**, not draconian rejection.
- **Stays in sync:** CI regenerates the schema and fails if the committed copy is stale
  (`git diff --exit-code`), so it can't drift from the code.
- **Deeper CI validation** (mention briefly): a separate typed catalog + Python validator
  (`scripts/build/parameterValidate.py`) does precise per-implementation checks — expands
  `xi:include`, checks `idRef` references, verifies selectors and parameter names against the
  selected implementation, and checks type/enum/range constraints. This is the strict CI gate;
  the XSD is the live-editor layer.
- **VSCode integration (this is the demo):**
  - The repo ships **`.vscode/settings.json`** that wires the schema up out of the box once the
    Red Hat **XML extension** (`redhat.vscode-xml`) is installed:
    ```json
    {
      "xml.fileAssociations": [
        { "pattern": "**/parameters/**/*.xml",
          "systemId": "${workspaceFolder}/schema/parameters.xsd" }
      ]
    }
    ```
  - Documented in `docs/manuals/user-guide/advanced.rst` → "Editor validation and autocompletion."
  - Single-file alternative for any XSD-aware editor: point the root element at the schema via
    `xsi:noNamespaceSchemaLocation` (resolves relative to the file's own directory). The repo's
    root `parameters.xml` already carries this as a demo.

**Slide plan (3 slides + demo):**
- **S4.1 — "The problem."** Old world: a typo'd parameter name or an invalid implementation label
  fails only at run time, minutes in. Show a small wrong XML snippet.
- **S4.2 — "The schema is generated from the code."** Diagram: `source/*.F90` directive comments
  → `parameterSchema.py` → `schema/parameters.xsd` → editor. Emphasize: nobody hand-maintains it;
  CI keeps it in sync.
- **S4.3 — LIVE DEMO SLIDE.** A near-empty slide that just says **"Live demo — VSCode"** with a
  short scripted checklist for the presenter (Claude Design: render these as on-screen speaker
  cues):
  - Open a parameter file under `parameters/` in VSCode (extension already installed).
  - Type a `functionClass` selector and trigger autocomplete → show the enumerated valid
    implementation labels appear.
  - Type an invalid implementation label → show the red squiggle / diagnostic instantly.
  - Hover / show the parameter name completion.
  > **[PRESENTER: VSCode runs live on your laptop for this. Have the repo open and
  > `redhat.vscode-xml` installed beforehand; open a file under a `parameters/` directory so the
  > shipped file-association matches.]**
- **S4.4 (optional) — "Also guards CI."** One line + small badge: same parameter definitions feed
  the strict Python validator that gates CI. Belt and suspenders.

---

## 7. Section 5 — The parameter `changes`-file resolver

**Headline:** Compose parameter files instead of copy-pasting them — patch a base file with
separate **change-files**, resolved the same way Galacticus resolves them internally.

**Key facts (accurate):**
- **New (June 2026, PR #1207):** a standalone **Python parameter resolver**,
  `python/Galacticus/Parameters/resolve.py`, that reproduces Galacticus' load-time processing
  **without running the Fortran binary**: **XInclude → apply change-files → validate `id`/`idRef`
  references → prune `active="…"` conditionals**. Wired into the launcher as `galacticus resolve`
  and `galacticus validate`.
- **The change-file itself:** a separate XML file with a `<changes>` root and a list of
  `<change>` elements that patch a base parameter file at load time — swap models, add/remove
  parameters, tweak values — **without editing the base file.** Heavily used in constraint
  pipelines (e.g. run one base model with sharp-k vs top-hat window functions).
- **Real example** (`constraints/pipelines/darkMatter/changesDespali2015.xml`):
  ```xml
  <changes>
    <change type="replaceWith" path="cosmologicalMassVariance"
            target="cosmologicalMassVariance/cosmologicalMassVariance"/>
    <change type="replace" path="haloMassFunction">
      <haloMassFunction value="despali2015"/>
    </change>
    <change type="update" path="outputFileName" append="true" value="_Despali2015"/>
  </changes>
  ```
- **Change types supported:** `remove`, `update` (optionally `append="true"`), `append`,
  `insertBefore`, `insertAfter`, `replace`, `replaceWith` (needs a `target` XPath), `encapsulate`
  (wrap the target), and `replaceOrAppend`. Each `<change>` takes a `type` and an XPath `path`.
- **Complementary older tool worth one mention:** `scripts/aux/parametersMigrate.py` +
  `scripts/aux/migrations.xml` **migrate old parameter files forward across versions** — every
  parameter rename/removal/restructure is keyed to the git commit that made it, and the tool
  applies only the migrations between your file's stamped revision and HEAD. Categories: rename,
  remove, value-transform, move-into-host, and complex restructurings via handler functions.
  Frame this as: **change-files = compose variants now; migrations = upgrade an old file.**

**Slide plan (2–3 slides):**
- **S5.1 — "Stop copy-pasting parameter files."** The pain: N near-identical parameter files that
  drift. The fix: one base + small change-files. Show the Despali2015 example verbatim (it's
  compact and real).
- **S5.2 — "One resolver, the real semantics."** Diagram the pipeline:
  `base.xml → XInclude → changes → reference check → conditional prune → fully-resolved.xml`.
  Emphasize it matches what the binary does internally, and you can now produce the resolved file
  with `galacticus resolve` for inspection/debugging without a run. List the change types as a
  compact grid.
- **S5.3 (optional) — "…and migrate old files forward."** Brief: `parametersMigrate.py` +
  `migrations.xml` for upgrading legacy files across renames. Keep it to one slide — it's the
  complementary point, not the headline.

---

## 8. Section 6 — PyPI install (LIVE / QUASI-LIVE DEMO)

**Headline:** `pip install galacticus` — a ready-to-run model, no compilation.

**Key facts (accurate):**
- **Package name:** **`galacticus`** on PyPI → `pip install galacticus`. README:
  *"gets you a ready-to-run model with no compilation — pre-built binaries, datasets, and tools
  are downloaded automatically for Linux and macOS."*
- **What it actually is:** NOT the Fortran source, NOT Python bindings — a pip-installable
  **launcher/installer** (`python/galacticus_launcher/`). It installs a `galacticus` console
  command. On first use it downloads, per platform, the **pre-built `Galacticus.exe`**, the
  **datasets** archive, and **tools** from GitHub Releases. No compiler needed.
- **Supported platforms:** Linux x86-64, macOS Intel, macOS Apple Silicon (M1) — one pure-Python
  wheel; the native binary is fetched at runtime, not shipped in the wheel.
- **Sub-commands:** `install`/`update`, `run <params.xml>` (`galacticus params.xml` is shorthand),
  `validate`, `resolve` (ties back to Section 5's resolver), `clean`, `info`.
- **Packaging/publish:** single PEP 621 `pyproject.toml`; version derived from git tags via
  `setuptools_scm`. Publishing is automated: `.github/workflows/publish.yml` triggers on a
  `vX.Y.Z` tag, builds sdist + wheel, and publishes to PyPI via **Trusted Publishing (OIDC — no
  stored token)**, then mirrors the binary/tool assets onto the version tag and pins the datasets
  commit.
- **Related (one mention):** analysis/plotting of outputs lives in a **separate** package,
  **`dendros`** (`pip install dendros`) — distinct from the `galacticus` launcher.

**Slide plan (2–3 slides + demo):**
- **S6.1 — "Installing used to mean building."** Contrast: clone + compile + fetch datasets +
  configure vs **`pip install galacticus`**. Big monospaced payoff line.
- **S6.2 — "What you get."** The launcher model: `pip install` → `galacticus install` pulls the
  right pre-built binary + datasets + tools for your platform → `galacticus run params.xml`. Small
  flow diagram. List the sub-commands, highlighting `run`/`validate`/`resolve` (callback to
  Sections 4–5). Note the three supported platforms.
- **S6.3 — DEMO SLIDE.** **"Live demo — pip install"** with speaker cues:
  - Show `pip install galacticus` (pre-run in a fresh venv, or show a recording — a real install
    pulls binaries and may be slow; a **pre-baked terminal or asciinema recording is safer than a
    cold live install**).
  - `galacticus info` to show the resolved platform/binary.
  - `galacticus run <small params.xml>` on a tiny model, or at least `galacticus validate` /
    `galacticus resolve` for something fast.
  > **[PRESENTER: a truly cold `pip install` downloads large binaries/datasets — for a reliable
  > demo either pre-install into a venv and show the commands, or play a short recording. Link out
  > to https://pypi.org/project/galacticus/ and the pip install docs page.]**
- **S6.4 (optional) — publish pipeline.** One line for the developers in the room: tag `vX.Y.Z` →
  GitHub Actions builds and publishes to PyPI via OIDC trusted publishing, no manual upload.

---

## 9. Closing

**Slide — "Recap."** Six one-line takeaways, matching the agenda, each with its "why you care":
- Docs → **galacticus.readthedocs.io** (one site, wiki folded in, can't go stale).
- `source/` → a real physics-organized tree (~1,965 flat files → 62 dirs).
- Faster model → four surgical optimizations, bit-identical output, real hotspots hit.
- Parameter files validate **as you type** (generated XSD + VSCode).
- Compose & migrate parameter files with **change-files** / the resolver.
- **`pip install galacticus`** — a running model in one command.

**Slide — "Links & try it."** A clean link board:
- Docs: https://galacticus.readthedocs.io/
- PyPI: https://pypi.org/project/galacticus/  → `pip install galacticus`
- Repo: github.com/galacticusorg/galacticus
- Analysis tools: `pip install dendros`
- VSCode: install `redhat.vscode-xml`, open the repo, edit any file under `parameters/`.

**Optional closing slide — "What's next / questions."**

---

## 10. Notes for Claude Design (do / don't)

- **Do** render all file paths, flags, and XML in monospaced, syntax-tinted blocks — they're
  load-bearing for this audience.
- **Do** reuse the before→after motif and the June-2026 timeline strip for cohesion.
- **Do** leave obvious, styled placeholders for: (a) the presenter's **profiling plot** (Section
  3), (b) the two **live-demo** slides (Sections 4 and 6) with the speaker-cue checklists.
- **Don't** invent profiling numbers or a run-time plot — the presenter supplies that. The only
  measured figures you may state are the ones quoted in this brief (e.g. `history_increment`
  ~7.4%→~1%; +7.5% avoided by pooling) and they should be attributed as measured on the lightcone
  benchmark.
- **Don't** over-explain Galacticus basics; this is an internal group that knows the code.
- Keep to ~30–38 slides for a 45-minute talk with two live demos.
