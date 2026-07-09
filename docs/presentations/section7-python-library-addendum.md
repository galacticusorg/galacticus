# Design Brief — ADDENDUM for Claude Design

## New Section 7 — Call Galacticus from Python (the `libgalacticus` library interface)

**Add this as a seventh content section.** Place it **immediately after the PyPI section** — the
two form a natural "Python story": PyPI gets you a *running model* with one command; this gets you
Galacticus's *internal physics routines* callable directly in Python. Match the visual language,
motif, and pacing of the rest of the deck.

> **Status note for the presenter (not a slide):** this work is currently on the feature branch
> `feat/library-interface-tree-access`, not yet on `master`. Frame it to the group as "landing
> now / here's what's coming," and pitch it as **June–July 2026** rather than folding it into the
> earlier June timeline. (~55 commits, +10,083 / −1,688 lines across 79 files.)

**Headline:** `libgalacticus.so` now exposes much of Galacticus to Python — `import galacticus`
and call the real Fortran routines directly. You can even build a merger tree, walk its nodes, and
read halo properties, all from Python. Ships with **8 executed tutorial notebooks**, published on
Read the Docs.

---

### Key facts (accurate — use as source of truth)

**What the user actually gets.** Building the shared library
(`make GALACTICUS_BUILD_OPTION=lib libgalacticus.so`) *also emits* an auto-generated Python module,
**`galacticus.py`**, at the repo root. Then in Python:
```python
import galacticus
galacticus.nodesInitialize()        # once, before creating trees/nodes
galacticus.verbositySet('silent')   # quiet the Fortran progress bars in a notebook
```
Every registered Galacticus class becomes a Python-callable object. You **construct objects from
parameter values and from other objects** — the same dependency-injection object graph Galacticus
builds from a parameter file, but assembled in Python. numpy arrays are marshalled in and out.

**How it works (the generation pipeline — this is the clever part).** Nobody hand-writes the
bindings; they're generated from the Fortran source:
- **Registry:** `source/libraryClasses.xml` lists which `functionClass`es are built into the
  library (e.g. `cosmologyFunctions`, `powerSpectrum`, `surveyGeometry`, `nodePropertyExtractor`,
  `mergerTreeConstructor`, …). Method signatures come from the existing Fortran `<method>` /
  `functionClass` directives in `source/`.
- **Generator:** `scripts/build/libraryInterfaces.py` emits **(a)** Fortran `bind(c)` wrapper
  subroutines (into `libgalacticus.inc`, included by the thin `source/libgalacticus.F90`) and
  **(b)** the ctypes glue `galacticus.py`.
- **Reusable generator package:** `python/LibraryInterfaces/` — `ArgSpec` (IR for one argument),
  `Pipeline` (four ordered enrichment stages: C types → C attributes → Python reassignments →
  Fortran reassignments), `Emitters` (render ArgSpec → C/Fortran/Python code), `Classification`
  (the single source of truth for "which signatures are translatable"), `Hierarchy` (derived-type
  hierarchy scan). A companion **audit tool** (`libraryInterfacesAudit.py`) reports how much of the
  API is exposed — per the commit history, method-level readiness rose to **~90%** on this branch
  *(figure from commit messages; presenter can re-run the audit for a live number).*
- So the pipeline is: **Fortran source (directives) → C `bind(c)` wrappers → ctypes
  `galacticus.py` → `import galacticus`.**

**Why it's a big deal — the interface now handles *real* Fortran signatures, not just scalar
in/out.** This branch is mostly about teaching the generator to translate the messy,
real-world argument patterns that Galacticus routines actually use:
- **Output arrays** — 1D, 2D, and logical `intent(out)`/`intent(inout)` allocatable arrays come
  back as numpy arrays; scalar `intent(out)` companions and in-place args supported alongside them.
- **Python-callable callbacks** — Galacticus can call *your* Python function: e.g. quadrature of a
  Python integrand over a 3D computational domain, or a Python-defined cross-section in a radiation
  field. (Registered/gated; Fortran has no closures, so a callback slot is process-global.)
- **Opaque handles for shared-type pointers** — routines that return a `pointer` to a shared type
  (like a merger tree) hand Python an opaque handle; Fortran keeps ownership. This is what unblocks
  `mergerTreeConstructor.construct` from Python.
- **Pointer write-back protocol** — enables the tree-walker idiom
  `node = c_void_p(None); while walker.next(node): …`.
- **Complex sized-output buffers** — `intent(out)` arrays whose dimensions are named by integer
  input args (e.g. an N×N×N grid), including the first `complex128` support — used by
  `surveyGeometry.windowFunctions`.

**"Tree access" — the loop is now closed.** Combining opaque handles + the write-back walker + a
new **`extractScalar`** convenience method on the `nodePropertyExtractor` base class
(`source/nodes/property_extractor/_class.F90`), from pure Python you can now: **build a
Cole et al. (2000) merger tree → hold its handle → walk its nodes → read per-node halo properties**
(masses, cosmic times, virial radii/velocities). That end-to-end capability is the branch's namesake
and its best "wow" moment.

**Eight executed tutorial notebooks** (`tutorials/*.ipynb`), each demonstrating real usage:
1. **Cosmology calculator** — expansion history, ages, distances; introduces the object-construction pattern.
2. **Matter power spectrum** — transfer function → growth → P(k,z), σ(M); warm-dark-matter free-streaming.
3. **Halo mass function** — δc(z), Δvir(z); Sheth–Tormen vs Tinker; collapsed fractions.
4. **Python callbacks** — Galacticus calling your Python integrands and cross-sections.
5. **Merger-tree building blocks** — M*(z), Parkinson–Cole–Helly branching rates, build cost vs resolution.
6. **Stellar populations** — Salpeter/Kroupa/Chabrier IMFs, SN II counts, lifetimes/yields, SN Ia delay-time distribution. *(needs the datasets repo.)*
7. **Survey geometries** — SDSS angular mask from mangle polygons, Vmax, FFTW 3D window functions validated against numpy. *(needs the datasets repo.)*
8. **Merger trees (capstone)** — build Cole 2000 trees, walk nodes, extract properties, main-branch assembly histories. Explicitly exercises all three new capabilities.

These notebooks are **executed in CI** (a "Python-Interface" job runs them against the freshly built
library — any cell error fails the build) and **published on Read the Docs** under a
**"Python library tutorials"** section (rendered via nbsphinx with the committed outputs; CI keeps
the outputs truthful). So the tutorials can't silently rot.

---

### Slide plan (4–5 slides + a demo)

- **S7.1 — "Galacticus, from Python."** The hook. Show the three-line snippet (`import galacticus`
  / `nodesInitialize` / construct-and-call). One sentence: *the same physics routines the model
  runs internally, now callable interactively.* Contrast with the old world (write a parameter
  file, run the whole executable, parse HDF5 output) → now: call one routine, get a numpy array.

- **S7.2 — "The bindings write themselves."** The pipeline diagram:
  `source/*.F90 directives` + `libraryClasses.xml` → `libraryInterfaces.py` →
  (`bind(c)` wrappers `+` ctypes `galacticus.py`) → `import galacticus`. Emphasize: generated, not
  hand-maintained; an audit tool tracks coverage (~90% of methods on this branch). This is the
  "why it's maintainable" slide for the developers in the room.

- **S7.3 — "It now speaks real Fortran signatures."** A compact capability grid (icon + one line
  each): output arrays (1D/2D/logical) · scalar companions & in-place args · **Python callbacks**
  (Galacticus calls *your* function) · **opaque handles** for shared objects · **pointer
  write-back** walker · **complex** sized-output buffers. Theme line: *"from scalar-in/scalar-out to
  the messy signatures Galacticus actually uses."*

- **S7.4 — "Closing the loop: merger trees in Python."** The capstone concept. A little flow:
  `mergerTreeConstructor.construct` → opaque tree handle → `walker.next(node)` loop →
  `nodePropertyExtractor.extractScalar(node)`. Show the walker idiom
  (`node = c_void_p(None); while walker.next(node): …`) and a one-line list of properties you can
  pull (mass, cosmic time, virial radius/velocity). Note this is the branch's namesake feature.

- **S7.5 — LIVE NOTEBOOK DEMO** (strongly recommended — this is the payoff). A near-empty
  **"Live demo — Jupyter"** slide with speaker cues (render as on-screen bullets):
  - Open notebook **08 (merger trees)** or **01 (cosmology calculator)** in Jupyter.
  - Build a tree / call a routine live; walk nodes; plot a main-branch assembly history or a P(k).
  - Mention: these exact notebooks run in CI on every change and are published on Read the Docs.
  > **[PRESENTER: this needs `libgalacticus.so` built on your laptop and `galacticus.py` on the
  > path. Since the notebooks are already executed and published, a safe fallback is to walk through
  > the published Read-the-Docs version, or run a pre-built notebook — a cold library build is slow.
  > Notebooks 06/07 need the datasets repo (`GALACTICUS_DATA_PATH`); notebooks 01–05/08 don't.]**
  > **[LINK-OUT: the "Python library tutorials" section on galacticus.readthedocs.io.]**

---

### Small updates to the rest of the deck

- **Agenda slide:** add a 7th item — *"7. Call Galacticus from Python (`import galacticus`)."*
- **Recap slide:** add a bullet — *"Galacticus's physics routines are now callable from Python —
  build & walk merger trees, 8 published tutorial notebooks."*
- **Links & try it board:** add — *"Python tutorials: the 'Python library tutorials' section on
  galacticus.readthedocs.io."*
- If you're showing the June-2026 timeline strip, extend it to **early July 2026** for this item and
  optionally tag it "on branch — landing now."
