#!/usr/bin/env python3
"""Assemble the report artifact HTML from figure fragments and data."""
import json

PERF = "/home/user/perf"

fig1 = open(f"{PERF}/fig1.svg.html").read()
fig2 = open(f"{PERF}/fig2.svg.html").read()
fig3 = open(f"{PERF}/fig3.svg.html").read()
fig5 = open(f"{PERF}/fig5.svg.html").read()
figdata = json.load(open(f"{PERF}/figures.json"))


def series_rows(key, label):
    rows = []
    for r in figdata[key]:
        rows.append(f"<tr><td>{r['massResolution']:.3g}</td><td>{r['ratio']:.3g}</td>"
                    f"<td>{r['nodes']:,.0f}</td><td>{r['tConstruct']:.2f}</td>"
                    f"<td>{r['tEvolve']:.2f}</td><td>{r['walks']:,}</td></tr>")
    return "\n".join(rows)


html = """<title>Galacticus Tree-Evolution Scaling</title>
<style>
  :root {
    --bg:#fafbfc; --surface:#ffffff; --ink:#14181c; --ink2:#4b5560; --muted:#8b939b;
    --grid:#e6e9ec; --axis:#c6ccd2; --accent:#2a78d6; --accent-ink:#1c5cab;
    --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a; --stall:rgba(227,73,72,0.08);
    --border:rgba(20,24,28,0.09); --codebg:#f1f4f7; --warnbg:#fdf6ec; --warnborder:#eda100;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --bg:#131619; --surface:#1a1e23; --ink:#f0f3f5; --ink2:#b9c2cb; --muted:#8b939b;
      --grid:#262b30; --axis:#3a4046; --accent:#3987e5; --accent-ink:#6da7ec;
      --s1:#3987e5; --s2:#d95926; --s3:#199e70; --stall:rgba(230,103,103,0.13);
      --border:rgba(255,255,255,0.09); --codebg:#20262c; --warnbg:#241d10; --warnborder:#c98500;
    }
  }
  :root[data-theme="dark"] {
    --bg:#131619; --surface:#1a1e23; --ink:#f0f3f5; --ink2:#b9c2cb; --muted:#8b939b;
    --grid:#262b30; --axis:#3a4046; --accent:#3987e5; --accent-ink:#6da7ec;
    --s1:#3987e5; --s2:#d95926; --s3:#199e70; --stall:rgba(230,103,103,0.13);
    --border:rgba(255,255,255,0.09); --codebg:#20262c; --warnbg:#241d10; --warnborder:#c98500;
  }
  body { background:var(--bg); color:var(--ink); font-family:"Source Sans 3",system-ui,sans-serif;
         font-size:16px; line-height:1.65; margin:0; }
  .wrap { max-width:880px; margin:0 auto; padding:48px 24px 96px; }
  h1,h2,h3 { font-family:"Spectral",Georgia,serif; line-height:1.25; text-wrap:balance; }
  h1 { font-size:34px; font-weight:600; margin:8px 0 10px; }
  h2 { font-size:23px; font-weight:600; margin:52px 0 12px; }
  h3 { font-size:17.5px; font-weight:600; margin:30px 0 8px; }
  p, li { max-width:72ch; }
  a { color:var(--accent-ink); }
  .eyebrow { font-size:12px; letter-spacing:0.12em; text-transform:uppercase; color:var(--muted);
             font-weight:600; }
  .lede { font-size:18.5px; color:var(--ink2); max-width:66ch; }
  code, .mono { font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.87em;
                background:var(--codebg); border-radius:4px; padding:1px 5px; }
  pre { background:var(--codebg); border:1px solid var(--border); border-radius:8px;
        padding:14px 16px; overflow-x:auto; }
  pre code { background:none; padding:0; }
  .tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:12px;
           margin:28px 0 8px; }
  .tile { background:var(--surface); border:1px solid var(--border); border-radius:10px;
          padding:14px 16px; }
  .tile .v { font-size:26px; font-weight:600; font-family:"Spectral",Georgia,serif; }
  .tile .k { font-size:13px; color:var(--ink2); margin-top:2px; line-height:1.4; }
  figure { margin:26px 0; background:var(--surface); border:1px solid var(--border);
           border-radius:12px; padding:18px 18px 10px; }
  figcaption { font-size:13.5px; color:var(--ink2); padding:10px 6px 8px; max-width:none; }
  .tick { font-size:12px; fill:var(--muted); }
  .axis-label { font-size:13px; fill:var(--ink2); }
  .series-label { font-size:13px; font-weight:600; }
  .guide-label { font-size:11.5px; fill:var(--muted); }
  .bar-label { font-size:13.5px; fill:var(--ink); }
  .bar-value { font-size:13.5px; fill:var(--ink2); }
  table { border-collapse:collapse; font-size:14px; margin:14px 0;
          font-variant-numeric:tabular-nums; }
  .tablewrap { overflow-x:auto; }
  th { text-align:right; font-weight:600; color:var(--ink2); padding:6px 12px;
       border-bottom:2px solid var(--axis); }
  th:first-child { text-align:left; }
  td { text-align:right; padding:5px 12px; border-bottom:1px solid var(--grid); }
  td:first-child { text-align:left; }
  .callout { background:var(--warnbg); border-left:3px solid var(--warnborder);
             border-radius:0 8px 8px 0; padding:12px 18px; margin:20px 0; max-width:72ch; }
  .strategies { counter-reset: strat; padding-left:0; }
  .strategies > li { list-style:none; counter-increment:strat; margin:18px 0; padding-left:44px;
                     position:relative; }
  .strategies > li::before { content:"S" counter(strat); position:absolute; left:0; top:1px;
                             font-family:"IBM Plex Mono",monospace; font-size:13px; font-weight:600;
                             color:var(--accent-ink); background:var(--codebg);
                             border-radius:6px; padding:2px 7px; }
  .footnote { color:var(--muted); font-size:13.5px; border-top:1px solid var(--grid);
              margin-top:56px; padding-top:16px; }
  b.hot { color:var(--accent-ink); }
</style>
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:wght@500;600&family=Source+Sans+3:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Mono:wght@400;600&display=swap">
<div class="wrap">
  <div class="eyebrow">Galacticus performance analysis &middot; master @ d6a011b</div>
  <h1>Why tree evolution steepens from 1/m to 1/m<sup>1.5</sup></h1>
  <p class="lede">The 1/m<sub>res</sub><sup>1.5</sup> regime is a crossover into a
  <em>quadratic</em> cost term driven by event-synchronized full-tree walks over the
  accumulated subhalo population &mdash; and where that crossover sits is controlled almost
  entirely by one evolver parameter.</p>

  <div class="tiles">
    <div class="tile"><div class="v">0.98&ndash;1.01</div><div class="k">measured slope, reference DMO config, to M/m<sub>res</sub> = 3&times;10<sup>5</sup></div></div>
    <div class="tile"><div class="v">(1/m)<sup>1.57</sup></div><div class="k">measured with the code-default evolver settings, M/m<sub>res</sub> = 10<sup>4</sup>&ndash;3&times;10<sup>5</sup></div></div>
    <div class="tile"><div class="v">8.3&times;</div><div class="k">evolve-time cost of the default satellite-deferral setting (0 instead of 0.75), at M/m<sub>res</sub> = 10<sup>5</sup></div></div>
    <div class="tile"><div class="v">M/m = 1.3&times;10<sup>7</sup></div><div class="k">where the reference config itself passes slope 1.5 (two-component fit)</div></div>
  </div>

  <h2>1 &middot; The measurement</h2>
  <p>Reference dark-matter-only subhalo model (orbiting satellites, Chandrasekhar dynamical
  friction, Zentner tidal stripping, Gnedin heating, merging at 0.01&thinsp;r<sub>vir</sub>,
  destruction below m<sub>res</sub>), fixed seed, single thread, prebuilt bleeding-edge binary.
  Two host masses sweep M/m<sub>res</sub> from 10<sup>3</sup> to 3&times;10<sup>6</sup>; a third
  series uses the evolver's code-default settings instead of the reference file's.</p>
  <figure>__FIG1__
  <figcaption><b>Evolve CPU time per tree.</b> Both host masses collapse onto one curve in
  M/m<sub>res</sub> (cost is host-mass invariant at fixed ratio). The reference configuration
  holds slope &asymp;1 to M/m<sub>res</sub>&nbsp;&asymp;&nbsp;3&times;10<sup>5</sup> and then bends
  (1.04 &rarr; 1.07 &rarr; 1.16); the defaults series runs at slope 1.45&ndash;1.7 through the
  practical range. Hover points for values.</figcaption></figure>

  <p>All twelve reference-series points fit
  <code>t = A(M/m) + B(M/m)&sup2;</code> to 0.019 dex rms, with
  A&nbsp;=&nbsp;1.71&times;10<sup>&minus;4</sup>&thinsp;s and
  B&nbsp;=&nbsp;1.34&times;10<sup>&minus;11</sup>&thinsp;s (B/A&nbsp;=&nbsp;7.9&times;10<sup>&minus;8</sup>).
  Local slope: 1.1 at M/m&nbsp;=&nbsp;1.4&times;10<sup>6</sup>, <b>1.5 at 1.3&times;10<sup>7</sup></b>,
  1.75 at 3.8&times;10<sup>7</sup> &mdash; a measured &ldquo;1/m<sup>1.5</sup>&rdquo; is this crossover in
  progress, and the true asymptote is quadratic. Node count stays linear
  (slope 0.95&ndash;0.97) and per-node <em>construction</em> cost stays flat (~30&thinsp;&micro;s)
  even at 3.4M nodes / 9&thinsp;GB, so the bend is algorithmic, not memory-hierarchy.</p>

  <div class="tablewrap"><table>
    <tr><th>m_res [M&#9737;]</th><th>M/m_res</th><th>nodes/tree</th><th>t_construct [s]</th><th>t_evolve [s]</th><th>walks</th></tr>
    __TABLE13__
  </table></div>
  <div class="tablewrap"><table>
    <tr><th>m_res (10&sup1;&sup2; M&#9737; host)</th><th>M/m_res</th><th>nodes/tree</th><th>t_construct [s]</th><th>t_evolve [s]</th><th>walks</th></tr>
    __TABLE12__
  </table></div>

  <h2>2 &middot; The mechanism</h2>
  <p>The standard evolver advances a tree by <b>repeated full walks</b>: every live node is
  visited, evolvable nodes compute their allowed step (<code>timeEvolveTo</code>) and are pushed
  forward until blocked; the outer loop repeats until a walk evolves nothing
  (<code>merger_trees/evolver/standard.F90:291</code>). Three properties of this loop produce the
  quadratic term:</p>
  <ul>
    <li><b>Events fragment the host's stepping.</b> A subhalo merger (triggered mid-ODE at
    r &lt; 0.01&thinsp;r<sub>vir</sub>) must synchronize with its merge target to within
    <code>timeOffsetMaximumAbsolute</code> = 10&thinsp;Myr; the target's step is capped at each
    mergee's merge time. Each of the O(S) merger/destruction events forces host-region steps
    of event size instead of the natural ~0.1/H, and <b>each fragment costs one full-tree
    walk</b>. Measured walk counts grow from 8 (M/m=10&sup3;) to 10,468 (3&times;10<sup>6</sup>).</li>
    <li><b>Every walk visits everything, and hosts re-scan their satellites.</b> Visits
    (walks &times; live nodes) scale as (1/m<sub>res</sub>)<sup>2.06</sup>; the host's
    <code>timeEvolveTo</code> scans its full satellite and mergee lists on every consideration
    (<code>standard.F90:963,990</code> &mdash; the stray comment at l.989 already suspected this),
    and each trunk promotion re-parents the entire O(S) satellite list
    (<code>node_evolver/standard.F90:1284</code>).</li>
    <li><b>Per-ODE-call overhead is ~130&thinsp;&micro;s regardless of step size</b> (locks,
    node-operator pre/post hooks, three serialization passes, scales, a freshly allocated GSL
    driver, initial-step-size search). A forced 10&thinsp;Myr micro-step costs the same as a
    1&thinsp;Gyr step.</li>
  </ul>
  <figure>__FIG3__
  <figcaption><b>Walk anatomy of one m<sub>res</sub>=10<sup>8</sup> tree</b> (throttled samples
  from the walk log). Healthy walks evolve ~80% of live nodes; in the shaded merger-sync
  cascades the evolver visits ~1,250 live nodes to advance ~20, sometimes 0. At this resolution
  these stalls are still cheap because deferred satellites skip their ODE call &mdash; remove the
  deferral and each visit pays the full ~130&thinsp;&micro;s (section 3).</figcaption></figure>
  <figure>__FIG2__
  <figcaption><b>Full-tree walk count per tree.</b> Superlinear growth (local slope up to ~1.5)
  driven by event-synchronization; the per-walk cost grows linearly with live nodes, giving the
  quadratic visit product.</figcaption></figure>

  <h2>3 &middot; The default that decides where the break happens</h2>
  <p>The reference configuration files carry
  <code>fractionTimestepSatelliteMinimum = 0.75</code>: a satellite that cannot take at least
  75% of its host-relative step defers instead of micro-stepping. <b>The code default is 0</b>,
  and stock configurations not built on the reference includes (<code>quickTest.xml</code>,
  <code>parametersProfile.xml</code>, &hellip;) run with it. With deferral off, every
  stalled-walk visit to a satellite becomes a micro ODE call with full per-call overhead:</p>
  <figure>__FIG5__
  <figcaption><b>Variants at m<sub>res</sub> = 10<sup>8</sup></b> (evolve time per tree,
  10<sup>13</sup> M&#9737; host). Identical trees and walk counts throughout &mdash; the 8.3&times;
  is pure per-visit overhead. Backtracking, the merger sync window (&times;10 either way), and
  ODE step-size reuse are all &le;4% effects here.</figcaption></figure>
  <p>Rerunning the resolution series with the default settings gives
  <b>t<sub>evolve</sub> &prop; (1/m<sub>res</sub>)<sup>1.57</sup></b> over
  M/m<sub>res</sub> = 10<sup>4</sup>&ndash;3&times;10<sup>5</sup> (local slopes 1.45 &rarr; 1.67 &rarr; 1.54, and
  15&times; the reference cost at M/m = 3&times;10<sup>5</sup>): the production-observed scaling,
  reproduced at practical resolutions. The mechanism is unchanged &mdash; the deferral simply
  multiplies the quadratic term's coefficient by ~10&sup2; and moves the slope-1.5 crossover
  from M/m ~ 10<sup>7</sup> down to ~3&times;10<sup>4</sup>.</p>
  <div class="callout"><b>Cheapest win:</b> set
  <code>&lt;fractionTimestepSatelliteMinimum value="0.75"/&gt;</code> in any subhalo-resolving
  run that doesn't already have it (or change the code default / warn on 0 with orbiting
  satellites). Measured: 8.3&times; at M/m = 10<sup>5</sup>, growing with resolution; identical
  tree structures; it is already the validated reference/validation setting.</div>

  <h2>4 &middot; Where the quadratic time goes (perf)</h2>
  <p>Self-time shares from <code>perf</code> on the release binary. The linear-term composition
  is resolution-independent (NFW velocity dispersion ~7%, log/pow ~8%, Brent root-finds in the
  heated-profile solve ~7%, GSL ODE core ~10%, malloc/free ~4%). The quadratic term's symbols
  appear only at high M/m<sub>res</sub> &mdash; and they are exactly the scan loops and list
  surgery identified in the code map:</p>
  <div class="tablewrap"><table>
    <tr><th>symbol (self-time)</th><th>M/m=10<sup>4</sup></th><th>10<sup>5</sup></th><th>10<sup>6</sup></th></tr>
    <tr><td><code>basicStandardTimeGet</code> &mdash; <code>basic%time()</code> in scan loops</td><td>&lt;0.5%</td><td>&lt;0.5%</td><td><b>4.9%</b></td></tr>
    <tr><td><code>satelliteOrbitingTimeOfMergingGet</code> &mdash; mergee scans</td><td>&lt;0.5%</td><td>&lt;0.5%</td><td><b>2.5%</b></td></tr>
    <tr><td><code>Satellite_Move_To_New_Host</code> + <code>standardPromote</code> &mdash; O(S) surgery</td><td>&lt;0.5%</td><td>&lt;0.5%</td><td><b>2.0%</b></td></tr>
    <tr><td><code>simpleTimeEvolveTo</code> &mdash; per-visit timestep evaluation</td><td>&lt;0.5%</td><td>&lt;0.5%</td><td><b>1.0%</b></td></tr>
    <tr><td>kernel page-zeroing &mdash; allocation churn</td><td>&lt;0.5%</td><td>0.8%</td><td><b>2.5%</b></td></tr>
  </table></div>
  <p>The first four sum to ~10.4% at M/m = 10<sup>6</sup> &mdash; matching the ~11% measured
  excess over the linear extrapolation there.</p>

  <h2>5 &middot; What is <em>not</em> responsible</h2>
  <ul>
    <li><b>Branch time-refinement in tree building</b>: node count is slightly <em>sub</em>linear
    in 1/m<sub>res</sub> (slope 0.95&ndash;0.97); nodes per branch don't grow.</li>
    <li><b>Construction</b>: linear throughout, a steady ~20% share.</li>
    <li><b>ODE tolerances / stiffness</b>: steps, RHS evaluations, and in-step CPU all linear to
    M/m = 3&times;10<sup>5</sup>; ~1.7 evaluations/step and ~130&thinsp;&micro;s/step throughout;
    step sizes limited by <code>tidalTensorPathIntegrated</code>.</li>
    <li><b>The merger sync window size itself</b>, backtracking, step-size reuse: &le;4%.</li>
    <li><b>Memory hierarchy</b>: flat per-node construction cost up to 9&thinsp;GB working sets.</li>
  </ul>

  <h2>6 &middot; Fixes, ranked</h2>
  <ol class="strategies">
    <li><b>Keep satellites deferring (config now; consider changing the default).</b>
    Measured 8.3&times; at M/m = 10<sup>5</sup>; moves the entire superlinear regime ~2.5 decades
    deeper in resolution. Zero code risk.</li>
    <li><b>Event-driven evolution (priority queue) &mdash; the structural fix.</b> Replace
    &ldquo;walk everything until nothing evolves&rdquo; with a time-ordered queue of evolvable
    nodes plus per-node waiter lists. The dependency information already exists: <code>timeEvolveTo</code>
    reports which node blocks each blocked node (<code>nodeLock</code>, today used only for deadlock
    reports) &mdash; invert it, and wake exactly the waiters when a node advances or an event fires.
    O((evolutions+wakes)&thinsp;log&thinsp;N) replaces O(walks&thinsp;&times;&thinsp;N<sub>live</sub>);
    the stalled-walk regime (1,250 visits per ~20 evolutions) disappears.</li>
    <li><b>O(1) satellite bookkeeping.</b> Cache the host's minimum-satellite time (lazy
    recompute when the argmin advances); keep a <code>lastSatellite</code> tail pointer or
    doubly-linked satellite list (O(1) insert/remove); splice satellite lists on trunk promotion
    instead of re-walking and re-parenting O(S) nodes. Directly removes the 10% perf excess
    measured at M/m = 10<sup>6</sup>, and its growth.</li>
    <li><b>Shrink the per-ODE-call constant.</b> Reuse a per-thread GSL driver keyed by system
    dimension instead of alloc/free per call; skip re-serialization scaffolding when the
    property count is unchanged (the common case for a satellite stepping repeatedly);
    <code>reuseODEStepSize=true</code> is a free ~4%. Scales down both the linear term and the
    residual micro-step cost.</li>
    <li><b>Batch merger synchronization.</b> Complete all mergers falling within one sync window
    in a single host step instead of one walk per merger &mdash; attacks the walk count where S1&ndash;S2
    are not enough (cluster-mass hosts).</li>
    <li><b>Cheapen the satellite RHS (linear-term constant).</b> Memoize/interpolate the
    heated-profile initial-radius solve (currently a Brent root-find per RHS evaluation) and NFW
    dispersion; avoid per-call mass-distribution copies (<code>memmove</code> in
    <code>darkMatterOnlyGet</code>). Plausibly 1.5&ndash;2&times; on subhalo-heavy models.</li>
    <li><b>Physics-level coarse-graining where acceptable.</b> The validation configuration's
    <code>mergerTreeBuildController=subsample</code> below a mass threshold is the established
    way to cut S itself.</li>
  </ol>

  <h2>7 &middot; Methods &amp; reproduction</h2>
  <p>Everything lives on branch
  <code>claude/galacticus-performance-analysis-giadol</code> under
  <code>performanceAnalysis/</code>: parameter generator (reference includes, E&amp;H transfer
  function per AB), run/collect/analysis scripts, and per-run data
  (<code>data/results.jsonl</code>). Timing instrumentation:
  <code>mergerTreeOperator=treeProcessingTimer</code> (per-tree construct/evolve CPU + node
  count in <code>metaData/treeTiming</code>), the <code>simple</code> evolve profiler
  (step/evaluation/CPU histograms), <code>verbosityLevel=warn</code> walk logs, and
  <code>perf record</code> on the static release binary. 4 trees per point (2 at the two
  deepest), single OpenMP thread, Xeon @ 2.1&thinsp;GHz.</p>
  <p class="footnote">Analysis and text generated by Claude (AI) from runs of the unmodified
  bleeding-edge binary at master d6a011b, 2026-08-29/30; per CONTRIBUTING.md this analysis has
  not yet been human-verified. Absolute timings are machine-specific; slopes, ratios, and
  B/A are the transferable results.</p>
</div>
"""

html = html.replace("__FIG1__", fig1).replace("__FIG2__", fig2)
html = html.replace("__FIG3__", fig3).replace("__FIG5__", fig5)
html = html.replace("__TABLE13__", series_rows("series1e13", "1e13"))
html = html.replace("__TABLE12__", series_rows("series1e12", "1e12"))

open(f"{PERF}/report.html", "w").write(html)
print(f"report.html written: {len(html)} bytes")
