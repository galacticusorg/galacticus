#!/usr/bin/env python3
"""Generate theme-aware inline-SVG figures for the report artifact.

Reads figures.json (+ perfshares.json when present) and writes fig*.svg
fragments (SVG elements meant to be inlined in HTML; colors reference CSS
custom properties so light/dark theming is handled by the page).
"""
import json
import math

PERF = "/home/user/perf"

# Geometry.
W, H = 760, 440
ML, MR, MT, MB = 64, 16, 20, 46


def sx(x, x0, x1):
    return ML + (x - x0) / (x1 - x0) * (W - ML - MR)


def sy(y, y0, y1):
    return H - MB - (y - y0) / (y1 - y0) * (H - MT - MB)


def fmt_pow(v):
    e = int(round(math.log10(v)))
    sup = str(e).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
    return f"10{sup}"


def log_axes(x0, x1, y0, y1, xlabel, ylabel):
    parts = []
    # Gridlines + ticks at decades.
    for e in range(math.ceil(x0), math.floor(x1) + 1):
        X = sx(e, x0, x1)
        parts.append(f'<line x1="{X:.1f}" y1="{MT}" x2="{X:.1f}" y2="{H-MB}" stroke="var(--grid)" stroke-width="1"/>')
        parts.append(f'<text x="{X:.1f}" y="{H-MB+18}" text-anchor="middle" class="tick">{fmt_pow(10**e)}</text>')
    for e in range(math.ceil(y0), math.floor(y1) + 1):
        Y = sy(e, y0, y1)
        parts.append(f'<line x1="{ML}" y1="{Y:.1f}" x2="{W-MR}" y2="{Y:.1f}" stroke="var(--grid)" stroke-width="1"/>')
        parts.append(f'<text x="{ML-8}" y="{Y+4:.1f}" text-anchor="end" class="tick">{fmt_pow(10**e)}</text>')
    parts.append(f'<line x1="{ML}" y1="{H-MB}" x2="{W-MR}" y2="{H-MB}" stroke="var(--axis)" stroke-width="1"/>')
    parts.append(f'<text x="{(ML+W-MR)/2}" y="{H-8}" text-anchor="middle" class="axis-label">{xlabel}</text>')
    parts.append(f'<text x="16" y="{(MT+H-MB)/2}" text-anchor="middle" class="axis-label" transform="rotate(-90 16 {(MT+H-MB)/2})">{ylabel}</text>')
    return parts


def polyline(pts, x0, x1, y0, y1, color, dash=None, width=2):
    d = " ".join(f"{sx(x,x0,x1):.1f},{sy(y,y0,y1):.1f}" for x, y in pts)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"{dash_attr}/>'


def markers(pts, x0, x1, y0, y1, color, r=4, titles=None):
    out = []
    for i, (x, y) in enumerate(pts):
        tip = f"<title>{titles[i]}</title>" if titles else ""
        out.append(f'<circle cx="{sx(x,x0,x1):.1f}" cy="{sy(y,y0,y1):.1f}" r="{r}" fill="{color}" stroke="var(--surface)" stroke-width="2" pointer-events="all">{tip}</circle>')
    return "".join(out)


def svg_wrap(name, parts, title):
    body = "\n".join(parts)
    return (f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{title}" '
            f'style="width:100%;height:auto;display:block">{body}</svg>')


def guide(x_anchor, y_anchor, slope, x0, x1, y0, y1, label, labelAt=1.0):
    """Dashed power-law guide through (x_anchor,y_anchor), clipped to the plot box."""
    xa, xb = x0 + 0.1, x1 - 0.1
    ya0, yb0 = y0 + 0.1, y1 - 0.1

    def yat(x):
        return y_anchor + slope * (x - x_anchor)

    # Clip to the y range.
    if yat(xa) < ya0:
        xa = x_anchor + (ya0 - y_anchor) / slope
    if yat(xb) > yb0:
        xb = x_anchor + (yb0 - y_anchor) / slope
    pts = [(xa, yat(xa)), (xb, yat(xb))]
    xL = xa + labelAt * (xb - xa)
    lx, ly = sx(xL, x0, x1) - 4, sy(yat(xL), y0, y1) - 8
    return (polyline(pts, x0, x1, y0, y1, "var(--muted)", dash="5 5", width=1.5)
            + f'<text x="{lx}" y="{ly}" text-anchor="end" class="guide-label">{label}</text>')


def fig_scaling(fig):
    """F1: t_evolve vs M/m for the three series + slope guides."""
    series = [("series1e13", "var(--s1)", "10¹³ M☉ host"),
              ("series1e12", "var(--s2)", "10¹² M☉ host"),
              ("seriesDefaults", "var(--s3)", "default evolver settings")]
    x0, x1 = 2.8, 6.8
    y0, y1 = -1.2, 3.4
    parts = log_axes(x0, x1, y0, y1, "M_host / m_res", "evolve time per tree [s]")
    parts.append(guide(3.0, math.log10(0.18), 1.0, x0, x1, y0, y1, "slope 1"))
    parts.append(guide(5.0, math.log10(17.1), 1.5, x0, x1, y0, y1, "slope 1.5"))
    for key, color, label in series:
        pts = [(math.log10(r["ratio"]), math.log10(r["tEvolve"])) for r in fig.get(key, []) if r["tEvolve"] > 0]
        if not pts:
            continue
        parts.append(polyline(pts, x0, x1, y0, y1, color))
        tips = [f"M/m = {r['ratio']:.3g}: {r['tEvolve']:.3g} s/tree" for r in fig.get(key, []) if r["tEvolve"] > 0]
        parts.append(markers(pts, x0, x1, y0, y1, color, titles=tips))
    parts.append(f'<text x="{ML+10}" y="{MT+16}" class="series-label" fill="var(--s1)">10¹³ M☉ host</text>')
    parts.append(f'<text x="{ML+10}" y="{MT+34}" class="series-label" fill="var(--s2)">10¹² M☉ host</text>')
    parts.append(f'<text x="{ML+10}" y="{MT+52}" class="series-label" fill="var(--s3)">default evolver settings</text>')
    return svg_wrap("fig1", parts, "Evolve time per tree versus resolution ratio")


def fig_walks(fig):
    """F2: walk count vs M/m for the series + slope guides."""
    x0, x1 = 2.8, 6.8
    y0, y1 = 0.5, 4.2
    parts = log_axes(x0, x1, y0, y1, "M_host / m_res", "full-tree walks per tree")
    parts.append(guide(3.2, 1.0, 1.0, x0, x1, y0, y1, "slope 1", labelAt=0.55))
    parts.append(guide(5.3, math.log10(195.0), 1.5, x0, x1, y0, y1, "slope 1.5"))
    for key, color, label in [("series1e13", "var(--s1)", "10¹³ M☉ host"),
                              ("series1e12", "var(--s2)", "10¹² M☉ host"),
                              ("seriesDefaults", "var(--s3)", "default settings")]:
        pts = [(math.log10(r["ratio"]), math.log10(r["walks"])) for r in fig.get(key, []) if r.get("walks", 0) > 0]
        if not pts:
            continue
        parts.append(polyline(pts, x0, x1, y0, y1, color))
        tips = [f"M/m = {r['ratio']:.3g}: {r['walks']} walks" for r in fig.get(key, []) if r.get("walks", 0) > 0]
        parts.append(markers(pts, x0, x1, y0, y1, color, titles=tips))
    parts.append(f'<text x="{ML+10}" y="{MT+16}" class="series-label" fill="var(--s1)">10¹³ M☉ host</text>')
    parts.append(f'<text x="{ML+10}" y="{MT+34}" class="series-label" fill="var(--s2)">10¹² M☉ host</text>')
    parts.append(f'<text x="{ML+10}" y="{MT+52}" class="series-label" fill="var(--s3)">default settings</text>')
    return svg_wrap("fig2", parts, "Tree-walk count versus resolution ratio")


def fig_stall(fig):
    """F3: live nodes and nodes evolved per sampled walk, one m_res=1e8 tree."""
    st = fig.get("stall")
    if not st:
        return ""
    walks, live, evolved = st["walks"], st["live"], st["evolved"]
    x0, x1 = 0, max(walks) * 1.02
    y0, y1 = 0, max(live) * 1.1

    def SX(x):
        return ML + (x - x0) / (x1 - x0) * (W - ML - MR)

    def SY(y):
        return H - MB - (y - y0) / (y1 - y0) * (H - MT - MB)

    parts = []
    for gy in range(0, int(y1), 400):
        parts.append(f'<line x1="{ML}" y1="{SY(gy):.1f}" x2="{W-MR}" y2="{SY(gy):.1f}" stroke="var(--grid)" stroke-width="1"/>')
        parts.append(f'<text x="{ML-8}" y="{SY(gy)+4:.1f}" text-anchor="end" class="tick">{gy}</text>')
    for gx in range(0, int(x1) + 1, 50):
        parts.append(f'<text x="{SX(gx):.1f}" y="{H-MB+18}" text-anchor="middle" class="tick">{gx}</text>')
    parts.append(f'<line x1="{ML}" y1="{H-MB}" x2="{W-MR}" y2="{H-MB}" stroke="var(--axis)" stroke-width="1"/>')
    parts.append(f'<text x="{(ML+W-MR)/2}" y="{H-8}" text-anchor="middle" class="axis-label">tree walk number</text>')
    parts.append(f'<text x="16" y="{(MT+H-MB)/2}" text-anchor="middle" class="axis-label" transform="rotate(-90 16 {(MT+H-MB)/2})">nodes</text>')
    # Shade stalled stretches (evolved < 10% of live).
    for i, w in enumerate(walks):
        if evolved[i] < 0.10 * max(live[i], 1):
            wPrev = walks[i - 1] if i > 0 else 0
            parts.append(f'<rect x="{SX(wPrev):.1f}" y="{MT}" width="{max(SX(w)-SX(wPrev),1):.1f}" height="{H-MT-MB}" fill="var(--stall)" />')
    for vals, color in [(live, "var(--s1)"), (evolved, "var(--s2)")]:
        d = " ".join(f"{SX(w):.1f},{SY(v):.1f}" for w, v in zip(walks, vals))
        parts.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round"/>')
    parts.append(f'<text x="{W-MR-10}" y="{MT+16}" text-anchor="end" class="series-label" fill="var(--s1)">live nodes in tree</text>')
    parts.append(f'<text x="{W-MR-10}" y="{MT+34}" text-anchor="end" class="series-label" fill="var(--s2)">nodes evolved in walk</text>')
    return svg_wrap("fig3", parts, "Per-walk live and evolved node counts, one tree at m_res = 1e8")


def fig_variants(fig):
    """F5: variant evolve times as horizontal bars."""
    order = [("base_m1.0e8", "baseline (reference config)"),
             ("varNoDefer_m1.0e8", "satellite deferral off (code default)"),
             ("varNoBacktrack_m1.0e8", "backtrackToSatellites = false"),
             ("varReuseStep_m1.0e8", "reuseODEStepSize = true"),
             ("varSyncLoose_m1.0e8", "merger sync window ×10"),
             ("varSyncTight_m1.0e8", "merger sync window ÷10")]
    v = fig.get("variants", {})
    rows = [(label, v[k]["tEvolve"]) for k, label in order if k in v]
    if not rows:
        return ""
    bw, gap, top = 26, 14, 16
    height = top + len(rows) * (bw + gap) + 30
    vmax = max(t for _, t in rows) * 1.15
    LW = 340
    parts = []
    for i, (label, t) in enumerate(rows):
        y = top + i * (bw + gap)
        w = (t / vmax) * (W - LW - MR)
        parts.append(f'<text x="{LW-10}" y="{y+bw/2+4}" text-anchor="end" class="bar-label">{label}</text>')
        parts.append(f'<rect x="{LW}" y="{y}" width="{w:.1f}" height="{bw}" rx="4" fill="var(--s1)"/>')
        parts.append(f'<text x="{LW+w+8:.1f}" y="{y+bw/2+4}" class="bar-value">{t:.1f} s</text>')
    return (f'<svg viewBox="0 0 {W} {height}" role="img" aria-label="Variant comparison at m_res = 1e8" '
            f'style="width:100%;height:auto;display:block">{"".join(parts)}</svg>')


def main():
    fig = json.load(open(f"{PERF}/figures.json"))
    out = {
        "fig1": fig_scaling(fig),
        "fig2": fig_walks(fig),
        "fig3": fig_stall(fig),
        "fig5": fig_variants(fig),
    }
    for name, svg in out.items():
        with open(f"{PERF}/{name}.svg.html", "w") as f:
            f.write(svg)
        print(name, len(svg), "bytes")


if __name__ == "__main__":
    main()
