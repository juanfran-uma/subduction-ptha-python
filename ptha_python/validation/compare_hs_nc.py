"""HTML report comparing two PTHA18-style HS scenario tables (netCDF).

Reference = PTHA18's published all_stochastic_slip_earthquake_events_<zone>.nc,
ours = the one step 7d writes. The report keeps four questions apart:
  A  format       is every variable the same name, order, type, dimension, text convention?
  B  invariants   PTHA18's own consistency checks (check_event_netcdf_files.R), on both files
  C  exactness    the engine on PTHA18's own inputs vs the published numbers
  D  distributions how the two tables differ in value, and why (different mesh, catalogue, random draws)

    python ptha_python/validation/compare_hs_nc.py

Any zone: --zone <name>. It needs the zone's reference table in nc_files/, the json/npz that
validate_weight_with_nonzero_rate.py writes in html/data/ (<zone>_validation_official.json and, for a
zone PTHA18 publishes with segments, <zone>_validation_plain.json), and the example folder of that name.
"""

import argparse
import html
import json
import os

import netCDF4
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
ZONE = "kermadectonga2"
EX = os.path.join(ROOT, "examples", "ptha18_compare", ZONE)
BLUE, ORANGE = "#2a78d6", "#eb6834"          # categorical slots 1 and 2 of the reference palette
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e4de", "#fcfcfb"

MEANING = {
    "max_nchar": "Dimension variable: width of the text matrices (characters of the longest string)",
    "table_rows": "Dimension variable: one entry per scenario (row of the table)",
    "event_index_string": "Unit sources that break, 1-based numbers of the unit-source table, joined by '-'",
    "event_slip_string": "Slip (m) on each of those unit sources, same order, joined by '_' (4 significant digits)",
    "Mw": "Moment magnitude with the constant rigidity (30 GPa), on the 7.2 to 9.8 grid of 0.1",
    "target_lon": "Longitude of the unit source the random slip field was centred on",
    "target_lat": "Latitude of that unit source",
    "peak_slip_downdip_ind": "Down-dip row of the cell where the slip field was asked to peak",
    "peak_slip_alongstrike_ind": "Along-strike column of that cell",
    "physical_corner_wavenumber_x": "Corner wavenumber (1/km, along strike) of the slip field's spectrum",
    "physical_corner_wavenumber_y": "Corner wavenumber (1/km, down dip) of the slip field's spectrum",
    "sourcename": "Source zone",
    "uniform_event_row": "Row of the uniform-slip (FAUS) scenario this slip field was drawn from, its parent",
    "rate_annual": "Annual rate of the scenario, constant rigidity (logic-tree mean)",
    "rate_annual_lower_ci": "Lower credible limit of that rate",
    "rate_annual_upper_ci": "Upper credible limit of that rate",
    "variable_mu_Mw": "Magnitude when each unit source uses the depth-dependent rigidity",
    "variable_mu_rate_annual": "Annual rate with the variable rigidity (logic-tree mean)",
    "variable_mu_rate_annual_lower_ci": "Lower credible limit, variable rigidity",
    "variable_mu_rate_annual_upper_ci": "Upper credible limit, variable rigidity",
    "variable_mu_rate_annual_median": "Median rate, variable rigidity",
    "variable_mu_rate_annual_16pc": "16th percentile of the rate, variable rigidity",
    "variable_mu_rate_annual_84pc": "84th percentile of the rate, variable rigidity",
    "variable_mu_weight_with_nonzero_rate": "Weight of the logic-tree branches giving this magnitude a non-zero rate, variable rigidity",
    "weight_with_nonzero_rate": "Weight of the logic-tree branches giving this magnitude a non-zero rate (0 = impossible)",
    "rate_annual_16pc": "16th percentile of the rate, constant rigidity",
    "rate_annual_84pc": "84th percentile of the rate, constant rigidity",
    "rate_annual_median": "Median rate, constant rigidity",
}
SOURCE = {
    "event_index_string": "step 7b event_index_string, unit sources sorted ascending",
    "event_slip_string": "step 7b event_slip_string, reordered with the indices",
    "Mw": "step 7b Mw, on R's seq(7.2, 9.8, 0.1) grid", "target_lon": "step 7b target_lon, mod 360",
    "target_lat": "step 7b target_lat", "peak_slip_downdip_ind": "step 7b", "peak_slip_alongstrike_ind": "step 7b",
    "physical_corner_wavenumber_x": "step 7b", "physical_corner_wavenumber_y": "step 7b",
    "sourcename": "the zone name", "uniform_event_row": "step 7b faus_event_id",
    "rate_annual": "step 7b rate_mean", "rate_annual_lower_ci": "step 7b rate_p0.025",
    "rate_annual_upper_ci": "step 7b rate_p0.975", "variable_mu_Mw": "step 7c variable_mu_Mw",
    "variable_mu_rate_annual": "step 7c variable_mu_rate_mean",
    "variable_mu_rate_annual_lower_ci": "step 7c variable_mu_rate_p0.025",
    "variable_mu_rate_annual_upper_ci": "step 7c variable_mu_rate_p0.975",
    "variable_mu_rate_annual_median": "step 7c variable_mu_rate_p0.5",
    "variable_mu_rate_annual_16pc": "step 7c variable_mu_rate_p0.16",
    "variable_mu_rate_annual_84pc": "step 7c variable_mu_rate_p0.84",
    "variable_mu_weight_with_nonzero_rate": "step 7 engine (variable-rigidity branch weights), via 7c",
    "weight_with_nonzero_rate": "step 7 engine (branch weights), via 7b",
    "rate_annual_16pc": "step 7b rate_p0.16", "rate_annual_84pc": "step 7b rate_p0.84",
    "rate_annual_median": "step 7b rate_p0.5",
}


# ---------------------------------------------------------------- data
def strings(ds, name):
    return np.array([str(s).strip("\x00 ") for s in netCDF4.chartostring(ds.variables[name][:])])


def load(path):
    ds = netCDF4.Dataset(path)
    ds.set_auto_mask(False)
    t = {n: np.asarray(ds[n][:]) for n, v in ds.variables.items() if v.dtype != "S1"}
    idx, slp = strings(ds, "event_index_string"), strings(ds, "event_slip_string")
    t["n_units"] = np.array([sum(1 for x in s.split("-") if x) for s in idx])
    t["peak_slip"] = np.array([max(float(x) for x in s.split("_") if x) for s in slp])
    t["_idx"], t["_slp"] = idx, slp
    t["mw_bin"] = np.round(t["Mw"], 1)
    t["_ds"] = ds
    return t


def per_mw(t):
    out = {}
    for m in np.unique(t["mw_bin"]):
        k = t["mw_bin"] == m
        pos = k & (t["rate_annual"] > 0)
        out[float(m)] = {
            "rows": int(k.sum()), "parents": int(len(np.unique(t["uniform_event_row"][k]))),
            "rate": float(t["rate_annual"][k].sum()), "rate_mu": float(t["variable_mu_rate_annual"][k].sum()),
            "w": float(t["weight_with_nonzero_rate"][k].mean()),
            "w_nz": float(t["weight_with_nonzero_rate"][k & (t["weight_with_nonzero_rate"] > 0)].mean()) if (k & (t["weight_with_nonzero_rate"] > 0)).any() else 0.0,
            "w_nz_min": float(t["weight_with_nonzero_rate"][k & (t["weight_with_nonzero_rate"] > 0)].min()) if (k & (t["weight_with_nonzero_rate"] > 0)).any() else 0.0,
            "w_nz_max": float(t["weight_with_nonzero_rate"][k].max()),
            "wmu": float(t["variable_mu_weight_with_nonzero_rate"][k].mean()),
            "dmu": float((t["variable_mu_Mw"][k] - t["Mw"][k]).mean()),
            "dmu_min": float(np.percentile(t["variable_mu_Mw"][k] - t["Mw"][k], 5)),
            "dmu_max": float(np.percentile(t["variable_mu_Mw"][k] - t["Mw"][k], 95)),
            "ps": [float(x) for x in np.percentile(t["peak_slip"][k], [5, 50, 95])],
            "nu": [float(x) for x in np.percentile(t["n_units"][k], [50, 95])],
            "kcx": float(np.median(t["physical_corner_wavenumber_x"][k])),
            "ci": float(np.median(t["rate_annual_upper_ci"][pos] / t["rate_annual"][pos])) if pos.any() else float("nan"),
        }
    return out


# ---------------------------------------------------------------- svg helpers
def nice_ticks(lo, hi, n=5):
    span = hi - lo
    step = 10 ** np.floor(np.log10(span / n))
    for f in (1, 2, 2.5, 5, 10):
        if span / (step * f) <= n:
            step *= f
            break
    a = np.ceil(lo / step) * step
    return [a + i * step for i in range(int((hi - a) / step + 1e-9) + 1)]


def fnum(v):
    if v == 0:
        return "0"
    if abs(v) >= 1000 or abs(v) < 0.01:
        return f"{v:.0e}".replace("e-0", "e-").replace("e+0", "e")
    return f"{v:g}"


def line_chart(series, xlabel, ylabel, ylog=False, w=540, h=300, ylim=None, hlines=(), note=""):
    L, R, T, B = 64, 14, 12, 46
    xs = np.concatenate([s["x"] for s in series])
    ally = np.concatenate([np.concatenate([s["y"]] + ([s["lo"], s["hi"]] if "lo" in s else [])) for s in series])
    if ylog:
        ally = ally[ally > 0]
    lo, hi = (ylim if ylim else (ally.min(), ally.max()))
    if ylog:
        lo, hi = 10 ** np.floor(np.log10(lo)), 10 ** np.ceil(np.log10(hi))
        f = lambda v: np.log10(np.maximum(v, lo))
        a, b = np.log10(lo), np.log10(hi)
        yt = [10.0 ** k for k in range(int(a), int(b) + 1)]
    else:
        pad = 0 if ylim else (hi - lo) * 0.06
        lo, hi = lo - pad, hi + pad
        f = lambda v: v
        a, b = lo, hi
        yt = nice_ticks(lo, hi)
    x0, x1 = xs.min(), xs.max()
    px = lambda v: L + (v - x0) / (x1 - x0) * (w - L - R)
    py = lambda v: T + (1 - (f(v) - a) / (b - a)) * (h - T - B)
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{html.escape(ylabel)} by {html.escape(xlabel)}">']
    for v in yt:
        o.append(f'<line x1="{L}" x2="{w-R}" y1="{py(v):.1f}" y2="{py(v):.1f}" stroke="{GRID}"/>'
                 f'<text x="{L-6}" y="{py(v)+4:.1f}" text-anchor="end" class="ax">{fnum(v)}</text>')
    for v in nice_ticks(x0, x1, 7):
        o.append(f'<text x="{px(v):.1f}" y="{h-B+16}" text-anchor="middle" class="ax">{fnum(v)}</text>')
    for v, lab in hlines:
        o.append(f'<line x1="{L}" x2="{w-R}" y1="{py(v):.1f}" y2="{py(v):.1f}" stroke="{INK2}" stroke-dasharray="4 4"/>'
                 f'<text x="{w-R-4}" y="{py(v)-4:.1f}" text-anchor="end" class="ax">{html.escape(lab)}</text>')
    o.append(f'<text x="{(L+w-R)/2}" y="{h-6}" text-anchor="middle" class="axl">{html.escape(xlabel)}</text>')
    o.append(f'<text transform="translate(14 {(T+h-B)/2}) rotate(-90)" text-anchor="middle" class="axl">{html.escape(ylabel)}</text>')
    for s in series:
        c = s["color"]
        if "lo" in s:
            up = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in zip(s["x"], s["hi"]))
            dn = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in list(zip(s["x"], s["lo"]))[::-1])
            o.append(f'<polygon points="{up} {dn}" fill="{c}" opacity="0.15"/>')
        pts = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in zip(s["x"], s["y"]) if (y > 0 or not ylog))
        o.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="2" stroke-linejoin="round"/>')
        for x, y in zip(s["x"], s["y"]):
            if y > 0 or not ylog:
                o.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="3" fill="{c}" stroke="{SURF}" stroke-width="1.5">'
                         f'<title>{html.escape(s["name"])}, {html.escape(xlabel)} {fnum(x)}: {y:.4g}</title></circle>')
    o.append("</svg>")
    legend = "".join(f'<span><i style="background:{s["color"]}"></i>{html.escape(s["name"])}</span>' for s in series)
    return f'<figure class="chart"><div class="legend">{legend}</div>{"".join(o)}<figcaption>{note}</figcaption></figure>'


def scatter_map(sets, w=540, h=380):
    L, R, T, B = 52, 12, 10, 40
    lon = np.concatenate([s["lon"] for s in sets])
    lat = np.concatenate([s["lat"] for s in sets])
    x0, x1, y0, y1 = lon.min() - 0.5, lon.max() + 0.5, lat.min() - 1, lat.max() + 1
    px = lambda v: L + (v - x0) / (x1 - x0) * (w - L - R)
    py = lambda v: T + (1 - (v - y0) / (y1 - y0)) * (h - T - B)
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="unit-source centres of both meshes">']
    for v in nice_ticks(y0, y1, 6):
        o.append(f'<line x1="{L}" x2="{w-R}" y1="{py(v):.1f}" y2="{py(v):.1f}" stroke="{GRID}"/>'
                 f'<text x="{L-6}" y="{py(v)+4:.1f}" text-anchor="end" class="ax">{fnum(v)}</text>')
    for v in nice_ticks(x0, x1, 6):
        o.append(f'<text x="{px(v):.1f}" y="{h-B+16}" text-anchor="middle" class="ax">{fnum(v)}</text>')
    o.append(f'<text x="{(L+w-R)/2}" y="{h-6}" text-anchor="middle" class="axl">longitude (degrees east)</text>')
    o.append(f'<text transform="translate(14 {(T+h-B)/2}) rotate(-90)" text-anchor="middle" class="axl">latitude</text>')
    for s in sets:
        for x, y in zip(s["lon"], s["lat"]):
            if s["hollow"]:
                o.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="4.5" fill="none" stroke="{s["color"]}" stroke-width="1.6"/>')
            else:
                o.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="2.6" fill="{s["color"]}"/>')
    o.append("</svg>")
    legend = "".join(f'<span><i style="background:{s["color"]}"></i>{html.escape(s["name"])}</span>' for s in sets)
    return f'<figure class="chart"><div class="legend">{legend}</div>{"".join(o)}</figure>'


def table(head, rows, cls=""):
    th = "".join(f"<th>{h}</th>" for h in head)
    tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'


def badge(ok, yes="PASS", no="DIFF"):
    return f'<span class="b {"ok" if ok else "no"}">{yes if ok else no}</span>'


# ---------------------------------------------------------------- diagrams
def pipeline_diagram(c):
    def box(x, y, w, h, col, t1, t2):
        dash = " stroke-dasharray='7 5'" if col == BLUE else ""
        return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{col}" fill-opacity="0.12" stroke="{col}" stroke-width="2"{dash}/>'
                f'<text x="{x+w/2}" y="{y+h/2-3}" text-anchor="middle" class="bt">{t1}</text>'
                f'<text x="{x+w/2}" y="{y+h/2+14}" text-anchor="middle" class="bs">{t2}</text>')
    def arrow(x1, y1, x2, y2):
        return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{INK2}" stroke-width="2" marker-end="url(#ar)"/>'
    o = ['<svg viewBox="0 0 980 296" role="img" aria-label="How both tables are made">',
         f'<defs><marker id="ar" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
         f'<path d="M0 0L10 5L0 10z" fill="{INK2}"/></marker></defs>',
         f'<text x="10" y="22" class="lane" fill="{BLUE}">nc from the Drive (the PTHA18 table): only READ, not re-run here</text>',
         f'<text x="10" y="192" class="lane" fill="{ORANGE}">This pipeline (Python): SLAB2.0 from scratch, --ptha false</text>']
    xs = [10, 205, 400, 595, 790]
    top = [("Unit-source table", f"{c['units_D']} unit sources"), ("Uniform-slip scenarios", f"FAUS, {c['par_D']:,} ruptures"),
           ("Logic tree rates", "constant and variable rigidity"), ("Random-slip scenarios", "HS, 15+ per FAUS rupture"),
           ("HS table (.nc)", f"{c['rows_D']:,} rows")]
    bot = [("Steps 1-2: SLAB2.0 mesh", f"{c['units_O']} unit sources"), ("Step 7 ruptures", f"FAUS, {c['par_O']:,} ruptures"),
           ("Step 7 + 7c rates", "same recipe, validated"), ("Step 7b slip fields", "HS, 15+ per FAUS rupture"),
           ("Step 7d export", f"{c['rows_O']:,} rows")]
    pills = [("differs: other mesh", 0), ("differs: other mesh", 0),
             ("same recipe: checked" if c["recipe_ok"] else "recipe DIFFERS", 1 if c["recipe_ok"] else 0),
             ("differs: random draws", 0), ("same format" if c["fmt_ok"] else "format DIFFERS", 1 if c["fmt_ok"] else 0)]
    for i, ((a, b), (cc, d)) in enumerate(zip(top, bot)):
        o.append(box(xs[i], 36, 180, 70, BLUE, a, b))
        o.append(box(xs[i], 206, 180, 70, ORANGE, cc, d))
        if i < 4:
            o.append(arrow(xs[i] + 182, 71, xs[i + 1] - 2, 71))
            o.append(arrow(xs[i] + 182, 241, xs[i + 1] - 2, 241))
    for i, (txt, ok) in enumerate(pills):
        col = "#1a7f37" if ok else "#b45309"
        o.append(f'<line x1="{xs[i]+90}" y1="108" x2="{xs[i]+90}" y2="204" stroke="{GRID}" stroke-width="2"/>')
        o.append(f'<rect x="{xs[i]+10}" y="128" width="160" height="30" rx="15" fill="{SURF}" stroke="{col}"/>'
                 f'<rect x="{xs[i]+10}" y="128" width="160" height="30" rx="15" fill="{col}" fill-opacity="0.12"/>'
                 f'<text x="{xs[i]+90}" y="148" text-anchor="middle" class="pill" fill="{col}">{txt}</text>')
    o.append("</svg>")
    return "".join(o)


ZOOM = """
<div id="zoom" hidden><div class="zbar"><b>Zoom</b><button id="zreset">Reset</button><button id="zclose">Close</button><span>scroll: zoom, drag: move, Esc: close</span></div><div id="zbody"></div></div>
<script>
(function () {
  var z = document.getElementById('zoom'), body = document.getElementById('zbody'), svg, vb, base, drag = null;
  function show() { svg.setAttribute('viewBox', vb.join(' ')); }
  function openFig(fig) {
    var s = fig.querySelector('svg'); if (!s) return;
    body.innerHTML = '';
    var leg = fig.querySelector('.legend'); if (leg) body.appendChild(leg.cloneNode(true));
    svg = s.cloneNode(true); body.appendChild(svg);
    base = svg.getAttribute('viewBox').split(' ').map(Number); vb = base.slice();
    z.hidden = false; document.body.style.overflow = 'hidden';
    svg.addEventListener('wheel', function (e) {
      e.preventDefault();
      var r = svg.getBoundingClientRect(), f = e.deltaY < 0 ? 0.85 : 1 / 0.85;
      var nw = Math.min(base[2], Math.max(base[2] / 30, vb[2] * f)), k = nw / vb[2];
      var px = vb[0] + (e.clientX - r.left) / r.width * vb[2], py = vb[1] + (e.clientY - r.top) / r.height * vb[3];
      vb = [px - (px - vb[0]) * k, py - (py - vb[1]) * k, nw, vb[3] * k]; show();
    }, { passive: false });
    svg.addEventListener('pointerdown', function (e) { drag = { x: e.clientX, y: e.clientY, vb: vb.slice() }; svg.setPointerCapture(e.pointerId); });
    svg.addEventListener('pointermove', function (e) {
      if (!drag) return; var r = svg.getBoundingClientRect();
      vb = [drag.vb[0] - (e.clientX - drag.x) / r.width * drag.vb[2], drag.vb[1] - (e.clientY - drag.y) / r.height * drag.vb[3], drag.vb[2], drag.vb[3]]; show();
    });
    svg.addEventListener('pointerup', function () { drag = null; });
  }
  function close() { z.hidden = true; body.innerHTML = ''; document.body.style.overflow = ''; }
  document.querySelectorAll('figure.chart').forEach(function (f) {
    f.tabIndex = 0; f.setAttribute('role', 'button'); f.title = 'Click to enlarge';
    f.addEventListener('click', function (e) { if (!e.target.closest('.legend')) openFig(f); });
    f.addEventListener('keydown', function (e) { if (e.key === 'Enter') openFig(f); });
  });
  document.getElementById('zclose').onclick = close;
  document.getElementById('zreset').onclick = function () { vb = base.slice(); show(); };
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
})();
</script>
"""


# ---------------------------------------------------------------- report
def main():
    global ZONE, EX
    _pre = argparse.ArgumentParser(add_help=False)
    _pre.add_argument("--zone", default=ZONE)
    ZONE = _pre.parse_known_args()[0].zone
    EX = os.path.join(ROOT, "examples", "ptha18_compare", ZONE)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zone", default=ZONE, help="source zone (default kermadectonga2)")
    ap.add_argument("--reference", default=os.path.join(ROOT, "nc_files", f"all_stochastic_slip_earthquake_events_{ZONE}.nc"))
    ap.add_argument("--ours", default=os.path.join(EX, "outputs", "nc", f"all_stochastic_slip_earthquake_events_{ZONE}.nc"))
    ap.add_argument("--validation-dir", default=os.path.join(ROOT, "html", "data"))
    ap.add_argument("--official-mesh", default=os.path.join(EX, "outputs_official", "ptha18_reference", f"unit_source_statistics_{ZONE}.nc"))
    ap.add_argument("--our-faus", default=os.path.join(EX, "outputs", f"scenario_rates_{ZONE}_source_zone.csv"),
                    help="our parent (uniform-slip) rates of the zone: the 0.5/0.5 mix on a segmented zone")
    ap.add_argument("--our-unseg-faus", default=os.path.join(EX, "outputs", f"scenario_rates_{ZONE}.csv"),
                    help="our unsegmented model alone (the first half of the mix)")
    ap.add_argument("--our-input", default=os.path.join(EX, "inputs", f"input_{ZONE}_scratch.json"))
    ap.add_argument("--our-mesh", default=os.path.join(EX, "outputs", "exports", "mesh_cells.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "html", f"{ZONE}_nc_comparison.html"))
    a = ap.parse_args()

    D, O = load(a.reference), load(a.ours)
    dD, dO = D["_ds"], O["_ds"]
    pD, pO = per_mw(D), per_mw(O)
    mws = sorted(pD)
    # ---- A: format
    names = list(dD.variables)
    same_order = names == list(dO.variables)
    same_type = all(dD[n].dtype == dO[n].dtype and dD[n].dimensions == dO[n].dimensions for n in names if n in dO.variables)
    mw_grid_same = np.array_equal(np.unique(D["Mw"]), np.unique(O["Mw"]))
    asc = lambda s: (lambda v: all(y > x for x, y in zip(v, v[1:])))([int(x) for x in s.split("-") if x])
    asc_D, asc_O = int(sum(map(asc, D["_idx"]))), int(sum(map(asc, O["_idx"])))
    var_rows = []
    fmt_ok_list = []
    for n in names:
        x = dD[n]
        kind = "text (char matrix)" if x.dtype == "S1" else {"float64": "float64", "int32": "int32"}[str(x.dtype)]
        if x.dtype == "S1" or n in ("max_nchar", "table_rows"):
            rng = ("", "")
        else:
            rng = (fnum(float(D[n].min())) + " to " + fnum(float(D[n].max())),
                   fnum(float(O[n].min())) + " to " + fnum(float(O[n].max())) if n in O else "missing")
        ok = n in dO.variables and dD[n].dtype == dO[n].dtype and dD[n].dimensions == dO[n].dimensions
        fmt_ok_list.append(ok)
        var_rows.append([f"<code>{n}</code>", kind, html.escape(MEANING.get(n, "")), rng[0], rng[1],
                         html.escape(SOURCE.get(n, "")), badge(ok, "same", "DIFF")])
    def sig_digits(tok):
        mant = tok.lower().split("e")[0].replace(".", "").lstrip("0")
        return len(mant.rstrip("0")) if mant else 0
    sig_max = lambda t: max(sig_digits(x) for sl in t["_slp"] for x in sl.split("_") if x)
    sep_ok = lambda t: bool(all(x.endswith("-") for x in t["_idx"]) and all(x.endswith("_") for x in t["_slp"]))
    sD, sO = sig_max(D), sig_max(O)
    _uq = np.unique(D["Mw"])
    _noisy = next((repr(float(x_)) for x_ in _uq if len(repr(float(x_))) > 5), None)
    grid_txt = f"{_uq[0]:g}, {_uq[1]:g}, ..., " + (f"{_noisy}, ..., " if _noisy else "") + f"{_uq[-1]:g}"
    attr_rows = [
        ["File format", dD.file_format, dO.file_format, badge(dD.file_format == dO.file_format)],
        ["table_rows dimension", "unlimited" if dD.dimensions["table_rows"].isunlimited() else "fixed",
         "unlimited" if dO.dimensions["table_rows"].isunlimited() else "fixed",
         badge(dD.dimensions["table_rows"].isunlimited() == dO.dimensions["table_rows"].isunlimited())],
        ["Variables (count, same set)", f"{len(names)}", f"{len(dO.variables)}", badge(set(names) == set(dO.variables) and len(names) == len(dO.variables))],
        ["Order of the variables", "this file's order", "same order" if same_order else "different order", badge(True, "same", "") if same_order else '<span class="b mid">INFO</span>'],
        ["Names, types and dimensions of every variable", "-", "all identical" if same_type else "differences", badge(same_type)],
        [f"Magnitude grid ({len(np.unique(D['Mw']))} values, float64 bit by bit)", grid_txt, "bit-identical" if mw_grid_same else "different", badge(mw_grid_same)],
        ["Unit sources of a scenario listed ascending", f"{asc_D:,} of {len(D['_idx']):,} rows", f"{asc_O:,} of {len(O['_idx']):,} rows", badge(asc_D == len(D['_idx']) and asc_O == len(O['_idx']))],
        ["Text separators ('-' indexes, '_' slips, trailing)", "yes" if sep_ok(D) else "no", "yes" if sep_ok(O) else "no", badge(sep_ok(D) and sep_ok(O))],
        ["Slip text precision (most significant digits in any slip)", f"{sD}", f"{sO}", badge(sD <= 4 and sO <= 4)],
        ["Numeric variables", "float64 and int32", "float64 and int32 (the type of every variable is compared above)", badge(same_type)],
    ]
    trunc = lambda v: html.escape(str(v) if len(str(v)) < 90 else str(v)[:87] + "...")
    attr_info = {
        "R_session_info": "Version of R and of its packages that wrote the file",
        "python_session_info": "Version of Python and of its packages that wrote ours (replaces R_session_info)",
        "unit_source_statistics_file": "Table that gives the unit-source numbers a position on the fault",
        "corresponding_uniform_slip_events_file": "File with the uniform-slip (parent) scenarios",
        "discretized_sources_file": "The unit-source mesh the scenarios were built on",
        "parent_script_name": "Script that wrote the file",
        "slip_type": "Kind of scenarios in the table", "shear_modulus": "Constant rigidity (Pa) of the Mw column",
        "source_summary_file": "Our summary.csv the table was reshaped from (extra)",
    }
    glob_rows = [[f"<code>{k}</code>", html.escape(attr_info.get(k, "")), trunc(dD.getncattr(k)) if k in dD.ncattrs() else "-",
                  trunc(dO.getncattr(k)) if k in dO.ncattrs() else "-"] for k in sorted(set(dD.ncattrs()) | set(dO.ncattrs()))]
    import glob
    other_orders = []
    for _f in sorted(glob.glob(os.path.join(os.path.dirname(a.reference), "all_stochastic_slip_earthquake_events_*.nc"))):
        if os.path.abspath(_f) != os.path.abspath(a.reference):
            with netCDF4.Dataset(_f) as _d:
                other_orders.append((os.path.basename(_f).replace("all_stochastic_slip_earthquake_events_", "").replace(".nc", ""),
                                     list(_d.variables) != names and set(_d.variables) == set(names)))
    differing = [z for z, d in other_orders if d]
    order_note_html = ((f'<div class="note">The order of the variables is not part of the format: PTHA18\'s own files do not agree on it. '
                       f'The files of {", ".join(differing)} have the same variables as this one but in a different order. Readers find a variable by its name.</div>')
                       if differing and not same_order else
                       ('<div class="note">The order of the variables differs from ours. It is not part of the format: readers find a variable by its name.</div>' if not same_order else ""))
    # ---- B: invariants
    def checks(t, label):
        w, r, mw = t["weight_with_nonzero_rate"], t["rate_annual"], t["Mw"]
        wm, rm = t["variable_mu_weight_with_nonzero_rate"], t["variable_mu_rate_annual"]
        k_hi, k_lo = np.round(mw, 3) > 9.6, (np.round(mw, 3) < 7.35) & (r > 0)
        u = t["uniform_event_row"]
        cnt = np.bincount(u)[1:]
        return {
            "weight == 0 implies rate == 0 (constant)": bool(np.all(r[w == 0] == 0)),
            "weight == 0 implies rate == 0 (variable)": bool(np.all(rm[wm == 0] == 0)),
            "Mw above 9.6 (thrust limit) has weight 0": bool(np.all(w[k_hi] == 0)) and bool(np.all(wm[k_hi] == 0)),
            "Mw below 7.35 with a rate has weight 1": bool(np.all(np.abs(w[k_lo] - 1) < 1e-6)),
            "|variable_mu_Mw - Mw| never above 0.4": bool(np.all(np.abs(t["variable_mu_Mw"] - mw) <= 0.4)),
            "all rates and weights are >= 0": bool(min(r.min(), rm.min(), w.min(), wm.min()) >= 0),
            "weights never above 1": bool(max(w.max(), wm.max()) <= 1 + 1e-9),
            "uniform_event_row is sorted and has no gaps": bool(np.all(np.diff(u) >= 0) and len(np.unique(u)) == u.max()),
            "every parent has at least 15 slip fields": bool(cnt.min() >= 15),
            "index and slip strings have the same number of unit sources": all(
                len([x for x in a_.split("-") if x]) == len([x for x in b_.split("_") if x]) for a_, b_ in zip(t["_idx"], t["_slp"])),
        }
    cD, cO = checks(D, "D"), checks(O, "O")
    chk_rows = [[k, badge(cD[k]), badge(cO[k])] for k in cD]
    # soft stats
    def ci_stats(t):
        r, lo, hi = t["rate_annual"], t["rate_annual_lower_ci"], t["rate_annual_upper_ci"]
        pos = r > 0
        inside = float(np.mean(((r >= lo) & (r <= hi))[pos]))
        order = np.vstack([t[n] for n in ("rate_annual_lower_ci", "rate_annual_16pc", "rate_annual_median", "rate_annual_84pc", "rate_annual_upper_ci")])
        viol = float(np.mean((np.diff(order, axis=0) < -1e-18).any(axis=0)))
        return inside, viol
    iD, vD = ci_stats(D)
    iO, vO = ci_stats(O)
    # ---- C: exactness
    val_dir = a.validation_dir
    vs = json.load(open(os.path.join(val_dir, f"{ZONE}_validation_official.json")))
    segn = np.load(os.path.join(val_dir, f"{ZONE}_validation_official.npz"))
    # how PTHA18 publishes this zone: segments in sourcezone_parameters.csv mean a mix with a segmented model
    szp = os.path.join(ROOT, "rptha", "R", "examples", "austptha_template", "DATA", "SOURCEZONE_PARAMETERS", "sourcezone_parameters.csv")
    pt_segs, pt_seg_names = [], []
    if os.path.exists(szp):
        import csv as _csv
        with open(szp, newline="") as _fh:
            for _r in _csv.DictReader(_fh):
                _r = {k.strip(): (v or "").strip() for k, v in _r.items() if k}
                if _r.get("sourcename") == ZONE and _r.get("segment_name"):
                    pt_segs.append(f"{_r['segment_boundary_alongstrike_index_lower']}-{_r['segment_boundary_alongstrike_index_upper']}")
                    pt_seg_names.append(_r["segment_name"].strip("_").capitalize())
    ptha_seg = bool(pt_segs)
    _plain = os.path.join(val_dir, f"{ZONE}_validation_plain.npz")
    unsn = np.load(_plain) if (ptha_seg and os.path.exists(_plain)) else None
    n_par = int(vs.get('n_parents', vs.get('n_parents_engine')))
    n_zero = int((segn['w_pub'] == 0).sum())
    _first = {}
    for _i, _k in enumerate(D["uniform_event_row"]):
        _first.setdefault(int(_k), _i)
    n_zero_all = int(sum(D["weight_with_nonzero_rate"][_i] == 0 for _i in _first.values()))
    _zr = vs["rate_mean (constant)"].get("engine_nonzero_where_published_zero")
    zero_rate_txt = (f" Where PTHA18's parent rate is 0 ({vs['rate_mean (constant)']['n_published_zero']} parent{'s' if vs['rate_mean (constant)']['n_published_zero'] != 1 else ''}) the engine's is 0 too." if _zr == 0 else
                     (f" The engine gives a non-zero rate to {_zr} parent{'s' if _zr != 1 else ''} whose published rate is 0." if _zr else ""))
    _zm = vs["weight_with_nonzero_rate"].get("zero_mismatch")
    zero_w_txt = f" The engine puts {"it" if n_zero == 1 else "them"} at 0 as well." if _zm == 0 else f" The engine disagrees with PTHA18 on {_zm} zero-weight patterns."
    beyond = vs.get("published_parents_beyond_engine_range", {})
    beyond_txt = (f" The published table also has {len(beyond['Mw'])} parent{'s' if len(beyond['Mw']) != 1 else ''} beyond the engine's magnitude range (Mw {', '.join(f'{m:.1f}' for m in beyond['Mw'])}), {'all with weight 0' if all(w_ == 0 for w_ in beyond['weight_with_nonzero_rate']) else 'with non-zero weights'}, which are not part of this count."
                  if beyond.get("Mw") else "")
    cval = []
    for lab, key, unit in (("weight_with_nonzero_rate", "weight_with_nonzero_rate", "abs"),
                           ("variable_mu_weight_with_nonzero_rate", "variable_mu_weight_with_nonzero_rate", "abs"),
                           ("rate_annual of the parent (constant rigidity)", "rate_mean (constant)", "rel"),
                           ("rate_annual of the parent (variable rigidity)", "variable_mu rate_mean", "rel")):
        r = vs[key]
        n = r.get("n", r.get("n_compared"))
        mx = r["max_abs_diff"] if unit == "abs" else r["max_rel_diff"]
        cval.append([f"<code>{lab}</code>", f"{r['n_within_1e-9']:,} of {n:,}", f"{mx:.1e} ({'absolute' if unit == 'abs' else 'relative'})",
                     badge(r["n_within_1e-9"] == n)])
    mwu = np.round(segn["mw"], 1)
    rng_txt = f"Mw {mwu.min():.1f} to {mwu.max():.1f}"
    match_upto = p96 = u96 = s96 = None
    unseg_rows = []
    if unsn is not None:
        mwu = np.round(unsn["mw"], 1)
        bins_sorted = sorted(set(mwu))
        match_upto = None
        for _m in bins_sorted:
            _k = (mwu == _m) & (unsn["w_pub"] > 0)
            if _k.any() and np.all(np.abs(unsn["w_eng"][_k] - unsn["w_pub"][_k]) < 1e-9):
                match_upto = _m
            else:
                break
        _k96 = mwu == 9.6
        p96 = float(unsn["w_pub"][_k96].mean()) if _k96.any() else float("nan")
        u96 = float(unsn["w_eng"][_k96].mean()) if _k96.any() else float("nan")
        s96 = 2 * p96 - u96
        def rng(arr):
            nz = arr[arr > 0]
            txt = f"{nz.min():.4f} to {nz.max():.4f}" if nz.size else "none above 0"
            z0 = int((arr == 0).sum())
            return txt + (f" ({z0} parents at 0)" if z0 else "")
        for m in (7.2, 8.0, 8.2, 8.6, 9.0, 9.4, 9.6):
            k = mwu == m
            unseg_rows.append([f"{m:.1f}", rng(unsn["w_pub"][k]), rng(unsn["w_eng"][k]), rng(segn["w_eng"][k])])
    # ---- D: charts    # ---- D: charts
    xs = np.array(mws)
    S = lambda key, p, scale=1.0: np.array([p[m][key] for m in mws]) * scale
    ex = lambda t: np.array([t["rate_annual"][t["Mw"] >= m - 1e-9].sum() for m in xs])
    exm = lambda t: np.array([t["variable_mu_rate_annual"][t["Mw"] >= m - 1e-9].sum() for m in xs])
    c1 = line_chart([{"name": "PTHA18", "color": BLUE, "x": xs, "y": S("rate", pD)}, {"name": "This pipeline", "color": ORANGE, "x": xs, "y": S("rate", pO)}],
                    "Mw (constant rigidity)", "sum of scenario rates (events per year)", ylog=True,
                    note="Rate of each magnitude bin, summed over the scenarios of the bin. Hover a point for its value.")
    c2 = line_chart([{"name": "PTHA18", "color": BLUE, "x": xs, "y": ex(D)}, {"name": "This pipeline", "color": ORANGE, "x": xs, "y": ex(O)}],
                    "Mw", "annual rate of Mw or larger", ylog=True,
                    note="Exceedance curve of the HS table, constant rigidity. It is the quantity a hazard study reads from the table.")
    c3 = line_chart([{"name": "PTHA18", "color": BLUE, "x": xs, "y": exm(D)}, {"name": "This pipeline", "color": ORANGE, "x": xs, "y": exm(O)}],
                    "Mw", "annual rate of Mw or larger", ylog=True,
                    note="The same curve with the variable rigidity (variable_mu_rate_annual, each scenario at its variable_mu_Mw is not re-binned: PTHA18 keeps the constant-rigidity bin).")
    zero_mws = [m for m in mws if pD[m]["w_nz"] == 0]
    zero_mw_txt = ("Mw " + ", ".join(f"{m:.1f}" for m in zero_mws)) if zero_mws else "none"
    c4_tail = ("The curves differ where the logic trees differ: both are the half-and-half mix of an unsegmented and a segmented model, with different segment definitions and parameters."
               if ptha_seg else "The curves differ where the logic trees differ (our Mw_max branches and parameters against PTHA18's).")
    mD_nz = np.array([pD[m]["w_nz"] > 0 for m in mws])
    mO_nz = np.array([pO[m]["w_nz"] > 0 for m in mws])
    c4 = line_chart([{"name": "PTHA18 (mean, min to max of the non-zero weights)", "color": BLUE, "x": xs[mD_nz], "y": S("w_nz", pD)[mD_nz],
                      "lo": S("w_nz_min", pD)[mD_nz], "hi": S("w_nz_max", pD)[mD_nz]},
                     {"name": "This pipeline (same)", "color": ORANGE, "x": xs[mO_nz], "y": S("w_nz", pO)[mO_nz],
                      "lo": S("w_nz_min", pO)[mO_nz], "hi": S("w_nz_max", pO)[mO_nz]}],
                    "Mw", "weight_with_nonzero_rate", ylim=(0, 1.02),
                    note=f"Share of the logic tree that allows a scenario of this magnitude. Parents with weight 0 are left out, and so are magnitudes where every parent has weight 0 ({zero_mw_txt} in PTHA18's table). PTHA18's table has {n_zero_all} parent{'s' if n_zero_all != 1 else ''} at weight 0 ({n_zero} inside the engine's range, see section C). {c4_tail}")
    import sys as _sys
    _sys.path.insert(0, os.path.join(ROOT, "ptha_python"))
    from pyptha import variable_mu as vmu
    with netCDF4.Dataset(a.official_mesh) as _m:
        _area = np.asarray(_m["length"][:], float) * np.asarray(_m["width"][:], float)
        _depth = np.asarray(_m["depth"][:], float)
    vmu_check = float(np.max(np.abs(vmu.variable_mu_Mw(list(D["_idx"]), list(D["_slp"]), _area, _depth) - D["variable_mu_Mw"])))
    c5 = line_chart([{"name": "PTHA18 (mean, 5 to 95%)", "color": BLUE, "x": xs, "y": S("dmu", pD), "lo": S("dmu_min", pD), "hi": S("dmu_max", pD)},
                     {"name": "This pipeline (mean, 5 to 95%)", "color": ORANGE, "x": xs, "y": S("dmu", pO), "lo": S("dmu_min", pO), "hi": S("dmu_max", pO)}],
                    "Mw (constant rigidity)", "variable_mu_Mw minus Mw", hlines=[(0, "no change")],
                    note=f"Soft shallow rock lowers the magnitude, stiff deep rock raises it. Overall range: PTHA18 {float((D['variable_mu_Mw'] - D['Mw']).min()):+.3f} to {float((D['variable_mu_Mw'] - D['Mw']).max()):+.3f}, ours {float((O['variable_mu_Mw'] - O['Mw']).min()):+.3f} to {float((O['variable_mu_Mw'] - O['Mw']).max()):+.3f}. Recomputing PTHA18's variable_mu_Mw from its own slip strings and its unit-source table with pyptha/variable_mu.py gives a largest difference of {vmu_check:.1e}.")
    ps = lambda p, i: np.array([p[m]["ps"][i] for m in mws])
    r50, r95 = ps(pO, 1) / ps(pD, 1), ps(pO, 2) / ps(pD, 2)
    ps_text = (f"Median peak slip, ours divided by PTHA18: {r50.min():.2f} to {r50.max():.2f}; 95th percentile: {r95.min():.2f} to {r95.max():.2f} (ours is below PTHA18's at the 95th percentile in {int((r95 < 1).sum())} of {len(r95)} magnitude bins). The slip fields are random draws, so this is a distribution check, not a row-by-row one.")
    c6 = line_chart([{"name": "PTHA18 (median, 5 to 95%)", "color": BLUE, "x": xs, "y": ps(pD, 1), "lo": ps(pD, 0), "hi": ps(pD, 2)},
                     {"name": "This pipeline (median, 5 to 95%)", "color": ORANGE, "x": xs, "y": ps(pO, 1), "lo": ps(pO, 0), "hi": ps(pO, 2)}],
                    "Mw", "peak slip of the scenario (m)", ylog=True,
                    note=ps_text)
    nu = lambda p, i: np.array([p[m]["nu"][i] for m in mws])
    c7 = line_chart([{"name": "PTHA18 (median, 95th percentile)", "color": BLUE, "x": xs, "y": nu(pD, 0), "lo": nu(pD, 0), "hi": nu(pD, 1)},
                     {"name": "This pipeline (median, 95th percentile)", "color": ORANGE, "x": xs, "y": nu(pO, 0), "lo": nu(pO, 0), "hi": nu(pO, 1)}],
                    "Mw", "unit sources in the scenario", ylog=True,
                    note="Size of the rupture in unit sources. It jumps where the rule that sets the number of cells per magnitude changes; the jump falls at slightly different magnitudes in the two tables because the cells differ.")
    kc = line_chart([{"name": "PTHA18 (median)", "color": BLUE, "x": xs, "y": S("kcx", pD)}, {"name": "This pipeline (median)", "color": ORANGE, "x": xs, "y": S("kcx", pO)}],
                    "Mw", "corner wavenumber along strike (1/km)", ylog=True,
                    note="Larger ruptures get smoother slip (smaller corner wavenumber), in both tables.")
    cw = line_chart([{"name": "PTHA18 (median)", "color": BLUE, "x": xs, "y": S("ci", pD)}, {"name": "This pipeline (median)", "color": ORANGE, "x": xs, "y": S("ci", pO)}],
                    "Mw", "upper_ci divided by rate_annual", ylog=True,
                    note="How far the upper credible limit sits above the mean rate. It grows with magnitude as the logic-tree branches disagree more (the lower limit is often 0, so the ratio to it is not used).")
    # ---- mesh
    import csv
    with netCDF4.Dataset(a.official_mesh) as m:
        olon, olat = np.mod(np.asarray(m["lon_c"][:], float), 360), np.asarray(m["lat_c"][:], float)
        oarea = float(np.sum(np.asarray(m["length"][:], float) * np.asarray(m["width"][:], float)))
        o_as, o_dd = int(np.max(m["alongstrike_number"][:])), int(np.max(m["downdip_number"][:]))
    with open(a.our_mesh, newline="") as fh:
        mr = list(csv.DictReader(fh))
    mlon, mlat = np.array([float(r["lon"]) for r in mr]) % 360, np.array([float(r["lat"]) for r in mr])
    m_as, m_dd = int(max(int(r["alongstrike_number"]) for r in mr)), int(max(int(r["downdip_number"]) for r in mr))
    mesh_chart = scatter_map([{"name": f"PTHA18 unit sources ({len(olon)})", "color": BLUE, "lon": olon, "lat": olat, "hollow": True},
                              {"name": f"This pipeline ({len(mr)})", "color": ORANGE, "lon": mlon, "lat": mlat, "hollow": False}])
    lat_span_o, lat_span_m = float(olat.max() - olat.min()), float(mlat.max() - mlat.min())
    mesh_extent_txt = (f"PTHA18's mesh spans latitudes {olat.min():.1f} to {olat.max():.1f} and longitudes {olon.min():.1f} to {olon.max():.1f}; "
                       f"ours spans {mlat.min():.1f} to {mlat.max():.1f} and {mlon.min():.1f} to {mlon.max():.1f}, "
                       f"{100 * lat_span_m / lat_span_o:.0f}% of PTHA18's latitude span and {100 * len(mr) / len(olon):.0f}% of its unit sources.")
    # latitude share of the rate
    bins = np.linspace(min(D["target_lat"].min(), O["target_lat"].min()), max(D["target_lat"].max(), O["target_lat"].max()), 31)
    cen = (bins[:-1] + bins[1:]) / 2
    shr = lambda t: np.histogram(t["target_lat"], bins=bins, weights=t["rate_annual"])[0] / t["rate_annual"].sum() * 100
    c_lat = line_chart([{"name": "PTHA18", "color": BLUE, "x": cen, "y": shr(D)}, {"name": "This pipeline", "color": ORANGE, "x": cen, "y": shr(O)}],
                       "target latitude of the scenario", "share of the total rate (%)",
                       note="Where along the trench the rate sits. The shape follows the plate convergence profile of each pipeline.")
    # per-Mw table
    pm_rows = []
    for m in mws:
        d, o = pD[m], pO[m]
        pm_rows.append([f"{m:.1f}", f"{d['parents']} / {o['parents']}", f"{d['rows']:,} / {o['rows']:,}", f"{d['rate']:.3e}", f"{o['rate']:.3e}",
                        f"{o['rate']/d['rate']:.2f}" if d["rate"] > 0 else "-", f"{d['w']:.3f}", f"{o['w']:.3f}"])
    tot = lambda t, k: float(t[k].sum())
    # ---- what kind of zone is this, and what did we run?
    import csv as _csv2
    ours_seg = os.path.exists(a.our_faus)
    our_faus_path = a.our_faus if ours_seg else a.our_unseg_faus
    run_flags = "--ptha false" + (" --segmented true" if ours_seg else "")
    def parent_rates(path):
        with open(path, newline="") as fh:
            fr = list(_csv2.DictReader(fh))
        return np.round([float(r["Mw"]) for r in fr], 1), np.array([float(r["rate_mean"]) for r in fr])
    ours_mw, ours_rate = parent_rates(our_faus_path)
    mwp = np.round(segn["mw"], 1)
    m96 = [m for m in mws if m in set(mwp)]
    pub_s = np.array([segn["r_pub"][mwp == m].sum() for m in m96])
    seg_s = np.array([ours_rate[ours_mw == m].sum() for m in m96])
    xm = np.array(m96)
    ix = lambda m: m96.index(m)
    gap_total = abs(1 - seg_s.sum() / pub_s.sum()) * 100
    with open(a.our_input, encoding="utf-8") as fh:
        our_segs = json.load(fh).get("segments", {})
    n_cols = int(max(int(r["alongstrike_number"]) for r in mr))
    ours_seg_txt = ", ".join(f"{v['alongstrike_slice'][0]}-{v['alongstrike_slice'][1]}" for v in our_segs.values()) or "none"
    pt_seg_txt = ", ".join(pt_segs) or "none"
    seg_names_txt = (", ".join(pt_seg_names[:-1]) + " and " + pt_seg_names[-1]) if len(pt_seg_names) > 1 else (pt_seg_names[0] if pt_seg_names else "")
    full_d1 = ptha_seg and ours_seg and unsn is not None
    gap_seg = gap_uns = share = 0.0
    if full_d1:
        # PTHA18's published rates; the engine on PTHA18's plain unsegmented input; our unsegmented half; our mix
        mwq = np.round(unsn["mw"], 1)
        off_s = np.array([unsn["r_eng"][mwq == m].sum() for m in m96])
        uns_mw, uns_rate = parent_rates(a.our_unseg_faus)
        uns_s = np.array([uns_rate[uns_mw == m].sum() for m in m96])
        okp = (pub_s > 0) & (off_s > 0)
        mm = [m for i_, m in enumerate(m96) if okp[i_]]
        sel = [m for m in (7.2, 7.6, 8.0, 8.4, 8.8, 9.0, 9.2, 9.4, 9.6) if m in mm]
        d1_rows = [[f"{m:.1f}", f"{pub_s[ix(m)]:.3e}", f"{off_s[ix(m)]:.3e}", f"{uns_s[ix(m)]:.3e}", f"{seg_s[ix(m)]:.3e}",
                    f"{uns_s[ix(m)]/pub_s[ix(m)]:.2f}", f"{seg_s[ix(m)]/pub_s[ix(m)]:.2f}"] for m in sel]
        d1_rows.append(["<b>total</b>", f"<b>{pub_s.sum():.4f}</b>", f"<b>{off_s.sum():.4f}</b>", f"<b>{uns_s.sum():.4f}</b>", f"<b>{seg_s.sum():.4f}</b>",
                        f"<b>{uns_s.sum()/pub_s.sum():.2f}</b>", f"<b>{seg_s.sum()/pub_s.sum():.2f}</b>"])
        gap_uns, gap_seg = abs(1 - uns_s.sum() / pub_s.sum()) * 100, abs(1 - seg_s.sum() / pub_s.sum()) * 100
        mesh_eff = abs(uns_s.sum() / off_s.sum() - 1) * 100
        hi_m = max(mm)
        hi_ratio = float(seg_s[ix(hi_m)] / pub_s[ix(hi_m)])
        closer = [abs(seg_s[ix(m)] / pub_s[ix(m)] - 1) < abs(uns_s[ix(m)] / pub_s[ix(m)] - 1) for m in mm]
        n_closer, n_bins = int(sum(closer)), len(closer)
        notc = [m for m, c in zip(mm, closer) if not c]
        not_closer_txt = (" (not closer at Mw " + ", ".join(f"{m:.1f}" for m in notc) + ")") if notc else ""
        others = [(m, float(uns_s[ix(m)] / off_s[ix(m)])) for m in (8.8, 9.4, 9.6) if m in mm]
        others_txt = ", ".join(f"{r_:.2f} at Mw {m:.1f}" for m, r_ in others) or "not available"
        share = (gap_uns - gap_seg) / gap_uns * 100 if gap_uns > 0 else 0.0
        low_m = [m for m in (7.2, 7.6, 8.0) if m in mm]
        lo_unseg = float(np.mean([uns_s[ix(m)] / pub_s[ix(m)] for m in low_m])) if low_m else float("nan")
        lo_seg = float(np.mean([seg_s[ix(m)] / pub_s[ix(m)] for m in low_m])) if low_m else float("nan")
        near = [abs(uns_s[ix(m)] / off_s[ix(m)] - 1) for m in mm if m <= 8.5]
        near_txt = f"{100 * max(near):.0f}%" if near else "n/a"
        allr = np.concatenate([off_s[okp] / pub_s[okp], uns_s[okp] / pub_s[okp], seg_s[okp] / pub_s[okp]])
        xk = np.array(mm)
        d1_chart = line_chart([{"name": "Engine on PTHA18's mesh, unsegmented", "color": INK2, "x": xk, "y": np.array([off_s[ix(m)] / pub_s[ix(m)] for m in mm])},
                               {"name": "Ours, own mesh, unsegmented half only", "color": "#1baf7a", "x": xk, "y": np.array([uns_s[ix(m)] / pub_s[ix(m)] for m in mm])},
                               {"name": "Ours, own mesh, segmented mix (this example)", "color": ORANGE, "x": xk, "y": np.array([seg_s[ix(m)] / pub_s[ix(m)] for m in mm])}],
                              "Mw", "parent rate divided by PTHA18's published rate", hlines=[(1.0, "PTHA18 = 1")],
                              ylim=(min(0.8, float(allr.min()) * 0.97), max(1.1, float(allr.max()) * 1.05)), w=1000, h=340,
                              note=f"The two unsegmented versions (gray and green) differ from each other by at most {near_txt} up to Mw 8.5: there the other inputs (PTHA18's mesh and parameters against ours) matter little. The segmented mix (orange) is closer to 1 than the unsegmented half in {n_closer} of {n_bins} magnitudes.")
        d1_html = f"""<h3>D.1 What explains the rates: the segmentation or the other inputs?</h3>
<p>Parent rates summed by magnitude, in four versions: PTHA18's published ones; the engine run on <b>PTHA18's own plain unsegmented input</b> (its mesh and its parameters); our own mesh, <b>unsegmented half only</b> (what the example's unsegmented model gives alone); and our own mesh with the <b>0.5 / 0.5 mix</b> of the unsegmented and segmented models, which is what the example's table holds (as PTHA18's does).</p>
{d1_chart}
{table(["Mw", "PTHA18 published", "Engine, PTHA18 mesh, unsegmented", "Ours, unsegmented half only", "Ours, segmented mix", "Unsegmented / PTHA18", "Segmented mix / PTHA18"], d1_rows)}
<p><b>Segmentation (same mesh, same parameters):</b> the unsegmented half alone is {gap_uns:.0f}% away from PTHA18's total and the mix is {gap_seg:.0f}% away, so the mix removes {share:.0f}% of the gap. At Mw 7.2 to 8.0 the ratio to PTHA18 goes from {lo_unseg:.2f} to {lo_seg:.2f}. The mix is closer to PTHA18 than the unsegmented half in {n_closer} of {n_bins} magnitude bins{not_closer_txt}. At the highest magnitudes ours is still higher (Mw {hi_m:.1f}: {hi_ratio:.2f} times). <b>All other inputs together (PTHA18's mesh and parameters against ours, both unsegmented):</b> the totals differ by {mesh_eff:.1f}%, but only because low magnitudes dominate the total; the ratio of ours to the engine on PTHA18's input is {others_txt}. These two inputs differ in more than the mesh (PTHA18's has its own conditional-probability model, b values, convergence and maximum magnitude), so this is not the effect of the mesh alone. Our segments are {ours_seg_txt} of {n_cols} columns (Berryman et al. 2015 end points on our mesh); PTHA18's are {pt_seg_txt} of {o_as} columns, so the boundaries are not directly comparable by index. Their parameters per segment also differ. This report does not separate those causes from each other.</p>"""
    else:
        okp = pub_s > 0
        ratio_s = seg_s[okp] / pub_s[okp]
        xr = xm[okp]
        y_top = 1.6
        clipped = [(m_, r_, pub_s[ix(m_)]) for m_, r_ in zip(xr, ratio_s) if r_ > y_top]
        clip_txt = ("" if not clipped else " Not drawn because they are above the plotted range: " + "; ".join(
            f"Mw {m_:.1f} ({r_:.2f} times, where PTHA18's rate is {pr_:.1e} events per year)" for m_, r_, pr_ in clipped) + ".")
        d1_chart = line_chart([{"name": "Ours / PTHA18", "color": ORANGE, "x": xr[ratio_s <= y_top], "y": ratio_s[ratio_s <= y_top]}],
                              "Mw", "parent rate divided by PTHA18's published rate", hlines=[(1.0, "PTHA18 = 1")],
                              ylim=(min(0.9, float(ratio_s.min()) * 0.97), max(1.1, min(float(ratio_s.max()), y_top) * 1.03)), w=1000, h=340,
                              note="Ours divided by PTHA18's, magnitude by magnitude (magnitudes where PTHA18's rate is 0 are left out). A value of 1 is a match." + clip_txt)
        sel_m = [m for m in (7.2, 7.6, 8.0, 8.4, 8.8, 9.0, 9.2, 9.4, 9.6) if m in m96 and pub_s[ix(m)] > 0]
        _worst_m = float(xr[int(np.argmax(np.abs(ratio_s - 1)))])
        if _worst_m not in sel_m:
            sel_m = sorted(sel_m + [_worst_m])
        d1_rows = [[f"{m:.1f}", f"{pub_s[ix(m)]:.3e}", f"{seg_s[ix(m)]:.3e}", f"{seg_s[ix(m)]/pub_s[ix(m)]:.2f}"] for m in sel_m]
        d1_rows.append(["<b>total</b>", f"<b>{pub_s.sum():.4f}</b>", f"<b>{seg_s.sum():.4f}</b>", f"<b>{seg_s.sum()/pub_s.sum():.2f}</b>"])
        lowm = [m for m in m96 if m <= 8.0 and pub_s[ix(m)] > 0]
        lo_r = [seg_s[ix(m)] / pub_s[ix(m)] for m in lowm]
        worst = int(np.argmax(np.abs(ratio_s - 1)))
        if ptha_seg:
            kind_txt = ("PTHA18 publishes this zone as a mix with a segmented model and our example is unsegmented, so part of the gap is the model itself "
                        "(on kermadectonga2 the segmented mix removed most of it).")
        else:
            kind_txt = "This zone has no segments in PTHA18, so both are one unsegmented model and the segmentation plays no part."
        d1_html = f"""<h3>D.1 How close are the rates?</h3>
<p>Parent (uniform-slip) rates summed by magnitude: PTHA18's published ones against ours. {kind_txt}</p>
{d1_chart}
{table(["Mw", "PTHA18 published", "Ours", "Ours / PTHA18"], d1_rows)}
<p>The totals are {seg_s.sum():.4f} (ours) and {pub_s.sum():.4f} (PTHA18), {gap_total:.0f}% apart. At Mw 7.2 to 8.0 ours is between {min(lo_r):.2f} and {max(lo_r):.2f} times PTHA18's; the largest difference is at Mw {xr[worst]:.1f} ({ratio_s[worst]:.2f} times). Differences come from our own mesh ({len(mr)} unit sources against PTHA18's {len(olon)}), our parameters (literature and the GCMT catalogue); this report does not separate them. The random slip fields of the HS table do not enter these parent rates.</p>"""
    ratio = {m: pO[m]["rate"] / pD[m]["rate"] for m in mws if pD[m]["rate"] > 0}
    cross = next((m for m in mws if m in ratio and ratio[m] > 1), None)
    tail_rate = ("Section D.1 separates the effect of the segmentation from that of the other inputs." if full_d1 else
                 "Both are the same kind of model, so these differences come from the inputs (mesh, parameters, catalogue), not from the model structure." if not ptha_seg else
                 "Section D.1 shows the gap to PTHA18's mix.")
    if cross is None:
        rate_text = "Our rates are lower than PTHA18's at every magnitude. " + tail_rate
    else:
        low = [ratio[m] for m in mws if m in ratio and m < cross] or [1.0]
        lo_txt = (f"Up to Mw {cross - 0.1:.1f} our rates are {abs(100*(1-max(low))):.0f} to {abs(100*(1-min(low))):.0f}% {'lower' if max(low) < 1 else 'different'} than PTHA18's; from Mw {cross:.1f} they are higher"
                  if cross > mws[0] else f"Our rates are higher than PTHA18's from Mw {cross:.1f}")
        mmax = max(ratio, key=ratio.get)
        rate_text = (lo_txt + f", up to {max(ratio.values()):.1f} times at Mw {mmax:.1f} (PTHA18's rate there is {pD[mmax]['rate']:.1e} events per year, which is {100 * pD[mmax]['rate'] / tot(D, 'rate_annual'):.2g}% of its total). " + tail_rate)
    # ---- the texts that depend on the zone kind
    total_card = ('<div class="card"><div class="k">Total annual rate</div><div class="v">' + f"{tot(O, 'rate_annual'):.3f} vs {tot(D, 'rate_annual'):.3f}" + '</div><div class="s">' + (
        f"ours vs PTHA18, {gap_seg:.0f}% apart. The unsegmented half alone would be {gap_uns:.0f}% apart: PTHA18 publishes a mix with a segmented model and so does this example (section D)" if full_d1 else
        f"ours vs PTHA18, {gap_total:.0f}% apart. This zone has no segments in PTHA18, so both are one unsegmented model (section D)" if not ptha_seg else
        f"ours vs PTHA18, {gap_total:.0f}% apart. PTHA18 publishes a mix with a segmented model and this example is unsegmented (section D)") + "</div></div>")
    max_w = max(vs["weight_with_nonzero_rate"]["max_abs_diff"], vs["variable_mu_weight_with_nonzero_rate"]["max_abs_diff"])
    max_r = max(vs["rate_mean (constant)"]["max_rel_diff"], vs["variable_mu rate_mean"]["max_rel_diff"])
    exact_txt = (f"Run on PTHA18's own inputs, our engine reproduces the published weights to {max_w:.0e} and the parent rates to {max_r:.0e} (percentiles are not compared). "
                 f"Our example, built from SLAB2.0 alone, ends {gap_total:.0f}% from PTHA18's total rate. ")
    tail_one = ("The remaining differences come from the inputs, not from the file: our mesh (the number of ruptures and rows), our own segment boundaries and parameters, and the random slip fields." if full_d1 else
                "This zone has no segments, so both are one model; the remaining differences come from the inputs, not from the file: our mesh (the number of ruptures and rows), our parameters and catalogue, and the random slip fields." if not ptha_seg else
                "PTHA18 publishes this zone as a mix with a segmented model and our example is unsegmented; the remaining differences come from that and from our mesh, parameters and the random slip fields.")
    one_sentence = '<div class="note good"><b>In one sentence.</b> The file has the same layout as PTHA18\'s. ' + exact_txt + tail_one + "</div>"
    how_html = (f"<p><b>What this report does.</b> (1) It reads PTHA18's published table; nothing of PTHA18 is recomputed from scratch (dashed blue boxes). (2) It builds our table from scratch with this pipeline (orange boxes): SLAB2.0 geometry, our own mesh and literature parameters (Berryman et al. 2015, Bird 2003), run with <code>{run_flags}</code>"
                + (" (segments at the trench end points of Berryman et al. 2015)" if ours_seg else "")
                + ", so steps 1 to 7 and 9 read no PTHA18 file. One input is PTHA18-derived: the table of peak-slip weights used by step 7b was recovered from PTHA18's published HS rates. (3) To test that our recipe is the same as PTHA18's, it also runs our engine on PTHA18's own input files and compares with the published numbers (section C). Both chains below have the same steps; what differs is what goes in at the start.</p>")
    d_note = (f'<div class="note"><b>What differs, and why.</b> This pipeline was run with <code>{run_flags}</code>: it builds its own mesh from SLAB2.0 ({len(mr)} unit sources of {m_as} along strike by {m_dd} down dip; PTHA18 has {len(olon)} of {o_as} by {o_dd}), '
              + ("places its segments at Berryman et al. (2015) end points, " if ours_seg else "") + "uses literature parameters, and draws its own random slip fields. "
              + (f"The mesh and the other inputs explain the counts (parents, rows) and a smaller part of the rates; the segmentation removes {share:.0f}% of the gap in the total rate: see D.1. " if full_d1 else
                 "The mesh explains the counts (parents, rows); the rates differ through our mesh, parameters and catalogue together, which this report does not separate: see D.1. ")
              + "Row-by-row equality is not expected because the slip fields are random; distributions are what to compare.</div>")
    c_intro = ((f"<p>Run the engine on PTHA18's own segmented input of {ZONE} (its official unit-source table, convergence, coupling and catalogue: unsegmented plus {seg_names_txt} segments, 0.5 / 0.5), with PTHA18's HS table as the variable-rigidity catalogue. Then compare, parent by parent, with what the published table says. Script: <code>ptha_python/validation/validate_weight_with_nonzero_rate.py</code>.</p>")
               if ptha_seg else
               (f"<p>Run the engine on PTHA18's own input of {ZONE} (its official unit-source table, Bird convergence, coupling and catalogue; this zone has no segments, so PTHA18 publishes it as one model), with PTHA18's HS table as the variable-rigidity catalogue. Then compare, parent by parent, with what the published table says. Script: <code>ptha_python/validation/validate_weight_with_nonzero_rate.py</code>.</p>"))
    if unsn is not None:
        match_txt = (f"the unsegmented run matches up to Mw {match_upto:.1f}, where every branch of both models allows the scenario"
                     if match_upto is not None else "the unsegmented run does not match at any magnitude")
        implied_txt = (f"at Mw 9.6 the published {p96:.3f} with an unsegmented half of {u96:.3f} implies a segmented half of about {s96:.3f} (both halves weighted 0.5)"
                       if p96 == p96 and p96 > 0 else "at the highest magnitudes the published weight is 0")
        c_seg_html = f"""<h3>Why the segmented run matches and the unsegmented one does not</h3>
<p>PTHA18 does not publish {ZONE} as one model. The published numbers are a <b>half-and-half mix</b>: 0.5 x an unsegmented model of the whole zone + 0.5 x a segmented model ({seg_names_txt}, each with its own convergence and maximum magnitude). <code>weight_with_nonzero_rate</code> is mixed with the same weights as the rates. Running the engine on a single unsegmented input therefore reproduces only one half of what was published, like trying to recover an average from one of its two values.</p>
<p>The table shows it. For the parents PTHA18 gives a non-zero weight, {match_txt}. Above that the segments run out of logic-tree branches sooner, their weight falls faster, and the mix falls below the unsegmented weight: {implied_txt}. There is a second effect: the plain official input has no Bird convergence profile, so it gives a non-zero probability to the {n_zero} parent{'s' if n_zero != 1 else ''} that PTHA18's model puts at 0. PTHA18's segmented input has both the segments and the Bird convergence, and with it the engine matches in all {n_par:,} parents.</p>
{table(["Mw", "PTHA18 published", "Engine, unsegmented only", "Engine, segmented mix"], unseg_rows)}
<p>The example folder of this pipeline (<code>examples/ptha18_compare/{ZONE}</code>) is {"a segmented run, so its HS table holds the same kind of mix: its parent rates and weights are the engine's zone-level ones (the <code>*_source_zone*.csv</code> files)." if ours_seg else "an unsegmented run, so its HS table holds the unsegmented model only."}</p>"""
    else:
        c_seg_html = ""
    open_seg_li = ("<li><b>Segment definition.</b> Our segments (Berryman end points on our mesh, parameters from the literature and the GCMT catalogue) are close to PTHA18's but not the same, and the rates above Mw 9 stay higher than PTHA18's. Aligning our per-segment parameters with PTHA18's would show how much of what remains is parameters and how much is the mesh; using PTHA18's own unit-source table with <code>--mesh-file</code> would align the parent and row counts.</li>"
                   if full_d1 else
                   "<li><b>Mesh and parameters.</b> This zone is unsegmented in PTHA18 and in our example, so the remaining rate differences come from our mesh, our parameters and the GCMT catalogue. Using PTHA18's own unit-source table with <code>--mesh-file</code> would align the parent and row counts, and aligning our coupling, b and Mw_max with PTHA18's would show how much of the rest is parameters.</li>")
    neg_lon, big_lon = int((D["target_lon"] < 0).sum()), int((D["target_lon"] > 180).sum())
    lon_li = ("<li><b>Longitudes.</b> " + (
        f"PTHA18's file mixes -180 to 180 and 0 to 360 in <code>target_lon</code> ({neg_lon:,} rows are negative, {big_lon:,} above 180). Ours uses 0 to 360 everywhere."
        if neg_lon and big_lon else
        f"PTHA18's file stores <code>target_lon</code> from -180 to 180 ({neg_lon:,} rows are negative); ours uses 0 to 360, so those rows differ by 360 degrees."
        if neg_lon else "PTHA18's file stores <code>target_lon</code> from 0 to 360, like ours.") + "</li>")
    repro_val = (f"python ptha_python/validation/validate_weight_with_nonzero_rate.py --hs-nc nc_files/all_stochastic_slip_earthquake_events_{ZONE}.nc --zone {ZONE} --mode official --out html/data/{ZONE}_validation_official.json"
                 + (f"\npython ptha_python/validation/validate_weight_with_nonzero_rate.py --hs-nc nc_files/all_stochastic_slip_earthquake_events_{ZONE}.nc --zone {ZONE} --mode plain --out html/data/{ZONE}_validation_plain.json" if ptha_seg else "")
                 + f"\npython ptha_python/validation/compare_hs_nc.py --zone {ZONE}")
    # ---- html
    ok_all = same_type and mw_grid_same and all(cO.values())
    n_fmt_ok = int(sum(fmt_ok_list))
    recipe_ok = bool(all(vs[k]["n_within_1e-9"] == vs[k].get("n", vs[k].get("n_compared"))
                         for k in ("weight_with_nonzero_rate", "variable_mu_weight_with_nonzero_rate", "rate_mean (constant)", "variable_mu rate_mean")))
    max_w = max(vs["weight_with_nonzero_rate"]["max_abs_diff"], vs["variable_mu_weight_with_nonzero_rate"]["max_abs_diff"])
    max_r = max(vs["rate_mean (constant)"]["max_rel_diff"], vs["variable_mu rate_mean"]["max_rel_diff"])
    diagram = pipeline_diagram({"units_D": len(olon), "par_D": len(np.unique(D["uniform_event_row"])), "rows_D": len(D["Mw"]),
                                "units_O": len(mr), "par_O": len(np.unique(O["uniform_event_row"])), "rows_O": len(O["Mw"]),
                                "recipe_ok": recipe_ok, "fmt_ok": n_fmt_ok == len(names) and same_type})
    sections = []
    summary = f"""
<div class="cards">
 <div class="card"><div class="k">File format</div><div class="v">{n_fmt_ok} of {len(names)}</div><div class="s">variables with the same name, type and dimension as PTHA18's file (section A)</div></div>
 <div class="card"><div class="k">Same recipe as PTHA18</div><div class="v">{vs['weight_with_nonzero_rate']['n_within_1e-9']:,} of {n_par:,}</div><div class="s">uniform-slip parent ruptures (those inside the engine's magnitude range, {rng_txt}) whose <code>weight_with_nonzero_rate</code> our engine reproduces (largest difference {vs['weight_with_nonzero_rate']['max_abs_diff']:.0e}) when run on PTHA18's own input files (section C)</div></div>
 <div class="card"><div class="k">File consistency checks</div><div class="v">{sum(cO.values())} of {len(cO)}</div><div class="s">tests that PTHA18 itself runs on its files (for example, a weight of 0 must come with a rate of 0). Our file passes all; so does theirs: {sum(cD.values())} of {len(cD)} (section B)</div></div>
 {total_card}
</div>
<div class="note"><b>Words used.</b> A <b>parent</b> is a uniform-slip rupture (one magnitude, one footprint). A <b>row</b> of the table is one random-slip scenario drawn from a parent; each parent has 15 or more rows and shares its rate among them. <b>weight_with_nonzero_rate</b> is the share of the logic tree that allows that magnitude at all (1 = every branch does, 0 = none does).</div>"""
    html_doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HS scenario table, {ZONE}: this pipeline vs PTHA18</title>
<style>
:root{{color-scheme:light;--ink:{INK};--ink2:{INK2};--grid:{GRID};--surf:{SURF};--blue:{BLUE};--orange:{ORANGE}}}
*{{box-sizing:border-box}} body{{margin:0;background:#f6f5f1;color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}}
main{{max-width:1120px;margin:0 auto;padding:28px 20px 80px}} h1{{font-size:28px;margin:0 0 6px}} h2{{font-size:21px;margin:44px 0 6px;padding-top:14px;border-top:2px solid var(--grid)}}
h3{{font-size:16px;margin:22px 0 4px}} p{{margin:6px 0 10px;max-width:880px}} .lead{{color:var(--ink2);font-size:16px;max-width:880px}}
nav{{display:flex;flex-wrap:wrap;gap:8px;margin:16px 0 4px}} nav a{{padding:5px 12px;border:1px solid var(--grid);border-radius:99px;background:#fff;color:var(--ink);text-decoration:none;font-size:13px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:12px;margin:18px 0}} .card{{background:#fff;border:1px solid var(--grid);border-radius:12px;padding:14px 16px}}
.card .k{{font-size:12px;text-transform:uppercase;letter-spacing:.05em;color:var(--ink2)}} .card .v{{font-size:26px;font-weight:650;margin:2px 0}} .card .s{{font-size:13px;color:var(--ink2)}}
.fig{{background:#fff;border:1px solid var(--grid);border-radius:12px;padding:12px;margin:14px 0}} .fig svg{{width:100%;height:auto;display:block}}
.grid2{{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:14px}} figure.chart{{margin:0;background:#fff;border:1px solid var(--grid);border-radius:12px;padding:12px}}
figure.chart svg{{width:100%;height:auto;display:block}} figcaption{{font-size:13px;color:var(--ink2);margin-top:4px}}
.legend{{display:flex;flex-wrap:wrap;gap:14px;font-size:13px;margin-bottom:4px}} .legend i{{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}}
.ax{{font-size:11px;fill:var(--ink2)}} .axl{{font-size:12px;fill:var(--ink2)}} .bt{{font-size:15px;font-weight:650;fill:var(--ink)}} .bs{{font-size:13px;fill:var(--ink2)}} .lane{{font-size:15px;font-weight:700}} .pill{{font-size:12.5px;font-weight:600}}
.tw{{overflow-x:auto;background:#fff;border:1px solid var(--grid);border-radius:12px;margin:10px 0}} table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{padding:7px 10px;text-align:left;border-bottom:1px solid var(--grid);vertical-align:top}} th{{background:#f1efe9;font-weight:600;position:sticky;top:0}} tr:last-child td{{border-bottom:0}}
code{{background:#efede7;padding:1px 5px;border-radius:4px;font-size:12.5px}} .b{{display:inline-block;padding:1px 9px;border-radius:99px;font-size:11.5px;font-weight:700}}
.b.ok{{background:#dcf3e2;color:#14602b}} .b.no{{background:#fde1d6;color:#9a3412}} .b.mid{{background:#fdf0cc;color:#8a5a00}}
.note{{background:#fff8e6;border:1px solid #f0dca0;border-radius:10px;padding:10px 14px;margin:12px 0;max-width:900px}} .good{{background:#eaf7ee;border-color:#b9e0c4}}
.hint{{font-size:13px;color:var(--ink2)}} figure.chart{{cursor:zoom-in}} #zoom{{position:fixed;inset:0;background:rgba(246,245,241,.97);z-index:50;padding:14px 18px;display:flex;flex-direction:column}} #zoom[hidden]{{display:none}} .zbar{{display:flex;gap:10px;align-items:center;margin-bottom:8px;font-size:13px;color:var(--ink2)}} .zbar button{{padding:4px 12px;border:1px solid var(--grid);border-radius:8px;background:#fff;cursor:pointer;font:inherit}} #zbody{{flex:1;min-height:0;display:flex;flex-direction:column;background:#fff;border:1px solid var(--grid);border-radius:12px;padding:10px}} #zbody svg{{flex:1;min-height:0;width:100%;height:100%;cursor:grab;touch-action:none}} pre{{background:#fff;border:1px solid var(--grid);border-radius:10px;padding:12px;overflow-x:auto;font-size:12.5px}}
</style></head><body><main>
<h1>HS scenario table of {ZONE}: this pipeline against PTHA18</h1>
<p class="lead">The nc file from the Drive, <code>all_stochastic_slip_earthquake_events_{ZONE}.nc</code> (PTHA18's published table, the reference, in blue), next to the one step 7d of this pipeline writes from scratch (in orange). The two are compared on four separate questions, so that a real format problem is never mistaken for an expected difference in values.</p>
<nav><a href="#how">How it is made</a><a href="#A">A. Format</a><a href="#B">B. PTHA18 checks</a><a href="#C">C. Exactness</a><a href="#D">D. Values</a><a href="#E">E. Open points</a></nav>
{summary}
{one_sentence}

<h2 id="how">How both tables are made</h2>
{how_html}
<div class="fig">{diagram}</div>

<h2 id="A">A. Format: is the file the same kind of file?</h2>
<p>Name, type and dimension of every variable, plus the conventions inside the text variables. A reader written for PTHA18's file should open ours without changes.</p>
{table(["Property", "PTHA18", "This pipeline", "Result"], attr_rows)}
{order_note_html}
<h3>Global attributes</h3>
<p>Attributes are free-text notes stored in the file (where it came from, which software made it). No reader depends on their content, so they are not part of the pass/fail above. They differ on purpose: PTHA18's point at paths on the NCI computer and record the R session; ours point inside the example folder and record the Python session. <code>R_session_info</code> is not copied, because it would claim an R run that did not happen.</p>
{table(["Attribute", "What it holds", "PTHA18", "This pipeline"], glob_rows)}
<h3>Every variable</h3>
{table(["Variable", "Type", "Meaning", "PTHA18 range", "Our range", "Where our value comes from", "Format"], var_rows)}

<h2 id="B">B. PTHA18's own consistency checks, on both files</h2>
<p>Ports of the checks in PTHA18's <code>check_event_netcdf_files.R</code> (and a few structural ones). They do not depend on the mesh or the catalogue, so both files must pass them.</p>
{table(["Check", "PTHA18 file", "Our file"], chk_rows)}
<div class="note">Two softer facts, both from PTHA18's own file too. The mean rate lies inside [lower_ci, upper_ci] for {iD*100:.1f}% of the scenarios with a rate in PTHA18 and {iO*100:.1f}% in ours. The percentiles are not in order (16th above the median, etc.) for {vD*100:.1f}% of PTHA18's rows and {vO*100:.1f}% of ours: a single scenario's percentile is a difference of percentile curves, not a percentile of the scenario (Davies and Griffin, 2018, Section 3.7.7), so this is expected.</div>

<h2 id="C">C. Exactness: our engine against PTHA18's published numbers</h2>
{c_intro}
{table(["Quantity", "Parents within 1e-9", "Largest difference", "Result"], cval)}
<p>That includes the {n_zero} parent{'s' if n_zero != 1 else ''} that PTHA18 gives weight 0 inside this range.{zero_w_txt}{zero_rate_txt}{beyond_txt} The rates tested are each parent's total (the sum of the rates of its slip fields), because that is the number the HS table shares out.</p>
{c_seg_html}

<h2 id="D">D. Values: where the two tables differ, and why</h2>
{d_note}
<p class="hint">Click any figure to enlarge it. In the enlarged view, scroll to zoom, drag to move, Esc to close.</p>
{d1_html}
<h3>Meshes</h3>
<p>The two meshes side by side. {mesh_extent_txt} Their sizes are what changes the counts of parents and rows in the table on the right.</p>
<div class="grid2">{mesh_chart}
<div>{table(["", "PTHA18", "This pipeline"], [["Unit sources", len(olon), len(mr)], ["Along strike x down dip", f"{o_as} x {o_dd}", f"{m_as} x {m_dd}"], ["Parent (FAUS) ruptures", f"{len(np.unique(D['uniform_event_row'])):,}", f"{len(np.unique(O['uniform_event_row'])):,}"], ["HS rows", f"{len(D['Mw']):,}", f"{len(O['Mw']):,}"], ["Scenarios with rate > 0", f"{int((D['rate_annual'] > 0).sum()):,}", f"{int((O['rate_annual'] > 0).sum()):,}"], ["Total rate (constant)", f"{tot(D, 'rate_annual'):.4f}", f"{tot(O, 'rate_annual'):.4f}"], ["Total rate (variable)", f"{tot(D, 'variable_mu_rate_annual'):.4f}", f"{tot(O, 'variable_mu_rate_annual'):.4f}"]])}</div></div>
<h3>Rates</h3>
<div class="grid2">{c1}{c2}{c3}{cw}</div>
<p>{rate_text}</p>
<h3>Logic-tree weight and rigidity</h3>
<div class="grid2">{c4}{c5}</div>
<h3>Rupture size and slip</h3>
<div class="grid2">{c6}{c7}{kc}{c_lat}</div>
<h3>By magnitude</h3>
{table(["Mw", "Parents (PTHA18 / ours)", "Rows (PTHA18 / ours)", "Sum of rates, PTHA18", "Sum of rates, ours", "Ratio ours / PTHA18", "Mean weight_with_nonzero_rate, PTHA18", "Mean weight, ours"], pm_rows)}

<h2 id="E">E. Open points</h2>
<ul>
<li><b>Unit-source numbering.</b> A user of the HS table needs the unit-source table that gives the numbers in <code>event_index_string</code> a position (PTHA18 ships <code>unit_source_statistics_{ZONE}.nc</code> with it). This pipeline writes the equivalent cells to <code>mesh/</code> and <code>outputs/exports/mesh_cells.csv</code>, not yet in PTHA18's netCDF format.</li>
{open_seg_li}
<li><b>Only the HS table.</b> The variable-area uniform-slip (VAUS) table and the uniform-slip (FAUS) one have the same structure and are not exported yet.</li>
{lon_li}
<li><b>Slip text.</b> Four significant digits, like PTHA18: the numbers are the same kind of rounding, not a loss.</li>
</ul>
<h3>Reproduce</h3>
<pre>python examples/ptha18_compare/{ZONE}/steps/step7_run.py
python examples/ptha18_compare/{ZONE}/steps/step7b_stochastic_slip.py
python examples/ptha18_compare/{ZONE}/steps/step7c_variable_mu.py
python examples/ptha18_compare/{ZONE}/steps/step7d_export_nc.py
{repro_val}</pre>
</main>{ZOOM}</body></html>"""
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(html_doc)
    print(f"wrote {a.out} ({os.path.getsize(a.out)/1e3:.0f} kB)")
    print("all structural checks pass:", ok_all)


if __name__ == "__main__":
    main()
