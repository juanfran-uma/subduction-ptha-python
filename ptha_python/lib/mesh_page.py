"""The 3D page of an example's mesh: <example>/mesh/mesh_3d.html.

One self-contained page (Plotly from a CDN, English, light) with the mesh in 3D and
in plan, the SLAB contours step 1 cut it from (so you can see the cut), the trench,
and the L x W of every cell as numbers. Built from files the example already has:

    data/slab2/unit_source_grid.npy          the mesh (step 2)
    data/convergence_per_cell.npy            Bird convergence per cell (step 3), optional
    data/slab2/unit_source_grid.meta.json    the target cell size, optional
    inputs/geometry/*_slab2_contours.shp     the SLAB contours (step 1), optional
                                             (absent for an external --mesh-file)

export_mesh.write_mesh_folder() calls write_page(); it can also be run alone:

    python ptha_python/lib/mesh_page.py examples/new_zones_40km/calabria2 [...]
"""

import glob
import json
import os
import sys

import numpy as np
from pyproj import Proj

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from pyptha import unit_sources as us  # noqa: E402

TEMPLATE = os.path.join(os.path.dirname(HERE), "templates", "mesh_3d.html.tmpl")


def _round(a, nd=2):
    return np.round(np.asarray(a, dtype=float), nd).tolist()


def _contours(example_dir, proj, lon0):
    """SLAB contours as [{level, lines: [[(x, y), ...], ...]}], x/y in km; [] if none."""
    found = sorted(glob.glob(os.path.join(example_dir, "inputs", "geometry", "*_slab2_contours.shp")))
    if not found:
        return []
    import geopandas as gpd
    out = []
    for _, row in gpd.read_file(found[0]).sort_values("level").iterrows():
        geom = row.geometry
        parts = list(geom.geoms) if geom.geom_type.startswith("Multi") else [geom]
        lines = []
        for g in parts:
            xy = np.array(g.coords)[:, :2]
            lon = (xy[:, 0] - lon0 + 180) % 360 - 180 + lon0
            x, y = proj(lon, xy[:, 1])
            lines.append([list(p) for p in zip(_round(x / 1000.0), _round(y / 1000.0))])
        out.append({"level": float(row["level"]), "lines": lines})
    return out


def page_data(example_dir):
    grid = np.load(os.path.join(example_dir, "data", "slab2", "unit_source_grid.npy"))   # [row, (lon, lat, depth), col]
    n_rows, n_cols = grid.shape[0] - 1, grid.shape[2] - 1
    lon = grid[:, 0, :]
    lon0 = float(np.degrees(np.arctan2(np.sin(np.radians(lon)).mean(), np.cos(np.radians(lon)).mean())))
    lon = (lon - lon0 + 180) % 360 - 180 + lon0
    proj = Proj(proj="tmerc", lat_0=0, lon_0=lon0, k_0=0.9996, datum="WGS84")
    x, y = proj(lon, grid[:, 1, :])
    P = np.stack([x / 1000.0, y / 1000.0, -grid[:, 2, :]], axis=-1)   # [row, col, (x, y, z)] in km, z down negative

    def dist(a, b):
        return np.linalg.norm(a - b, axis=-1)

    L = (dist(P[:-1, 1:], P[:-1, :-1]) + dist(P[1:, 1:], P[1:, :-1])) / 2     # top and bottom edges
    W = (dist(P[1:, :-1], P[:-1, :-1]) + dist(P[1:, 1:], P[:-1, 1:])) / 2     # left and right edges
    a, b, c, d = P[:-1, :-1], P[:-1, 1:], P[1:, 1:], P[1:, :-1]
    A = (np.linalg.norm(np.cross(b - a, c - a), axis=-1) + np.linalg.norm(np.cross(c - a, d - a), axis=-1)) / 2
    depth = (grid[:-1, 2, :-1] + grid[:-1, 2, 1:] + grid[1:, 2, 1:] + grid[1:, 2, :-1]) / 4

    # the engine's own cell order is along strike outer, down dip inner
    st = us.discretized_source_approximate_summary_statistics(grid)

    def per_cell(v):
        return np.asarray(v, dtype=float).reshape(n_cols, n_rows).T

    data = {
        "name": os.path.basename(os.path.normpath(example_dir)),
        "R": n_rows, "C": n_cols,
        "x": _round(P[..., 0], 1), "y": _round(P[..., 1], 1), "z": _round(P[..., 2], 2),
        "L": _round(L, 1), "W": _round(W, 1), "A": _round(A, 0), "D": _round(depth, 1),
        "strike": _round(per_cell(st["strike"]), 0), "dip": _round(per_cell(st["dip"]), 0),
        "lon": _round(per_cell(st["lon_c"]), 3), "lat": _round(per_cell(st["lat_c"]), 3),
        "trench_km": float(dist(P[0, 1:], P[0, :-1]).sum()),
        "deep_km": float(dist(P[-1, 1:], P[-1, :-1]).sum()),
        "depth_max": float(grid[:, 2, :].max()),
        "tL": None, "tW": None, "conv": None,
        "contours": _contours(example_dir, proj, lon0),
    }
    meta = os.path.join(example_dir, "data", "slab2", "unit_source_grid.meta.json")
    if os.path.exists(meta):
        with open(meta, encoding="utf-8") as f:
            m = json.load(f)
        data["tL"] = m.get("desired_unit_source_length")
        data["tW"] = m.get("desired_unit_source_width")
    conv = os.path.join(example_dir, "data", "convergence_per_cell.npy")
    if os.path.exists(conv):
        v = np.load(conv)
        if len(v) == n_rows * n_cols:
            data["conv"] = _round(per_cell(v), 2)
    return data


def write_page(example_dir, out_dir=None):
    """Write <example_dir>/mesh/mesh_3d.html and return its path."""
    data = page_data(example_dir)
    n = f"{data['C']} x {data['R']}"
    target = f" The target cell size was {data['tL']:g} x {data['tW']:g} km." if data["tL"] and data["tW"] else ""
    title = f"{data['name']}: mesh in 3D"
    sub = (f"{n} = {data['C'] * data['R']} unit sources (columns x rows) of {data['name']}, with the SLAB contours they were cut from." if data["contours"]
           else f"{n} = {data['C'] * data['R']} unit sources (columns x rows) of {data['name']}.") + target
    with open(TEMPLATE, encoding="utf-8") as f:
        html = f.read()
    html = (html.replace("__TITLE__", title).replace("__SUBTITLE__", sub)
            .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    out_dir = out_dir or os.path.join(example_dir, "mesh")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "mesh_3d.html")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(html)
    return path


if __name__ == "__main__":
    for d in sys.argv[1:]:
        print(write_page(d))
