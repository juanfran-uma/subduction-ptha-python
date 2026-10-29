"""Mesh edit (--columns *-edited): redesign a trench (rptha) mesh whose cells misread the slab's dip.

rptha's mesh puts the rows by an along-strike rule. Where the arc bends (an S shaped zone such as cyprus2),
a few rows cross the depth contours almost along strike: those cells are slivers (engine width of 6 to 12 km
against 40 km), and the dip the engine reads from them (5 to 14 deg) is far from the dip of the slab under
them (28 to 47 deg). The moment balance uses that dip, so the defect reaches the rates.

``edit_if_ugly`` keeps the trench mesh untouched unless some cell's engine dip is more than DIP_TOL_DEG from
the dip of its own plane. Then it designs a replacement with the same number of columns and rows:

  1. the slab surface is step 1's contours (the first and the last row stay on the first and the last one);
  2. every column is a down-dip line crossing the contours; where it crosses each is a free variable
     (the two end columns stay on the ends of the contours);
  3. the rows are at equal slope distance along every column;
  4. the variables are fitted (scipy least_squares) so the corners are close to 90 deg, the strike edges
     close to the desired length, the dip edges close to the desired width, the columns do not cross and
     each column is smooth.

The replacement is used only if it passes ``accept`` (fewer cells over the tolerance, no more cells narrower than
half the desired width, no folded or overlapping cell, no crossing columns, true 3D area within AREA_TOL); otherwise the original mesh is kept and the reason is logged.
"""

from __future__ import annotations

import numpy as np
from pyproj import Proj
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix

from . import contour_discretisation as cd
from . import unit_sources as us

VERSION = "1"            # enters the mesh fingerprint (grid_cache): bump it when the design or the checks below change
DIP_TOL_DEG = 8.0
AREA_TOL = 0.05
# strike length, dip length, right angles, smooth columns, columns not crossing (tuned on cyprus2)
WEIGHTS = dict(L=3.0, W=0.5, A=1.5, S=0.3, G=3.0)


def _projection(lon, lat):
    """Local transverse Mercator in km, centred on the points. Longitudes may be 0-360 or -180-180."""
    lon = np.asarray(lon, float)
    lon0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(lon))))))
    proj = Proj(proj="tmerc", lat_0=0, lon_0=lon0, k_0=0.9996, datum="WGS84")
    wrap = 360.0 if float(np.max(lon)) > 180.0 else 0.0

    def fwd(lo, la):
        x, y = proj((np.asarray(lo, float) - lon0 + 180) % 360 - 180 + lon0, np.asarray(la, float))
        return np.stack([x, y], axis=-1) / 1000.0

    def inv(xy):
        lo, la = proj(xy[..., 0] * 1000.0, xy[..., 1] * 1000.0, inverse=True)
        return (lo % wrap if wrap else lo), la

    return fwd, inv


def _engine_dip(grid):
    R, C = grid.shape[0] - 1, grid.shape[2] - 1
    st = us.discretized_source_approximate_summary_statistics(grid)
    return np.asarray(st["dip"], float).reshape(C, R).T      # along strike outer, down dip inner


def _plane(grid):
    """(dip deg, area km2) of every cell's own plane, from its two diagonals; shape (rows, cols)."""
    fwd, _ = _projection(grid[:, 0, :], grid[:, 1, :])
    xy = fwd(grid[:, 0, :], grid[:, 1, :])
    N = np.concatenate([xy, -grid[:, 2, :, None]], axis=-1)              # (R+1, C+1, 3): x, y, z up (km)
    n = np.cross(N[1:, 1:] - N[:-1, :-1], N[1:, :-1] - N[:-1, 1:])
    norm = np.linalg.norm(n, axis=-1)
    return np.degrees(np.arccos(np.clip(np.abs(n[..., 2]) / norm, 0, 1))), norm / 2


def dip_mismatch(grid):
    """abs(dip the engine reads - dip of the cell's own plane), deg, shape (rows, cols)."""
    return np.abs(_engine_dip(grid) - _plane(grid)[0])


class _Surface:
    """The slab surface as contour lines in km, each ordered like the mesh's columns."""

    def __init__(self, levels, lines):
        self.levels = np.asarray(levels, float)
        self.lines = lines
        self.cum = [np.r_[0, np.cumsum(np.hypot(*np.diff(k, axis=0).T))] for k in lines]

    def point(self, k, u):
        s = np.clip(u, 0, 1) * self.cum[k][-1]
        return np.stack([np.interp(s, self.cum[k], self.lines[k][:, 0]), np.interp(s, self.cum[k], self.lines[k][:, 1])], axis=-1)

    def foot(self, k, p, n=4000):
        uu = np.linspace(0, 1, n)
        pts = self.point(k, uu)
        return uu[np.argmin(np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1]))]

    def start_fractions(self, C, alpha=0.5):
        """First guess: perpendicular feet from the trench blended with uniform fractions."""
        u = np.zeros((len(self.levels), C + 1))
        u[0] = np.linspace(0, 1, C + 1)
        for k in range(1, len(self.levels)):
            for j in range(1, C):
                u[k, j] = self.foot(k, self.point(k - 1, u[k - 1, j]))
            u[k, C] = 1.0
        return (1 - alpha) * u + alpha * np.linspace(0, 1, C + 1)[None, :]

    def nodes(self, u, R):
        """(R+1, C+1, 3): x, y (km), depth; R+1 rows at equal slope distance along every column."""
        P = np.stack([self.point(k, u[k]) for k in range(len(self.levels))])
        out = np.zeros((R + 1, P.shape[1], 3))
        for j in range(P.shape[1]):
            xyz = np.c_[P[:, j, 0], P[:, j, 1], self.levels]
            s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]
            for a in range(3):
                out[:, j, a] = np.interp(np.linspace(0, s[-1], R + 1), s, xyz[:, a])
        return out


def _surface(contours, grid, fwd):
    """The contours between the trench row's depth and the last row's depth, as a _Surface; None if fewer than two."""
    cs = sorted(((float(d), np.asarray(c, float)[:, :2]) for d, c in contours), key=lambda t: t[0])
    levels = np.array([d for d, _ in cs])
    first = levels[np.argmin(np.abs(levels - np.median(grid[0, 2, :])))]
    last = levels[np.argmin(np.abs(levels - np.median(grid[-1, 2, :])))]
    ref0 = fwd(grid[:, 0, 0], grid[:, 1, 0]).mean(axis=0)
    ref1 = fwd(grid[:, 0, -1], grid[:, 1, -1]).mean(axis=0)
    keep, lines = [], []
    for d, c in cs:
        if not first <= d <= last:
            continue
        xy = fwd(c[:, 0], c[:, 1])
        xy = xy[np.r_[True, np.any(np.diff(xy, axis=0) != 0, axis=1)]]
        if np.hypot(*(xy[0] - ref0)) + np.hypot(*(xy[-1] - ref1)) > np.hypot(*(xy[-1] - ref0)) + np.hypot(*(xy[0] - ref1)):
            xy = xy[::-1]
        keep.append(d)
        lines.append(xy)
    return _Surface(keep, lines) if len(keep) >= 2 else None


def _residuals(v, surf, C, R, length, width, w):
    K = len(surf.levels)
    u = np.zeros((K, C + 1))
    u[:, C] = 1.0
    u[:, 1:C] = v.reshape(K, C - 1)
    X = surf.nodes(u, R)
    es, ed = X[:, 1:] - X[:, :-1], X[1:] - X[:-1]
    ls, ld = np.linalg.norm(es, axis=2), np.linalg.norm(ed, axis=2)
    res = [w["L"] * (ls / length - 1).ravel(), w["W"] * (ld / width - 1).ravel()]
    a, b = es[:-1] / ls[:-1, :, None], es[1:] / ls[1:, :, None]
    c, d = ed[:, :-1] / ld[:, :-1, None], ed[:, 1:] / ld[:, 1:, None]
    res += [w["A"] * np.sum(p * q, axis=2).ravel() for p, q in ((a, c), (a, d), (b, c), (b, d))]
    vv = u[:, 1:C]
    res.append(w["S"] * 100 * (vv[2:] - 2 * vv[1:-1] + vv[:-2]).ravel())
    res.append(w["G"] * 100 * np.minimum(np.diff(u, axis=1) - 0.2 / C, 0).ravel())
    return np.concatenate(res)


def _sparsity(K, C, R):
    """Which residual of _residuals depends on which variable (a scipy sparse matrix).

    Node column m (1..C-1) owns the variables k * (C - 1) + m - 1, k = 0..K-1. Every residual involves one column
    (dip edges, smoothness) or two neighbours (strike edges, corners, crossing), in the order _residuals returns them.
    Declared to least_squares so the numerical Jacobian needs a handful of evaluations instead of one per variable.
    """
    pair = lambda shape: np.broadcast_to(np.arange(shape[1]), shape)            # residual (i, j): node columns j, j + 1  # noqa: E731
    blocks = [(pair((R + 1, C)), 2), (np.broadcast_to(np.arange(C + 1), (R, C + 1)), 1)] + [(pair((R, C)), 2)] * 4 \
        + [(np.broadcast_to(np.arange(1, C), (K - 2, C - 1)), 1), (pair((K, C)), 2)]
    rows, cols, r0 = [], [], 0
    for first, width in blocks:
        first = np.ravel(first)
        for off in range(width):
            m = first + off
            hit = np.nonzero((m >= 1) & (m <= C - 1))[0]
            for k in range(K):
                rows.append(r0 + hit)
                cols.append(k * (C - 1) + m[hit] - 1)
        r0 += first.size
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    return coo_matrix((np.ones(rows.size), (rows, cols)), shape=(r0, K * (C - 1)))


def redesign(grid, contours, desired_length_km, desired_width_km, weights=WEIGHTS):
    """A new grid with the shape of ``grid`` (rows+1, 3, cols+1) on the same contours, or None if it cannot be built."""
    R, C = grid.shape[0] - 1, grid.shape[2] - 1
    if C < 3:
        return None
    fwd, inv = _projection(grid[:, 0, :], grid[:, 1, :])
    surf = _surface(contours, grid, fwd)
    if surf is None:
        return None
    K = len(surf.levels)
    v0 = surf.start_fractions(C)[:, 1:C].ravel()
    sol = least_squares(_residuals, v0, args=(surf, C, R, desired_length_km, desired_width_km, weights),
                        bounds=(0, 1), max_nfev=200, x_scale=0.02, jac_sparsity=_sparsity(K, C, R))
    u = np.zeros((K, C + 1))
    u[:, C] = 1.0
    u[:, 1:C] = sol.x.reshape(K, C - 1)
    if not (np.diff(u, axis=1) > 0).all():
        return None                                         # columns cross
    X = surf.nodes(u, R)
    lon, lat = inv(X[..., :2])
    return np.stack([lon, lat, X[..., 2]], axis=1)


def _thin(grid, desired_width_km):
    """Cells whose engine width (area / length) is less than half the desired width."""
    width = np.asarray(us.discretized_source_approximate_summary_statistics(grid)["width"], float)
    return int((width < desired_width_km / 2).sum())


def accept(old, new, desired_width_km):
    """(True, '') if ``new`` is a safe replacement of ``old``, else (False, why)."""
    if new is None:
        return False, "no valid design (too few contours, too few columns or crossing columns)"
    defects = cd.mesh_defects(new)
    if defects["total"]:
        return False, f"the design has defects {defects}"
    n_old, n_new = int((dip_mismatch(old) > DIP_TOL_DEG).sum()), int((dip_mismatch(new) > DIP_TOL_DEG).sum())
    if n_new >= n_old:
        return False, f"cells over {DIP_TOL_DEG:g} deg did not decrease ({n_old} -> {n_new})"
    thin_old, thin_new = _thin(old, desired_width_km), _thin(new, desired_width_km)
    if thin_new > thin_old:
        return False, f"cells narrower than half the desired width increased ({thin_old} -> {thin_new})"
    change = float(_plane(new)[1].sum() / _plane(old)[1].sum() - 1)       # true 3D area: the engine's own area is what slivers distort
    if abs(change) > AREA_TOL:
        return False, f"the area changes by {100 * change:+.1f}% (limit {100 * AREA_TOL:g}%)"
    return True, ""


def edit_if_ugly(grid, contours, desired_length_km, desired_width_km, log=print):
    """``grid`` itself unless some cell's dip is misread; then the redesign if it is safe, else ``grid`` and the reason."""
    mis = dip_mismatch(grid)
    bad = int((mis > DIP_TOL_DEG).sum())
    if bad == 0:
        log(f"  mesh edit (--columns *-edited): no cell's dip is more than {DIP_TOL_DEG:g} deg from its plane "
            f"(worst {mis.max():.1f} deg); the original mesh is kept")
        return grid
    log(f"  mesh edit (--columns *-edited): {bad} of {mis.size} cells read a dip more than {DIP_TOL_DEG:g} deg from their plane "
        f"(worst {mis.max():.1f} deg); designing a replacement on the same contours ...")
    try:
        new = redesign(grid, contours, desired_length_km, desired_width_km)
    except Exception as exc:                                # a failed design must not stop the run: the original mesh is valid
        log(f"  mesh edit (--columns *-edited): the design failed ({type(exc).__name__}: {exc}); the original mesh is kept")
        return grid
    ok, why = accept(grid, new, desired_width_km)
    if not ok:
        log(f"  mesh edit (--columns *-edited): replacement rejected, the original mesh is kept: {why}")
        return grid
    after = dip_mismatch(new)
    log(f"  mesh edit (--columns *-edited): replacement used; cells over {DIP_TOL_DEG:g} deg {bad} -> {int((after > DIP_TOL_DEG).sum())}, "
        f"worst {mis.max():.1f} -> {after.max():.1f} deg")
    return new
