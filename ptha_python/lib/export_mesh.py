"""Write an example's mesh as plain files, in <example>/mesh/.

    mesh.dat               one line per cell: lon lat depth x 4 corners, depth in km
                           NEGATIVE down (the --mesh-file format; pass it to generate.py
                           with --mesh-depth-units km or auto, never m)
    mesh_cell_corners.csv  the same cells with their (downdip_number, alongstrike_number),
                           the 4 corners as lon, lat, depth_km (positive down) and the
                           Bird convergent slip of the cell (mm/yr), as step 3 computed
                           it (not the constant override, if CONVERGENCE_OVERRIDE is set)
    unit_source_grid.npy   the engine's grid, [row, (lon, lat, depth_km), column]
    mesh_3d.html           the mesh in 3D and in plan with the SLAB contours and the L x W
                           of every cell (lib/mesh_page.py)

Everything comes from files step 2 and step 3 already wrote
(data/slab2/unit_source_grid.npy, data/convergence_per_cell.npy), so it runs
in a second. generate.py's step 3 calls write_mesh_folder(); for examples
generated before that existed:

    python ptha_python/lib/export_mesh.py examples/new_zones_40km/calabria2 [...]
    python ptha_python/lib/export_mesh.py --all examples          # every example under it
"""

import argparse
import contextlib
import glob
import io
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from pyptha import mesh_file  # noqa: E402
import mesh_page  # noqa: E402

HEADER = (["downdip_number", "alongstrike_number"]
          + [f"{a}{k}{b}" for k in range(1, 5) for a, b in (("lon", ""), ("lat", ""), ("depth", "_km"))]
          + ["convergent_slip_mm_per_yr"])


def cells(grid):
    """Yield (downdip_number, alongstrike_number, 4 corners) in the engine's unit-source
    order (along strike outer, down dip inner). Corners go trench-left, trench-right,
    deep-right, deep-left, as an external mesh file lists them."""
    for i in range(grid.shape[2] - 1):
        for j in range(grid.shape[0] - 1):
            yield j + 1, i + 1, [grid[j, :, i], grid[j, :, i + 1], grid[j + 1, :, i + 1], grid[j + 1, :, i]]


def write_mesh_folder(example_dir, check=True):
    """Write <example_dir>/mesh/ and return its path."""
    grid = np.load(os.path.join(example_dir, "data", "slab2", "unit_source_grid.npy"))
    conv = np.load(os.path.join(example_dir, "data", "convergence_per_cell.npy"))
    n_cells = (grid.shape[0] - 1) * (grid.shape[2] - 1)
    if len(conv) != n_cells:
        raise ValueError(f"{example_dir}: {len(conv)} convergence values for {n_cells} cells; "
                         "re-run steps 2 and 3 so both describe the same mesh")

    dat, lines = [], [",".join(HEADER)]
    for n, (dd, al, corners) in enumerate(cells(grid)):
        dat.append(np.concatenate([[lon, lat, -dep] for lon, lat, dep in corners]))
        row = [str(dd), str(al)]
        for lon, lat, dep in corners:
            row += [f"{lon:.6f}", f"{lat:.6f}", f"{dep:.4f}"]
        lines.append(",".join(row + [f"{float(conv[n]):.6f}"]))

    # nothing is written to mesh/ unless mesh.dat reads back as the engine's grid
    # (the reader needs two or more columns to tell left from right)
    if check and grid.shape[2] > 2:
        with tempfile.TemporaryDirectory() as tmp:
            probe = os.path.join(tmp, "mesh.dat")
            np.savetxt(probe, np.array(dat), fmt="%.8f", newline="\n")
            with contextlib.redirect_stdout(io.StringIO()):  # the reader logs as if it were an external mesh
                back = mesh_file.read_quadrilateral_mesh(probe, depth_units="km")
        if back.shape != grid.shape or not np.allclose(back, grid, atol=1e-4):
            raise ValueError("mesh.dat does not read back as the example's unit_source_grid.npy")

    out = os.path.join(example_dir, "mesh")
    os.makedirs(out, exist_ok=True)
    np.savetxt(os.path.join(out, "mesh.dat"), np.array(dat), fmt="%.8f", newline="\n")
    with open(os.path.join(out, "mesh_cell_corners.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    shutil.copyfile(os.path.join(example_dir, "data", "slab2", "unit_source_grid.npy"),
                    os.path.join(out, "unit_source_grid.npy"))
    try:  # the page is a view of the files above: its failure must not undo them
        mesh_page.write_page(example_dir, out)
    except Exception as e:  # noqa: BLE001
        print(f"  WARNING: mesh_3d.html not written ({e})")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="example folders, or with --all the folders that contain examples")
    ap.add_argument("--all", action="store_true",
                    help="treat each path as a parent: every <path>/*/<zone> with a convergence_per_cell.npy")
    args = ap.parse_args()

    folders = []
    for p in args.paths:
        if args.all:
            folders += sorted(os.path.dirname(os.path.dirname(f)) for f in glob.glob(
                os.path.join(p, "*", "*", "data", "convergence_per_cell.npy")))
        else:
            folders.append(p)
    if not folders:
        print("no example found (with --all the path is the folder that holds the groups, e.g. examples)")
        return 1
    failed = 0
    for d in folders:
        try:
            print(f"  {d} -> {write_mesh_folder(d)}")
        except (OSError, ValueError, IndexError, KeyError) as e:
            failed += 1
            print(f"  {d}: SKIPPED ({e})")
    print(f"{len(folders) - failed} written, {failed} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
