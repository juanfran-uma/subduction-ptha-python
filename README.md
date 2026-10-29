# PTHA in Python

A Python port of the rupture and rate parts of **rptha**, the code behind the 2018 Australian probabilistic tsunami hazard assessment (PTHA18, Davies and Griffin, 2018). For a subduction zone covered by the SLAB2.0 or SLAB1.0 geometry models, it builds the earthquake sources on the plate interface, the logic tree of recurrence models, and the exceedance rate of earthquakes per year. It uses public data only.

The code stops at earthquake rates. It does not simulate tsunamis or compute hazard curves.

## What you get

Each run creates a folder that holds:

- `steps/`: nine steps as plain Python scripts, run in order. `step_total.py` runs them all.
- `inputs/`, `data/`, `outputs/`: what each step reads and writes.
- `report.html`: one page (its figures are embedded; only the map loads from the internet) that explains every step, with the numbers it produced and the exceedance-rate curve. When step 8 has run, it also compares the result with PTHA18's own run.
- `RUN.html`: the exact commands that created the folder.

Steps 7b (heterogeneous and variable-area slip), 7c (variable shear modulus, only with `--variable-mu on`) and 8 (PTHA18's official run, which needs R) are optional.

## Quick start

You need Python 3.13 (the version the code was run with) and the packages in `requirements.txt`. Every command is run from the repository root, the folder that `git clone` creates.

```
git clone <repository URL>
cd <the cloned folder>
python -m venv .venv
.venv\Scripts\Activate.ps1             (Windows PowerShell; cmd: .venv\Scripts\activate.bat; macOS and Linux: source .venv/bin/activate)
pip install -r requirements.txt
python ptha_python/list_regions.py     (the SLAB regions and the zone names)
python ptha_python/generate.py kermadec --zone kermadectonga2 --folder kermadectonga2
python kermadectonga2/steps/step_total.py --skip-official
```

Then open `kermadectonga2/report.html` in a browser. On macOS and Linux the Python command may be `python3` until the environment is active. The SLAB2.0 raster for this zone ships in `ptha_python/data/slab2/`, so step 1 downloads nothing. Step 4 downloads the GCMT catalogue (about 30 MB), so the first run needs a network connection. In the reference runs, a complete run with step 8 skipped took 2.1 to 4.6 minutes.

`--skip-official` leaves out step 8, the PTHA18 comparison, which needs R with the rptha package and a 1.34 GB download. Everything else works without R. If you leave the flag out and R is missing, step 8 stops with a message, the run goes on, and the report has no comparison with PTHA18. The folder `kermadectonga2/` holds everything the run writes: `report.html`, `outputs/` (the CSV files of the rates), `inputs/`, `data/` and `figures/`.

## Options of `generate.py`

`python ptha_python/generate.py NAME [options]` creates the folder. Run `python ptha_python/generate.py --help` for the full text.

| Option | Needed? | Default | What it does |
|---|---|---|---|
| `NAME` | **NEEDED** | (required) | Zone name, for example `cascadia` or `'south america'`. It selects the SLAB region. |
| `--zone` | OPTIONAL (guessed from the name) | guessed from `NAME` | The PTHA18 zone name. With `--ptha false`, a zone listed in `lib/berryman_params.py`. |
| `--folder` | OPTIONAL | `NAME` in lower case, with other characters replaced by `_` | The folder to create. It can be nested, for example `runs/calabria2`. |
| `--id` | OPTIONAL (only when the region name is ambiguous) | from `NAME` | The ScienceBase item id of the SLAB2.0 region, when the name is ambiguous. |
| `--ptha` | OPTIONAL | `false` | `true` lets the steps read PTHA18's own cell width, scaling and row count, and needs a zone PTHA18 modelled. `false` uses public data only. |
| `--discretizer` | OPTIONAL | `optimal` | How contours become unit sources: `optimal` (rptha's mesh, repaired only if broken), `lm` (rptha's method as PTHA18 runs it) or `mid` (a method that rptha has deprecated and PTHA18 does not use). |
| `--columns` | OPTIONAL | `trench` | Unit sources along strike: `trench` (trench length / 50 km, PTHA18's rule), `average`, or `trench-edited` / `average-edited`. |
| `--cell-size` | OPTIONAL | `mean` | `mean` (about 50 x 50 km on average, PTHA18's rule) or `strasser` (no cell larger than the smallest rupture divided by `--cell-k`). |
| `--cell-k` | OPTIONAL (only with --cell-size strasser) | `1.5` | With `--cell-size strasser`, the divisor applied to the smallest rupture. |
| `--clip` | OPTIONAL | none | `LON_MIN,LON_MAX,LAT_MIN,LAT_MAX`: use only that window of the SLAB raster. |
| `--mesh-file` | OPTIONAL | none | Use an external quadrilateral mesh instead of building one from SLAB. |
| `--mesh-depth-units` | OPTIONAL (only with --mesh-file) | `auto` | Depth units of `--mesh-file`: `auto` (metres if any depth exceeds 200, else km), `m` or `km`. |
| `--variable-mu` | OPTIONAL | `off` | `on` also runs step 7c: rates with a variable shear modulus next to the constant one. |
| `--rupture-size` | OPTIONAL | `rptha` | `rptha` (one block of cells per magnitude, PTHA18's rule) or `local` (block sized by the real area of the cells at each place). |
| `--convergence` | OPTIONAL | `bird` | Plate-boundary table of step 3: `bird` (Bird, 2003) or `bird-griffin` (the table PTHA18 used). |
| `--segmented` | OPTIONAL | `false` | `true` also runs the zone as along-strike segments, with the LEVEL 0 weights. |
| `--segment-boundaries` | OPTIONAL (only with --segmented true) | `berryman` | Where segments end: `berryman` (end points of Berryman et al., 2015, Table 3.1) or `bird` (where Bird's plate pair changes). |

## External mesh

If you already have a quadrilateral mesh, step 1 and step 2 can read it instead of building one from SLAB.

1. Write a text file with one line per unit source. Each line holds 12 numbers: the four corners of the cell as `lon lat depth`, going round the cell. The file `inputs_meshes/alaskaaleutians_quadrilateral_coors.dat` is an example.
2. The mesh must be structured: rows down dip, columns along strike, and neighbouring cells share their corners.
3. Run:

```
python ptha_python/generate.py alaska --zone alaskaaleutians --folder alaskaaleutians --mesh-file inputs_meshes/alaskaaleutians_quadrilateral_coors.dat
python alaskaaleutians/steps/step_total.py
```

The file is copied into the folder's `inputs/geometry/`. `--mesh-file` cannot be combined with `--cell-size strasser` or with the `-edited` values of `--columns`.

## Repository layout

```
ptha_python/          the package
  generate.py         creates a zone folder from the step templates
  list_regions.py     lists the SLAB regions and the zone names
  lib/                zone parameters and the geometry and convergence helpers the steps use
  pyptha/             the rupture and rate engine
  logic_tree/         the logic-tree driver that the steps call
  templates/          the step scripts and the report that generate.py fills in
  data/               public inputs shipped with the package (Bird tables, SLAB2.0 rasters)
  html/docs/          the documentation (start at html/docs/index.html)
  validation/         checks against PTHA18's published files
rptha/                data files copied from the rptha R package (its R code is not included)
official_ptha_data/   downloads PTHA18's saved session and extracts its trees (needs R)
inputs_meshes/        example external meshes for --mesh-file
requirements.txt, LICENSE
```

## What needs the network

- **Step 1** downloads a SLAB raster only when it is not already in the zone folder (`data/slab2/` or `data/slab1/`) or in `ptha_python/data/slab2/` or `ptha_python/data/slab1/`. The SLAB2.0 rasters for `cal`, `car`, `hel`, `ker`, `mak`, `png`, `puy` and `sum` are cached; no SLAB1.0 raster is, so a zone whose name does not end in 2 (for example `cascadia`) always downloads one the first time. SLAB2.0 comes from ScienceBase and SLAB1.0 from USGS. A run with `--mesh-file` downloads nothing in step 1. Both sites sometimes block scripts. If step 1 reports that, download the file in a browser and put it in the folder named in the message.
- **Step 4** downloads the Global CMT catalogue (NDK files from the Lamont-Doherty archive, the GCMT project): one cumulative file of about 23 MB and one file per month since January 2021, about 30 MB in all. Files already in `data/gcmt/` are not downloaded again, and months not yet published are skipped. It keeps the window 1976-01-01 to 2017-03-01 that PTHA18 used.
- **The validation checks** download PTHA18's published files from NCI the first time they need them (see `ptha_python/html/docs/validation.html`).
- **Step 8** downloads PTHA18's saved session (1.34 GB, once) from NCI's THREDDS server. It also needs R with the rptha package, which is not in this repository. See the rptha GitHub repository (GeoscienceAustralia/ptha) for the R package.
- **The map in the step 2 section of the report** loads Leaflet (from cdnjs.cloudflare.com) and Esri tiles when the page is opened in a browser.
- Everything else runs offline, including the Bird tables, which are in `ptha_python/data/bird/`.

## Validation

The comparisons with PTHA18 are in `ptha_python/validation/`. Start with `validate_against_ptha18.py`. The page `ptha_python/html/docs/validation.html` lists every check.

## Documentation

Open `ptha_python/html/docs/index.html` in a browser. It explains the method, the nine steps, the options and the report, and the troubleshooting messages.

## Licence and credits

- The code is under the BSD 3-Clause licence (see `LICENSE`).
- rptha is Copyright (c) 2015, Geoscience Australia, under the BSD 3-Clause licence (see `rptha/LICENSE`). The original code is at https://github.com/GeoscienceAustralia/ptha.
- PTHA18: Davies, G. and Griffin, J. (2018), *The 2018 Australian probabilistic tsunami hazard assessment: hazard from earthquake generated tsunamis*, Geoscience Australia Record 2018/41.
- Segment end points and literature values: Berryman, K. et al. (2015), *The GEM Faulted Earth Subduction Interface Characterisation Project, Version 2.0*.
- Plate boundaries: Bird, P. (2003), *An updated digital model of plate boundaries*, Geochemistry, Geophysics, Geosystems 4(3).
- Geometry: SLAB2.0 (USGS, ScienceBase) and SLAB1.0 (USGS). Earthquakes: the Global CMT catalogue.
