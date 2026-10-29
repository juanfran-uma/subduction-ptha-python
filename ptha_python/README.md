# ptha_python: package reference

This folder holds the package that builds a source-zone folder, the engine that computes the rates, and the documentation. The repository's front page is `../README.md`. The documentation pages are in `html/docs/`, starting with `html/docs/index.html`.

## Layout

```
ptha_python/
  generate.py          builds a zone folder from the templates; the command-line options
  list_regions.py      lists the 27 SLAB2.0 regions and the SLAB1.0 zones
  lib/                 zone parameters and helpers used by the steps
  pyptha/              the engine: scaling, unit sources, events, rates, slip, variable shear modulus
  logic_tree/          run_logic_tree.py, the logic-tree driver that steps 7, 7c and 8 call
  templates/           the step scripts, RUN.html and the report, filled in by generate.py
  data/                Bird tables, PTHA18's slip quantile weights, cached SLAB2.0 rasters
  html/docs/           the documentation pages, a shared stylesheet (assets/) and figures (img/)
  validation/          checks against PTHA18's published files
```

What each part holds:

- **`lib/`**: `berryman_params.py` (per-zone literature inputs: segment end points, cutoffs, coupling, b values, largest magnitude); `official_geometry_params.py` (raster lookup and per-zone geometry); `slab_contours.py` (SLAB raster to depth contours); `slab2_catalog.py` (the 27 regions); `bird_convergence.py` (plate convergence per along-strike column); `segmentation.py` (LEVEL 0 boundaries); `plot_meshes.py` (mesh figures); `hs_official_compare.py` (step 8's comparison of slip fields).
- **`pyptha/`**: `contour_discretisation.py` (contours to unit sources), `unit_sources.py` (cell statistics), `mesh_edit.py` (the `-edited` column rules), `mesh_file.py` (external meshes), `scaling.py` (rupture scaling), `events.py` (uniform-slip ruptures), `rates.py` (Gutenberg-Richter rates and the logic-tree machinery), `logic_tree.py` (LEVELs 0, 3, 4, 5), `moment_balance.py` (moment balance and the coupling prior), `stochastic_slip.py` (heterogeneous slip), `hs_vaus_rates.py` (rates of those scenarios), `variable_mu.py` (depth-dependent rigidity), `okada.py` (surface deformation), `grid_cache.py` (one saved mesh per run). `pyptha/examples/` holds four short scripts that call the engine directly.
- **`logic_tree/`**: `run_logic_tree.py` reads one input JSON and writes the CSV outputs. `gr_curve_and_convergence.py` summarises the posterior curve and the convergence per zone.
- **`templates/`**: `step1_fetch_slab2.py.tmpl` to `step9_report.py.tmpl`, `step_total.py.tmpl`, `RUN.html.tmpl`. `generate.py` replaces the `{{...}}` placeholders and writes them into the zone folder.
- **`data/`**: `bird/` (the plate-boundary tables, `ptha18_convergence_sources.csv` and their README), `ptha18_peak_slip_quantile_weights.csv`, `slab2/` (eight SLAB2.0 rasters: `cal`, `car`, `hel`, `ker`, `mak`, `png`, `puy`, `sum`), `slab1/` (empty; SLAB1.0 rasters are downloaded on demand).
- **`validation/`**: `validate_against_ptha18.py`, `validate_variable_mu.py`, `validate_mesh_file.py`, `recover_peak_slip_weights.py`, and `ptha18_reference/` (PTHA18's contours, unit-source grids and statistics files).

## The nine steps

`generate.py` writes the step scripts into `<folder>/steps/`. Run them from the repository root, one at a time, or all at once with `step_total.py`. Paths below are relative to the zone folder.

| Step | Script | Reads | Writes |
|---|---|---|---|
| 1 | `step1_fetch_slab2.py` | SLAB raster (the zone folder's `data/slab2/` or `data/slab1/`, then `ptha_python/data/`; downloaded if in neither) | `inputs/geometry/<zone>_slab2_contours.shp`, `data/step1_contours_info.json`, `figures/step1_trench_and_contours.png`. With `--mesh-file` only the JSON |
| 2 | `step2_build_grid.py` | the contours (or the external mesh file) | `data/slab2/unit_source_grid.npy` and `unit_source_grid.meta.json` (the saved mesh), `figures/scratch_mesh_map.png` |
| 3 | `step3_convergence.py` | the mesh; `data/bird/` | `data/convergence.txt`, `data/convergence_per_cell.npy`, `data/convergence_per_column.json`, `data/plate_boundary_change.json` (only if it finds a change) |
| 4 | `step4_fetch_gcmt.py` | the GCMT NDK files (downloaded) | `data/gcmt/*.ndk`, `data/gcmt/gcmt_catalogue_1976_present.csv` |
| 5 | `step5_subset_gcmt.py` | the NDK files of step 4; the mesh | `data/gcmt/<zone>_gcmt_subset.csv`; with `--segmented true` also one file per segment and `data/gcmt/segments.json` |
| 6 | `step6_write_input.py` | steps 1 to 5; the constants of `step2_build_grid.py`; `lib/berryman_params.py` | `inputs/input_<zone>_scratch.json` (kept if it exists; `--force` rebuilds it) |
| 7 | `step7_run.py` | the input JSON | `outputs/`, replaced as a whole (summary, rate curves, percentiles, branches, scenario rates, `figures/`) |
| 7b | `step7b_stochastic_slip.py` | the input JSON; step 7's scenarios | `outputs/hs_slip_fields/`, `outputs/vaus_slip_fields/` |
| 7c | `step7c_variable_mu.py` | the input JSON; steps 7 and 7b | `inputs/input_<zone>_scratch_variable_mu.json`, the `*_variable_mu*` files in `outputs/`, extra columns in the 7b summaries |
| 8 | `step8_official.py` | `official_ptha_data/` and the rptha data (needs R) | `inputs/input_<zone>_official.json`, `outputs_official/`, `outputs_official_status.txt` |
| 9 | `step9_report.py` | everything above | `report.html`, `outputs/exports/` (the CSV files behind it), mesh figures in `figures/`, and `report_unsegmented.html` on a segmented run |

Steps 7b, 7c and 8 are optional. `step_total.py` reports a failed optional step and continues, and it stops at the first failure among the required steps. It runs step 7c only for a folder made with `--variable-mu on`, and step 8 fails with a message when R is missing.

Flags of the scripts:

| Script | Flag | Effect |
|---|---|---|
| `step_total.py` | `--skip-official` | do not run step 8 |
| `step_total.py` | `--skip-hs` | do not run steps 7b and 7c |
| `step_total.py` | `--auto-clip` | pass `--auto-clip` to step 1 and `--force` to step 6 (needs one earlier run) |
| `step_total.py` | `--download-official-hs-vaus` | pass `--download-hs` to step 8 |
| `step1_fetch_slab2.py` | `--auto-clip` | cut the columns that step 3 found on another plate boundary |
| `step6_write_input.py` | `--force` | overwrite an existing (possibly hand-edited) input JSON |
| `step7b_stochastic_slip.py` | `--max-placements-per-mw N`, `--max-fields-per-mw N` | cap the number of heterogeneous-slip fields (faster, less exact) |
| `step7c_variable_mu.py` | `--force` | run in a folder made with `--variable-mu off` |
| `step8_official.py` | `--download-hs` | download PTHA18's HS and VAUS catalogues and compare them with step 7b |

Options that change a step's behaviour are documented in `html/docs/reference.html`. The generator options are listed in `../README.md`.

## How a zone folder is organised

```
<folder>/
  README.md          the zone's own README, written by generate.py
  RUN.html           the commands that created the folder
  report.html        the report (step 9)
  steps/             the step scripts, and step_total.py
  inputs/            the input JSON files and inputs/geometry/ (contours, or the copied external mesh)
  data/              intermediate data: slab2/ (mesh), gcmt/, convergence files
  figures/           the figures of steps 1 and 2
  outputs/           the engine's results and exports/ (the CSV files behind the report)
  outputs_official/  PTHA18's official run (step 8), only if it ran
```

`<folder>` is the name given by `--folder`. If that folder exists already, `generate.py` adds a suffix such as `_2`, so an existing folder is never overwritten.

## Adding a zone that is not in `lib/berryman_params.py`

The zone must be one that `--ptha false` can build. Those are the zones with a Berryman et al. (2015) mapping.

1. Check that its SLAB2.0 region is in `list_regions.py`, or find the ScienceBase id to pass with `--id`.
2. In `lib/berryman_params.py`, add the zone to `ZONE_TO_BERRYMAN_SEGMENTS`, listing the Table 3.1 rows that make up the zone. Add it to `ZONE_DOMINANT_SEGMENT` with the row that sets its seismogenic cutoff. Those rows come from `BERRYMAN_ROWS`, so the names must match exactly.
3. Optional: add entries to `ZONE_WHOLE_MARGIN` (the zone's whole-margin row, if Table 3.1 has one), `ZONE_TRENCH_ENDS` (where a part of a SLAB region is cut) and `ZONE_DEFAULT_CLIP`.
4. Run `python ptha_python/generate.py <name> --zone <zone> --folder <zone>`. Step 1 prints the geometry inputs it recovered. Check them before running the rest.

Without the mapping, step 1 stops with a message that the zone "has no Berryman et al. (2015) segment mapping, so there is no seismogenic cutoff": the seismogenic cutoff comes only from the mapping. `generate.py` prints the geometry inputs it recovered, and says which are missing, before the folder is built.

## Engine and validation

- **Engine**: `pyptha/` has no dependency on the zone folders. `logic_tree/run_logic_tree.py` is the entry point the steps call with one input JSON. The output files are listed in `html/docs/reference.html`.
- **Validation**: `validation/validate_against_ptha18.py` compares the mesh, engine, event rates, GCMT selection and segments with PTHA18's files. `validate_variable_mu.py` checks step 7c, `validate_mesh_file.py` checks `--mesh-file`. `recover_peak_slip_weights.py` recovers PTHA18's slip weights. The checks and their commands are in `html/docs/validation.html`.
- **Code correspondence**: `html/docs/code_map.html` maps each Python module to the rptha file it ports.
