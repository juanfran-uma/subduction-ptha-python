# Example results

Four finished runs, two around the Pacific and two in the Mediterranean, so that you can run the same zone yourself and compare your numbers with these. Open each `report.html` in a browser to see the whole run, and `mesh/mesh_3d.html` to see its mesh in 3D with the SLAB contours it was cut from.

| Example | Zone | Mesh (rows x columns) | How it was generated |
|---|---|---|---|
| `ptha18_compare/kermadectonga2` | Kermadec-Tonga (Pacific), modelled by PTHA18 | 3 x 73 | `generate.py kermadec --zone kermadectonga2 --ptha false` |
| `ptha18_compare/cascadia` | Cascadia (Pacific), modelled by PTHA18 | 3 x 23 | `generate.py cascadia --zone cascadia --ptha false` |
| `new_zones/calabria2` | Calabrian arc (Mediterranean), not modelled by PTHA18 | 6 x 15 | `generate.py calabria --zone calabria2 --ptha false --rupture-size local` |
| `new_zones/hellenic2` | Hellenic arc (Mediterranean), not modelled by PTHA18 | 7 x 30 | `generate.py hellenic --zone hellenic2 --ptha false --rupture-size local` |

The two `ptha18_compare` zones also have `outputs_official/`: PTHA18's own inputs for the zone, run through the same engine (step 8, which needs R), so that the report can put the two results side by side. The two Mediterranean zones are ones PTHA18 never modelled, so there is no official run to compare with.

## Compare your run with these

Generate the same zone in the repository root (the folder name is free, here it is the zone name) and run it. For example, for Calabria:

```
python ptha_python/generate.py calabria --zone calabria2 --folder calabria2 --ptha false --rupture-size local
python calabria2/steps/step_total.py --skip-official
```

Then compare the files in `calabria2/outputs/` with the ones in `examples/new_zones/calabria2/outputs/`. We did this from a fresh clone of the repository (the other three examples were not repeated): `rate_curves.csv`, `exceedance_rate_percentiles.csv`, `scenario_rates_calabria2.csv`, `logic_tree_summary.csv` and `mesh/mesh.dat` came out with **the same contents**. They are not byte for byte equal on Windows, because the repository stores line endings as LF and a Windows run writes CRLF, so compare ignoring them (for example `diff --strip-trailing-cr a b`, or load both CSV files in pandas). If the numbers differ, something in the setup is different (the options of `generate.py`, a package version, or a newer GCMT month in step 4).

The `generate.py` options of every example are also written in its `RUN.html`.

## What to expect

Exceedance rate per year of earthquakes of at least Mw 7.5, 8.0, 8.5 and 9.0 (the `exceedance_rate_unsegmented` column of `outputs/rate_curves.csv`):

| Example | Mw 7.5 | Mw 8.0 | Mw 8.5 | Mw 9.0 |
|---|---|---|---|---|
| kermadectonga2 | 0.118 | 0.0402 | 0.0127 | 0.00285 |
| kermadectonga2, PTHA18 inputs (`outputs_official/`) | 0.118 | 0.0398 | 0.0123 | 0.00259 |
| cascadia | 0.0276 | 0.00958 | 0.00318 | 0.000753 |
| cascadia, PTHA18 inputs (`outputs_official/`) | 0.0268 | 0.00916 | 0.00297 | 0.000641 |
| calabria2 | 0.00700 | 0.00158 | 0.000264 | 0.0000205 |
| hellenic2 | 0.0208 | 0.00767 | 0.00271 | 0.000802 |

## What is in each folder, and what is left out

Each folder is what `step_total.py` writes: `report.html`, `RUN.html`, `README.md`, `steps/` (the nine scripts), `inputs/`, `data/`, `outputs/`, `figures/` and `mesh/`. To keep the repository small, three heavy kinds of files are **not** included. The steps write them again when you run the folder:

- `data/gcmt/*.ndk` and `data/gcmt/gcmt_catalogue_1976_present.csv`: the GCMT catalogue (about 33 MB), downloaded and parsed again by step 4. The subset of earthquakes that step 5 kept, `data/gcmt/<zone>_gcmt_subset.csv`, is included.
- `outputs/hs_slip_fields/` and `outputs/vaus_slip_fields/`: the slip fields of step 7b.
- `outputs/logic_tree_branches_<zone>.csv` (and the same file in `outputs_official/`): one row per logic-tree branch, about 7 MB each.

The rest of the 26 runs behind these examples are not in the repository. To get another zone, see `ptha_python/html/docs/recipes.html`.
