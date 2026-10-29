"""Generate a complete, self-contained PTHA18-style source-zone example for
any subduction zone covered by SLAB2.0 / SLAB1.0 (or for an external mesh).

The generated folder holds the nine steps as plain scripts (steps/), the
inputs and outputs of each, a RUN.html with the exact commands used, and a
README.md. Running steps/step_total.py builds the geometry, convergence,
seismicity, logic tree, scenario ruptures, rates and the HTML report.

Usage
-----
    python ptha_python/generate.py <name> [--zone Z] [--folder F]
        [--ptha true|false] [--segmented true|false]
        [--segment-boundaries berryman|bird]
        [--columns trench|average|trench-edited|average-edited]
        [--cell-size mean|strasser] [--cell-k K]
        [--rupture-size rptha|local] [--variable-mu on|off]
        [--convergence bird|bird-griffin] [--discretizer optimal|lm|mid]
        [--clip LON_MIN,LON_MAX,LAT_MIN,LAT_MAX]
        [--mesh-file FILE [--mesh-depth-units auto|m|km]]

Examples
--------
    python ptha_python/generate.py calabria --zone calabria2 --folder calabria2 --rupture-size local
    python ptha_python/generate.py kermadec --zone kermadectonga2 --folder kermadectonga2 --segmented true
    python ptha_python/generate.py caribbean --zone antilles2 --folder antilles2 --clip 297,307,5,18.5
    python ptha_python/generate.py alaska --zone alaskaaleutians --folder alaskaaleutians --mesh-file inputs_meshes/alaskaaleutians_quadrilateral_coors.dat
    python calabria2/steps/step_total.py

`python ptha_python/list_regions.py` lists the SLAB regions and
`html/docs/zones.html` the zone names. Every option is described in
`html/docs/reference.html`.
"""

import argparse
import csv
import html
import os
import re
import shlex
import sys
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
TEMPLATES_DIR = os.path.join(HERE, "templates")

# This package's own names, substituted into every generated step so an
# example runs the engine of the package that generated it (see the mapping
# in main()).
PKG = os.path.basename(HERE)
PYPKG = "pyptha"
LTPKG = "logic_tree"
# rptha/ lives at the repository root
_RPTHA_ROOT = ROOT
SZP = os.path.join(_RPTHA_ROOT, "rptha", "R", "examples", "austptha_template", "DATA",
                   "SOURCEZONE_PARAMETERS", "sourcezone_parameters.csv")

sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))  # helper modules
from slab2_catalog import by_id, find_region, guess_ptha18_zones  # noqa: E402
from official_geometry_params import width_km, official_area_and_dip  # noqa: E402
from berryman_params import (berryman_cutoff_km, berryman_coupling,  # noqa: E402
                             berryman_b_anchor, berryman_mw_max_observed,
                             ZONE_TO_BERRYMAN_SEGMENTS, ZONE_DEFAULT_CLIP,
                             ZONE_TRENCH_ENDS)

TEMPLATE_FILES = [
    "step1_fetch_slab2.py",
    "step2_build_grid.py",
    "step3_convergence.py",
    "step4_fetch_gcmt.py",
    "step5_subset_gcmt.py",
    "step6_write_input.py",
    "step7_run.py",
    "step7b_stochastic_slip.py",
    "step7c_variable_mu.py",
    "step8_official.py",
    "step9_report.py",
    "step_total.py",
]


def unsegmented_zone_names():
    names = []
    with open(SZP, newline="") as f:
        for row in csv.DictReader(f):
            if not row.get("segment_name", "").strip():
                names.append(row["sourcename"].strip())
    return names


def zone_exists(name):
    return name in unsegmented_zone_names()


def zone_row_weight(name):
    """row_weight for `name`'s unsegmented row, or None if not found."""
    with open(SZP, newline="") as f:
        for raw_row in csv.DictReader(f):
            row = {k.strip(): (v.strip() if v else "") for k, v in raw_row.items()}
            if row["sourcename"] == name and not row.get("segment_name"):
                try:
                    return float(row["row_weight"])
                except (KeyError, ValueError):
                    return None
    return None


def resolve_region(name, explicit_id):
    if explicit_id:
        region = by_id(explicit_id)
        if region is None:
            raise SystemExit(
                f"'{explicit_id}' is not one of the 27 known SLAB2.0 "
                f"ScienceBase item ids.\nRun "
                f"python ptha_python/list_regions.py "
                f"to see the full list.")
        return region

    hits = find_region(name)
    if len(hits) == 1:
        return hits[0]
    if len(hits) == 0:
        raise SystemExit(
            f"no SLAB2.0 region matches '{name}'.\nRun "
            f"python ptha_python/list_regions.py "
            f"to see the 27 available regions, then either use one of "
            f"their names/prefixes or pass --id directly.")
    lines = "\n".join(f"  {r['title']}  (prefix '{r['prefix']}', "
                      f"id {r['id']})" for r in hits)
    raise SystemExit(
        f"'{name}' matches more than one SLAB2.0 region:\n{lines}\n"
        f"Be more specific, or pass --id <scibase_item_id> directly.")


def warn_if_retired(zone_name):
    if zone_row_weight(zone_name) == 0.0:
        print(f"note: '{zone_name}' is a retired PTHA18 zone "
             f"(row_weight = 0), step 8 will find no logic tree to compare "
             f"against.")


def resolve_zone_from_berryman(explicit_zone, name, region_title):
    """--ptha false: the zone name, from this package's own Berryman et al.
    (2015) mapping only (lib/berryman_params.py); sourcezone_parameters.csv
    is never read."""
    if explicit_zone:
        if explicit_zone not in ZONE_TO_BERRYMAN_SEGMENTS:
            print(f"  note: '{explicit_zone}' has no Berryman et al. (2015) "
                  f"mapping in lib/berryman_params.py, so step 1 will stop. "
                  f"Add the zone to ZONE_TO_BERRYMAN_SEGMENTS and "
                  f"ZONE_DOMINANT_SEGMENT there first (all valid names: "
                  f"html/docs/zones.html).")
        return explicit_zone
    zone_names = sorted(ZONE_TO_BERRYMAN_SEGMENTS)
    squashed = name.lower().replace(" ", "").replace("-", "").replace("_", "")
    if squashed in zone_names:
        return squashed
    guesses = guess_ptha18_zones(region_title, zone_names)
    if len(guesses) == 1:
        return guesses[0]
    if not guesses:
        raise SystemExit(
            f"could not guess a zone name for '{region_title}' from "
            f"lib/berryman_params.py.\nPass --zone <name>.")
    raise SystemExit(
        f"'{region_title}' matches more than one zone:\n"
        + "\n".join(f"  {z}" for z in guesses) + "\nPass --zone <name> to pick one.")


def resolve_ptha18_zone(explicit_zone, name, region_title, allow_new_zone=False):
    if allow_new_zone:  # --ptha false
        return resolve_zone_from_berryman(explicit_zone, name, region_title)
    zone_names = unsegmented_zone_names()

    if explicit_zone:
        if not zone_exists(explicit_zone):
            if not allow_new_zone:
                raise SystemExit(
                    f"'{explicit_zone}' is not an unsegmented sourcename in "
                    f"sourcezone_parameters.csv.\nIf this is a real "
                    f"subduction zone SLAB2.0/SLAB1.0 covers but PTHA18 "
                    f"never modelled, pass --ptha false to build a run for "
                    f"it anyway (falls back to generic defaults for "
                    f"literature values PTHA18 does not have for this "
                    f"zone, see report_geometry_inputs' output).")
            berry = ("it has a Berryman et al. (2015) mapping"
                     if explicit_zone in ZONE_TO_BERRYMAN_SEGMENTS else
                     "and it has NO Berryman et al. (2015) mapping either")
            print(f"  note: '{explicit_zone}' is not a PTHA18 zone (no row "
                  f"in sourcezone_parameters.csv; {berry}), proceeding "
                  f"with --ptha false. step 8 "
                  f"(official comparison) will not be available for this "
                  f"zone; see the geometry-inputs report below for which "
                  f"other numbers use a generic default vs. need to be set "
                  f"by hand.")
            return explicit_zone
        warn_if_retired(explicit_zone)
        return explicit_zone

    guesses = guess_ptha18_zones(region_title, zone_names)
    active = [z for z in guesses if zone_row_weight(z) != 0.0]

    squashed = name.lower().replace(" ", "").replace("-", "").replace("_", "")
    if zone_exists(squashed):
        if zone_row_weight(squashed) == 0.0 and len(active) == 1:
            print(f"note: '{squashed}' is retired, using '{active[0]}'.")
            return active[0]
        warn_if_retired(squashed)
        return squashed

    if len(guesses) == 1:
        warn_if_retired(guesses[0])
        return guesses[0]
    if len(guesses) == 0:
        raise SystemExit(
            f"could not guess a PTHA18 zone name for '{region_title}'.\n"
            f"Pass --zone <exact_sourcename>.")
    if len(active) == 1:
        return active[0]

    lines = "\n".join(
        f"  {z}" + ("  (retired)" if zone_row_weight(z) == 0.0 else "")
        for z in guesses)
    raise SystemExit(
        f"'{region_title}' matches more than one PTHA18 zone:\n{lines}\n"
        f"Pass --zone <exact_sourcename> to pick one.")


CONVERGENCE_SOURCES_CSV = os.path.join(HERE, "data", "bird", "ptha18_convergence_sources.csv")


def convergence_source_note(zone, choice):
    """A note on where PTHA18 took `zone`'s convergence from, and what the
    two --convergence tables give on PTHA18's own mesh (measured by
    html/docs/figures_src/make_convergence_sources.py). Empty if the zone is
    not a PTHA18 subduction zone."""
    if not os.path.exists(CONVERGENCE_SOURCES_CSV):
        return ""
    with open(CONVERGENCE_SOURCES_CSV, newline="") as f:
        row = next((r for r in csv.DictReader(f) if r["zone"] == zone), None)
    if row is None:
        return ""
    both = (f"on PTHA18's mesh: bird {row['conv_bird']} mm/yr, bird-griffin "
            f"{row['conv_bird_griffin']} mm/yr")
    if row["use_bird_convergence"] == "0":
        return (f"  note: PTHA18 did not use Bird on '{zone}' but a constant "
                f"{row['constant_mm_per_yr']} mm/yr (set CONVERGENCE_OVERRIDE_MM_PER_YR "
                f"in step 3 for that); {both}")
    if row["source"] == "bird":
        return f"  note: PTHA18 took {zone}'s convergence from Bird's own steps; {both}"
    who = "Griffin's traces" if row["source"] == "griffin" else "Bird plus some of Griffin's or Davies' traces"
    hint = ("" if choice == "bird-griffin" else
            "; --convergence bird-griffin reproduces PTHA18's value")
    return (f"  NOTE: PTHA18 took {zone}'s convergence from {who} "
            f"({row['traces']}), not from Bird alone; {both}{hint}")


def pick_folder_name(requested):
    candidate = requested
    n = 1
    while os.path.exists(os.path.join(ROOT, candidate)):
        n += 1
        candidate = f"{requested}_{n}"
    return candidate


def render(text, mapping):
    for key, val in mapping.items():
        text = text.replace("{{" + key + "}}", val)
    leftover = re.findall(r"\{\{[A-Z0-9_]+\}\}", text)
    if leftover:
        raise SystemExit(
            f"template still has unresolved placeholders: {sorted(set(leftover))}")
    return text


def report_geometry_inputs(zone, use_ptha=True):
    """Print the recovered numbers before generating, so they are visible.

    If any of these looks wrong, it is better to find out now than after
    running the whole pipeline on it. With use_ptha False nothing of PTHA18
    is read.
    """
    print("\n  Geometry inputs recovered for this zone")
    print("  " + "-" * 62)
    cutoff, dominant_segment = berryman_cutoff_km(zone)
    if cutoff is None:
        print(f"  seismogenic cutoff   NOT AVAILABLE, '{zone}' has no "
              f"Berryman et al. (2015) segment mapping: add the zone to ZONE_TO_BERRYMAN_SEGMENTS and ZONE_DOMINANT_SEGMENT in lib/berryman_params.py (step 1 stops without it).")
    else:
        print(f"  seismogenic cutoff   {cutoff:6.1f} km below trench")
        print(f"                       (Berryman et al. 2015, dominant segment "
              f"'{dominant_segment}')")
    if not use_ptha:
        print("  unit-source width      50.0 km  (--ptha false: the general "
              "default, see step2_build_grid.py)")
    else:
        try:
            width = width_km(zone)
            print(f"  unit-source width    {width:6.1f} km  "
                  f"(sourcezone_parameters.csv, a resolution parameter, not a "
                  f"physical 'declared' value)")
        except SystemExit:
            print(f"  unit-source width    NOT AVAILABLE, '{zone}' has no row "
                  f"in sourcezone_parameters.csv; falling back to the general "
                  f"default of 50 km (see step2_build_grid.py's own docstring)")
    has_berryman = zone in ZONE_TO_BERRYMAN_SEGMENTS
    coupling = berryman_coupling(zone)
    b_anchor = berryman_b_anchor(zone)
    mw_max_observed = berryman_mw_max_observed(zone)
    coupling_src = "Berryman et al. 2015" if has_berryman else \
        "generic default, PTHA18's own uniform prior, no Berryman mapping for this zone"
    b_src = "Berryman et al. 2015" if has_berryman else \
        "generic default, no Berryman mapping for this zone"
    print(f"  coupling (min/pref/max)  {coupling}  ({coupling_src})")
    print(f"  b (min/pref/max)         {b_anchor}  ({b_src})")
    if mw_max_observed is None:
        print(f"  mw_max_observed          NOT AVAILABLE, '{zone}' has no "
              f"Berryman et al. (2015) segment mapping: add the zone to ZONE_TO_BERRYMAN_SEGMENTS and ZONE_DOMINANT_SEGMENT in lib/berryman_params.py (step 1 stops without it).")
    else:
        print(f"  mw_max_observed          {mw_max_observed:g}  (Berryman et al. 2015)")
    if not use_ptha:
        return  # no PTHA18 file is read with --ptha false
    off = official_area_and_dip(zone)
    if off:
        print(f"  official area        {off['area_km2']:10,.1f} km2   "
              f"<- step 2 will compare against this")
        print(f"  official mean dip    {off['mean_dip_deg']:6.3f} deg")
        print(f"  official mesh        {off['n_unit_sources']} unit sources, "
              f"depth {off['depth_min_km']:.1f}-{off['depth_max_km']:.1f} km")
    else:
        print(f"  official table not downloaded; step 2 will skip the "
              f"area/dip comparison.")
        print(f"  get it with: official_ptha_data/fetch_official_inputs.py "
              f"{zone}")


def main():
    ap = argparse.ArgumentParser(
        description="Generate a PTHA18-style source-zone example for any "
                    "SLAB2.0/SLAB1.0 zone.")
    ap.add_argument("name", help="zone name, e.g. cascadia, 'south america'")
    ap.add_argument("--zone", default=None,
                    help="zone name: with --ptha true an exact "
                         "sourcezone_parameters.csv sourcename, with --ptha "
                         "false a name in lib/berryman_params.py")
    ap.add_argument("--folder", default=None,
                    help="folder to create, relative to the repository root; "
                         "may be nested (e.g. runs/calabria2)")
    ap.add_argument("--id", default=None, help="ScienceBase item id")
    ap.add_argument("--ptha", default="false", choices=["true", "false"],
                    help="true: steps may read PTHA18's files, "
                         "per-zone unit-source width and scaling "
                         "(sourcezone_parameters.csv) and official down-dip "
                         "row count (the .nc file), and requires the zone to "
                         "be one PTHA18 modelled. false (default): steps 1-7 and 9 read "
                         "no PTHA18 or rptha file (50 km width, "
                         "geometry-derived row count, Strasser/30 GPa); only "
                         "step 8, the official run, does. Also allows --zone "
                         "to be a zone PTHA18 never modelled, provided it "
                         "has a berryman_params.py mapping. The Bird table is "
                         "chosen by --convergence, not by this flag.")
    ap.add_argument("--discretizer", default="optimal",
                    choices=["optimal", "lm", "mid"],
                    help="optimal (default): rptha's own mesh, kept unchanged "
                         "whenever it has no folded, collapsed or overlapping "
                         "cell; only a broken mesh is repaired (nearest-point "
                         "start, then end-aligned contours). "
                         "lm: rptha's method exactly as PTHA18 runs it, no "
                         "repairs. mid: rptha's older, deprecated method "
                         "(eq_spacing cuts); PTHA18 does not use it.")
    ap.add_argument("--columns", default="trench",
                    choices=["trench", "average", "trench-edited", "average-edited"],
                    help="how many unit sources along strike the optimal "
                         "mesh gets. trench (default): rptha's and PTHA18's "
                         "rule, trench length / 50 km. average: the average row "
                         "length / 50 km, for a zone whose mesh tapers far more "
                         "than any PTHA18 mesh (step 2 prints the taper; "
                         "PTHA18: at most 1.17), e.g. calabria2. "
                         "Not PTHA18's procedure. "
                         "trench-edited and average-edited: the trench "
                         "or the average mesh, kept unless "
                         "some cell reads a dip more than 8 degrees from the "
                         "dip of its own plane (slivers where an S shaped arc "
                         "makes rows run along strike, e.g. cyprus2); then "
                         "step 2 designs a replacement with the same columns "
                         "and rows on the same contours, with right-angled "
                         "cells, and uses it only if it is safe (fewer such "
                         "cells, no more cells under half the desired width, "
                         "no fold or overlap, area within 5%%). Not PTHA18's "
                         "procedure. Needs --discretizer optimal and "
                         "--cell-size mean, not --mesh-file. The per-column "
                         "Bird convergence is matched to the new mesh, so "
                         "rates are not comparable with the trench mesh's.")
    ap.add_argument("--cell-size", default="mean", choices=["mean", "strasser"],
                    help="how big the unit sources may be. "
                         "mean (default): rptha's and PTHA18's rule, cells of "
                         "~50 x 50 km ON AVERAGE (columns from the trench "
                         "length, rows from the mean down-dip length), so on a "
                         "fan or a zone of varying width some cells are much "
                         "larger (calabria2). "
                         "strasser: NO cell larger than the smallest rupture "
                         "the run builds (Strasser at Mmin = 7.2: 54 x 44 km) "
                         "divided by --cell-k; same mesher, more rows and "
                         "columns. Not PTHA18's procedure; 3-7 times more unit "
                         "sources. Needs --discretizer optimal.")
    ap.add_argument("--cell-k", type=float, default=1.5,
                    help="with --cell-size strasser: the cap is the "
                         "smallest rupture divided by this. 1.5 (default): "
                         "cells of at most 36 x 29 km, i.e. no larger than an "
                         "Mw 7.2 rupture one sigma smaller than Strasser's mean "
                         "in length and width; 2: at most 27 x 22 km.")
    ap.add_argument("--clip", default=None,
                    help="optional LON_MIN,LON_MAX,LAT_MIN,LAT_MAX window: "
                         "step 1 uses only this part of the SLAB raster. For "
                         "a SLAB region that covers more than the source zone "
                         "(e.g. only the Lesser Antilles part of 'Caribbean'). "
                         "Longitudes in the raster's own convention.")
    ap.add_argument("--mesh-file", default=None,
                    help="use this quadrilateral mesh instead of building "
                         "one from SLAB (steps 1-2). A text file with one line "
                         "per unit source and 12 numbers: its 4 corners as "
                         "lon lat depth, going round the cell (e.g. "
                         "inputs_meshes/alaskaaleutians_quadrilateral_coors.dat). It must be "
                         "a structured mesh (rows down dip x columns along "
                         "strike, cells sharing their corners). Copied into "
                         "the example's inputs/geometry/.")
    ap.add_argument("--mesh-depth-units", default="auto", choices=["auto", "m", "km"],
                    help="with --mesh-file: depth units of the file. auto "
                         "(default): metres if any |depth| > 200, else km; the "
                         "sign (negative down or positive down) is detected.")
    ap.add_argument("--variable-mu", default="off", choices=["on", "off"],
                    help="off (default): every rate uses the constant "
                         "rigidity of 30 GPa (PTHA18's rate_annual). on: "
                         "step_total also runs step 7c, which adds PTHA18's "
                         "variable shear modulus rates (variable_mu_*, the "
                         "ones its headline maps use with HS) next to the "
                         "constant ones. Needs step 7b.")
    ap.add_argument("--rupture-size", default="rptha", choices=["rptha", "local"],
                    help="how many cells each uniform-slip rupture gets. "
                         "rptha (default): rptha's and PTHA18's rule, one "
                         "block of cells per magnitude from the zone's mean "
                         "cell size. local: a block per placement, from the "
                         "real km of the cells there, so the rupture area "
                         "follows the scaling relation on a mesh whose cells "
                         "differ a lot (a fan such as calabria2); each "
                         "rupture's conditional probability is also "
                         "multiplied by q, the rupture-overlap correction. "
                         "Not PTHA18's procedure.")
    ap.add_argument("--convergence", default="bird", choices=["bird", "bird-griffin"],
                    help="which plate-boundary table step 3 reads, both in "
                         "ptha_python/data/bird/. bird (default): Bird "
                         "(2003)'s public catalogue. bird-griffin: the table "
                         "PTHA18 used, Bird plus Jonathan Griffin's source-zone "
                         "traces with their own plate rates. They agree on the "
                         "zones where PTHA18 used Bird's own steps and can "
                         "differ a lot where it used Griffin's traces (eastern "
                         "Indonesia, New Guinea): see "
                         "html/docs/convergence_sources.html.")
    ap.add_argument("--segmented", default="false", choices=["true", "false"],
                    help="false (default): run the zone unsegmented (LEVEL 0 "
                         "weight 1). true: also split the zone into "
                         "along-strike segments and run the full LEVEL 0 tree "
                         "(unsegmented branch 0.5, union of segments 0.5, the "
                         "segments summing within that branch). Boundaries "
                         "from --segment-boundaries, public data either way, "
                         "so this stays public-data-only under --ptha false. "
                         "See lib/segmentation.py.")
    ap.add_argument("--segment-boundaries", default="berryman",
                    choices=["berryman", "bird"],
                    help="with --segmented true: where the segments end. "
                         "berryman (default): at the trench end points "
                         "Berryman et al. (2015) Table 3.1 gives every "
                         "segment, the segments PTHA18 used (on PTHA18's own "
                         "meshes the boundaries land within 0-2 columns of "
                         "PTHA18's). bird: where Bird (2003)'s plate pair "
                         "changes along the trench (misses a segment Bird has "
                         "no plate for, e.g. Hikurangi).")
    args = ap.parse_args()

    # an external mesh (checked here, before anything is written)
    mesh_name = None
    if args.mesh_file:
        if not os.path.exists(args.mesh_file):
            raise SystemExit(f"--mesh-file: {args.mesh_file} not found")
        from pyptha import mesh_file as _mesh_file
        _mesh_file.read_quadrilateral_mesh(args.mesh_file, depth_units=args.mesh_depth_units)
        mesh_name = os.path.basename(args.mesh_file)
    if args.cell_size == "strasser":
        if args.mesh_file:
            raise SystemExit("--cell-size strasser builds the mesh from SLAB; "
                             "it cannot be combined with --mesh-file")
        if args.discretizer != "optimal":
            raise SystemExit("--cell-size strasser needs --discretizer optimal")
        if not args.cell_k > 0:
            raise SystemExit("--cell-k must be positive")
    trench_edited = args.columns.endswith("-edited")          # the mesh, redesigned only if it misreads the dip
    column_rule = args.columns.replace("-edited", "")
    if trench_edited and (args.mesh_file or args.discretizer != "optimal" or args.cell_size != "mean"):
        raise SystemExit(f"--columns {args.columns} needs --discretizer optimal and --cell-size mean, "
                         "and it cannot be combined with --mesh-file")
    segmented = args.segmented == "true"
    clip = None
    if args.clip:
        try:
            clip = tuple(float(v) for v in args.clip.split(","))
        except ValueError:
            clip = ()
        if len(clip) != 4 or clip[0] >= clip[1] or clip[2] >= clip[3]:
            raise SystemExit("--clip must be LON_MIN,LON_MAX,LAT_MIN,LAT_MAX "
                             "with LON_MIN < LON_MAX and LAT_MIN < LAT_MAX")

    use_ptha = args.ptha == "true"

    region = resolve_region(args.name, args.id)
    ptha18_zone = resolve_ptha18_zone(args.zone, args.name, region["title"],
                                      allow_new_zone=not use_ptha)
    if clip is None and ptha18_zone in ZONE_DEFAULT_CLIP:
        clip = ZONE_DEFAULT_CLIP[ptha18_zone]
    # a zone that is one part of a SLAB region is cut at Berryman's end points
    trench_ends = ZONE_TRENCH_ENDS.get(ptha18_zone)

    requested = args.folder or re.sub(r"[^a-z0-9_]+", "_", args.name.lower()).strip("_")
    # --folder may be nested (e.g. runs/calabria2): the steps then climb
    # one more level to find the package (ROOT_UP), and the run and figure
    # labels use the last component only (NAME)
    requested = requested.replace("\\", "/").strip("/")
    folder = pick_folder_name(requested)
    example_dir = os.path.join(ROOT, folder)

    print("=" * 70)
    print(f"Generating '{folder}/' for {region['title']}")
    print(f"  SLAB2.0 prefix '{region['prefix']}', PTHA18 zone '{ptha18_zone}'")
    print(f"  --ptha={args.ptha}  ({'may read PTHA18 files (width, row count, scaling, '
          'comparisons)' if use_ptha else 'steps 1-7 and 9 read no PTHA18 or rptha file; only step 8 does'})")
    discretizer_note = {
        "optimal": "rptha's own mesh; repaired only if it has folded, "
                   "collapsed or overlapping cells",
        "lm": "rptha's multigrid Levenberg-Marquardt method exactly as "
              "PTHA18 ran it, no repairs",
        "mid": "rptha's older, deprecated eq_spacing method",
    }[args.discretizer]
    print(f"  --discretizer={args.discretizer}  ({discretizer_note})")
    if args.cell_size == "strasser":
        print(f"  --cell-size=strasser --cell-k={args.cell_k:g}  (no unit source larger "
              f"than the Mmin rupture / {args.cell_k:g}; not PTHA18's procedure)")
    if trench_edited:
        print(f"  --columns={args.columns}  (the {column_rule} mesh is kept unless a cell misreads the dip by more "
              "than 8 deg; then step 2 designs a safe replacement; not PTHA18's procedure)")
    print(f"  --convergence={args.convergence}  ("
          + ("Bird (2003)'s public catalogue" if args.convergence == "bird" else
             "Bird + Griffin's traces, the table PTHA18 used") + ")")
    note = convergence_source_note(ptha18_zone, args.convergence)
    if note:
        print(note)
    if clip:
        print(f"  --clip={clip}  (step 1 uses only this lon/lat window of "
              f"the SLAB raster)")
    if trench_ends and not mesh_name:
        print(f"  zone ends: step 1 cuts the zone along the down-dip lines "
              f"through Berryman et al. 2015's end points {trench_ends}")
    print("=" * 70)

    report_geometry_inputs(ptha18_zone, use_ptha)

    mapping = {
        "FOLDER": folder,
        "NAME": folder.rsplit("/", 1)[-1],
        # the exact command line, shown at the top of RUN.html
        "GENERATE_CMD": html.escape(
            f"python {PKG}/generate.py "
            + " ".join(shlex.quote(a) for a in sys.argv[1:])),
        "ROOT_UP": ", ".join(['".."'] * (folder.count("/") + 1)),
        "ZONE": ptha18_zone,
        "PTHA18_ZONE": ptha18_zone,
        # the generated steps import the engine from the package that
        # generated them
        "PKG": PKG,
        "PYPKG": PYPKG,
        "LTPKG": LTPKG,
        "SEGMENTED": "True" if segmented else "False",
        "SEGMENT_BOUNDARIES": args.segment_boundaries,
        "TITLE": region["title"],
        "PREFIX": region["prefix"],
        "SCIBASE_ITEM": region["id"],
        "USE_PTHA": "True" if use_ptha else "False",
        "DISCRETIZER": repr(args.discretizer),
        "COLUMN_RULE": repr(column_rule),
        "CELL_SIZE": repr(args.cell_size),
        "CELL_K": repr(args.cell_k),
        "TRENCH_EDITED": "True" if trench_edited else "False",
        "RUPTURE_SIZE": repr(args.rupture_size),
        "CONVERGENCE": repr(args.convergence),
        # Manual-download fallback links (see step1_fetch_slab2.py.tmpl's
        # "If the automatic download fails" section). region["title"] is
        # this region's full ScienceBase title (e.g. "Kamchatka-Kuril
        # Islands-Japan"), the same string download_depth_grid() itself
        # searches on, so the link actually finds this zone's item.
        "SLAB2_REGION_PHRASE": urllib.parse.quote(region["title"]),
        "SLAB1_REGION_PREFIX": region["prefix"],
        "CLIP_BBOX": repr(clip) if clip else "None",
        "TRENCH_ENDS": repr(trench_ends) if trench_ends else "None",
        "VARIABLE_MU": "True" if args.variable_mu == "on" else "False",
        "MESH_FILE": repr(mesh_name) if mesh_name else "None",
        "MESH_DEPTH_UNITS": repr(args.mesh_depth_units),
    }

    for sub in ("steps", "data/slab1", "data/slab2", "data/gcmt", "inputs/geometry", "outputs"):
        os.makedirs(os.path.join(example_dir, sub), exist_ok=True)

    if mesh_name:
        # the external mesh, copied so the example is self-contained
        import shutil
        shutil.copyfile(args.mesh_file, os.path.join(example_dir, "inputs", "geometry", mesh_name))
        print(f"  copied the mesh -> inputs/geometry/{mesh_name}")

    print()
    for fname in TEMPLATE_FILES:
        with open(os.path.join(TEMPLATES_DIR, fname + ".tmpl"),
                  "r", encoding="utf-8") as f:
            text = f.read()
        with open(os.path.join(example_dir, "steps", fname),
                  "w", encoding="utf-8") as f:
            f.write(render(text, mapping))
        print(f"  wrote steps/{fname}")

    with open(os.path.join(TEMPLATES_DIR, "RUN.html.tmpl"),
              "r", encoding="utf-8") as f:
        run_html = render(f.read(), mapping)
    with open(os.path.join(example_dir, "RUN.html"), "w", encoding="utf-8") as f:
        f.write(run_html)
    print("  wrote RUN.html")

    cutoff, dominant_segment = berryman_cutoff_km(ptha18_zone)
    coupling = berryman_coupling(ptha18_zone)
    b_anchor = berryman_b_anchor(ptha18_zone)
    mw_max_observed = berryman_mw_max_observed(ptha18_zone)
    cutoff_str = (f"{cutoff:g} km below the trench (Berryman et al. 2015 Table 3.1, "
                  f"dominant segment '{dominant_segment}')" if cutoff is not None else
                  "**NOT AVAILABLE**: no Berryman et al. (2015) mapping "
                  "for this zone; add it to `lib/berryman_params.py` "
                  "(step 1 stops without it)")
    mw_max_str = (f"{mw_max_observed:g}" if mw_max_observed is not None else
                  "**NOT AVAILABLE**: no Berryman et al. (2015) mapping "
                  "for this zone; add it to `lib/berryman_params.py` "
                  "(step 1 stops without it)")
    if not use_ptha:
        width_str = "50 km (the general default)"
    else:
        try:
            width_str = f"{width_km(ptha18_zone):g} km (sourcezone_parameters.csv)"
        except SystemExit:
            width_str = "50 km (generic default, no row in sourcezone_parameters.csv for this zone)"

    if mesh_name:
        geometry_line = (f"external mesh: the unit sources are the cells of "
                         f"`inputs/geometry/{mesh_name}`, used as given "
                         f"(`--mesh-file`); steps 1-2 do not use SLAB")
        mesh_line = "not applicable (external mesh)"
    else:
        geometry_line = ("SLAB: the trench is the up-dip edge of the raster, "
                         "contours every 5 km below the nearby trench"
                         + (f"; zone cut at Berryman et al. (2015)'s end points {trench_ends}"
                            if trench_ends else "")
                         + (f"; clipped to the lon/lat window {clip}" if clip else ""))
        if args.cell_size == "strasser":
            mesh_line = (f"`--cell-size strasser --cell-k {args.cell_k:g}`: no unit source larger "
                         f"than the Mmin rupture / {args.cell_k:g}")
        else:
            mesh_line = (f"`--columns {args.columns}`, `--discretizer {args.discretizer}`, "
                         f"cells of about 50 km")
    rupture_line = ("local: each rupture is sized from the real km of the cells where it sits, "
                    "with the overlap correction q" if args.rupture_size == "local"
                    else "rptha's rule: one block of cells per magnitude")
    if segmented:
        segment_line = ("on (LEVEL 0), boundaries from "
                        + ("Berryman et al. (2015) Table 3.1" if args.segment_boundaries == "berryman"
                           else "Bird (2003) plate pairs"))
    else:
        segment_line = "off (LEVEL 0 weight 1)"
    mu_line = ("constant 30 GPa, plus PTHA18's variable shear modulus rates (step 7c)"
               if args.variable_mu == "on" else "constant 30 GPa (step 7c only runs with --force)")
    conv_line = ("Bird (2003)" if args.convergence == "bird"
                 else "Bird (2003) plus Griffin's traces (the table PTHA18 used)")
    ptha_line = ("may be read (`--ptha true`)" if use_ptha
                 else "steps 1-7 and 9 read none (`--ptha false`); only step 8 does")
    readme = f"""# {region['title']}

Generated by `{PKG}/generate.py` for the SLAB2.0/SLAB1.0 region
"{region['title']}" (ScienceBase item {region['id']}, prefix
"{region['prefix']}"), matched to the source zone `{ptha18_zone}`.

## Settings of this run

| | |
|---|---|
| Geometry | {geometry_line} |
| Mesh | {mesh_line} |
| Rupture size | {rupture_line} |
| Segmentation | {segment_line} |
| Shear modulus | {mu_line} |
| Plate convergence | {conv_line} |
| PTHA18 files | {ptha_line} |

## Inputs taken from the literature

| | |
|---|---|
| Seismogenic cutoff | {cutoff_str} |
| Coupling (min/pref/max) | {coupling} |
| b-value (min/pref/max) | {b_anchor} |
| mw_max_observed | {mw_max_str} |
| Unit-source width | {width_str} |

## Run

One step at a time, from the repository root:

    python {folder}/steps/step1_fetch_slab2.py       # SLAB contours (or reads the external mesh)
    python {folder}/steps/step2_build_grid.py        # unit-source mesh
    python {folder}/steps/step3_convergence.py       # plate convergence per column
    python {folder}/steps/step4_fetch_gcmt.py        # GCMT catalogue
    python {folder}/steps/step5_subset_gcmt.py       # events of this zone
    python {folder}/steps/step6_write_input.py       # the input JSON
    python {folder}/steps/step7_run.py               # logic tree, ruptures and rates
    python {folder}/steps/step7b_stochastic_slip.py  # HS and VAUS variants (optional)
    python {folder}/steps/step7c_variable_mu.py      # variable shear modulus rates (after 7b)
    python {folder}/steps/step8_official.py          # PTHA18's official run (optional, needs R)
    python {folder}/steps/step9_report.py            # report.html

Or all in one go:

    python {folder}/steps/step_total.py

`RUN.html` lists the exact commands this folder was made with.
"""
    with open(os.path.join(example_dir, "README.md"), "w", encoding="utf-8") as f:
        f.write(readme)
    print("  wrote README.md")

    print("\nDone. Start with:")
    print(f"  python {folder}/steps/step1_fetch_slab2.py")


if __name__ == "__main__":
    sys.exit(main())
