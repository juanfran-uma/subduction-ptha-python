"""weight_with_nonzero_rate and rates of a zone against PTHA18's published HS table.

Runs the engine on PTHA18's own input of a zone, with PTHA18's published HS table as the
variable-rigidity catalogue, then compares, parent by parent, the engine's weight_with_nonzero_rate,
variable_mu_weight_with_nonzero_rate and rates (constant and variable rigidity) with the values in
that table, through its uniform_event_row. Each HS row carries its parent's weights, and the rates
of a parent's HS rows add up to the parent's rate.

  --mode official   (default) PTHA18's full official input (official_segmented_input of
                    validate_against_ptha18.py): its Bird convergence and, for a zone with segments,
                    the unsegmented model plus every segment, mixed as PTHA18 publishes it
                    (0.5 / 0.5 for kermadectonga2). A zone without segments is just its one model.
                    This is what the published table contains.
  --mode plain      the plain official input inputs/input_<zone>.json (one unsegmented model, no
                    Bird convergence profile): on a zone with segments it can only match where the
                    segments agree with it, which shows why the mix matters.

Needs, from the repo root: official_ptha_data/trees (extracted once with
official_ptha_data/fetch_official_inputs.py), inputs/ and PTHA18's HS table of the zone.

    python ptha_python/validation/validate_weight_with_nonzero_rate.py \\
        --hs-nc nc_files/all_stochastic_slip_earthquake_events_kermadectonga2.nc \\
        --zone kermadectonga2 --mode official --out html/data/kermadectonga2_validation_official.json
"""

import argparse
import json
import os
import shutil
import subprocess
import sys

import netCDF4
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(PKG, os.pardir))
sys.path.insert(0, HERE)
sys.path.insert(0, PKG)
sys.path.insert(0, os.path.join(PKG, "lib"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hs-nc", required=True, help="PTHA18's all_stochastic_slip_earthquake_events_<zone>.nc")
    ap.add_argument("--zone", default="kermadectonga2")
    ap.add_argument("--mode", choices=["official", "plain"], default="official")
    ap.add_argument("--out", required=True, help="json with the results (an .npz with the arrays goes next to it)")
    a = ap.parse_args()
    z = a.zone
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    if a.mode == "official":
        import validate_against_ptha18 as v
        cfg = v.official_segmented_input(z)
    else:
        cfg = json.load(open(os.path.join(ROOT, "inputs", f"input_{z}.json")))
    cfg["run_name"] = f"{z}_weight_validation_{a.mode}"
    cfg["percentiles"]["N"] = 2000          # LEVEL 5 is not compared here
    cfg["events"].update({"scaling_relation": "Strasser", "shear_modulus_Pa": 3.0e10})
    cfg["geometry"]["official_statistics_nc"] = os.path.join(ROOT, cfg["geometry"]["official_statistics_nc"])
    cfg["variable_shear_modulus"] = {"hs_events": os.path.abspath(a.hs_nc)}
    tmp = os.path.join(HERE, f"_input_{z}_weights_{a.mode}.json")
    json.dump(cfg, open(tmp, "w"), indent=1)
    out = os.path.join(ROOT, "runs", "python", cfg["run_name"])
    # the zone-level files hold the 0.5/0.5 mix when there are segments, the plain ones otherwise
    suffix = "_source_zone" if cfg["segments"] else ""
    try:
        subprocess.run([sys.executable, os.path.join(PKG, "logic_tree", "run_logic_tree.py"), tmp],
                       cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        pz = pd.read_csv(os.path.join(out, f"scenario_rates_{z}{suffix}.csv"))
        pzm = pd.read_csv(os.path.join(out, f"scenario_rates_{z}{suffix}_variable_mu.csv"))
    finally:
        os.remove(tmp)
        shutil.rmtree(out, ignore_errors=True)

    with netCDF4.Dataset(a.hs_nc) as ds:
        ds.set_auto_mask(False)
        g = lambda n: np.asarray(ds[n][:], float)
        u = g("uniform_event_row").astype(int)
        n_published = int(u.max())
        first = {}
        for i, k in enumerate(u):
            first.setdefault(int(k), i)
        n = len(pz)                      # the engine's Mw range stops at its Mmax
        rows = np.array([first[k] for k in range(1, n + 1)])
        pub = {k: g(k)[rows] for k in ("weight_with_nonzero_rate", "variable_mu_weight_with_nonzero_rate", "Mw")}
        fam = {}
        for name in ("rate_annual", "variable_mu_rate_annual"):
            s = np.zeros(n_published + 2)
            np.add.at(s, u, g(name))
            fam[name] = s[1:n + 1]
        out_of_range = {k: g(k)[[first[j] for j in range(n + 1, n_published + 1)]] for k in
                        ("Mw", "weight_with_nonzero_rate")}

    np.savez(a.out.replace(".json", ".npz"), mw=pub["Mw"], w_eng=pz.weight_with_nonzero_rate.to_numpy(),
             w_pub=pub["weight_with_nonzero_rate"], wmu_eng=pzm.variable_mu_weight_with_nonzero_rate.to_numpy(),
             wmu_pub=pub["variable_mu_weight_with_nonzero_rate"], r_eng=pz.rate_mean.to_numpy(), r_pub=fam["rate_annual"],
             rmu_eng=pzm.rate_mean.to_numpy(), rmu_pub=fam["variable_mu_rate_annual"])
    res = {"zone": z, "mode": a.mode, "n_parents": int(n), "n_parents_published": n_published,
           "published_parents_beyond_engine_range": {"Mw": [round(float(m), 1) for m in out_of_range["Mw"]],
                                                    "weight_with_nonzero_rate": [float(w) for w in out_of_range["weight_with_nonzero_rate"]]},
           "segments": list(cfg["segments"])}
    for lab, e, r in (("weight_with_nonzero_rate", pz.weight_with_nonzero_rate.to_numpy(), pub["weight_with_nonzero_rate"]),
                      ("variable_mu_weight_with_nonzero_rate", pzm.variable_mu_weight_with_nonzero_rate.to_numpy(),
                       pub["variable_mu_weight_with_nonzero_rate"])):
        d = np.abs(e - r)
        res[lab] = {"max_abs_diff": float(d.max()), "n_within_1e-9": int((d < 1e-9).sum()), "n": int(n),
                    "zero_mismatch": int(((e == 0) != (r == 0)).sum()), "n_zero_published": int((r == 0).sum())}
    for lab, e, r in (("rate_mean (constant)", pz.rate_mean.to_numpy(), fam["rate_annual"]),
                      ("variable_mu rate_mean", pzm.rate_mean.to_numpy(), fam["variable_mu_rate_annual"])):
        ok = r > 0
        rel = np.abs(e - r)[ok] / r[ok]
        res[lab] = {"max_rel_diff": float(rel.max()), "n_within_1e-9": int((rel < 1e-9).sum()), "n_compared": int(ok.sum()),
                    "n_published_zero": int((~ok).sum()), "engine_nonzero_where_published_zero": int((e[~ok] != 0).sum())}
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
    if a.mode == "plain":
        return 0                          # no pass/fail: on a segmented zone it is expected to differ
    good = all(res[k]["n_within_1e-9"] == res[k].get("n", res[k].get("n_compared")) for k in
               ("weight_with_nonzero_rate", "variable_mu_weight_with_nonzero_rate", "rate_mean (constant)", "variable_mu rate_mean"))
    good &= all(res[k]["zero_mismatch"] == 0 for k in ("weight_with_nonzero_rate", "variable_mu_weight_with_nonzero_rate"))
    good &= all(res[k]["engine_nonzero_where_published_zero"] == 0 for k in ("rate_mean (constant)", "variable_mu rate_mean"))
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
