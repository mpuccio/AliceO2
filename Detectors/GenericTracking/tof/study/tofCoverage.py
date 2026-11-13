#!/usr/bin/env python3
# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# SPDX-License-Identifier: GPL-3.0-only
"""Truth-hit availability among reconstructed tracks, before testing a matcher.

One QA timeframe and its source's TF3 transport hits. This measures availability,
not matching efficiency or digitization efficiency. Source/event/track identity
is mandatory: a hit in another event does not make a track eligible.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess

import numpy as np
import uproot


def interval(numerator, denominator):
    if denominator == 0:
        return {"fraction": None, "low": None, "high": None}
    value = numerator / denominator
    # Wilson interval with z=1 (approximately 68% coverage).
    centre = (value + 0.5 / denominator) / (1.0 + 1.0 / denominator)
    half = np.sqrt(value * (1.0 - value) / denominator + 0.25 / denominator**2) / (1.0 + 1.0 / denominator)
    return {"fraction": value, "low": float(centre - half), "high": float(centre + half)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa", type=Path, required=True)
    parser.add_argument("--tof-hits", type=Path, required=True)
    parser.add_argument("--first-outer-chip", type=int, required=True, help="from this simulation's TOF GeometryTGeo")
    parser.add_argument("--total-chips", type=int, required=True)
    parser.add_argument("--source", type=int, default=0, help="MC source corresponding to --tof-hits")
    parser.add_argument("--min-pt", type=float, default=0.1)
    parser.add_argument("--max-abs-eta", type=float, default=1.0)
    parser.add_argument("--include-secondaries", action="store_true")
    parser.add_argument("--output", type=Path, required=True, help="output filename prefix")
    args = parser.parse_args()
    if not 0 < args.first_outer_chip < args.total_chips or not 0 <= args.source < 256:
        parser.error("invalid chip partition or source")
    with uproot.open(args.qa) as qa:
        if qa["run"].num_entries != 1:
            parser.error("provide exactly one QA timeframe: MC identities repeat between timeframes")
        tracks = qa["tracks"].arrays(library="np")
    with uproot.open(args.tof_hits) as tof:
        tree = tof["o2sim"]
        track_ids = tree["TF3Hit.mTrackID"].array(library="np")
        chips = tree["TF3Hit.mDetectorID"].array(library="np")
    inner, outer = set(), set()
    hit_rows = 0
    for event, (ids, sensors) in enumerate(zip(track_ids, chips)):
        if len(ids) != len(sensors):
            raise ValueError("inconsistent hit columns")
        if np.any(sensors < 0) or np.any(sensors >= args.total_chips):
            raise ValueError("hit sensor outside the declared barrel geometry")
        for track_id, chip in zip(ids, sensors):
            if track_id < 0:
                continue
            identity = int(track_id) | (event << 31) | (args.source << 50)
            (inner if chip < args.first_outer_chip else outer).add(identity)
            hit_rows += 1
    ids = tracks["mcParticle"]
    selected = (~tracks["fake"].astype(bool) & ~tracks["clone"].astype(bool) & tracks["fitValid"].astype(bool) &
                (tracks["ptGen"] >= args.min_pt) & (np.abs(tracks["etaGen"]) < args.max_abs_eta) &
                (((ids >> np.uint64(50)) & np.uint64(255)) == args.source))
    if not args.include_secondaries:
        selected &= tracks["origin"] == 0  # TrackingQA::ParticleOrigin::kPhysicalPrimary
    identities = ids[selected]
    if len(np.unique(identities)) != len(identities):
        raise ValueError("duplicate identities remain after excluding clones")
    selected_pt = tracks["ptGen"][selected]
    selected_eta = tracks["etaGen"][selected]
    has_inner = np.array([int(identity) in inner for identity in identities], dtype=bool)
    has_outer = np.array([int(identity) in outer for identity in identities], dtype=bool)
    counts = {"selected_tracks": len(identities), "inner": int(has_inner.sum()), "outer": int(has_outer.sum()),
              "both": int((has_inner & has_outer).sum()), "inner_only": int((has_inner & ~has_outer).sum()),
              "outer_only": int((~has_inner & has_outer).sum()), "neither": int((~has_inner & ~has_outer).sum())}
    pt_edges = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0, np.inf])
    bins = []
    for low, high in zip(pt_edges[:-1], pt_edges[1:]):
        mask = (selected_pt >= low) & (selected_pt < high)
        den = int(mask.sum())
        row = {"pt_low": float(low), "pt_high": None if np.isinf(high) else float(high), "denominator": den}
        for name, available in [("inner", has_inner), ("outer", has_outer), ("both", has_inner & has_outer)]:
            num = int((mask & available).sum())
            row[name] = {"numerator": num, **interval(num, den)}
        bins.append(row)
    report = {"kind": "TOF transport-hit availability, NOT matching efficiency",
              "inputs": {"qa": str(args.qa.resolve()), "tof_hits": str(args.tof_hits.resolve())},
              "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
              "selection": {"min_generated_pt_GeV": args.min_pt, "max_abs_generated_eta": args.max_abs_eta,
                            "physical_primaries_only": not args.include_secondaries,
                            "tracking": "good, non-clone, valid refit", "source": args.source},
              "geometry": {"first_outer_chip": args.first_outer_chip, "total_barrel_chips": args.total_chips},
              "events": len(track_ids), "hit_rows": hit_rows, "input_reconstructed_tracks": len(ids),
              "counts": counts, "pt_bins": bins}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    Path(str(args.output) + ".json").write_text(json.dumps(report, indent=2) + "\n")
    with Path(str(args.output) + ".csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mc_particle", "pt_gen_GeV", "eta_gen", "inner_hit", "outer_hit"])
        writer.writerows(zip(identities, selected_pt, selected_eta, has_inner, has_outer))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name, label in [("inner", "Inner TOF"), ("outer", "Outer TOF"), ("both", "Both")]:
        finite = [b for b in bins if b["pt_high"] is not None and b["denominator"]]
        x = [np.sqrt(b["pt_low"] * b["pt_high"]) for b in finite]
        y = [b[name]["fraction"] for b in finite]
        error = [[b[name]["fraction"] - b[name]["low"] for b in finite], [b[name]["high"] - b[name]["fraction"] for b in finite]]
        ax.errorbar(x, y, yerr=error, marker="o", capsize=2, label=label)
    ax.set(xscale="log", ylim=(0, 1.05), xlabel="Generated pT [GeV/c]", ylabel="Fraction with a TOF transport hit",
           title="TOF hit availability among reconstructed tracks")
    ax.text(0.02, 0.03, f"Good unique refits; |eta(gen)| < {args.max_abs_eta}; transport hits\nNot a matching-efficiency measurement", transform=ax.transAxes, fontsize=9)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(args.output) + ".png", dpi=160)
    fig.savefig(str(args.output) + ".pdf")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
