#!/usr/bin/env python3
# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""Run digit clustering, CPU tracking and hybrid tracking in an O2 environment.

Input batches and calibration snapshots are read-only. All products and logs go
under --output. A mismatch, empty output, or failed workflow exits nonzero.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess


def link(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        destination.unlink()
    elif destination.exists():
        raise RuntimeError(f"Refusing to replace {destination}")
    destination.symlink_to(source.resolve())


def run(command, directory, log, upstream=None):
    with (directory / log).open("w") as output:
        # A DPL workflow whose stdin is not a terminal waits for an upstream
        # workflow description on it: the first one gets an empty stdin.
        if upstream is None:
            subprocess.run(command, cwd=directory, stdin=subprocess.DEVNULL, stdout=output, stderr=output, check=True)
        else:
            producer = subprocess.Popen(upstream, cwd=directory, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=output)
            try:
                result = subprocess.run(command, cwd=directory, stdin=producer.stdout,
                                        stdout=output, stderr=output)
                producer.stdout.close()
                producer_code = producer.wait()
                result.check_returncode()
                if producer_code:
                    raise RuntimeError(f"Cluster reader failed: {producer_code}")
            finally:
                if producer.poll() is None:
                    producer.terminate()
                    producer.wait()


def compare(cpu, hybrid):
    import ROOT
    files = [ROOT.TFile.Open(str(path)) for path in (cpu, hybrid)]
    if any(not file or file.IsZombie() for file in files):
        raise RuntimeError("Cannot open tracking outputs")
    try:
        trees = [file.Get("o2sim") for file in files]
        if any(not tree for tree in trees) or trees[0].GetEntries() <= 0 or trees[0].GetEntries() != trees[1].GetEntries():
            raise RuntimeError("Empty or unequal tracking entry counts")
        tracks = 0
        branches = ("ITSTrack", "ITSTrackClusIdx", "ITSTracksROF", "ITSTrackMCTruth")
        for entry in range(trees[0].GetEntries()):
            for tree in trees:
                if tree.GetEntry(entry) <= 0:
                    raise RuntimeError(f"Cannot read tracking entry {entry}")
            for branch in branches:
                if not all(tree.GetBranch(branch) for tree in trees):
                    raise RuntimeError(f"Missing branch {branch}")
            if not ROOT.equivalentTrackingEvent(*(getattr(tree, branch) for tree in trees for branch in branches)):
                raise RuntimeError(f"CPU/hybrid mismatch: entry {entry}")
            tracks += len(trees[0].ITSTrack)
        if not tracks:
            raise RuntimeError("Comparison produced no tracks")
        return tracks
    finally:
        for file in files:
            file.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="directory containing batch_000, ...")
    parser.add_argument("--ccdb", type=Path, required=True, help="local calibration snapshot directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tracker", default="o2-its-ca-tracker-workflow")
    parser.add_argument("--backend", choices=("cuda", "hip"), default="cuda")
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--last", type=int, default=19)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--mode", choices=("sync", "async"), default="sync")
    parser.add_argument("--vertex-source", choices=("diamond", "truth"), default="diamond")
    parser.add_argument("--reuse-clusters", action="store_true")
    parser.add_argument("--keep-going", action="store_true", help="check all batches and report failures; still exit nonzero")
    args = parser.parse_args()
    if args.first < 0 or args.last < args.first or args.threads < 1:
        parser.error("invalid batch range or thread count")
    args.output = args.output.resolve()
    args.input = args.input.resolve()
    args.ccdb = args.ccdb.resolve()
    if args.output == args.input or args.input in args.output.parents:
        parser.error("--output must be outside the input data directory")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if ROOT.gROOT.LoadMacro(str(Path(__file__).with_name("compareITSTracking.C"))) < 0:
        raise RuntimeError("Cannot load comparison helper")
    report = {"backend": args.backend, "mode": args.mode, "vertices": args.vertex_source,
              "threads": args.threads, "comparison": "exact per-ROF track identities and metadata; float absolute/relative tolerances",
              "tolerances": {"reference": [1e-5, 1e-4], "angle": [1e-6, 1e-4],
                             "parameters": [1e-5, 2e-3], "covariance": [1e-8, 2e-2], "chi2": [1e-2, 5e-3]},
              "batches": [], "failures": []}
    for number in range(args.first, args.last + 1):
        name = f"batch_{number:03d}"
        try:
            source = args.input / name
            if not (source / "itsdigits.root").is_file():
                raise RuntimeError(f"Missing ITS digits in {source}")
            output = args.output / name
            snapshot = output / "ccdb"
            for original in args.ccdb.rglob("snapshot.root"):
                if original.exists():
                    link(original, snapshot / original.relative_to(args.ccdb))
            for condition, filename in (("GRPECS", "o2sim_grpecs.root"),
                                        ("GRPMagField", "o2sim_grpMagField.root"),
                                        ("GeometryAligned", "o2sim_geometry-aligned.root")):
                link(source / filename, snapshot / "GLO/Config" / condition / "snapshot.root")
            for mode in ("cluster", "cpu", args.backend):
                directory = output / mode
                directory.mkdir(parents=True, exist_ok=True)
                for original in list(source.glob("*.root")) + [source / "o2simdigitizerworkflow_configuration.ini"]:
                    link(original, directory / original.name)
            condition = ["--condition-backend", f"file://{snapshot}"]
            timing = "ITSAlpideParam.roFrameLengthInBC=198;ITSAlpideParam.roFrameBiasInBC=64"
            clusters = output / "cluster/o2clus_its.root"
            if not args.reuse_clusters or not clusters.exists():
                run(["o2-its-reco-workflow", "--disable-tracking", *condition,
                     "--configKeyValues", timing, "-b", "--run"], output / "cluster", "cluster.log")
            for backend in ("cpu", args.backend):
                directory = output / backend
                link(clusters, directory / clusters.name)
                command = [args.tracker, "--tracklet-backend", backend, "--vertex-source", args.vertex_source,
                           "--tracking-mode", args.mode, "--use-geom", *condition, "--configKeyValues",
                           timing + f";ITSCommonCATrackerParam.nThreads={args.threads}", "-b", "--run"]
                if backend != "cpu":
                    command.extend(["--validate-tracklets", "--validate-cells", "--validate-neighbours", "--validate-roads", "--validate-refit"])
                run(command, directory, "tracking.log", ["o2-its-cluster-reader-workflow", "--with-mc", "-b"])
            tracks = compare(output / "cpu/o2trac_its_ca.root", output / args.backend / "o2trac_its_ca.root")
            log = (output / args.backend / "tracking.log").read_text()
            candidates = sum(map(int, re.findall(r"CPU/GPU tracklet validation:.*candidates=(\d+) matched", log)))
            cells = sum(map(int, re.findall(r"CPU/GPU cell validation:.*cells=(\d+) matched", log)))
            neighbours = sum(map(int, re.findall(r"GPU neighbour validation: (\d+) neighbours", log)))
            roads = sum(map(int, re.findall(r"GPU road validation: (\d+) extensions", log)))
            refits = sum(map(int, re.findall(r"Validated (\d+) GPU refits against CPU", log)))
            if not refits or not candidates or not cells or not neighbours or not roads or "tracking dropped this TimeFrame" in log:
                raise RuntimeError(f"No validated GPU candidates or a dropped timeframe in {name}")
            report["batches"].append({"batch": name, "tracks": tracks, "validated_candidates": candidates, "validated_cells": cells, "validated_neighbours": neighbours, "validated_road_extensions": roads, "validated_refits": refits})
            (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            print(f"{name}: {tracks} equivalent tracks, {candidates} validated candidates, {cells} validated cells, {neighbours} validated neighbours, {roads} validated road extensions, {refits} validated refits", flush=True)
        except Exception as error:
            report["failures"].append({"batch": name, "error": str(error)})
            (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            print(f"{name}: FAILED: {error}", flush=True)
            if not args.keep_going:
                raise
    if report["failures"]:
        raise RuntimeError(f"{len(report['failures'])} batch comparisons failed; see report.json")



if __name__ == "__main__":
    main()
