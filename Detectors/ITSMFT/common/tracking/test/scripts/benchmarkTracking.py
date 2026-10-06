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

"""Benchmark CPU/GPU tracking using fixtures prepared by validateHybridTracking.py.

Validation is disabled. Runs are sequential, use fresh workflow processes, and
alternate backend order. The tracker timer includes host preparation, GPU
allocation/transfers/synchronization and MC labeling; workflow wall time also
includes startup, conditions, cluster input and ROOT output. This is a first-TF
benchmark, not a measurement of steady-state kernel throughput.
"""
import argparse
import json
from pathlib import Path
import re
import statistics
import time

from validateHybridTracking import link, run


def summarize(samples):
    batches = sorted({sample['batch'] for sample in samples})
    result = {}
    for backend in sorted({sample['backend'] for sample in samples}):
        rows = [sample for sample in samples if sample['backend'] == backend]
        per_batch = []
        for batch in batches:
            selected = [sample for sample in rows if sample['batch'] == batch]
            if selected:
                per_batch.append({'batch': batch,
                                  'tracking_ms': statistics.median(x['tracking_ms'] for x in selected),
                                  'wall_s': statistics.median(x['wall_s'] for x in selected),
                                  'tracks': selected[0]['tracks']})
        times = [x['tracking_ms'] for x in per_batch]
        result[backend] = {'batches': per_batch,
                           'sum_batch_medians_ms': sum(times),
                           'median_batch_ms': statistics.median(times),
                           'min_batch_ms': min(times), 'max_batch_ms': max(times),
                           'mean_workflow_wall_s': statistics.mean(x['wall_s'] for x in per_batch)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tracker', required=True)
    parser.add_argument('--backend', choices=('cuda', 'hip'), default='cuda')
    parser.add_argument('--first', type=int, default=0)
    parser.add_argument('--last', type=int, default=19)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--mode', choices=('sync', 'async'), default='sync')
    parser.add_argument('--vertex-source', choices=('diamond', 'truth'), default='diamond')
    args = parser.parse_args()
    if args.first < 0 or args.last < args.first or args.repetitions < 1 or args.threads < 1:
        parser.error('invalid batch range, repetition count or thread count')
    args.fixtures = args.fixtures.resolve()
    args.output = args.output.resolve()
    if args.output == args.fixtures or args.fixtures in args.output.parents:
        parser.error('--output must be outside the fixture directory')
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output / 'report.json').exists():
        parser.error('use a new output directory to preserve previous measurements')
    report = {'configuration': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              'timing_scope': __doc__, 'samples': []}
    pattern = re.compile(r'ITS CA tracking produced (\d+) tracks in ([\d.]+) ms')
    def measure(number, backend, repetition):
        name = f'batch_{number:03d}'
        fixture = args.fixtures / name
        directory = args.output / name / backend
        directory.mkdir(parents=True, exist_ok=True)
        for original in (fixture / 'cpu').glob('*.root'):
            if original.name != 'o2trac_its_ca.root':
                link(original, directory / original.name)
        link(fixture / 'cpu/o2simdigitizerworkflow_configuration.ini', directory / 'o2simdigitizerworkflow_configuration.ini')
        if not (directory / 'o2clus_its.root').is_file() or not (fixture / 'ccdb').is_dir():
            raise RuntimeError(f'Missing prepared fixture {fixture}')
        command = [args.tracker, '--tracklet-backend', backend, '--vertex-source', args.vertex_source,
                   '--tracking-mode', args.mode, '--use-geom', '--condition-backend', f'file://{fixture / "ccdb"}',
                   '--configKeyValues', 'ITSAlpideParam.roFrameLengthInBC=198;ITSAlpideParam.roFrameBiasInBC=64;'
                   f'ITSCommonCATrackerParam.nThreads={args.threads}', '-b', '--run']
        log = f'tracking-{repetition}.log'
        start = time.perf_counter()
        run(command, directory, log, ['o2-its-cluster-reader-workflow', '--with-mc', '-b'])
        wall = time.perf_counter() - start
        text = (directory / log).read_text()
        matches = pattern.findall(text)
        if not matches or 'tracking dropped this TimeFrame' in text:
            raise RuntimeError(f'Missing tracking results in {directory / log}')
        tracks = sum(int(x[0]) for x in matches)
        if not tracks:
            raise RuntimeError('Empty tracking output')
        return {'batch': name, 'backend': backend, 'repetition': repetition,
                'tracks': tracks, 'timeframes': len(matches),
                'frame_tracks': [int(x[0]) for x in matches],
                'frame_timings_ms': [float(x[1]) for x in matches],
                'tracking_ms': sum(float(x[1]) for x in matches), 'wall_s': wall,
                'command': command, 'log': str(directory / log)}

    # Prime executable/filesystem caches. Each measured run still has a fresh
    # CUDA context and processes its first timeframe.
    for backend in ('cpu', args.backend):
        measure(args.first, backend, 'warmup')
    for repetition in range(args.repetitions):
        for number in range(args.first, args.last + 1):
            order = ('cpu', args.backend) if (number + repetition) % 2 == 0 else (args.backend, 'cpu')
            for backend in order:
                sample = measure(number, backend, repetition)
                report['samples'].append(sample)
                report['summary'] = summarize(report['samples'])
                (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                print(f'{sample["batch"]} repeat={repetition} {backend}: '
                      f'{sample["tracking_ms"]:.2f} ms tracking, {sample["wall_s"]:.3f} s workflow, '
                      f'{sample["tracks"]} tracks', flush=True)
    for backend, summary in report['summary'].items():
        print(f'{backend}: {summary["sum_batch_medians_ms"]:.2f} ms total of batch medians', flush=True)


if __name__ == '__main__':
    main()
