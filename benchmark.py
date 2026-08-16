#!/usr/bin/env python3
"""Timing harness for noblegasmd's numba kernel: reports wall-clock for a few
(N, n_steps) points, so performance on e.g. Colab can be measured directly
rather than guessed. Run with: python benchmark.py
"""
from __future__ import annotations

import time

from noblegasmd import run
from noblegasmd.api import _ensure_warm

POINTS = [
    (108, 2000),
    (216, 20000),   # default production point
    (216, 50000),   # He-equivalent step count
    (512, 5000),
]


def main() -> None:
    print("Warming up JIT (not timed)...")
    t0 = time.perf_counter()
    _ensure_warm()
    print(f"  warmup: {time.perf_counter() - t0:.2f}s\n")

    print(f"{'N':>6}  {'n_steps':>8}  {'wall (s)':>10}  {'steps/s':>10}")
    for n_particles, n_steps in POINTS:
        t0 = time.perf_counter()
        run(gas="Ar", T=300.0, rho=40.0, n_particles=n_particles,
            n_steps=n_steps, seed=0)
        elapsed = time.perf_counter() - t0
        print(f"{n_particles:>6}  {n_steps:>8}  {elapsed:>10.3f}  {n_steps / elapsed:>10.1f}")


if __name__ == "__main__":
    main()
