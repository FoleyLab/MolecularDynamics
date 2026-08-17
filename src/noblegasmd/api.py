"""Public API: run() for a single simulation, sweep() for a parameter grid."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from . import engine
from .units import KB_SI, NA, get_gas_constants, steps_for_gas

_warmed_up = False


def _ensure_warm() -> None:
    global _warmed_up
    if not _warmed_up:
        engine.warmup()
        _warmed_up = True


@dataclass
class RunResult:
    """Result of a single :func:`run` call.

    Averages (``T_avg``, ``P_avg``, ``Z``, ``gc``) follow MD.cpp's convention:
    with ``n_equil=0`` they are means over all recorded iterations, matching
    the published methodology. ``n_equil > 0`` discards that many leading
    iterations before averaging, which is a documented departure from the
    paper, useful only for pedagogical experiments.
    """

    gas: str
    T_set: float
    rho_set: float
    N: int
    n_steps: int
    n_equil: int
    seed: Optional[int]
    T_avg: float          # K
    P_avg: float          # Pa
    Z: float               # compressibility factor, PV/(N kB T)
    gc: float               # PV/(nT), J/(mol K)
    V: float                 # m^3
    energy_drift: float       # max |E(t)-E(0)| / |E(0)|, natural units
    instantaneous_T: np.ndarray = field(repr=False)  # K, per iteration
    instantaneous_P: np.ndarray = field(repr=False)  # Pa, per iteration
    kinetic_energy: np.ndarray = field(repr=False)   # natural units, per iteration
    potential_energy: np.ndarray = field(repr=False)  # natural units, per iteration
    trajectory: Optional[np.ndarray] = field(default=None, repr=False)  # (n_records, N, 3)


def run(
    gas: str = "Ar",
    T: float = 300.0,
    rho: float = 40.0,
    n_particles: int = 216,
    n_steps: Optional[int] = None,
    n_equil: int = 0,
    seed: Optional[int] = None,
    record_trajectory: bool = False,
) -> RunResult:
    """Run a single NVE Lennard-Jones simulation and return thermodynamic averages.

    Parameters
    ----------
    gas : one of "He", "Ne", "Ar", "Kr", "Xe".
    T : initial temperature, Kelvin. Must be positive.
    rho : number density, mol/m^3. Must be positive.
    n_particles : number of particles (simple-cubic lattice init).
    n_steps : number of Verlet steps; defaults to MD.cpp's convention
        (50000 for He, 20000 otherwise).
    n_equil : leading iterations discarded before averaging. Default 0
        matches the published methodology (no discard); values > 0 are a
        documented departure, for pedagogical experiments only.
    seed : seed for the velocity-initialization RNG (numpy PCG64). The
        reference C uses unseeded (fixed-seed) ``rand()``; this port uses an
        independent RNG by design, so trajectories will not match the C
        step-for-step — validate against statistical averages only.
    record_trajectory : if True, also return the full position trajectory
        (memory: ~24 bytes * N * (n_steps+1)).
    """
    if T <= 0:
        raise ValueError(f"T must be positive, got {T}")
    if rho <= 0:
        raise ValueError(f"rho must be positive, got {rho}")

    gc_const = get_gas_constants(gas)

    Vol = n_particles / (rho * NA)
    Vol /= gc_const.vol_fac
    if Vol < n_particles:
        raise ValueError(
            f"Density too high: N={n_particles} particles but only {Vol:.4f} "
            "natural-unit volume available. Simulations with density greater "
            "than 1 particle/(natural unit of volume) may diverge."
        )
    L = Vol ** (1.0 / 3.0)

    dt_numerator, default_n_steps = steps_for_gas(gas)
    dt = dt_numerator / gc_const.timefac
    if n_steps is None:
        n_steps = default_n_steps
    if n_equil < 0 or n_equil > n_steps:
        raise ValueError(f"n_equil must be in [0, n_steps], got {n_equil}")

    T_init_nu = T / gc_const.temp_fac

    rng = np.random.Generator(np.random.PCG64(seed))
    pos0 = engine.initialize_lattice(n_particles, L)
    vel0 = engine.initialize_velocities(n_particles, T_init_nu, rng)

    _ensure_warm()
    trajectory = None
    if record_trajectory:
        T_nu, P_nu, KE, PE, _pos_f, _vel_f, trajectory = engine.simulate_with_trajectory(
            pos0, vel0, L, dt, n_steps
        )
    else:
        T_nu, P_nu, KE, PE, _pos_f, _vel_f = engine.simulate(pos0, vel0, L, dt, n_steps)

    T_si = T_nu * gc_const.temp_fac
    P_si = P_nu * gc_const.press_fac

    # Match MD.cpp exactly at n_equil=0: sum over all NumTime+1 recorded
    # iterations, divide by NumTime (== n_steps), not by the record count.
    if n_equil == 0:
        T_avg = float(np.sum(T_si) / n_steps)
        P_avg = float(np.sum(P_si) / n_steps)
    else:
        T_avg = float(np.mean(T_si[n_equil:]))
        P_avg = float(np.mean(P_si[n_equil:]))

    V_si = Vol * gc_const.vol_fac
    Z = P_avg * V_si / (n_particles * KB_SI * T_avg)
    gc_val = NA * P_avg * V_si / (n_particles * T_avg)

    E_nu = KE + PE
    energy_drift = float(np.max(np.abs(E_nu - E_nu[0])) / np.abs(E_nu[0]))

    return RunResult(
        gas=gas, T_set=T, rho_set=rho, N=n_particles, n_steps=n_steps,
        n_equil=n_equil, seed=seed,
        T_avg=T_avg, P_avg=P_avg, Z=Z, gc=gc_val, V=V_si,
        energy_drift=energy_drift,
        instantaneous_T=T_si, instantaneous_P=P_si,
        kinetic_energy=KE, potential_energy=PE,
        trajectory=trajectory,
    )


def sweep(
    gas: str = "Ar",
    T: Iterable[float] = (100.0, 200.0, 300.0, 400.0),
    rho: Iterable[float] = (40.0,),
    n_particles: int = 216,
    n_steps: Optional[int] = None,
    n_equil: int = 0,
    n_replicates: int = 1,
    seed: Optional[int] = 0,
) -> pd.DataFrame:
    """Run :func:`run` over a grid of (T, rho) state points x replicate seeds.

    Returns a tidy DataFrame with one row per (state point, replicate):
    gas, T_set, rho_set, seed, T_avg, P_avg, Z, gc, V, N, n_steps.
    """
    rows = []
    base_rng = np.random.default_rng(seed)
    for T_val in T:
        for rho_val in rho:
            for rep in range(n_replicates):
                run_seed = int(base_rng.integers(0, 2**32 - 1))
                result = run(
                    gas=gas, T=T_val, rho=rho_val, n_particles=n_particles,
                    n_steps=n_steps, n_equil=n_equil, seed=run_seed,
                )
                rows.append({
                    "gas": result.gas, "T_set": T_val, "rho_set": rho_val,
                    "seed": run_seed, "T_avg": result.T_avg, "P_avg": result.P_avg,
                    "Z": result.Z, "gc": result.gc, "V": result.V,
                    "N": result.N, "n_steps": result.n_steps,
                    "energy_drift": result.energy_drift,
                })
    return pd.DataFrame(rows)
