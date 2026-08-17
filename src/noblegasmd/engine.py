"""Numba kernels: Lennard-Jones forces, velocity Verlet, elastic walls.

Physics matches ``MD.cpp`` (Foley, Sweet, Akinfenwa, GPLv3) exactly: natural
units (sigma = epsilon = m = kB = 1), full O(N^2) pairwise LJ with no cutoff,
hard elastic walls on [0, L)^3 that reverse velocity without repositioning.

Two step variants are provided:

- ``vv_step_naive`` mirrors the reference's redundant double force-evaluation
  per step exactly, for equivalence testing only.
- ``vv_step_cached`` is the production path: one force evaluation per step,
  reusing the previous step's post-kick acceleration as this step's initial
  acceleration (numerically identical to the naive version, ~2x faster).
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def compute_forces(pos: np.ndarray, N: int) -> np.ndarray:
    """O(N^2) pairwise Lennard-Jones acceleration (m=1 => a=F), no cutoff."""
    acc = np.zeros((N, 3))
    for i in range(N - 1):
        for j in range(i + 1, N):
            rij0 = pos[i, 0] - pos[j, 0]
            rij1 = pos[i, 1] - pos[j, 1]
            rij2 = pos[i, 2] - pos[j, 2]
            rSqd = rij0 * rij0 + rij1 * rij1 + rij2 * rij2
            f = 24.0 * (2.0 * rSqd**-7 - rSqd**-4)
            acc[i, 0] += rij0 * f
            acc[i, 1] += rij1 * f
            acc[i, 2] += rij2 * f
            acc[j, 0] -= rij0 * f
            acc[j, 1] -= rij1 * f
            acc[j, 2] -= rij2 * f
    return acc


@njit(cache=True, fastmath=True)
def compute_potential(pos: np.ndarray, N: int) -> float:
    """Total LJ potential energy, summed once per distinct pair i<j."""
    pot = 0.0
    for i in range(N - 1):
        for j in range(i + 1, N):
            rij0 = pos[i, 0] - pos[j, 0]
            rij1 = pos[i, 1] - pos[j, 1]
            rij2 = pos[i, 2] - pos[j, 2]
            rSqd = rij0 * rij0 + rij1 * rij1 + rij2 * rij2
            rnorm = np.sqrt(rSqd)
            quot = 1.0 / rnorm
            term1 = quot**12
            term2 = quot**6
            pot += 4.0 * (term1 - term2)
    return pot


@njit(cache=True, fastmath=True)
def kinetic_energy(vel: np.ndarray, N: int) -> float:
    kin = 0.0
    for i in range(N):
        v2 = vel[i, 0] ** 2 + vel[i, 1] ** 2 + vel[i, 2] ** 2
        kin += 0.5 * v2
    return kin


@njit(cache=True, fastmath=True)
def mean_squared_velocity(vel: np.ndarray, N: int) -> float:
    vx2 = 0.0
    vy2 = 0.0
    vz2 = 0.0
    for i in range(N):
        vx2 += vel[i, 0] * vel[i, 0]
        vy2 += vel[i, 1] * vel[i, 1]
        vz2 += vel[i, 2] * vel[i, 2]
    return (vx2 + vy2 + vz2) / N


@njit(cache=True, fastmath=True)
def _apply_walls(pos: np.ndarray, vel: np.ndarray, N: int, L: float, dt: float) -> float:
    """Elastic walls: reverse velocity (not position) on out-of-box coords.

    Returns psum, the summed momentum-transfer contribution to pressure.
    """
    psum = 0.0
    for i in range(N):
        for k in range(3):
            if pos[i, k] < 0.0:
                vel[i, k] *= -1.0
                psum += 2.0 * abs(vel[i, k]) / dt
            if pos[i, k] >= L:
                vel[i, k] *= -1.0
                psum += 2.0 * abs(vel[i, k]) / dt
    return psum


@njit(cache=True, fastmath=True)
def vv_step_naive(pos: np.ndarray, vel: np.ndarray, N: int, L: float, dt: float) -> float:
    """One velocity-Verlet step, recomputing forces twice (matches MD.cpp's
    ``VelocityVerlet`` exactly, including its redundant leading force call).
    Mutates pos/vel in place. Returns instantaneous pressure (natural units).
    """
    acc = compute_forces(pos, N)  # redundant recompute, as in the reference
    for i in range(N):
        for k in range(3):
            pos[i, k] += vel[i, k] * dt + 0.5 * acc[i, k] * dt * dt
            vel[i, k] += 0.5 * acc[i, k] * dt
    acc = compute_forces(pos, N)
    for i in range(N):
        for k in range(3):
            vel[i, k] += 0.5 * acc[i, k] * dt
    psum = _apply_walls(pos, vel, N, L, dt)
    return psum / (6.0 * L * L)


@njit(cache=True, fastmath=True)
def vv_step_cached(pos: np.ndarray, vel: np.ndarray, acc: np.ndarray, N: int,
                    L: float, dt: float) -> float:
    """One velocity-Verlet step given the acceleration already evaluated at
    the current position (from the previous step's second force call, or the
    pre-loop initial force call). Mutates pos/vel/acc in place; acc is left
    holding the force at the *new* position for reuse by the next step.
    Returns instantaneous pressure (natural units).
    """
    for i in range(N):
        for k in range(3):
            pos[i, k] += vel[i, k] * dt + 0.5 * acc[i, k] * dt * dt
            vel[i, k] += 0.5 * acc[i, k] * dt
    new_acc = compute_forces(pos, N)
    for i in range(N):
        for k in range(3):
            vel[i, k] += 0.5 * new_acc[i, k] * dt
            acc[i, k] = new_acc[i, k]
    psum = _apply_walls(pos, vel, N, L, dt)
    return psum / (6.0 * L * L)


@njit(cache=True, fastmath=True)
def simulate(pos0: np.ndarray, vel0: np.ndarray, L: float, dt: float,
             n_steps: int):
    """Run n_steps+1 velocity-Verlet iterations (matching MD.cpp's ``i<NumTime+1``
    loop bound), recording per-iteration instantaneous T, P, KE, PE (all
    natural units). Returns (T_nu, P_nu, KE, PE, pos_final, vel_final).
    """
    N = pos0.shape[0]
    pos = pos0.copy()
    vel = vel0.copy()
    acc = compute_forces(pos, N)

    n_records = n_steps + 1
    T_nu = np.empty(n_records)
    P_nu = np.empty(n_records)
    KE = np.empty(n_records)
    PE = np.empty(n_records)

    for i in range(n_records):
        p_inst = vv_step_cached(pos, vel, acc, N, L, dt)
        mvs = mean_squared_velocity(vel, N)
        T_nu[i] = mvs / 3.0
        P_nu[i] = p_inst
        KE[i] = kinetic_energy(vel, N)
        PE[i] = compute_potential(pos, N)

    return T_nu, P_nu, KE, PE, pos, vel


@njit(cache=True, fastmath=True)
def simulate_with_trajectory(pos0: np.ndarray, vel0: np.ndarray, L: float, dt: float,
                              n_steps: int):
    """Same as ``simulate`` but additionally records positions at every
    iteration (memory-heavy; only used when the caller requests a trajectory).
    """
    N = pos0.shape[0]
    pos = pos0.copy()
    vel = vel0.copy()
    acc = compute_forces(pos, N)

    n_records = n_steps + 1
    T_nu = np.empty(n_records)
    P_nu = np.empty(n_records)
    KE = np.empty(n_records)
    PE = np.empty(n_records)
    traj = np.empty((n_records, N, 3))

    for i in range(n_records):
        p_inst = vv_step_cached(pos, vel, acc, N, L, dt)
        mvs = mean_squared_velocity(vel, N)
        T_nu[i] = mvs / 3.0
        P_nu[i] = p_inst
        KE[i] = kinetic_energy(vel, N)
        PE[i] = compute_potential(pos, N)
        traj[i] = pos

    return T_nu, P_nu, KE, PE, pos, vel, traj


def initialize_lattice(N: int, L: float) -> np.ndarray:
    """Simple-cubic lattice positions, matching MD.cpp's ``initialize()``."""
    n = int(np.ceil(N ** (1.0 / 3.0)))
    pos = np.empty((N, 3))
    p = 0
    spacing = L / n
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if p < N:
                    pos[p, 0] = (i + 0.5) * spacing
                    pos[p, 1] = (j + 0.5) * spacing
                    pos[p, 2] = (k + 0.5) * spacing
                p += 1
    return pos


def initialize_velocities(N: int, T_init_nu: float, rng: np.random.Generator) -> np.ndarray:
    """Gaussian velocities, COM removed, rescaled to T_init_nu (natural units),
    matching MD.cpp's ``initializeVelocities()`` (including its (N-1) factor).
    """
    vel = rng.normal(0.0, 1.0, size=(N, 3))
    vcm = vel.mean(axis=0)
    vel -= vcm
    vsqd_sum = np.sum(vel ** 2)
    lam = np.sqrt(3.0 * (N - 1) * T_init_nu / vsqd_sum)
    vel *= lam
    return vel


def warmup() -> None:
    """Trigger JIT compilation on a throwaway tiny system so it doesn't land
    on the user's first timed simulation."""
    pos = initialize_lattice(4, 4.0)
    vel = np.zeros((4, 3))
    simulate(pos, vel, 4.0, 0.001, 1)
    # also warm the naive path used by the equivalence test
    p2 = pos.copy()
    v2 = np.zeros((4, 3))
    vv_step_naive(p2, v2, 4, 4.0, 0.001)
