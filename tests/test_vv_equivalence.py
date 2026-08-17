"""The cached single-force-per-step velocity Verlet must be numerically
identical to the naive reference implementation that recomputes forces twice
per step (matching MD.cpp's VelocityVerlet exactly)."""
import numpy as np

from noblegasmd import engine


def test_cached_matches_naive_trajectory():
    rng = np.random.default_rng(0)
    N, L, dt, n_steps = 20, 6.0, 1e-3, 50

    pos0 = engine.initialize_lattice(N, L)
    vel0 = engine.initialize_velocities(N, 1.0, rng)

    pos_naive = pos0.copy()
    vel_naive = vel0.copy()
    for _ in range(n_steps):
        engine.vv_step_naive(pos_naive, vel_naive, N, L, dt)

    pos_cached = pos0.copy()
    vel_cached = vel0.copy()
    acc_cached = engine.compute_forces(pos_cached, N)
    for _ in range(n_steps):
        engine.vv_step_cached(pos_cached, vel_cached, acc_cached, N, L, dt)

    np.testing.assert_allclose(pos_cached, pos_naive, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(vel_cached, vel_naive, rtol=1e-9, atol=1e-12)
