"""Oracle-free: at low density, the compressibility factor Z -> 1 (ideal gas).

The instantaneous pressure estimator (momentum transfer at wall collisions)
is intrinsically noisy for a single seed -- collisions are rare events, so a
single trajectory's Z can swing several percent even at low density. We
therefore average over a handful of replicate seeds, exactly as the published
methodology and the oracle-agreement test do. rho=40 mol/m^3 is MD.cpp's own
reference value for "ideal gas at STP" (see the density prompt in MD.cpp).
"""
from noblegasmd import sweep

IDEAL_GAS_Z_TOLERANCE = 0.02
IDEAL_GAS_RHO = 40.0
N_REPLICATES = 5


def test_z_approaches_one_at_low_density():
    df = sweep(gas="Ar", T=[300.0], rho=[IDEAL_GAS_RHO], n_replicates=N_REPLICATES, seed=42)
    z_mean = df["Z"].mean()
    assert abs(z_mean - 1.0) < IDEAL_GAS_Z_TOLERANCE, (
        f"Expected ensemble-mean Z near 1 at low density, got Z={z_mean}"
    )


def test_z_approaches_one_across_low_density_temperatures():
    for T in (100.0, 200.0, 300.0, 400.0):
        df = sweep(gas="Ar", T=[T], rho=[IDEAL_GAS_RHO], n_replicates=N_REPLICATES, seed=int(T))
        z_mean = df["Z"].mean()
        assert abs(z_mean - 1.0) < IDEAL_GAS_Z_TOLERANCE, (
            f"T={T}: expected ensemble-mean Z near 1 at low density, got Z={z_mean}"
        )
