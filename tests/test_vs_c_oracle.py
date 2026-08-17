"""Python ensemble-mean Z must agree with the deterministic C oracle (MD.cpp)
within combined statistical tolerance, across an Ar grid spanning ideal to
strongly non-ideal densities. MD is chaotic and the Python port uses a
different RNG, so we validate ensemble-mean thermodynamic averages only,
never trajectory identity.
"""
import numpy as np
import pytest
from oracle import run_oracle

from noblegasmd import sweep

Z_AGREEMENT_TOLERANCE = 0.03  # 3% relative, per domain-expert sign-off

TEMPERATURES = [100.0, 200.0, 300.0, 400.0]
DENSITIES = [500.0, 2000.0]  # dense enough that the oracle's single deterministic sample is low-noise;
# rho=40 (sparse wall collisions) is covered instead by the oracle-free, ensemble-averaged
# check in test_ideal_gas_limit.py, since at that density even the C oracle's own single
# run swings several percent from its own long-run mean.
N_REPLICATES = 3


@pytest.mark.parametrize("T_set", TEMPERATURES)
@pytest.mark.parametrize("rho_set", DENSITIES)
def test_ensemble_mean_z_agrees_with_c_oracle(T_set, rho_set):
    df = sweep(
        gas="Ar", T=[T_set], rho=[rho_set],
        n_replicates=N_REPLICATES, seed=int(T_set * 100 + rho_set),
    )
    python_z_mean = df["Z"].mean()

    oracle_result = run_oracle("Ar", T_set, rho_set, title=f"oracle_{T_set:.0f}_{rho_set:.0f}")

    rel_err = abs(python_z_mean - oracle_result.Z) / abs(oracle_result.Z)
    assert rel_err < Z_AGREEMENT_TOLERANCE, (
        f"T={T_set} rho={rho_set}: python Z_mean={python_z_mean:.5f} "
        f"(std={df['Z'].std():.5f}) vs oracle Z={oracle_result.Z:.5f}, "
        f"rel_err={rel_err:.4f}"
    )
