"""Oracle-free: elastic walls + conservative forces => total energy is conserved.

Do NOT test momentum conservation here: walls inject impulse by design, so
linear momentum is not conserved and is not a physics bug.
"""
from noblegasmd import run

ENERGY_DRIFT_TOLERANCE = 1e-3


def test_energy_conserved_over_full_run():
    result = run(gas="Ar", T=300.0, rho=40.0, seed=1)
    assert result.energy_drift < ENERGY_DRIFT_TOLERANCE, (
        f"Relative energy drift {result.energy_drift} exceeds "
        f"{ENERGY_DRIFT_TOLERANCE}"
    )


def test_energy_conserved_at_higher_density():
    # Denser system => more frequent close encounters, a harder case for
    # energy conservation.
    result = run(gas="Ar", T=300.0, rho=5000.0, seed=2)
    assert result.energy_drift < ENERGY_DRIFT_TOLERANCE, (
        f"Relative energy drift {result.energy_drift} exceeds "
        f"{ENERGY_DRIFT_TOLERANCE}"
    )
