"""Physical constants and per-gas unit-conversion factors.

Transcribed verbatim from ``MD.cpp`` (Foley, Sweet, Akinfenwa, GPLv3) — the C
source is the authoritative reference; do not "clean up" these values, they
must match to full precision (see tests/test_units.py).
"""
from __future__ import annotations

from dataclasses import dataclass

# Avogadro's number, mol^-1
NA = 6.022140857e23
# Boltzmann constant, m^2 kg / (s^2 K)
KB_SI = 1.38064852e-23


@dataclass(frozen=True)
class GasConstants:
    """Natural-unit <-> SI conversion factors for one noble gas."""

    name: str
    vol_fac: float    # m^3 per natural unit of volume
    press_fac: float  # Pa per natural unit of pressure
    temp_fac: float   # K per natural unit of temperature
    timefac: float    # s per natural unit of time


GAS_CONSTANTS: dict[str, GasConstants] = {
    "He": GasConstants("He", 1.8399744000000005e-29, 8152287.336171632,
                        10.864459551225972, 1.7572698825166272e-12),
    "Ne": GasConstants("Ne", 2.0570823999999997e-29, 27223022.27659913,
                        40.560648991243625, 2.1192341945685407e-12),
    "Ar": GasConstants("Ar", 3.7949992920124995e-29, 51695201.06691862,
                        142.0950000000000, 2.09618e-12),
    "Kr": GasConstants("Kr", 4.5882712000000004e-29, 59935428.40275003,
                        199.1817584391428, 8.051563913585078e-13),
    "Xe": GasConstants("Xe", 5.4872e-29, 70527773.72794868,
                        280.30305642163006, 9.018957925790732e-13),
}


def get_gas_constants(gas: str) -> GasConstants:
    """Look up per-gas unit factors, matching MD.cpp's fallback-to-Ar behavior
    only for the empty/unset case; unlike the C, an unrecognized gas string is
    an error here rather than a silent fallback."""
    try:
        return GAS_CONSTANTS[gas]
    except KeyError as exc:
        valid = ", ".join(GAS_CONSTANTS)
        raise ValueError(f"Unknown gas {gas!r}; must be one of: {valid}") from exc


def steps_for_gas(gas: str) -> tuple[float, int]:
    """Return (dt_fraction_seconds, NumTime) matching MD.cpp's timestep rule.

    dt is returned as the SI seconds numerator (0.2e-14 or 0.5e-14); divide by
    timefac to get dt in natural units, exactly as MD.cpp does.
    """
    if gas == "He":
        return 0.2e-14, 50000
    return 0.5e-14, 20000
