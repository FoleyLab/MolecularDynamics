"""Constants round-trip: per-gas factors must match MD.cpp to full precision."""
from noblegasmd.units import GAS_CONSTANTS, KB_SI, NA, get_gas_constants

# Transcribed directly from the literal values in MD.cpp's main().
EXPECTED = {
    "He": dict(vol_fac=1.8399744000000005e-29, press_fac=8152287.336171632,
               temp_fac=10.864459551225972, timefac=1.7572698825166272e-12),
    "Ne": dict(vol_fac=2.0570823999999997e-29, press_fac=27223022.27659913,
               temp_fac=40.560648991243625, timefac=2.1192341945685407e-12),
    "Ar": dict(vol_fac=3.7949992920124995e-29, press_fac=51695201.06691862,
               temp_fac=142.0950000000000, timefac=2.09618e-12),
    "Kr": dict(vol_fac=4.5882712000000004e-29, press_fac=59935428.40275003,
               temp_fac=199.1817584391428, timefac=8.051563913585078e-13),
    "Xe": dict(vol_fac=5.4872e-29, press_fac=70527773.72794868,
               temp_fac=280.30305642163006, timefac=9.018957925790732e-13),
}


def test_global_constants_match_md_cpp():
    assert NA == 6.022140857e23
    assert KB_SI == 1.38064852e-23


def test_all_five_gases_present():
    assert set(GAS_CONSTANTS) == set(EXPECTED)


def test_per_gas_constants_match_md_cpp_exactly():
    for gas, expected in EXPECTED.items():
        gc = get_gas_constants(gas)
        assert gc.vol_fac == expected["vol_fac"], gas
        assert gc.press_fac == expected["press_fac"], gas
        assert gc.temp_fac == expected["temp_fac"], gas
        assert gc.timefac == expected["timefac"], gas


def test_unknown_gas_raises():
    import pytest
    with pytest.raises(ValueError):
        get_gas_constants("Rn")
