"""Golden-oracle harness: drives the compiled MD.cpp reference binary via stdin
and parses its <title>_average.txt output. Used only by tests, never by the package.
"""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
ORACLE_SRC = REPO_ROOT / "MD.cpp"
ORACLE_BIN = Path(__file__).resolve().parent / "MD_oracle.exe"


@dataclass(frozen=True)
class OracleResult:
    T: float   # average temperature, K
    P: float   # average pressure, Pa
    gc: float  # PV/nT, J/(mol K)
    Z: float   # compressibility factor, unitless
    V: float   # volume, m^3
    N: int


def build_oracle(force: bool = False) -> Path:
    """Compile MD.cpp into ORACLE_BIN if not already built (or force rebuild)."""
    if ORACLE_BIN.exists() and not force:
        return ORACLE_BIN
    subprocess.run(
        ["g++", "-O3", "-I", str(REPO_ROOT), "-L", str(REPO_ROOT),
         str(ORACLE_SRC), "-o", str(ORACLE_BIN)],
        check=True,
    )
    return ORACLE_BIN


def run_oracle(gas: str, T: float, rho: float, title: str = "oracle_run") -> OracleResult:
    """Run the C reference MD.cpp for one state point and return its averages.

    MD.cpp prompts on stdin for: title, gas, T (K), rho (mol/m^3) and writes
    averages to '<title>_average.txt' in the current working directory. The
    binary's RNG is unseeded `rand()`, which glibc/libstdc++ seed identically
    to seed 1 on every run, so results are deterministic.
    """
    build_oracle()
    stdin_text = f"{title}\n{gas}\n{T}\n{rho}\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(
            [str(ORACLE_BIN)],
            input=stdin_text,
            text=True,
            cwd=tmpdir,
            check=True,
            capture_output=True,
            timeout=600,
        )
        avg_path = Path(tmpdir) / f"{title}_average.txt"
        return _parse_average_file(avg_path)


def _parse_average_file(path: Path) -> OracleResult:
    lines = path.read_text().splitlines()
    # Line 0: header, line 1: dashes, line 2: data row.
    data_line = lines[2]
    tokens = data_line.split()
    # time, T, P, gc, Z, V, N
    _time, T, P, gc, Z, V, N = tokens
    return OracleResult(T=float(T), P=float(P), gc=float(gc), Z=float(Z), V=float(V), N=int(N))
