# Handoff Spec: Pure-Python (numba) Port of the NVE Lennard-Jones MD Engine

**Purpose.** Reimplement the real-gas molecular-dynamics engine currently written in C
(`MD.cpp`, Foley/Sweet/Akinfenwa, GPLv3) as a small, `pip install`-able Python package.
The package must reproduce the published quasi-isotherms (Figure 4 of the *J. Chem. Educ.*
article) and give students a clean programmatic API for parameter sweeps and numerical
experiments — runnable in Google Colab with zero build step.

**Source of truth.** `MD.cpp` is authoritative for all physics and all numeric constants.
Where this spec and the C disagree, **the C wins** — transcribe constants directly from the
source rather than from this document, and cross-check against the reference values quoted here.

**The golden rule of validation.** MD is chaotic and the port uses a different RNG, so
trajectories will *not* match step-for-step and are not expected to. Validate on
**statistical/thermodynamic averages** (Z, P, T) and on invariants (energy conservation),
never on trajectory identity.

---

## 1. Physics specification (must match `MD.cpp` exactly)

### Ensemble & integrator
- Pure **NVE** (microcanonical). No thermostat or barostat anywhere in the main loop.
- **Velocity Verlet**, timestep `dt`:
  1. `r += v·dt + 0.5·a·dt²`
  2. `v += 0.5·a·dt`  (first half-kick)
  3. recompute accelerations at new `r`
  4. `v += 0.5·a·dt`  (second half-kick)
  5. apply walls (below)
- **Optimization note:** the reference recomputes forces *twice* per step (once redundantly at
  the top of `VelocityVerlet`, once mid-step). A standard velocity Verlet that caches the
  previous step's second force computation is numerically identical and ~2× faster. Implement
  the cached single-force-per-step version, but assert equivalence against the reference during
  testing.

### Forces (no cutoff)
- Full **O(N²)** pairwise Lennard-Jones over all distinct pairs `i<j`. **No cutoff, no shift,
  no tail correction, no neighbor list, no periodic minimum image.**
- In natural units (σ = ε = m = k_B = 1), with `rSqd = |r_i − r_j|²`, the scalar factor is
  `f = 24·(2·rSqd⁻⁷ − rSqd⁻⁴)` and the force vector contribution is `f · (r_i − r_j)`,
  applied with Newton's third law: `a_i += f·rij`, `a_j −= f·rij`.
- Potential energy: `U = Σ_{i≠j, i<j} 4·(r⁻¹² − r⁻⁶)` (report the sum over distinct pairs; note
  the reference `Potential()` double-counts by looping all `i≠j`, so if you reuse its convention
  divide appropriately — verify which value the reference actually reports and match it).

### Boundary conditions (hard elastic walls)
- Cubic box, side `L` (natural units), spanning `[0, L]` per axis.
- After each full VV step, for every particle and axis `k`:
  - if `r[k] < 0`  → `v[k] *= −1`, accumulate `psum += 2·m·|v[k]| / dt`
  - if `r[k] ≥ L`  → `v[k] *= −1`, accumulate `psum += 2·m·|v[k]| / dt`
- **Do not reflect the position back into the box.** Only the velocity is reversed; the position
  may sit marginally outside `[0, L]` for a step. Match this exactly.
- This wall-collision momentum transfer *is* the pressure (kinetic-theory definition of
  infinitesimal pressure), which is why linear momentum is **not** conserved (walls inject
  impulse) but total energy **is** (collisions are elastic).

### Per-step observables
- Mean-square velocity `mvs = (Σ_i |v_i|²)/N`.
- Instantaneous temperature (K): `T = m·mvs/(3·k_B) · TempFac` — note the **3** d.o.f. per
  particle here (not `3(N−1)`; that factor appears only in the init rescale, below).
- Instantaneous pressure (Pa): `P = (psum/(6·L²)) · PressFac`.
- Kinetic and potential energy in natural units.

### Averaging & reported quantities  ⚠️ subtle
- Accumulate `Tavg` and `Pavg` over **all** `NumTime` steps starting at step 0 — **no
  equilibration discard** — then divide by `NumTime`.
- Compressibility is computed from the **averaged** P and T, not from averaged instantaneous Z:
  `Z = Pavg · (Vol·VolFac) / (N · kBSI · Tavg)`
  `gc = NA · Pavg · (Vol·VolFac) / (N · Tavg)`   (this is PV/nT; → R ≈ 8.314 in the ideal limit)
- To reproduce Figure 4, the default must average over all steps. Expose an optional
  `n_equil` (default **0**) for pedagogical experiments, and document that any value > 0 departs
  from the published methodology.

### Initialization
- **N = 216** by default. Simple-cubic lattice: `n = ceil(N^(1/3))` per axis (= 6 for 216),
  spacing `pos = L/n`, positions at `((i+0.5)·pos, (j+0.5)·pos, (k+0.5)·pos)`.
- Velocities: draw each component from a Gaussian, remove center-of-mass velocity, then rescale
  to the target temperature with `lambda = sqrt(3·(N−1)·Tinit / vSqdSum)` (the `(N−1)` accounts
  for the removed COM d.o.f.). `Tinit` is in natural units (input Kelvin ÷ `TempFac`).
- RNG: use `numpy.random.Generator(PCG64(seed))`. Bit-matching C's `rand()` is neither possible
  nor required.

### Timestep & step count (per gas)
- Non-helium: `dt = 0.5e-14 / timefac`, `NumTime = 20000`.
- Helium: `dt = 0.2e-14 / timefac`, `NumTime = 50000`.

---

## 2. Units and per-gas constants

Natural units per noble gas; SI is used only for reporting. **Transcribe these from `MD.cpp`**;
values below are for cross-checking only.

Global: `NA = 6.022140857e23`, `kBSI = 1.38064852e-23`.

| Gas | VolFac (m³)              | PressFac (Pa)         | TempFac (K)          | timefac (s)             |
|-----|-------------------------|-----------------------|----------------------|-------------------------|
| He  | 1.8399744000000005e-29  | 8152287.336171632     | 10.864459551225972   | 1.7572698825166272e-12  |
| Ne  | 2.0570823999999997e-29  | 27223022.27659913     | 40.560648991243625   | 2.1192341945685407e-12  |
| Ar  | 3.7949992920124995e-29  | 51695201.06691862     | 142.0950000000000    | 2.09618e-12             |
| Kr  | 4.5882712000000004e-29  | 59935428.40275003     | 199.1817584391428    | 8.051563913585078e-13   |
| Xe  | 5.4872e-29              | 70527773.72794868     | 280.30305642163006   | 9.018957925790732e-13   |

Volume from density: input `rho` in **mol/m³**, then `Vol = N/(rho·NA)`, convert to natural
units `Vol /= VolFac`, box side `L = Vol^(1/3)`. Preserve the reference's guard rails
(reject `rho ≤ 0`, warn/της on `Vol < N`, i.e. density above ~1 particle per natural-unit volume).

---

## 3. Target API

Ergonomic for both students and an inline LLM assistant. Rich docstrings and type hints
directly improve autocomplete quality (and reduce wasted AI calls), so treat them as part of
the deliverable.

```python
from noblegasmd import run, sweep   # package name TBD; check PyPI availability

# Single simulation
result = run(
    gas: str = "Ar",            # He | Ne | Ar | Kr | Xe
    T: float = 300.0,           # initial temperature, Kelvin
    rho: float = 40.0,          # number density, mol/m^3
    n_particles: int = 216,
    n_steps: int | None = None, # default 20000 (50000 for He)
    n_equil: int = 0,           # steps discarded before averaging; 0 = match reference
    seed: int | None = None,
    record_trajectory: bool = False,
)
# result exposes at least: T_avg (K), P_avg (Pa), Z, gc (PV/nT), V (m^3), N,
# energy drift diagnostic, and optionally per-step arrays / trajectory.

# Parameter sweep -> tidy DataFrame, one row per (state point, replicate)
df = sweep(
    gas="Ar",
    T=[100, 200, 300, 400],
    rho=np.linspace(1.0, 5000.0, 25),
    n_replicates=3,             # independent seeds per state point for error bars
    seed=0,
)
# Columns: gas, T_set, rho_set, seed, T_avg, P_avg, Z, gc, V, N, n_steps
# Figure 4 = df.groupby("T_set") plotted as Z (or gc) vs rho_set (or P_avg).
```

Keep the heavy loop inside numba; `run`/`sweep` are thin Python wrappers. **Never** place an
LLM/API call inside a sweep loop — that reintroduces the per-run token cost the whole design
avoids.

---

## 4. Implementation guidance (numba)

- Put the integrator + force kernel in `@njit(fastmath=True, cache=True)` functions operating on
  preallocated `(N,3)` float64 arrays. Write the force loop as **explicit nested loops with
  Newton's third law** — do *not* build N×N numpy temporaries (slower and memory-heavy).
- Consider `parallel=True` + `prange` on the outer force loop for a free ~2× on Colab's second
  vCPU; benchmark with and without.
- **Warm up the JIT** with a throwaway tiny run at import or first call so compile time doesn't
  land on the user's first timed simulation.
- No per-step Python callbacks or printing inside the kernel; collect diagnostics into arrays and
  process after the loop.
- Ship a small timing harness (`benchmark.py`) that reports wall-clock for a couple of
  `(N, n_steps)` points, so the team can measure Colab performance directly rather than guess.

---

## 5. Validation protocol

Build the C as a **golden oracle** and compare averages. `MD.cpp` reads its inputs from stdin
(title, gas, T in K, rho in mol/m³) and writes averages to `<title>_average.txt`; drive it by
piping inputs, parse `Z`, `P`, `T` from that file.

Required tests (pytest):
1. **Ideal-gas limit** (oracle-free): at low density, `Z → 1`. Assert `|Z − 1| < ~0.02` for the
   lowest-density state points. Tune threshold with the domain expert.
2. **Energy conservation** (oracle-free): elastic walls + conservative forces ⇒ total energy is
   conserved. Assert relative drift `|E(t) − E(0)| / |E(0)|` stays below a small bound over a full
   run. (Do **not** test momentum conservation — walls break it by design.)
3. **Agreement with the C oracle**: for a grid of Ar state points (e.g. T ∈ {100, 200, 300, 400} K
   across a density range spanning ideal → strongly non-ideal), the Python **ensemble mean** of Z
   over `n_replicates` seeds must agree with the C value within combined statistical error
   (start with ~2–3 %, or within 2σ of the ensemble — final tolerance is the domain expert's call).
   The C is deterministic (unseeded `rand()` ⇒ fixed seed 1), so it yields one value per input;
   compare it against the Python ensemble.
4. **Constants round-trip**: unit-test that the per-gas factors and derived quantities match
   `MD.cpp` to full precision.

Acceptance test / **definition of done**: a `notebooks/reproduce_figure4.ipynb` that calls
`sweep()` and regenerates the published quasi-isotherms within agreed tolerance, plus all tests
green and a clean `pip install` in a fresh Colab runtime.

---

## 6. Packaging & publishing

- Pure Python + numba ⇒ a **pure-Python wheel** (`py3-none-any`), no compiled extensions, no
  `cibuildwheel`/manylinux machinery. Dependencies: `numpy`, `numba`, `pandas` (and `matplotlib`
  for the example notebook only).
- `pyproject.toml` with a standard backend (hatchling or setuptools). Target the Python versions
  Colab currently ships and that numba supports; pin a floor.
- **License must be GPLv3** (derivative of GPLv3 `MD.cpp`). Include `LICENSE` and proper
  attribution to the original authors and the *J. Chem. Educ.* article.
- Verify the package name is free on PyPI before committing to it.
- **Publishing / credential boundary:** Claude Code can build the sdist/wheel, do a **TestPyPI**
  dry-run install, and set up a GitHub Actions release workflow using **PyPI Trusted Publishing
  (OIDC)** so no long-lived token exists. Do **not** hand it a PyPI API token; the human performs
  the final publish / authorizes the trusted publisher.

### Suggested repo layout
```
<pkg>/
  pyproject.toml
  README.md
  LICENSE                       # GPLv3
  src/<pkg>/__init__.py
  src/<pkg>/units.py            # per-gas constants (transcribed from MD.cpp)
  src/<pkg>/engine.py           # numba kernels: forces, velocity Verlet, walls
  src/<pkg>/api.py              # run(), sweep(), result object
  benchmark.py
  tests/
    oracle/                     # MD.cpp + build script + stdin driver
    test_units.py
    test_energy_conservation.py
    test_ideal_gas_limit.py
    test_vs_c_oracle.py
  notebooks/
    reproduce_figure4.ipynb     # the acceptance test
```

---

## 7. Decisions to confirm with the domain expert

- Final numerical tolerances for oracle agreement and the ideal-gas limit (§5).
- Whether `n_equil` should default to 0 (faithful to the paper) — recommended yes.
- Package name.
- Whether the shipped trajectory format needs to stay VMD-compatible or can be simplified
  (the browser app already covers live visualization, so the Python package can focus on
  thermodynamics).

