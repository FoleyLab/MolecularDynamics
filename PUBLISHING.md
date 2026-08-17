# Publishing `noblegasmd`

Claude Code will not touch PyPI credentials or trigger a real publish. This
document covers what's already prepared and what a human needs to do.

## Already done

- `pyproject.toml` builds a pure-Python wheel (`py3-none-any`) via hatchling.
- `python -m build` produces `dist/noblegasmd-<version>-py3-none-any.whl` and
  the sdist; both pass `twine check`.
- A clean-venv install of the built wheel was verified locally (`pip install
  dist/*.whl` in a fresh venv, then `import noblegasmd; noblegasmd.run(...)`).
- `.github/workflows/release.yml` builds the package and publishes via
  **PyPI Trusted Publishing (OIDC)** — no API token is stored in this repo.

## What a human needs to do

1. **Register the trusted publisher** (one-time, per target):
   - On [test.pypi.org](https://test.pypi.org) and/or
     [pypi.org](https://pypi.org), create (or claim) the `noblegasmd` project,
     then under its "Publishing" settings add a trusted publisher pointing at:
     - Repository: `jayfoleyiv/MolecularDynamics` (or wherever this ends up)
     - Workflow: `.github/workflows/release.yml`
     - Environment: `testpypi` or `pypi` (must match the workflow's
       `environment:` field for that job)
   - Create matching GitHub Environments named `testpypi` and `pypi` under
     the repo's Settings > Environments (optionally with required reviewers
     for `pypi`, as a manual publish gate).
2. **TestPyPI dry run**: Actions tab -> "Publish to PyPI" -> "Run workflow" ->
   target = `testpypi`. Then verify with:
   ```bash
   pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ noblegasmd
   ```
   (`--extra-index-url` is needed because `numpy`/`numba`/`pandas` aren't on
   TestPyPI.)
3. **Real publish**: push a tag `vX.Y.Z` (or run the workflow manually with
   target = `pypi`) once you're satisfied.

## Local dry run (no GitHub Actions, no upload)

```bash
python -m pip install --upgrade build twine
python -m build                      # writes dist/*.whl and dist/*.tar.gz
python -m twine check dist/*         # metadata sanity check, no network
python -m venv /tmp/noblegasmd_check
/tmp/noblegasmd_check/bin/pip install dist/noblegasmd-*-py3-none-any.whl
/tmp/noblegasmd_check/bin/python -c "from noblegasmd import run; print(run(gas='Ar', T=300, rho=40, n_steps=500, seed=1).Z)"
```
