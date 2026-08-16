## Python package (`noblegasmd`)

A `pip install`-able, numba-accelerated Python port of the same NVE Lennard-Jones engine,
runnable with zero build step (e.g. in Google Colab). `MD.cpp` below remains the authoritative
physics reference; the port is validated against it (see `tests/`).

```bash
pip install noblegasmd
```

```python
from noblegasmd import run, sweep

result = run(gas="Ar", T=300.0, rho=40.0)
print(result.Z, result.P_avg, result.T_avg)

df = sweep(gas="Ar", T=[100, 200, 300, 400], rho=[1, 40, 500, 2000], n_replicates=3)
```

See [`notebooks/reproduce_figure4.ipynb`](notebooks/reproduce_figure4.ipynb) to regenerate the
published quasi-isotherms, and [`MD_python_port_spec.md`](MD_python_port_spec.md) for the full
port specification. Development install: `pip install -e ".[dev]"`, then `pytest`.

---

## C++ reference implementation (`MD.cpp`)

- Source code for Molecular Dynamics Program - can compile and run on Linux, Windows, and Mac OSX

- More information about this program, including detailed instructions for its use, can be found [here for instructions](https://pubs.acs.org/doi/suppl/10.1021/acs.jchemed.7b00747) and [here for discussion of its use in an undergraduate laboratory setting](https://pubs.acs.org/doi/pdf/10.1021/acs.jchemed.7b00747)

- This folder should containt the source code (MD.cpp) and a makefile which can be used to compile the source to a machine-executable file.

- A copy of the gnu public license (LICENSE.md) should be included.

- To compile this code using a gnu C/C++ compiler and create an executable called 'MD.exe' in a Linux/Unix environment, type
  
  `make`
  
- To run the program in a Linux/Unix environment, type
  
  `./MD.exe` 

- The program will run interactively.  Follow the prompts to customize your simulation

*Note for Windows users on Cygwin installation:*  Our J. Chem. Ed. article suggests installing **all** Cygwin packages, which is quite large and time consuming.  Success with a much lighter installation has been reported 
by selecting the default package installation plus 3 additional packages under the "devel" sub-heading.
As accessed on 01/08/2020, these packages and version numbers are as follows:

- gcc-g++ (7.4.0-1) 
- git (2.21.0-1)
- make (4.2.1-2)



