# Annular Quantum Well Backend Testing and Validation

## 1. Overview

This package implements and validates a two-dimensional particle in an
infinitely deep annular well, the ring-shaped counterpart to the triangular
well project. The allowed region is `r_inner < r < r_outer`; the wavefunction
vanishes on both circular boundaries.

The backend calculates radial eigenvalues, normalized wavefunctions,
probability densities, and energy levels. Automated tests check inputs,
boundary conditions, energy scaling, normalization, and plot output.

## 2. Files

- `annulusbackend.py` — Main annular infinite-well calculations and plotting.
- `test_annulus.py` — Automated tests for inputs, calculations, and output files.
- `verify_annulus_normalization.py` — Independent Gauss-Legendre check of
  wavefunction normalization.
- `README.md` — Model, setup, testing procedure, and limitations.
- `requirements.txt` — Python package dependencies.

## 3. Model and user inputs

The inputs are inner radius `r_inner`, outer radius `r_outer`, angular quantum
number `l` (nonnegative integer), and radial quantum number `q` (positive
integer). The radii must satisfy `0 < r_inner < r_outer`.

Separation in polar coordinates gives the radial Bessel solution

```text
R(r) = A [J_l(k r) Y_l(k r_inner) - Y_l(k r) J_l(k r_inner)]
```

This construction makes `R(r_inner)=0`; the allowed radial wavenumbers `k`
are roots that also make `R(r_outer)=0`. The energy is

```text
E = hbar^2 k^2 / (2 mass)
```

The wavefunction normalization includes the polar area element
`r dr dtheta`, so the area integral of `|psi|^2` over the annulus equals 1.
For `l > 0`, cosine and sine angular states have the same energy (twofold
degeneracy). `run_annulus` plots the real cosine state. The `wavefunction`
helper also allows the sine and complex angular states.

The backend returns `energy` using the `mass` and `hbar` arguments. Their
defaults are both 1, so the default result is in dimensionless units unless
consistent physical units are supplied.

## 4. Requirements and installation

Python 3 and the packages in `requirements.txt` are required. Install them
from the project directory with:

```bash
python -m pip install -r requirements.txt
```

The backend uses NumPy, SciPy, and Matplotlib. The verifier and tests also use
NumPy. `math`, `pathlib`, and `unittest` are part of Python's standard library.

## 5. Generate the plots

Run the example from the project directory:

```bash
python annulusbackend.py
```

It generates three PNG images in `generated_plots/` for `r_inner=0.35`,
`r_outer=1.0`, `l=1`, and `q=1`.

To supply different values from Python:

```python
from annulusbackend import run_annulus

result = run_annulus(
    r_inner=0.35,
    r_outer=1.0,
    angular_momentum=1,
    radial_index=2,
    output_dir="generated_plots",
)
print(result["energy"])
print(result["plots"])
```

The returned dictionary includes the inputs, radial wavenumber, energy,
normalization description, angular degeneracy, the energy-level record, and
absolute paths to the plots. It also contains JSON-serializable `wavefunction`
and `probability_density` objects. Each object has 2D `x`, `y`, and `values`
arrays; the arrays are sampled on the plotting grid. Values outside the
annulus are set to zero, and `plot_grid.annulus_mask` identifies those cells.
The returned wavefunction is the normalized real cosine state used in the
plots, and probability density is its squared magnitude. For a readable
formula, `normalized_wavefunction_equation` returns the numerically evaluated
coefficient and wave number in the Bessel-function expression;
`normalization_constant_A` is the coefficient for the radial bracket exactly
as written in `radial_solution`.

## 6. Run the tests

Run automated tests:

```bash
python test_annulus.py
```

They cover valid and invalid input, positive and ordered radial roots,
zero-valued boundary conditions, wavefunction components, numerical
normalization for sample states, energy consistency and inverse-square size
scaling, and expected plot filenames.

Run the independent normalization check:

```bash
python verify_annulus_normalization.py
```

It accepts `--r-inner`, `--r-outer`, `--l`, `--q`, and `--order`. For example:

```bash
python verify_annulus_normalization.py --r-inner 0.3 --r-outer 1 --l 2 --q 1
```

The result should be close to 1. Small differences can arise from root-finding
and quadrature precision.

## 7. Expected plot files

For the state `l=1`, `q=1`, the backend writes:

- `wavefunction_l1_q1_cos.png` — signed real wavefunction amplitude.
- `probability_density_l1_q1_cos.png` — nonnegative probability density
  `|psi|^2`.
- `energy_level_l1_q1.png` — selected state energy level.

The inner hole and the area outside the outer radius are masked in the spatial
plots. The energy plot shows the selected energy level.

## 8. Validation scope and limitations

The tests check that the code follows its implemented equations and selected
mathematical properties. The boundary check samples both circular walls; the
normalization check independently performs quadrature over the annular domain.
The energy test verifies the formula and expected inverse-square scaling with
size. These are implementation checks, not a complete independent derivation
of the quantum model.

The backend uses numerical root finding for the Bessel boundary condition.
The wavefunction radial normalization is evaluated numerically. Results are
intended for the nonrelativistic, two-dimensional, hard-wall model with
positive finite radii and integer quantum indices.

## 9. Integration with the shared project

The UI can collect the two radii and the `(l,q)` state, call `run_annulus`,
then display the returned numerical results and plot paths. Keep the function
arguments and returned dictionary stable while integrating this branch with
the rest of the team's shape backends.
