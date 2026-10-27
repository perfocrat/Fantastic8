# Fantastic8 Physics Verification Report

**Audit scope:** read-only verification of the square/rectangular, circular,
and equilateral-triangle modules against the supplied paper PDF.

**Primary source:** `reference paper for physics project.pdf`, a 46-page PDF
whose first page identifies arXiv `quant-ph/0307070v1`, dated 9 July 2003.
The supplied SciSpace Markdown was treated only as a secondary guide and
checked against the PDF.

**Status vocabulary**

- **PAPER-SUPPORTED:** stated in the supplied PDF.
- **INDEPENDENTLY DERIVED:** follows from the paper or standard separation of
  variables, but is not presented there as a solved result.
- **NUMERICALLY VERIFIED:** calculated in the isolated audit checks.
- **NOT VERIFIED:** no safe test was run or the supplied source cannot support
  the claim.
- **NOT APPLICABLE:** not part of this module's implementation.

## 1. Executive summary

The three modules were identified unambiguously:

1. [`shapes/rectangle.py`](../shapes/rectangle.py)
2. [`shapes/circle.py`](../shapes/circle.py)
3. [`shapes/equilateral_triangle.py`](../shapes/equilateral_triangle.py)

The square, circular, and triangular modules implement **stationary
eigenstates**, not the paper's complete Gaussian wave-packet construction,
time evolution, revival analysis, or classical orbit analysis. The
implementation uses analytic formulas and special functions; it does not
numerically solve the 2-D time-independent Schrödinger PDE on a grid.

Ten selected numerical cases were executed in an isolated audit script:
four rectangle cases (three square and one unequal rectangle), three
equilateral-triangle cases, and three circle cases. **All 10 passed the
stated comparison criteria; zero failed.** The criterion is defined in
Section 7 and does not imply full scientific certification of every code path.
Numerical normalization integrals were within approximately
`1.4e-14` of unity or better. The largest sampled wall residual among these
cases was below `9e-15` for the triangle and below `5e-16` for the rectangle;
circle Bessel-root residuals were below `4e-16`.

The most important qualification is the triangle's indexing convention. The
paper uses `(M,N)` with `M >= 2N`; the module uses `(m,n)` with `m >= n`.
The states correspond under `M=m+n`, `N=n`. With this mapping, tested energies
and wavefunctions match the paper, including the special non-degenerate state
and the two-dimensional degenerate eigenspaces. Without documenting this
mapping, the module/UI labels can mislead readers into thinking they are
paper-convention quantum numbers.

The circle module numbers radial roots from `n=1`, while the paper's radial
node label starts at `n_r=0`; thus `n_code=n_r+1`. The tested circular energies
and normalized real states agree with independently bracketed Bessel roots
under that mapping.

The rectangle implementation is an analytical **general rectangle**
extension. The paper's worked result is for a **square**; it explicitly poses
the general rectangular extension as Problem 1. No claim is made here that the
paper directly derives the unequal-side module result.

### Overall module assessment

| Module | Assessment |
|---|---|
| Rectangle | Formula, wavefunction, normalization, density, and wall behavior are consistent with the paper's square case and an independently derived unequal-rectangle extension for the tested states. Source uses fixed SI constants and has no direct configurable mass/`hbar` interface. |
| Circle | Tested Bessel roots, energy, normalization, real angular modes, and boundary behavior agree. Source's one-based radial index differs from the paper's zero-based radial-node label. |
| Equilateral triangle | Tested mode formulas, normalization, boundary behavior, and eigenspace membership agree after the documented index and coordinate transformations. Source emits energy metadata/factors, not a numeric physical energy value. |
| Project as a paper reproduction | Incomplete: no Gaussian expansion coefficients, time propagation, autocorrelation, revival analysis, or paper-specific wave-packet calculations are implemented in these three modules. |

## 2. Research source and extraction audit

### 2.1 What the PDF itself verifies

The PDF title page (PDF p. 1) shows the title *Wave packet construction in
two-dimensional quantum billiards: Blueprints for the square, equilateral
triangle, and circular cases*, authors M. A. Doncheski, S. Heppelmann,
R. W. Robinett, and D. C. Tussey, and the arXiv identifier/date above. The
abstract describes Gaussian wave-packet calculations and square, triangle,
and circular infinite-well geometries.

The supplied PDF does **not** show a journal citation, volume/pages, or DOI in
its extracted pages; searching the extracted PDF text found no `10.1119` DOI
or journal metadata string. The journal citation and DOI stated in the audit
instructions therefore could not be independently confirmed from this file.
No external source was used to fill that gap.

### 2.2 Verified paper equations and locations

| Paper result | PDF location | Verification |
|---|---|---|
| 1-D infinite-well energies and normalized sine eigenfunction | Eqs. (19)–(20), p. 7 | Verified: `E_n = hbar^2 pi^2 n^2/(2 mu a^2)` and `u_n = sqrt(2/a) sin(n pi x/a)`. |
| Numerical convention `hbar = 2 mu = a = 1`, example widths, and use of 40 lowest states | Eq. (25) and text, p. 8 | Verified. The 40-state truncation is for the 1-D packet example, not a universal truncation for all geometries. |
| Square energy and product eigenstate | Eqs. (35a)–(35b), p. 12 | Verified. The geometry is explicitly `Lx x Ly = a x a`, a square. |
| General rectangle extension | Problem 1, p. 31; incommensurate-side discussion p. 12 | Verified as an extension/problem, not as the paper's worked stationary spectrum. |
| Triangle vertices, energy, and allowed-label inequality | Eq. (51), p. 16 | Verified: vertices `(0,0)`, `(a/2,sqrt(3)a/2)`, `(-a/2,sqrt(3)a/2)`; integral labels with `m >= 2n`. |
| Triangle degenerate eigenfunctions | Eqs. (52)–(53), pp. 16–17 | Verified; two symmetry states for `m>2n`. |
| Triangle special state and symmetry | Eqs. (54)–(56), p. 17 | Verified; `m=2n` has one non-degenerate state; the symmetric limiting form is related by `sqrt(2)` and the antisymmetric one vanishes. |
| Circular potential, Bessel states, and angular momentum | Eqs. (68)–(69), p. 20 | Verified: `V=0` for `r<R`, `infinity` for `r>=R`; order `|m|`; `m` is signed. |
| Circular root condition, spectrum, radial-node label, and asymptotic estimate | Eqs. (70)–(72), p. 21 | Verified: the paper's `n_r` starts at zero; Eq. (72) is explicitly an approximation, not an exact root formula. |
| Normalized circular eigenstate and radial integral | Eqs. (81)–(84), p. 23 | Verified. The angular factor is `exp(i m theta)/sqrt(2 pi)` and the radial normalization uses `r dr`. |
| Circular packet parameter example | Text associated with Fig. 5, p. 25 | Verified: `y0=0`, `p0x=0`, `p0y=100`, `R=1`, with `x0/R=0.0,-0.5,0.70`. |
| Square packet example values | Eq. (25), Fig. 4 text, p. 14 | Verified: `p0=400 pi`, center `(a/2,a/2)`, widths `0.05`, and stated time scales approximately `0.8e-3`, `2.5e-3`, and `0.64`. |
| Square-orbit Table I | p. 32 | Verified: `T_cl^(po)/tau=sqrt(p^2+q^2)`, `tau=2a/v0`, and `tan(theta)=q/p`. |
| Circular-orbit Table II | p. 33 | Verified: `(p,q)=(2,1)` gives `L/R=4.00` and `Rmin/R=0`; `(3,1)` gives `L/R=5.20` and `Rmin/R=0.50`. |

The paper discusses the Schrödinger equation and says that the triangle
functions satisfy the 2-D equation, but the inspected PDF does not print the
full 2-D time-independent operator in the common Cartesian or polar PDE form.
The extraction's qualification on that point is supported.

### 2.3 SciSpace Markdown findings

The extraction's mathematical claims checked above have **no definite
equation-number, page-number, sign, coefficient, or arithmetic error detected**
in this audit. Its square/triangle/circle scope distinction is especially
important and is correct.

| Extraction claim/provenance | Audit finding |
|---|---|
| Paper studies square, triangle, and circle packet dynamics, not a grid-based stationary PDE solver | **Supported** by the abstract and the paper's methods and examples. |
| Paper's square result is not a complete arbitrary-rectangle result | **Supported**: p. 12 presents `a x a`; p. 31 poses a rectangular extension. |
| Triangle restriction `m>=2n`, with special case at equality | **Supported** by Eq. (51) p. 16 and Eqs. (52)–(56) pp. 16–17. |
| Circle radial label begins at zero, and Bessel asymptotic is not exact | **Supported** by Eqs. (70)–(72), pp. 21–22. |
| The exact Bessel-root values in the requested test cases | The extraction does not attribute exact decimal roots to the paper; it correctly says the paper has no raw exact-root table. The values used in this audit were independently bracketed and solved, not taken as paper tabulations. |
| Formal journal citation/DOI | **Not verifiable from the supplied PDF**; the file identifies the arXiv preprint but contains no journal/DOI text found by inspection. |
| Extraction header's `/home/sandbox/...` source path | **Provenance not verifiable** from the PDF; it is not the local path of the supplied file. This does not invalidate its equations, but the extraction history cannot be checked from file content. |

## 3. Source modules and actual implementation

### 3.0 Dependencies and helper-file inventory

The three shape modules do not import a project-local Python helper module.
Their numerical and plotting helpers are defined in the same file. The
separate project adapter [`shapes/adapters.py`](../shapes/adapters.py) calls
these modules; it is not imported by any of the three audited shape files.

| Module | Imports relevant to calculation |
|---|---|
| `shapes/rectangle.py` | NumPy, Matplotlib, `pathlib.Path`, and `math.pi/sqrt`. |
| `shapes/circle.py` | NumPy, Matplotlib, `pathlib.Path`, SciPy `jv` and `jn_zeros`. |
| `shapes/equilateral_triangle.py` | `math`, NumPy, Matplotlib, `matplotlib.tri`, and `pathlib.Path`. |

### 3.1 Square / rectangle

**File:** [`shapes/rectangle.py`](../shapes/rectangle.py)  
**Relevant functions:** `solve_schrodinger`, `energy`, `lowest_states`,
`psi`, `probability_density`, `rectangle_mesh`, `verify`,
`fields_on_mesh`, `plot_energy_levels`, `plot_grid`, `run_rectangle`.

- **Inputs and coordinates:** `run_rectangle(a,b,num_states=6,grid=80,
  output_dir="generated_plots")` accepts positive side lengths, a positive
  number of states, and a mesh size. It uses `x in [0,a]`, `y in [0,b]`.
  The selected-state indices are not arguments to `run_rectangle`; `energy`
  and `psi` accept `m,n` directly. `lowest_states` enumerates both indices
  from 1 through 20. The module does not validate `m,n` in `energy`/`psi`.
- **Potential and boundary:** `solve_schrodinger()` prints the zero-inside,
  infinite-outside potential and Dirichlet walls; it performs no numerical
  solve. The sine factors make `psi` zero on the four edges.
- **Energy:** `energy()` analytically evaluates
  `pi^2 hbar^2/(2 MASS) * (m^2/a^2+n^2/b^2)`. The hard-coded constants are
  SI `HBAR` and electron `MASS`; the code has no mass or `hbar` arguments.
- **Wavefunction and density:** `psi()` returns
  `sqrt(4/(a*b))*sin(m*pi*x/a)*sin(n*pi*y/b)`. Its normalization has
  dimensions of inverse length. `probability_density()` returns
  `abs(psi)**2`, dimension inverse area.
- **Numerical/plot behavior:** `verify()` uses an inclusive uniform grid and
  plain `sum(probability)*dx*dy`, then prints the estimate. `run_rectangle()`
  calculates a selected collection of lowest-state arrays and saves energy,
  wavefunction, and density PNGs. It returns arrays and file paths; it does
  not return a potential grid. The `run_rectangle()` plot/print entry point
  was not executed in this audit because it writes into its configured
  directory; its pure formula functions were tested instead.

**Scientific comparison:** for `a=b`, the formula agrees with paper Eq. (35a)
and the product of Eq. (20). For `a!=b`, the same result follows by applying
1-D separation independently to widths `a` and `b`; this unequal-side formula
is **INDEPENDENTLY DERIVED**, not directly given as the paper's solved
rectangle case. The code calls the geometry “rectangle” even in the square
benchmark, which is fine as an extension, but paper comparisons must identify
the equal-side square case.

### 3.2 Circle

**File:** [`shapes/circle.py`](../shapes/circle.py)  
**Relevant functions:** `_validate`, `_lowest_states`, `_grid`,
`_normalization`, `_psi`, `_plot_shape`, `_plot_energy`, `_plot_states`,
`run_circle`.

- **Inputs and coordinates:** `run_circle(R=1,n_states=6,grid_size=120,
  mass=1,hbar=1,output_dir=...,make_plots=True,show_plots=False)` uses
  center-origin polar coordinates internally and a Cartesian grid for outputs.
  `R`, mass, and `hbar` must be positive finite; grid size must be >=30.
  The energy inputs' physical units are caller-defined and must be consistent.
- **Indices:** `_lowest_states` uses angular order `m>=0`, radial index `n>=1`,
  and C/S real modes for nonzero `m`. The paper uses signed angular
  `m=0,+/-1,...` and radial-node index `n_r=0,1,...`. Correspondence is
  `|m_paper|=m_code`, `n_code=n_r+1`. For nonzero `m`, code's cosine/sine
  modes are real linear combinations of the paper's degenerate `+m` and `-m`
  complex states.
- **Potential and boundary:** the circular infinite wall is described by the
  model. `_grid()` marks points `r<=R`; `_psi()` masks exterior points to NaN.
  At `r=R`, the Bessel root enforces zero amplitude.
- **Energy and state calculation:** `_lowest_states()` obtains zeros through
  `scipy.special.jn_zeros`, sorts by `alpha^2`, and retains the first
  `n_states`, counting C/S degeneracy. Energy is
  `hbar^2*alpha^2/(2*mass*R^2)`.
- **Wavefunction and normalization:** `_psi()` evaluates `J_m(alpha*r/R)`
  times `1`, `cos(m theta)`, or `sin(m theta)`. `_normalization()` uses the
  analytic Bessel-root radial integral and the angular integral (`2 pi` for
  `m=0`, `pi` for real nonzero modes), so the result is normalized over
  `r dr dtheta`. `run_circle()` uses `abs(psi)**2`.
- **Returns/plots:** `run_circle()` returns state dictionaries, NumPy arrays
  for coordinates/fields, energy lists, optional paths, and optional
  Matplotlib figures. With `make_plots=False`, as used in this audit, it does
  not create the output directory or plot files. This module has no
  time-dependent packet coefficients or a potential-value array.

**Scientific comparison:** equations and real-state conventions are consistent
with Eqs. (68)–(70) and (81)–(82) after the radial-index shift. The paper's
exact Bessel roots are not tabulated; roots below were independently bracketed
from `J_m(z)=0` using `brentq`, distinct from the module's `jn_zeros` call.

### 3.3 Equilateral triangle

**File:** [`shapes/equilateral_triangle.py`](../shapes/equilateral_triangle.py)  
**Relevant functions:** `validate_inputs`, `triangle_mesh`, `k_vectors`,
`wavefunction`, `_save_state_plot`, `_save_energy_level_plot`,
`run_triangle`.

- **Inputs and coordinates:** `run_triangle(a,m,n,output_dir=...,divisions=120)`
  validates a positive finite side, integer indices, and its own index rule
  `1<=n<=m<=5`. It evaluates in a triangle with vertices
  `(0,0)`, `(a,0)`, `(a/2,sqrt(3)a/2)`. This is congruent to the paper triangle.
- **Paper-to-module index mapping:** paper labels `(M,N)` map to module labels
  `(m,n)=(M-N,N)`, equivalently `M=m+n`, `N=n`. Then
  `M^2+N^2-MN = m^2+mn+n^2`, and `M>=2N` is equivalent to `m>=n`.
  The module rule is therefore a transformed convention, not the paper's
  literal `m>=2n`. The module's bound `m<=5` is also not a paper upper bound.
- **Potential/boundary:** the module documents `potential_inside=0` in its
  result metadata, but does not calculate a spatial potential array or
  explicitly define the outside potential. The six-plane-wave combination
  vanishes at the edges of its triangular domain for tested modes.
- **Energy:** `run_triangle()` forms `m^2+mn+n^2` and a factor
  `8*pi^2/(9*a^2)`. Its dictionary returns a dimensionless coefficient,
  symbolic energy strings, and `energy_factor_without_hbar2_over_M`; it does
  **not** return a standalone numeric energy in joules. The full equation is
  equivalent to paper Eq. (51) under the label map.
- **Wavefunction and normalization:** `wavefunction()` evaluates six
  sinusoidal/cosinusoidal plane-wave terms with fixed signs. `run_triangle()`
  currently generates only the `"S"` family; for module `m==n` it applies a
  special normalization factor and rejects the unsupported `"C"` family.
  Under `M=m+n,N=n`, module `m==n` is the paper's special `M=2N` case.
- **Plots and output:** `run_triangle()` writes selected-state wavefunction,
  probability, and energy PNGs. It returns phase/equation/energy metadata but
  no numerical field arrays and omits the locally constructed plot-path
  dictionary. The plotted probability is `values**2`, equivalent to
  `|psi|^2` for this real wavefunction. Numerical test plots were directed
  only to `physics_audit_report/audit_artifacts/triangle_plots/`.

**Coordinate transformation for comparison:** source coordinates map to the
paper convention by `x_paper=x_module-a/2`,
`y_paper=sqrt(3)*a/2-y_module`. The Jacobian magnitude is one. The source
triangle's orientation therefore does not itself indicate a physics error.

**Scientific comparison:** the tested source sine state for `(module m,n)=(1,1)`
is the negative of the paper Eq. (54) state after the coordinate map; a global
minus sign is physically irrelevant. For module `(2,1)` and `(3,1)`, the
paper-convention modes are `(M,N)=(3,1)` and `(4,1)`. The source states
project into the paper's respective two-state degenerate subspaces with
numerical projection residual zero at the displayed precision. Those results
do not establish equivalence for every allowed input.

## 4. Independently derived benchmark formulas

### Square / rectangle

Paper-supported square (`a x a`) formula, Eq. (35a):

`E = hbar^2*pi^2*(nx^2+ny^2)/(2*mu*a^2)`.

Under the paper convention `hbar=1`, `2*mu=1`, `a=1`, the derived values
are:

| `(nx,ny)` | Reference energy |
|---|---:|
| `(1,1)` | `2*pi^2 = 19.739208802178716` |
| `(1,2)` | `5*pi^2 = 49.348022005446786` |
| `(2,2)` | `8*pi^2 = 78.95683520871486` |

These are **INDEPENDENTLY DERIVED** from Eq. (35a), not separately tabulated
numbers in the paper.

For unequal sides, the independently derived result is
`E = hbar^2*pi^2/(2*mu)*(nx^2/Lx^2+ny^2/Ly^2)`.
For `Lx=1`, `Ly=2`, `(nx,ny)=(1,1)` under the same dimensionless convention,
this is `1.25*pi^2 = 12.337005501361697`. This is not a paper-stated solved
rectangle result.

The module fixes SI `hbar` and electron mass. Therefore, its actual numeric
outputs for geometry values 1 and 2 are SI-scaled values, not the paper's
dimensionless numbers. The audit compared actual module energy with the same
analytic equation using the module's constants; the dimensionless equation
ratios agree with the paper values above.

### Triangle

Paper Eq. (51):

`E = hbar^2/(2*mu*a^2)*(4*pi/3)^2*(M^2+N^2-M*N)`, with integer labels and
`M>=2N`.

For `a=1`, `hbar=1`, `2*mu=1`, the three test cases `(M,N)=(2,1),(3,1),(4,1)`
give `16*pi^2/9 * (3,7,13)` or approximately
`52.637890139143245`, `122.8217436580009`, and `228.09752393628736`.

The module return values tested were the reduced factors
`8*pi^2/9 * (m^2+mn+n^2)` with the `hbar^2/mass` factor excluded. These are
`26.318945069571622`, `61.41087182900045`, and `114.04876196814368`.
They are exactly half the paper-convention numeric energies for `mu=0.5`;
that factor is expected because the module's returned key explicitly omits
`hbar^2/M`. The full-energy string uses that factor, but the module does not
return a complete computed scalar energy.

### Circle

Paper Eqs. (68)–(70) imply `E=hbar^2*z^2/(2*mu*R^2)`, where `z` is an
ordered positive zero of `J_|m|`. For `R=1`, `hbar=1`, `mu=0.5`, independent
`brentq` roots and reference energies are:

| Paper `(m,n_r)` | Code `(m,n)` | Root `z` | Reference/computed energy |
|---|---|---:|---:|
| `(0,0)` | `(0,1)` | `2.4048255576958` | `5.7831859629468` |
| `(1,0)` | `(1,1)` | `3.8317059702075` | `14.681970642124` |
| `(0,1)` | `(0,2)` | `5.5200781102863` | `30.471262343662` |

These roots/energies are **NUMERICALLY VERIFIED** with an independent bracketed
root solve and the paper's energy formula; the PDF does not tabulate these
exact decimals.

## 5. Actual test results

The full machine-readable output, inputs, and notes are in
[`test_results.csv`](test_results.csv). The independent test source is
[`audit_artifacts/independent_checks.py`](audit_artifacts/independent_checks.py).

### 5.1 Results summary

| Group | Cases | Status | Summary |
|---|---:|---|---|
| Rectangle square/unequal rectangle | 4 | 4 PASS | Relative energy error 0; normalization integrals approximately `0.999999999999989`; tested boundary residuals below `5e-16`. |
| Triangle mapped paper modes | 3 | 3 PASS | Reduced-factor relative errors `0` to `1.25e-16`; integrals approximately `1.000000000000003`; wall residuals below `9e-15`. |
| Circle | 3 | 3 PASS | Energy relative error 0 at displayed double precision; independent radial quadrature integral approximately `1.000000000000002`; Bessel root residuals below `4e-16`. |
| **Total** | **10** | **10 PASS, 0 FAIL** | No overall “physics correctness percentage” is claimed. |

Audit acceptance criterion: energy/scalar relative error <=`1e-10` when
comparable, normalization error <=`1e-10`, sampled wall/root residual
<=`1e-10`, and normalized wavefunction comparison error <=`1e-10`. These
are explicit audit thresholds for these analytic formulas and high-order
quadratures, not tolerances specified by the paper.

### 5.2 Numerical method and comparison details

- **Rectangle:** 100-point Gauss–Legendre product quadrature over the actual
  rectangle domain; source `psi` and `probability_density` compared with
  independently written sine-product references. Boundary samples covered
  each of the four walls. The module's own coarse `verify()` function was not
  run.
- **Triangle:** 90-point Gauss–Legendre quadrature in `y` and across the
  geometry-dependent `x` interval; paper Eqs. (52)–(54) were transcribed as an
  independent reference. The special state used sign-aligned L2 and density
  comparisons. Degenerate cases were compared to the full paper eigenspace,
  not incorrectly required to equal either chosen basis function. Sampled
  points on all three boundaries.
- **Circle:** Bessel roots bracketed and solved with `brentq` using brackets
  around the benchmark zeros, rather than calling the module's `jn_zeros`
  root routine for the reference. The radial integral was independently
  computed by Gauss–Legendre quadrature and the angular integral by a
  2048-point periodic sum. Compared the code's cosine state with the
  normalized real combination of the paper's degenerate complex `+m/-m`
  eigenstates.
- Triangle image outputs were generated by the audited module into the
  dedicated audit artifact subdirectory. No comparison figure was generated.
- Circle testing used `run_circle(..., make_plots=False)`; it made no plot
  files. Rectangle `run_rectangle()` was not executed because it writes plot
  files and prints reports; its pure energy, wavefunction, and density
  functions were tested.

### 5.3 Scope of PASS statuses

The pass count refers only to these 10 parameterized benchmark cases.
Triangle scalar tests check the source's reduced energy factor, **not** a full
numeric physical energy return. Triangle wavefunction comparisons verify one
non-degenerate state and two tested degenerate eigenspaces. They do not prove
every possible mode is correct. The circle comparison checks three states
and the `C` real mode; the code's `S` mode was not separately compared.

## 6. Discrepancies and limitations ranked by severity

### High: paper's wave-packet physics is absent

The project modules calculate stationary eigenstates and energies. The paper's
main focus includes Gaussian expansion coefficients, time-dependent phases
`exp(-iEt/hbar)`, autocorrelation/revival behavior, and semiclassical
periodicities. These are not implemented in the three audited shape modules.
If the product is described as reproducing the paper's packet results, that
claim is unsupported by the current implementation.

### Medium: quantum-number conventions need explicit mapping

- Triangle source labels are transformed: `(M,N)_paper=(m+n,n)_module`.
  The UI labels source values `m,n`, without inherently making them paper
  values. A raw comparison using same-named pairs is invalid.
- Circle source radial index `n` begins at 1; paper `n_r` begins at 0.
  A raw comparison without `n_code=n_r+1` compares different roots.

### Medium: energy return and units are inconsistent across shapes

- Rectangle returns SI joules using fixed electron constants.
- Circle permits arbitrary consistent `mass`/`hbar` inputs and defaults to
  `1.0`, dimensionless unless caller defines otherwise.
- Triangle source returns symbolic formulas and a reduced factor, not a numeric
  energy. It has no `mass` or `hbar` function arguments.
- The project adapter converts values, but this report audits the modules'
  direct function behavior as well. Do not compare raw outputs across modules
  without normalizing units.

### Medium: direct solver output and plots have limitations

- Rectangle `run_rectangle()` returns fields for the lowest-state collection,
  not a selected-state parameter. Its plot can therefore be unrelated to an
  arbitrary UI-selected state unless the selected state is among those
  returned.
- Triangle `run_triangle()` does not return its numerical wavefunction,
  density, coordinates, or plot paths despite saving selected-state PNGs.
- Circle returns NumPy arrays and, when plots are enabled, Matplotlib Figure
  objects; neither is directly JSON-serializable.
- No module returns a numerical potential grid.

### Low: verification coverage is finite

No grid-convergence study, broad state-space sweep, orthogonality matrix,
independent derivative-based Schrödinger residual, or full degeneracy-state
catalogue was run. The tested analytic formulas do not use a grid-based
eigenvalue solver, so finite-difference eigenvalue convergence is not
applicable; the plotted grids are visualization grids only.

## 7. Untested items and assumptions

- No full application/UI run or HTTP test was performed.
- No generated plot was visually inspected or compared pixel-by-pixel.
- No triangle states beyond the three mapped paper benchmarks were tested.
- No rectangle state indices above 2, arbitrary aspect ratios beyond `1x2`, or
  `lowest_states()` ordering/degeneracy sweep were tested.
- No circle angular modes above `m=1`, higher radial roots beyond the second
  `m=0` root, or Cartesian-grid integration error sweep was tested.
- Orthogonality between different states was not numerically tested.
- Potential output is NOT APPLICABLE to current field comparisons because the
  modules do not calculate/return a potential grid; only model descriptions
  or metadata are present.
- The independent triangle paper formulas were transcribed from the PDF
  equations and evaluated over the mapped domain. The test only covers the
  selected low modes.
- The formal journal citation/DOI cannot be confirmed from the supplied PDF.

## 8. Testing artifacts and reproducibility

Created audit artifacts:

- [`physics_verification_report.md`](physics_verification_report.md)
- [`human_explanation.txt`](human_explanation.txt)
- [`test_results.csv`](test_results.csv)
- [`audit_artifacts/independent_checks.py`](audit_artifacts/independent_checks.py)
- Nine triangle PNG outputs in
  [`audit_artifacts/triangle_plots/`](audit_artifacts/triangle_plots/)

The script has a guard that refuses to overwrite an existing CSV or triangle
plot output directory. Its initial version completed the rectangle, triangle,
and circle test cases with the existing `.venv` and
`PYTHONDONTWRITEBYTECODE=1`. After aligning the triangle special state's
irrelevant global sign and switching the independent circle reference from
`jn_zeros` to bracketed `brentq` plus quadrature, those affected metrics were
recomputed in targeted no-write checks and the CSV was updated. The finalized
script was not rerun as a whole because its overwrite guard correctly refuses
to replace the outputs already created by the initial run. No packages were
installed and no existing source was edited.

**Safety disclosure:** in the course of confirming that the PDF was
text-extractable, a temporary page-text file was created at
`/tmp/fantastic8-paper.txt`; Matplotlib also created its configured cache under
`/tmp/f8-matplotlib-config` during isolated plot generation. These are outside
the requested audit folder and were not part of the project source. They were
left untouched; no cleanup/deletion operation was run.

The pre-existing project working tree already contained untracked
`UI_BACKEND_INTEGRATION_REPORT.md`; the supplied PDF also appears as an
untracked workspace file. Both were left unchanged. No tracked project source
was modified, and no Git history/state-changing operation or package
installation was performed.
