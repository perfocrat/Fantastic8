# Fantastic8
physics sem 1 project 
# Fantastic8 — Quantum Well Explorer

Fantastic8 is a Python project for exploring quantum states in two-dimensional
infinite-potential wells with different geometries. The Streamlit UI lets a user
choose a shape and geometry, select quantum numbers, and view calculated
wavefunctions, probability densities, and energy information.

> **Integration status:** The UI and shape solvers are present, but they do not
> yet share a consistent interface. The UI currently contains illustrative
> placeholder fields, and the ellipse and sector solvers rely on module-level
> settings. Do not treat the placeholder plots as physical results.

## Features and shape modules

| Shape shown in UI | Backend module | Current entry point / status |
|---|---|---|
| Circle | `shapes/circle.py` | `run_circle(...)` |
| Rectangle | `shapes/rectangle.py` | `run_rectangle(...)` |
| Equilateral Triangle | `shapes/equilateral_triangle.py` | `run_triangle(a, m, n, ...)` |
| Ellipse | `shapes/ellipse.py` | Needs a callable adapter that accepts geometry and quantum-state selection |
| Circular Sector | `shapes/sector.py` | Needs a callable adapter that accepts geometry and quantum-state selection |
| Annulus | `shapes/annulus.py` | `run_annulus(r_inner, r_outer, angular_momentum, radial_index, ...)` |

## Project structure

```text
Fantastic8/
├── config.py
├── requirements.txt
├── README.md
├── ui/
│   ├── __init__.py
│   └── interface.py
├── shapes/
│   ├── __init__.py
│   ├── circle.py
│   ├── rectangle.py
│   ├── equilateral_triangle.py
│   ├── ellipse.py
│   ├── sector.py
│   └── annulus.py
└── generated_plots/       # created at runtime
```

Keep virtual environments, `.git/`, caches, `.pyc` files, and generated plot files
out of the source repository unless there is a specific reason to commit them.

## Requirements

- Python 3.10 or newer recommended
- Dependencies are listed in `requirements.txt`

Install dependencies from the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run the UI

From the project root (the folder containing `config.py` and `requirements.txt`):

```bash
streamlit run ui/interface.py
```

Streamlit will print a local URL, usually `http://localhost:8501`.

## How the application should work

1. The user selects exactly one shape in the UI.
2. The UI validates the geometry and quantum numbers for that shape.
3. The UI calls only that shape's adapter/solver.
4. The solver returns a consistent result object containing numerical arrays,
   energy data, and optionally plot paths or Matplotlib figures.
5. The UI renders those results. It must not generate substitute placeholder
   physics fields when a solver fails.

## Standard result contract to use during integration

For easier integration, each adapter should return the same top-level keys where
possible:

```python
{
    "shape": "circle",
    "inputs": {"m": 1, "n": 1},
    "energy": 0.0,                 # selected state energy, in agreed units
    "wavefunction": array_or_data,  # selected state's normalized psi(x, y)
    "probability_density": array_or_data,  # abs(psi)**2
    "coordinates": {
        "x": array_or_data,
        "y": array_or_data,
        "mask": array_or_data
    },
    "plots": {
        "wavefunction": optional_path_or_figure,
        "probability_density": optional_path_or_figure,
        "energy_levels": optional_path_or_figure
    }
}
```

The values above are an interface example, not computed results. Each shape adapter
must map the shape-specific return format to this contract. Use `None` for an
unavailable optional plot and raise a clear `ValueError` for invalid input.

## Important physics consistency notes

- **Quantum-number conventions differ.** For example, the triangle model requires
  `1 <= n <= m <= 5`; the annulus uses angular momentum `l >= 0` and radial index
  `q >= 1`. Do not silently reinterpret an invalid pair.
- **Energy units must match.** `rectangle.py` currently defines SI values for
  Planck's reduced constant and electron mass, while circle and sector use
  dimensionless defaults. Choose one convention and document it across all solvers
  before comparing energies from different shapes.
- **Geometry parameters must be passed at runtime.** The ellipse and sector modules
  currently use module-level settings; refactor their calculation functions to
  accept the geometry selected in the UI rather than mutating globals per request.
- **One solver per request.** Selecting one shape should call only its solver.
- **Numerical results and placeholder visuals are not interchangeable.** Keep
  placeholder generation disabled once the actual backend is connected.

## Integration and testing checklist

- [ ] Confirm that every shape module imports without error.
- [ ] Define one adapter function per shape with a consistent input/output contract.
- [ ] Test each adapter directly before connecting it to Streamlit.
- [ ] Add input validation for geometry and shape-specific quantum numbers.
- [ ] Check normalization and boundary conditions where supported.
- [ ] Confirm `probability_density` is calculated from `abs(wavefunction)**2`.
- [ ] Confirm plots use the selected state and selected geometry.
- [ ] Test every shape in the UI, including invalid inputs.
- [ ] Confirm only the selected solver runs.
- [ ] Check `git status` and avoid committing `.venv`, `.git`, caches, and generated outputs.

## Team contribution guide

- Shape-module contributors: maintain the mathematical solver, input validation,
  numerical outputs, and solver-specific tests.
- UI contributor: maintain controls, geometry previews, and rendering of returned
  solver results.
- Integration lead: maintain the common adapter/result contract, route each shape
  selection to exactly one solver, resolve units and quantum-number mappings, and
  test the complete application.
