"""
Comprehensive tests for a rectangular infinite potential well backend.

Tests:
1. Energy values against the analytical infinite-well equation
2. Wavefunction against the analytical solution
3. Probability density equals |psi|^2
4. Normalization over the rectangular region
5. Boundary conditions
6. Orthogonality of distinct states
7. Time-independent Schrödinger equation residual inside the well
8. Full backend run and plot generation

Usage:
    python3 test_rectangle_backend_complete.py "/path/to/backend_file"

This script assumes the backend provides:
    energy(a, b, m, n)
    psi(a, b, m, n, x, y)
    probability_density(a, b, m, n, x, y)
    run_rectangle(a, b, num_states=..., grid=..., output_dir=...)

Important:
- These tests verify the rectangular INFINITE-well model.
- They cannot prove that the physical model is appropriate for a different system.
- A numerical test cannot independently verify the derivation of Schrödinger's equation.
  It can check whether the implemented wavefunction satisfies the time-independent
  Schrödinger equation for the assumed potential inside the well.
"""

import importlib.util
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys

import numpy as np


# ---------- Test settings ----------
A = 5e-9                    # Rectangle width, metres
B = 4e-9                    # Rectangle height, metres
GRID_POINTS = 401
QUANTUM_NUMBERS = range(1, 6)
HBAR = 1.054571817e-34       # J s
ELECTRON_MASS = 9.1093837015e-31  # kg
ENERGY_RTOL = 1e-12
WAVEFUNCTION_RTOL = 1e-10
NORMALIZATION_TOL = 1e-4
ORTHOGONALITY_TOL = 1e-4
SCHRODINGER_RELATIVE_TOL = 2e-3


def load_backend(file_path):
    """Import a Python backend from the supplied file path."""
    file_path = Path(file_path)
    if not file_path.is_file():
        raise FileNotFoundError(f"Backend file not found: {file_path}")

    spec = importlib.util.spec_from_file_location(
        "rectangle_backend_under_test", file_path
    )
    if spec is None or spec.loader is None:
        loader = SourceFileLoader(
            "rectangle_backend_under_test", str(file_path)
        )
        spec = importlib.util.spec_from_loader(
            "rectangle_backend_under_test", loader
        )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def report(name, passed, details=""):
    status = "PASS" if passed else "FAIL"
    suffix = f" — {details}" if details else ""
    print(f"[{status}] {name}{suffix}")
    return bool(passed)


def analytical_energy(a, b, m, n):
    return (
        np.pi**2 * HBAR**2 / (2 * ELECTRON_MASS)
    ) * (m**2 / a**2 + n**2 / b**2)


def analytical_psi(a, b, m, n, x, y):
    return (
        2.0 / np.sqrt(a * b)
        * np.sin(m * np.pi * x / a)
        * np.sin(n * np.pi * y / b)
    )


def integrate_2d(values, x, y):
    """Integrate a grid shaped (len(y), len(x)) over x and y."""
    integrate_x = np.trapezoid(values, x=x, axis=1)
    return np.trapezoid(integrate_x, x=y, axis=0)


def run_tests(backend):
    a, b = A, B
    x = np.linspace(0.0, a, GRID_POINTS)
    y = np.linspace(0.0, b, GRID_POINTS)
    X, Y = np.meshgrid(x, y)

    total_pass = 0
    total_fail = 0

    def check(name, condition, details=""):
        nonlocal total_pass, total_fail
        if report(name, condition, details):
            total_pass += 1
        else:
            total_fail += 1

    print("\n" + "=" * 76)
    print("RECTANGULAR INFINITE POTENTIAL WELL — COMPREHENSIVE TESTS")
    print(f"Rectangle: a={a:.3e} m, b={b:.3e} m")
    print(f"Grid: {GRID_POINTS} x {GRID_POINTS}")
    print("=" * 76)

    # 1–4. Energy, wavefunction, density and normalization for all 25 states.
    print("\n1–4. Energy, wavefunction, probability density and normalization")
    state_functions = {}

    for m in QUANTUM_NUMBERS:
        for n in QUANTUM_NUMBERS:
            label = f"(m={m}, n={n})"
            expected_e = analytical_energy(a, b, m, n)
            actual_e = backend.energy(a, b, m, n)
            check(
                f"Energy {label}",
                np.isclose(actual_e, expected_e, rtol=ENERGY_RTOL, atol=0.0),
                f"backend={actual_e:.8e} J, theory={expected_e:.8e} J",
            )

            expected = analytical_psi(a, b, m, n, X, Y)
            actual = np.asarray(backend.psi(a, b, m, n, X, Y))
            state_functions[(m, n)] = actual

            check(
                f"Wavefunction {label}",
                actual.shape == expected.shape
                and np.allclose(
                    actual, expected, rtol=WAVEFUNCTION_RTOL, atol=1e-8
                ),
            )

            density = np.asarray(
                backend.probability_density(a, b, m, n, X, Y)
            )
            check(
                f"Probability density {label}",
                density.shape == actual.shape
                and np.allclose(
                    density, np.abs(actual) ** 2,
                    rtol=1e-10, atol=1e-8,
                ),
            )

            norm = float(integrate_2d(np.abs(actual) ** 2, x, y))
            check(
                f"Normalization {label}",
                np.isclose(norm, 1.0, rtol=NORMALIZATION_TOL,
                           atol=NORMALIZATION_TOL),
                f"integral={norm:.8f}",
            )

    # 5. Boundary conditions: psi should vanish on all four edges.
    print("\n5. Boundary conditions")
    for state, psi_grid in state_functions.items():
        edge_max = max(
            float(np.max(np.abs(psi_grid[0, :]))),
            float(np.max(np.abs(psi_grid[-1, :]))),
            float(np.max(np.abs(psi_grid[:, 0]))),
            float(np.max(np.abs(psi_grid[:, -1]))),
        )
        wavefunction_scale = max(float(np.max(np.abs(psi_grid))), 1e-300)
        boundary_relative_error = edge_max / wavefunction_scale
        check(
            f"Boundary conditions state {state}",
            boundary_relative_error < 1e-12,
            f"max boundary |psi|={edge_max:.3e}, "
            f"relative boundary error={boundary_relative_error:.3e}",
        )

    # 6. Orthogonality: integral psi_mn * psi_pq = 0 for distinct states.
    print("\n6. Orthogonality of distinct states")
    states = list(state_functions.keys())
    for i, state_1 in enumerate(states):
        psi_1 = state_functions[state_1]
        for state_2 in states[i + 1:]:
            psi_2 = state_functions[state_2]
            overlap = float(integrate_2d(psi_1 * psi_2, x, y))
            check(
                f"Orthogonality {state_1} vs {state_2}",
                abs(overlap) < ORTHOGONALITY_TOL,
                f"overlap={overlap:.3e}",
            )

    # 7. Check the time-independent Schrödinger equation inside the well:
    # -(hbar^2/2m) * Laplacian(psi) = E*psi where V=0.
    # Finite differences approximate second derivatives, so use a tolerance.
    print("\n7. Time-independent Schrödinger equation inside the well")
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    for m in QUANTUM_NUMBERS:
        for n in QUANTUM_NUMBERS:
            psi_grid = state_functions[(m, n)]
            energy_value = float(backend.energy(a, b, m, n))

            d2x = (
                psi_grid[1:-1, 2:] - 2 * psi_grid[1:-1, 1:-1]
                + psi_grid[1:-1, :-2]
            ) / dx**2
            d2y = (
                psi_grid[2:, 1:-1] - 2 * psi_grid[1:-1, 1:-1]
                + psi_grid[:-2, 1:-1]
            ) / dy**2
            laplacian = d2x + d2y
            psi_inner = psi_grid[1:-1, 1:-1]

            lhs = -(HBAR**2 / (2 * ELECTRON_MASS)) * laplacian
            rhs = energy_value * psi_inner
            denominator = max(float(np.max(np.abs(rhs))), 1e-300)
            relative_residual = float(np.max(np.abs(lhs - rhs)) / denominator)

            check(
                f"Schrödinger equation state ({m}, {n})",
                relative_residual < SCHRODINGER_RELATIVE_TOL,
                f"relative finite-difference residual={relative_residual:.3e}",
            )
    # 8. Full backend run and output plots.
    print("\n8. Full run and plot generation")
    output_dir = Path("rectangle_test_plots")
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = backend.run_rectangle(
            a, b, num_states=6, grid=80, output_dir=str(output_dir)
        )
        expected_keys = {
            "a", "b", "states", "X", "Y", "fields",
            "energy_plot", "psi_plot", "rho_plot",
        }
        missing = expected_keys - set(result.keys())
        check(
            "run_rectangle returns expected fields",
            not missing,
            f"missing keys={sorted(missing)}" if missing else "",
        )

        for key in ("energy_plot", "psi_plot", "rho_plot"):
            plot_path = Path(result[key])
            check(
                f"Plot generated: {key}",
                plot_path.is_file(),
                str(plot_path),
            )
    except Exception as exc:
        check(
            "run_rectangle completes",
            False,
            f"{type(exc).__name__}: {exc}",
        )

    print("\n" + "=" * 76)
    print(f"TEST SUMMARY: {total_pass} passed, {total_fail} failed")
    print("=" * 76)
    print(
        "\nScope note: These tests check implementation against the standard "
        "rectangular infinite-well model. They cannot prove the derivation "
        "of Schrödinger's equation or establish that this model suits every "
        "physical system."
    )
    return total_fail == 0


def main():
    if len(sys.argv) != 2:
        print(
            'Usage: python3 test_rectangle_backend_complete.py '
            '"/path/to/rectangle_backend.py"'
        )
        raise SystemExit(2)

    try:
        backend = load_backend(sys.argv[1])
    except Exception as exc:
        print(f"Could not load backend: {type(exc).__name__}: {exc}")
        raise SystemExit(2)

    success = run_tests(backend)
    raise SystemExit(0 if success else 1)


if __name__ == "__main__":
    main()
