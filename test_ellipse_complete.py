"""
Complete black-box test suite for Fantastic8's ellipse solver.

This suite tests the solver without editing its source file. It tests the
valid UI quantum-number combinations:
    m = 0..5
    n = 0..5, mapped to radial root index j = n + 1
    C family for m = 0..5
    S family for m = 1..5 (S, m=0 is not a valid family)

IMPORTANT:
- Keep this file separate from the teammate's solver.
- The solver must be importable as a normal Python module.
- If the solver is named ellipse_testcopy.py, use ELLIPSE_MODULE=ellipse_testcopy.
- The test suite adaptively increases the root-search range; it does not assume
  every state has a root below k*a = 6.
- Tests verify the implementation's observable behavior; passing tests do not
  prove that every possible mathematical defect is impossible.

Install dependencies in the project's active virtual environment:
    python -m pip install numpy scipy sympy matplotlib pytest

Run all tests:
    MPLBACKEND=Agg ELLIPSE_MODULE=ellipse_testcopy python -m pytest -v test_ellipse_complete.py

Save a test log:
    MPLBACKEND=Agg ELLIPSE_MODULE=ellipse_testcopy python -m pytest -v --tb=short test_ellipse_complete.py 2>&1 | tee ellipse_test_log.txt
"""
import importlib
import os
from functools import lru_cache

import numpy as np
import pytest

MODULE_NAME = os.environ.get("ELLIPSE_MODULE", "ellipse")
solver = importlib.import_module(MODULE_NAME)

M_VALUES = range(6)
N_VALUES = range(6)


def finite(values):
    return bool(np.all(np.isfinite(np.asarray(values))))


def close(actual, expected, atol=1e-8, rtol=1e-7):
    assert np.allclose(actual, expected, atol=atol, rtol=rtol), (
        f"actual={actual!r}, expected={expected!r}, "
        f"atol={atol}, rtol={rtol}"
    )


def valid_families(m):
    return ("C",) if m == 0 else ("C", "S")


@lru_cache(maxsize=None)
def roots_for_state_family(kind, m):
    """
    Find enough roots for all six UI radial indices in one call per (kind, m).
    Increase the search ceiling if fewer than six roots are found.
    """
    kmax = 6.0
    last_roots = []
    while kmax <= 100.0:
        last_roots = list(solver.roots_for(kind, m, kmax))
        if len(last_roots) >= 6:
            return tuple(float(q) for q in last_roots)
        kmax *= 1.5
    raise AssertionError(
        f"Found only {len(last_roots)} roots for kind={kind}, m={m} "
        f"after searching up to k*a={kmax / 1.5:.3f}; need 6"
    )


@pytest.fixture(scope="module")
def state_catalog():
    """Map each valid (family, m, UI n) to its q root and energy."""
    catalog = {}
    for m in M_VALUES:
        for kind in valid_families(m):
            roots = roots_for_state_family(kind, m)
            assert finite(roots), f"Non-finite roots for {kind}, m={m}"
            assert all(q > 0 for q in roots), f"Non-positive roots for {kind}, m={m}"
            assert np.all(np.diff(roots) > 0), (
                f"Roots are not strictly increasing for {kind}, m={m}: {roots}"
            )
            for n in N_VALUES:
                j = n + 1
                q = roots[j - 1]
                catalog[(kind, m, n)] = {
                    "kind": kind,
                    "m": m,
                    "n": n,
                    "j": j,
                    "q": q,
                    "energy": 4.0 * q / solver.ECC**2,
                }
    return catalog


# ---------- Module/API and fixed geometry ----------

def test_module_imports():
    assert solver is not None


def test_required_public_api_exists():
    required = (
        "solve_schrodinger", "lowest_states", "roots_for", "char_value",
        "fourier", "angular", "radial", "norm_const", "to_elliptic",
        "psi", "verify", "ellipse_mesh", "fields_on_mesh",
        "plot_energy_levels", "plot_grid", "print_table",
    )
    missing = [name for name in required
               if not callable(getattr(solver, name, None))]
    assert not missing, f"Missing required callable(s): {missing}"


def test_geometry_constants_and_ellipse_relations():
    required = ("SEMI_MAJOR", "SEMI_MINOR", "B_OVER_A", "ECC", "C_FOC", "XI0")
    missing = [name for name in required if not hasattr(solver, name)]
    assert not missing, f"Missing geometry constant(s): {missing}"

    for name in required:
        assert np.isfinite(getattr(solver, name)), f"{name} is not finite"

    assert solver.SEMI_MAJOR > solver.SEMI_MINOR > 0
    close(solver.B_OVER_A, solver.SEMI_MINOR / solver.SEMI_MAJOR)
    close(solver.ECC**2, 1.0 - solver.B_OVER_A**2)
    close(solver.C_FOC, solver.ECC)
    close(np.tanh(solver.XI0), solver.B_OVER_A)


def test_symbolic_schrodinger_and_coordinate_checks():
    h2, k2 = solver.solve_schrodinger()
    assert h2 is not None
    assert k2 is not None


# ---------- Coordinate conversion ----------

def test_coordinate_conversion_scalar_input():
    xi, eta = solver.to_elliptic(0.1, 0.1)
    assert np.isfinite(xi) and np.isfinite(eta)
    assert xi >= -1e-10


def test_coordinate_conversion_array_input_and_reconstruction():
    ratio = solver.B_OVER_A
    x = np.array([0.0, 0.05, 0.15, -0.2, 0.3, -0.1])
    y = np.array([0.0, 0.1, -0.05, 0.12, -0.1, -0.2 * ratio])
    assert np.all(x*x + (y/ratio)**2 < 1.0)

    xi, eta = solver.to_elliptic(x, y)
    assert np.shape(xi) == x.shape
    assert np.shape(eta) == y.shape
    assert finite(xi) and finite(eta)
    assert np.all(np.asarray(xi) >= -1e-10)

    xr = solver.C_FOC * np.cosh(xi) * np.cos(eta)
    yr = solver.C_FOC * np.sinh(xi) * np.sin(eta)
    close(xr, x, atol=2e-8, rtol=2e-8)
    close(yr, y, atol=2e-8, rtol=2e-8)


# ---------- Angular and radial separated solutions ----------

@pytest.mark.parametrize(
    "kind,m",
    [(kind, m) for m in M_VALUES for kind in valid_families(m)],
    ids=[f"{kind}-m{m}" for m in M_VALUES for kind in valid_families(m)],
)
def test_fourier_coefficients_and_angular_periodicity(kind, m):
    q = 0.35
    harmonics, coeffs = solver.fourier(kind, m, q)
    assert len(harmonics) > 0
    assert len(harmonics) == len(coeffs)
    assert finite(harmonics) and finite(coeffs)
    assert np.linalg.norm(coeffs) > 0

    eta = np.linspace(-np.pi, np.pi, 121)
    values = np.asarray(solver.angular(kind, m, q, eta))
    shifted = np.asarray(solver.angular(kind, m, q, eta + 2*np.pi))
    assert values.shape == eta.shape
    assert finite(values)
    close(values, shifted, atol=2e-7, rtol=2e-7)


@pytest.mark.parametrize("kind,m", [("C", 0), ("C", 1), ("S", 1), ("C", 2), ("S", 2)])
def test_radial_initial_conditions(kind, m):
    h = 1e-6
    r0 = float(np.asarray(solver.radial(kind, m, 0.2, np.array([0.0]))).ravel()[0])
    r1 = float(np.asarray(solver.radial(kind, m, 0.2, np.array([h]))).ravel()[0])
    if kind == "C":
        assert abs(r0 - 1.0) < 1e-8
        assert abs((r1 - r0) / h) < 2e-3
    else:
        assert abs(r0) < 1e-8
        assert abs((r1 - r0) / h - 1.0) < 2e-3


# ---------- Every valid UI quantum-number combination ----------
# There are 36 C-family states and 30 S-family states: 66 total.
# Each parameterized failure identifies the precise family, m, and UI n.

STATE_CASES = [
    (kind, m, n)
    for m in M_VALUES
    for kind in valid_families(m)
    for n in N_VALUES
]
STATE_IDS = [f"{kind}-m{m}-n{n}-j{n+1}" for kind, m, n in STATE_CASES]


@pytest.mark.parametrize("kind,m,n", STATE_CASES, ids=STATE_IDS)
def test_every_quantum_state_has_a_valid_root_and_energy(kind, m, n, state_catalog):
    state = state_catalog[(kind, m, n)]
    q = state["q"]
    energy = state["energy"]

    assert np.isfinite(q) and q > 0
    assert np.isfinite(energy) and energy > 0
    close(energy, 4.0*q/solver.ECC**2, atol=2e-7, rtol=2e-7)


@pytest.mark.parametrize("kind,m,n", STATE_CASES, ids=STATE_IDS)
def test_every_quantum_state_satisfies_radial_boundary(kind, m, n, state_catalog):
    state = state_catalog[(kind, m, n)]
    q = state["q"]
    boundary_value = np.asarray(
        solver.radial(kind, m, q, np.array([solver.XI0]))
    ).ravel()[0]
    assert np.isfinite(boundary_value)
    assert abs(float(boundary_value)) < 2e-5, (
        f"Radial boundary residual={boundary_value} for "
        f"{kind}, m={m}, UI n={n}, radial j={n+1}, q={q}"
    )


@pytest.mark.parametrize("kind,m,n", STATE_CASES, ids=STATE_IDS)
def test_every_quantum_state_has_valid_normalization_constant(kind, m, n, state_catalog):
    state = state_catalog[(kind, m, n)]
    norm = solver.norm_const(m, kind, state["q"])
    assert np.isfinite(norm) and norm > 0, (
        f"Invalid normalization constant for {kind}, m={m}, "
        f"UI n={n}, radial j={n+1}: {norm}"
    )


@pytest.mark.parametrize("kind,m,n", STATE_CASES, ids=STATE_IDS)
def test_every_quantum_state_wavefunction_is_finite_and_zero_on_boundary(
    kind, m, n, state_catalog
):
    state = state_catalog[(kind, m, n)]
    theta = np.linspace(0.0, 2.0*np.pi, 61)
    xb = np.cos(theta)
    yb = solver.B_OVER_A * np.sin(theta)
    values = np.asarray(solver.psi(m, kind, state["q"], xb, yb))
    assert values.shape == xb.shape
    assert finite(values)
    residual = float(np.max(np.abs(values)))
    assert residual < 2e-5, (
        f"Wavefunction boundary residual={residual} for "
        f"{kind}, m={m}, UI n={n}, radial j={n+1}"
    )


# ---------- Lowest states, normalization and theory verification ----------

def test_lowest_states_are_sorted_and_consistent():
    states = solver.lowest_states(6)
    assert len(states) == 6
    energies = np.asarray([s[0] for s in states], dtype=float)
    assert finite(energies)
    assert np.all(energies > 0)
    assert np.all(np.diff(energies) >= -1e-9)

    for energy, m, j, kind, q in states:
        assert m >= 0 and j >= 1
        assert kind in ("C", "S")
        assert not (kind == "S" and m == 0)
        assert q > 0 and np.isfinite(q)
        close(energy, 4.0*q/solver.ECC**2, atol=2e-7, rtol=2e-7)



def test_lowest_state_normalization_integral():
    _, m, j, kind, q = solver.lowest_states(1)[0]

    radial_values = solver.radial(
        kind, m, q, solver.XI_NODES
    )
    angular_values = solver.angular(
        kind, m, q, solver.ETA_NODES
    )

    area = solver.C_FOC**2 * (
        np.sinh(solver.XI_NODES)[:, None]**2
        + np.sin(solver.ETA_NODES)[None, :]**2
    )

    wave = (
        solver.norm_const(m, kind, q)
        * radial_values[:, None]
        * angular_values[None, :]
    )

    integral = np.sum(
        solver.XI_W[:, None]
        * solver.ETA_W
        * np.abs(wave)**2
        * area
    )

    assert np.isfinite(integral)
    assert abs(float(integral) - 1.0) < 2e-8, (
        f"Normalization integral={integral}"
    )



def test_builtin_verification_for_six_lowest_states():
    solver.verify(solver.lowest_states(6))


# ---------- Domain, mesh and plotting ----------

def test_wavefunction_and_probability_density_inside_domain():
    ratio = solver.B_OVER_A
    x = np.array([0.0, 0.1, -0.2, 0.3])
    y = np.array([0.0, 0.1, 0.1, -0.1])
    assert np.all(x*x + (y/ratio)**2 < 1.0)

    for _, m, j, kind, q in solver.lowest_states(6):
        wave = np.asarray(solver.psi(m, kind, q, x, y))
        density = np.abs(wave)**2
        assert wave.shape == x.shape
        assert finite(wave) and finite(density)
        assert np.all(density >= 0)


@pytest.mark.parametrize("grid_size", [5, 10, 20])
def test_mesh_and_fields_are_valid(grid_size):
    states = solver.lowest_states(3)
    mesh = solver.ellipse_mesh(grid_size)
    assert len(mesh.x) == len(mesh.y)
    assert len(mesh.x) > 1
    assert finite(mesh.x) and finite(mesh.y)
    assert np.all(np.asarray(mesh.x)**2 + (np.asarray(mesh.y)/solver.B_OVER_A)**2 <= 1 + 1e-9)

    fields = solver.fields_on_mesh(states, mesh)
    assert len(fields) == len(states)
    for field in fields:
        assert np.shape(field) == np.shape(mesh.x)
        assert finite(field)


def test_plotting_and_table_functions_run_without_display():
    import matplotlib.pyplot as plt

    states = solver.lowest_states(2)
    mesh = solver.ellipse_mesh(10)
    fields = solver.fields_on_mesh(states, mesh)

    solver.plot_energy_levels(states)
    fig_psi = solver.plot_grid(states, mesh, fields, "psi")
    fig_rho = solver.plot_grid(states, mesh, fields, "rho")
    assert fig_psi is not None and fig_rho is not None
    solver.print_table(states)
    plt.close("all")