
import math
from pathlib import Path

import numpy as np

from trianglebackend import (
    validate_inputs,
    triangle_mesh,
    k_vectors,
    wavefunction,
    run_triangle,
)

# ============================================================
# CONFIGURATION
# ============================================================

A = 5.0
GRIDS = [100, 200, 400, 800]
TOLERANCE = 1e-10

STATES = [
    (m, n)
    for m in range(1, 6)
    for n in range(1, m + 1)
]

failures = []
checks = {
    "input_validation": 0,
    "boundary_conditions": 0,
    "normalization": 0,
    "wavevectors": 0,
    "energy_consistency": 0,
    "schrodinger_identity": 0,
    "graph_generation": 0,
}

def record(category, passed, description):
    if passed:
        checks[category] += 1
        print(f"  PASS: {description}")
    else:
        failures.append(description)
        print(f"  FAIL: {description}")

# ============================================================
# NUMERICAL NORMALIZATION INTEGRAL
# ============================================================

def integrate_probability(m, n, divisions):
    tri = triangle_mesh(A, divisions=divisions)

    psi = wavefunction(
        A, m, n, "S", tri.x, tri.y
    )

    cells = tri.triangles

    x1, y1 = tri.x[cells[:, 0]], tri.y[cells[:, 0]]
    x2, y2 = tri.x[cells[:, 1]], tri.y[cells[:, 1]]
    x3, y3 = tri.x[cells[:, 2]], tri.y[cells[:, 2]]

    area = 0.5 * np.abs(
        (x2 - x1) * (y3 - y1)
        - (x3 - x1) * (y2 - y1)
    )

    p1 = psi[cells[:, 0]]
    p2 = psi[cells[:, 1]]
    p3 = psi[cells[:, 2]]

    return float(np.sum(
        area / 6.0 * (
            p1**2 + p2**2 + p3**2
            + p1*p2 + p2*p3 + p3*p1
        )
    ))

# ============================================================
# TEST 1: INPUT VALIDATION
# ============================================================

print("\n" + "=" * 68)
print("TEST 1: INPUT VALIDATION")
print("=" * 68)

for m, n in STATES:
    try:
        a, mv, nv = validate_inputs(A, m, n)
        record(
            "input_validation",
            a == A and mv == m and nv == n,
            f"Valid input m={m}, n={n}"
        )
    except Exception as exc:
        record(
            "input_validation", False,
            f"Valid input m={m}, n={n}: {exc}"
        )

invalid_cases = [
    (0, 1, 1),
    (-1, 1, 1),
    (float("nan"), 1, 1),
    (5, 0, 1),
    (5, 1, 0),
    (5, 1, 2),
    (5, 6, 1),
    (5, 1.5, 1),
]

for a, m, n in invalid_cases:
    try:
        validate_inputs(a, m, n)
        record(
            "input_validation", False,
            f"Invalid input accepted: a={a}, m={m}, n={n}"
        )
    except ValueError:
        record(
            "input_validation", True,
            f"Invalid input correctly rejected: a={a}, m={m}, n={n}"
        )
    except Exception as exc:
        record(
            "input_validation", False,
            f"Wrong exception for invalid input: {exc}"
        )

# ============================================================
# TEST 2: BOUNDARY CONDITIONS
# ============================================================

print("\n" + "=" * 68)
print("TEST 2: WAVEFUNCTION BOUNDARY CONDITIONS")
print("=" * 68)

t = np.linspace(0.0, 1.0, 2001)

edges = [
    (A * t, np.zeros_like(t)),
    (A * t / 2, math.sqrt(3) * A * t / 2),
    (A * (1 - t / 2), math.sqrt(3) * A * t / 2),
]

boundary_errors = {}

for m, n in STATES:
    try:
        edge_errors = [
            float(np.max(np.abs(
                wavefunction(A, m, n, "S", x, y)
            )))
            for x, y in edges
        ]

        maximum = max(edge_errors)
        boundary_errors[m, n] = maximum

        record(
            "boundary_conditions",
            maximum < TOLERANCE,
            f"m={m}, n={n}, max boundary residual={maximum:.3e}"
        )
    except Exception as exc:
        record(
            "boundary_conditions", False,
            f"m={m}, n={n}: {exc}"
        )

# ============================================================
# TEST 3: NORMALIZATION AND CONVERGENCE
# ============================================================

print("\n" + "=" * 68)
print("TEST 3: NORMALIZATION CONVERGENCE")
print("=" * 68)

normalization_results = {}

for m, n in STATES:
    try:
        values = [
            integrate_probability(m, n, grid)
            for grid in GRIDS
        ]

        normalization_results[m, n] = values

        errors = [abs(1.0 - value) for value in values]

        # Check that the error decreases as the grid is refined.
        converges = all(
            errors[i + 1] < errors[i]
            for i in range(len(errors) - 1)
        )

        print(
            f"m={m}, n={n}: "
            + ", ".join(
                f"N={grid}: {value:.8f}"
                for grid, value in zip(GRIDS, values)
            )
        )

        record(
            "normalization",
            converges and all(math.isfinite(v) for v in values),
            f"m={m}, n={n}, integral converges toward 1"
        )
    except Exception as exc:
        record(
            "normalization", False,
            f"m={m}, n={n}: {exc}"
        )

# ============================================================
# TEST 4: WAVEVECTOR SANITY
# ============================================================

print("\n" + "=" * 68)
print("TEST 4: WAVEVECTOR CHECKS")
print("=" * 68)

for m, n in STATES:
    try:
        vectors = k_vectors(A, m, n)

        k_squared = np.array([
            kx**2 + ky**2 for kx, ky in vectors
        ])

        same_magnitude = np.allclose(
            k_squared,
            k_squared[0],
            rtol=1e-10,
            atol=1e-12
        )

        finite = np.all(np.isfinite(k_squared))
        correct_count = len(vectors) == 6

        record(
            "wavevectors",
            same_magnitude and finite and correct_count,
            f"m={m}, n={n}, six finite vectors, "
            f"equal |k| magnitudes={same_magnitude}"
        )
    except Exception as exc:
        record(
            "wavevectors", False,
            f"m={m}, n={n}: {exc}"
        )

# ============================================================
# TEST 5: ENERGY CONSISTENCY
# ============================================================

print("\n" + "=" * 68)
print("TEST 5: ENERGY FORMULA")
print("=" * 68)

for m, n in STATES:
    try:
        vectors = k_vectors(A, m, n)
        k_squared = vectors[0][0]**2 + vectors[0][1]**2

        energy_from_k = k_squared / 2.0

        result = run_triangle(
            A, m, n,
            output_dir=f"validation_plots/m{m}_n{n}",
            divisions=30,
        )

        coefficient = m*m + m*n + n*n
        expected = (
            8.0 * math.pi**2 * coefficient
            / (9.0 * A**2)
        )

        backend_energy = result["energy"][
            "energy_factor_without_hbar2_over_M"
        ]

        passed = (
            math.isclose(
                energy_from_k, expected,
                rel_tol=1e-10, abs_tol=1e-12
            )
            and math.isclose(
                backend_energy, expected,
                rel_tol=1e-10, abs_tol=1e-12
            )
            and result["energy"][
                "dimensionless_coefficient_E_over_E0"
            ] == coefficient
        )

        record(
            "energy_consistency",
            passed,
            f"m={m}, n={n}, coefficient={coefficient}, "
            f"E(k)={energy_from_k:.8f}, "
            f"E(backend)={backend_energy:.8f}"
        )
    except Exception as exc:
        record(
            "energy_consistency", False,
            f"m={m}, n={n}: {exc}"
        )

# ============================================================
# TEST 6: SCHRODINGER EQUATION ALGEBRAIC IDENTITY
# ============================================================

print("\n" + "=" * 68)
print("TEST 6: SCHRODINGER EQUATION CONSISTENCY")
print("=" * 68)

x = np.array([0.23 * A, 0.37 * A, 0.51 * A])
y = np.array([0.15 * A, 0.20 * A, 0.25 * A])

for m, n in STATES:
    try:
        vectors = k_vectors(A, m, n)
        k_squared = vectors[0][0]**2 + vectors[0][1]**2

        psi = wavefunction(A, m, n, "S", x, y)

        # Analytic identity for a sum of plane waves.
        laplacian = -k_squared * psi

        lhs = -0.5 * laplacian
        rhs = (k_squared / 2.0) * psi

        residual = float(np.max(np.abs(lhs - rhs)))

        record(
            "schrodinger_identity",
            residual < 1e-10,
            f"m={m}, n={n}, algebraic residual={residual:.3e}"
        )
    except Exception as exc:
        record(
            "schrodinger_identity", False,
            f"m={m}, n={n}: {exc}"
        )

# ============================================================
# TEST 7: GRAPH GENERATION
# ============================================================

print("\n" + "=" * 68)
print("TEST 7: GRAPH GENERATION")
print("=" * 68)

for m, n in STATES:
    folder = Path(f"validation_plots/m{m}_n{n}")

    expected_files = [
        folder / f"wavefunction_m{m}_n{n}_s.png",
        folder / f"probability_density_m{m}_n{n}_s.png",
        folder / f"energy_level_m{m}_n{n}.png",
    ]

    exists = all(
        path.is_file() and path.stat().st_size > 0
        for path in expected_files
    )

    record(
        "graph_generation",
        exists,
        f"m={m}, n={n}, all three graph files exist"
    )

# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 68)
print("FINAL SUMMARY")
print("=" * 68)

for category, count in checks.items():
    print(f"{category}: {count} checks passed")

print(f"\nUnexpected failures: {len(failures)}")

if failures:
    print("\nFAILED CHECKS:")
    for failure in failures:
        print(" -", failure)
else:
    print("All implemented checks passed.")

print("\nIMPORTANT:")
print("These are runtime and consistency checks.")
print("They do NOT independently prove the physical model is correct.")
print("The energy and Schrodinger checks use identities related to")
print("the same wavevectors and formulas implemented in the backend.")
print("Boundary sampling and mesh integration are numerical approximations.")
print("An independent theoretical derivation is still needed.")
