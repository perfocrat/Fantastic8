"""
Test harness for circle_backend.py in the Fantastic8 project.

Run from the Fantastic8 project root:
    python3 test_circle_backend.py

Tests the backend's actual API:
    run_circle(R, n_states, grid_size, mass, hbar, make_plots, ...)

The backend returns Cartesian 2D grids and energy-sorted states, not a
one-dimensional radial grid or a solve(m, n) function.
"""

import math
import traceback
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros

import circle_backend


# Keep these aligned with the project's chosen units and backend defaults.
RADIUS = 1.0
MASS = 1.0
HBAR = 1.0
GRID_SIZE = 160
TARGET_M = range(1, 6)
TARGET_N = range(1, 6)

# The backend sorts states by energy and counts C/S angular degeneracy
# separately. A larger count is needed to include all (m,n) targets through 5.
N_STATES = 120

# Tolerances account for numerical integration on a finite Cartesian grid.
ENERGY_RTOL = 1e-9
NORMALIZATION_RTOL = 0.06
DENSITY_RTOL = 1e-8


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def expected_energy(m, n):
    alpha = float(jn_zeros(m, n)[-1])
    return HBAR**2 * alpha**2 / (2.0 * MASS * RADIUS**2)


def integrate_cartesian(values, mask, dx, dy):
    """Approximate integral over the circular domain using the Cartesian grid."""
    safe_values = np.where(mask, values, 0.0)
    return float(np.sum(safe_values) * dx * dy)


def run_tests():
    report = []
    passed_cases = 0
    failed_cases = 0
    individual_passes = 0
    individual_failures = 0

    report.extend([
        "CIRCULAR INFINITE-WELL BACKEND TEST REPORT",
        "=" * 48,
        f"Backend: {Path(circle_backend.__file__).name}",
        f"Assumptions: 2D circular infinite well; R={RADIUS}, mass={MASS}, hbar={HBAR}",
        f"Grid: {GRID_SIZE} x {GRID_SIZE}; requested states: {N_STATES}",
        "Inputs tested: m=1..5, n=1..5 (25 unique quantum-number combinations)",
        "Energy theory: E_mn = hbar^2 * alpha_mn^2 / (2 M R^2), J_m(alpha_mn)=0",
        "Normalization is checked numerically over the Cartesian grid; a finite-grid tolerance is used.",
        "",
    ])

    # Call the backend once so every state uses the same state list and grid.
    try:
        result = circle_backend.run_circle(
            R=RADIUS,
            n_states=N_STATES,
            grid_size=GRID_SIZE,
            mass=MASS,
            hbar=HBAR,
            make_plots=False,
            show_plots=False,
        )
    except Exception:
        report.append("[FATAL] circle_backend.run_circle(...) failed:")
        report.append(traceback.format_exc())
        output = "\n".join(report)
        Path("circle_test_results.txt").write_text(output, encoding="utf-8")
        print(output)
        return 1

    required_top_keys = {"coordinates", "states", "plots", "probability_densities"}
    for key in required_top_keys:
        try:
            check(key in result, f"Backend result is missing top-level key '{key}'")
            individual_passes += 1
        except AssertionError as exc:
            individual_failures += 1
            report.append(f"[FAIL] API: {exc}")

    coords = result.get("coordinates", {})
    for key in ("x", "y", "r", "theta", "mask"):
        try:
            check(key in coords, f"coordinates is missing '{key}'")
            check(np.asarray(coords[key]).shape == (GRID_SIZE, GRID_SIZE),
                  f"coordinates['{key}'] has unexpected shape {np.asarray(coords[key]).shape}")
            individual_passes += 1
        except AssertionError as exc:
            individual_failures += 1
            report.append(f"[FAIL] Grid: {exc}")

    x = np.asarray(coords.get("x"))
    y = np.asarray(coords.get("y"))
    mask = np.asarray(coords.get("mask"), dtype=bool)
    dx = float(x[0, 1] - x[0, 0])
    dy = float(y[1, 0] - y[0, 0])
    states = result.get("states", [])

    # A state list contains separate cosine/sine angular states for m > 0.
    # Check each requested (m,n) pair once; either angular member is valid.
    for m in TARGET_M:
        for n in TARGET_N:
            report.append(f"[CASE] m={m}, n={n}")
            try:
                matches = [
                    state for state in states
                    if int(state.get("m", -1)) == m and int(state.get("n", -1)) == n
                ]
                check(bool(matches),
                      f"State (m={m}, n={n}) not returned; increase N_STATES if needed")
                state = matches[0]

                energy = float(state["energy"])
                theory = expected_energy(m, n)
                check(np.isfinite(energy), f"Energy is not finite: {energy}")
                check(math.isclose(energy, theory, rel_tol=ENERGY_RTOL, abs_tol=1e-10),
                      f"Energy mismatch: backend={energy:.12g}, theory={theory:.12g}")
                individual_passes += 1
                report.append(f"  [PASS] Energy: backend={energy:.8g}, theory={theory:.8g}")

                psi = np.asarray(state["wavefunction"])
                rho = np.asarray(state["probability_density"])
                check(psi.shape == mask.shape,
                      f"Wavefunction shape {psi.shape} does not match grid {mask.shape}")
                check(rho.shape == mask.shape,
                      f"Density shape {rho.shape} does not match grid {mask.shape}")
                check(np.all(np.isfinite(psi[mask])),
                      "Wavefunction contains non-finite values inside the circle")
                check(np.all(np.isfinite(rho[mask])),
                      "Probability density contains non-finite values inside the circle")
                individual_passes += 1
                report.append("  [PASS] Wavefunction/density arrays match the 2D grid and are finite inside")

                expected_rho = np.abs(psi) ** 2
                density_error = float(np.nanmax(np.abs(rho - expected_rho)))
                density_scale = max(1.0, float(np.nanmax(expected_rho)))
                check(density_error <= DENSITY_RTOL * density_scale,
                      f"Probability density does not match |psi|^2 (max error {density_error:.3g})")
                individual_passes += 1
                report.append("  [PASS] Probability density equals |wavefunction|^2")

                norm = integrate_cartesian(np.abs(psi) ** 2, mask, dx, dy)
                check(np.isfinite(norm), f"Normalization integral is not finite: {norm}")
                check(abs(norm - 1.0) <= NORMALIZATION_RTOL,
                      f"Numerical normalization integral={norm:.6f}; expected approximately 1.0")
                individual_passes += 1
                report.append(f"  [PASS] Numerical normalization integral={norm:.6f}")

                # Backend represents the infinite-well boundary by mask and NaNs
                # outside the circle; it does not return a potential array.
                check(np.all(np.isnan(psi[~mask])),
                      "Wavefunction should be NaN outside the circular mask")
                check(np.all(np.isnan(rho[~mask])),
                      "Probability density should be NaN outside the circular mask")
                individual_passes += 1
                report.append("  [PASS] Wavefunction and density are masked outside the well")

                report.append(f"  [INFO] Angular type checked: {state.get('angular_type', 'not provided')}")
                passed_cases += 1
            except Exception as exc:
                failed_cases += 1
                report.append(f"  [FAIL] {type(exc).__name__}: {exc}")
            report.append("")

    # Check plot generation separately, without cluttering the numerical test run.
    report.append("[PLOT TEST]")
    try:
        plot_result = circle_backend.run_circle(
            R=RADIUS,
            n_states=6,
            grid_size=GRID_SIZE,
            mass=MASS,
            hbar=HBAR,
            output_dir="generated_plots_test",
            make_plots=True,
            show_plots=False,
        )
        plots = plot_result.get("plots", {})
        expected_plot_keys = {"shape", "energy_levels", "wavefunctions", "probability_density"}
        missing = expected_plot_keys - set(plots)
        check(not missing, f"Missing plot entries: {sorted(missing)}")
        missing_files = [name for name, filename in plots.items() if not Path(filename).is_file()]
        check(not missing_files, f"Plot files were not created for: {missing_files}")
        individual_passes += 1
        report.append(f"  [PASS] All expected plot files were generated: {', '.join(sorted(expected_plot_keys))}")
    except Exception as exc:
        individual_failures += 1
        report.append(f"  [FAIL] {type(exc).__name__}: {exc}")
        report.append(traceback.format_exc())

    report.extend([
        "",
        "=" * 48,
        f"Unique (m,n) cases passed: {passed_cases}/25",
        f"Unique (m,n) cases failed: {failed_cases}/25",
        f"Individual checks passed: {individual_passes}",
        f"Individual checks failed: {individual_failures}",
        "",
        "Notes:",
        "- This harness uses the backend's actual run_circle() API and Cartesian grids.",
        "- The backend returns cosine/sine angular variants for m > 0; each (m,n) pair is tested once.",
        "- Normalization is approximated by a Cartesian-grid integral, so a finite-grid tolerance is used.",
        "- This script checks the infinite-well mask; the backend does not return a numeric potential array.",
        "- Review any failed checks before merging the backend into main.",
    ])

    output = "\n".join(report)
    Path("circle_test_results.txt").write_text(output, encoding="utf-8")
    print(output)
    print("\nSaved report to: " + str(Path("circle_test_results.txt").resolve()))
    return 0 if failed_cases == 0 and individual_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run_tests())
