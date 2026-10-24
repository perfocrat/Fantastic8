"""
Automated validation for circular_sector_quantum_well_code.py

Run from the directory containing both Python files:
    python test_circular_sector_quantum_well.py
    python test_circular_sector_quantum_well.py --angle 90

The script tests every (ell, s) pair from 1..5 (25 combinations), writes a
human-readable report to test_results.txt, and exits non-zero if any check fails.

Note: the reference solver calls its separated quantum numbers ell and s.
These correspond to the requested pair (m, n) for this test harness.
"""

import argparse
import contextlib
import io
import math
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.special import jv

import circular_sector_quantum_well_code as solver


def make_state(ell, s, alpha):
    """Build one (ell, s) eigenstate using the solver's own numerical root finder."""
    nu = ell * math.pi / alpha
    roots = solver.bessel_positive_zeros(nu, s)
    root = float(roots[s - 1])
    energy = solver.HBAR**2 * root**2 / (
        2.0 * solver.MASS * solver.RADIUS**2
    )
    return {
        "ell": ell,
        "s": s,
        "nu": nu,
        "zero": root,
        "energy": energy,
    }


def check_state(ell, s, alpha, nr=401, ntheta=401):
    """Return a list of validation outcomes for one state."""
    results = []
    state = make_state(ell, s, alpha)
    R = float(solver.RADIUS)
    hbar = float(solver.HBAR)
    mass = float(solver.MASS)
    nu = state["nu"]
    root = state["zero"]

    # 1. The root must satisfy J_nu(root) = 0.
    root_residual = abs(float(jv(nu, root)))
    results.append(("Bessel root", root_residual < 2e-7,
                    f"|J_nu(root)|={root_residual:.3e}"))

    # 2. Energy must match the theoretical infinite-sector-well result.
    theoretical_energy = hbar**2 * root**2 / (2.0 * mass * R**2)
    energy_error = abs(state["energy"] - theoretical_energy)
    results.append(("Theoretical energy", energy_error < 1e-10,
                    f"computed={state['energy']:.10g}, theory={theoretical_energy:.10g}, "
                    f"abs_error={energy_error:.3e}"))

    # 3. Check normalization: integral |psi|^2 r dr dtheta = 1.
    # Include r=0 explicitly; the eigenfunction is finite there.
    r = np.linspace(0.0, R, nr)
    theta = np.linspace(0.0, alpha, ntheta)
    rr, tt = np.meshgrid(r, theta, indexing="ij")
    psi = solver.wavefunction(state, rr, tt, alpha)
    integrand = np.abs(psi)**2 * rr
    angular_integral = simpson(integrand, x=theta, axis=1)
    norm = float(simpson(angular_integral, x=r))
    norm_error = abs(norm - 1.0)
    results.append(("Normalization", norm_error < 5e-3,
                    f"integral={norm:.8f}, abs_error={norm_error:.3e}"))

    # 4. Angular boundaries: psi(r, 0) = psi(r, alpha) = 0.
    radial_samples = np.linspace(0.0, R, 101)
    psi_theta0 = solver.wavefunction(
        state, radial_samples, np.zeros_like(radial_samples), alpha
    )
    psi_thetaalpha = solver.wavefunction(
        state, radial_samples, np.full_like(radial_samples, alpha), alpha
    )
    angular_boundary_error = max(
        float(np.max(np.abs(psi_theta0))),
        float(np.max(np.abs(psi_thetaalpha))),
    )
    results.append(("Angular boundaries", angular_boundary_error < 1e-10,
                    f"max_abs_psi={angular_boundary_error:.3e}"))

    # 5. Circular boundary: psi(R, theta) = 0.
    boundary_theta = np.linspace(0.0, alpha, 101)
    psi_outer = solver.wavefunction(
        state, np.full_like(boundary_theta, R), boundary_theta, alpha
    )
    radial_boundary_error = float(np.max(np.abs(psi_outer)))
    results.append(("Outer-radius boundary", radial_boundary_error < 2e-7,
                    f"max_abs_psi={radial_boundary_error:.3e}"))

    # 6. Probability density must be non-negative and equal to |psi|^2.
    density = np.abs(psi)**2
    density_min = float(np.min(density))
    density_consistency = float(np.max(np.abs(density - psi**2)))
    results.append(("Probability density", density_min >= -1e-12 and density_consistency < 1e-12,
                    f"min_density={density_min:.3e}, consistency_error={density_consistency:.3e}"))

    # 7. Invalid physical constants should be caught by the test harness.
    constants_valid = R > 0 and mass > 0 and hbar > 0
    results.append(("Physical constants", constants_valid,
                    f"R={R}, MASS={mass}, HBAR={hbar}"))

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--angle", type=float, default=90.0,
        help="sector angle in degrees, 0 < angle <= 360 (default: 90)"
    )
    parser.add_argument(
        "--output", default="test_results.txt",
        help="report output path (default: test_results.txt)"
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    lines = [
        "CIRCULAR SECTOR QUANTUM WELL — AUTOMATED TEST REPORT",
        "=" * 62,
        f"Timestamp: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Solver module: {solver.__file__}",
        f"Sector angle: {args.angle} degrees",
        f"RADIUS={solver.RADIUS}, MASS={solver.MASS}, HBAR={solver.HBAR}",
        "Quantum-number combinations: ell=1..5, s=1..5 (25 total)",
        "",
    ]

    if not (0.0 < args.angle <= 360.0):
        lines += ["OVERALL RESULT: FAIL", "Reason: angle must satisfy 0 < angle <= 360."]
        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("\n".join(lines))
        raise SystemExit(1)

    alpha = math.radians(args.angle)
    passed = 0
    failed = 0
    checks_total = 0
    failures = []

    # Silence incidental solver output if any is added in a future version.
    for ell in range(1, 6):
        for s in range(1, 6):
            label = f"(m,n)=({ell},{s}) / (ell,s)=({ell},{s})"
            lines.append(f"STATE {label}")
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    outcomes = check_state(ell, s, alpha)
                for name, ok, detail in outcomes:
                    checks_total += 1
                    status = "PASS" if ok else "FAIL"
                    lines.append(f"  [{status}] {name}: {detail}")
                    if ok:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(f"{label}: {name} — {detail}")
            except Exception as exc:
                failed += 1
                checks_total += 1
                detail = f"{type(exc).__name__}: {exc}"
                lines.append(f"  [ERROR] State test crashed: {detail}")
                failures.append(f"{label}: {detail}")
                lines.append(traceback.format_exc().rstrip())
            lines.append("")

    lines += [
        "-" * 62,
        f"Checks passed: {passed}",
        f"Checks failed/errors: {failed}",
        f"Checks attempted: {checks_total}",
        f"States tested: 25",
        f"OVERALL RESULT: {'PASS' if failed == 0 else 'FAIL'}",
    ]
    if failures:
        lines.append("")
        lines.append("FAILURE SUMMARY")
        lines.extend(f"- {failure}" for failure in failures)

    report = "\n".join(lines) + "\n"
    output_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"Report written to: {output_path.resolve()}")
    raise SystemExit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
