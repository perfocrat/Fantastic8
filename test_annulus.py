"""Automated tests for the annular infinite-well backend."""
from pathlib import Path
import math
import unittest

import numpy as np

from annulusbackend import (
    boundary_residual,
    energy,
    radial_wavenumber,
    run_annulus,
    validate_inputs,
    wavefunction,
)


class AnnulusBackendTests(unittest.TestCase):
    def test_valid_inputs(self):
        self.assertEqual(validate_inputs(0.5, 1.0, 2, 3), (0.5, 1.0, 2, 3))

    def test_invalid_radii_and_indices(self):
        invalid_cases = [
            (0, 1, 0, 1), (-0.1, 1, 0, 1), (1, 1, 0, 1),
            (1.1, 1, 0, 1), (0.2, math.inf, 0, 1),
            (0.2, 1, -1, 1), (0.2, 1, 0, 0), (0.2, 1, True, 1),
            (0.2, 1, 1, 1.5),
        ]
        for args in invalid_cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                validate_inputs(*args)

    def test_radial_wavenumbers_are_positive_and_ordered(self):
        for ell in range(4):
            roots = [radial_wavenumber(0.3, 1.0, ell, q) for q in (1, 2, 3)]
            self.assertTrue(all(k > 0 for k in roots))
            self.assertTrue(roots[0] < roots[1] < roots[2])

    def test_boundary_conditions(self):
        for ell in range(4):
            for q in range(1, 4):
                residual = boundary_residual(0.3, 1.0, ell, q)
                self.assertLess(residual, 1e-9)

    def test_wavefunction_shape_and_angular_states(self):
        r = np.array([0.5, 0.7, 0.9])
        theta = np.array([0.0, 0.5, 1.0])
        for ell in range(4):
            real = wavefunction(0.3, 1.0, ell, 1, r, theta, "real")
            self.assertEqual(real.shape, r.shape)
            self.assertTrue(np.isfinite(real).all())
        with self.assertRaises(ValueError):
            wavefunction(0.3, 1.0, 0, 1, r, theta, "imag")
        self.assertTrue(np.isfinite(wavefunction(0.3, 1.0, 2, 1, r, theta, "complex")).all())

    def test_normalization_for_supported_sample_states(self):
        # Independent Gauss-Legendre quadrature in r and theta.
        nodes, weights = np.polynomial.legendre.leggauss(100)
        ri, ro = 0.3, 1.0
        r = (ro - ri) * (nodes + 1) / 2 + ri
        wr = weights * (ro - ri) / 2
        theta = math.pi * (nodes + 1)
        wt = math.pi * weights
        rr, tt = np.meshgrid(r, theta, indexing="ij")
        for ell, q in ((0, 1), (0, 2), (1, 1), (2, 1), (3, 2)):
            psi = wavefunction(ri, ro, ell, q, rr, tt)
            integral = np.sum(np.abs(psi) ** 2 * rr * wr[:, None] * wt[None, :])
            self.assertAlmostEqual(float(integral), 1.0, places=9)

    def test_energy_consistency_and_inverse_square_scaling(self):
        e1 = energy(0.3, 1.0, 1, 1)
        e2 = energy(0.6, 2.0, 1, 1)
        k = radial_wavenumber(0.3, 1.0, 1, 1)
        self.assertAlmostEqual(e1, k * k / 2.0, places=12)
        self.assertAlmostEqual(e1 / e2, 4.0, places=8)

    def test_graph_generation_and_expected_names(self):
        folder = Path(__file__).parent / "work" / "test_annulus_output"
        folder.mkdir(parents=True, exist_ok=True)
        result = run_annulus(0.3, 1.0, 0, 1, output_dir=folder, divisions=60)
        self.assertEqual(result["shape"], "annulus")
        expected = {
            "wavefunction": "wavefunction_l0_q1_cos.png",
            "probability_density": "probability_density_l0_q1_cos.png",
            "energy_level": "energy_level_l0_q1.png",
        }
        for kind, name in expected.items():
            path = Path(result["plots"][kind])
            self.assertEqual(path.name, name)
            self.assertTrue(path.is_file())

    def test_return_value_contains_ui_ready_data(self):
        folder = Path(__file__).parent / "work" / "test_annulus_ui_data"
        result = run_annulus(0.3, 1.0, 1, 1, output_dir=folder, divisions=40)
        self.assertIn("wavefunction", result)
        self.assertIn("probability_density", result)
        self.assertIn("energy_levels", result)
        psi = np.asarray(result["wavefunction"]["values"])
        probability = np.asarray(result["probability_density"]["values"])
        self.assertEqual(psi.shape, (40, 40))
        self.assertEqual(probability.shape, (40, 40))
        self.assertTrue(np.isfinite(psi).all())
        self.assertTrue(np.isfinite(probability).all())
        np.testing.assert_allclose(probability, psi ** 2, atol=1e-12)
        self.assertAlmostEqual(result["energy_levels"][0]["energy"], result["energy"])
        self.assertEqual(np.asarray(result["wavefunction"]["x"]).shape, (40, 40))
        self.assertEqual(np.asarray(result["wavefunction"]["y"]).shape, (40, 40))
        self.assertIn("J_1", result["normalized_wavefunction_equation"])
        self.assertGreater(result["normalization_constant_A"], 0)

    def test_normalization_constant_from_backend_matches_independent_integral(self):
        ri, ro, ell, q = 0.35, 1.0, 1, 1
        nodes, weights = np.polynomial.legendre.leggauss(160)
        r = (ro - ri) * (nodes + 1.0) / 2.0 + ri
        wr = weights * (ro - ri) / 2.0
        theta = math.pi * (nodes + 1.0)
        wt = math.pi * weights
        rr, tt = np.meshgrid(r, theta, indexing="ij")
        psi = wavefunction(ri, ro, ell, q, rr, tt)
        integral = np.sum(np.abs(psi) ** 2 * rr * wr[:, None] * wt[None, :])
        self.assertAlmostEqual(float(integral), 1.0, places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
