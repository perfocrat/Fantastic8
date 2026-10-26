import math
import json
import unittest
from unittest.mock import patch

import numpy as np

import server


class CalculateRequestTests(unittest.TestCase):
    def test_dispatches_one_shape_with_converted_inputs(self):
        result = {
            "shape": "Rectangle",
            "inputs": {"width": 8e-9, "height": 5e-9, "m": 2, "n": 1},
            "energy": 1.0e-21,
            "states": [{"m": 2, "n": 1, "E": 1.0e-21, "selected": True}],
            "wavefunction": np.ones((2, 2)),
            "probability_density": np.ones((2, 2)),
            "coordinates": {
                "x": np.ones((2, 2)),
                "y": np.ones((2, 2)),
                "mask": np.ones((2, 2), dtype=bool),
                "units": "m",
            },
        }
        payload = {
            "geometry": "rectangle",
            "dimensions": {"width": 8, "height": 5},
            "quantumNumbers": {"m": 2, "n": 1},
        }

        with patch("server.calculate_shape", return_value=result) as calculate:
            response = server.calculate_request(payload)

        calculate.assert_called_once_with(
            "rectangle", width=8e-9, height=5e-9, m=2, n=1, grid=server.GRID_SIZE
        )
        self.assertAlmostEqual(response["energy_mev"], 1.0e-21 / server.MEV_IN_JOULES)
        self.assertEqual(response["potential"]["outside_domain"], "infinite")

    def test_calculates_every_shape(self):
        requests = [
            ("rectangle", {"width": 8, "height": 5}, {"m": 2, "n": 1}),
            ("triangle", {"side": 7}, {"m": 2, "n": 1}),
            ("circle", {"radius": 4}, {"m": 0, "n": 1}),
            ("sector", {"radius": 4, "angle": 90, "startAngle": 30}, {"m": 1, "n": 1}),
            ("annulus", {"innerRadius": 2, "outerRadius": 5}, {"m": 0, "n": 1}),
            ("ellipse", {"semiMajor": 5, "semiMinor": 3}, {"m": 0, "n": 1}),
        ]
        for shape, dimensions, quantum_numbers in requests:
            with self.subTest(shape=shape):
                response = server.calculate_request(
                    {
                        "geometry": shape,
                        "dimensions": dimensions,
                        "quantumNumbers": quantum_numbers,
                    }
                )
                wavefunction = np.asarray(response["wavefunction"], dtype=float)
                density = np.asarray(response["probability_density"], dtype=float)
                mask = np.asarray(response["coordinates"]["mask"], dtype=bool)
                x = np.asarray(response["coordinates"]["x"], dtype=float)
                y = np.asarray(response["coordinates"]["y"], dtype=float)

                expected_shape = {
                    "triangle": "equilateral_triangle",
                    "sector": "circular_sector",
                }.get(shape, shape)
                self.assertEqual(response["shape"].lower().replace(" ", "_"), expected_shape)
                self.assertEqual(wavefunction.shape, (server.GRID_SIZE, server.GRID_SIZE))
                self.assertEqual(density.shape, wavefunction.shape)
                self.assertEqual(mask.shape, wavefunction.shape)
                self.assertTrue(np.isfinite(response["energy"]))
                self.assertGreater(response["energy"], 0)
                self.assertTrue(np.isfinite(response["energy_mev"]))
                self.assertTrue(np.isfinite(wavefunction[mask]).all())
                self.assertTrue(np.isfinite(density[mask]).all())
                self.assertTrue(math.isfinite(response["inputs"].get("m", response["inputs"].get("angular_momentum", 1))))
                dx = np.nanmean(np.abs(np.diff(x, axis=1)))
                dy = np.nanmean(np.abs(np.diff(y, axis=0)))
                probability = np.sum(density[mask]) * dx * dy
                self.assertAlmostEqual(probability, 1.0, delta=0.002)
                json.dumps(server._json_safe(response), allow_nan=False)

    def test_rejects_invalid_shape_inputs(self):
        invalid_requests = [
            {"geometry": "unknown", "dimensions": {}, "quantumNumbers": {"m": 1, "n": 1}},
            {"geometry": "rectangle", "dimensions": {"width": 8}, "quantumNumbers": {"m": 1, "n": 1}},
            {"geometry": "triangle", "dimensions": {"side": 7}, "quantumNumbers": {"m": 1, "n": 2}},
            {"geometry": "annulus", "dimensions": {"innerRadius": 5, "outerRadius": 2}, "quantumNumbers": {"m": 0, "n": 1}},
            {"geometry": "circle", "dimensions": {"radius": 4}, "quantumNumbers": {"m": -1, "n": 1}},
        ]
        for payload in invalid_requests:
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                server.calculate_request(payload)


if __name__ == "__main__":
    unittest.main()