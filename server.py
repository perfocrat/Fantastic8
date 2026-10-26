"""Serve the HTML application and its local physics calculation API."""

from __future__ import annotations

import argparse
import json
import math
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import Any

import numpy as np

from config import FRONTEND_DIR, GRID_SIZE, SERVER_HOST, SERVER_PORT
from shapes.adapters import calculate_shape


FRONTEND_ROOT = FRONTEND_DIR
MEV_IN_JOULES = 1.602176634e-22
MAX_REQUEST_BYTES = 16_384

SHAPE_INPUTS = {
    "rectangle": ("width", "height"),
    "triangle": ("side",),
    "circle": ("radius",),
    "sector": ("radius", "angle", "startAngle"),
    "annulus": ("innerRadius", "outerRadius"),
    "ellipse": ("semiMajor", "semiMinor"),
}


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number.")
    return number


def _quantum_number(value: Any, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return value


def calculate_request(payload: Any) -> dict[str, Any]:
    """Validate a browser request, dispatch one solver, and prepare its result."""
    if not isinstance(payload, dict):
        raise ValueError("Request body must be a JSON object.")

    shape = payload.get("geometry")
    if not isinstance(shape, str) or shape not in SHAPE_INPUTS:
        raise ValueError("Choose one of the available well geometries.")
    dimensions = payload.get("dimensions")
    quantum_numbers = payload.get("quantumNumbers")
    if not isinstance(dimensions, dict) or not isinstance(quantum_numbers, dict):
        raise ValueError("Dimensions and quantumNumbers must be objects.")

    required = SHAPE_INPUTS[shape]
    if set(dimensions) != set(required):
        raise ValueError(f"{shape} requires these dimensions: {', '.join(required)}.")

    minimum_m = 0 if shape in {"circle", "annulus", "ellipse"} else 1
    m = _quantum_number(quantum_numbers.get("m"), "m", minimum_m)
    n = _quantum_number(quantum_numbers.get("n"), "n", 1)

    solver_inputs: dict[str, Any] = {"m": m, "n": n, "grid": GRID_SIZE}
    for key in required:
        value = _finite_number(dimensions.get(key), key)
        if key in {"angle", "startAngle"}:
            solver_inputs["start_angle" if key == "startAngle" else key] = value
        else:
            solver_key = {
                "innerRadius": "r_inner",
                "outerRadius": "r_outer",
                "semiMajor": "semi_major",
                "semiMinor": "semi_minor",
            }.get(key, key)
            solver_inputs[solver_key] = value * 1e-9

    if shape == "annulus":
        solver_inputs["angular_momentum"] = solver_inputs.pop("m")
        solver_inputs["radial_index"] = solver_inputs.pop("n")

    result = calculate_shape(shape, **solver_inputs)
    coordinates = result["coordinates"]
    mask = np.asarray(coordinates["mask"], dtype=bool)
    states = []
    for state in result["states"]:
        state_data = dict(state)
        state_data["energy_mev"] = float(state["E"]) / MEV_IN_JOULES
        states.append(state_data)

    return {
        "shape": result["shape"],
        "inputs": result["inputs"],
        "energy": result["energy"],
        "energy_units": "J",
        "energy_mev": result["energy"] / MEV_IN_JOULES,
        "states": states,
        "wavefunction": result["wavefunction"],
        "probability_density": result["probability_density"],
        "coordinates": coordinates,
        "potential": {
            "values": np.where(mask, 0.0, np.nan),
            "units": "J",
            "outside_domain": "infinite",
        },
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


class RequestHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, directory=str(FRONTEND_ROOT), **kwargs)

    def _send_json(self, status: int, body: dict[str, Any]) -> None:
        encoded = json.dumps(_json_safe(body), allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(encoded)

    def do_GET(self) -> None:
        if self.path == "/api/health":
            self._send_json(200, {"status": "ok"})
            return
        super().do_GET()

    def do_POST(self) -> None:
        if self.path != "/api/calculate":
            self._send_json(404, {"error": "Not found."})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > MAX_REQUEST_BYTES:
                raise ValueError("Request body is empty or too large.")
            payload = json.loads(self.rfile.read(length))
            result = calculate_request(payload)
        except (json.JSONDecodeError, UnicodeDecodeError, ValueError) as exc:
            self._send_json(400, {"error": str(exc)})
            return
        except Exception:
            self.log_error("Physics calculation failed")
            self._send_json(500, {"error": "The physics calculation failed."})
            return
        self._send_json(200, result)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Quantum Well Explorer.")
    parser.add_argument("--host", default=SERVER_HOST)
    parser.add_argument("--port", type=int, default=SERVER_PORT)
    args = parser.parse_args()
    server = HTTPServer((args.host, args.port), RequestHandler)
    print(f"Quantum Well Explorer running at http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()