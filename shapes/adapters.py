"""Translate UI inputs and solver results into one SI-based interface."""

from __future__ import annotations

import io
import math
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros

from shapes import circle as circle_module
from shapes import ellipse as ellipse_module
from shapes import sector as sector_module
from shapes.annulus import energy as annulus_energy
from shapes.annulus import run_annulus
from shapes.equilateral_triangle import run_triangle
from shapes.equilateral_triangle import wavefunction as triangle_wavefunction
from shapes.rectangle import energy as rectangle_energy
from shapes.rectangle import psi as rectangle_psi
from shapes.rectangle import run_rectangle


HBAR_SI = 1.054571817e-34
ELECTRON_MASS_SI = 9.1093837015e-31
PLOT_ROOT = Path("generated_plots")


def _positive_number(value, name):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive finite number.") from exc
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValueError(f"{name} must be a positive finite number.")
    return parsed


def _integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    parsed = int(value)
    if parsed < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return parsed


def _grid_size(grid):
    return _integer(grid, "grid", minimum=30)


def _quiet_backend_call(function, *args, **kwargs):
    output = io.StringIO()
    with redirect_stdout(output), redirect_stderr(output):
        return function(*args, **kwargs)


def _result(shape, inputs, energy, wavefunction, probability_density, X, Y,
            domain_mask, states, plots=None, extras=None):
    wave = np.asarray(wavefunction)
    density = np.asarray(probability_density, dtype=float)
    x = np.asarray(X, dtype=float)
    y = np.asarray(Y, dtype=float)
    mask = np.asarray(domain_mask, dtype=bool)
    if wave.shape != density.shape or wave.shape != x.shape or wave.shape != y.shape or wave.shape != mask.shape:
        raise ValueError(f"{shape} solver returned incompatible coordinate and field arrays.")
    if np.iscomplexobj(wave):
        density = np.abs(wave) ** 2
        wave = np.real(wave)
    result = {
        "shape": shape,
        "shape_id": shape.lower().replace(" ", "_"),
        "inputs": inputs,
        "energy": float(energy),
        "energy_units": "J",
        "wavefunction": np.asarray(wave, dtype=float),
        "probability_density": density,
        "X": x,
        "Y": y,
        "coordinates": {"x": x, "y": y, "mask": mask, "units": "m"},
        "states": states,
        "plots": {
            "wavefunction": None,
            "probability_density": None,
            "energy_levels": None,
            **(plots or {}),
        },
    }
    if extras:
        result.update(extras)
    return result


def _state(m, n, energy, **metadata):
    return {"m": int(m), "n": int(n), "E": float(energy), **metadata}


def _select_state(states, m, n, energy, **metadata):
    for state in states:
        state.pop("selected", None)
    selected = next(
        (
            state
            for state in states
            if state["m"] == m
            and state["n"] == n
            and all(state.get(key) == value for key, value in metadata.items())
        ),
        None,
    )
    if selected is None:
        selected = _state(m, n, energy, **metadata)
        states.append(selected)
    selected["selected"] = True
    return states


def calculate_rectangle(width, height, m, n, grid=70):
    """Calculate a rectangle state; geometry and coordinates are in metres."""
    width = _positive_number(width, "width")
    height = _positive_number(height, "height")
    m = _integer(m, "m")
    n = _integer(n, "n")
    grid = _grid_size(grid)
    result = _quiet_backend_call(
        run_rectangle,
        a=width,
        b=height,
        num_states=6,
        grid=grid,
        output_dir=str(PLOT_ROOT / "rectangle"),
    )
    solver_x = np.asarray(result["X"], dtype=float)
    solver_y = np.asarray(result["Y"], dtype=float)
    field = rectangle_psi(width, height, m, n, solver_x, solver_y)
    X = solver_x - width / 2
    Y = solver_y - height / 2
    selected_energy = rectangle_energy(width, height, m, n)
    states = [
        _state(item["m"], item["n"], item["E"])
        for item in result["states"]
    ]
    states = _select_state(states, m, n, selected_energy)
    return _result(
        "Rectangle",
        {"width": width, "height": height, "m": m, "n": n},
        selected_energy,
        field,
        np.abs(field) ** 2,
        X,
        Y,
        np.ones_like(X, dtype=bool),
        states,
        {
            "wavefunction": result.get("psi_plot"),
            "probability_density": result.get("rho_plot"),
            "energy_levels": result.get("energy_plot"),
        },
    )


def calculate_triangle(side, m, n, grid=70):
    """Calculate the supported sine-family equilateral-triangle state."""
    side = _positive_number(side, "side")
    m = _integer(m, "m")
    n = _integer(n, "n")
    if not (1 <= n <= m <= 5):
        raise ValueError("Triangle quantum numbers must satisfy 1 <= n <= m <= 5.")
    grid = _grid_size(grid)
    out = PLOT_ROOT / "triangle"
    result = _quiet_backend_call(
        run_triangle, a=side, m=m, n=n, output_dir=str(out), divisions=grid
    )
    x = np.linspace(-side / 2, side / 2, grid)
    y = np.linspace(0, math.sqrt(3) * side / 2, grid)
    X, Y = np.meshgrid(x, y)
    solver_x = X + side / 2
    mask = (
        (Y >= 0)
        & (Y <= math.sqrt(3) * (X + side / 2))
        & (Y <= math.sqrt(3) * (side / 2 - X))
    )
    field = np.where(
        mask,
        triangle_wavefunction(side, m, n, "S", solver_x, Y),
        np.nan,
    )
    coefficient = m**2 + m * n + n**2

    def energy_for(mi, ni):
        return (
            (mi**2 + mi * ni + ni**2)
            * 8
            * math.pi**2
            * HBAR_SI**2
            / (9 * ELECTRON_MASS_SI * side**2)
        )

    energy = energy_for(m, n)
    states = _select_state(
        [_state(1, 1, energy_for(1, 1))],
        m,
        n,
        energy,
    )
    plots = {
        "wavefunction": str((out / f"wavefunction_m{m}_n{n}_s.png").resolve()),
        "probability_density": str((out / f"probability_density_m{m}_n{n}_s.png").resolve()),
        "energy_levels": str((out / f"energy_level_m{m}_n{n}.png").resolve()),
    }
    return _result(
        "Equilateral Triangle",
        {"side": side, "m": m, "n": n},
        energy,
        field,
        np.abs(field) ** 2,
        X,
        Y,
        mask,
        states,
        plots,
        {"triangle_details": result, "energy_coefficient": coefficient},
    )


def calculate_circle(radius, m, n, grid=70):
    """Calculate one circular-well mode; m is angular and n radial."""
    radius = _positive_number(radius, "radius")
    m = _integer(m, "m", minimum=0)
    n = _integer(n, "n")
    grid = _grid_size(grid)
    result = _quiet_backend_call(
        circle_module.run_circle,
        R=radius,
        n_states=6,
        grid_size=grid,
        mass=1.0,
        hbar=1.0,
        output_dir=str(PLOT_ROOT / "circle"),
        make_plots=False,
    )
    X = np.asarray(result["coordinates"]["x"], dtype=float)
    Y = np.asarray(result["coordinates"]["y"], dtype=float)
    r = np.asarray(result["coordinates"]["r"], dtype=float)
    theta = np.asarray(result["coordinates"]["theta"], dtype=float)
    mask = np.asarray(result["coordinates"]["mask"], dtype=bool)
    zero = float(jn_zeros(m, n)[-1])
    field = circle_module._psi(radius, m, "C", zero, r, theta, mask)
    energy = HBAR_SI**2 * zero**2 / (2 * ELECTRON_MASS_SI * radius**2)
    states = [
        _state(
            item["m"],
            item["n"],
            HBAR_SI**2 * item["bessel_zero"]**2
            / (2 * ELECTRON_MASS_SI * radius**2),
            angular_type=item["angular_type"],
        )
        for item in result["states"]
    ]
    states = _select_state(
        states,
        m,
        n,
        energy,
        angular_type="C",
    )
    return _result(
        "Circle",
        {"radius": radius, "m": m, "n": n},
        energy,
        field,
        np.abs(field) ** 2,
        X,
        Y,
        mask,
        states,
        extras={"angular_type": "C", "bessel_zero": zero},
    )


def calculate_annulus(r_inner, r_outer, angular_momentum, radial_index, grid=70):
    """Calculate an annular state; radii and coordinates are in metres."""
    inner = _positive_number(r_inner, "r_inner")
    outer = _positive_number(r_outer, "r_outer")
    if inner >= outer:
        raise ValueError("Choose radii with 0 < r_inner < r_outer.")
    angular_momentum = _integer(angular_momentum, "angular_momentum", minimum=0)
    radial_index = _integer(radial_index, "radial_index")
    grid = _grid_size(grid)
    result = _quiet_backend_call(
        run_annulus,
        r_inner=inner,
        r_outer=outer,
        angular_momentum=angular_momentum,
        radial_index=radial_index,
        output_dir=str(PLOT_ROOT / "annulus"),
        divisions=grid,
        mass=1.0,
        hbar=1.0,
    )
    X = np.asarray(result["wavefunction"]["x"], dtype=float)
    Y = np.asarray(result["wavefunction"]["y"], dtype=float)
    mask = ~np.asarray(result["plot_grid"]["annulus_mask"], dtype=bool)
    field = np.where(
        mask,
        np.asarray(result["wavefunction"]["values"], dtype=float),
        np.nan,
    )
    density = np.where(
        mask,
        np.asarray(result["probability_density"]["values"], dtype=float),
        np.nan,
    )
    scale = HBAR_SI**2 / ELECTRON_MASS_SI
    energy = float(result["energy"]) * scale
    ground_energy = annulus_energy(inner, outer, 0, 1) * scale
    states = _select_state(
        [_state(0, 1, ground_energy)],
        angular_momentum,
        radial_index,
        energy,
    )
    return _result(
        "Annulus",
        {
            "r_inner": inner,
            "r_outer": outer,
            "angular_momentum": angular_momentum,
            "radial_index": radial_index,
        },
        energy,
        field,
        density,
        X,
        Y,
        mask,
        states,
        {
            "wavefunction": result["plots"]["wavefunction"],
            "probability_density": result["plots"]["probability_density"],
            "energy_levels": result["plots"]["energy_level"],
        },
        {"annulus_details": result},
    )


def _ellipse_mode_state(kind, m, radial_index, minimum_k):
    k_max = max(6.0, minimum_k)
    while k_max <= 1000:
        roots = ellipse_module.roots_for(kind, m, k_max)
        if len(roots) >= radial_index:
            return roots[radial_index - 1]
        k_max *= 1.5
    raise ValueError("Could not find the requested ellipse radial state.")


def calculate_ellipse(semi_major, semi_minor, m, n, grid=70):
    """Calculate an ellipse state; m is angular order and n radial index."""
    semi_major = _positive_number(semi_major, "semi_major")
    semi_minor = _positive_number(semi_minor, "semi_minor")
    if semi_minor >= semi_major:
        raise ValueError("Ellipse geometry requires 0 < semi_minor < semi_major.")
    m = _integer(m, "m", minimum=0)
    n = _integer(n, "n")
    grid = _grid_size(grid)

    names = (
        "SEMI_MAJOR", "SEMI_MINOR", "B_OVER_A", "ECC", "C_FOC", "XI0",
        "XI_NODES", "XI_W", "_cache", "_norm_cache",
    )
    previous = {name: getattr(ellipse_module, name) for name in names}
    try:
        ratio = semi_minor / semi_major
        eccentricity = math.sqrt(1.0 - ratio**2)
        xi0 = math.atanh(ratio)
        gl_x, gl_w = np.polynomial.legendre.leggauss(200)
        ellipse_module.SEMI_MAJOR = semi_major
        ellipse_module.SEMI_MINOR = semi_minor
        ellipse_module.B_OVER_A = ratio
        ellipse_module.ECC = eccentricity
        ellipse_module.C_FOC = eccentricity
        ellipse_module.XI0 = xi0
        ellipse_module.XI_NODES = 0.5 * xi0 * (gl_x + 1)
        ellipse_module.XI_W = 0.5 * xi0 * gl_w
        ellipse_module._cache = {}
        ellipse_module._norm_cache = {}

        kind = "C"
        q_value = _ellipse_mode_state(kind, m, n, n * math.pi + m + 2)
        selected_scaled_energy = 4 * q_value / eccentricity**2
        x = np.linspace(-semi_major, semi_major, grid)
        y = np.linspace(-semi_minor, semi_minor, grid)
        X, Y = np.meshgrid(x, y)
        mask = (X / semi_major) ** 2 + (Y / semi_minor) ** 2 <= 1
        field = np.asarray(
            ellipse_module.psi(m, kind, q_value, X / semi_major, Y / semi_major),
            dtype=float,
        ) / semi_major
        field = np.where(mask, field, np.nan)
        energy_scale = HBAR_SI**2 / (2 * ELECTRON_MASS_SI * semi_major**2)
        selected_energy = selected_scaled_energy * energy_scale

        ground_q = _ellipse_mode_state("C", 0, 1, 6.0)
        ground_scaled = 4 * ground_q / eccentricity**2
        ground_energy = ground_scaled * energy_scale
        states = _select_state(
            [_state(0, 1, ground_energy, kind="C")],
            m,
            n,
            selected_energy,
            kind=kind,
        )
        return _result(
            "Ellipse",
            {"semi_major": semi_major, "semi_minor": semi_minor, "m": m, "n": n},
            selected_energy,
            field,
            np.abs(field) ** 2,
            X,
            Y,
            mask,
            states,
            extras={
                "ellipse_state": {
                    "kind": kind,
                    "q": q_value,
                    "energy_scaled": selected_scaled_energy,
                }
            },
        )
    finally:
        for name, value in previous.items():
            setattr(ellipse_module, name, value)


def calculate_sector(radius, angle, m, n, start_angle=0.0, grid=70):
    """Calculate a sector state; m=ell and n=s, both one-based."""
    radius = _positive_number(radius, "radius")
    try:
        angle_degrees = float(angle)
        start_degrees = float(start_angle)
    except (TypeError, ValueError) as exc:
        raise ValueError("Sector angles must be finite numbers in degrees.") from exc
    if not math.isfinite(angle_degrees) or not 0 < angle_degrees <= 360:
        raise ValueError("angle must be in degrees and satisfy 0 < angle <= 360.")
    if not math.isfinite(start_degrees):
        raise ValueError("start_angle must be finite.")
    m = _integer(m, "m")
    n = _integer(n, "n")
    grid = _grid_size(grid)

    previous = {
        "RADIUS": sector_module.RADIUS,
        "HBAR": sector_module.HBAR,
        "MASS": sector_module.MASS,
    }
    try:
        sector_module.RADIUS = radius
        sector_module.HBAR = 1.0
        sector_module.MASS = 1.0
        alpha = math.radians(angle_degrees)
        nu = m * math.pi / alpha
        root = float(sector_module.bessel_positive_zeros(nu, n)[-1])
        selected_state = {
            "ell": m,
            "s": n,
            "nu": nu,
            "zero": root,
            "energy": root**2 / (2 * radius**2),
        }
        x = np.linspace(-radius, radius, grid)
        y = np.linspace(-radius, radius, grid)
        X, Y = np.meshgrid(x, y)
        r = np.hypot(X, Y)
        theta = (np.arctan2(Y, X) - math.radians(start_degrees)) % (2 * math.pi)
        mask = (r <= radius) & (theta <= alpha)
        field = np.asarray(
            sector_module.wavefunction(selected_state, r, theta, alpha), dtype=float
        )
        field = np.where(mask, field, np.nan)
        energy = selected_state["energy"] * HBAR_SI**2 / ELECTRON_MASS_SI
        ground_nu = math.pi / alpha
        ground_root = float(sector_module.bessel_positive_zeros(ground_nu, 1)[0])
        ground_energy = (
            ground_root**2 / (2 * radius**2) * HBAR_SI**2 / ELECTRON_MASS_SI
        )
        states = _select_state(
            [_state(1, 1, ground_energy)],
            m,
            n,
            energy,
        )
        return _result(
            "Circular Sector",
            {
                "radius": radius,
                "angle_deg": angle_degrees,
                "start_angle_deg": start_degrees,
                "m": m,
                "n": n,
            },
            energy,
            field,
            np.abs(field) ** 2,
            X,
            Y,
            mask,
            states,
            extras={"sector_state": selected_state},
        )
    finally:
        for name, value in previous.items():
            setattr(sector_module, name, value)


def calculate_shape(shape, **kwargs):
    """Dispatch a UI shape label to exactly one shape adapter."""
    name = str(shape).strip().lower()
    adapters = {
        "rectangle": calculate_rectangle,
        "circle": calculate_circle,
        "equilateral triangle": calculate_triangle,
        "triangle": calculate_triangle,
        "annulus": calculate_annulus,
        "ellipse": calculate_ellipse,
        "circular sector": calculate_sector,
        "sector": calculate_sector,
    }
    if name not in adapters:
        raise ValueError(f"Unsupported shape '{shape}'.")
    return adapters[name](**kwargs)


__all__ = [
    "calculate_rectangle",
    "calculate_triangle",
    "calculate_circle",
    "calculate_annulus",
    "calculate_ellipse",
    "calculate_sector",
    "calculate_shape",
]
