"""Visualize a particle in a two-dimensional, infinitely deep annular well.

The particle is confined to r_inner < r < r_outer. The wavefunction is
separated into angular modes exp(i*l*theta) and radial Bessel solutions.
Dirichlet boundary conditions (psi=0) are imposed at both circular walls.

Requires numpy, scipy, and matplotlib. Run this file directly for an example.
"""
from __future__ import annotations

from pathlib import Path
from functools import lru_cache
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, yv


DEFAULT_GRID = 240


def validate_inputs(r_inner: float, r_outer: float, angular_momentum: int,
                    radial_index: int) -> tuple[float, float, int, int]:
    """Validate radii and mode indices (both indices start at 1/0)."""
    try:
        inner, outer = float(r_inner), float(r_outer)
    except (TypeError, ValueError) as exc:
        raise ValueError("Radii must be positive finite numbers.") from exc
    if not math.isfinite(inner) or not math.isfinite(outer) or inner <= 0 or outer <= inner:
        raise ValueError("Choose finite radii with 0 < r_inner < r_outer.")
    for name, value, minimum in (("angular_momentum", angular_momentum, 0),
                                 ("radial_index", radial_index, 1)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}.")
    return inner, outer, int(angular_momentum), int(radial_index)


def _radial_boundary_function(order: int, inner: float, outer: float, k: float) -> float:
    """Determinant for a Bessel combination that vanishes at both walls."""
    return float(jv(order, k * inner) * yv(order, k * outer)
                 - yv(order, k * inner) * jv(order, k * outer))


@lru_cache(maxsize=256)
def _cached_radial_wavenumber(inner: float, outer: float, order: int, index: int) -> float:
    """Cached numerical root search for repeated wavefunction and energy calls."""
    width = outer - inner
    k_max = (index + order + 12) * math.pi / width
    samples = np.linspace(1e-8 / outer, k_max, max(12000, 2000 * (index + order + 1)))
    vals = np.array([_radial_boundary_function(order, inner, outer, k) for k in samples])
    brackets = np.flatnonzero(vals[:-1] * vals[1:] < 0)
    if len(brackets) < index:
        raise RuntimeError("Could not bracket the requested radial eigenvalue.")
    lo, hi = samples[brackets[index - 1]], samples[brackets[index - 1] + 1]
    function = lambda k: _radial_boundary_function(order, inner, outer, k)
    return float(brentq(function, lo, hi, xtol=1e-12, rtol=1e-12))


def radial_wavenumber(r_inner: float, r_outer: float, angular_momentum: int,
                      radial_index: int) -> float:
    """Find the requested positive radial eigenvalue k in inverse length units."""
    inner, outer, order, index = validate_inputs(
        r_inner, r_outer, angular_momentum, radial_index
    )
    # Scan densely, bracket each sign-changing root, and cache the result.
    return _cached_radial_wavenumber(inner, outer, order, index)


def _radial_shape(r: np.ndarray, inner: float, k: float, order: int) -> np.ndarray:
    """Unnormalized radial solution, zero at r_inner and at the outer wall."""
    return (jv(order, k * r) * yv(order, k * inner)
            - yv(order, k * r) * jv(order, k * inner))


def _normalization_coefficient(inner: float, outer: float, k: float,
                               order: int, component: str) -> float:
    """Return the radial coefficient for a unit-normalized 2D eigenstate."""
    # Gauss-Legendre quadrature accurately integrates the oscillatory radial
    # factor, including the polar Jacobian r. The angular factor is normalized
    # separately: integral cos(l theta)^2 dtheta = pi for l>0, or 2*pi for l=0.
    nodes, weights = np.polynomial.legendre.leggauss(256)
    radii = (outer - inner) * (nodes + 1.0) / 2.0 + inner
    radial_weights = weights * (outer - inner) / 2.0
    raw = _radial_shape(radii, inner, k, order)
    radial_integral = float(np.sum(np.abs(raw) ** 2 * radii * radial_weights))
    angular_integral = 2.0 * math.pi if order == 0 or component == "complex" else math.pi
    return 1.0 / math.sqrt(radial_integral * angular_integral)


def wavefunction(r_inner: float, r_outer: float, angular_momentum: int,
                 radial_index: int, r: np.ndarray, theta: np.ndarray,
                 component: str = "real") -> np.ndarray:
    """Evaluate a normalized eigenstate on polar-coordinate arrays.

    component can be 'real', 'imag', or 'complex'. A nonzero angular mode has
    two degenerate real states proportional to cos(l*theta) and sin(l*theta);
    'real' selects the cosine member, while 'imag' selects the sine member.
    """
    inner, outer, order, index = validate_inputs(
        r_inner, r_outer, angular_momentum, radial_index
    )
    component = component.lower()
    if component not in ("real", "imag", "complex"):
        raise ValueError("component must be 'real', 'imag', or 'complex'.")
    if order == 0 and component == "imag":
        raise ValueError("The l=0 state has no sine angular member.")
    k = radial_wavenumber(inner, outer, order, index)
    radial = _radial_shape(np.asarray(r, dtype=float), inner, k, order)
    coefficient = _normalization_coefficient(inner, outer, k, order, component)
    if order == 0:
        angular = np.ones_like(np.asarray(theta, dtype=float))
    elif component == "complex":
        angular = np.exp(1j * order * np.asarray(theta, dtype=float))
    else:
        trig = np.sin if component == "imag" else np.cos
        angular = trig(order * np.asarray(theta, dtype=float))
    return coefficient * radial * angular


def energy(r_inner: float, r_outer: float, angular_momentum: int,
           radial_index: int, mass: float = 1.0, hbar: float = 1.0) -> float:
    """Return E = hbar^2 k^2/(2 mass), in the units supplied by the caller."""
    inner, outer, order, index = validate_inputs(
        r_inner, r_outer, angular_momentum, radial_index
    )
    if not math.isfinite(mass) or mass <= 0 or not math.isfinite(hbar) or hbar <= 0:
        raise ValueError("mass and hbar must be positive finite numbers.")
    k = radial_wavenumber(inner, outer, order, index)
    return hbar * hbar * k * k / (2.0 * mass)


def _polar_grid(r_inner: float, r_outer: float, divisions: int):
    if divisions < 20:
        raise ValueError("divisions must be at least 20.")
    x = np.linspace(-r_outer, r_outer, divisions)
    y = np.linspace(-r_outer, r_outer, divisions)
    xx, yy = np.meshgrid(x, y)
    rr = np.hypot(xx, yy)
    theta = np.arctan2(yy, xx)
    mask = (rr < r_inner) | (rr > r_outer)
    return xx, yy, rr, theta, mask


def run_annulus(r_inner: float, r_outer: float, angular_momentum: int,
                radial_index: int, output_dir: str | Path = "generated_plots",
                divisions: int = DEFAULT_GRID, mass: float = 1.0,
                hbar: float = 1.0) -> dict:
    """Create wavefunction, probability, and energy plots and return metadata."""
    inner, outer, order, index = validate_inputs(
        r_inner, r_outer, angular_momentum, radial_index
    )
    xx, yy, rr, theta, mask = _polar_grid(inner, outer, divisions)
    psi = wavefunction(inner, outer, order, index, rr, theta, "real")
    probability_density = np.abs(psi) ** 2
    psi_plot = np.ma.array(np.real(psi), mask=mask)
    prob_plot = np.ma.array(probability_density, mask=mask)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stem = f"l{order}_q{index}"

    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    limit = max(float(np.max(np.abs(psi_plot.compressed()))), 1e-12)
    artist = ax.pcolormesh(xx, yy, psi_plot, shading="auto", cmap="RdBu_r",
                           vmin=-limit, vmax=limit)
    fig.colorbar(artist, ax=ax, label=r"$\psi(x,y)$ (real part)")
    ax.set_title(f"Annular-well wavefunction — l={order}, q={index}")
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout()
    psi_path = output / f"wavefunction_{stem}_cos.png"
    fig.savefig(psi_path, dpi=160, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    artist = ax.pcolormesh(xx, yy, prob_plot, shading="auto", cmap="viridis")
    fig.colorbar(artist, ax=ax, label=r"$|\psi(x,y)|^2$")
    ax.set_title(f"Probability density — l={order}, q={index}")
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout()
    prob_path = output / f"probability_density_{stem}_cos.png"
    fig.savefig(prob_path, dpi=160, bbox_inches="tight"); plt.close(fig)

    k = radial_wavenumber(inner, outer, order, index)
    e = energy(inner, outer, order, index, mass, hbar)
    normalization_constant = _normalization_coefficient(inner, outer, k, order, "real")
    radial_equation = (
        f"{normalization_constant:.12g} * [J_{order}({k:.12g}*r) "
        f"* Y_{order}({k*inner:.12g}) - Y_{order}({k:.12g}*r) "
        f"* J_{order}({k*inner:.12g})]"
    )
    angular_equation = "1" if order == 0 else f"cos({order}*theta)"
    normalized_equation = f"psi(r,theta) = ({radial_equation}) * {angular_equation}"
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    level = e
    ax.hlines(level, 0.15, 0.75, linewidth=3)
    ax.plot([0.45], [level], marker="o", markersize=6)
    ax.text(0.82, level, f"Selected state (l,q)=({order},{index})\nE = {e:.6g}",
            va="center", ha="left")
    ax.set_xlim(0, 1.8); ax.set_ylim(0, max(1.0, level * 1.35))
    ax.set_xticks([]); ax.set_ylabel("Energy")
    ax.set_title(f"Selected annular-well energy — radii {inner:g} to {outer:g}")
    ax.grid(axis="y", alpha=0.3); fig.tight_layout()
    energy_path = output / f"energy_level_{stem}.png"
    fig.savefig(energy_path, dpi=160, bbox_inches="tight"); plt.close(fig)

    # Return plain nested lists so the result can be serialized to JSON and
    # passed directly through a UI/API boundary without NumPy-specific types.
    wavefunction_data = np.where(mask, 0.0, np.real(psi)).tolist()
    probability_data = np.where(mask, 0.0, probability_density).tolist()
    coordinate_x = xx.tolist()
    coordinate_y = yy.tolist()
    energy_levels = [{
        "angular_momentum_l": order,
        "radial_index_q": index,
        "radial_wavenumber_k": k,
        "energy": e,
        "degeneracy": 1 if order == 0 else 2,
    }]
    return {
        "shape": "annulus",
        "inputs": {"r_inner": inner, "r_outer": outer,
                   "angular_momentum_l": order, "radial_index_q": index,
                   "mass": mass, "hbar": hbar},
        "boundary_condition": "psi(r_inner)=psi(r_outer)=0",
        "radial_solution": "A*[J_l(k*r)*Y_l(k*r_inner)-Y_l(k*r)*J_l(k*r_inner)]",
        "normalization_constant_A": normalization_constant,
        "normalized_wavefunction_equation": normalized_equation,
        "radial_wavenumber_k": k,
        "energy_expression": "E = hbar^2*k^2/(2*mass)",
        "energy": e,
        "energy_levels": energy_levels,
        "angular_degeneracy": 1 if order == 0 else 2,
        "normalization": "integral_annulus |psi|^2 dA = 1",
        "wavefunction": {
            "representation": "real cosine angular state",
            "x": coordinate_x,
            "y": coordinate_y,
            "values": wavefunction_data,
        },
        "probability_density": {
            "definition": "|psi|^2",
            "x": coordinate_x,
            "y": coordinate_y,
            "values": probability_data,
        },
        "plot_grid": {
            "x": xx[0, :].tolist(),
            "y": yy[:, 0].tolist(),
            "annulus_mask": mask.tolist(),
            "outside_region_value": 0.0,
        },
        "plots": {"wavefunction": str(psi_path.resolve()),
                  "probability_density": str(prob_path.resolve()),
                  "energy_level": str(energy_path.resolve())},
    }


def boundary_residual(r_inner: float, r_outer: float, angular_momentum: int,
                      radial_index: int, samples: int = 16) -> float:
    """Return the largest absolute wavefunction value sampled on either wall."""
    inner, outer, order, index = validate_inputs(
        r_inner, r_outer, angular_momentum, radial_index
    )
    if isinstance(samples, bool) or not isinstance(samples, (int, np.integer)) or samples < 2:
        raise ValueError("samples must be an integer of at least 2.")
    theta = np.linspace(0.0, 2.0 * math.pi, int(samples), endpoint=False)
    inner_values = wavefunction(inner, outer, order, index,
                                np.full_like(theta, inner), theta)
    outer_values = wavefunction(inner, outer, order, index,
                                np.full_like(theta, outer), theta)
    return float(max(np.max(np.abs(inner_values)), np.max(np.abs(outer_values))))


if __name__ == "__main__":
    result = run_annulus(r_inner=0.35, r_outer=1.0,
                         angular_momentum=1, radial_index=1)
    for key, value in result.items():
        print(f"{key}: {value}")
