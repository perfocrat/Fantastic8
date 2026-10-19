"""Backend for a 2D particle in a circular infinite potential well.

UI entry point: run_circle(R=1.0, n_states=6, grid_size=120)
Dependencies: numpy, scipy, matplotlib
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros

NUM_STATES = 6
GRID = 120
DEFAULT_MASS = 1.0
DEFAULT_HBAR = 1.0
PLOT_DIR = Path("generated_plots")


def _validate(R, n_states, grid_size, mass, hbar):
    if not np.isfinite(R) or R <= 0:
        raise ValueError("R must be a positive finite number.")
    if not isinstance(n_states, (int, np.integer)) or n_states < 1:
        raise ValueError("n_states must be a positive integer.")
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 30:
        raise ValueError("grid_size must be an integer >= 30.")
    if not np.isfinite(mass) or mass <= 0:
        raise ValueError("mass must be positive and finite.")
    if not np.isfinite(hbar) or hbar <= 0:
        raise ValueError("hbar must be positive and finite.")


def _lowest_states(count):
    """States sorted by energy. m=0 has one angular state; m>0 has C/S pair."""
    candidates = []
    max_m = max(8, count + 3)
    for m in range(max_m + 1):
        for n, alpha in enumerate(jn_zeros(m, count), start=1):
            candidates.append((float(alpha**2), m, n, "C", float(alpha)))
            if m > 0:
                candidates.append((float(alpha**2), m, n, "S", float(alpha)))
    candidates.sort(key=lambda s: (s[0], s[1], s[2], s[3]))
    return candidates[:count]


def _grid(R, size):
    coords = np.linspace(-R, R, size)
    X, Y = np.meshgrid(coords, coords)
    r = np.hypot(X, Y)
    theta = np.arctan2(Y, X)
    mask = r <= R
    return X, Y, r, theta, mask


def _normalization(R, m, alpha):
    """Analytical normalization over area element r dr dtheta."""
    angular_integral = 2.0 * np.pi if m == 0 else np.pi
    radial_integral = 0.5 * R**2 * jv(m + 1, alpha)**2
    return 1.0 / np.sqrt(angular_integral * radial_integral)


def _psi(R, m, angular_type, alpha, r, theta, mask):
    N = _normalization(R, m, alpha)
    radial = jv(m, alpha * r / R)
    if m == 0:
        angular = np.ones_like(theta)
    elif angular_type == "C":
        angular = np.cos(m * theta)
    else:
        angular = np.sin(m * theta)
    return np.where(mask, N * radial * angular, np.nan)


def _plot_shape(R, out):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.add_patch(plt.Circle((0, 0), R, fill=False, linewidth=2))
    ax.annotate("R", xy=(R / np.sqrt(2), R / np.sqrt(2)),
                xytext=(0.25 * R, 0.75 * R), arrowprops={"arrowstyle": "->"})
    ax.set(xlim=(-1.15 * R, 1.15 * R), ylim=(-1.15 * R, 1.15 * R),
           xlabel="x", ylabel="y", title="2D circular infinite potential well")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = out / "circle_shape.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    return fig, str(path)


def _plot_energy(states, R, mass, hbar, out):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    energies = []
    for i, (_, m, n, kind, alpha) in enumerate(states, start=1):
        E = hbar**2 * alpha**2 / (2.0 * mass * R**2)
        energies.append(float(E))
        ax.hlines(E, 0.15, 0.75, linewidth=2)
        label = f"State {i}: (m,n)=({m},{n})"
        if m > 0:
            label += f", {kind}"
        ax.text(0.82, E, f"{label}\nE = {E:.5g}", va="center", fontsize=8)
    ax.set_xlim(0, 2.8)
    ax.set_xticks([])
    ax.set_ylabel("Energy E")
    ax.set_title(r"Circular well: $E=\hbar^2\alpha_{mn}^2/(2MR^2)$")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    path = out / "circle_energy_levels.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    return fig, str(path), energies


def _plot_states(states, X, Y, wavefunctions, densities, R, out):
    count = len(states)
    cols = 3
    rows = int(np.ceil(count / cols))
    fig_psi, axes_psi = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.8 * rows), squeeze=False)
    fig_rho, axes_rho = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.8 * rows), squeeze=False)
    for ax in axes_psi.ravel():
        ax.axis("off")
    for ax in axes_rho.ravel():
        ax.axis("off")

    for i, (state, psi, rho) in enumerate(zip(states, wavefunctions, densities)):
        _, m, n, kind, _ = state
        axp, axr = axes_psi.ravel()[i], axes_rho.ravel()[i]
        axp.axis("on")
        axr.axis("on")
        vmax = np.nanmax(np.abs(psi))
        im1 = axp.pcolormesh(X, Y, psi, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        im2 = axr.pcolormesh(X, Y, rho, shading="auto", cmap="viridis")
        axp.add_patch(plt.Circle((0, 0), R, fill=False, color="black", linewidth=0.8))
        axr.add_patch(plt.Circle((0, 0), R, fill=False, color="black", linewidth=0.8))
        axp.set_title(f"State {i+1}: (m,n)=({m},{n}) {kind}")
        axr.set_title(f"State {i+1}: (m,n)=({m},{n}) {kind}")
        for ax in (axp, axr):
            ax.set_aspect("equal")
            ax.set(xlim=(-R, R), ylim=(-R, R), xlabel="x", ylabel="y")
        fig_psi.colorbar(im1, ax=axp, fraction=0.046, pad=0.04)
        fig_rho.colorbar(im2, ax=axr, fraction=0.046, pad=0.04)
    fig_psi.suptitle("Normalized wavefunctions")
    fig_rho.suptitle(r"Probability density $|\psi(x,y)|^2$")
    fig_psi.tight_layout()
    fig_rho.tight_layout()
    psi_path = out / "circle_wavefunctions.png"
    rho_path = out / "circle_probability_density.png"
    fig_psi.savefig(psi_path, dpi=160, bbox_inches="tight")
    fig_rho.savefig(rho_path, dpi=160, bbox_inches="tight")
    return fig_psi, fig_rho, str(psi_path), str(rho_path)


def run_circle(R=1.0, n_states=NUM_STATES, grid_size=GRID,
               mass=DEFAULT_MASS, hbar=DEFAULT_HBAR,
               output_dir=PLOT_DIR, make_plots=True, show_plots=False):
    """Compute states and optionally save plots. UI should call this function.

    R: circle radius (the UI slider parameter)
    n_states: number of lowest states to return, counting degenerate C/S states
    grid_size: Cartesian plotting grid size
    mass, hbar: physical constants in a consistent unit system
    output_dir: folder for generated PNGs
    make_plots: set False for numerical-only calculations
    show_plots: display figures (normally False when called by a UI)
    """
    _validate(R, n_states, grid_size, mass, hbar)
    out = Path(output_dir)
    if make_plots:
        out.mkdir(parents=True, exist_ok=True)

    states = _lowest_states(n_states)
    X, Y, r, theta, mask = _grid(R, grid_size)
    wavefunctions, densities, energies, state_results = [], [], [], []

    for index, (_, m, n, kind, alpha) in enumerate(states, start=1):
        psi = _psi(R, m, kind, alpha, r, theta, mask)
        rho = np.where(mask, np.abs(psi)**2, np.nan)
        energy = float(hbar**2 * alpha**2 / (2.0 * mass * R**2))
        wavefunctions.append(psi)
        densities.append(rho)
        energies.append(energy)
        state_results.append({
            "state": index, "m": int(m), "n": int(n),
            "angular_type": kind, "bessel_zero": float(alpha),
            "energy": energy, "energy_factor": float(alpha**2),
            "wavefunction": psi, "probability_density": rho,
        })

    plots, figures = {}, []
    if make_plots:
        fig, path = _plot_shape(R, out)
        figures.append(fig)
        plots["shape"] = path
        fig, path, _ = _plot_energy(states, R, mass, hbar, out)
        figures.append(fig)
        plots["energy_levels"] = path
        figp, figr, ppath, rpath = _plot_states(
            states, X, Y, wavefunctions, densities, R, out
        )
        figures.extend([figp, figr])
        plots["wavefunctions"] = ppath
        plots["probability_density"] = rpath

    result = {
        "shape": "circle",
        "model": "2D circular infinite potential well",
        "parameters": {"R": float(R), "n_states": int(n_states),
                       "grid_size": int(grid_size), "mass": float(mass), "hbar": float(hbar)},
        "coordinates": {"x": X, "y": Y, "r": r, "theta": theta, "mask": mask},
        "energies": energies,
        "wavefunctions": wavefunctions,
        "probability_densities": densities,
        "states": state_results,
        "plots": plots,
        "figures": figures,
    }
    if show_plots and make_plots:
        plt.show()
    return result


if __name__ == "__main__":
    # Local test only; the UI should call run_circle(...) directly.
    result = run_circle(R=1.0, n_states=NUM_STATES, grid_size=GRID,
                        make_plots=True, show_plots=True)
    print("Shape:", result["shape"])
    print("Parameters:", result["parameters"])
    print("\nEnergy levels:")
    for state in result["states"]:
        print(f"State {state['state']}: (m,n)=({state['m']},{state['n']}), "
              f"type={state['angular_type']}, E={state['energy']:.8g}")
    print("\nPlots saved to:")
    for name, path in result["plots"].items():
        print(f"  {name}: {path}")
