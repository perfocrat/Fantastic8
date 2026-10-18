import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from math import pi, sqrt


# ============================================================
# CONSTANTS
# ============================================================

HBAR = 1.054571817e-34       # J s
MASS = 9.1093837015e-31      # electron mass


# ============================================================
# 1. SCHRODINGER EQUATION
# ============================================================

def solve_schrodinger():

    print("\n" + "=" * 70)
    print("RECTANGULAR INFINITE POTENTIAL WELL")
    print("=" * 70)

    print("\nPotential:")
    print("V(x,y) = 0 inside the rectangle")
    print("V(x,y) = infinity outside the rectangle")

    print("\nTime-independent Schrodinger equation:")
    print("-hbar^2/(2M) * (d²psi/dx² + d²psi/dy²) = E psi")

    print("\nBoundary conditions:")
    print("psi(0,y) = psi(a,y) = 0")
    print("psi(x,0) = psi(x,b) = 0")

    print("\nNormalized wavefunction:")
    print("psi_mn(x,y) = sqrt(4/(ab))")
    print("              * sin(m*pi*x/a) * sin(n*pi*y/b)")

    print("\nEnergy:")
    print("E_mn = (pi² hbar² / 2M)")
    print("       * (m²/a² + n²/b²)")


# ============================================================
# 2. ENERGY
# ============================================================

def energy(a, b, m, n):

    return (
        (pi ** 2 * HBAR ** 2) / (2 * MASS)
        * (
            (m ** 2 / a ** 2)
            + (n ** 2 / b ** 2)
        )
    )


# ============================================================
# 3. LOWEST ENERGY STATES
# ============================================================

def lowest_states(a, b, count=6):

    states = []

    for m in range(1, 21):
        for n in range(1, 21):

            E = energy(a, b, m, n)

            states.append({
                "m": m,
                "n": n,
                "E": E
            })

    states.sort(key=lambda state: state["E"])

    return states[:count]


# ============================================================
# 4. NORMALIZED WAVEFUNCTION
# ============================================================

def psi(a, b, m, n, x, y):

    normalization = sqrt(4 / (a * b))

    return (
        normalization
        * np.sin(m * pi * x / a)
        * np.sin(n * pi * y / b)
    )


# ============================================================
# 5. PROBABILITY DENSITY
# ============================================================

def probability_density(a, b, m, n, x, y):

    wave = psi(a, b, m, n, x, y)

    return np.abs(wave) ** 2


# ============================================================
# 6. RECTANGLE MESH
# ============================================================

def rectangle_mesh(a, b, N=80):

    x = np.linspace(0, a, N)
    y = np.linspace(0, b, N)

    X, Y = np.meshgrid(x, y)

    return X, Y


# ============================================================
# 7. NORMALIZATION VERIFICATION
# ============================================================

def verify(a, b, states, N=200):

    x = np.linspace(0, a, N)
    y = np.linspace(0, b, N)

    X, Y = np.meshgrid(x, y)

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    print("\nNormalization check:")
    print("-" * 50)

    for state in states:

        m = state["m"]
        n = state["n"]

        wave = psi(a, b, m, n, X, Y)

        probability = np.abs(wave) ** 2

        integral = np.sum(probability) * dx * dy

        print(
            f"(m,n)=({m},{n})    "
            f"integral |psi|² dxdy = {integral:.6f}"
        )


# ============================================================
# 8. ENERGY LEVEL PLOT
# ============================================================

def plot_energy_levels(a, b, states, output_dir):

    fig, ax = plt.subplots(figsize=(8, 6))

    energies = [state["E"] for state in states]

    # Ground-state energy for normalization
    E0 = energies[0]

    for i, state in enumerate(states):

        normalized_energy = state["E"] / E0

        ax.hlines(
            normalized_energy,
            i - 0.35,
            i + 0.35
        )

        ax.text(
            i,
            normalized_energy,
            f"({state['m']},{state['n']})",
            ha="center",
            va="bottom"
        )

    ax.set_xlabel("State")
    ax.set_ylabel("E / E₁")

    ax.set_title(
        "Energy Levels of Rectangular Quantum Box\n"
        f"a = {a:.3e} m, b = {b:.3e} m"
    )

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    path = output_dir / "rectangle_levels.png"

    fig.savefig(path, dpi=200)

    plt.close(fig)

    return path


# ============================================================
# 9. CALCULATE FIELDS
# ============================================================

def fields_on_mesh(a, b, states, X, Y):

    fields = []

    for state in states:

        m = state["m"]
        n = state["n"]

        wave = psi(a, b, m, n, X, Y)

        probability = np.abs(wave) ** 2

        fields.append({
            "m": m,
            "n": n,
            "psi": wave,
            "rho": probability
        })

    return fields


# ============================================================
# 10. WAVEFUNCTION / PROBABILITY PLOTS
# ============================================================

def plot_grid(
    a,
    b,
    states,
    X,
    Y,
    fields,
    what,
    output_dir
):

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(12, 7)
    )

    axes = axes.ravel()

    for i, field in enumerate(fields):

        if what == "psi":

            data = field["psi"]
            label = "ψ"

        else:

            data = field["rho"]
            label = "|ψ|²"

        ax = axes[i]

        contour = ax.contourf(
            X,
            Y,
            data,
            levels=30
        )

        ax.set_aspect("equal")

        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")

        ax.set_title(
            f"{label}, "
            f"(m,n)=({field['m']},{field['n']})"
        )

        fig.colorbar(
            contour,
            ax=ax
        )

    fig.suptitle(
        f"Rectangular Quantum Box - {label}\n"
        f"a = {a:.3e} m, b = {b:.3e} m"
    )

    plt.tight_layout()

    if what == "psi":

        filename = "rectangle_psi.png"

    else:

        filename = "rectangle_rho.png"

    path = output_dir / filename

    fig.savefig(
        path,
        dpi=200
    )

    plt.close(fig)

    return path


# ============================================================
# 11. PRINT ENERGY TABLE
# ============================================================

def print_table(states):

    print("\n" + "=" * 60)
    print("LOWEST ENERGY STATES")
    print("=" * 60)

    print(
        f"{'State':<8}"
        f"{'m':<6}"
        f"{'n':<6}"
        f"{'Energy (J)':<20}"
    )

    print("-" * 60)

    for i, state in enumerate(states, start=1):

        print(
            f"{i:<8}"
            f"{state['m']:<6}"
            f"{state['n']:<6}"
            f"{state['E']:<20.6e}"
        )


# ============================================================
# 12. MAIN FUNCTION FOR FRONTEND
# ============================================================

def run_rectangle(
    a,
    b,
    num_states=6,
    grid=80,
    output_dir="generated_plots"
):

    # --------------------------------------------------------
    # Validate slider values
    # --------------------------------------------------------

    if a <= 0:
        raise ValueError("a must be greater than zero.")

    if b <= 0:
        raise ValueError("b must be greater than zero.")

    if num_states <= 0:
        raise ValueError(
            "num_states must be greater than zero."
        )

    # --------------------------------------------------------
    # Create output folder
    # --------------------------------------------------------

    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Display equation
    # --------------------------------------------------------

    solve_schrodinger()

    # --------------------------------------------------------
    # Calculate lowest states
    # --------------------------------------------------------

    states = lowest_states(
        a,
        b,
        count=num_states
    )

    print_table(states)

    # --------------------------------------------------------
    # Verify normalization
    # --------------------------------------------------------

    verify(
        a,
        b,
        states
    )

    # --------------------------------------------------------
    # Create mesh
    # --------------------------------------------------------

    X, Y = rectangle_mesh(
        a,
        b,
        grid
    )

    # --------------------------------------------------------
    # Calculate wavefunction and probability
    # --------------------------------------------------------

    fields = fields_on_mesh(
        a,
        b,
        states,
        X,
        Y
    )

    # --------------------------------------------------------
    # Generate plots
    # --------------------------------------------------------

    energy_plot = plot_energy_levels(
        a,
        b,
        states,
        output_dir
    )

    psi_plot = plot_grid(
        a,
        b,
        states,
        X,
        Y,
        fields,
        "psi",
        output_dir
    )

    rho_plot = plot_grid(
        a,
        b,
        states,
        X,
        Y,
        fields,
        "rho",
        output_dir
    )

    # --------------------------------------------------------
    # Return everything to the integrator/frontend
    # --------------------------------------------------------

    return {
        "a": a,
        "b": b,
        "states": states,
        "X": X,
        "Y": Y,
        "fields": fields,
        "energy_plot": str(energy_plot),
        "psi_plot": str(psi_plot),
        "rho_plot": str(rho_plot)
    }
