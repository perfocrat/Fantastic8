
"""Adapter between the Streamlit UI and the rectangle solver."""

from pathlib import Path

from shapes.rectangle import run_rectangle, energy, psi


def calculate_rectangle(width, height, m, n, grid=70):
    """
    width and height must be in metres.
    m and n are positive integer quantum numbers.
    """

    if width <= 0 or height <= 0:
        raise ValueError("Width and height must be positive.")

    if not isinstance(m, int) or not isinstance(n, int) or m < 1 or n < 1:
        raise ValueError("Quantum numbers m and n must be positive integers.")

    # Generate the six lowest states for the energy-level display.
    result = run_rectangle(
        a=width,
        b=height,
        num_states=6,
        grid=grid,
        output_dir=str(Path("generated_plots") / "rectangle"),
    )

    # Calculate the requested state directly, even if it is not
    # among the six lowest-energy states.
    selected_energy = energy(width, height, m, n)
    selected_wavefunction = psi(
        width, height, m, n, result["X"], result["Y"]
    )
    selected_density = abs(selected_wavefunction) ** 2

    return {
        "shape": "Rectangle",
        "inputs": {
            "width": width,
            "height": height,
            "m": m,
            "n": n,
        },
        "energy": selected_energy,
        "wavefunction": selected_wavefunction,
        "probability_density": selected_density,
        "X": result["X"],
        "Y": result["Y"],
        "states": result["states"],
        "energy_plot": result["energy_plot"],
        "psi_plot": result["psi_plot"],
        "rho_plot": result["rho_plot"],
    }
