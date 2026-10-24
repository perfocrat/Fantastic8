
"""
Quantum Well Explorer
---------------------
Streamlit UI for exploring two-dimensional infinite quantum wells.

This module collects geometry and quantum-state inputs, dispatches the
selected shape solver, and renders its wavefunction, probability density,
and energy results.
"""

import sys
from pathlib import Path

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from matplotlib.patches import (
    Polygon,
    Rectangle,
    Circle,
    Wedge,
    Ellipse,
)

from shapes.adapters import calculate_shape

# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Quantum Well Explorer",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ============================================================
# CUSTOM STYLING
# ============================================================

st.markdown(
    """
    <style>
    .stApp {
        background-color: #0b1020;
        color: #e8edf7;
    }

    [data-testid="stSidebar"] {
        background-color: #111a2e;
        border-right: 1px solid #25334d;
    }

    h1, h2, h3 {
        color: #e8edf7;
    }

    .main-title {
        font-size: 2.5rem;
        font-weight: 750;
        color: #e8edf7;
        margin-bottom: 0.1rem;
    }

    .subtitle {
        color: #9eacc5;
        font-size: 1rem;
        margin-bottom: 1.5rem;
    }

    .info-card {
        background: #141f35;
        border: 1px solid #293954;
        border-radius: 12px;
        padding: 16px;
        margin-bottom: 12px;
    }

    .info-label {
        color: #9eacc5;
        font-size: 0.85rem;
    }

    .info-value {
        color: #f0f4ff;
        font-size: 1.5rem;
        font-weight: 700;
    }

    .small-note {
        color: #9eacc5;
        font-size: 0.82rem;
    }

    div[data-testid="stMetric"] {
        background-color: #141f35;
        border: 1px solid #293954;
        padding: 14px;
        border-radius: 10px;
    }

    </style>
    """,
    unsafe_allow_html=True,
)

# ============================================================
# CONSTANTS
# ============================================================

SHAPES = [
    "Equilateral Triangle",
    "Rectangle",
    "Circle",
    "Circular Sector",
    "Annulus",
    "Ellipse",
]

SHAPE_DESCRIPTIONS = {
    "Equilateral Triangle": "A particle confined inside an equilateral triangular boundary.",
    "Rectangle": "A particle confined inside a rectangular boundary.",
    "Circle": "A particle confined inside a circular boundary.",
    "Circular Sector": "A particle confined inside a sector of a circle.",
    "Annulus": "A particle confined between two concentric circular boundaries.",
    "Ellipse": "A particle confined inside an elliptical boundary.",
}

# ============================================================
# SIDEBAR: GEOMETRY INPUTS
# ============================================================

st.sidebar.title("Well Configuration")
st.sidebar.caption("Configure your quantum well.")

shape = st.sidebar.selectbox(
    "Well geometry",
    SHAPES,
    help="Choose the shape of the two-dimensional quantum well.",
)

st.sidebar.divider()
st.sidebar.subheader("Geometric Parameters")

geometry = {}

if shape == "Equilateral Triangle":

    geometry["side"] = st.sidebar.slider(
        "Side length, a (nm)",
        min_value=1.0,
        max_value=20.0,
        value=5.0,
        step=0.5,
        help="Length of each side of the equilateral triangle.",
    )

elif shape == "Rectangle":

    geometry["width"] = st.sidebar.slider(
        "Width, W (nm)",
        min_value=1.0,
        max_value=20.0,
        value=8.0,
        step=0.5,
        help="Horizontal dimension of the rectangular well.",
    )

    geometry["height"] = st.sidebar.slider(
        "Height, H (nm)",
        min_value=1.0,
        max_value=20.0,
        value=5.0,
        step=0.5,
        help="Vertical dimension of the rectangular well.",
    )

elif shape == "Circle":

    geometry["radius"] = st.sidebar.slider(
        "Radius, R (nm)",
        min_value=1.0,
        max_value=10.0,
        value=5.0,
        step=0.5,
        help="Radius of the circular well.",
    )

elif shape == "Circular Sector":

    geometry["radius"] = st.sidebar.slider(
        "Radius, R (nm)",
        min_value=1.0,
        max_value=10.0,
        value=5.0,
        step=0.5,
        help="Radius of the circular sector.",
    )

    geometry["angle"] = st.sidebar.slider(
        "Sector angle, θ (degrees)",
        min_value=15,
        max_value=345,
        value=90,
        step=15,
        help="Angular opening of the sector.",
    )

    geometry["start_angle"] = st.sidebar.slider(
        "Starting angle, φ (degrees)",
        min_value=0,
        max_value=345,
        value=0,
        step=15,
        help="Orientation of the first radial boundary.",
    )

elif shape == "Annulus":

    geometry["outer_radius"] = st.sidebar.slider(
        "Outer radius, Ro (nm)",
        min_value=2.0,
        max_value=10.0,
        value=5.0,
        step=0.5,
        help="Radius of the outer circular boundary.",
    )

    geometry["inner_radius"] = st.sidebar.slider(
        "Inner radius, Ri (nm)",
        min_value=0.5,
        max_value=min(9.0, geometry["outer_radius"] - 0.5),
        value=min(2.0, geometry["outer_radius"] - 0.5),
        step=0.5,
        help="Radius of the inner circular boundary. It must be smaller than Ro.",
    )

elif shape == "Ellipse":

    geometry["semi_major"] = st.sidebar.slider(
        "Semi-major axis, A (nm)",
        min_value=1.0,
        max_value=10.0,
        value=6.0,
        step=0.5,
        help="Half of the ellipse's longer diameter.",
    )

    geometry["semi_minor"] = st.sidebar.slider(
        "Semi-minor axis, B (nm)",
        min_value=0.1,
        max_value=geometry["semi_major"] - 0.1,
        value=min(3.0, geometry["semi_major"] - 0.1),
        step=0.1,
        help="Half of the ellipse's shorter diameter; must be smaller than A.",
    )

# ============================================================
# SIDEBAR: QUANTUM NUMBERS
# ============================================================

st.sidebar.divider()
st.sidebar.subheader("Quantum State")

if shape == "Equilateral Triangle":
    m_label = "Quantum number, m"
    n_label = "Quantum number, n"
    m_minimum = n_minimum = 1
    m_maximum = 5
    n_maximum = None
    state_help = "The triangle solver requires 1 <= n <= m <= 5."
elif shape == "Circle":
    m_label = "Angular mode, m"
    n_label = "Radial index, n"
    m_minimum = 0
    n_minimum = 1
    m_maximum = n_maximum = 10
    state_help = "The circular solver uses m >= 0 and n >= 1."
elif shape == "Annulus":
    m_label = "Angular momentum, l"
    n_label = "Radial index, q"
    m_minimum = 0
    n_minimum = 1
    m_maximum = n_maximum = 10
    state_help = "The annulus solver uses l >= 0 and q >= 1."
elif shape == "Circular Sector":
    m_label = "Angular index, ell"
    n_label = "Radial index, s"
    m_minimum = n_minimum = 1
    m_maximum = n_maximum = 10
    state_help = "The sector solver uses ell >= 1 and s >= 1."
elif shape == "Ellipse":
    m_label = "Angular mode, m"
    n_label = "Radial index, j"
    m_minimum = 0
    n_minimum = 1
    m_maximum = n_maximum = 10
    state_help = "The ellipse solver uses m >= 0 and j >= 1 (cosine family)."
else:
    m_label = "Quantum number, m"
    n_label = "Quantum number, n"
    m_minimum = n_minimum = 1
    m_maximum = n_maximum = 10
    state_help = "Rectangle states use positive integer indices."

m = st.sidebar.slider(
    m_label,
    min_value=m_minimum,
    max_value=m_maximum,
    value=1,
    step=1,
    help=state_help,
    key="triangle_m" if shape == "Equilateral Triangle" else None,
)

if shape == "Equilateral Triangle":
    st.session_state.setdefault("triangle_n", n_minimum)
    st.session_state["triangle_n"] = min(
        max(n_minimum, st.session_state["triangle_n"]),
        m,
    )
    n = st.sidebar.slider(
        n_label,
        min_value=n_minimum,
        max_value=max(m, n_minimum + 1),
        step=1,
        help=f"{state_help} When m=1, n is fixed at 1.",
        key="triangle_n",
        disabled=m == n_minimum,
    )
else:
    n = st.sidebar.slider(
        n_label,
        min_value=n_minimum,
        max_value=n_maximum,
        value=1,
        step=1,
        help=state_help,
    )

st.sidebar.divider()
st.sidebar.subheader("Visualisation Settings")

resolution = st.sidebar.slider(
    "Grid resolution",
    min_value=30,
    max_value=150,
    value=70,
    step=10,
    help="Number of points used along each axis of the solver display grid.",
)

display_mode = st.sidebar.radio(
    "Visualisation",
    [
        "Both plots",
        "Wavefunction only",
        "Probability density only",
    ],
)

# ============================================================
# GEOMETRY HELPERS
# ============================================================

def get_domain(shape, geometry):
    """
    Returns the coordinate limits and a mask describing the
    selected geometry.

    The mask is used only for visualisation.
    """

    if shape == "Equilateral Triangle":
        a = geometry["side"]
        x_min, x_max = -a / 2, a / 2
        y_min, y_max = 0, np.sqrt(3) * a / 2

        def mask(X, Y):
            return (
                (Y >= 0)
                & (Y <= np.sqrt(3) * (X + a / 2))
                & (Y <= np.sqrt(3) * (a / 2 - X))
            )

    elif shape == "Rectangle":
        W = geometry["width"]
        H = geometry["height"]

        x_min, x_max = -W / 2, W / 2
        y_min, y_max = -H / 2, H / 2

        def mask(X, Y):
            return (
                (np.abs(X) <= W / 2)
                & (np.abs(Y) <= H / 2)
            )

    elif shape == "Circle":
        R = geometry["radius"]

        x_min, x_max = -R, R
        y_min, y_max = -R, R

        def mask(X, Y):
            return X**2 + Y**2 <= R**2

    elif shape == "Circular Sector":
        R = geometry["radius"]
        angle = geometry["angle"]
        start = geometry["start_angle"]

        x_min, x_max = -R, R
        y_min, y_max = -R, R

        def mask(X, Y):
            theta = np.degrees(np.arctan2(Y, X))
            relative_angle = (theta - start) % 360

            return (
                (X**2 + Y**2 <= R**2)
                & (relative_angle <= angle)
            )

    elif shape == "Annulus":
        Ri = geometry["inner_radius"]
        Ro = geometry["outer_radius"]

        x_min, x_max = -Ro, Ro
        y_min, y_max = -Ro, Ro

        def mask(X, Y):
            radius_squared = X**2 + Y**2

            return (
                (radius_squared >= Ri**2)
                & (radius_squared <= Ro**2)
            )

    else:  # Ellipse
        A = geometry["semi_major"]
        B = geometry["semi_minor"]

        x_min, x_max = -A, A
        y_min, y_max = -B, B

        def mask(X, Y):
            return (
                (X / A)**2 + (Y / B)**2 <= 1
            )

    return x_min, x_max, y_min, y_max, mask


# ============================================================
# GEOMETRY PREVIEW
# ============================================================

def draw_geometry(ax, shape, geometry):
    """Draws the selected well boundary."""

    ax.set_aspect("equal")

    if shape == "Equilateral Triangle":
        a = geometry["side"]

        vertices = [
            (-a / 2, 0),
            (a / 2, 0),
            (0, np.sqrt(3) * a / 2),
        ]

        patch = Polygon(
            vertices,
            closed=True,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        ax.add_patch(patch)

    elif shape == "Rectangle":
        W = geometry["width"]
        H = geometry["height"]

        patch = Rectangle(
            (-W / 2, -H / 2),
            W,
            H,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        ax.add_patch(patch)

    elif shape == "Circle":
        R = geometry["radius"]

        patch = Circle(
            (0, 0),
            R,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        ax.add_patch(patch)

    elif shape == "Circular Sector":
        R = geometry["radius"]
        start = geometry["start_angle"]
        end = start + geometry["angle"]

        patch = Wedge(
            (0, 0),
            R,
            start,
            end,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        ax.add_patch(patch)

    elif shape == "Annulus":
        Ri = geometry["inner_radius"]
        Ro = geometry["outer_radius"]

        outer = Circle(
            (0, 0),
            Ro,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        inner = Circle(
            (0, 0),
            Ri,
            fill=False,
            edgecolor="#ffb454",
            linewidth=2.5,
        )

        ax.add_patch(outer)
        ax.add_patch(inner)

    elif shape == "Ellipse":
        A = geometry["semi_major"]
        B = geometry["semi_minor"]

        patch = Ellipse(
            (0, 0),
            2 * A,
            2 * B,
            fill=False,
            edgecolor="#53b9ff",
            linewidth=2.5,
        )

        ax.add_patch(patch)

    ax.set_title(
        "Well Geometry",
        color="white",
        fontsize=13,
        pad=12,
    )

    ax.set_xlabel("x (nm)", color="#c7d3e8")
    ax.set_ylabel("y (nm)", color="#c7d3e8")

    ax.tick_params(colors="#9eacc5")
    ax.grid(alpha=0.15)

    for spine in ax.spines.values():
        spine.set_color("#34435c")

    ax.set_facecolor("#111a2e")


# ============================================================
# MAIN HEADER
# ============================================================

st.markdown(
    '<div class="main-title">Quantum Well Explorer</div>',
    unsafe_allow_html=True,
)

st.markdown(
    '<div class="subtitle">'
    'An interactive visualisation interface for two-dimensional '
    'infinite quantum wells.'
    '</div>',
    unsafe_allow_html=True,
)

st.info("Each selection is calculated by its corresponding shape solver.")

# ============================================================
# SELECTED PARAMETERS
# ============================================================

st.subheader("Selected Configuration")

col1, col2, col3 = st.columns(3)

col1.metric("Well geometry", shape)
first_index = m_label.rsplit(",", maxsplit=1)[-1].strip()
second_index = n_label.rsplit(",", maxsplit=1)[-1].strip()
col2.metric("Quantum state", f"{first_index}={m}, {second_index}={n}")
col3.metric("Grid resolution", f"{resolution} × {resolution}")

st.caption(SHAPE_DESCRIPTIONS[shape])

# ============================================================
# GEOMETRY SUMMARY
# ============================================================

with st.expander("View geometric parameters", expanded=True):

    parameter_columns = st.columns(3)

    parameter_items = []

    if shape == "Equilateral Triangle":
        parameter_items = [
            ("Side length, a", f"{geometry['side']:.1f} nm"),
        ]

    elif shape == "Rectangle":
        parameter_items = [
            ("Width, W", f"{geometry['width']:.1f} nm"),
            ("Height, H", f"{geometry['height']:.1f} nm"),
        ]

    elif shape == "Circle":
        parameter_items = [
            ("Radius, R", f"{geometry['radius']:.1f} nm"),
        ]

    elif shape == "Circular Sector":
        parameter_items = [
            ("Radius, R", f"{geometry['radius']:.1f} nm"),
            ("Sector angle, θ", f"{geometry['angle']}°"),
            ("Starting angle, φ", f"{geometry['start_angle']}°"),
        ]

    elif shape == "Annulus":
        parameter_items = [
            ("Outer radius, Ro", f"{geometry['outer_radius']:.1f} nm"),
            ("Inner radius, Ri", f"{geometry['inner_radius']:.1f} nm"),
            (
                "Ri / Ro",
                f"{geometry['inner_radius'] / geometry['outer_radius']:.2f}",
            ),
        ]

    elif shape == "Ellipse":
        parameter_items = [
            ("Semi-major axis, A", f"{geometry['semi_major']:.1f} nm"),
            ("Semi-minor axis, B", f"{geometry['semi_minor']:.1f} nm"),
        ]

    for i, (label, value) in enumerate(parameter_items):
        with parameter_columns[i % 3]:
            st.metric(label, value)

# ============================================================
# CALCULATE THE SELECTED STATE
# ============================================================

x_min, x_max, y_min, y_max, _ = get_domain(shape, geometry)
adapter_inputs = {"m": m, "n": n, "grid": resolution}
if shape == "Equilateral Triangle":
    adapter_inputs.update(side=geometry["side"] * 1e-9)
elif shape == "Rectangle":
    adapter_inputs.update(
        width=geometry["width"] * 1e-9,
        height=geometry["height"] * 1e-9,
    )
elif shape == "Circle":
    adapter_inputs.update(radius=geometry["radius"] * 1e-9)
elif shape == "Circular Sector":
    adapter_inputs.update(
        radius=geometry["radius"] * 1e-9,
        angle=geometry["angle"],
        start_angle=geometry["start_angle"],
    )
elif shape == "Annulus":
    adapter_inputs.update(
        r_inner=geometry["inner_radius"] * 1e-9,
        r_outer=geometry["outer_radius"] * 1e-9,
        angular_momentum=m,
        radial_index=n,
    )
    adapter_inputs.pop("m")
    adapter_inputs.pop("n")
elif shape == "Ellipse":
    adapter_inputs.update(
        semi_major=geometry["semi_major"] * 1e-9,
        semi_minor=geometry["semi_minor"] * 1e-9,
    )

solver_result = None
solver_error = None
try:
    solver_result = calculate_shape(shape, **adapter_inputs)
    X = np.asarray(solver_result["X"], dtype=float) * 1e9
    Y = np.asarray(solver_result["Y"], dtype=float) * 1e9
    field = np.asarray(solver_result["wavefunction"], dtype=float)
    probability_density = np.asarray(solver_result["probability_density"], dtype=float)
except Exception as exc:
    solver_result = None
    solver_error = f"{type(exc).__name__}: {exc}"
    st.error(f"{shape} solver could not calculate this state: {solver_error}")

# ============================================================
# VISUALISATION SECTION
# ============================================================

st.divider()
st.subheader("Well Visualisation")

geometry_column, wavefunction_column = st.columns(2)

with geometry_column:

    fig_geometry, ax_geometry = plt.subplots(
        figsize=(5, 4),
    )

    fig_geometry.patch.set_facecolor("#0b1020")

    draw_geometry(
        ax_geometry,
        shape,
        geometry,
    )

    ax_geometry.set_xlim(x_min * 1.15, x_max * 1.15)
    ax_geometry.set_ylim(y_min - 0.15 * (y_max-y_min),
                         y_max + 0.15 * (y_max-y_min))

    st.pyplot(fig_geometry)
    plt.close(fig_geometry)

with wavefunction_column:

    st.markdown("**Selected quantum state**")
    st.markdown(
        f"Selected solver state: `{solver_result['inputs']}`"
        if solver_result is not None
        else f"Requested state: `{m}, {n}`"
    )

    if solver_result is not None:
        st.markdown(f"Calculated using the {solver_result['shape']} solver.")
    else:
        st.markdown("No field is shown because the solver did not return a result.")

# ============================================================
# WAVEFUNCTION AND PROBABILITY PLOTS
# ============================================================

if solver_result is not None and display_mode in ["Both plots", "Wavefunction only"]:

    st.markdown("### Wavefunction preview")

    fig_wave, ax_wave = plt.subplots(
        figsize=(8, 4.5),
    )

    fig_wave.patch.set_facecolor("#0b1020")
    ax_wave.set_facecolor("#111a2e")

    contour = ax_wave.contourf(
        X,
        Y,
        field,
        levels=30,
        cmap="coolwarm",
    )

    fig_wave.colorbar(
        contour,
        ax=ax_wave,
        label="Wavefunction amplitude",
    )

    ax_wave.set_title(f"Wavefunction ψ — {shape}", color="white")

    ax_wave.set_xlabel("x (nm)", color="white")
    ax_wave.set_ylabel("y (nm)", color="white")
    ax_wave.tick_params(colors="white")
    ax_wave.set_aspect("equal")

    st.pyplot(fig_wave)
    plt.close(fig_wave)

if solver_result is not None and display_mode in ["Both plots", "Probability density only"]:

    st.markdown("### Probability density preview")

    fig_probability, ax_probability = plt.subplots(
        figsize=(8, 4.5),
    )

    fig_probability.patch.set_facecolor("#0b1020")
    ax_probability.set_facecolor("#111a2e")

    contour = ax_probability.contourf(
        X,
        Y,
        probability_density,
        levels=30,
        cmap="viridis",
    )

    fig_probability.colorbar(
        contour,
        ax=ax_probability,
        label="Probability density |ψ|²",
    )

    ax_probability.set_title(f"Probability density — {shape}", color="white")

    ax_probability.set_xlabel("x (nm)", color="white")
    ax_probability.set_ylabel("y (nm)", color="white")
    ax_probability.tick_params(colors="white")
    ax_probability.set_aspect("equal")

    st.pyplot(fig_probability)
    plt.close(fig_probability)

# ============================================================
# ENERGY LEVEL SECTION
# ============================================================

st.divider()
st.subheader("Energy Information")

energy_col1, energy_col2, energy_col3 = st.columns(3)

if solver_result is not None:
    states = solver_result["states"]
    ground_state = min(states, key=lambda state: state["E"])
    selected_energy = float(solver_result["energy"])
    ground_energy = float(ground_state["E"])
    energy_difference = selected_energy - ground_energy

    energy_col1.metric("Ground-state energy", f"{ground_energy:.4e} J")
    energy_col2.metric("Selected-state energy", f"{selected_energy:.4e} J")
    energy_col3.metric("Energy difference", f"{energy_difference:.4e} J")
    st.success(f"Energy and fields are from the {solver_result['shape']} solver.")
    fig_energy, ax_energy = plt.subplots(figsize=(7, 3.5))
    fig_energy.patch.set_facecolor("#0b1020")
    for index, state in enumerate(states, start=1):
        ax_energy.hlines(state["E"], index - 0.35, index + 0.35, color="#53b9ff", linewidth=2)
        ax_energy.text(
            index,
            state["E"],
            f"({state['m']},{state['n']})",
            color="#e8edf7",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax_energy.set_xlim(0.5, len(states) + 0.5)
    ax_energy.set_xlabel("Returned state")
    ax_energy.set_ylabel("Energy (J)")
    ax_energy.set_title(f"Energy levels — {shape}", color="white")
    ax_energy.tick_params(colors="white")
    ax_energy.grid(axis="y", alpha=0.25)
    st.pyplot(fig_energy)
    plt.close(fig_energy)
else:
    energy_col1.metric("Ground-state energy", "Not calculated")
    energy_col2.metric("Selected-state energy", "Not calculated")
    energy_col3.metric("Energy difference", "Not calculated")
    st.warning(
        "Energy and fields are unavailable because the selected solver failed."
    )

# ============================================================
# QUANTUM NUMBER INFORMATION
# ============================================================

with st.expander("Understanding the quantum numbers"):

    st.markdown(
        """
        The quantum numbers identify a particular state of the
        confined particle.

        - **m:** First state index.
        - **n:** Second state index.

        Their mathematical interpretation depends on the geometry.
        Different geometries can use different quantum-number conventions.

        Each shape uses its own solver and quantum-number convention.
        """
    )

# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "Quantum Well Explorer | Six shape-specific solvers"
)
