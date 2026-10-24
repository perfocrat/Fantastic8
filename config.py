"""Shared configuration for Fantastic8 / Quantum Well Explorer.

These values are the intended project-wide defaults. The current shape modules
do not all import this file yet; during integration, refactor their settings
to use these constants or pass values through solver functions.
"""
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent
SHAPES_DIR = PROJECT_ROOT / "shapes"
UI_DIR = PROJECT_ROOT / "ui"
OUTPUT_DIR = PROJECT_ROOT / "generated_plots"

# Streamlit UI
APP_TITLE = "Quantum Well Explorer"
PAGE_ICON = "⚛️"
PAGE_LAYOUT = "wide"

# Numerical defaults
GRID_SIZE = 120
NUM_STATES = 6

# Physics constants:
# Use dimensionless units (hbar = mass = 1) consistently across every solver,
# OR convert every solver to SI units before comparing energies.
# Do not mix these defaults with the SI constants currently in rectangle.py.
HBAR = 1.0
PARTICLE_MASS = 1.0

# Shape labels used by the UI -> backend module / entry point
SHAPE_SOLVERS = {
    "Circle": ("shapes.circle", "run_circle"),
    "Rectangle": ("shapes.rectangle", "run_rectangle"),
    "Equilateral Triangle": ("shapes.equilateral_triangle", "run_triangle"),
    "Ellipse": ("shapes.ellipse", None),  # needs a callable API accepting geometry
    "Circular Sector": ("shapes.sector", None),  # needs a callable API accepting geometry
    "Annulus": ("shapes.annulus", "run_annulus"),
}

# UI geometry field names and intended backend meanings
GEOMETRY_FIELDS = {
    "Circle": ("radius",),
    "Rectangle": ("width", "height"),
    "Equilateral Triangle": ("side",),
    "Ellipse": ("semi_major", "semi_minor"),
    "Circular Sector": ("radius", "angle", "start_angle"),
    "Annulus": ("inner_radius", "outer_radius"),
}

# Quantum-number constraints vary by model; enforce them in the UI/adapter.
TRIANGLE_QUANTUM_NUMBER_MIN = 1
TRIANGLE_QUANTUM_NUMBER_MAX = 5
TRIANGLE_REQUIRES_M_GE_N = True

# Plot generation
SAVE_PLOTS = True
PLOT_DPI = 160

