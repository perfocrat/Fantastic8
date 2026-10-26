"""Shared application and physics-adapter configuration."""
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent
SHAPES_DIR = PROJECT_ROOT / "shapes"
FRONTEND_DIR = PROJECT_ROOT / "frontend"
OUTPUT_DIR = PROJECT_ROOT / "generated_plots"

APP_TITLE = "Quantum Well Explorer"
SERVER_HOST = "127.0.0.1"
SERVER_PORT = 8000

# Numerical defaults used by the HTML API.
GRID_SIZE = 70
NUM_STATES = 6

# Shape labels and their existing adapter entry points.
SHAPE_SOLVERS = {
    "Circle": ("shapes.adapters", "calculate_circle"),
    "Rectangle": ("shapes.adapters", "calculate_rectangle"),
    "Equilateral Triangle": ("shapes.adapters", "calculate_triangle"),
    "Ellipse": ("shapes.adapters", "calculate_ellipse"),
    "Circular Sector": ("shapes.adapters", "calculate_sector"),
    "Annulus": ("shapes.adapters", "calculate_annulus"),
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

# Quantum-number constraints vary by model; enforce them in the adapter.
TRIANGLE_QUANTUM_NUMBER_MIN = 1
TRIANGLE_QUANTUM_NUMBER_MAX = 5
TRIANGLE_REQUIRES_M_GE_N = True

# Plot generation
SAVE_PLOTS = True
PLOT_DPI = 160

