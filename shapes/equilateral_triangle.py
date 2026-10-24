from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
SQ3 = math.sqrt(3.0)
SIGNS = (1.0, 1.0, 1.0, -1.0, -1.0, -1.0)
DEFAULT_GRID = 120

def validate_inputs(a: float, m: int, n: int) -> tuple[float, int, int]:
    try:
        a_value = float(a)
    except (TypeError, ValueError) as exc:
        raise ValueError('Triangle side a must be a positive number.') from exc
    if not math.isfinite(a_value) or a_value <= 0:
        raise ValueError('Triangle side a must be a positive finite number.')
    for name, value in (('m', m), ('n', n)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError(f'{name} must be an integer from 1 to 5.')
        if not 0 <= int(value) <= 5:
            raise ValueError(f'{name} must be between 1 and 5.')
    m_value, n_value = (int(m), int(n))
    if m_value == 0 or n_value == 0:
        raise ValueError('The current triangular-well equations support positive quantum numbers only. Please choose 1 <= n <= m <= 5; zero is not a valid state in this model.')
    if m_value < n_value:
        raise ValueError('The current derivation uses m >= n. Please choose values with 1 <= n <= m <= 5; the backend will not silently swap them.')
    return (a_value, m_value, n_value)

def triangle_mesh(side: float, divisions: int=DEFAULT_GRID):
    if divisions < 10:
        raise ValueError('divisions must be at least 10.')
    h = side / divisions
    points, ids = ([], {})
    for j in range(divisions + 1):
        for i in range(divisions + 1 - j):
            ids[i, j] = len(points)
            points.append((i * h + j * h / 2.0, j * h * SQ3 / 2.0))
    triangles = []
    for j in range(divisions):
        for i in range(divisions - j):
            triangles.append((ids[i, j], ids[i + 1, j], ids[i, j + 1]))
            if i + j + 2 <= divisions:
                triangles.append((ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]))
    points = np.asarray(points, dtype=float)
    triangulation = mtri.Triangulation(points[:, 0], points[:, 1], np.asarray(triangles, dtype=int))
    return triangulation

def k_vectors(a: float, m: int, n: int):
    kx = 2.0 * math.pi * (2 * m + n) / (3.0 * a)
    ky = 2.0 * math.pi * n / (SQ3 * a)
    return ((kx, ky), (-kx / 2.0 - SQ3 * ky / 2.0, SQ3 * kx / 2.0 - ky / 2.0), (-kx / 2.0 + SQ3 * ky / 2.0, -SQ3 * kx / 2.0 - ky / 2.0), (kx, -ky), (-kx / 2.0 + SQ3 * ky / 2.0, SQ3 * kx / 2.0 + ky / 2.0), (-kx / 2.0 - SQ3 * ky / 2.0, -SQ3 * kx / 2.0 + ky / 2.0))

def wavefunction(a: float, m: int, n: int, kind: str, x, y):
    kind = kind.upper()
    if kind not in ('S', 'C'):
        raise ValueError("kind must be 'S' (sine) or 'C' (cosine).")
    if m == n and kind == 'C':
        raise ValueError("For m == n, the supplied model supports only kind='S'.")
    norm = 2.0 / (3.0 ** 0.75 * a)
    if m == n:
        norm /= math.sqrt(2.0)
    trig = np.sin if kind == 'S' else np.cos
    total = np.zeros_like(np.asarray(x, dtype=float), dtype=float)
    for (kx, ky), sign in zip(k_vectors(a, m, n), SIGNS):
        total += sign * trig(kx * x + ky * y)
    return norm * total

def _phase_strings(a: float, m: int, n: int):
    return [f't{i} = ({kx:.8g})*x + ({ky:.8g})*y' for i, (kx, ky) in enumerate(k_vectors(a, m, n), start=1)]

def _save_state_plot(tri, a, m, n, kind, output_path, probability=False):
    values = wavefunction(a, m, n, kind, tri.x, tri.y)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    if probability:
        plotted = values ** 2
        artist = ax.tripcolor(tri, plotted, shading='gouraud', cmap='viridis')
        title = f'Probability density |psi|^2 — m={m}, n={n}, {kind}'
        colorbar_label = '$|\\psi(x,y)|^2$'
    else:
        plotted = values
        limit = max(float(np.max(np.abs(plotted))), 1e-12)
        artist = ax.tripcolor(tri, plotted, shading='gouraud', cmap='RdBu_r', vmin=-limit, vmax=limit)
        title = f'Normalized wave function psi — m={m}, n={n}, {kind}'
        colorbar_label = '$\\psi(x,y)$'
    fig.colorbar(artist, ax=ax, label=colorbar_label)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, a)
    ax.set_ylim(0, SQ3 * a / 2.0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches='tight')
    plt.close(fig)

def _save_energy_level_plot(a: float, m: int, n: int, energy_coefficient: int, output_path):
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    level = float(energy_coefficient)
    ax.hlines(level, 0.15, 0.75, linewidth=3)
    ax.plot([0.45], [level], marker='o', markersize=6)
    ax.text(0.82, level, f'Selected state (m,n)=({m},{n})\nE = {energy_coefficient} E0', va='center', ha='left')
    ax.set_xlim(0, 1.8)
    ax.set_ylim(0, max(1.0, level * 1.35))
    ax.set_xticks([])
    ax.set_ylabel('E / E0')
    ax.set_title(f'Selected energy level — a={a:g}, m={m}, n={n}')
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches='tight')
    plt.close(fig)

def run_triangle(a: float, m: int, n: int, output_dir='generated_plots', divisions: int=DEFAULT_GRID) -> dict:
    a, m, n = validate_inputs(a, m, n)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    energy_coefficient = m * m + m * n + n * n
    e0_factor = 8.0 * math.pi ** 2 / (9.0 * a ** 2)
    energy_factor = e0_factor * energy_coefficient
    kinds = ['S']
    tri = triangle_mesh(a, divisions)
    plot_paths = {}
    for kind in kinds:
        suffix = f'm{m}_n{n}_{kind.lower()}'
        psi_path = out / f'wavefunction_{suffix}.png'
        prob_path = out / f'probability_density_{suffix}.png'
        _save_state_plot(tri, a, m, n, kind, psi_path, probability=False)
        _save_state_plot(tri, a, m, n, kind, prob_path, probability=True)
        plot_paths[kind] = {'wavefunction_plot': str(psi_path.resolve()), 'probability_density_plot': str(prob_path.resolve())}
    energy_plot_path = out / f'energy_level_m{m}_n{n}.png'
    _save_energy_level_plot(a, m, n, energy_coefficient, energy_plot_path)
    trig_name = 'sin'
    wavefunction_equation = (
        f'psi_S = N * [{trig_name}(t1) + {trig_name}(t2) + '
        f'{trig_name}(t3) - {trig_name}(t4) - '
        f'{trig_name}(t5) - {trig_name}(t6)]'
    )
    normalization = f'2/(3^(3/4)*{a:g})'
    if m == n:
        normalization += '/sqrt(2)'
    return {
        'shape': 'equilateral_triangle',
        'inputs': {
            'a': a,
            'm': m,
            'n': n,
        },
        'potential_inside': 0,
        'valid_quantum_number_range_for_current_model': '1 <= n <= m <= 5',
        'energy': {
            'dimensionless_coefficient_E_over_E0': energy_coefficient,
            'E0_expression': '8*pi^2*hbar^2/(9*M*a^2)',
            'selected_energy_expression': (
                f'E = {energy_coefficient} * 8*pi^2*hbar^2/'
                f'(9*M*({a:g})^2)'
            ),
            'energy_factor_without_hbar2_over_M': energy_factor,
        },
        'normalization_constant': normalization,
        'available_combinations': kinds,
        'phase_terms': _phase_strings(a, m, n),
        'wavefunction_equations': {'S': wavefunction_equation},
    }
