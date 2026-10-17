import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

NUM_STATES = 6
GRID = 120

SQ3 = np.sqrt(3.0)
SIGNS = (+1, +1, +1, -1, -1, -1)

x, y, a, M, hbar = sp.symbols("x y a M hbar", positive=True)
m, n = sp.symbols("m n", positive=True, integer=True)
s3 = sp.sqrt(3)


def show(label, expr):
    print(label)
    print(sp.pretty(expr, use_unicode=False))
    print()


def solve_schrodinger():
    kx = 2 * sp.pi * (2 * m + n) / (3 * a)
    ky = 2 * sp.pi * n / (s3 * a)
    vecs = [
        (kx, ky),
        (-kx / 2 - s3 * ky / 2, s3 * kx / 2 - ky / 2),
        (-kx / 2 + s3 * ky / 2, -s3 * kx / 2 - ky / 2),
        (kx, -ky),
        (-kx / 2 + s3 * ky / 2, s3 * kx / 2 + ky / 2),
        (-kx / 2 - s3 * ky / 2, -s3 * kx / 2 + ky / 2),
    ]
    phases = [sp.factor(sp.simplify(p * x + q * y)) for p, q in vecs]
    energy = sp.factor(sp.simplify(hbar ** 2 * (vecs[0][0] ** 2 + vecs[0][1] ** 2) / (2 * M)))
    for f in (sp.sin, sp.cos):
        psi_sym = sum(s * f(t) for s, t in zip(SIGNS, phases))
        residual = sp.simplify(-hbar ** 2 / (2 * M) * (psi_sym.diff(x, 2) + psi_sym.diff(y, 2))
                               - energy * psi_sym)
        assert residual == 0
        for line in (0, s3 * x, s3 * (a - x)):
            assert sp.simplify(psi_sym.subs(y, line)) == 0
    area = s3 * a ** 2 / 4
    N_diff = sp.simplify(sp.sqrt(1 / (3 * area)))
    N_same = sp.simplify(N_diff / sp.sqrt(2))
    return phases, energy, N_diff, N_same


def print_results():
    phases, energy, N_diff, N_same = solve_schrodinger()
    print("NORMALIZED ENERGY EIGENFUNCTIONS")
    print()
    print("psi_S(m,n)(x,y) = N [ sin(t1) + sin(t2) + sin(t3) - sin(t4) - sin(t5) - sin(t6) ]")
    print("psi_C(m,n)(x,y) = N [ cos(t1) + cos(t2) + cos(t3) - cos(t4) - cos(t5) - cos(t6) ]")
    print()
    for j, t in enumerate(phases, 1):
        show(f"t{j} =", t)
    show("N (m != n) =", N_diff)
    show("N (m == n, only psi_S exists) =", N_same)
    print("ENERGY EIGENVALUE")
    print()
    show("E(m, n) =  (m >= n >= 1)", energy)
    E0 = 8 * sp.pi ** 2 * hbar ** 2 / (9 * M * a ** 2)
    show("E(m, n) = E0 (m^2 + m n + n^2),  where E0 =", E0)


def k_vectors(m, n):
    kx = 2 * np.pi * (2 * m + n) / 3.0
    ky = 2 * np.pi * n / SQ3
    return [
        (kx, ky),
        (-kx / 2 - SQ3 * ky / 2, SQ3 * kx / 2 - ky / 2),
        (-kx / 2 + SQ3 * ky / 2, -SQ3 * kx / 2 - ky / 2),
        (kx, -ky),
        (-kx / 2 + SQ3 * ky / 2, SQ3 * kx / 2 + ky / 2),
        (-kx / 2 - SQ3 * ky / 2, -SQ3 * kx / 2 + ky / 2),
    ]


def psi(m, n, kind, x, y):
    N = 2.0 / 3.0 ** 0.75
    if m == n:
        N /= np.sqrt(2.0)
    f = np.sin if kind == "S" else np.cos
    total = 0.0
    for (kx, ky), s in zip(k_vectors(m, n), SIGNS):
        total = total + s * f(kx * x + ky * y)
    return N * total


def lowest_states(count):
    states = []
    K = int(np.ceil(np.sqrt(count))) + 8
    for m in range(1, K):
        for n in range(1, m + 1):
            states.append((m * m + m * n + n * n, m, n, "S"))
            if m != n:
                states.append((m * m + m * n + n * n, m, n, "C"))
    states.sort()
    return states[:count]


def triangle_mesh(N):
    h = 1.0 / N
    pts, ids = [], {}
    for j in range(N + 1):
        for i in range(N + 1 - j):
            ids[(i, j)] = len(pts)
            pts.append((i * h + j * h / 2, j * h * SQ3 / 2))
    tris = []
    for j in range(N):
        for i in range(N - j):
            tris.append((ids[(i, j)], ids[(i + 1, j)], ids[(i, j + 1)]))
            if i + j + 2 <= N:
                tris.append((ids[(i + 1, j)], ids[(i + 1, j + 1)], ids[(i, j + 1)]))
    pts = np.array(pts)
    return mtri.Triangulation(pts[:, 0], pts[:, 1], np.array(tris))


def plot_energy_levels(states):
    groups = []
    for i, (v, m, n, kind) in enumerate(states, 1):
        if groups and groups[-1][0] == v:
            groups[-1][1].append((i, m, n, kind))
        else:
            groups.append([v, [(i, m, n, kind)]])
    plt.figure(figsize=(8, 5.5))
    for v, mem in groups:
        plt.hlines(v, 0.15, 0.75, linewidth=2)
        names = ", ".join(str(a[0]) for a in mem)
        tag = "State" if len(mem) == 1 else f"States ({len(mem)}-fold)"
        mn = f"(m,n)=({mem[0][1]},{mem[0][2]})"
        plt.text(0.82, v, f"{tag} {names}: {mn}\nE = {v} E0", va="center")
    plt.xlim(0, 2.2)
    plt.xticks([])
    plt.ylabel("E / E0     (E0 = 8 pi^2 hbar^2 / 9 M a^2)")
    plt.title("Energy levels: E(m,n) = E0 (m^2 + m n + n^2)")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()


def plot_grid(states, tri, what):
    cols = 3
    rows = int(np.ceil(len(states) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.9 * rows),
                             squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (v, m, n, kind) in zip(axes.ravel(), states):
        ax.axis("on")
        f = psi(m, n, kind, tri.x, tri.y)
        if f[np.argmax(np.abs(f))] < 0:
            f = -f
        if what == "psi":
            vmax = np.max(np.abs(f))
            art = ax.tripcolor(tri, f, shading="gouraud", cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax)
        else:
            art = ax.tripcolor(tri, f ** 2, shading="gouraud", cmap="viridis")
        ax.set_title(f"({m},{n}) {kind},  E = {v} E0", fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlabel("x / a")
        ax.set_ylabel("y / a")
        fig.colorbar(art, ax=ax, fraction=0.046, pad=0.04)
    if what == "psi":
        fig.suptitle(r"Normalized wave functions  $\psi \cdot a$", fontsize=13)
    else:
        fig.suptitle(r"Probability density  $|\psi|^2 \cdot a^2$", fontsize=13)
    fig.tight_layout()


def main():
    print_results()
    states = lowest_states(NUM_STATES)
    tri = triangle_mesh(GRID)
    plot_energy_levels(states)
    plot_grid(states, tri, "psi")
    plot_grid(states, tri, "rho")
    plt.show()


if __name__ == "__main__":
    main()
