"""Check annular-state normalization using Gauss-Legendre quadrature.

This independently integrates the normalized wavefunction returned by the
backend over radius and angle, including the polar area element r dr dtheta.
"""
import argparse
import math

import numpy as np

from annulusbackend import wavefunction


def normalization_integral(r_inner: float, r_outer: float, angular_momentum: int,
                           radial_index: int, order: int = 160) -> float:
    if order < 16:
        raise ValueError("Quadrature order must be at least 16.")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    r = (r_outer - r_inner) * (nodes + 1.0) / 2.0 + r_inner
    wr = weights * (r_outer - r_inner) / 2.0
    theta = math.pi * (nodes + 1.0)
    wt = math.pi * weights
    rr, tt = np.meshgrid(r, theta, indexing="ij")
    psi = wavefunction(r_inner, r_outer, angular_momentum, radial_index, rr, tt)
    return float(np.sum(np.abs(psi) ** 2 * rr * wr[:, None] * wt[None, :]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-inner", type=float, default=0.35)
    parser.add_argument("--r-outer", type=float, default=1.0)
    parser.add_argument("--l", type=int, default=1)
    parser.add_argument("--q", type=int, default=1)
    parser.add_argument("--order", type=int, default=160)
    args = parser.parse_args()
    result = normalization_integral(args.r_inner, args.r_outer, args.l, args.q, args.order)
    print(f"State (l,q)=({args.l},{args.q}), annulus={args.r_inner:g}..{args.r_outer:g}")
    print(f"Integral of |psi|^2 over the annulus: {result:.12g}")
    print(f"Absolute difference from 1: {abs(result - 1.0):.3g}")


if __name__ == "__main__":
    main()
