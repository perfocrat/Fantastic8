
import math
import numpy as np

from trianglebackend import k_vectors, SIGNS, SQ3


def raw_wavefunction(a, m, n, x, y):
    """Wavefunction sum before applying any normalization constant."""
    result = np.zeros(np.broadcast_shapes(np.shape(x), np.shape(y)), dtype=float)
    for (kx, ky), sign in zip(k_vectors(a, m, n), SIGNS):
        result += sign * np.sin(kx * x + ky * y)

    return result


def integrate_raw_probability(a, m, n, order=160):
    """
    Numerically integrate |raw_wavefunction|^2 over the triangle.
    Does not use the normalization constant from the backend.
    """
    nodes, weights = np.polynomial.legendre.leggauss(order)

    # Transform quadrature nodes to t in [0, 1].
    t = (nodes + 1.0) / 2.0
    wt = weights / 2.0

    total_integral = 0.0

    # Split at x = a/2, where the triangle's upper boundary changes.
    for x_start, x_end, left_half in [
        (0.0, a / 2.0, True),
        (a / 2.0, a, False),
    ]:
        x = x_start + (nodes + 1.0) * (x_end - x_start) / 2.0
        wx = weights * (x_end - x_start) / 2.0

        if left_half:
            ymax = SQ3 * x
        else:
            ymax = SQ3 * (a - x)

        # Map t in [0,1] to y in [0, ymax(x)].
        xx = x[:, None]
        yy = ymax[:, None] * t[None, :]

        raw = raw_wavefunction(a, m, n, xx, yy)

        # dy = ymax(x) dt: include the Jacobian.
        integrand = raw**2 * ymax[:, None]

        total_integral += np.sum(
            wx[:, None] * wt[None, :] * integrand
        )

    return float(total_integral)


def verify_state(a, m, n, order=160):
    integral = integrate_raw_probability(a, m, n, order)

    # Infer the normalization constant required to make the integral 1.
    inferred_N = 1.0 / math.sqrt(integral)

    expected_N = 2.0 / (3.0**0.75 * a)
    if m == n:
        expected_N /= math.sqrt(2.0)

    relative_error = abs(inferred_N - expected_N) / expected_N

    print(f"(m,n)=({m},{n})")
    print(f"  Raw integral       : {integral:.12g}")
    print(f"  Inferred N         : {inferred_N:.12g}")
    print(f"  Backend N formula  : {expected_N:.12g}")
    print(f"  Relative difference: {relative_error:.3e}")

    return relative_error


if __name__ == "__main__":
    a = 5.0
    errors = []

    for m in range(1, 6):
        for n in range(1, m + 1):
            errors.append(verify_state(a, m, n))
    
    print("\nLargest relative difference:", max(errors))
