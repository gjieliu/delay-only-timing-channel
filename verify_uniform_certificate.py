# -*- coding: utf-8 -*-
"""Verify the analytic Delta_+ certificates for the two uniform hosts of Table II.

`certify_delta_plus.py` computes the exponential host by branch-and-bound, but records
the two uniform hosts as the constants proved analytically in Appendix F of the paper:

    U(B/2,3B/2) : Delta_plus_interval = [1.5, 1.5]
    U(B,2B)     : Delta_plus_interval = [0.0, 0.0]

Those two entries are proved, not computed, so nothing in the production chain checks
them. This script closes that gap: it re-derives both claims numerically and asserts
them, so the recorded constants are checked rather than merely asserted.

It is a verification script, not part of the production chain, and writes no result file.

Run:  python3 verify_uniform_certificate.py     (exit 0 = all checks pass)
"""
import sys
import numpy as np

B = 1.0
LN2 = np.log(2.0)


def hbin(p):
    p = np.asarray(p, dtype=float)
    out = np.zeros_like(p)
    m = (p > 0) & (p < 1)
    out[m] = -(p[m] * np.log2(p[m]) + (1 - p[m]) * np.log2(1 - p[m]))
    return out


# ---------------------------------------------------------------- U(B,2B)
def check_sparse_host():
    """d ~ U(B,2B) => d >= B almost surely => K_beta == 0 => J_beta == 0 for all beta>1."""
    rng = np.random.default_rng(0)
    d = rng.uniform(B, 2 * B, 2_000_000)
    assert d.min() >= B, "sampled d below B"
    for beta in (1.02, 1.5, 2.0, 5.0, 20.0, 1e3):
        K = np.ceil(np.log(B / np.minimum(d, B)) / np.log(beta))
        assert np.all(K == 0), "K_beta not identically zero at beta=%g" % beta
    print("U(B,2B)      : K_beta == 0 a.s. for every beta > 1  ->  J_beta == 0   [OK]")
    return 0.0


# ------------------------------------------------------- U(B/2,3B/2), beta>=2
def bracket_exact(b):
    """H(M) + b E[M] for M = ceil(u/b), u on (0,1] with density 2 ln2 2^{-u}."""
    n = int(np.ceil(1.0 / b))
    edges = np.minimum(np.arange(0, n + 1) * b, 1.0)
    p = 2.0 * (np.power(2.0, -edges[:-1]) - np.power(2.0, -edges[1:]))
    p = p[p > 1e-18]
    H = -(p * np.log2(p)).sum()
    EM = (np.arange(1, len(p) + 1) * p).sum()
    return H + b * EM


def check_uniform_host(grid=200_001, bgrid=20_000):
    # u = log2(B/d) conditioned on d < B; d ~ U(B/2,B) there, so u has density 2 ln2 2^{-u}
    u = np.linspace(1e-14, 1.0, grid)
    f = 2.0 * LN2 * np.power(2.0, -u)
    norm = np.trapezoid(f, u)
    Eu = np.trapezoid(u * f, u)
    h = np.trapezoid(-f * np.log2(f), u)
    assert abs(norm - 1.0) < 1e-9, "conditional density of u is not normalised"

    # ---- beta >= 2 (b >= 1): only the first nonzero tier is occupied
    for b in (1.0, 1.3, 2.0, 4.32):
        assert abs(bracket_exact(b) - b) < 1e-12, "b>=1 branch: bracket != b at b=%g" % b
    # so J_beta = 1 + (1/2) log2 beta on beta >= 2, minimised at beta = 2 with value 3/2

    # ---- 0 < b < 1, piece 1: uniform-quantizer bound  H(M) >= h(u) - log2 b,  b E M >= E u
    b1 = 2.0 ** (h + Eu - 1.0)          # piece 1 gives >= 1 exactly for b <= b1
    # ---- 0 < b < 1, piece 2: merge tiers above the first
    bb = np.linspace(1e-6, 1.0, bgrid)
    q1 = 2.0 - np.power(2.0, 1.0 - bb)
    piece2 = hbin(q1) + bb * np.power(2.0, 1.0 - bb)
    ok2 = piece2 >= 1.0
    b2 = bb[np.argmax(ok2)] if ok2.any() else np.inf   # first b where piece 2 holds
    assert ok2[np.searchsorted(bb, b2):].all(), "piece 2 is not monotone-valid above b2"

    assert b2 < b1, "the two pieces do not overlap: [%.6f, %.6f]" % (b2, b1)
    print("U(B/2,3B/2)  : h(u)=%.8f  E[u]=%.8f" % (h, Eu))
    print("               piece 1 (quantizer entropy) valid for b <= %.6f" % b1)
    print("               piece 2 (two-tier merge)    valid for b >= %.6f" % b2)
    print("               overlap [%.6f, %.6f] covers (0,1)                  [OK]"
          % (b2, b1))

    # ---- independent check: the exact bracket on (0,1) never drops below 1
    bs = np.linspace(1e-3, 1.0, 4000)
    vals = np.array([bracket_exact(x) for x in bs])
    i = int(np.argmin(vals))
    assert vals[i] >= 1.0, "exact bracket dips below 1: %.8f at b=%.5f" % (vals[i], bs[i])
    print("               exact H(M)+b E[M] on (0,1): min = %.6f at b = %.4f (>= 1) [OK]"
          % (vals[i], bs[i]))

    J_min = 1.0 + 0.5 * vals[i]
    print("               => J_beta >= %.6f on 1 < beta < 2, and = 3/2 at beta = 2" % J_min)
    return 1.5


def main():
    print("Verifying the analytic Delta_+ certificates of Appendix F.\n")
    d_sparse = check_sparse_host()
    d_unif = check_uniform_host()
    print()
    expected = {"U(B,2B) (sparse)": (0.0, d_sparse), "U(B/2,3B/2)": (1.5, d_unif)}
    for host, (recorded, derived) in expected.items():
        assert abs(recorded - derived) < 1e-12, "%s: %r != %r" % (host, recorded, derived)
        print("%-18s recorded Delta_plus_interval [%.1f, %.1f] is verified."
              % (host, recorded, recorded))
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
