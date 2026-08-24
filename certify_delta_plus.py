#!/usr/bin/env python3
"""表 2 独立宿主的 Delta_plus 全局证书。

指数宿主采用三段论证：解析排除 beta 接近 1 与 beta 足够大的两端，
再对紧区间上的 t=log(beta) 作高精度区间分支定界。区间下界把所有
未展开的尺度层合并为一个尾类，因此不会把截断后的熵误当作原熵。
两个有界宿主使用解析证书。
"""

from __future__ import annotations

import heapq
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import scipy
from scipy.optimize import minimize_scalar


OUT = Path(__file__).resolve().parent / "result" / "delta_plus_certificate.json"
MP_DPS = 80
LAMBDA = mp.mpf("0.5")  # Exp(mean=2B), B=1
BETA_LEFT = mp.mpf("1.05")
BETA_RIGHT = mp.mpf("17")
INITIAL_CELLS = 4096
VALUE_TOL = mp.mpf("5e-5")
ROUND_GUARD = mp.mpf("1e-12")
TAIL_EXPONENT = mp.mpf("38")


def hbit(p: mp.mpf) -> mp.mpf:
    if p <= 0 or p >= 1:
        return mp.mpf("0")
    return -p * mp.log(p, 2)


def exp_cdf(x: mp.mpf) -> mp.mpf:
    return -mp.expm1(-LAMBDA * x)


def exp_point_upper(beta: mp.mpf) -> mp.mpf:
    """在给定 beta 上返回 J_beta 的严格上包络。"""
    t = mp.log(beta)
    r = 1 / beta
    n = max(16, int(mp.ceil(mp.mpf("55") / t)))

    p0 = mp.e ** (-LAMBDA)
    entropy = hbit(p0)
    for k in range(1, n + 1):
        upper = r ** (k - 1)
        lower = r**k
        pk = mp.e ** (-LAMBDA * lower) - mp.e ** (-LAMBDA * upper)
        entropy += hbit(pk)

    # f_D(d) <= lambda 给出 p_k <= lambda(1-r)r^(k-1)。
    c = LAMBDA * (1 - r)
    rn = r**n
    entropy_tail_nats = c * rn * (
        -mp.log(c) / (1 - r)
        - mp.log(r) * (mp.mpf(n) / (1 - r) + r / (1 - r) ** 2)
    )
    entropy += entropy_tail_nats / mp.log(2)

    # E K = sum_{j>=0} F_D(r^j)，余项用 F_D(x) <= lambda x。
    mean_k = mp.fsum(exp_cdf(r**j) for j in range(n + 1))
    mean_k += LAMBDA * r ** (n + 1) / (1 - r)
    return entropy + mean_k * mp.log(beta, 2) + ROUND_GUARD


def exp_cell_lower(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """对 t in [a,b] 给出 J_exp(t) 的保守下界。"""
    n = max(12, int(mp.ceil(TAIL_EXPONENT / a)))
    p0 = mp.e ** (-LAMBDA)
    entropy_lower = hbit(p0)
    mean_lower = mp.mpf("0")

    for k in range(1, n + 1):
        # F(u)-F(l)；利用阈值随 t 单调变化分别包住上下端。
        p_lo = exp_cdf(mp.e ** (-(k - 1) * b)) - exp_cdf(mp.e ** (-k * a))
        p_hi = exp_cdf(mp.e ** (-(k - 1) * a)) - exp_cdf(mp.e ** (-k * b))
        p_lo = max(mp.mpf("0"), p_lo)
        p_hi = min(mp.mpf("1"), max(mp.mpf("0"), p_hi))
        entropy_lower += min(hbit(p_lo), hbit(p_hi))
        mean_lower += k * p_lo

    # K>n 的全部层合并成一个尾类；合并只会降低熵。
    tail_lo = exp_cdf(mp.e ** (-n * b))
    tail_hi = exp_cdf(mp.e ** (-n * a))
    entropy_lower += min(hbit(tail_lo), hbit(tail_hi))
    mean_lower += (n + 1) * tail_lo

    lower = entropy_lower + mean_lower * a / mp.log(2)
    return max(mp.mpf("0"), lower - ROUND_GUARD)


def exp_small_beta_lower(beta_max: mp.mpf) -> mp.mpf:
    """对 1<beta<=beta_max 的解析下界。"""
    q = exp_cdf(mp.mpf("1"))
    pmax = LAMBDA * (1 - 1 / beta_max)
    # H(K)=h_2(q)+qH(K|K>=1)，且 H >= -log p_max。
    return hbit(q) + hbit(1 - q) + q * mp.log(q / pmax, 2)


def exp_large_beta_lower(beta_min: mp.mpf) -> mp.mpf:
    """对 beta>=beta_min 的解析下界。"""
    q = exp_cdf(mp.mpf("1"))
    return hbit(q) + hbit(1 - q) + q * mp.log(beta_min, 2)


def certify_exponential() -> dict:
    with mp.workdps(MP_DPS):
        # 浮点优化只提供一个可行见证；全局性由后续下界证书给出。
        fit = minimize_scalar(
            lambda z: float(exp_point_upper(mp.e ** mp.mpf(str(z)))),
            bounds=(math.log(float(BETA_LEFT)), math.log(float(BETA_RIGHT))),
            method="bounded",
            options={"xatol": 1e-13},
        )
        witness_t = mp.mpf(str(float(fit.x)))
        witness_beta = mp.e**witness_t
        global_upper = exp_point_upper(witness_beta)

        small_lower = exp_small_beta_lower(BETA_LEFT)
        large_lower = exp_large_beta_lower(BETA_RIGHT)
        if min(small_lower, large_lower) <= global_upper:
            raise AssertionError("analytic tail exclusion does not clear witness upper bound")

        left = mp.log(BETA_LEFT)
        right = mp.log(BETA_RIGHT)
        width = (right - left) / INITIAL_CELLS
        heap: list[tuple[mp.mpf, int, mp.mpf, mp.mpf, mp.mpf]] = []
        serial = 0
        for i in range(INITIAL_CELLS):
            a = left + i * width
            b = left + (i + 1) * width
            lb = exp_cell_lower(a, b)
            heapq.heappush(heap, (lb, serial, a, b, lb))
            serial += 1

        splits = 0
        while global_upper - heap[0][4] > VALUE_TOL:
            _, _, a, b, _ = heapq.heappop(heap)
            mid = (a + b) / 2
            for x, y in ((a, mid), (mid, b)):
                lb = exp_cell_lower(x, y)
                heapq.heappush(heap, (lb, serial, x, y, lb))
                serial += 1
            splits += 1
            if splits > 200000:
                raise RuntimeError("branch-and-bound did not converge")

        compact_lower = heap[0][4]
        global_lower = min(compact_lower, small_lower, large_lower)

        # 任何仍可能优于见证上界的单元都必须包含全局实现者。
        possible = [(a, b) for _, _, a, b, lb in heap if lb <= global_upper]
        beta_lo = mp.e ** min(a for a, _ in possible)
        beta_hi = mp.e ** max(b for _, b in possible)

        return {
            "status": "certified",
            "method": "analytic endpoint exclusion plus probability-enclosure branch-and-bound",
            "beta_domain": [1.0, "inf"],
            "compact_beta_interval": [float(BETA_LEFT), float(BETA_RIGHT)],
            "endpoint_lower_bounds": {
                "beta_le_1.05": float(small_lower),
                "beta_ge_17": float(large_lower),
            },
            "beta_star_interval": [float(beta_lo), float(beta_hi)],
            "beta_witness": float(witness_beta),
            "Delta_plus_interval": [float(global_lower), float(global_upper)],
            "certificate_width": float(global_upper - global_lower),
            "initial_cells": INITIAL_CELLS,
            "adaptive_splits": splits,
            "mp_dps": MP_DPS,
            "round_guard": float(ROUND_GUARD),
            "target_width": float(VALUE_TOL),
            "tail_exponent": float(TAIL_EXPONENT),
        }


def main() -> None:
    exp_cert = certify_exponential()
    result = {
        "schema_version": 1,
        "quantity": "Delta_plus=inf_{beta>1} J_beta",
        "scope": "Table 2 independent-host exemplars only",
        "software": {
            "python": platform.python_version(),
            "mpmath": mp.__version__,
            "scipy": scipy.__version__,
        },
        "hosts": {
            "Exp(mean=2B)": exp_cert,
            "U(B/2,3B/2)": {
                "status": "certified",
                "method": "analytic partition-entropy inequality",
                "beta_star_interval": [2.0, 2.0],
                "Delta_plus_interval": [1.5, 1.5],
            },
            "U(B,2B) (sparse)": {
                "status": "certified",
                "method": "analytic K_beta=0 almost surely",
                "beta_star": "any beta>1",
                "Delta_plus_interval": [0.0, 0.0],
            },
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    lo, hi = exp_cert["Delta_plus_interval"]
    blo, bhi = exp_cert["beta_star_interval"]
    print(f"Exp(mean=2B): Delta_plus in [{lo:.8f}, {hi:.8f}]")
    print(f"               beta_star in [{blo:.8f}, {bhi:.8f}]")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
