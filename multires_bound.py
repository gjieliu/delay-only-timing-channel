"""平稳遍历随机宿主的多分辨率模格可达界。

归一化 B=1。对宿主间隔 D 定义几何分层

    K=0,                         D >= B,
    K=k>=1,   B beta^{-k} <= D < B beta^{-(k-1)},

并在第 k 层使用周期 m_k=B beta^{-k}。逐包可达发送区间的长度满足
L_i >= min(B,D_i) >= m_{K_i}，故编码端总能把发送 IPD 放入指定的
模 m_{K_i} 陪集。一般平稳遍历层过程的状态移除代价是熵率
Hbar(K)。本脚本中的边缘分布基线取独立宿主，因而 Hbar(K)=H(K_1)：

    C >= sum_k p_k C_circle(m_k/sigma) - Hbar(K).

这里 C_circle(r) 是周长与高斯标准差之比为 r 的模环加性高斯信道容量。
脚本同时计算几何分层常数

    J_beta = Hbar(K) + E[K] log2(beta),

以及旧式单阈值界 q C_circle(b/sigma)-h_2(q)，供正文比较。
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "result")
B = 1.0


def binary_entropy(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def entropy_bits(p: np.ndarray) -> float:
    q = np.asarray(p, dtype=float)
    q = q[q > 0.0]
    return float(-np.sum(q * np.log2(q)))


def symmetric_markov_entropy_rate(stay_probability: float) -> float:
    """平稳对称二状态 Markov 链的比特熵率。

    两个状态的平稳概率均为 1/2，转移矩阵对角元为 stay_probability。
    """
    s = float(stay_probability)
    if not 0.0 <= s <= 1.0:
        raise ValueError("stay_probability must lie in [0,1]")
    return binary_entropy(1.0 - s)


def circle_capacity(ratio: float) -> float:
    """模环加性标准高斯信道容量，ratio=period/sigma。

    均匀输入最优，容量为 log2(ratio)-h_2(Z mod ratio)。大 ratio 时
    折叠尾概率可忽略，直接使用未折叠高斯熵；中小 ratio 数值积分。
    """
    r = float(ratio)
    if r <= 0.0:
        return 0.0
    if r >= 18.0:
        return max(0.0, math.log2(r / math.sqrt(2.0 * math.pi * math.e)))

    ngrid = 8192
    z = (np.arange(ngrid) + 0.5) * (r / ngrid) - r / 2.0
    if r < 1.0:
        # Poisson 求和形式；小周期时只需极少 Fourier 项。
        density = np.ones_like(z) / r
        for n in range(1, 16):
            coef = math.exp(-2.0 * math.pi**2 * n * n / (r * r))
            if coef < 1e-16:
                break
            density += (2.0 / r) * coef * np.cos(2.0 * math.pi * n * z / r)
    else:
        # 直接叠加标准高斯副本。
        kmax = int(math.ceil(9.0 / r)) + 1
        density = np.zeros_like(z)
        for k in range(-kmax, kmax + 1):
            t = z + k * r
            density += np.exp(-0.5 * t * t) / math.sqrt(2.0 * math.pi)

    dz = r / ngrid
    density /= float(np.sum(density) * dz)
    h = float(-np.sum(density * np.log2(np.maximum(density, 1e-300))) * dz)
    return max(0.0, math.log2(r) - h)


@dataclass(frozen=True)
class Host:
    name: str
    cdf: Callable[[float], float]


HOSTS = [
    Host("Exp(mean=2B)", lambda x: 0.0 if x <= 0.0 else 1.0 - math.exp(-x / 2.0)),
    Host("U(0,2B)", lambda x: min(1.0, max(0.0, x / 2.0))),
    Host("Exp(mean=B)", lambda x: 0.0 if x <= 0.0 else 1.0 - math.exp(-x)),
    Host("Exp(mean=B/2)", lambda x: 0.0 if x <= 0.0 else 1.0 - math.exp(-2.0 * x)),
    Host("U(B/2,3B/2)", lambda x: min(1.0, max(0.0, x - 0.5))),
    Host("U(B,2B) (sparse)", lambda x: min(1.0, max(0.0, x - 1.0))),
]


def geometric_probabilities(cdf: Callable[[float], float], beta: float,
                            tail_tol: float = 1e-14) -> np.ndarray:
    if beta <= 1.0:
        raise ValueError("beta must exceed one")
    probs = [max(0.0, 1.0 - cdf(B))]
    for k in range(1, 100000):
        upper = B * beta ** (-(k - 1))
        lower = B * beta ** (-k)
        probs.append(max(0.0, cdf(upper) - cdf(lower)))
        if cdf(lower) <= tail_tol:
            break
    p = np.asarray(probs, dtype=float)
    p /= float(np.sum(p))
    return p


def partition_penalty(cdf: Callable[[float], float], beta: float) -> dict:
    p = geometric_probabilities(cdf, beta)
    k = np.arange(len(p), dtype=float)
    label = entropy_bits(p)
    resolution = float(np.sum(k * p) * math.log2(beta))
    return {
        "beta": float(beta),
        "label_entropy": label,
        "resolution_loss": resolution,
        "J_beta": label + resolution,
        "tiers": int(len(p)),
    }


def select_beta_bounded(cdf: Callable[[float], float]) -> dict:
    """在固定区间内数值选择分层比。

    返回的 J_beta_hat 是候选点 beta_hat 上的可证固定设计代价，
    不是理论全局下确界 Delta_plus=inf_{beta>1} J_beta 的数值认证。
    """
    if cdf(B) <= 1e-15:
        penalty = partition_penalty(cdf, 2.0)
    else:
        # beta 过近 1 会制造过多层，过大则量化损失变大。
        fit = minimize_scalar(
            lambda z: partition_penalty(cdf, math.exp(z))["J_beta"],
            bounds=(math.log(1.02), math.log(20.0)),
            method="bounded",
            options={"xatol": 1e-8},
        )
        penalty = partition_penalty(cdf, math.exp(float(fit.x)))
    return {
        "beta_hat": penalty["beta"],
        "label_entropy": penalty["label_entropy"],
        "resolution_loss": penalty["resolution_loss"],
        "J_beta_hat": penalty["J_beta"],
        "tiers": penalty["tiers"],
    }


def multires_rate(cdf: Callable[[float], float], gamma: float, beta: float) -> float:
    p = geometric_probabilities(cdf, beta)
    caps = np.array([circle_capacity(gamma * beta ** (-k)) for k in range(len(p))])
    return max(0.0, float(np.dot(p, caps)) - entropy_bits(p))


def single_threshold_rate(cdf: Callable[[float], float], gamma: float) -> tuple[float, float]:
    """同一状态移除引理下的最优单阈值界及其 b/B。"""
    # b/B 在 [1/gamma^2,1] 上搜索；极小周期的模环容量为零。
    lo = max(1e-12, gamma ** -2)

    def objective(logb: float) -> float:
        b = math.exp(logb)
        q = max(0.0, 1.0 - cdf(b))
        return -(q * circle_capacity(b * gamma) - binary_entropy(q))

    grid = np.linspace(math.log(lo), 0.0, 241)
    vals = np.array([objective(x) for x in grid])
    j = int(np.argmin(vals))
    left = grid[max(0, j - 1)]
    right = grid[min(len(grid) - 1, j + 1)]
    if right == left:
        b = math.exp(float(grid[j]))
        return max(0.0, -float(vals[j])), b
    fit = minimize_scalar(objective, bounds=(left, right), method="bounded",
                          options={"xatol": 1e-8})
    b = math.exp(float(fit.x))
    return max(0.0, -float(fit.fun)), b


def full_width_modulo_benchmark(gamma: float) -> float:
    """满宽度模环基准。

    该量与幅度受限高斯信道共享高信噪比体积渐近，但在有限噪声下
    不等同于 Smith 幅度约束容量。
    """
    return circle_capacity(gamma)


def residue_representative(lo: float, hi: float, period: float, residue: float) -> float:
    """返回闭区间 [lo,hi] 中同余于 residue (mod period) 的最小代表点。"""
    if period <= 0.0 or not (0.0 <= residue < period):
        raise ValueError("invalid period or residue")
    point = residue + math.ceil((lo - residue) / period - 1e-14) * period
    if point < lo - 1e-10 or point > hi + 1e-10:
        raise ValueError("interval does not contain the requested residue")
    return min(max(point, lo), hi)


def main() -> dict:
    rows = {
        "_meta": {
            "schema_version": 4,
            "paper_symbols": {
                "fixed_design_penalty": "J_beta",
                "reported_numerical_bound": "J_beta_hat",
                "theoretical_optimized_endpoint": "Delta_plus=inf_{beta>1} J_beta",
            },
            "beta_selection": {
                "method": "scipy.optimize.minimize_scalar bounded in log(beta)",
                "interval": [1.02, 20.0],
                "xatol": 1e-8,
                "global_optimum_certified": False,
            },
            "buffer_B": B,
            "gamma_definition": "B/sigma_xi",
            "rate_unit": "bits/packet",
            "full_width_curve": "modulo-circle benchmark, not exact finite-noise amplitude capacity",
        }
    }
    print("=" * 78)
    print("多分辨率模格常数间隙（B=1）")
    print("=" * 78)
    print(f"{'host':<23}{'beta_hat':>9}{'Hbar(K)':>10}{'分辨率损失':>13}{'J_beta_hat':>12}")
    for host in HOSTS:
        opt = select_beta_bounded(host.cdf)
        rows[host.name] = {"partition": {k: round(v, 6) if isinstance(v, float) else v
                                         for k, v in opt.items()}}
        print(f"{host.name:<23}{opt['beta_hat']:>9.3f}{opt['label_entropy']:>10.3f}"
              f"{opt['resolution_loss']:>13.3f}{opt['J_beta_hat']:>12.3f}")

    gammas = [2**k for k in range(3, 19)]
    for host in HOSTS:
        beta = rows[host.name]["partition"]["beta_hat"]
        curves = {}
        for gamma in gammas:
            multi = multires_rate(host.cdf, gamma, beta)
            single, threshold = single_threshold_rate(host.cdf, gamma)
            curves[str(gamma)] = {
                "full_width_modulo": round(full_width_modulo_benchmark(gamma), 6),
                "multires": round(multi, 6),
                "single_threshold": round(single, 6),
                "single_b_over_B": round(threshold, 8),
            }
        rows[host.name]["curves"] = curves

    os.makedirs(RES, exist_ok=True)
    out = os.path.join(RES, "multires.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print(f"\nwrote {os.path.relpath(out, os.path.join(HERE, '..'))}")
    return rows


if __name__ == "__main__":
    main()
