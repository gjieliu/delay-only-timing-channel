# -*- coding: utf-8 -*-
"""表 tab:reset-control（丁族）的生产脚本。

条件于 d<B 时，宿主边缘恒为 Unif(B/4,3B/4)，关联由对称二状态
持续链承载；独立地以概率 p 插入 Unif(B,3B/2) 满宽度步。同一
(p, seed) 下两个 theta 共用均匀分位、满宽度位置与状态转移随机流。

用法：python3 run_ding.py           # 正式档，写 result/reset_control.json
      python3 run_ding.py --quick   # 冒烟档，写 result/reset_control_quick.json
"""
import json
import math
import os
import sys

import numpy as np

import corridor as C

B = 1.0
PS = (0.0, 0.10, 0.25, 0.50)
THETAS = (0.50, 0.99)


def gen_ding(theta, p, n, seed):
    """生成丁族路径；随机流只由 (p, seed) 决定，使不同 theta 可配对。"""
    ru = np.random.default_rng(1_000 + seed)
    rz = np.random.default_rng(2_000 + seed)
    rs = np.random.default_rng(3_000 + seed)
    u = ru.random(n)
    full = rz.random(n) < p
    stay_u = rs.random(n)
    state = 0
    out = np.empty(n)
    for i in range(n):
        if i and stay_u[i] >= theta:
            state ^= 1
        if full[i]:
            out[i] = 1.00 + 0.50 * u[i]
        elif state == 0:
            out[i] = 0.25 + 0.25 * u[i]
        else:
            out[i] = 0.50 + 0.25 * u[i]
    return out


def width_window_uniform(a=0.25, b=0.75, order=256):
    """Unif(a,b) 的单字母窗口；Gauss--Legendre 求积。"""
    x, w = np.polynomial.legendre.leggauss(order)
    d = (b - a) * x / 2.0 + (a + b) / 2.0
    s = np.minimum(d, B) / B
    lo = np.sum(w * np.log2(2.0 / (2.0 - (1.0 - s) ** 2))) / 2.0
    hi = np.sum(w * (-np.log2(s))) / 2.0
    return float(lo), float(hi)


def mean_se(values):
    mean = float(np.mean(values))
    se = (float(np.std(values, ddof=1) / math.sqrt(len(values)))
          if len(values) > 1 else 0.0)
    return mean, se


def main():
    quick = "--quick" in sys.argv
    if quick:
        grid_points, path_length, seeds = 300, 20_000, (1,)
    else:
        grid_points, path_length, seeds = 1200, 400_000, (1, 2, 3, 4)

    base_lo, base_hi = width_window_uniform()
    rows = []
    max_se = 0.0
    print("--- %s丁族：N=%d, n=%d, seeds=%s ---"
          % ("快速" if quick else "正式", grid_points, path_length, seeds),
          flush=True)
    for p in PS:
        by_theta = {}
        for theta in THETAS:
            values = [float(C.delta_minus(gen_ding(theta, p, path_length, seed),
                                       N=grid_points)) for seed in seeds]
            mean, se = mean_se(values)
            max_se = max(max_se, se)
            by_theta[str(theta)] = {"mean": mean, "se": se, "by_seed": values}
            print("p=%.2f theta=%.2f  Delta_-=%.4f  se=%.5f"
                  % (p, theta, mean, se), flush=True)

        low = by_theta[str(THETAS[0])]
        high = by_theta[str(THETAS[1])]
        paired = [b - a for a, b in zip(low["by_seed"], high["by_seed"])]
        sensitivity, sensitivity_se = mean_se(paired)
        rows.append({
            "p": p,
            "theta": by_theta,
            "sensitivity": sensitivity,
            "sensitivity_se": sensitivity_se,
            "mass_dilution_prediction": None,
            "ratio": None,
            "window_lo": (1.0 - p) * base_lo,
            "window_hi": (1.0 - p) * base_hi,
        })

    s0 = rows[0]["sensitivity"]
    for row in rows:
        pred = s0 * (1.0 - row["p"])
        row["mass_dilution_prediction"] = pred
        row["ratio"] = row["sensitivity"] / pred

    out = {
        "meta": {
            "quick": quick,
            "schema_version": 2,
            "script": "run_ding.py",
            "paper_symbol": "Delta_minus",
            "B": B,
            "grid_points": grid_points,
            "path_length": path_length,
            "repetitions": len(seeds),
            "seeds": list(seeds),
            "rule": "midpoint",
            "p_values": list(PS),
            "theta_values": list(THETAS),
            "sub_buffer_conditional_law": "Unif(B/4,3B/4)",
            "full_width_law": "Unif(B,3B/2)",
            "paired_random_streams_across_theta": True,
            "units": "bits/packet",
        },
        "base_window": {"lo": base_lo, "hi": base_hi},
        "rows": rows,
        "max_standard_error": max_se,
    }
    res_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result")
    os.makedirs(res_dir, exist_ok=True)
    filename = "reset_control_quick.json" if quick else "reset_control.json"
    path = os.path.join(res_dir, filename)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("已写 %s；最大标准误 %.5f" % (os.path.normpath(path), max_se), flush=True)


if __name__ == "__main__":
    main()
