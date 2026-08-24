# -*- coding: utf-8 -*-
"""表2（tab:correlation，相关宿主 9 行）的生产脚本。

甲/乙/丙三族：两状态对称马氏尺度层（持续概率 theta），条件均匀间隔；
乙 = 甲整体 ×1/2，丙 = 甲整体 ×1/4，B=1 不变。
走廊指数 Lambda 由正算子迭代（中点矩形规则）估计，Delta_minus = -Lambda（以 log2 B=0 为基准）。

正式参数（与论文表 2 一致）：N=1600 网格，n=2,000,000 步，种子 (11,12,13)，burn=3000；
网格自检：丙族 theta=0.5，n=200,000，seed=11，N∈{800,1600,3200}。
用法：python3 correlation_full.py            # 全量（约 30–60 分钟），写 result/correlation_full.json
      python3 correlation_full.py --quick    # 快速冒烟（单种子短路径），写 result/correlation_full_quick.json
"""
import json
import math
import os
import sys

import numpy as np

B = 1.0


def lam(ds, N=2000, burn=3000):
    """路径条件走廊指数 Lambda 的中点矩形离散：
    (T_d phi)(x) = int_0^{min(B, x+d)} phi(t) dt，phi_1 = 1，
    Lambda = lim (1/n) sum log2(归一化常数)。返回 (1/cnt)*sum log2 c。"""
    h = B / N
    x = (np.arange(N) + 0.5) * h
    phi = np.ones(N) / B
    tot, cnt = 0.0, 0
    for i, d in enumerate(ds):
        cum = np.empty(N + 1)
        cum[0] = 0.0
        np.cumsum(phi, out=cum[1:])
        cum[1:] *= h
        u = np.minimum(B, x + d)
        k = np.minimum((u / h).astype(np.int64), N - 1)
        F = cum[k] + phi[k] * (u - k * h)
        c = F.sum() * h
        phi = F / c
        if i >= burn:
            tot += np.log2(c)
            cnt += 1
    return tot / cnt


def markov2(sub, theta, n, seed=1):
    """两状态持续概率 theta 的对称马氏链，平稳分布 (1/2,1/2)；
    每个状态内 d 均匀落在各自子区间——单时刻边缘与 theta 无关。"""
    r = np.random.default_rng(seed)
    stay = r.random(n) < theta
    v = r.random(n)
    s = 0
    out = np.empty(n)
    for i in range(n):
        if i and not stay[i]:
            s ^= 1
        a, b = sub[s]
        out[i] = a + (b - a) * v[i]
    return out


def h2(p):
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def est(sub, th, n, N, seeds):
    v = [-lam(markov2(sub, th, n, seed=sd), N=N) for sd in seeds]
    return float(np.mean(v)), float(np.max(v) - np.min(v)), [float(x) for x in v]


JIA = [(1.0, 1.5), (0.5, 1.0)]
YI = [(0.5, 0.75), (0.25, 0.5)]
BING = [(0.25, 0.375), (0.125, 0.25)]

FAMILIES = [
    ("jia", "甲", JIA, 1.0, 0.5, 0.5),      # (键, 名, 子区间, 尺度, P(d>=B), E K log2 beta)
    ("yi", "乙", YI, 0.5, 0.0, 1.5),
    ("bing", "丙", BING, 0.25, 0.0, 2.5),
]
THETAS = (0.50, 0.90, 0.99)


def main():
    quick = "--quick" in sys.argv
    if quick:
        n_main, seeds, grid_n = 200_000, (11,), 50_000
        path_lens = (50_000, 100_000, 200_000)
    else:
        n_main, seeds, grid_n = 2_000_000, (11, 12, 13), 200_000
        path_lens = (500_000, 1_000_000, 2_000_000)
    N = 1600

    grid_checks = []
    for th_g in (0.5, 0.99):
        print("--- 网格自检（丙, theta=%.2f, n=%d, seed=11）---" % (th_g, grid_n),
              flush=True)
        gc = {}
        for Ng in (800, 1600, 3200):
            m, _, _ = est(BING, th_g, grid_n, Ng, (11,))
            gc[str(Ng)] = m
            print("   N=%-5d Delta_minus=%.4f" % (Ng, m), flush=True)
        grid_checks.append({"family": "bing", "theta": th_g, "path_length": grid_n,
                            "seed": 11, "Delta_minus_by_grid": gc})

    print("--- 路径长度自检（丙, theta=0.99, N=%d, seed=11）---" % N, flush=True)
    path_check = {"family": "bing", "theta": 0.99, "grid_points": N, "seed": 11,
                  "Delta_minus_by_length": {}}
    for nl in path_lens:
        m, _, _ = est(BING, 0.99, nl, N, (11,))
        path_check["Delta_minus_by_length"][str(nl)] = m
        print("   n=%-8d Delta_minus=%.4f" % (nl, m), flush=True)

    print("\n--- %s (N=%d, n=%s, seeds=%s) ---"
          % ("快速" if quick else "正式", N, n_main, seeds), flush=True)
    fams_out = []
    max_spread = 0.0
    for key, name, sub, scale, p_full, ek in FAMILIES:
        print("\n%s  尺度=%s  P(d>=B)=%.1f  E[K]log2b=%.3f" % (name, scale, p_full, ek), flush=True)
        rows = []
        for th in THETAS:
            m, sp, per_seed = est(sub, th, n_main, N, seeds)
            H = h2(th)
            j_beta = ek + H
            width = j_beta - m
            max_spread = max(max_spread, sp)
            rows.append({
                "theta": th, "H_rate": H, "EK_log2beta": ek, "J_beta": j_beta,
                "Delta_minus_mean": m, "Delta_minus_spread": sp,
                "Delta_minus_by_seed": per_seed,
                "J_beta_minus_Delta_minus": width,
                "tex_row": "$%.2f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.4f$ & $%.3f$"
                           % (th, H, ek, j_beta, m, width),
            })
            print("   theta=%.2f  Hbar=%.3f  J_beta=%.3f  Delta_minus=%.4f  (spread %.4f)  width %.3f"
                  % (th, H, j_beta, m, sp, width), flush=True)
        # 配对差（同种子同底层随机流，theta=0.99 与 0.50 之差）及配对标准误
        d99, d50 = rows[-1]["Delta_minus_by_seed"], rows[0]["Delta_minus_by_seed"]
        diffs = [a - b for a, b in zip(d99, d50)]
        paired_se = (float(np.std(diffs, ddof=1)) / math.sqrt(len(diffs))
                     if len(diffs) > 1 else 0.0)
        fams_out.append({
            "key": key, "name": name, "sub_intervals": sub, "scale": scale,
            "P_d_geq_B": p_full, "EK_log2beta": ek, "rows": rows,
            "Delta_minus_change": rows[-1]["Delta_minus_mean"] - rows[0]["Delta_minus_mean"],
            "Delta_minus_change_by_seed": diffs,
            "Delta_minus_change_paired_se": paired_se,
        })
        print("   Δ_-(0.99−0.50) 差值=%.4f  配对标准误=%.5f" %
              (rows[-1]["Delta_minus_mean"] - rows[0]["Delta_minus_mean"], paired_se), flush=True)

    out = {
        "meta": {
            "quick": quick, "schema_version": 3, "B": B, "beta": 2, "grid_points": N,
            "path_length": n_main, "repetitions": len(seeds), "seeds": list(seeds),
            "burn_in": 3000, "thetas": list(THETAS),
            "rule": "midpoint", "script": "correlation_full.py",
            "paper_symbols": {
                "fixed_design_penalty": "J_beta",
                "corridor_lower_endpoint": "Delta_minus",
            },
            "units": "bits/packet",
        },
        "grid_checks": grid_checks,
        "path_check": path_check,
        "families": fams_out,
        "max_spread": max_spread,
    }
    res_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result")
    os.makedirs(res_dir, exist_ok=True)
    fn = os.path.join(res_dir, "correlation_full_quick.json" if quick else "correlation_full.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n最大实现间极差 %.4f；已写 %s" % (max_spread, os.path.normpath(fn)), flush=True)


if __name__ == "__main__":
    main()
