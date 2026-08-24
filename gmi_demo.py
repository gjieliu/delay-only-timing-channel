# -*- coding: utf-8 -*-
"""命题"小噪声损失的显式常数上界"中常数 c1 = E log2 m + H(X1|U1) − h(X1) 的数值估计。

三个设计 × 3 种子：
  Unif(B/2,3B/2)（甲边缘，beta=2，层 {0,1}，E K log2 beta = 0.5）；
  Unif(B/4,3B/4)（乙边缘，beta=2，层 {1,2}，E K log2 beta = 1.5）；
  Unif(B/4,3B/4)（乙边缘，beta=4，单层 {1}，E K log2 beta = 2）。

方法学要点（认证化设计）：
* 层端点取精确值：d~Unif[lo,hi) 半开抽样，胞元对 [lo,hi) 解析积分，无端点扰动；
* 平稳滞留池取单一长轨迹的确定性等距稀释（无子采样随机性），对全部评估点共用；
* 评估用策略字母 u 由独立随机流新鲜抽取，与轨迹解耦；
* H(X1|U1) 的标准误由 u 间波动直接给出；池规模、密度网格、评估点数、
  密度池的半倍/双倍扰动为共用 u 流的配对比较，隔离各近似环节的移动量；
* c1 为 Monte Carlo 估计，报告种子间极差与配对自检移动量，不构成确定性认证区间。

用法：python3 gmi_demo.py           # 全量（约 15 分钟），写 result/gmi_demo.json
      python3 gmi_demo.py --quick   # 快速冒烟，写 result/gmi_demo_quick.json
"""
import json
import math
import os
import sys

import numpy as np


def cells(delta, v, m, lo, hi):
    out = []
    b1, b2 = lo, min(hi, max(delta, lo))
    if b2 > b1:
        out.append((b2 - b1, 0))
    glo = max(lo, delta)
    if hi > glo:
        a0, a1 = glo - delta, hi - delta
        for l in range(max(math.ceil((a0 - v) / m), 0), math.ceil((a1 - v) / m) + 1):
            l_, h_ = max(a0, v + (l - 1) * m), min(a1, v + l * m)
            if h_ > l_:
                out.append((h_ - l_, max(l, 0)))
    return out


def thin(arr, size):
    """确定性等距稀释：无放回、无随机性，覆盖整段轨迹。"""
    idx = np.linspace(0, len(arr) - 1, size).astype(np.int64)
    return arr[np.unique(idx)]


def run(tiers, dlo_all, dhi_all, seed, n=300_000, n_ev=1800, n_pool=4400,
        n_fx=5600, n_grid=4400):
    def simulate(n):
        r = np.random.default_rng(seed)
        d = r.uniform(dlo_all, dhi_all, n)
        vs = {k: r.random(n) * m for k, (m, lo, hi) in tiers.items()}
        dl = np.empty(n + 1)
        dl[0] = 0.0
        for i in range(n):
            k = [kk for kk, (m, lo, hi) in tiers.items() if lo <= d[i] < hi][0]
            m = tiers[k][0]
            v = vs[k][i]
            a = max(d[i] - dl[i], 0.0)
            dl[i + 1] = dl[i] + v + math.ceil((a - v) / m) * m - d[i]
        return dl

    DL = simulate(n)
    burn = 2000
    pool_H = thin(DL[burn:n], n_pool)     # H(X|U) 边缘化用固定池
    pool_fx = thin(DL[burn:n], n_fx)      # h(X) 密度混合用固定池

    # H(X1|U1)：u 由独立流抽取；对固定池作解析胞元枚举
    rng_u = np.random.default_rng(seed + 2000)
    Hs = []
    for _ in range(n_ev):
        u = {k: rng_u.random() * m for k, (m, lo, hi) in tiers.items()}
        agg = {}
        for dl_ in pool_H:
            for k, (m, lo, hi) in tiers.items():
                for w, l in cells(dl_, u[k], m, lo, hi):
                    agg[(k, l)] = agg.get((k, l), 0.0) + w
        Zt = sum(agg.values())
        Hs.append(-sum(w / Zt * math.log2(w / Zt) for w in agg.values() if w > 1e-15))
    H_XU = float(np.mean(Hs))
    H_se = float(np.std(Hs, ddof=1) / math.sqrt(len(Hs)))

    # h(X1)：解析条件密度对固定池平均
    xg = np.linspace(0, dhi_all + 1.05, n_grid)
    dx = xg[1] - xg[0]

    def fx(delta):
        f = np.zeros_like(xg)
        for k, (m, lo, hi) in tiers.items():
            p = (hi - lo) / (dhi_all - dlo_all)
            pb = max(min(hi, delta) - lo, 0.0) / (hi - lo)
            if pb > 0:
                f += p * pb * ((xg >= 0) & (xg < m)) / m
            glo, ghi = max(0.0, lo - delta), hi - delta
            if ghi > glo:
                f += p * (1 - pb) * np.clip(
                    (np.minimum(xg, ghi) - np.maximum(xg - m, glo)) / (ghi - glo), 0, 1) / m
        return f

    fs = np.zeros_like(xg)
    for dl_ in pool_fx:
        fs += fx(dl_)
    fs /= len(pool_fx)
    mk = fs > 1e-14
    h_X = float(-np.sum(fs[mk] * np.log2(fs[mk])) * dx)
    ElogM = sum((hi - lo) / (dhi_all - dlo_all) * math.log2(m) for m, lo, hi in tiers.values())
    return ElogM + H_XU - h_X, H_XU, h_X, ElogM, H_se


# 层端点精确（旧版 1.5001/0.7501 已弃用）
FAMS = {
    "Unif(B/2,3B/2), beta=2": (
        {0: (1.0, 1.0, 1.5), 1: (0.5, 0.5, 1.0)}, 0.5, 1.5, 2, 0.5, 1.0),
    "Unif(B/4,3B/4), beta=2": (
        {1: (0.5, 0.5, 0.75), 2: (0.25, 0.25, 0.5)}, 0.25, 0.75, 2, 1.5, 1.0),
    "Unif(B/4,3B/4), beta=4": (
        {1: (0.25, 0.25, 0.75)}, 0.25, 0.75, 4, 2.0, 0.0),
}


def main():
    quick = "--quick" in sys.argv
    if quick:
        seeds, kw = (31,), dict(n=100_000, n_ev=300, n_pool=1200, n_fx=1500)
    else:
        seeds, kw = (31, 32, 33), dict(n=300_000, n_ev=1800, n_pool=4400, n_fx=5600)

    fams_out = {}
    for name, (tiers, lo, hi, beta, ek_log2b, h_layer) in FAMS.items():
        cs, per_seed = [], []
        for seed in seeds:
            c1, HXU, hX, ElogM, Hse = run(tiers, lo, hi, seed, **kw)
            cs.append(c1)
            per_seed.append({"seed": seed, "H_XU": HXU, "H_XU_se": Hse,
                             "h_X": hX, "Elog2m": ElogM, "c1": c1})
            print("%-15s seed=%d  H(X|U)=%.4f (se %.4f)  h(X)=%.4f  c1=%.4f"
                  % (name, seed, HXU, Hse, hX, c1), flush=True)
        mean = float(np.mean(cs))
        rng_ = float(max(cs) - min(cs))
        fams_out[name] = {
            "tiers": {str(k): list(v) for k, v in tiers.items()},
            "beta": beta, "EK_log2beta": ek_log2b, "H_layer": h_layer,
            "J_beta": ek_log2b + h_layer,
            "per_seed": per_seed, "c1_mean": mean, "c1_range": rng_,
            "tightened_bound": ek_log2b + mean,
        }
        print("%-15s c1 均值=%.4f  极差=%.4f  上界=%.3f"
              % (name, mean, rng_, ek_log2b + mean), flush=True)

    # 配对收敛自检（正式档）：甲边缘 seed=31，共用同一 u 随机流，
    # 只改动池规模/网格/评估点数/密度池，移动量即各近似环节的贡献。
    conv = None
    if not quick:
        tiers_a, lo_a, hi_a, _, _, _ = FAMS["Unif(B/2,3B/2), beta=2"]
        base = run(tiers_a, lo_a, hi_a, 31, **kw)[0]
        variants = {
            "pool_half": dict(kw, n_pool=kw["n_pool"] // 2),
            "pool_double": dict(kw, n_pool=kw["n_pool"] * 2),
            "grid_half": dict(kw, n_grid=2200),
            "grid_double": dict(kw, n_grid=8800),
            "ev_half": dict(kw, n_ev=kw["n_ev"] // 2),
            "fx_half": dict(kw, n_fx=kw["n_fx"] // 2),
            "fx_double": dict(kw, n_fx=kw["n_fx"] * 2),
        }
        vals = {}
        for vname, kv in variants.items():
            vals[vname] = run(tiers_a, lo_a, hi_a, 31, **kv)[0]
            print("配对自检 %-12s c1=%.4f (基准 %.4f, Δ=%+.5f)"
                  % (vname, vals[vname], base, vals[vname] - base), flush=True)
        conv = {"base_seed": 31, "base_c1": base, "variants": vals,
                "max_shift": max(abs(v - base) for v in vals.values())}
        print("配对自检最大移动 %.5f" % conv["max_shift"], flush=True)

    out = {
        "meta": {"quick": quick, "schema_version": 3,
                 "paper_symbol": "J_beta",
                 "seeds": list(seeds), "burn_in": 2000,
                 "B": 1.0, "design_count": len(FAMS), "n_grid": 4400,
                 "pooling": "deterministic-thin",
                 "script": "gmi_demo.py", **kw},
        "families": fams_out,
        "convergence_check": conv,
    }
    res_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result")
    os.makedirs(res_dir, exist_ok=True)
    fn = os.path.join(res_dir, "gmi_demo_quick.json" if quick else "gmi_demo.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("已写 %s" % os.path.normpath(fn), flush=True)


if __name__ == "__main__":
    main()
