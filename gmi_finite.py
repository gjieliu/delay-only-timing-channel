# -*- coding: utf-8 -*-
"""命题"小噪声损失的显式常数上界"的有限噪声数值检查。

对甲边缘 Unif(B/2,3B/2)（beta=2，层 {0,1}）在若干有限噪声水平 sigma 下，
按命题原样取匹配设计：p_r(y|u) = E_{delta~mu, d}[phi_sigma(y - X)]，
alpha0 = sigma/(B+sigma)，Lambda0 = 8 (x_max/sigma)^2 log2(e)，
以 Monte Carlo 估计已证下界 E_mu[Z]，并报告相应损失上界的数值估计
G_est(sigma) = log2(B/(sqrt(2 pi e) sigma)) - E_mu[Z]，
与零噪声显式上界常数 E K log2 beta + c1（gmi_demo.py 的估计）对照。

方法学要点：层端点取精确值（d~Unif[lo,hi) 半开抽样 + [lo,hi) 解析积分，无端点扰动）；
评估滞留取自独立轨迹（seed+7000），与参考密度所用平稳池（seed 轨迹）解耦；正式档
另对参考密度的平稳池与混合样本规模作半倍/双倍扰动，报告嵌套近似的移动量。
所报标准误仅为外层抽样波动，整组数字是命题下界的 Monte Carlo 估计，
不构成确定性认证区间。
用法：python3 gmi_finite.py           # 全量（约 20–30 分钟），写 result/gmi_finite.json
      python3 gmi_finite.py --quick   # 快速冒烟，写 result/gmi_finite_quick.json
依赖：无（自含）。
"""
import json
import math
import os
import sys

import numpy as np

B = 1.0
TIERS = {0: (1.0, 1.0, 1.5), 1: (0.5, 0.5, 1.0)}   # k: (m, dlo, dhi)，端点精确
DLO, DHI = 0.5, 1.5
EK_LOG2B = 0.5
XMAX = DHI + B     # 输入支撑上界 d_max + B


def simulate(n, seed):
    r = np.random.default_rng(seed)
    d = r.uniform(DLO, DHI, n)
    vs = {k: r.random(n) * m for k, (m, lo, hi) in TIERS.items()}
    dl = np.empty(n + 1)
    dl[0] = 0.0
    for i in range(n):
        k = 0 if d[i] >= 1.0 else 1
        m = TIERS[k][0]
        v = vs[k][i]
        a = max(d[i] - dl[i], 0.0)
        dl[i + 1] = dl[i] + v + math.ceil((a - v) / m) * m - d[i]
    return dl


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


def atom_law(u, pool_dl):
    """给定 u=(v0,v1)，对滞留池与 d 解析平均的原子集合与权重。"""
    agg = {}
    for dl_ in pool_dl:
        for k, (m, lo, hi) in TIERS.items():
            for w, l in cells(dl_, u[k], m, lo, hi):
                agg[(k, l)] = agg.get((k, l), 0.0) + w
    tot = sum(agg.values())
    xs = np.array([u[k] + l * TIERS[k][0] for (k, l) in agg])
    ws = np.array([w / tot for w in agg.values()])
    return xs, ws


def gauss(y, xs, sig):
    return np.exp(-((y - xs) ** 2) / (2 * sig * sig)) / (math.sqrt(2 * math.pi) * sig)


def evaluate(sig, pool_marg, pool_eval, rng, n_pool, n_ubar, n_ev):
    """按命题设计估计 E_mu[Z]：参考密度用 pool_marg，评估滞留用 pool_eval。"""
    ubar_atoms = []
    for _ in range(n_ubar):
        u = {k: rng.random() * TIERS[k][0] for k in TIERS}
        ubar_atoms.append(atom_law(u, pool_marg[rng.integers(0, len(pool_marg),
                                                             n_pool)]))
    bar_xs = np.concatenate([xs for xs, _ in ubar_atoms])
    bar_ws = np.concatenate([ws for _, ws in ubar_atoms]) / n_ubar

    alpha0 = sig / (B + sig)
    lam0 = 8.0 * (XMAX / sig) ** 2 * math.log2(math.e)
    Zs = []
    for _ in range(n_ev):
        dl0 = pool_eval[rng.integers(0, len(pool_eval))]
        u = {k: rng.random() * TIERS[k][0] for k in TIERS}
        dd = rng.uniform(DLO, DHI)
        k = 0 if dd >= 1.0 else 1
        m = TIERS[k][0]
        a = max(dd - dl0, 0.0)
        x = u[k] + math.ceil((a - u[k]) / m) * m
        y = x + sig * rng.standard_normal()
        xs, ws = atom_law(u, pool_marg[rng.integers(0, len(pool_marg), n_pool)])
        pr = float(np.dot(ws, gauss(y, xs, sig)))
        pbar = float(np.dot(bar_ws, gauss(y, bar_xs, sig)))
        ratio = ((1 - alpha0) * pr + alpha0 * pbar) / pbar
        Zs.append(min(math.log2(ratio), lam0))
    Zs = np.array(Zs)
    return float(Zs.mean()), float(Zs.std(ddof=1) / math.sqrt(len(Zs))), alpha0, lam0


def main():
    quick = "--quick" in sys.argv
    if quick:
        n_path, n_pool, n_ubar, n_ev, sigmas, seeds = (
            60_000, 400, 150, 1_500, (0.02,), (41,))
    else:
        n_path, n_pool, n_ubar, n_ev, sigmas, seeds = (
            300_000, 1_200, 600, 12_000, (0.04, 0.02, 0.01), (41, 42))
    burn = 2000

    out_pts = []
    pools = {}
    for seed in seeds:
        pool_marg = simulate(n_path, seed)[burn:]
        pool_eval = simulate(n_path, seed + 7000)[burn:]
        pools[seed] = (pool_marg, pool_eval)
        for sig in sigmas:
            rng = np.random.default_rng([seed, int(round(sig * 1000))])
            gmi, se, alpha0, lam0 = evaluate(sig, pool_marg, pool_eval, rng,
                                             n_pool, n_ubar, n_ev)
            A = math.log2(B / (math.sqrt(2 * math.pi * math.e) * sig))
            pt = {"seed": seed, "sigma": sig, "alpha0": alpha0, "Lambda0": lam0,
                  "A_sigma": A, "gmi": gmi, "gmi_se": se,
                  "G_bound": A - gmi, "G_bound_se": se,
                  "n_eval": n_ev, "n_pool": n_pool, "n_ubar": n_ubar}
            out_pts.append(pt)
            print("seed=%d sigma=%.3f  A=%.4f  E[Z]=%.4f (se %.4f)  "
                  "G_est=%.4f  对照零噪声上界常数 EKlog2b+c1" %
                  (seed, sig, A, gmi, se, A - gmi), flush=True)

    # 嵌套误差自检（正式档）：seed=41、sigma=0.02，参考密度平稳池与混合样本半倍/双倍
    nest = None
    if not quick:
        base = next(p["G_bound"] for p in out_pts
                    if p["seed"] == 41 and p["sigma"] == 0.02)
        pool_marg, pool_eval = pools[41]
        vals = {}
        for tag, (np_, nu_) in (("half", (n_pool // 2, n_ubar // 2)),
                                ("double", (n_pool * 2, n_ubar * 2))):
            rng = np.random.default_rng([41, 20, np_, nu_])
            gmi, se, _, _ = evaluate(0.02, pool_marg, pool_eval, rng,
                                     np_, nu_, 6000)
            A = math.log2(B / (math.sqrt(2 * math.pi * math.e) * 0.02))
            vals[tag] = {"G_bound": A - gmi, "se": se,
                         "n_pool": np_, "n_ubar": nu_}
            print("嵌套自检 %-6s G_est=%.4f (se %.4f)" % (tag, A - gmi, se),
                  flush=True)
        nest = {"seed": 41, "sigma": 0.02, "base_G": base, "variants": vals,
                "max_shift": max(abs(v["G_bound"] - base) for v in vals.values())}
        print("嵌套自检最大移动 %.4f" % nest["max_shift"], flush=True)

    out = {
        "meta": {"quick": quick, "family": "Unif(B/2,3B/2)", "beta": 2, "B": B,
                 "x_max": XMAX, "EK_log2beta": EK_LOG2B,
                 "path_length": n_path, "burn_in": burn, "seeds": list(seeds),
                 "eval_path_seed_offset": 7000,
                 "sigmas": list(sigmas), "script": "gmi_finite.py"},
        "points": out_pts,
        "nesting_check": nest,
    }
    res_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result")
    os.makedirs(res_dir, exist_ok=True)
    fn = os.path.join(res_dir, "gmi_finite_quick.json" if quick else "gmi_finite.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("已写 %s" % os.path.normpath(fn), flush=True)


if __name__ == "__main__":
    main()
