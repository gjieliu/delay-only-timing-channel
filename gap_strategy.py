"""复现 tab:gap —— 宿主未知与已知容量之差 Delta_b（状态增广 Shannon 策略，记忆松弛）。

因果状态在编码端已知时，无宿主边信息的可达是 Shannon-1958 策略容量
  C_ng = max_{p(u), x=f(u,s)} I(U;Y),  U ⊥ S   （非非因果 GP 的 I(U;Y)-I(U;S)=genie）。
宿主已知（genie）= sum_s w_s A(区间_s, sigma)。差 = genie - nogenie。

三种诊断（均为同心区间：中心 0，最大重叠 = 对译码端最坏），精确 Blahut-Arimoto、细网格：
  纯偏移对照：宽度恒定、仅中心起伏 -> 差应 ->0（偏移可由因果模格抵消，方法校验）。
  纯宽度对照：同心宽度 {2B,0.8B} 等概 -> 差 -> 正常数（反射壁宽度不可抵消）。
  Exp(均值 2B)：由 delta~U[0,B]（高信噪比达容量输入的滞留边缘）与 d~Exp(2B) 导出的可达宽度
                L=B-(delta-d)^+ 的分布，离散为同心宽度混合 -> 差（量级同 rho）。
稀疏宿主差恒为 0（可达宽度恒为 B，无宽度混淆；见定理 1）。
"""
import itertools
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "result")
B = 1.0


def ba(P, tol=1e-8, maxit=400):
    logP = np.log(np.clip(P, 1e-300, None))
    q = np.full(P.shape[0], 1.0 / P.shape[0])
    C = 0.0
    for _ in range(maxit):
        r = q @ P
        logr = np.log(np.clip(r, 1e-300, None))
        D = (P * (logP - logr[None, :])).sum(axis=1)
        m = D.max()
        w = q * np.exp(D - m)
        Z = w.sum()
        Cn = (m + math.log(Z)) / math.log(2)
        if abs(Cn - C) < tol:
            C = Cn
            break
        C, q = Cn, w / Z
    return C


def gauss_rows(xvals, y, sigma):
    W = np.exp(-0.5 * ((y[None, :] - np.asarray(xvals)[:, None]) / sigma) ** 2)
    return W / W.sum(axis=1, keepdims=True)


def _grids(states, sigma, kcap):
    dx = sigma / 2.5
    return [min(kcap, max(6, int(round((hi - lo) / dx)) + 1)) for (_w, lo, hi) in states]


def _ygrid(states, sigma, cap=900):
    lo = min(l for (_w, l, _h) in states) - 0.6
    hi = max(h for (_w, _l, h) in states) + 0.6
    ny = min(cap, int((hi - lo) / (sigma / 4)) + 1)
    return np.linspace(lo, hi, ny)


def genie(states, sigma, kcap=110):
    y = _ygrid(states, sigma)
    ks = _grids(states, sigma, kcap)
    tw = sum(w for (w, _l, _h) in states)
    c = 0.0
    for (w, lo, hi), k in zip(states, ks):
        c += (w / tw) * ba(gauss_rows(np.linspace(lo, hi, k), y, sigma))
    return c


def nogenie(states, sigma, kcap=110, combo_cap=18000):
    y = _ygrid(states, sigma)
    ks = _grids(states, sigma, kcap)
    # 控制组合总数（kx^S），必要时压 kcap
    while np.prod(ks) > combo_cap and max(ks) > 6:
        ks = [max(6, k - 1) for k in ks]
    tw = sum(w for (w, _l, _h) in states)
    ws = [w / tw for (w, _l, _h) in states]
    Wl = [gauss_rows(np.linspace(lo, hi, k), y, sigma) for (w, lo, hi), k in zip(states, ks)]
    combos = list(itertools.product(*[range(k) for k in ks]))
    P = np.empty((len(combos), len(y)))
    for ci, ch in enumerate(combos):
        acc = np.zeros(len(y))
        for s, c in enumerate(ch):
            acc += ws[s] * Wl[s][c]
        P[ci] = acc
    return ba(P)


def concentric(widths, weights):
    return [(w, -L / 2.0, L / 2.0) for L, w in zip(widths, weights)]


def gap(states, gamma, kcap=110):
    sigma = B / gamma
    return genie(states, sigma, kcap) - nogenie(states, sigma, kcap)


# ---- Exp(均值 2B) 可达宽度混合 --------------------------------------------
def exp_width_states(mean=2.0, nbin=4, Nmc=300000, seed=0):
    rng = np.random.default_rng(seed)
    delta = rng.uniform(0.0, B, Nmc)          # 高信噪比达容量输入的滞留边缘 ~ U[0,B]
    d = rng.exponential(mean, Nmc)
    L = B - np.maximum(delta - d, 0.0)        # 可达宽度 in (0,B]
    edges = np.linspace(L.min(), B, nbin + 1)
    widths, weights = [], []
    for k in range(nbin):
        hi_incl = (L <= edges[k + 1]) if k == nbin - 1 else (L < edges[k + 1])
        mask = (L >= edges[k]) & hi_incl
        w = float(mask.mean())
        if w > 5e-3:
            widths.append(float(0.5 * (edges[k] + edges[k + 1])))
            weights.append(w)
    return concentric(widths, weights)


GAMMAS = [8, 16, 32, 48]


def main():
    print("=" * 60)
    print("tab:gap 复现 —— Delta_b = genie - nogenie（状态增广 Shannon 策略）")
    print("=" * 60)

    offset = [(0.5, -0.7, 0.7), (0.5, 0.3, 1.7)]     # 纯偏移（宽度恒 1.4）
    width = concentric([2.0, 0.8], [0.5, 0.5])        # 纯宽度（同心）
    expst = exp_width_states(mean=2.0)
    print("Exp(2B) 可达宽度混合状态：", [(round(w, 3), round(hi - lo, 3)) for (w, lo, hi) in expst])

    rows = {"pure_offset": {}, "exp2B": {}, "pure_width": {}}
    print(f"\n{'gamma':>6}{'纯偏移':>10}{'Exp(2B)':>10}{'纯宽度':>10}")
    for g in GAMMAS:
        go = gap(offset, g)
        ge = gap(expst, g, kcap=42)
        gw = gap(width, g)
        rows["pure_offset"][g] = round(go, 3)
        rows["exp2B"][g] = round(ge, 3)
        rows["pure_width"][g] = round(gw, 3)
        print(f"{g:>6}{go:>10.3f}{ge:>10.3f}{gw:>10.3f}", flush=True)

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "gap.json"), "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print("\nwrote result/gap.json")
    return rows


if __name__ == "__main__":
    main()
