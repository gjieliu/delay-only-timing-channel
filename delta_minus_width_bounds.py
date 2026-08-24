# -*- coding: utf-8 -*-
"""命题 prop:geom-width 单字母窗口（§5）的生产脚本：
  Delta_minus_lo = E log2( 2/(2-(1-s)^2) ),
  Delta_minus_hi = E log2(1/s),  s=(d∧B)/B。

(1) 甲/乙/丙三族窗口用高精度 Simpson 求解析积分；
(2) 表 2（result/correlation_full.json）九个 Delta_minus 的窗口包含性与左端相对紧度；
(3) 六个附加 i.i.d. 边缘的抽查（窗口与中点规则 Delta_minus 直接估计）。
论文 §5 引用：窗口 [0.031,0.221]、[0.215,1.065]、[0.486,2.065]，
紧度 95%–98%、76%–82%、56%–59%，15 项包含性检查全部通过。
用法：python3 delta_minus_width_bounds.py           # 全量，写 result/delta_minus_width_bounds.json
      python3 delta_minus_width_bounds.py --quick   # 冒烟，写 result/delta_minus_width_bounds_quick.json
依赖：result/correlation_full.json（先跑 correlation_full.py；quick 模式亦可用 *_quick.json）
"""
import json
import math
import os
import sys

import numpy as np

B = 1.0


def lam(ds, N=2000, burn=3000):
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


def simpson_uniform(y, a, b):
    """均匀网格复合 Simpson（点数为奇数）。"""
    n = len(y) - 1
    assert n % 2 == 0
    h = (b - a) / n
    return h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


def delta_minus_lo_hi_uniform(a, b, M=4_000_001):
    """d ~ Unif(a,b) 的两个解析界（解析积分的高精度数值实现）。"""
    x = np.linspace(a, b, M)
    s = np.minimum(x, B) / B
    flo = np.log2(2.0 / (2.0 - (1.0 - s) ** 2))
    fhi = -np.log2(s)
    return (float(simpson_uniform(flo, a, b) / (b - a)),
            float(simpson_uniform(fhi, a, b) / (b - a)))


def delta_minus_lo_hi_from_samples(d):
    s = np.minimum(d, B) / B
    lo = float(np.mean(np.log2(2.0 / (2.0 - (1.0 - s) ** 2))))
    hi = float(np.mean(-np.log2(s)))
    return lo, hi


def main():
    quick = "--quick" in sys.argv
    M = 400_001 if quick else 4_000_001
    n_iid = 100_000 if quick else 400_000

    here = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.join(here, "result")
    corr_fn = os.path.join(res_dir, "correlation_full_quick.json" if quick
                           else "correlation_full.json")
    if not os.path.exists(corr_fn) and quick:
        corr_fn = os.path.join(res_dir, "correlation_full.json")
    if not os.path.exists(corr_fn):
        sys.exit("缺少 %s：请先运行 correlation_full.py" % os.path.normpath(corr_fn))
    with open(corr_fn, encoding="utf-8") as f:
        corr = json.load(f)

    FAM = {"jia": ("甲 Unif(0.5,1.5)", 0.5, 1.5),
           "yi": ("乙 Unif(0.25,0.75)", 0.25, 0.75),
           "bing": ("丙 Unif(0.125,0.375)", 0.125, 0.375)}

    checks_total, checks_ok = 0, 0
    fams_out = []
    print("=== 三族窗口（Simpson, M=%d）与表 2 九个 Delta_minus 的包含性 ===" % M, flush=True)
    for fam in corr["families"]:
        key = fam["key"]
        label, a, b = FAM[key]
        lo, hi = delta_minus_lo_hi_uniform(a, b, M=M)
        print("%-22s 窗口 [%.4f, %.4f]  窗宽 %.4f" % (label, lo, hi, hi - lo), flush=True)
        rows = []
        tights = []
        for row in fam["rows"]:
            r = row["Delta_minus_mean"]
            ok = bool(lo - 1e-3 <= r <= hi + 1e-3)
            checks_total += 1
            checks_ok += ok
            t = 100.0 * lo / r
            tights.append(t)
            rows.append({"theta": row["theta"], "Delta_minus": r, "contained": bool(ok),
                         "tightness_pct": t})
            print("   theta=%.2f  Delta_minus=%.4f  %s  紧度 %.1f%%"
                  % (row["theta"], r, "OK" if ok else "!!! 违反", t), flush=True)
        fams_out.append({"key": key, "label": label, "a": a, "b": b,
                         "Delta_minus_lo": lo, "Delta_minus_hi": hi,
                         "window_width": hi - lo,
                         "rows": rows,
                         "tightness_min": min(tights), "tightness_max": max(tights)})

    print("\n=== 附加 i.i.d. 边缘抽查（n=%d, N=1600, 中点规则）===" % n_iid, flush=True)
    rng = np.random.default_rng(7)
    cases = {
        "Unif(0,2B)": rng.uniform(0.0, 2.0, n_iid),
        "Unif(0,B)": rng.uniform(0.0, 1.0, n_iid),
        "两点 {0.2B,2B} 等概": np.where(rng.random(n_iid) < 0.5, 0.2, 2.0),
        "Exp(mean=B)": rng.exponential(1.0, n_iid),
        "Exp(mean=3B)": rng.exponential(3.0, n_iid),
        "Pareto a=2 (min .3B)": 0.3 / rng.random(n_iid) ** 0.5,
    }
    iid_out = []
    for name, d in cases.items():
        lo, hi = delta_minus_lo_hi_from_samples(d)
        r = float(-lam(d, N=1600))
        ok = bool(lo - 2e-3 <= r <= hi + 2e-3)
        checks_total += 1
        checks_ok += ok
        iid_out.append({"name": name, "Delta_minus_lo": lo, "Delta_minus": r,
                        "Delta_minus_hi": hi,
                        "contained": bool(ok)})
        print("%-22s  Delta_minus_lo=%.4f  Delta_minus=%.4f  Delta_minus_hi=%.4f   %s"
              % (name, lo, r, hi, "OK" if ok else "!!! 违反"), flush=True)

    out = {
        "meta": {"quick": quick, "schema_version": 3,
                 "simpson_points": M, "iid_path_length": n_iid,
                 "grid_points": 1600, "burn_in": 3000, "B": B,
                 "correlation_source": os.path.basename(corr_fn),
                 "script": "delta_minus_width_bounds.py",
                 "paper_symbol": "Delta_minus", "units": "bits/packet"},
        "families": fams_out,
        "iid_checks": iid_out,
        "checks_total": int(checks_total), "checks_ok": int(checks_ok),
        "all_contained": bool(checks_ok == checks_total),
    }
    os.makedirs(res_dir, exist_ok=True)
    fn = os.path.join(res_dir, "delta_minus_width_bounds_quick.json" if quick
                      else "delta_minus_width_bounds.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n包含性 %d/%d；已写 %s" % (checks_ok, checks_total, os.path.normpath(fn)),
          flush=True)


if __name__ == "__main__":
    main()
