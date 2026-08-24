"""复现表 tab:penalty —— 反射走廊对数体积指数 Lambda(P_process,B) 与损失下端 Delta_minus。

固定起点走廊为
    A_n(d^n) = { u : delta_0=0, delta_i=sum_{l<=i}u_l,
                (delta_{i-1}-d_i)^+ <= delta_i <= B }.
程序采用增广走廊的前向递推（[0,B] 网格，B=1）：
    g_0(delta_0) == 1,
    g_i(delta_i) = int_0^{min(B,delta_i+d_i)} g_{i-1}(delta_{i-1}) d delta_{i-1},
    Q_n = int_0^B g_n.
这里 Q_n 还对块首滞留 delta_0 积分。正文的确定性夹逼证明
Q_n 与固定 delta_0=0 的 Vol(A_n) 只相差次指数因子，故二者具有同一
Lambda = lim (1/n)log2 volume。数值递推因此计算该共同的路径条件指数。
退火：平均核 (T g)(d') = int_0^B g(d) Sbar((d-d')^+) dd 的顶特征值对数，Sbar(t)=P(D>=t)。

一般理论允许平稳遍历相关宿主；主表在独立宿主基线上计算，相关宿主由
correlation_full.py 调用同型的路径递推。约定 B=1，故
Delta_minus=-Lambda_q。结构化结果同时记录离散参数、路径重复数、Monte Carlo
标准误和代表性网格细化检查。
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "result")
B = 1.0
GRID_POINTS = 600
PATH_LENGTH = 24000
SEED = 0
REPETITIONS = 20
GRID_CHECK_POINTS = (300, 600, 1200)
GRID_CHECK_LENGTH = 12000
GRID_CHECK_SEED = 19
GRID_CHECK_REPETITIONS = 2


# ---- 路径条件指数：增广走廊体积的随机递推 -----------------------------
def lambda_quenched_path(path_sampler, G=600, n=6000, seed=0, reps=3):
    x = np.linspace(0.0, B, G)
    h = x[1] - x[0]
    out = []
    for r in range(reps):
        rng = np.random.default_rng(seed + r)
        path = np.asarray(path_sampler(rng, n), dtype=float)
        if path.shape != (n,) or np.any(path <= 0.0):
            raise ValueError("path_sampler must return n positive host intervals")
        g = np.ones(G)              # g_1 于 [0,B]
        log2vol = 0.0               # 累计被剥离的归一化因子之 log2
        for d in path:
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * h)])  # cum[j]=int_0^{x_j} g
            upper = np.minimum(B, x + d)
            gnew = np.interp(upper, x, cum)   # g_i(delta)=int_0^{min(B,delta+d)} g_{i-1}
            s = max(gnew.max(), 1e-300)
            gnew /= s
            log2vol += math.log2(s)
            g = gnew
        Lam = (log2vol + math.log2(max(np.trapezoid(g, x), 1e-300))) / n
        out.append(Lam)
    sample_sd = float(np.std(out, ddof=1)) if len(out) > 1 else 0.0
    return float(np.mean(out)), sample_sd


def lambda_quenched(sampler, G=600, n=6000, seed=0, reps=3):
    """兼容原独立宿主接口；一般相关宿主应调用 lambda_quenched_path。"""
    return lambda_quenched_path(
        lambda rng, length: np.asarray([sampler(rng) for _ in range(length)]),
        G=G, n=n, seed=seed, reps=reps,
    )


# ---- 退火指数：平均转移核顶特征值 -------------------------------------------
def lambda_annealed(survival, G=600, iters=6000, tol=1e-12):
    x = np.linspace(0.0, B, G)
    h = x[1] - x[0]
    weights = np.full(G, h)
    weights[[0, -1]] *= 0.5
    dd = x[None, :] - x[:, None]              # x_k - x_j
    Kmat = np.where(dd > 0, survival(np.maximum(dd, 0.0)), 1.0)
    g = np.ones(G)
    lam = 1.0
    for _ in range(iters):
        gn = Kmat @ (g * weights)
        s = np.trapezoid(gn, x)
        gn = gn / max(s, 1e-300)
        lam_new = s / np.trapezoid(g, x)
        if abs(lam_new - lam) < tol:
            lam, g = lam_new, gn
            break
        lam, g = lam_new, gn
    return math.log2(max(lam, 1e-300))


# ---- 宿主：采样器 + 生存函数 Sbar(t)=P(D>=t) --------------------------------
def exp_host(mean):
    return (lambda rng: rng.exponential(mean),
            lambda t: np.exp(-t / mean))


def unif_host(a, b):
    # D ~ U(a,b);  Sbar(t)=1 (t<=a), (b-t)/(b-a) (a<t<b), 0 (t>=b)
    def surv(t):
        t = np.asarray(t, dtype=float)
        return np.clip((b - t) / (b - a), 0.0, 1.0)
    return (lambda rng: rng.uniform(a, b), surv)


def det_host(val):
    return (lambda rng: val,
            lambda t: (np.asarray(t) <= val).astype(float))


# 与 tab:penalty 六行一一对应
HOSTS = [
    ("Exp(mean=2B)",          exp_host(2.0)),
    ("U(0,2B)",               unif_host(0.0, 2.0)),
    ("Exp(mean=B)",           exp_host(1.0)),
    ("Exp(mean=B/2)",         exp_host(0.5)),
    ("U(B/2,3B/2)",           unif_host(0.5, 1.5)),
    ("U(B,2B) (sparse)",      unif_host(1.0, 2.0)),
]


def main(quick=False):
    grid_points = 300 if quick else GRID_POINTS
    path_length = 6000 if quick else PATH_LENGTH
    repetitions = 3 if quick else REPETITIONS
    grid_check_points = (150, 300, 600) if quick else GRID_CHECK_POINTS
    grid_check_length = 3000 if quick else GRID_CHECK_LENGTH
    grid_check_repetitions = 1 if quick else GRID_CHECK_REPETITIONS
    print("=" * 66)
    print("tab:penalty 复现 —— Lambda_q / Delta_minus / (1+Delta_minus)   (B=1, log2 B=0)")
    print("=" * 66)
    print(f"{'host':<22}{'Lam_q':>10}{'Lam_ann':>10}{'Delta_-':>11}{'1+Delta_-':>12}")
    grid_values = {}
    exp2_sampler = HOSTS[0][1][0]
    for grid in grid_check_points:
        lam_grid, _ = lambda_quenched(
            exp2_sampler,
            G=grid,
            n=grid_check_length,
            seed=GRID_CHECK_SEED,
            reps=grid_check_repetitions,
        )
        grid_values[str(grid)] = round(-lam_grid, 8)
    grid_span = max(grid_values.values()) - min(grid_values.values())

    rows = {
        "_meta": {
            "quick": quick,
            "schema_version": 3,
            "script": "lambda_delta_minus.py",
            "paper_symbol": "Delta_minus",
            "buffer_B": B,
            "log_unit": "bits/packet",
            "estimator": (
                "path-conditioned augmented-corridor volume exponent; "
                "equal to fixed-start exponent by deterministic sandwich"
            ),
            "grid_points": grid_points,
            "path_length": path_length,
            "seed": SEED,
            "repetitions": repetitions,
            "reported_uncertainty": "Monte Carlo standard error across paths",
            "grid_convergence_check": {
                "host": "Exp(mean=2B)",
                "path_length": grid_check_length,
                "seed": GRID_CHECK_SEED,
                "repetitions": grid_check_repetitions,
                "Delta_minus_by_grid": grid_values,
                "max_span": round(grid_span, 8),
            },
        }
    }
    for name, (smp, surv) in HOSTS:
        lq, sd = lambda_quenched(
            smp,
            G=grid_points,
            n=path_length,
            seed=SEED,
            reps=repetitions,
        )
        la = lambda_annealed(surv, G=GRID_POINTS)
        delta_minus = -lq                 # log2 B - Lambda, B=1
        se = sd / math.sqrt(repetitions)
        print(f"{name:<22}{lq:>10.4f}{la:>10.4f}{delta_minus:>11.3f}{1 + delta_minus:>12.3f}"
              f"  se={se:.4f}")
        rows[name] = dict(lambda_quenched=round(lq, 4), lambda_quenched_std=round(sd, 4),
                          Delta_minus_mc_se=round(se, 5),
                          lambda_annealed=round(la, 4), Delta_minus=round(delta_minus, 3),
                          one_plus_Delta_minus=round(1 + delta_minus, 3))
    os.makedirs(RES, exist_ok=True)
    filename = "lambda_delta_minus_quick.json" if quick else "lambda_delta_minus.json"
    with open(os.path.join(RES, filename), "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print("\nwrote result/%s" % filename)
    return rows


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
