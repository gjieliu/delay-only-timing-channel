"""走廊体积指数 Lambda 与 Delta_minus = log2(B) - Lambda 的正算子迭代数值求解 (B=1)。
(T_d f)(x) = int_0^{min(1, x+d)} f(t) dt ;  phi_1 = 1 ;  V_n = int_0^1 phi_n .
"""
import numpy as np

B = 1.0

def delta_minus(dseq, N=1600):
    dt = 1.0 / N
    t = (np.arange(N) + 0.5) * dt          # 中点网格
    edges = (np.arange(N + 1)) * dt        # 累积积分的节点
    phi = np.ones(N)
    logsum = 0.0
    for d in dseq:
        C = np.concatenate(([0.0], np.cumsum(phi) * dt))   # C[k] = int_0^{edges[k]} phi
        y = np.minimum(1.0, t + d)
        phi = np.interp(y, edges, C)
        V = phi.sum() * dt
        if V <= 0:
            raise RuntimeError('underflow')
        logsum += np.log2(V)
        phi = phi / V
    Lam = logsum / len(dseq)
    return np.log2(B) - Lam


# ---------- 宿主生成器 ----------
def host_iid(rng, n, sampler):
    return sampler(rng, n)

def host_2state(rng, n, theta, lo0, hi0, lo1, hi1):
    """对称二状态持续链，P(S_i=S_{i-1})=theta；状态0->Unif(lo0,hi0)，状态1->Unif(lo1,hi1)"""
    s = rng.integers(0, 2)
    out = np.empty(n)
    stay = rng.random(n) < theta
    u = rng.random(n)
    for i in range(n):
        if i > 0 and not stay[i]:
            s = 1 - s
        if s == 0:
            out[i] = lo0 + (hi0 - lo0) * u[i]
        else:
            out[i] = lo1 + (hi1 - lo1) * u[i]
    return out

def host_control(rng, n, theta, p):
    """丁族对照：子缓冲条件律固定为 Unif(B/4,3B/4)（由持续二状态 L/H 承载相关性），
    独立地以概率 p 插入满宽度步 Unif(B,3B/2)。P(d>=B)=p，而 d<B 的条件律与 p 无关。"""
    s = rng.integers(0, 2)
    stay = rng.random(n) < theta
    z = rng.random(n) < p
    u = rng.random(n)
    out = np.empty(n)
    for i in range(n):
        if i > 0 and not stay[i]:
            s = 1 - s
        if z[i]:
            out[i] = 1.0 + 0.5 * u[i]              # Unif(B, 3B/2)
        elif s == 0:
            out[i] = 0.25 + 0.25 * u[i]            # Unif(B/4, B/2)
        else:
            out[i] = 0.50 + 0.25 * u[i]            # Unif(B/2, 3B/4)
    return out
