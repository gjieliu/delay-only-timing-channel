"""多分辨率模格定理的独立数值与几何测试。"""

import json
import math
from pathlib import Path

import numpy as np

from multires_bound import (
    HOSTS,
    circle_capacity,
    full_width_modulo_benchmark,
    geometric_probabilities,
    multires_rate,
    partition_penalty,
    residue_representative,
    select_beta_bounded,
    symmetric_markov_entropy_rate,
)


def test_circle_asymptotic():
    for ratio in (18.0, 32.0, 128.0):
        exact = circle_capacity(ratio)
        asym = math.log2(ratio / math.sqrt(2.0 * math.pi * math.e))
        assert abs(exact - asym) < 1e-12
    vals = [circle_capacity(r) for r in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)]
    assert all(a <= b + 1e-10 for a, b in zip(vals, vals[1:]))
    assert vals[0] >= 0.0


def test_circle_global_envelope_and_infinite_tier_tail():
    """核验无穷层极限所用的全局容量包络与有限均值尾控制。"""
    ratios = np.geomspace(0.125, 256.0, 80)
    for ratio in ratios:
        upper = 0.5 * math.log2(1.0 + ratio * ratio / 12.0)
        assert circle_capacity(float(ratio)) <= upper + 3e-9

    beta = 2.0
    tail_ratio = 0.65
    k = np.arange(200, dtype=float)
    p = (1.0 - tail_ratio) * tail_ratio**k
    mean_k = float(np.dot(p, k))
    errors = []
    for gamma in (2.0**10, 2.0**16, 2.0**22, 2.0**28):
        exact = float(np.dot(
            p,
            [circle_capacity(gamma * beta ** (-int(j))) for j in k],
        ))
        asym = (
            math.log2(gamma / math.sqrt(2.0 * math.pi * math.e))
            - mean_k * math.log2(beta)
        )
        errors.append(abs(exact - asym))
    assert all(b < a for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 5e-5


def test_partition_and_sparse_endpoint():
    for host in HOSTS:
        opt = select_beta_bounded(host.cdf)
        p = geometric_probabilities(host.cdf, opt["beta_hat"])
        penalty = partition_penalty(host.cdf, opt["beta_hat"])
        assert abs(float(np.sum(p)) - 1.0) < 1e-12
        assert abs(penalty["J_beta"] - opt["J_beta_hat"]) < 1e-12
        assert opt["J_beta_hat"] >= -1e-12
    sparse = HOSTS[-1]
    opt = select_beta_bounded(sparse.cdf)
    assert abs(opt["J_beta_hat"]) < 1e-12
    for gamma in (8.0, 64.0, 1024.0):
        assert abs(multires_rate(sparse.cdf, gamma, opt["beta_hat"])
                   - full_width_modulo_benchmark(gamma)) < 1e-10


def test_delta_plus_global_certificate():
    """独立抽查表 2 全局证书的解析端点、见证与可解特例。"""
    path = Path(__file__).resolve().parent / "result" / "delta_plus_certificate.json"
    cert = json.loads(path.read_text(encoding="utf-8"))
    assert cert["schema_version"] == 1
    hosts = cert["hosts"]

    exp_cert = hosts["Exp(mean=2B)"]
    lo, hi = exp_cert["Delta_plus_interval"]
    assert exp_cert["status"] == "certified"
    assert 0 <= hi - lo <= 5.01e-5
    assert round(lo, 3) == round(hi, 3) == 2.548
    assert min(exp_cert["endpoint_lower_bounds"].values()) > hi
    exp_host = next(h for h in HOSTS if h.name == "Exp(mean=2B)")
    witness = partition_penalty(exp_host.cdf, exp_cert["beta_witness"])["J_beta"]
    assert lo <= witness <= hi + 2e-7

    uniform = next(h for h in HOSTS if h.name == "U(B/2,3B/2)")
    assert hosts["U(B/2,3B/2)"]["Delta_plus_interval"] == [1.5, 1.5]
    for beta in np.geomspace(1.001, 32.0, 600):
        assert partition_penalty(uniform.cdf, float(beta))["J_beta"] >= 1.5 - 2e-12

    sparse = hosts["U(B,2B) (sparse)"]
    assert sparse["Delta_plus_interval"] == [0.0, 0.0]


def test_reachable_interval_residue():
    rng = np.random.default_rng(17)
    beta = 3.7
    for d in (1e-8, 0.2, 1.0, 3.0):
        worst_delta = 1.0
        lo = max(d - worst_delta, 0.0)
        hi = d + 1.0 - worst_delta
        assert abs((hi - lo) - min(d, 1.0)) < 1e-12

    for _ in range(20000):
        delta = rng.uniform(0.0, 1.0)
        d = rng.exponential(1.3)
        clipped = min(1.0, d)
        k = 0 if clipped >= 1.0 else int(math.ceil(math.log(1.0 / clipped, beta)))
        period = beta ** (-k)
        assert period <= clipped + 1e-12
        lo = max(d - delta, 0.0)
        hi = d + 1.0 - delta
        assert abs((hi - lo) - (1.0 - max(delta - d, 0.0))) < 1e-12
        assert hi - lo + 1e-12 >= clipped
        assert hi - lo + 1e-12 >= period
        residue = rng.uniform(0.0, period)
        sent_gap = residue_representative(lo, hi, period, residue)
        next_delta = delta + sent_gap - d
        injection = next_delta - max(delta - d, 0.0)
        representative_index = (sent_gap - residue) / period
        assert -1e-10 <= next_delta <= 1.0 + 1e-10
        assert injection >= -1e-10
        assert representative_index >= -1e-10
        assert abs(representative_index - round(representative_index)) < 1e-8


def test_asymptotic_constant_gap():
    host = HOSTS[0]  # Exp(mean=2B)
    opt = select_beta_bounded(host.cdf)
    gamma = 2.0**18
    gap = (full_width_modulo_benchmark(gamma)
           - multires_rate(host.cdf, gamma, opt["beta_hat"]))
    assert abs(gap - opt["J_beta_hat"]) < 2e-3


def test_correlated_label_entropy_rate():
    assert abs(symmetric_markov_entropy_rate(0.5) - 1.0) < 1e-12
    assert abs(symmetric_markov_entropy_rate(0.9) - 0.4689955936) < 1e-9
    assert abs(symmetric_markov_entropy_rate(0.99) - 0.0807931359) < 1e-9
    assert symmetric_markov_entropy_rate(0.99) < symmetric_markov_entropy_rate(0.9)


def test_joint_metric_wrong_message_normalization():
    """独立均匀错误码字下，联合译码似然比的条件均值应为一。"""
    for period in (0.7, 1.3, 2.5):
        ngrid = 40000
        residues = (np.arange(ngrid) + 0.5) * period / ngrid
        for observation in (0.17, 1.91):
            z = np.mod(observation - residues, period)
            density = np.zeros_like(z)
            kmax = int(math.ceil(10.0 / period)) + 2
            for k in range(-kmax, kmax + 1):
                t = z + k * period
                density += np.exp(-0.5 * t * t) / math.sqrt(2.0 * math.pi)
            likelihood_ratio = period * density
            assert abs(float(np.mean(likelihood_ratio)) - 1.0) < 2e-6


def test_free_corridor_scaling_step():
    """独立检查自由起点体积夹逼所用的缩放映射。"""
    rng = np.random.default_rng(29)
    for q in range(1, 10):
        alpha = q / (q + 1.0)
        for _ in range(200):
            partial = np.concatenate([[0.0], rng.uniform(0.0, 1.0, size=q)])
            increments = np.diff(partial)
            d = np.maximum(1e-3, -increments + rng.uniform(0.0, 1.0, size=q))
            assert np.all(increments >= -d - 1e-12)
            span = float(np.max(partial) - np.min(partial))
            assert span <= 1.0 + 1e-12

            scaled = alpha * increments
            scaled_partial = np.concatenate([[0.0], np.cumsum(scaled)])
            scaled_span = float(np.max(scaled_partial) - np.min(scaled_partial))
            assert np.all(scaled >= -d - 1e-12)
            assert scaled_span <= alpha * span + 1e-12
            assert 1.0 - scaled_span >= 1.0 / (q + 1.0) - 1e-12


def _corridor_volume_grid(path, augmented, grid_points=2401):
    """网格积分固定起点体积 V_q 或增广体积 Q_q（B=1）。"""
    x = np.linspace(0.0, 1.0, grid_points)
    h = x[1] - x[0]
    g = np.ones(grid_points)
    steps = path if augmented else path[1:]
    for d in steps:
        cumulative = np.concatenate(
            [[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * h)]
        )
        upper = np.minimum(1.0, x + float(d))
        g = np.interp(upper, x, cumulative)
    return float(np.trapezoid(g, x))


def test_fixed_corridor_subadditivity_and_augmented_bridge():
    """独立检查 Kingman 所需次可加性及 Q_q--V_q 确定性夹逼。"""
    rng = np.random.default_rng(41)
    for q in range(1, 8):
        path = rng.uniform(0.05, 1.4, size=q)
        fixed = _corridor_volume_grid(path, augmented=False)
        augmented = _corridor_volume_grid(path, augmented=True)
        alpha = q / (q + 1.0)
        lower = alpha**q * fixed / (q + 1.0)
        assert lower <= augmented + 2e-3
        assert augmented <= fixed + 2e-3
        envelope = float(np.prod(np.minimum(path, 1.0)))
        assert envelope <= fixed + 2e-3
        assert fixed <= 1.0 + 2e-3

    for _ in range(20):
        left = rng.uniform(0.05, 1.4, size=3)
        right = rng.uniform(0.05, 1.4, size=4)
        whole = _corridor_volume_grid(
            np.concatenate([left, right]), augmented=False
        )
        product = (
            _corridor_volume_grid(left, augmented=False)
            * _corridor_volume_grid(right, augmented=False)
        )
        assert whole <= product + 3e-3

    sparse = np.full(7, 1.2)
    assert abs(_corridor_volume_grid(sparse, augmented=False) - 1.0) < 2e-3


def main():
    test_circle_asymptotic()
    test_circle_global_envelope_and_infinite_tier_tail()
    test_partition_and_sparse_endpoint()
    test_delta_plus_global_certificate()
    test_reachable_interval_residue()
    test_asymptotic_constant_gap()
    test_correlated_label_entropy_rate()
    test_joint_metric_wrong_message_normalization()
    test_free_corridor_scaling_step()
    test_fixed_corridor_subadditivity_and_augmented_bridge()
    print("multires theory: ok")


if __name__ == "__main__":
    main()
