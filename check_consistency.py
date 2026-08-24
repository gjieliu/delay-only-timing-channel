# -*- coding: utf-8 -*-
"""投稿正文、技术证明附录、完整结构化结果与图文件的一致性校验。

完整 JSON 继续覆盖独立基线、相关宿主、单字母窗口、清零对照、GMI 常数及有限噪声检查；
正文只要求出现投稿叙事所需的两张紧凑表、一张差距图与 GMI 正反例。这样既防止主文数字
漂移，也不迫使补充验证层重新膨胀回正文。
"""

import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def load(name):
    with open(os.path.join(HERE, "result", name), encoding="utf-8") as f:
        return json.load(f)


def norm(s):
    return re.sub(r"[ \t\n]+", " ", s)


def require(tex, fragment):
    assert norm(fragment) in tex, "稿件缺失片段: %r" % fragment


def r4(x):
    return round(x, 4)


def main():
    with open(os.path.join(ROOT, "manuscript", "main.tex"), encoding="utf-8") as f:
        raw = f.read()
    tex = norm(raw)
    with open(os.path.join(ROOT, "manuscript", "appendix_proofs.tex"),
              encoding="utf-8") as f:
        appendix_raw = f.read()
    proof_tex = norm(raw + "\n" + appendix_raw)
    with open(os.path.join(ROOT, "sketelon.md"), encoding="utf-8") as f:
        skeleton_raw = f.read()
    skeleton = norm(skeleton_raw)
    with open(os.path.join(ROOT, "design.md"), encoding="utf-8") as f:
        design_raw = f.read()
    with open(os.path.join(ROOT, "manuscript", "supplement", "numerical_evidence.md"),
              encoding="utf-8") as f:
        supplement_raw = f.read()
    multi = load("multires.json")
    dplus_cert = load("delta_plus_certificate.json")
    dminus = load("lambda_delta_minus.json")
    corr = load("correlation_full.json")
    rw = load("delta_minus_width_bounds.json")
    gmi = load("gmi_demo.json")
    gfin = load("gmi_finite.json")
    reset = load("reset_control.json")

    # ---- 代码—论文符号接口：旧 schema 与旧产物名不得回流 ----
    symbol_results = (multi, dminus, corr, rw, gmi)
    serialized = "\n".join(json.dumps(x, ensure_ascii=False) for x in symbol_results)
    for obsolete in ('"Dlay"', '"delta_ref', '"total"'):
        assert obsolete not in serialized, "旧结果字段回流: %s" % obsolete
    for name, row in multi.items():
        if name == "_meta":
            continue
        assert "beta_star" not in row["partition"]
        assert "Delta_plus" not in row["partition"]
    for obsolete_name in (
        "lambda_delta_ref.py",
        "delta_ref_width_bounds.py",
        os.path.join("result", "lambda_delta_ref.json"),
        os.path.join("result", "delta_ref_width_bounds.json"),
    ):
        assert not os.path.exists(os.path.join(HERE, obsolete_name)), \
            "旧代码或结果文件回流: %s" % obsolete_name

    # ---- 全局措辞与结构 ----
    assert "无代价性二分律" not in raw
    assert "确定型宿主的精确高信噪比容量" not in raw
    assert "宿主间隔 $\\{d_i\\}$ 独立同分布" not in raw
    require(tex, "平稳遍历宿主下只延后时序信道容量界")
    require(tex, "多分辨率因果模格可达下界")
    require(tex, "宿主过程律诱导的容量夹逼")
    require(tex, "宿主过程律诱导的常数阶容量夹逼")
    require(tex, "\\overline H(\\mathbf K)")
    require(tex, "\\label{asm:alignment}")
    require(tex, "\\label{def:capacity}")
    require(proof_tex, "\\label{eq:joint-metric}")
    require(proof_tex, "\\label{eq:circle-global-bound}")
    require(proof_tex, "\\label{eq:corridor-subadditive}")
    require(proof_tex, "\\label{eq:free-volume-sandwich}")
    require(proof_tex, "\\label{eq:fixed-augmented-sandwich}")
    require(tex, "\\input{appendix_proofs}")
    for appendix_label in (
        "app:multires-proof", "app:circle-proofs", "app:gmi-proofs",
        "app:corridor-exponent", "app:outer-proof",
    ):
        require(proof_tex, "\\label{%s}" % appendix_label)
    assert appendix_raw.count("\\begin{proof}[") == 10, \
        "技术附录应恰含十个具名完整证明"
    for proof_label in (
        "prop:multires", "lem:circle-asymptotic", "lem:circle-tail",
        "lem:doeblin", "prop:gmi", "prop:gmi-const", "lem:lyapunov",
        "prop:ref-width", "lem:support-volume", "prop:outer",
    ):
        require(proof_tex, "\\ref{%s}的证明" % proof_label)
    require(tex, "正文对长证明保留可独立审阅的证明路线")
    require(tex, "完整技术推导统一列于文末附录")
    require(tex, "\\newcommand{\\snr}{\\mathsf{snr}}")
    require(tex, "\\newcommand{\\Rfw}{C_0}")
    require(tex, "\\newcommand{\\Dlayer}{\\Delta_{+}}")
    require(tex, "\\newcommand{\\Jlayer}{J}")
    require(tex, "\\newcommand{\\Dref}{\\Delta_{-}}")
    assert "R_{\\mathrm{fw}}" not in raw
    assert "\\Delta_{\\mathrm{lay}}" not in raw
    assert "\\Delta_{\\mathrm{ref}}" not in raw
    assert "\\Dlayer^\\star" not in raw
    assert "\\Dlayer(\\beta" not in raw
    require(tex, "\\Jlayer_\\beta(P_{\\mathbf d})")
    require(tex, "\\Dlayer(P_{\\mathbf d}):=\\inf_{\\beta>1}\\Jlayer_\\beta(P_{\\mathbf d})")
    require(tex, "\\widehat\\beta")
    require(tex, "\\Jlayer_{\\widehat\\beta}")
    require(tex, "给出 $\\Dlayer$ 的全局证书")
    assert "$\\beta^\\star=3.837$" not in raw
    assert "可证数值展示区间" not in raw
    assert "两端相向移动证明的是可证区间" not in raw
    require(tex, "\\label{eq:intro-capacity-sandwich}")
    require(tex, "\\label{eq:capacity-sandwich}")
    assert "O_{P_{\\mathbf d}}(1)" not in raw
    assert "\\label{eq:constant-order}" not in raw
    assert "\\gamma:=\\buf/\\jit" not in raw
    assert "\\Delta_{\\mathrm{geo}}" not in raw
    assert "\\rhogeom" not in raw and "\\rho_{\\mathrm{geom}}" not in raw
    assert "\\varpi" not in raw and "\\varrho" not in raw
    require(tex, "满宽度模环基准")
    require(proof_tex, "固定起点体积的次可加性")
    require(proof_tex, "无穷层模环容量的统一尾控制")
    # 骨架—设计—补充材料的契约一致性。
    assert "L_{\\mathrm{mean}}" not in skeleton_raw
    assert "L_{\\mathrm{avg}}" in skeleton_raw
    assert "\\inf_{\\beta>1}\\{H(K_{\\beta,1})" in skeleton_raw
    assert "`gap.json`" not in skeleton_raw
    assert "`multires.json`、`fig_multires_gap.pdf`" in skeleton_raw
    assert "可证数值展示区间" not in skeleton_raw
    assert "两端相向移动只证明" not in skeleton_raw
    assert "`beta_hat` 与 `J_beta_hat`" in skeleton_raw
    assert "`Delta_plus_interval`" in skeleton_raw
    assert "## 11. 2026-08-21 数值候选设计与理论端点分离及独立主例认证" in design_raw
    assert "## 13. 2026-08-21 正文证明路线与正式附录分层" in design_raw
    assert "`manuscript/appendix_proofs.tex`" in design_raw
    assert "`manuscript/appendix_proofs.tex`" in skeleton_raw
    assert "不能被归入外部补充材料" in skeleton_raw
    assert "2.54795143" in design_raw and "2.54800135" in design_raw
    assert "delta_plus_certificate.json" in supplement_raw
    assert "`fig_multires_rate.pdf` 列为补充材料图件" in supplement_raw
    # 译码端命题的软化表述（审稿修复轮）：越界措辞不得回归
    assert "只能是后一处" not in raw
    assert "这在互信息层面是唯一途径" not in raw
    require(tex, "尚不足以断言该放宽在整个码类上无损")
    # 小噪声桥接命题（fix33；fix35 起为"上界"表述）
    require(tex, "小噪声损失的显式常数上界")
    assert "小噪声端的显式费用常数" not in raw
    assert "极限常数" not in raw and "小噪声取值" not in raw
    require(tex, "\\label{prop:gmi-const}")
    require(tex, "\\label{eq:c1}")
    require(tex, "\\label{eq:gmi-const}")

    # ---- 独立宿主完整结果与投稿表 1 ----
    assert multi["_meta"]["schema_version"] == 4
    assert multi["_meta"]["paper_symbols"] == {
        "fixed_design_penalty": "J_beta",
        "reported_numerical_bound": "J_beta_hat",
        "theoretical_optimized_endpoint": "Delta_plus=inf_{beta>1} J_beta",
    }
    selection = multi["_meta"]["beta_selection"]
    assert selection["interval"] == [1.02, 20.0]
    assert selection["xatol"] == 1e-8
    # 局部候选搜索本身不冒充全局证书；证书由独立结果文件承担。
    assert selection["global_optimum_certified"] is False
    assert dplus_cert["schema_version"] == 1
    assert dplus_cert["quantity"] == "Delta_plus=inf_{beta>1} J_beta"
    assert set(dplus_cert["software"]) == {"python", "mpmath", "scipy"}
    chosts = dplus_cert["hosts"]
    assert set(chosts) == {
        "Exp(mean=2B)", "U(B/2,3B/2)", "U(B,2B) (sparse)"
    }
    exp_cert = chosts["Exp(mean=2B)"]
    dlo, dhi = exp_cert["Delta_plus_interval"]
    assert exp_cert["status"] == "certified"
    assert exp_cert["compact_beta_interval"] == [1.05, 17.0]
    assert exp_cert["certificate_width"] <= 5e-5
    assert exp_cert["mp_dps"] == 80 and exp_cert["tail_exponent"] == 38.0
    assert round(dlo, 3) == round(dhi, 3) == 2.548
    assert min(exp_cert["endpoint_lower_bounds"].values()) > dhi
    assert abs(exp_cert["beta_witness"]
               - multi["Exp(mean=2B)"]["partition"]["beta_hat"]) < 2e-6
    assert chosts["U(B/2,3B/2)"]["Delta_plus_interval"] == [1.5, 1.5]
    assert chosts["U(B,2B) (sparse)"]["Delta_plus_interval"] == [0.0, 0.0]
    assert multi["_meta"]["full_width_curve"].startswith("modulo-circle")
    assert dminus["_meta"]["schema_version"] == 3
    assert dminus["_meta"]["paper_symbol"] == "Delta_minus"
    assert dminus["_meta"]["estimator"].startswith(
        "path-conditioned augmented-corridor"
    )
    assert dminus["_meta"]["grid_convergence_check"]["max_span"] < 5e-4
    assert dminus["_meta"]["repetitions"] >= 12
    max_dminus_se = max(
        dminus[name]["Delta_minus_mc_se"] for name in dminus if name != "_meta"
    )
    require(tex, str(dminus["_meta"]["path_length"]))
    require(tex, str(dminus["_meta"]["repetitions"]))
    require(tex, "$%.4f$" % max_dminus_se)

    order = [
        "Exp(mean=2B)",
        "U(0,2B)",
        "Exp(mean=B)",
        "Exp(mean=B/2)",
        "U(B/2,3B/2)",
        "U(B,2B) (sparse)",
    ]
    for name in order:
        p = multi[name]["partition"]
        r = dminus[name]["Delta_minus"]
        assert dminus[name]["Delta_minus_mc_se"] < 0.002
        assert p["J_beta_hat"] >= -1e-12 and r >= -1e-12

    # 投稿表只保留一般非稀疏、GMI 正例和精确满宽度端点三行。
    for name in ("Exp(mean=2B)", "U(B/2,3B/2)"):
        p = multi[name]["partition"]
        r = dminus[name]["Delta_minus"]
        require(tex, "$%.3f$ & $%.3f$ & $%.3f$"
                % (p["beta_hat"], p["J_beta_hat"], r))
    sparse = multi["U(B,2B) (sparse)"]["partition"]
    assert abs(sparse["J_beta_hat"]) < 1e-12
    require(tex, "任意 & $0.000$ & $0.000$")
    require(tex, "$[0.143,2.548]$")
    require(tex, "容量损失数值外框")
    require(tex, "$\\Dlayer$ 已认证到表中精度")
    require(tex, "2.54795143\\le\\Dlayer\\le2.54800135")
    require(tex, "整个区间仍不是经连续算子误差认证的容量区间")
    require(tex, "这一认证不延伸到相关宿主表")

    # ---- 相关宿主完整结果与投稿表 2(a) ----
    cm = corr["meta"]
    assert cm["quick"] is False
    assert cm["schema_version"] == 3
    assert cm["paper_symbols"] == {
        "fixed_design_penalty": "J_beta",
        "corridor_lower_endpoint": "Delta_minus",
    }
    assert cm["grid_points"] == 1600 and cm["path_length"] == 2000000
    assert cm["repetitions"] == 3 and cm["rule"] == "midpoint"
    require(tex, "1600 点网格")
    require(tex, "$2\\times10^6$")
    require(tex, "作 3 次独立重复")

    fams = {f["key"]: f for f in corr["families"]}
    assert list(fams) == ["jia", "yi", "bing"]
    optimized_caps = {"jia": math.inf, "yi": 2.0, "bing": 3.0}
    for f in corr["families"]:
        assert len(f["rows"]) == 3
        for row in f["rows"]:
            assert abs(row["J_beta"] - row["EK_log2beta"] - row["H_rate"]) < 1e-12
            assert abs(row["J_beta_minus_Delta_minus"]
                       - row["J_beta"] + row["Delta_minus_mean"]) < 1e-12
            dopt = min(row["J_beta"], optimized_caps[f["key"]])
            width_opt = dopt - row["Delta_minus_mean"]
            assert width_opt >= 0

    # 主文只报告甲、乙两族的两个端点；中间点与丙族留在完整结果层。
    reset_prob = {"jia": "$1/2$", "yi": "$0$"}
    family_tex = {"jia": "甲", "yi": "乙"}
    for key in ("jia", "yi"):
        for row in (fams[key]["rows"][0], fams[key]["rows"][2]):
            dopt = min(row["J_beta"], optimized_caps[key])
            width_opt = dopt - row["Delta_minus_mean"]
            require(
                tex,
                "%s & %s & $%.2f$ & $%.3f$ & $%.4f$ & $%.3f$"
                % (family_tex[key], reset_prob[key], row["theta"], dopt,
                   row["Delta_minus_mean"], width_opt),
            )

    # 完整结果层仍核验网格、路径及三族配对稳定性。
    assert corr["max_spread"] <= 0.0026
    require(tex, "各行实现间极差不超过 $0.0026$")
    gcs = {g["theta"]: g["Delta_minus_by_grid"] for g in corr["grid_checks"]}
    assert set(gcs) == {0.5, 0.99}
    for g in gcs.values():
        assert set(g) == {"800", "1600", "3200"}
        assert len({"%.4f" % v for v in g.values()}) == 1, "网格细化四位小数应不变"
    pc = corr["path_check"]["Delta_minus_by_length"]
    assert sorted(int(k) for k in pc) == [500000, 1000000, 2000000]
    pl_bound = math.ceil((max(pc.values()) - min(pc.values())) * 1e4) / 1e4
    assert pl_bound <= 0.0010
    # 持续性效应：Delta_- 的差值与倍数（按表中四位小数导出）
    deltas = {k: r4(r4(fams[k]["rows"][2]["Delta_minus_mean"])
                    - r4(fams[k]["rows"][0]["Delta_minus_mean"])) for k in fams}
    require(tex, "而 $\\Dref$ 的离散估计分别上升 $%.4f$ 与 $%.4f$ 比特/包"
            % (deltas["jia"], deltas["yi"]))
    require(tex, "$J_{\\mathrm{sel}}$")
    jia_spread = max(row["Delta_minus_spread"] for row in fams["jia"]["rows"])
    assert round(deltas["jia"] / jia_spread) == 3
    # 配对标准误：同种子共用底层随机流；移动量须为其十倍以上
    pses = [fams[k]["Delta_minus_change_paired_se"] for k in ("jia", "yi", "bing")]
    require(tex, "配对标准误约为 $%.5f$ 与 $%.5f$" % tuple(pses[:2]))
    for k, se in zip(("jia", "yi", "bing"), pses):
        assert deltas[k] >= 10 * se, "配对显著性不足: %s" % k

    # ---- 单字母窗口：完整 15 点校验，正文只报告总通过数 ----
    assert rw["meta"]["quick"] is False
    assert rw["meta"]["schema_version"] == 3
    assert rw["meta"]["paper_symbol"] == "Delta_minus"
    assert rw["meta"]["correlation_source"] == "correlation_full.json"
    assert rw["all_contained"] and rw["checks_total"] == 15
    require(tex, "全部 15/15 个")
    for f in rw["families"]:
        assert f["Delta_minus_lo"] <= f["Delta_minus_hi"]
        assert 0 <= f["tightness_min"] <= f["tightness_max"] <= 100

    # ---- 清零对照完整结果与投稿表 2(b) ----
    rm = reset["meta"]
    assert rm["quick"] is False and rm["paper_symbol"] == "Delta_minus"
    assert rm["schema_version"] == 2
    assert rm["grid_points"] == 1200 and rm["path_length"] == 400000
    assert rm["repetitions"] == 4 and rm["seeds"] == [1, 2, 3, 4]
    assert rm["paired_random_streams_across_theta"] is True
    assert rm["rule"] == "midpoint"
    reset_rows = reset["rows"]
    assert [row["p"] for row in reset_rows] == [0.0, 0.1, 0.25, 0.5]
    ratios = []
    for row in reset_rows:
        lo = row["theta"]["0.5"]["mean"]
        hi = row["theta"]["0.99"]["mean"]
        assert row["window_lo"] <= lo <= row["window_hi"]
        assert row["window_lo"] <= hi <= row["window_hi"]
        ratios.append(row["ratio"])
    assert all(a > b for a, b in zip(ratios, ratios[1:]))
    assert all(row["sensitivity"] <= row["mass_dilution_prediction"] + 1e-12
               for row in reset_rows)
    for row in (reset_rows[0], reset_rows[-1]):
        lo = row["theta"]["0.5"]["mean"]
        hi = row["theta"]["0.99"]["mean"]
        require(
            tex,
            "$%.2f$ & $%.4f$ & $%.4f$ & $%.5f$ & $%.5f$ & $%.2f$"
            % (row["p"], lo, hi, row["sensitivity"],
               row["mass_dilution_prediction"], row["ratio"]),
        )
    max_se = math.ceil(reset["max_standard_error"] * 10000) / 10000
    low_theta_se = max(row["theta"]["0.5"]["se"] for row in reset_rows)
    low_theta_se = math.ceil(low_theta_se * 100000) / 100000
    require(tex, "最大标准误为 $%.4f$" % max_se)
    assert low_theta_se <= 0.00026
    yi = {row["theta"]: row["Delta_minus_mean"] for row in fams["yi"]["rows"]}
    cross_diff = max(abs(reset_rows[0]["theta"]["0.5"]["mean"] - yi[0.5]),
                     abs(reset_rows[0]["theta"]["0.99"]["mean"] - yi[0.99]))
    cross_diff = math.ceil(cross_diff * 10000) / 10000
    assert cross_diff <= 0.0017

    # ---- GMI 常数 c1（gmi_demo）与有限噪声检查（gmi_finite）----
    assert gmi["meta"]["quick"] is False
    assert gmi["meta"]["n"] == 300000 and len(gmi["meta"]["seeds"]) == 3
    assert gmi["meta"]["burn_in"] == 2000
    assert gmi["meta"]["n_pool"] == 4400 and gmi["meta"]["n_fx"] == 5600
    assert gmi["meta"]["n_grid"] == 4400 and gmi["meta"]["n_ev"] == 1800
    fa = gmi["families"]["Unif(B/2,3B/2), beta=2"]
    fb2 = gmi["families"]["Unif(B/4,3B/4), beta=2"]
    fb4 = gmi["families"]["Unif(B/4,3B/4), beta=4"]
    assert abs(fa["c1_mean"] - fb2["c1_mean"]) < 5e-4       # 同 beta 尺度自相似
    assert fa["c1_range"] < 0.02 and fb2["c1_range"] < 0.02
    assert fb4["c1_range"] < 0.02
    assert gmi["meta"]["schema_version"] == 3
    assert gmi["meta"]["paper_symbol"] == "J_beta"
    assert abs(fa["J_beta"] - 1.5) < 1e-12
    assert abs(fb2["J_beta"] - 2.5) < 1e-12
    assert abs(fb4["J_beta"] - 2.0) < 1e-12
    require(tex, "三个独立种子")
    require(tex, "$c_1=%.4f$" % fa["c1_mean"])
    require(tex, "实现间极差 $%.4f$" % fa["c1_range"])
    se_mean = sum(p["H_XU_se"] for p in fa["per_seed"]) / len(fa["per_seed"])
    assert se_mean < 0.004
    cc = gmi["convergence_check"]
    assert cc is not None and cc["max_shift"] < 0.005
    require(tex, "最大移动为 $%.4f$" % cc["max_shift"])
    require(tex, "Monte Carlo 点估计为 $%.3f$" % fa["tightened_bound"])
    require(tex, "这一比较不是确定性认证区间")
    require(tex, "$c_1=%.4f$ （实现间极差 $%.4f$）"
            % (fb4["c1_mean"], fb4["c1_range"]))
    require(tex, "点估计 $%.3f$ 与固定分层界 $\\Jlayer_4=%.3f$ 的差异不可分辨"
            % (fb4["tightened_bound"], fb4["J_beta"]))
    assert fa["c1_mean"] < 1.0 and fb2["c1_mean"] < 1.0      # 数值估计低于层熵率 1
    require(tex, "\\label{eq:single-threshold}")

    # 有限噪声 GMI 留在补充验证层；继续核验其与渐近上界的一致性。
    assert gfin["meta"]["quick"] is False
    xmax = gfin["meta"]["x_max"]
    limit_const = fa["tightened_bound"]
    by_sigma = {}
    for p in gfin["points"]:
        by_sigma.setdefault(p["sigma"], []).append(p)
        eps3 = p["sigma"] * math.log2(math.e) + (8 / math.log(2)) * 2 ** (
            -2 * (xmax / p["sigma"]) ** 2 * math.log2(math.e))
        assert p["G_bound"] <= limit_const + eps3 + 4 * p["G_bound_se"], \
            "有限噪声点超出命题界: sigma=%s" % p["sigma"]
    assert set(by_sigma) == {0.04, 0.02, 0.01}
    means = {s: sum(x["G_bound"] for x in ps) / len(ps)
             for s, ps in by_sigma.items()}
    assert means[0.04] < means[0.01] < limit_const     # 两端向上界常数移动且不越过
    nest = gfin["nesting_check"]
    assert nest is not None and nest["max_shift"] < 0.05
    nb = math.ceil(nest["max_shift"] * 100) / 100
    assert nb <= 0.04

    # ---- 图文件：三幅产物均保留，投稿正文只使用差距图 ----
    for name in ("fig_multires_rate", "fig_multires_gap", "fig_partition_penalty"):
        path = os.path.join(HERE, "fig", name + ".pdf")
        assert os.path.getsize(path) > 1000, "%s.pdf 缺失或过小" % name
    require(tex, "\\includegraphics[width=0.70\\linewidth]{../code/fig/fig_multires_gap.pdf}")
    require(tex, "\\label{fig:gap}")
    require(tex, "\\ref{fig:gap}")
    require(tex, "独立证书另保证 $\\Dlayer$ 在三位小数处同为 $2.548$")
    assert "fig_multires_rate.pdf" not in raw
    assert "fig_partition_penalty.pdf" not in raw
    assert "\\label{fig:rate}" not in raw and "\\label{fig:penalty}" not in raw
    assert "\\label{tab:reset-control}" not in raw

    # ---- Exp(mean=2B) 曲线的定性结论 ----
    c = multi["Exp(mean=2B)"]["curves"]
    assert c["1024"]["multires"] > c["1024"]["single_threshold"]
    gap = c["262144"]["full_width_modulo"] - c["262144"]["multires"]
    target = multi["Exp(mean=2B)"]["partition"]["J_beta_hat"]
    assert abs(gap - target) < 2e-3

    print("consistency: ok")


if __name__ == "__main__":
    main()
