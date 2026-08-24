# 只延后时序信道（有界滞留）容量界 —— 数值生产代码

本仓库是下述论文全部数值结果的生产代码：

> G. Liu, W. Liu and H. Ding, "Capacity Bounds for Delay-Only Timing Channels with
> Bounded Dwell," 投稿于 *IEEE Transactions on Information Theory*, 2026.

论文给出有界滞留预算下、只延后的分组时序信道容量的双侧界

```
C_0 - Delta_+  <=  C  <=  C_0 - Delta_- + o(1),   C_0 = log2( B / (sqrt(2*pi*e) * sigma) ),
```

其中 `Delta_+` 是多分辨率模格构造给出的可达端点，`Delta_-` 是走廊体积指数给出的逆定理
端点。本仓库的脚本对一组宿主过程计算这两个端点，对表 2 的三个独立宿主给出 `Delta_+`
的全局证书，并绘制图件。

**本仓库只收录源码。** 结果 JSON 与图件 PDF 都是产物，不入库，由 `.gitignore` 排除；
运行脚本即可重新生成。论文里每一个数字对应哪个脚本、写到哪个字段、应当得到什么值，
逐条列在 `MANIFEST.md`，读者不需要我们的产物副本就能核对自己的运行结果。

目录是扁平的：全部脚本位于仓库根目录，运行后在同级生成 `result/`（JSON）与
`fig/`（PDF），两者均不入库。

`make_fig2.py` 之外的脚本与生产链**逐字节相同**，注释与说明保持中文原样、未作翻译，
以保证本仓库的代码可证就是产出论文数字的代码。

## 依赖

Python ≥ 3.9、NumPy ≥ 1.22、SciPy（分层比的一维有界搜索）、mpmath（`Delta_+` 证书的
80 位区间算术）；`make_figs.py` 与 `make_fig2.py` 另需 Matplotlib。见 `requirements.txt`。

版本记录有两份，此处如实并列，不作事后统一：

* 蒙特卡洛与求积结果由 **Python 3.11.15 + NumPy 2.4.4** 产生；
* `Delta_+` 全局证书由 **Python 3.12.2 + mpmath 1.3.0 + SciPy 1.18.0** 产生（每次运行
  都会把当时的版本写入证书输出的 `software` 字段）。

全部随机流由显式种子经 `numpy.random.default_rng` 派生，同版本 NumPy 下逐位可复现；
证书与多分辨率基线是确定性的。

本仓库的扁平布局已在第三套配置（Python 3.11.15 + mpmath 1.4.1 + SciPy 1.17.1）上实跑
验证：

* `certify_delta_plus.py` 除记录运行版本的 `software` 字段外输出逐位相同——同一区间
  `[2.5479514302830033, 2.5480013506287813]`、同一 `beta_star_interval`、同样的
  `adaptive_splits` = 2116；
* `multires_bound.py` 输出与生产结果逐字段相同，**唯一例外**是正文未列表的宿主
  `Exp(mean=B)` 的 `beta_hat`（4.187529 对 4.187528）——差在第六位小数，属 SciPy 有界
  搜索容差范围内的版本差异；正文引用的全部数值不受影响；
* `test_multires.py` 通过（须先跑 `certify_delta_plus.py` 生成证书）；
* 图 2 重绘流程按下节命令跑通，切片与投稿包中的 `fig2_data.json` 完全一致。

## 运行

    python3 run_all.py          # 快速冒烟（分钟级），写 result/*_quick.json
    python3 run_all.py --full   # 正式全量，写 result/*.json 与 fig/*.pdf

正式档各步大致耗时（2 核参考值）：

| 步骤 | 耗时 |
|---|---|
| `certify_delta_plus.py` | 约 1 分钟 |
| `multires_bound.py` + `lambda_delta_minus.py` + `make_figs.py` | 约 10 分钟 |
| `correlation_full.py` | 约 40–60 分钟（九行，2×10⁶ 步 × 3 种子，含网格与路径长度自检）|
| `delta_minus_width_bounds.py` | 约 5 分钟（依赖 `result/correlation_full.json`）|
| `gmi_demo.py` | 约 15 分钟 |
| `gmi_finite.py` | 约 20–30 分钟 |
| `run_ding.py` | 约 5–10 分钟 |

## 脚本 → 结果 → 正文

| 脚本 | 结果 | 支撑正文 |
|---|---|---|
| `multires_bound.py` | `result/multires.json` | 独立宿主的固定设计代价 $J_\beta$、有界搜索候选点 $\widehat\beta$ 与所报界 $J_{\widehat\beta}$，以及三幅代码图件 |
| `certify_delta_plus.py` | `result/delta_plus_certificate.json` | 表 2 三个独立宿主的全局 $\Delta_+$ 证书；指数主例用解析端点排除加高精度区间分支定界，另两例用解析证书 |
| `lambda_delta_minus.py` | `result/lambda_delta_minus.json` | 表 1 走廊指数与下端点 $\Delta_-$（梯形规则）|
| `correlation_full.py` | `result/correlation_full.json` | 表 2 九行（中点规则）、配对标准误、网格与路径自检 |
| `delta_minus_width_bounds.py` | `result/delta_minus_width_bounds.json` | $\Delta_-$ 单字母窗口、包含性 15/15、紧度 |
| `gmi_demo.py` | `result/gmi_demo.json` | 三个“边缘—分层比”设计的常数 $c_1$ 与上界 |
| `gmi_finite.py` | `result/gmi_finite.json` | 有限噪声端 GMI 下界的数值估计 |
| `run_ding.py` | `result/reset_control.json` | 清零概率对相关性敏感度的对照表 |
| `make_figs.py` | `fig/*.pdf` | 中文母本的图 2–4 |
| `make_fig2.py` | `fig_multires_gap.pdf` | 英文稿实际排印的图 2（见下）|

## 论文图 2 的重绘

英文稿排印的图 2 由 `make_fig2.py` 绘制，它只改标签、字体与画布尺寸，不重算任何数值。
它读取 `fig2_data.json`——`result/multires.json` 中 `Exp(mean=2B)` 键的逐字切片。
跑完正式档后：

    python3 -c "import json;d=json.load(open('result/multires.json'));json.dump({'source':'result/multires.json key \"Exp(mean=2B)\"','Exp(mean=2B)':d['Exp(mean=2B)']},open('fig2_data.json','w'),indent=1)"
    python3 make_fig2.py        # -> fig_multires_gap.pdf

重绘结果与投稿包中的 PDF 除嵌入的 `/CreationDate` 时间戳外完全一致。

## 关于 `check_consistency.py`

该脚本是“正文每个数字与 JSON 一致”的机器校验，为完整性收录于此，但**在本仓库内无法
独立运行**：它要读作者内部的稿件树（`../manuscript/main.tex`、
`../manuscript/appendix_proofs.tex`、`../sketelon.md`、`../design.md`、
`../manuscript/supplement/numerical_evidence.md`），这些不随本仓库发布。核对论文数字
请用 `MANIFEST.md`，那里逐条写明了对应关系与应得的数值。

## 方法学说明与证据等级

Python 层与 JSON schema 用 `J_beta` 表示固定 $\beta$ 的构造代价，用 `beta_hat` 与
`J_beta_hat` 表示 $\beta\in[1.02,20]$ 上由 SciPy 有界一维搜索得到的候选点及其固定设计
代价，用 `Delta_minus` 表示走廊损失下端。该局部搜索**不**以 `beta_star` 或 `Delta_plus`
标记输出；理论端点仍为 $\Delta_+=\inf_{\beta>1}J_\beta$。

论文严格区分“已认证”与“估计”，代码保持同一口径：

* **已认证。** 表 2 三个独立宿主的 $\Delta_+$。指数宿主用解析端点排除加概率包络分支
  定界，在 80 位算术下覆盖全部 $\beta>1$，证书宽度 < 5×10⁻⁵；两个均匀宿主为零宽解析
  区间。
* **确定性估计，不是证书。** 走廊指数以及由它得到的 $\Delta_-$ 是体积指数的固定网格
  离散，不是连续算子证书；网格敏感度写入输出。
* **蒙特卡洛估计。** GMI 常数 $c_1$ 与有限噪声 GMI 数值是带标准误的蒙特卡洛估计，
  正文按此措辞引用。
* **仅固定设计界。** `J_beta_hat` 只是候选点上可证的固定设计界；相关宿主的 `J_sel`
  同样是固定设计界，不属于全局认证范围。

随机性与自检：`gmi_demo.py` 的平稳滞留池为轨迹的确定性等距稀释，评估用策略字母由独立
随机流抽取；`gmi_finite.py` 的评估滞留取自独立轨迹（seed+7000）。`correlation_full.py`
内置网格（800/1600/3200）与路径长度自检；`gmi_demo.py` 内置共用随机流的配对扰动自检；
`gmi_finite.py` 内置参考密度池与混合样本的嵌套误差自检。

## 许可

MIT，见 `LICENSE`。
