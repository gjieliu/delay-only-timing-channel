# 数值复现包

论文《平稳遍历宿主下只延后时序信道容量界》全部数值结果的生产代码。
正文引用的每一个数字都要求由 `result/` 下的结构化 JSON 经固定格式化后逐字出现，
由 `check_consistency.py` 机器校验。
本文件中的命令默认从 `ver8.1/code/` 运行；正式结果与图件分别写入同级的 `result/` 和 `fig/`，论文源码位于相邻的 `../manuscript/`。复现包保留全部正式计算；主文最终只引用 `sketelon.md` 固定的最小证据集，其余结果转入补充验证层。

## 依赖

需要 Python ≥ 3.9、NumPy ≥ 1.22、SciPy（用于分层比的一维优化）与 mpmath（用于表 2 的高精度全局证书）；
`make_figs.py` 另需 Matplotlib。
正式结果由 Python 3.11.15 + NumPy 2.4.4 生成；全部脚本使用
`numpy.random.default_rng` 的显式种子，同版本 NumPy 下逐位可复现。

## 快速冒烟（分钟级）

    python3 run_all.py

随机数值脚本以 `--quick` 档运行并写 `result/*_quick.json`，不覆盖正式 JSON；
确定性的多分辨率基线与图件会按相同配置重建。

## 正式复现（全量）

    python3 run_all.py --full        # 表 1、表 2、清零对照、GMI 与图件全量重算
    python3 check_consistency.py     # 应输出 consistency: ok

正式档各步大致耗时（2 核参考值）：`certify_delta_plus.py` 约 1 分钟，`multires_bound`/`lambda_delta_minus`/`make_figs`
共约 10 分钟；`correlation_full.py` 约 40–60 分钟（相关宿主表九行，2×10^6 步 × 3 种子，
含双 θ 网格自检与路径长度自检）；`delta_minus_width_bounds.py` 约 5 分钟（依赖
`correlation_full.json`）；`gmi_demo.py` 约 15 分钟（常数 c1，含配对收敛自检）；
`gmi_finite.py` 约 20–30 分钟（有限噪声 GMI 检查，含嵌套误差自检）；
`run_ding.py` 约 5–10 分钟（清零概率对照，四个独立种子）。

## 脚本 → 结果 → 正文对应

| 脚本 | 结果 | 支撑正文 |
|---|---|---|
| `multires_bound.py` | `multires.json` | 独立宿主表的固定设计代价 $J_\beta$、有界搜索候选点 $\widehat\beta$、所报界 $J_{\widehat\beta}$ 与三幅代码图件（字段 `J_beta`、`beta_hat` 与 `J_beta_hat`）|
| `certify_delta_plus.py` | `delta_plus_certificate.json` | 表 2 三个独立宿主的全局 $\Delta_+$ 证书；指数主例采用解析端点排除与高精度区间分支定界，另两个主例采用解析证书 |
| `lambda_delta_minus.py` | `lambda_delta_minus.json` | 表 1 走廊指数与论文下端点 $\Delta_-$（字段 `Delta_minus`，梯形规则）|
| `correlation_full.py` | `correlation_full.json` | 表 2 九行（中点规则）、配对标准误、网格/路径自检 |
| `delta_minus_width_bounds.py` | `delta_minus_width_bounds.json` | $\Delta_-$ 单字母窗口、包含性 15/15、紧度 |
| `gmi_demo.py` | `gmi_demo.json` | 三个“边缘—分层比”设计的常数 c1 与上界（含最优分层比对照）|
| `gmi_finite.py` | `gmi_finite.json` | 有限噪声端 GMI 下界的数值估计 |
| `run_ding.py` | `reset_control.json` | 清零概率对相关性敏感度的丁族对照表 |
| `make_figs.py` | `fig/*.pdf` | 图 2–4 |

## 方法学说明

Python 层与 JSON schema 使用 `J_beta` 表示固定 $\beta$ 的构造代价，使用
`beta_hat` 与 `J_beta_hat` 表示在 $\beta\in[1.02,20]$ 上由 SciPy 有界一维搜索得到的候选点及其固定设计代价，使用 `Delta_minus` 表示走廊损失下端。该局部搜索不以 `beta_star` 或 `Delta_plus` 标记输出；理论端点仍为 $\Delta_+=\inf_{\beta>1}J_\beta$。表 2 的三个独立主例由独立脚本输出 `Delta_plus_interval`：指数宿主证书覆盖全部 $\beta>1$ 且宽度小于 $5\times10^{-5}$，有界均匀宿主与满宽度宿主给出零宽解析区间。相关宿主的 `J_sel` 仍是固定设计界，不属于该全局认证范围。

* 随机性控制：所有随机流均由显式种子派生；`gmi_demo.py` 的平稳滞留池为轨迹的
  确定性等距稀释，评估用策略字母由独立随机流抽取；`gmi_finite.py` 的评估滞留
  取自独立轨迹（seed+7000）。
* 自检：`correlation_full.py` 内置网格（800/1600/3200）与路径长度自检；
  `gmi_demo.py` 内置共用随机流的配对扰动自检（池/网格/评估点/密度池半倍与双倍）；
  `gmi_finite.py` 内置参考密度池与混合样本的嵌套误差自检。
* 证据等级：c1 与有限噪声 GMI 数字为 Monte Carlo 估计（正文按此措辞引用），
  未构成确定性认证区间；走廊指数为固定网格离散估计，不是连续算子证书；
  `J_beta_hat` 单独看只是候选点上的可证固定设计界；只有 `delta_plus_certificate.json` 所列的三个独立宿主具有全局 $\Delta_+$ 证书。

已废弃的早期脚本归档在项目根目录 `_to_detelet/versions/ver1-pre-tit/code/`，不参与当前生产链。
