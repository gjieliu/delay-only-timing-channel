"""多分辨率模格主结果出图。"""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "fig")
RESULT = os.path.join(HERE, "result", "multires.json")

font_candidates = (
    glob.glob("/usr/local/texlive/*/texmf-dist/fonts/opentype/public/fandol/"
              "FandolSong-Regular.otf")
    + glob.glob("/usr/share/texlive/texmf-dist/fonts/opentype/public/fandol/"
                "FandolSong-Regular.otf")
)
if font_candidates:
    fm.fontManager.addfont(font_candidates[0])
    plt.rcParams["font.family"] = fm.FontProperties(fname=font_candidates[0]).get_name()
plt.rcParams["axes.unicode_minus"] = False


def load():
    with open(RESULT, encoding="utf-8") as f:
        return json.load(f)


def fig_rate(data):
    row = data["Exp(mean=2B)"]
    gam = np.array([int(g) for g in row["curves"]])
    full = np.array([row["curves"][str(g)]["full_width_modulo"] for g in gam])
    multi = np.array([row["curves"][str(g)]["multires"] for g in gam])
    single = np.array([row["curves"][str(g)]["single_threshold"] for g in gam])

    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    ax.plot(gam, full, "k--", lw=1.5, label="满宽度模环基准")
    ax.plot(gam, multi, "C0-o", ms=4.5, label="多分辨率模格下界")
    ax.plot(gam, single, "C3-s", ms=4.0, label="最优单阈值下界")
    ax.set_xscale("log", base=2)
    ax.set_xlabel(r"时序信噪比 $\mathsf{snr}=B/\sigma_\xi$")
    ax.set_ylabel("速率（比特/包）")
    ax.set_title(r"独立指数宿主 $\mathrm{Exp}(1/(2B))$：多分辨率构造恢复常数间隙")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_multires_rate.pdf"))
    plt.close(fig)


def fig_gap(data):
    row = data["Exp(mean=2B)"]
    gam = np.array([int(g) for g in row["curves"]])
    full = np.array([row["curves"][str(g)]["full_width_modulo"] for g in gam])
    multi_gap = full - np.array([row["curves"][str(g)]["multires"] for g in gam])
    single_gap = full - np.array([row["curves"][str(g)]["single_threshold"] for g in gam])
    j_beta_hat = row["partition"]["J_beta_hat"]

    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    ax.plot(gam, multi_gap, "C0-o", ms=4.5, label="多分辨率绝对差")
    ax.plot(gam, single_gap, "C3-s", ms=4.0, label="单阈值绝对差")
    ax.axhline(j_beta_hat, color="C0", ls="--", lw=1.2,
               label=rf"选定分层界 $J_{{\widehat{{\beta}}}}={j_beta_hat:.3f}$")
    ax.set_xscale("log", base=2)
    ax.set_xlabel(r"时序信噪比 $\mathsf{snr}=B/\sigma_\xi$")
    ax.set_ylabel("相对满宽度模环基准的差（比特/包）")
    ax.set_title("单阈值损失继续增长，多分辨率损失趋于有限上界")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_multires_gap.pdf"))
    plt.close(fig)


def fig_penalty(data):
    names = ["Exp\n均值 $2B$", "$U(0,2B)$", "Exp\n均值 $B$",
             "Exp\n均值 $B/2$", "$U(B/2,3B/2)$", "$U(B,2B)$"]
    keys = ["Exp(mean=2B)", "U(0,2B)", "Exp(mean=B)", "Exp(mean=B/2)",
            "U(B/2,3B/2)", "U(B,2B) (sparse)"]
    label = np.array([data[k]["partition"]["label_entropy"] for k in keys])
    resolution = np.array([data[k]["partition"]["resolution_loss"] for k in keys])
    x = np.arange(len(keys))

    fig, ax = plt.subplots(figsize=(6.5, 4.15))
    ax.bar(x, label, color="C0",
           label=r"尺度层标识损失 $\bar H(\mathbf{K})$（独立时为 $H(K_1)$）")
    ax.bar(x, resolution, bottom=label, color="C1",
           label=r"分辨率量化损失 $\mathrm{E}[K]\log_2\beta$")
    ax.set_xticks(x, names)
    ax.set_ylabel("分层损失界（比特/包）")
    ax.set_title("不同独立宿主基线的多分辨率损失界")
    ax.grid(True, axis="y", ls=":", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_partition_penalty.pdf"))
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    data = load()
    fig_rate(data)
    fig_gap(data)
    fig_penalty(data)
    print("wrote fig/fig_multires_rate.pdf, fig/fig_multires_gap.pdf, "
          "fig/fig_partition_penalty.pdf")


if __name__ == "__main__":
    main()
