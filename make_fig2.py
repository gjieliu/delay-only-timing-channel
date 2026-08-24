"""Fig. 2 (English) — gap to the full-width modulo benchmark.

Reads fig2_data.json, a verbatim slice of the official production result
ver8.1/code/result/multires.json (key "Exp(mean=2B)").  No numeric value is
recomputed or altered here: this script only re-renders the figure with
English labels and IEEE-style typography.

Differences from the Chinese producer ver8.1/code/make_figs.py::fig_gap
  * labels/legend in English;
  * the in-figure title is removed -- in IEEE style the caption carries it,
    and the interpretation now lives in the body text of Section V;
  * serif (STIX) fonts to match IEEEtran's Times;
  * canvas 5.0x3.3 in instead of 6.2x4.1 in, so that at width=0.70\\linewidth
    the rendered height is unchanged (3.00 in) while the labels come out
    about 12% larger.

Run:  python3 make_fig2.py      ->  fig_multires_gap.pdf
"""

import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
})

with open("fig2_data.json", encoding="utf-8") as f:
    row = json.load(f)["Exp(mean=2B)"]

snr = np.array(sorted(int(g) for g in row["curves"]))
full = np.array([row["curves"][str(g)]["full_width_modulo"] for g in snr])
multi_gap = full - np.array([row["curves"][str(g)]["multires"] for g in snr])
single_gap = full - np.array([row["curves"][str(g)]["single_threshold"] for g in snr])
j_hat = row["partition"]["J_beta_hat"]

fig, ax = plt.subplots(figsize=(5.0, 3.3))
ax.plot(snr, multi_gap, "C0-o", ms=4.0, lw=1.3, label="Multiresolution construction")
ax.plot(snr, single_gap, "C3-s", ms=3.6, lw=1.3, label="Optimal single threshold")
ax.axhline(j_hat, color="C0", ls="--", lw=1.1,
           label=rf"Fixed design cost $J_{{\widehat{{\beta}}}}={j_hat:.3f}$")
ax.set_xscale("log", base=2)
ax.set_xlabel(r"Timing signal-to-noise ratio $\mathrm{SNR}=B/\sigma_\xi$")
ax.set_ylabel("Gap to full-width benchmark (bits/packet)")
ax.grid(True, which="both", ls=":", alpha=0.5)
ax.legend(loc="upper left", framealpha=0.9)
fig.tight_layout(pad=0.4)
fig.savefig("fig_multires_gap.pdf")
plt.close(fig)
print("wrote fig_multires_gap.pdf  (J_beta_hat = %.6f)" % j_hat)
