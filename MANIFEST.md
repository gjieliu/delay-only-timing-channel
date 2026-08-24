# Manifest — every reported number, and the run that produces it

This repository tracks source only. The files named below (`result/*.json`, `fig/*.pdf`)
do not exist until you run `python3 run_all.py --full`; the **Expected value** column is
what that run should give, so this table doubles as a verification checklist.

Field paths use `/` for nesting inside the JSON, e.g. `Exp(mean=2B)/partition/J_beta_hat`.

Two naming conventions need stating once, because the code was written in Chinese and the
paper in English:

* **Host family labels.** Correlated-host families are keyed by Chinese ordinals:
  `jia` / 甲 = **family A** (marginal `Unif(B/2, 3B/2)`, contains full-width resets with
  probability 1/2), `yi` / 乙 = **family B** (marginal `Unif(B/4, 3B/4)`, no resets),
  `bing` / 丙 = a third family carried by the code that the paper does not tabulate.
* **Symbols.** `J_beta` is the fixed-design tiering cost at a given ratio; `beta_hat` and
  `J_beta_hat` are the candidate ratio from the bounded 1-D search and its provable
  fixed-design bound; `Delta_minus` is the corridor lower endpoint. The local search
  deliberately does **not** emit `beta_star` or `Delta_plus` — those names appear only in
  `delta_plus_certificate.json`, which is the globally certified quantity.

---

## Table II — Numerical bracket on the capacity loss, three independent hosts

Scripts: `multires_bound.py` (the `beta_hat` column), `certify_delta_plus.py` (the
certified `Delta_+`), `lambda_delta_minus.py` (the `Delta_-` column).

| Paper cell | Printed | Output file | Field | Expected value |
|---|---|---|---|---|
| Exp(mean 2B), `beta_hat` | 3.837 | `result/multires.json` | `Exp(mean=2B)/partition/beta_hat` | 3.837144 |
| Exp(mean 2B), `Delta_+` | 2.548 | `result/delta_plus_certificate.json` | `hosts/Exp(mean=2B)/Delta_plus_interval` | [2.5479514302830033, 2.5480013506287813] |
| Exp(mean 2B), `Delta_-` | 0.143 | `result/lambda_delta_minus.json` | `Exp(mean=2B)/Delta_minus` | 0.143 |
| Unif(B/2,3B/2), `beta_hat` | 2.000 | `result/delta_plus_certificate.json` | `hosts/U(B/2,3B/2)/beta_star_interval` | [2.0, 2.0] |
| Unif(B/2,3B/2), `Delta_+` | 1.500 | same | `hosts/U(B/2,3B/2)/Delta_plus_interval` | [1.5, 1.5] |
| Unif(B/2,3B/2), `Delta_-` | 0.032 | `result/lambda_delta_minus.json` | `U(B/2,3B/2)/Delta_minus` | 0.032 |
| Unif(B,2B), `beta_hat` | any | `result/delta_plus_certificate.json` | `hosts/U(B,2B) (sparse)/beta_star` | "any beta>1" |
| Unif(B,2B), `Delta_+` | 0.000 | same | `hosts/U(B,2B) (sparse)/Delta_plus_interval` | [0.0, 0.0] |
| Unif(B,2B), `Delta_-` | 0.000 | `result/lambda_delta_minus.json` | `U(B,2B) (sparse)/Delta_minus` | 0.0 |

In-text statements attached to the same table, all from
`result/delta_plus_certificate.json`, host `Exp(mean=2B)`:

| Paper statement | Field | Expected value |
|---|---|---|
| `2.54795143 <= Delta_+ <= 2.54800135` | `Delta_plus_interval` | [2.5479514302830033, 2.5480013506287813] |
| certificate width < 5e-5 | `certificate_width`, `target_width` | 4.992034577808031e-05, 5e-05 |
| minimizer localized to `beta* in [3.533, 4.142]` | `beta_star_interval` | [3.533245203248791, 4.142451740767473] |
| two grouping bounds confine the infimum to `beta in [1.05, 17]` | `compact_beta_interval`, `endpoint_lower_bounds` | [1.05, 17.0]; `beta_le_1.05` = 2.5592309664180948, `beta_ge_17` = 2.575293286002973 — both above the certified upper end 2.54800135, which is what excludes the two non-compact ends |
| 80-digit interval branch-and-bound (search parameters listed in Appendix F) | `mp_dps`, `initial_cells`, `adaptive_splits`, `round_guard`, `tail_exponent` | 80, 4096, 2116, 1e-12, 38.0 |
| certification route per host | `hosts/*/method` | Exp: "analytic endpoint exclusion plus probability-enclosure branch-and-bound"; U(B/2,3B/2): "analytic partition-entropy inequality"; U(B,2B): "analytic K_beta=0 almost surely" |
| loss bracket about [0.143, 2.548] | derived from the two columns above | — |

Grid insensitivity of the `Delta_-` column, quoted in the text as orders of magnitude
below the effects discussed: `result/lambda_delta_minus.json`,
`_meta/grid_convergence_check` — grids 300/600/1200 give 0.14711318 / 0.14711283 /
0.14711274, `max_span` = 4.4e-07.

## Figure 2 — achievable gap to the full-width modulo-circle benchmark

Data: `result/multires.json`, key `Exp(mean=2B)`, subkey `curves`, over `snr = B/sigma`
from 2^3 to 2^18; each entry carries `full_width_modulo`, `multires`, `single_threshold`,
`single_b_over_B`.

* `make_figs.py` draws the Chinese master's version into `fig/`.
* The figure printed in the English manuscript is drawn by `make_fig2.py` from a verbatim
  slice of the same key (the README gives the one-line extraction command). Only labels,
  fonts and canvas size differ; no value is recomputed. The result matches the PDF in the
  submission package except for the embedded `/CreationDate` timestamp.
* The asymptote quoted in the text, `J_{beta_hat} = 2.548` at `beta_hat = 3.837`:
  `Exp(mean=2B)/partition/J_beta_hat` = 2.548001, `.../beta_hat` = 3.837144.
* The decomposition quoted in Appendix F, `H = 1.459724` and `E[K] log2 beta = 1.088277`
  summing to 2.548001: `.../label_entropy` and `.../resolution_loss`.

## Section V — GMI comparison

Script `gmi_demo.py` → `result/gmi_demo.json`, key `families`.

| Paper statement | Field | Expected value |
|---|---|---|
| Unif(B/2,3B/2), beta=2: `c_1 = 0.4992` | `Unif(B/2,3B/2), beta=2/c1_mean` | 0.49922916646747856 |
| inter-replication range 0.0005 | `.../c1_range` | 0.0004766993565352706 |
| largest paired sensitivity shift 0.0038 | `convergence_check/max_shift` | 0.003804068928274784 |
| GMI point estimate 0.999 | `.../tightened_bound` | 0.9992291664674786 |
| compared against `J_2 = 1.500` | `.../J_beta` | 1.5 |
| negative control Unif(B/4,3B/4), beta=4: `c_1 = 0.0014` | `Unif(B/4,3B/4), beta=4/c1_mean` | 0.0014146730972927624 |
| range 0.0066 | `.../c1_range` | 0.006572394192160158 |
| auxiliary-channel estimate 2.001 | `.../tightened_bound` | 2.0014146730972926 |
| `J_4 = 2.000` | `.../J_beta` | 2.0 |
| three independent replications | `meta/seeds`, `meta/design_count` | seeds begin at 31; 3 designs |

Finite-noise GMI check (`gmi_finite.py` → `result/gmi_finite.json`): six points over seeds
{41, 42} and sigma in {0.04, 0.02, 0.01}, each with `gmi`, `gmi_se`, `G_bound`,
`G_bound_se`; nesting self-check `nesting_check/max_shift` = 0.0352.

## Table III(a) — Host memory and full-width resets, fixed marginal

Script `correlation_full.py` → `result/correlation_full.json`, key `families`; rows
selected by `key` and `theta`.

| Paper row | Printed `Delta_-` | Field | Expected value | `J_sel` |
|---|---|---|---|---|
| A, theta=0.50 | 0.0320 | `families[jia]/rows[theta=0.5]/Delta_minus_mean` | 0.03201007743085629 | `.../J_beta` = 1.5 |
| A, theta=0.99 | 0.0329 | `families[jia]/rows[theta=0.99]/Delta_minus_mean` | 0.03292839785069049 | `.../J_beta` = 0.5807931358959112 → 0.581 |
| B, theta=0.50 | 0.2626 | `families[yi]/rows[theta=0.5]/Delta_minus_mean` | 0.2625962533840501 | see note |
| B, theta=0.99 | 0.2821 | `families[yi]/rows[theta=0.99]/Delta_minus_mean` | 0.2820977276204561 | `.../J_beta` = 1.5807931358959113 → 1.581 |

**Note on `J_sel` in the family-B, theta=0.50 row.** The paper defines `J_sel` as the
smaller fixed tiering bound over the reported design set {beta=2, beta=4}. At theta=0.50
family B has `J_2 = 2.500` (`families[yi]/rows[theta=0.5]/J_beta`) but `J_4 = 2.000`,
which is smaller; `J_4` does not depend on the tier chain and comes from
`result/gmi_demo.json`, `families/Unif(B/4,3B/4), beta=4/J_beta` = 2.0. Hence the
tabulated 2.000 and width 2.000 − 0.2626 = 1.737. At theta=0.99, `J_2 = 1.581` is the
smaller of the two — which is why the design attaining the minimum changes with theta.

Movements quoted in the text: family A rises by 0.0009 with paired standard error 0.00008
(`families[jia]/Delta_minus_change` = 0.0009183204198341954,
`.../Delta_minus_change_paired_se` = 8.312547541962445e-05); family B rises by 0.0195 with
paired standard error about 0.00039 (same two fields under `families[yi]`:
0.01950147423640597 and 0.00038951669492009176).

Self-checks quoted in the text: `grid_checks` (800/1600/3200 grid points) and `path_check`
(500k/1M/2M steps); `max_spread` = 0.0025675784273097024.

## Table III(b) — Fixed constrained-step scale, varying reset probability

Script `run_ding.py` → `result/reset_control.json`, key `rows`, indexed by `p`. The run
covers `p` in {0.00, 0.10, 0.25, 0.50}; the paper tabulates the two extreme rows.

| Paper cell | Printed | Field | Expected value |
|---|---|---|---|
| p=0.00, `Delta_-(.50)` | 0.2626 | `rows[p=0.0]/theta/0.5/mean` | 0.2625833451756345 |
| p=0.00, `Delta_-(.99)` | 0.2804 | `rows[p=0.0]/theta/0.99/mean` | 0.280420441462732 |
| p=0.00, `s(p)` | 0.01784 | `rows[p=0.0]/sensitivity` | 0.01783709628709751 |
| p=0.00, ratio | 1.00 | `rows[p=0.0]/ratio` | 1.0 (equal by construction) |
| p=0.50, `Delta_-(.50)` | 0.1182 | `rows[p=0.5]/theta/0.5/mean` | 0.11820307071860486 |
| p=0.50, `Delta_-(.99)` | 0.1210 | `rows[p=0.5]/theta/0.99/mean` | 0.12103301822111909 |
| p=0.50, `s(p)` | 0.00283 | `rows[p=0.5]/sensitivity` | 0.0028299475025142275 |
| p=0.50, `s(0)(1-p)` | 0.00892 | `rows[p=0.5]/mass_dilution_prediction` | 0.008918548143548755 |
| p=0.50, ratio | 0.32 | `rows[p=0.5]/ratio` | 0.31731033537799247 |
| max standard error on the ratio | 0.0014 | `max_standard_error` | 0.0013491717429663956 |
| four independent seeds | — | `meta/seeds`, `meta/repetitions` | seeds begin at 1; 4 repetitions |

**Cross-run agreement quoted in the text.** Parts (a) and (b) come from two independent
production runs with different grids, path lengths and random streams. The paper reports
0.2804 (b) against 0.2821 (a) at theta=0.99, and 0.26258 against 0.26260 at theta=0.50 —
the latter pair being `rows[p=0.0]/theta/0.5/mean` = 0.2625833451756345 and
`families[yi]/rows[theta=0.5]/Delta_minus_mean` = 0.2625962533840501, which round to the
same four decimals. These are genuinely two separate runs, not one number printed twice.

## Single-letter window — "all 15 estimates fall inside the bracket"

Script `delta_minus_width_bounds.py` → `result/delta_minus_width_bounds.json`:
`checks_total` = 15, `checks_ok` = 15, `all_contained` = true. Per-family windows are in
`families[*]/Delta_minus_lo` and `.../Delta_minus_hi`, with per-row containment and
tightness percentages in `families[*]/rows`. This script consumes `correlation_full.json`
(`meta/correlation_source`) and must be run after it.

## Produced but not cited in the submitted paper

* `result/gap.json` (from `gap_strategy.py`) — earlier gap summary at SNR 8/16/32/48,
  superseded by `multires.json/curves`.
* `result/multires.json` keys `U(0,2B)`, `Exp(mean=B)`, `Exp(mean=B/2)` — additional
  independent hosts; the paper tabulates three. Note that `Exp(mean=B)/partition/beta_hat`
  moves in the sixth decimal across SciPy versions (4.187528 vs 4.187529), which is within
  the tolerance of the bounded 1-D search. No quantity cited in the paper is affected.
* `result/correlation_full.json` family `bing` / 丙 — a third correlated family.
* `fig/fig_multires_rate.pdf`, `fig/fig_partition_penalty.pdf` — figures of the Chinese
  master; the submitted English manuscript carries one figure.
* `corridor.py` (shared corridor kernel), `test_multires.py` (unit tests for the
  multiresolution construction).
