"""Makes the report/slide figures from results/ into figures/."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import venus_env as ve
from thermal_model import COATINGS, SIG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
R = json.load(open(os.path.join(RES, "results.json")))

NAVY, ORANGE, BLUE, GREEN, GRAY = "#1F3B57", "#C2571A", "#3C7DB8", "#4F8A3E", "#7A828C"
NODE_COLORS = {"Electronics": NAVY, "Radiator": BLUE, "Battery": ORANGE, "Propellant tank": GREEN, "Bus": GRAY}
plt.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.labelsize": 13, "xtick.labelsize": 12,
    "ytick.labelsize": 12, "legend.fontsize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.8, "axes.grid": True, "grid.color": "#E3E6EA",
    "grid.linewidth": 0.8, "legend.frameon": False, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(os.path.join(FIG, name))
    plt.close(fig)


# eclipse duration vs beta
b, e = np.array(R["eclipse_vs_beta"]).T
fig, ax = plt.subplots(figsize=(6.4, 3.9))
ax.plot(b, e, color=NAVY, lw=2.4, label="250 km circular (gravity orbit)")
b_star = np.degrees(np.arcsin(ve.R_V / (ve.R_V + 250e3)))
ax.axvspan(b_star, 90, color="#F3E3D8", lw=0)
ax.text((b_star + 90) / 2, 30, "No eclipse\n(full sun)", ha="center", color=ORANGE, fontsize=12)
ax.plot(0, R["eclipse_ellip_max_min"], "o", color=ORANGE, ms=8, label="294 x 8,540 km elliptical, worst case")
ax.annotate(f"{R['eclipse_ellip_max_min']:.1f} min", (0, R["eclipse_ellip_max_min"]), (8, 55), color=ORANGE, fontsize=12)
ax.annotate(f"{e[0]:.1f} min", (0, e[0]), (6, e[0] - 8), color=NAVY, fontsize=12)
ax.annotate(f"β* = {b_star:.1f}°", (b_star, 0), (b_star - 21, 6), color=GRAY, fontsize=12)
ax.set_xlabel("Beta angle β (deg)")
ax.set_ylabel("Eclipse duration (min)")
ax.set_xlim(-2, 90); ax.set_ylim(0, 65)
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02))
save(fig, "eclipse_vs_beta.png")

# absorbed flux on the -Y radiator over one orbit
fig, ax = plt.subplots(figsize=(6.6, 3.9))
a_eol, eps = COATINGS["OSR"][1], COATINGS["OSR"][2]
for beta, col in ((0, ORANGE), (30, NAVY), (90, BLUE)):
    f = ve.face_fluxes(250e3, 250e3, beta, n=360, S=ve.S_PERI, albedo=ve.ALB_HOT, q_ir=ve.IR_HOT)
    q = a_eol * (f["q_sun"]["anti"] + f["q_alb"]["anti"]) + eps * f["q_ir"]["anti"]
    ax.plot(f["t"] / 60, q, color=col, lw=2.2, label=f"β = {beta}°  (orbit avg {q.mean():.0f} W/m²)")
q_emit = eps * SIG * (273.15 + 15) ** 4
ax.axhline(q_emit, color=GRAY, ls="--", lw=1.2)
ax.text(46, q_emit + 14, f"OSR emission at 15 °C ≈ {q_emit:.0f} W/m²", color=GRAY, ha="center", fontsize=12)
ax.set_xlabel("Time from orbit noon (min)")
ax.set_ylabel("Absorbed flux on radiator (W/m²)")
ax.set_xlim(0, 92); ax.set_ylim(0, 360)
ax.legend(loc="center", ncol=1, bbox_to_anchor=(0.5, 0.58))
save(fig, "radiator_absorbed_flux.png")

# hot case vs beta
S = R["sizing"]
sw = np.array(S["sweep"])
fig, ax = plt.subplots(figsize=(6.6, 4.0))
for i, nm in ((0, "Electronics"), (3, "Propellant tank"), (2, "Battery"), (1, "Radiator")):
    ax.plot(S["betas"], sw[:, i], "o-", color=NODE_COLORS[nm], lw=2.2, ms=5, label=nm)
for val, col in ((39, NAVY), (29, GREEN), (24, ORANGE)):
    ax.axhline(val, color=col, ls="--", lw=1.2)
ax.plot([], [], ls="--", color=GRAY, lw=1.2, label="Max limit − 11 °C (matching color)")
ax.axvline(S["beta_hot"], color=GRAY, lw=1, ls=":")
ax.text(S["beta_hot"] + 1.5, 41, f"Hot case β = {S['beta_hot']}°", color="#333333", fontsize=11)
ax.set_xlabel("Beta angle β (deg)")
ax.set_ylabel("Orbit-average temperature (°C)")
ax.set_xlim(-3, 92); ax.set_ylim(0, 44)
ax.legend(loc="lower left", ncol=2, fontsize=10.5)
save(fig, "hot_case_beta_sweep.png")

# heater demand by cold case
names = list(R["cold"].keys())
Q = np.array([R["cold"][n] for n in names])
order = np.argsort(Q.sum(1))
fig, ax = plt.subplots(figsize=(6.8, 4.0))
left = np.zeros(len(names))
for j, (lab, col) in enumerate((("Electronics panel", NAVY), ("Battery", ORANGE), ("Propellant tank", GREEN))):
    ax.barh(np.arange(len(names)), Q[order, j], left=left, color=col, height=0.6, label=lab)
    left += Q[order, j]
for i, v in enumerate(left):
    ax.text(v + 3, i, f"{v:.0f} W", va="center", fontsize=11)
ax.set_yticks(np.arange(len(names)), [names[k] for k in order])
ax.set_xlabel("Required heater power, orbit average (W)")
ax.set_xlim(0, 240)
ax.grid(axis="y", visible=False)
ax.legend(loc="lower right")
save(fig, "heater_demand.png")


# transients
def transient_plot(key, fname):
    d = np.loadtxt(os.path.join(RES, "transients", f"{key}.csv"), delimiter=",", skiprows=1)
    t = d[:, 0]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[3.2, 1]))
    ecl = d[:, 9] > 0.5
    for a in (ax, ax2):
        a.fill_between(t, 0, 1, where=ecl, color="#DDE2E8", transform=a.get_xaxis_transform(), lw=0)
    for i, nm in enumerate(("Electronics", "Radiator", "Battery", "Propellant tank", "Bus")):
        ax.plot(t, d[:, 1 + i], color=NODE_COLORS[nm], lw=2.0, label=nm)
    ax.set_ylabel("Temperature (°C)")
    ax.legend(loc="upper center", ncol=5, fontsize=10, bbox_to_anchor=(0.5, 1.16), handlelength=1.5)
    for j, (lab, col) in enumerate((("EL heater", NAVY), ("BAT heater", ORANGE), ("PROP heater", GREEN))):
        on = d[:, 6 + j] > 0.5
        ax2.fill_between(t, 2 - j - 0.35, 2 - j + 0.35, where=on, color=col, lw=0, step="post")
    ax2.set_yticks([2, 1, 0], ["EL", "BAT", "PROP"])
    ax2.set_ylim(-0.6, 2.6)
    ax2.set_ylabel("Heater ON")
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel("Time (min, last two orbits; shaded = eclipse)")
    ax2.set_xlim(t[0], t[-1])
    save(fig, fname)


transient_plot("cold_op_eclipse", "transient_eclipse_cold_op.png")
transient_plot("safe_58min_eclipse", "transient_safe_58min_eclipse.png")
transient_plot("hot_op", "transient_hot_op.png")

# trade map
fig, ax = plt.subplots(figsize=(5.4, 4.1))
style = {"T1 radiator finish": "o", "T2 radiator location": "s", "T3 cold-case control": "^"}
short = {"Silvered Teflon": ("Silvered Teflon", 0.04, 0, "left"), "Z93 white paint": ("Z93 paint", 0.0, 20, "center"),
         "Zenith (-Z)": ("Zenith", 0.04, 0, "left"), "Ram (+X)": ("Ram", 0.05, 12, "left"),
         "Oversized (1.5x) radiator + heaters": ("Oversized", 0.05, -16, "left"),
         "Louvered radiator + heaters": ("Louvers", 0.04, 0, "left")}
for trade, rows in R["trades"].items():
    for r in rows:
        if r["concept"].startswith("A:"):
            continue
        lab, dx, dy, ha = short[r["concept"][3:]]
        ax.plot(r["A_rad"], r["heater_surv"], style[trade], ms=10, color=NAVY, mfc="white", mew=2)
        ax.text(r["A_rad"] + dx, r["heater_surv"] + dy, lab, fontsize=12, color="#333333", va="center", ha=ha)
sel = R["trades"]["T1 radiator finish"][0]
ax.plot(sel["A_rad"], sel["heater_surv"], "o", ms=12, color=ORANGE)
ax.text(sel["A_rad"] + 0.03, sel["heater_surv"] - 34, "Selected", fontsize=12, color=ORANGE, fontweight="bold")
for trade, m in style.items():
    ax.plot([], [], m, color=NAVY, mfc="white", mew=2, ms=9, label=trade.split(" ", 1)[1].capitalize())
ax.legend(loc="lower right", fontsize=11)
ax.set_xlabel("Radiator area (m²)")
ax.set_ylabel("Survival heater power (W)")
ax.set_xlim(1.1, 2.2); ax.set_ylim(0, 360)
save(fig, "trade_map.png")

# margins: allowable vs predicted ranges
M = R["margins"]
labels = {"EL": "Electronics", "BAT": "Battery", "PROP": "Propellant tank"}
fig, ax = plt.subplots(figsize=(6.8, 3.4))
for i, k in enumerate(("PROP", "BAT", "EL")):
    lo, hi, slo, shi = M[k]["limits"]
    ax.barh(i, shi - slo, left=slo, height=0.62, color="#E6E9ED", label="Survival range" if i == 0 else None)
    ax.barh(i, hi - lo, left=lo, height=0.62, color="#B9C6D3", label="Operating range" if i == 0 else None)
    ax.barh(i, M[k]["hot"] - M[k]["cold"], left=M[k]["cold"], height=0.26, color=ORANGE,
            label="Predicted (cold op to hot op)" if i == 0 else None)
    ax.text(hi + 1, i + 0.12, f"+{M[k]['m_hot']:.1f} °C", fontsize=11, color=NAVY, va="center")
ax.set_yticks([0, 1, 2], [labels[k] for k in ("PROP", "BAT", "EL")])
ax.set_xlabel("Temperature (°C)")
ax.set_xlim(-20, 72)
ax.grid(axis="y", visible=False)
ax.legend(loc="upper center", ncol=3, fontsize=10.5, bbox_to_anchor=(0.5, 1.2))
save(fig, "margin_ranges.png")
print("figures written to", FIG)
