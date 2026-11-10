"""Runs model checks, design cases, radiator/heater sizing, transients and trades.

Writes results/summary.md, results/results.json and results/transients/*.csv
"""
import json
import os

import numpy as np

import venus_env as ve
from thermal_model import Case, Design, Network, POWER, K0, SIG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
os.makedirs(os.path.join(RES, "transients"), exist_ok=True)

LIMITS = {"EL": (-5, 50, -15, 60), "BAT": (5, 35, -10, 45), "PROP": (10, 40, 5, 50)}
MARGIN = 11.0            # deg C, passive components (Gilmore 15.6)
HEATER_FACTOR = 1.25     # heater capacity margin
V_MIN, V_NOM, V_MAX = 24.0, 28.0, 32.0
BETAS = (0, 15, 30, 45, 60, 74, 90)
ELLIP = dict(alt_p=294e3, alt_a=8540e3, n_orbit_pts=1440)
HOT = dict(S=ve.S_PERI, albedo=ve.ALB_HOT, q_ir=ve.IR_HOT, eol=True)
COLD = dict(S=ve.S_APH, albedo=ve.ALB_COLD, q_ir=ve.IR_COLD, eol=False)

COLD_CASES = [
    Case("Ops: circular, β = 0°", beta=0, mode="gravity", **COLD),
    Case("Ops: circular, β = 90°", beta=90, mode="gravity", **COLD),
    Case("Ops: elliptical, β = 0°", beta=0, mode="sar", **ELLIP, **COLD),
    Case("Ops: elliptical, β = 90°", beta=90, mode="sar", **ELLIP, **COLD),
    Case("Safe: circular, β = 90°", beta=90, mode="safe", **COLD),
    Case("Safe: elliptical, 58-min eclipse", beta=0, mode="safe", **ELLIP, **COLD),
    Case("Safe: elliptical, β = 90°", beta=90, mode="safe", **ELLIP, **COLD),
    Case("Safe: cruise near 1 AU", S=1367.5, eol=False, mode="cruise", deep_space=True),
]

out, md = {}, []


def hot_case(beta, mode="gravity"):
    return Case(f"Hot op: circular, β = {beta}°", beta=beta, mode=mode, **HOT)


def hot_sweep(design):
    return np.array([Network(design, hot_case(b)).steady()[0] for b in BETAS])


def required_heaters(design, case):
    net = Network(design, case)
    T0 = net.steady()[0]
    on_temp = {"EL": (0, design.set_el[0]), "BAT": (2, design.set_bat[0]), "PROP": (3, design.set_prop[0])}
    active = [n for n, (i, Ts) in on_temp.items() if T0[i] < Ts]
    while active:
        T, Qh, _ = net.steady(tuple(active))
        idx = {"EL": 0, "BAT": 1, "PROP": 2}
        negative = [n for n in active if Qh[idx[n]] < -1e-6]
        if not negative:
            return T, np.clip(Qh, 0, None)
        active = [n for n in active if n not in negative]
    return T0, np.zeros(3)


def size_radiators(base, **changes):
    d = Design(**{**base.__dict__, **changes})
    for A in np.arange(0.60, 2.81, 0.05):
        d.A_rad, d.A_brad = round(float(A), 2), 0.05
        T = hot_sweep(d).max(0)
        if T[0] <= LIMITS["EL"][1] - MARGIN and T[3] <= LIMITS["PROP"][1] - MARGIN:
            break
    else:
        return None
    for Ab in np.arange(0.02, 0.30, 0.01):
        d.A_brad = round(float(Ab), 2)
        if hot_sweep(d).max(0)[2] <= LIMITS["BAT"][1] - MARGIN:
            break
    return d


def table(header, rows):
    md.append("| " + " | ".join(header) + " |")
    md.append("|" + "---|" * len(header))
    md.extend("| " + " | ".join(str(c) for c in r) + " |" for r in rows)
    md.append("")


# ---- 1. model checks --------------------------------------------------------
md.append("# Results summary\n")
md.append("## 1. Model checks\n")
t, r, _ = ve.orbit_states(250e3, 250e3, 0.0, n=4)
F, _ = ve.planet_factors(r, ve.face_normals(r, 0.0))
F_exact = (ve.R_V / (ve.R_V + 250e3))**2
P = ve.circular_period(250e3)
te_exact = P * np.arcsin(ve.R_V / (ve.R_V + 250e3)) / np.pi
te_num = ve.eclipse_duration(250e3, 250e3, 0.0)[0]

R_old, ve.R_V = ve.R_V, 6070e3          # Gilmore Table 2.11 uses R = 6070 km, h = 607 km
sphere = {}
for beta in (0, 90):
    for label, S in (("perihelion", 2759.0), ("aphelion", 2650.0)):
        _, r, _ = ve.orbit_states(607e3, 607e3, beta, n=720)
        _, _, Fs, Fas = ve.planet_factors(r, ve.face_normals(r, beta), R=6070e3, sphere=True)
        q = 0.25 * S * (~ve.in_eclipse(r)) + 0.8 * S * Fas + 153.0 * Fs
        Tk = (q / SIG)**0.25 - K0
        sphere[f"beta {beta} {label}"] = (round(float(Tk.max()), 1), round(float(Tk.min()), 1))
ve.R_V = R_old
gilmore = {"beta 0 perihelion": (122, -105), "beta 0 aphelion": (119, -105),
           "beta 90 perihelion": (67, 67), "beta 90 aphelion": (64, 64)}
table(["Check", "Model", "Reference"], [
    ["Nadir view factor, 250 km", f"{F['nadir'][0]:.5f}", f"(R/r)^2 = {F_exact:.5f}"],
    ["Side-face view factor, 250 km", f"{F['anti'][0]:.4f}", "-"],
    ["Eclipse, 250 km, beta 0 (min)", f"{te_num/60:.2f}", f"P asin(R/r)/pi = {te_exact/60:.2f}"],
    *[[f"Black sphere 607 km, {k}: max / min (C)", f"{v[0]} / {v[1]}", f"Gilmore T2.11: {gilmore[k][0]} / {gilmore[k][1]}"]
      for k, v in sphere.items()],
])
out["checks"] = dict(F_nadir=F["nadir"][0], F_exact=F_exact, F_side=F["anti"][0],
                     eclipse_num_min=te_num / 60, eclipse_exact_min=te_exact / 60,
                     period_min=P / 60, ref_sphere=sphere)

ecl_beta = [(b, ve.eclipse_duration(250e3, 250e3, b)[0] / 60) for b in range(0, 91, 2)]
ecl_ellip = max(ve.eclipse_duration(294e3, 8540e3, 0.0, w)[0] for w in range(0, 360, 5)) / 60
out["eclipse_vs_beta"] = ecl_beta
out["eclipse_ellip_max_min"] = ecl_ellip
md.append(f"Elliptical orbit (294 x 8540 km): period {ve.eclipse_duration(294e3, 8540e3, 0)[1]/3600:.2f} h, "
          f"longest eclipse {ecl_ellip:.1f} min.\n")

# ---- 2. radiator sizing ----------------------------------------------------
md.append("## 2. Radiator sizing (hot case: perihelion, EOL, max albedo/IR)\n")
d_min = size_radiators(Design())
design = Design(**{**d_min.__dict__, "A_rad": round(d_min.A_rad + 0.10, 2)})
sweep = hot_sweep(design)
table(["beta (deg)", "EL", "RAD", "BAT", "PROP", "BUS"],
      [[b, *(f"{v:.1f}" for v in row)] for b, row in zip(BETAS, sweep)])
beta_hot = BETAS[int(sweep[:, 0].argmax())]
net = Network(design, hot_case(beta_hot))
T_hot = net.steady()[0]
q_abs = net.qabs_rad.mean()
A_hand = POWER["gravity"]["el"] / (net.e_rad * SIG * (T_hot[1] + K0)**4 - q_abs)
T_rad_hand = ((POWER["gravity"]["el"] + design.A_rad * q_abs) / (net.e_rad * SIG * design.A_rad))**0.25 - K0
md.append(f"Minimum area {d_min.A_rad:.2f} m2, selected {design.A_rad:.2f} m2 (+0.10 m2 allowance); "
          f"battery radiator {design.A_brad:.2f} m2. Hot-case beta = {beta_hot} deg.\n")
md.append(f"Hand check: q_abs = {q_abs:.1f} W/m2, single-node radiator T = {T_rad_hand:.1f} C "
          f"(model {T_hot[1]:.1f} C); area to reject {POWER['gravity']['el']:.0f} W = {A_hand:.2f} m2.\n")
out["sizing"] = dict(A_rad_min=d_min.A_rad, A_rad=design.A_rad, A_brad=design.A_brad, beta_hot=beta_hot,
                     betas=list(BETAS), sweep=sweep.tolist(), q_abs=q_abs, A_hand=A_hand,
                     T_rad_hand=T_rad_hand, T_rad_model=T_hot[1], T_el_model=T_hot[0])

# ---- 3. cold cases and heater sizing ----------------------------------------
md.append("## 3. Cold cases and heater sizing (aphelion, BOL, min albedo/IR)\n")
req = {}
rows = []
for c in COLD_CASES:
    T0 = Network(design, c).steady()[0]
    _, Q = required_heaters(design, c)
    req[c.name] = Q
    rows.append([c.name, f"{T0[0]:.1f} / {T0[2]:.1f} / {T0[3]:.1f}", f"{Q[0]:.1f} / {Q[1]:.1f} / {Q[2]:.1f}"])
table(["Case", "Unheated EL / BAT / PROP (C)", "Required heater EL / BAT / PROP (W)"], rows)
worst = np.max(list(req.values()), axis=0)
heaters = {}
rows = []
for name, w, n_el in zip(("EL", "BAT", "PROP"), worst, (4, 1, 1)):
    R_each = np.floor(V_MIN**2 / (HEATER_FACTOR * w / n_el) * 2) / 2
    R = R_each / n_el
    heaters[name] = dict(required=w, design=HEATER_FACTOR * w, elements=n_el, R_element=R_each,
                         P24=V_MIN**2 / R, P28=V_NOM**2 / R, P32=V_MAX**2 / R)
    h = heaters[name]
    rows.append([name, f"{w:.1f}", f"{h['design']:.1f}", f"{n_el} x {R_each:.1f}",
                 f"{h['P24']:.1f} / {h['P28']:.1f} / {h['P32']:.1f}"])
table(["Heater", "Worst required (W)", "x1.25 (W)", "Elements x R (ohm)", "P at 24 / 28 / 32 V (W)"], rows)
design.P_inst = tuple(heaters[n]["P24"] for n in ("EL", "BAT", "PROP"))
out["cold"] = {k: v.tolist() for k, v in req.items()}
out["heaters"] = heaters

# ---- 4. transients --------------------------------------------------------
md.append("## 4. Orbital transients (RK4, dt = 10 s, 240 h run, stats over last 120 h)\n")
transients = {
    "hot_op": hot_case(beta_hot),
    "peak_payload_sar": Case("SAR pass, elliptical", beta=0, mode="sar", **ELLIP, **HOT),
    "cold_op_eclipse": COLD_CASES[0],
    "safe_58min_eclipse": COLD_CASES[5],
    "hot_survival": Case("Hot survival (safe)", beta=beta_hot, mode="safe", **HOT),
}
rows = []
out["transients"] = {}
for key, c in transients.items():
    net = Network(design, c)
    T0, _ = required_heaters(design, c)
    n_orb = int(np.ceil(240 * 3600 / net.T))
    tt, hist, hon, st, bal = net.transient(T0, n_orb, stat_window=120 * 3600)
    keep = tt >= tt[-1] - 2 * net.T
    idx = ((tt[keep] % net.T) / net.T * len(net.t)).astype(int) % len(net.t)
    data = np.column_stack([(tt[keep] - tt[keep][0]) / 60, hist[keep],
                            np.vstack([hon, hon[-1:]])[keep], net.ecl[idx]])
    np.savetxt(os.path.join(RES, "transients", f"{key}.csv"), data, delimiter=",", fmt="%.4f",
               header="time_min,T_EL,T_RAD,T_BAT,T_PROP,T_BUS,heater_EL,heater_BAT,heater_PROP,eclipse",
               comments="")
    out["transients"][key] = dict(name=c.name, Tmin=st["Tmin"].tolist(), Tmax=st["Tmax"].tolist(),
                                  duty=st["duty"].tolist(), Pavg=st["Pavg"].tolist(),
                                  balance=bal, period_min=net.T / 60)
    rows.append([c.name, *(f"{st['Tmin'][i]:.1f} / {st['Tmax'][i]:.1f}" for i in range(5)),
                 " / ".join(f"{100*x:.0f}%" for x in st["duty"]), f"{st['Pavg'].sum():.1f}",
                 f"{bal['residual_pct']:+.4f}%"])
table(["Case", "EL", "RAD", "BAT", "PROP", "BUS", "Heater duty EL/BAT/PROP", "Avg heater (W)", "Energy residual"], rows)
b = out["transients"]["cold_op_eclipse"]["balance"]
md.append(f"Energy balance over the last orbit of the eclipse case (MJ): internal + heaters {b['E_internal']/1e6:.3f}, "
          f"absorbed {b['E_absorbed']/1e6:.3f}, radiated {b['E_radiated']/1e6:.3f}, MLI loss {b['E_mli']/1e6:.3f}, "
          f"stored {b['dU']/1e6:.3f}, residual {b['residual']:.1f} J.\n")
net = Network(design, COLD_CASES[0])
T0, _ = required_heaters(design, COLD_CASES[0])
s10 = net.transient(T0, 4, dt=10.0)[3]
s5 = net.transient(T0, 4, dt=5.0)[3]
dt_err = float(np.max(np.abs(np.r_[s10["Tmin"] - s5["Tmin"], s10["Tmax"] - s5["Tmax"]])))
md.append(f"Time-step check: dt = 10 s vs 5 s, max difference {dt_err:.3f} K.\n")
out["dt_check_K"] = dt_err

# ---- 5. margins ------------------------------------------------------------
md.append("## 5. Temperature margins (M_hot = T_max,allow - T_pred ; M_cold = T_pred - T_min,allow)\n")
hot_max = out["transients"]["hot_op"]["Tmax"]
cold_min = np.min([out["transients"][k]["Tmin"] for k in ("cold_op_eclipse", "peak_payload_sar")], axis=0)
safe_min = out["transients"]["safe_58min_eclipse"]["Tmin"]
rows, margins = [], {}
for name, i in (("EL", 0), ("BAT", 2), ("PROP", 3)):
    lo, hi, slo, shi = LIMITS[name]
    m = dict(hot=hot_max[i], m_hot=hi - hot_max[i], cold=cold_min[i], m_cold=cold_min[i] - lo,
             safe=safe_min[i], m_safe=safe_min[i] - slo, limits=LIMITS[name])
    margins[name] = m
    status = "PASS" if m["m_hot"] >= MARGIN and m["m_cold"] > 0 and m["m_safe"] > 0 else "CHECK"
    rows.append([name, f"{m['hot']:.1f}", hi, f"{m['m_hot']:+.1f}", f"{m['cold']:.1f}", lo,
                 f"{m['m_cold']:+.1f}", f"{m['safe']:.1f}", slo, status])
table(["Node", "Hot pred", "Op max", "Hot margin", "Cold pred", "Op min", "Cold margin",
       "Safe-mode min", "Surv min", "Status"], rows)
out["margins"] = margins

# ---- 6. design trades -------------------------------------------------------
md.append("## 6. Design trades (each concept re-sized to the same hot-case margins)\n")
PANEL, FINISH, LOUVER = 3.0, {"OSR": 0.60, "AgTeflon": 0.35, "Z93": 0.25}, 4.5   # kg/m2


def evaluate(d, label):
    T_hot = hot_sweep(d).max(0)
    q = {c.name: required_heaters(d, c)[1] for c in COLD_CASES}
    op = np.max([q[c.name] for c in COLD_CASES[:4]], axis=0).sum()
    surv = np.max(list(q.values()), axis=0).sum()
    T_cold = Network(d, COLD_CASES[0]).steady()[0]
    mass = (d.A_rad + d.A_brad) * (PANEL + FINISH[d.rad_coating]) + (LOUVER * d.A_rad if d.louvers else 0.0)
    return dict(concept=label, A_rad=d.A_rad, A_brad=d.A_brad, EL_hot=T_hot[0], PROP_hot=T_hot[3],
                BAT_hot=T_hot[2], EL_cold_unheated=T_cold[0], heater_op=op, heater_surv=surv, mass=mass)


def resized(**changes):
    d = size_radiators(design, **changes)
    d.A_rad = round(d.A_rad + 0.10, 2)
    return d


trades = {
    "T1 radiator finish": [evaluate(design, "A: Quartz-mirror OSR"),
                           evaluate(resized(rad_coating="AgTeflon"), "B: Silvered Teflon"),
                           evaluate(resized(rad_coating="Z93"), "C: Z93 white paint")],
    "T2 radiator location": [evaluate(design, "A: Anti-sun side (-Y)"),
                             evaluate(resized(rad_face="zenith"), "B: Zenith (-Z)"),
                             evaluate(resized(rad_face="ram"), "C: Ram (+X)")],
    "T3 cold-case control": [evaluate(design, "A: Fixed radiator + heaters"),
                             evaluate(resized(louvers=True), "B: Louvered radiator + heaters"),
                             evaluate(Design(**{**design.__dict__, "A_rad": round(1.5 * design.A_rad, 2)}),
                                      "C: Oversized (1.5x) radiator + heaters")],
}
d_louver = resized(louvers=True)
d_stuck = Design(**{**d_louver.__dict__, "louver_eps_open": d_louver.louver_eps_closed})
T_stuck = Network(d_stuck, hot_case(beta_hot)).steady()[0]
T_fail = Network(design, COLD_CASES[-1]).steady()[0]
for name, rows_ in trades.items():
    md.append(f"### {name}\n")
    table(["Concept", "A_rad (m2)", "A_bat (m2)", "Hot EL", "Hot PROP", "Hot BAT", "Cold-op EL unheated",
           "Op heater (W)", "Survival heater (W)", "Radiator mass (kg)"],
          [[r["concept"], f"{r['A_rad']:.2f}", f"{r['A_brad']:.2f}", f"{r['EL_hot']:.1f}", f"{r['PROP_hot']:.1f}",
            f"{r['BAT_hot']:.1f}", f"{r['EL_cold_unheated']:.1f}", f"{r['heater_op']:.0f}",
            f"{r['heater_surv']:.0f}", f"{r['mass']:.1f}"] for r in rows_])
md.append(f"Failure checks: louvers stuck closed in the hot case -> EL {T_stuck[0]:.1f} C; "
          f"heaters failed off in cruise -> EL {T_fail[0]:.1f} C.\n")
out["trades"] = trades
out["failure_checks"] = dict(louver_stuck_closed_EL=T_stuck[0], heaters_off_cruise_EL=T_fail[0])
out["design"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in design.__dict__.items()}

with open(os.path.join(RES, "results.json"), "w") as f:
    json.dump(out, f, indent=1, default=float)
with open(os.path.join(RES, "summary.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(md))
print("results written to", RES)
