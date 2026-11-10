"""Five-node lumped thermal network of the Venus orbiter bus.

Nodes: EL (electronics), RAD (-Y radiator), BAT (battery + its radiator patch),
       PROP (hydrazine tank), BUS (structure behind MLI)

C_i dT_i/dt = sum of conduction, radiation and MLI exchange + internal/heater power
"""
import copy
from dataclasses import dataclass

import numpy as np

import venus_env as ve

SIG = ve.SIGMA
K0 = 273.15
NODES = ("EL", "RAD", "BAT", "PROP", "BUS")

LX, LY, LZ = 1.6, 1.4, 1.8          # bus dimensions, m
FACE_AREA = {"ram": LY * LZ, "wake": LY * LZ, "sun": LX * LZ, "anti": LX * LZ,
             "nadir": LX * LY, "zenith": LX * LY}

# (alpha_BOL, alpha_EOL at 4 yr, emittance) - Gilmore Table 4.1, Figs 4.7, 4.9, 4.10
COATINGS = {"OSR": (0.08, 0.20, 0.80), "AgTeflon": (0.09, 0.24, 0.78), "Z93": (0.20, 0.48, 0.92)}

# internal dissipation by mode: electronics node W, battery (sunlit, eclipse) W
POWER = {
    "gravity": dict(el=255.0, bat=(8.0, 15.0)),
    "sar_on": dict(el=475.0, bat=(8.0, 15.0)),
    "sar_off": dict(el=240.0, bat=(8.0, 15.0)),
    "safe": dict(el=130.0, bat=(5.0, 10.0)),
    "cruise": dict(el=130.0, bat=(3.0, 3.0)),
}
SAR_WINDOW_MIN = 40.0


@dataclass
class Design:
    A_rad: float = 1.20
    A_brad: float = 0.05
    rad_coating: str = "OSR"
    rad_face: str = "anti"
    louvers: bool = False
    louver_eps_closed: float = 0.15
    louver_eps_open: float = 0.70
    louver_T_closed: float = 5.0
    louver_T_open: float = 15.0
    mli_alpha_bol: float = 0.32       # Beta cloth
    mli_alpha_eol: float = 0.55
    mli_eps: float = 0.86
    e_star: float = 0.03
    C: tuple = (144e3, 13.5e3, 16.5e3, 257e3, 162e3)    # J/K
    G_el_rad: float = 20.0            # W/K
    G_rad_bus: float = 1.5
    G_bat_bus: float = 0.3
    G_prop_bus: float = 0.4
    epsA_el_bus: float = 0.865 * 3.6  # m^2
    A_box_footprint: float = 1.0
    eps_rad_in: float = 0.86
    estarA_bat_bus: float = 0.03 * 0.6
    estarA_prop_bus: float = 0.03 * 2.0
    set_el: tuple = (0.0, 3.0)        # heater ON/OFF, C
    set_bat: tuple = (10.0, 13.0)
    set_prop: tuple = (13.0, 16.0)
    P_inst: tuple = (0.0, 0.0, 0.0)   # installed heater power EL, BAT, PROP (W)

    def rad_props(self, eol):
        a_bol, a_eol, eps = COATINGS[self.rad_coating]
        return (a_eol if eol else a_bol), eps


@dataclass
class Case:
    name: str
    alt_p: float = 250e3
    alt_a: float = 250e3
    beta: float = 0.0
    omega: float = 0.0
    S: float = ve.S_MEAN
    albedo: float = ve.ALB_NOM
    q_ir: float = ve.IR_NOM
    eol: bool = True
    mode: str = "gravity"
    deep_space: bool = False
    n_orbit_pts: int = 360


_env_cache = {}


def environment(case):
    key = (case.alt_p, case.alt_a, case.beta, case.omega, case.S, case.albedo,
           case.q_ir, case.deep_space, case.n_orbit_pts)
    if key in _env_cache:
        return _env_cache[key]
    if case.deep_space:       # cruise: sun on +Y, no planet
        t = np.linspace(0, 5400, case.n_orbit_pts, endpoint=False)
        z = np.zeros_like(t)
        q = {f: dict(sun=z.copy(), alb=z.copy(), ir=z.copy()) for f in ve.FACES}
        q["sun"]["sun"] = np.full_like(t, case.S)
        out = (t, 5400.0, np.zeros_like(t, bool), q)
    else:
        e = ve.face_fluxes(case.alt_p, case.alt_a, case.beta, case.omega, case.n_orbit_pts,
                           case.S, case.albedo, case.q_ir)
        q = {f: dict(sun=e["q_sun"][f], alb=e["q_alb"][f], ir=e["q_ir"][f]) for f in ve.FACES}
        out = (e["t"], e["period"], e["eclipse"], q)
    _env_cache[key] = out
    return out


def internal_power(case, t, period, ecl):
    if case.mode == "sar":     # SAR on +/-20 min around periapsis (t = 0)
        on = np.minimum(t, period - t) <= SAR_WINDOW_MIN * 30
        el = np.where(on, POWER["sar_on"]["el"], POWER["sar_off"]["el"])
        pb = POWER["sar_on"]["bat"]
    else:
        el = np.full_like(t, POWER[case.mode]["el"])
        pb = POWER[case.mode]["bat"]
    return el, np.where(ecl, pb[1], pb[0])


class Network:
    def __init__(self, design, case):
        self.d = d = copy.deepcopy(design)
        self.c = case
        self.t, self.T, self.ecl, q = environment(case)
        self.Pel, self.Pbat = internal_power(case, self.t, self.T, self.ecl)
        self.a_rad, self.e_rad = d.rad_props(case.eol)
        a_mli = d.mli_alpha_eol if case.eol else d.mli_alpha_bol
        rf = d.rad_face
        self.qabs_rad = self.a_rad * (q[rf]["sun"] + q[rf]["alb"]) + self.e_rad * q[rf]["ir"]
        self.mli_faces = []
        for f in ve.FACES:
            A = FACE_AREA[f] - (d.A_rad + d.A_brad if f == rf else 0.0)
            if A < -1e-9:
                raise ValueError("radiator larger than face")
            qa = a_mli * (q[f]["sun"] + q[f]["alb"]) + d.mli_eps * q[f]["ir"]
            self.mli_faces.append((A, qa))
        self.A_rad_in = max(d.A_rad - d.A_box_footprint, 0.0)
        self.Qabs_rad = d.A_rad * self.qabs_rad
        self.Qabs_brad = d.A_brad * self.qabs_rad
        self.Cvec = np.array(d.C)

    def rad_eps(self, T_rad_C):
        d = self.d
        if not d.louvers:
            return self.e_rad, 1.0
        f = np.clip((T_rad_C - d.louver_T_closed) / (d.louver_T_open - d.louver_T_closed), 0, 1)
        e = d.louver_eps_closed + f * (d.louver_eps_open - d.louver_eps_closed)
        return e, e / self.e_rad

    def _mli_loss(self, T_bus, pick):
        d = self.d
        k = d.e_star / (d.mli_eps + d.e_star)
        return sum(k * A * (d.mli_eps * SIG * T_bus**4 - pick(qa)) for A, qa in self.mli_faces)

    def flows(self, TK, k, Qh, mean=False):
        """Net heat into each node (W) at orbit index k, or orbit-averaged."""
        d = self.d
        T_el, T_rad, T_bat, T_pr, T_bus = TK
        pick = (lambda a: a.mean()) if mean else (lambda a: a[k])
        e_r, scale = self.rad_eps(T_rad - K0)
        q_er = d.G_el_rad * (T_el - T_rad)
        q_eb = SIG * d.epsA_el_bus * (T_el**4 - T_bus**4)
        q_rb = d.G_rad_bus * (T_rad - T_bus) + SIG * d.eps_rad_in * self.A_rad_in * (T_rad**4 - T_bus**4)
        q_bb = d.G_bat_bus * (T_bat - T_bus) + SIG * d.estarA_bat_bus * (T_bat**4 - T_bus**4)
        q_pb = d.G_prop_bus * (T_pr - T_bus) + SIG * d.estarA_prop_bus * (T_pr**4 - T_bus**4)
        return np.array([
            pick(self.Pel) + Qh[0] - q_er - q_eb,
            q_er + scale * pick(self.Qabs_rad) - e_r * SIG * d.A_rad * T_rad**4 - q_rb,
            pick(self.Pbat) + Qh[1] + pick(self.Qabs_brad) - self.e_rad * SIG * d.A_brad * T_bat**4 - q_bb,
            Qh[2] - q_pb,
            q_eb + q_rb + q_bb + q_pb - self._mli_loss(T_bus, pick),
        ])

    def energy_terms(self, TK, k, Qh):
        d = self.d
        e_r, scale = self.rad_eps(TK[1] - K0)
        internal = self.Pel[k] + self.Pbat[k] + sum(Qh)
        absorbed = scale * self.Qabs_rad[k] + self.Qabs_brad[k]
        rejected = e_r * SIG * d.A_rad * TK[1]**4 + self.e_rad * SIG * d.A_brad * TK[2]**4
        return internal, absorbed, rejected, self._mli_loss(TK[4], lambda a: a[k])

    def steady(self, heater_nodes=()):
        """Orbit-average steady state. Nodes in heater_nodes are held at their
        thermostat ON temperature and the heater power becomes the unknown."""
        d = self.d
        hmap = {"EL": (0, 0, d.set_el[0]), "BAT": (2, 1, d.set_bat[0]), "PROP": (3, 2, d.set_prop[0])}
        active = [hmap[n] for n in heater_nodes]
        x = np.full(5, 290.0)
        for ni, _, _ in active:
            x[ni] = 10.0

        def unpack(x):
            TK, Qh = x.copy(), np.zeros(3)
            for ni, hi, Ts in active:
                Qh[hi], TK[ni] = x[ni], Ts + K0
            return TK, Qh

        def res(x):
            TK, Qh = unpack(x)
            return self.flows(TK, 0, Qh, mean=True)

        for _ in range(200):
            r = res(x)
            J = np.column_stack([(res(x + dx) - r) / 1e-3 for dx in np.eye(5) * 1e-3])
            step = np.linalg.solve(J, -r)
            x += np.clip(step, -50, 50)
            if np.abs(step).max() < 1e-7 and np.abs(r).max() < 1e-6:
                break
        TK, Qh = unpack(x)
        return TK - K0, Qh, self.flows(TK, 0, Qh, mean=True)

    def transient(self, T0_C, n_orbits, dt=10.0, stat_window=None):
        """RK4 with bang-bang thermostats. Stats over the last stat_window seconds."""
        d = self.d
        TK = np.array(T0_C, float) + K0
        P = np.array(d.P_inst)
        sets = [(0, d.set_el), (2, d.set_bat), (3, d.set_prop)]
        on = [False] * 3
        n_steps = int(round(n_orbits * self.T / dt))
        npts = len(self.t)
        hist = np.zeros((n_steps + 1, 5)); hist[0] = TK - K0
        hon = np.zeros((n_steps, 3))
        tt = np.arange(n_steps + 1) * dt
        E = np.zeros(4)
        t_last = (n_orbits - 1) * self.T
        U0 = None
        for s in range(n_steps):
            for h, (ni, (T_on, T_off)) in enumerate(sets):
                Tc = TK[ni] - K0
                if not on[h] and Tc <= T_on:
                    on[h] = True
                elif on[h] and Tc >= T_off:
                    on[h] = False
            Qh = P * np.array(on)
            k = int((s * dt % self.T) / self.T * npts) % npts
            f = lambda T: self.flows(T, k, Qh) / self.Cvec
            k1 = f(TK); k2 = f(TK + 0.5 * dt * k1); k3 = f(TK + 0.5 * dt * k2); k4 = f(TK + dt * k3)
            if s * dt >= t_last:
                if U0 is None:
                    U0 = (self.Cvec * TK).sum()
                E += np.array(self.energy_terms(TK, k, Qh)) * dt
            TK = TK + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            hist[s + 1] = TK - K0
            hon[s] = on
        dU = (self.Cvec * TK).sum() - U0
        win = self.T if stat_window is None else stat_window
        last = tt >= tt[-1] - win
        stats = dict(Tmin=hist[last].min(0), Tmax=hist[last].max(0),
                     duty=hon[last[:-1]].mean(0), Pavg=(hon[last[:-1]] * P).mean(0))
        resid = E[0] + E[1] - E[2] - E[3] - dU
        balance = dict(E_internal=E[0], E_absorbed=E[1], E_radiated=E[2], E_mli=E[3], dU=dU,
                       residual=resid, residual_pct=100 * resid / (E[0] + E[1]))
        return tt, hist, hon, stats, balance
