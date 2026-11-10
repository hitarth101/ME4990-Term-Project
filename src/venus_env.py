"""Venus orbital environment: orbit geometry, eclipse, and face fluxes.

Sun direction is +x. Orbit normal h = (sin b, 0, cos b) so beta = asin(s.h).
Spacecraft is nadir-pointing with a yaw flip so the sun is always on the +Y side.
Albedo is diffuse with cosine fall-off from the subsolar point (Gilmore Ch. 2).
"""
import numpy as np

SIGMA = 5.670374419e-8
R_V = 6051.8e3            # m   (NSSDC)
MU_V = 3.24859e14         # m^3/s^2 (NSSDC)

S_1AU = 1367.5            # W/m^2, Gilmore Eq. 2.9
S_PERI = S_1AU / 0.7184**2
S_MEAN = S_1AU / 0.7233**2
S_APH = S_1AU / 0.7282**2
ALB_NOM, ALB_HOT, ALB_COLD = 0.80, 0.82, 0.78    # Gilmore Table 2.9
IR_NOM, IR_HOT, IR_COLD = 153.0, 178.4, 138.5    # Gilmore Tables 2.9, 2.10

FACES = ("nadir", "zenith", "ram", "wake", "sun", "anti")


def orbit_basis(beta_deg):
    b = np.radians(beta_deg)
    h = np.array([np.sin(b), 0.0, np.cos(b)])
    e1 = np.array([np.cos(b), 0.0, -np.sin(b)])     # orbit noon direction
    e2 = np.array([0.0, 1.0, 0.0])
    return h, e1, e2


def circular_period(alt):
    return 2 * np.pi * np.sqrt((R_V + alt)**3 / MU_V)


def orbit_states(alt_p, alt_a, beta_deg, omega_deg=0.0, n=720):
    """Position vectors at n equal time steps. omega = periapsis angle from orbit noon."""
    rp, ra = R_V + alt_p, R_V + alt_a
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    T = 2 * np.pi * np.sqrt(a**3 / MU_V)
    t = np.linspace(0.0, T, n, endpoint=False)
    M = 2 * np.pi * t / T
    E = M.copy()
    for _ in range(50):
        E -= (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e * np.cos(E))
    u = np.radians(omega_deg) + nu
    _, e1, e2 = orbit_basis(beta_deg)
    r_vec = r[:, None] * (np.cos(u)[:, None] * e1 + np.sin(u)[:, None] * e2)
    return t, r_vec, T


def in_eclipse(r_vec):
    s = np.array([1.0, 0.0, 0.0])
    along = r_vec @ s
    perp = np.linalg.norm(r_vec - along[:, None] * s, axis=1)
    return (along < 0) & (perp < R_V)


def face_normals(r_vec, beta_deg):
    h, _, _ = orbit_basis(beta_deg)
    rh = r_vec / np.linalg.norm(r_vec, axis=1)[:, None]
    ram = np.cross(h, rh)
    ram /= np.linalg.norm(ram, axis=1)[:, None]
    hh = np.broadcast_to(h, rh.shape)
    return {"nadir": -rh, "zenith": rh, "ram": ram, "wake": -ram, "sun": hh.copy(), "anti": -hh}


def _cap_grid(r_mag, R, n_u=160, n_phi=160):
    # polar grid on the visible cap, clustered near the sub-spacecraft point
    lam_h = np.arccos(R / r_mag)
    u = (np.arange(n_u) + 0.5) / n_u
    lam = lam_h * u**2
    dlam = 2 * lam_h * u / n_u
    dphi = 2 * np.pi / n_phi
    phi = (np.arange(n_phi) + 0.5) * dphi
    L, P = np.meshgrid(lam, phi, indexing="ij")
    pts = np.stack([np.sin(L) * np.cos(P), np.sin(L) * np.sin(P), np.cos(L)], -1)
    dA = R**2 * np.sin(L) * dlam[:, None] * dphi
    return pts.reshape(-1, 3), dA.reshape(-1)


def _rotation_to(z):
    z = z / np.linalg.norm(z)
    a = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = np.cross(a, z); x /= np.linalg.norm(x)
    return np.stack([x, np.cross(z, x), z], axis=1)


def planet_factors(r_vec, normals, R=R_V, sphere=False):
    """IR view factor F and albedo factor Fa per face (and for a small sphere if asked)."""
    s = np.array([1.0, 0.0, 0.0])
    n = len(r_vec)
    F = {k: np.zeros(n) for k in normals}
    Fa = {k: np.zeros(n) for k in normals}
    Fs, Fas = np.zeros(n), np.zeros(n)
    grids = {}
    for i, rv in enumerate(r_vec):
        key = round(np.linalg.norm(rv), -2)
        if key not in grids:
            grids[key] = _cap_grid(np.linalg.norm(rv), R)
        local, dA = grids[key]
        p_hat = local @ _rotation_to(rv).T
        d = R * p_hat - rv
        dist = np.linalg.norm(d, axis=1)
        dh = d / dist[:, None]
        cos_s = np.clip(-(dh * p_hat).sum(1), 0, None)
        lit = np.clip(p_hat @ s, 0, None)
        w = cos_s * dA / (np.pi * dist**2)
        for k, nv in normals.items():
            cos_f = np.clip(dh @ nv[i], 0, None)
            F[k][i] = (w * cos_f).sum()
            Fa[k][i] = (w * cos_f * lit).sum()
        Fs[i], Fas[i] = 0.25 * w.sum(), 0.25 * (w * lit).sum()
    return (F, Fa, Fs, Fas) if sphere else (F, Fa)


def face_fluxes(alt_p, alt_a, beta_deg, omega_deg=0.0, n=720, S=S_MEAN, albedo=ALB_NOM, q_ir=IR_NOM):
    """Incident solar, albedo and IR flux (W/m^2) on each face over one orbit."""
    t, r_vec, T = orbit_states(alt_p, alt_a, beta_deg, omega_deg, n)
    ecl = in_eclipse(r_vec)
    nrm = face_normals(r_vec, beta_deg)
    F, Fa = planet_factors(r_vec, nrm)
    out = {"t": t, "period": T, "eclipse": ecl, "q_sun": {}, "q_alb": {}, "q_ir": {}}
    for k, nv in nrm.items():
        out["q_sun"][k] = S * np.clip(nv[:, 0], 0, None) * (~ecl)
        out["q_alb"][k] = albedo * S * Fa[k]
        out["q_ir"][k] = q_ir * F[k]
    return out


def eclipse_duration(alt_p, alt_a, beta_deg, omega_deg=0.0, n=4000):
    t, r_vec, T = orbit_states(alt_p, alt_a, beta_deg, omega_deg, n)
    return in_eclipse(r_vec).mean() * T, T
