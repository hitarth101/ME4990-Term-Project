# Results summary

## 1. Model checks

| Check | Model | Reference |
|---|---|---|
| Nadir view factor, 250 km | 0.92223 | (R/r)^2 = 0.92223 |
| Side-face view factor, 250 km | 0.3248 | - |
| Eclipse, 250 km, beta 0 (min) | 37.71 | P asin(R/r)/pi = 37.69 |
| Black sphere 607 km, beta 0 perihelion: max / min (C) | 121.2 / -105.7 | Gilmore T2.11: 122 / -105 |
| Black sphere 607 km, beta 0 aphelion: max / min (C) | 117.4 / -105.7 | Gilmore T2.11: 119 / -105 |
| Black sphere 607 km, beta 90 perihelion: max / min (C) | 67.1 / 67.1 | Gilmore T2.11: 67 / 67 |
| Black sphere 607 km, beta 90 aphelion: max / min (C) | 63.8 / 63.8 | Gilmore T2.11: 64 / 64 |

Elliptical orbit (294 x 8540 km): period 3.28 h, longest eclipse 58.1 min.

## 2. Radiator sizing (hot case: perihelion, EOL, max albedo/IR)

| beta (deg) | EL | RAD | BAT | PROP | BUS |
|---|---|---|---|---|---|
| 0 | 24.5 | 13.0 | 20.0 | 23.2 | 23.2 |
| 15 | 26.4 | 14.5 | 21.2 | 25.5 | 25.5 |
| 30 | 26.5 | 14.4 | 20.9 | 25.9 | 25.9 |
| 45 | 25.1 | 12.9 | 19.2 | 24.6 | 24.6 |
| 60 | 22.5 | 10.3 | 16.0 | 22.0 | 22.0 |
| 74 | 24.5 | 11.6 | 13.4 | 24.7 | 24.7 |
| 90 | 16.0 | 4.0 | 6.2 | 15.0 | 15.0 |

Minimum area 1.10 m2, selected 1.20 m2 (+0.10 m2 allowance); battery radiator 0.05 m2. Hot-case beta = 30 deg.

Hand check: q_abs = 83.1 W/m2, single-node radiator T = 11.0 C (model 14.4 C); area to reject 255 W = 1.12 m2.

## 3. Cold cases and heater sizing (aphelion, BOL, min albedo/IR)

| Case | Unheated EL / BAT / PROP (C) | Required heater EL / BAT / PROP (W) |
|---|---|---|
| Ops: circular, β = 0° | 7.0 / 4.1 / 3.9 | 0.0 / 2.9 / 5.4 |
| Ops: circular, β = 90° | 3.1 / -5.0 / 0.1 | 0.0 / 7.7 / 7.3 |
| Ops: elliptical, β = 0° | 1.2 / -5.4 / -3.0 | 0.0 / 7.7 / 9.1 |
| Ops: elliptical, β = 90° | 1.3 / -9.4 / -2.5 | 0.0 / 9.8 / 8.6 |
| Safe: circular, β = 90° | -24.5 / -28.2 / -24.5 | 94.7 / 12.3 / 10.3 |
| Safe: elliptical, 58-min eclipse | -36.2 / -36.6 / -36.6 | 136.9 / 12.5 / 11.2 |
| Safe: elliptical, β = 90° | -36.0 / -39.7 / -35.9 | 135.0 / 14.1 / 10.9 |
| Safe: cruise near 1 AU | -51.9 / -59.4 / -54.3 | 181.4 / 17.4 / 12.6 |

| Heater | Worst required (W) | x1.25 (W) | Elements x R (ohm) | P at 24 / 28 / 32 V (W) |
|---|---|---|---|---|
| EL | 181.4 | 226.8 | 4 x 10.0 | 230.4 / 313.6 / 409.6 |
| BAT | 17.4 | 21.8 | 1 x 26.0 | 22.2 / 30.2 / 39.4 |
| PROP | 12.6 | 15.7 | 1 x 36.5 | 15.8 / 21.5 / 28.1 |

## 4. Orbital transients (RK4, dt = 10 s, 240 h run, stats over last 120 h)

| Case | EL | RAD | BAT | PROP | BUS | Heater duty EL/BAT/PROP | Avg heater (W) | Energy residual |
|---|---|---|---|---|---|---|---|---|
| Hot op: circular, β = 30° | 26.1 / 26.9 | 12.6 / 17.1 | 20.8 / 21.1 | 25.9 / 25.9 | 24.7 / 27.1 | 0% / 0% / 0% | 0.0 | -0.0000% |
| SAR pass, elliptical | 11.7 / 15.8 | -1.1 / 6.2 | 10.0 / 13.0 | 13.0 / 16.0 | 10.0 / 13.2 | 0% / 16% / 13% | 5.5 | -0.0000% |
| Ops: circular, β = 0° | 8.1 / 8.6 | -3.1 / -0.9 | 10.0 / 13.0 | 13.0 / 16.0 | 4.9 / 6.3 | 0% / 17% / 43% | 10.4 | +0.0001% |
| Safe: elliptical, 58-min eclipse | -0.0 / 3.0 | -11.4 / -6.4 | 10.0 / 13.0 | 13.0 / 16.0 | -2.8 / -0.7 | 63% / 58% / 73% | 170.0 | +0.0007% |
| Hot survival (safe) | 5.7 / 6.6 | -4.1 / 0.7 | 10.0 / 13.0 | 13.0 / 16.0 | 7.3 / 9.9 | 0% / 22% / 31% | 9.8 | -0.0001% |

Energy balance over the last orbit of the eclipse case (MJ): internal + heaters 1.605, absorbed 0.362, radiated 1.697, MLI loss 0.184, stored 0.085, residual 2.1 J.

Time-step check: dt = 10 s vs 5 s, max difference 0.011 K.

## 5. Temperature margins (M_hot = T_max,allow - T_pred ; M_cold = T_pred - T_min,allow)

| Node | Hot pred | Op max | Hot margin | Cold pred | Op min | Cold margin | Safe-mode min | Surv min | Status |
|---|---|---|---|---|---|---|---|---|---|
| EL | 26.9 | 50 | +23.1 | 8.1 | -5 | +13.1 | -0.0 | -15 | PASS |
| BAT | 21.1 | 35 | +13.9 | 10.0 | 5 | +5.0 | 10.0 | -10 | PASS |
| PROP | 25.9 | 40 | +14.1 | 13.0 | 10 | +3.0 | 13.0 | 5 | PASS |

## 6. Design trades (each concept re-sized to the same hot-case margins)

### T1 radiator finish

| Concept | A_rad (m2) | A_bat (m2) | Hot EL | Hot PROP | Hot BAT | Cold-op EL unheated | Op heater (W) | Survival heater (W) | Radiator mass (kg) |
|---|---|---|---|---|---|---|---|---|---|
| A: Quartz-mirror OSR | 1.20 | 0.05 | 26.5 | 25.9 | 21.2 | 7.0 | 19 | 211 | 4.5 |
| B: Silvered Teflon | 1.30 | 0.05 | 26.2 | 25.4 | 22.3 | 5.7 | 20 | 222 | 4.5 |
| C: Z93 white paint | 1.30 | 0.05 | 26.8 | 25.7 | 22.5 | 3.5 | 50 | 264 | 4.4 |

### T2 radiator location

| Concept | A_rad (m2) | A_bat (m2) | Hot EL | Hot PROP | Hot BAT | Cold-op EL unheated | Op heater (W) | Survival heater (W) | Radiator mass (kg) |
|---|---|---|---|---|---|---|---|---|---|
| A: Anti-sun side (-Y) | 1.20 | 0.05 | 26.5 | 25.9 | 21.2 | 7.0 | 19 | 211 | 4.5 |
| B: Zenith (-Z) | 1.55 | 0.07 | 28.5 | 26.4 | 22.0 | -0.4 | 91 | 275 | 5.8 |
| C: Ram (+X) | 1.80 | 0.08 | 29.7 | 27.2 | 23.3 | 1.4 | 101 | 317 | 6.8 |

### T3 cold-case control

| Concept | A_rad (m2) | A_bat (m2) | Hot EL | Hot PROP | Hot BAT | Cold-op EL unheated | Op heater (W) | Survival heater (W) | Radiator mass (kg) |
|---|---|---|---|---|---|---|---|---|---|
| A: Fixed radiator + heaters | 1.20 | 0.05 | 26.5 | 25.9 | 21.2 | 7.0 | 19 | 211 | 4.5 |
| B: Louvered radiator + heaters | 1.35 | 0.05 | 27.1 | 26.6 | 21.4 | 20.2 | 4 | 26 | 11.1 |
| C: Oversized (1.5x) radiator + heaters | 1.80 | 0.05 | 13.1 | 12.6 | 12.8 | -6.7 | 94 | 309 | 6.7 |

Failure checks: louvers stuck closed in the hot case -> EL 74.0 C; heaters failed off in cruise -> EL -51.9 C.
